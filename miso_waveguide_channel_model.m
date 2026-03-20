function model = miso_waveguide_channel_model(config)
%MISO_WAVEGUIDE_CHANNEL_MODEL 复用几何与复合信道建模，供信号模型主函数调用。
% 输入：
%   config.params   : 系统参数结构体；缺省时使用默认参数
%   config.userXY   : [K,2] 用户平面坐标；缺省时随机生成
%   config.X        : [M,N] PA 纵向位置；缺省时自动生成
%   config.theta    : [M,N] PA 俯仰角；缺省时按最近用户自动生成
%   config.phi      : [M,N] PA 方位角；缺省时按最近用户自动生成
%
% 输出 model 结构体包含：
%   params, userXY, X, users, feedPoints, paPositions,
%   theta, phi, steeringTargets, H, details

    if nargin < 1
        config = struct();
    end

    params = getFieldOrDefault(config, 'params', defaultParameters());
    userXY = getFieldOrDefault(config, 'userXY', generateRandomUsers(params));
    X = getFieldOrDefault(config, 'X', generatePALayout(params));

    [users, feedPoints, paPositions] = buildGeometry(userXY, X, params);

    if isfield(config, 'theta') && isfield(config, 'phi')
        theta = config.theta;
        phi = config.phi;
        steeringTargets = [];
    else
        [theta, phi, steeringTargets] = generateBeamAngles(paPositions, users);
    end

    [H, details] = compositeChannel(users, X, theta, phi, params);

    model = struct();
    model.params = params;
    model.userXY = userXY;
    model.X = X;
    model.users = users;
    model.feedPoints = feedPoints;
    model.paPositions = paPositions;
    model.theta = theta;
    model.phi = phi;
    model.steeringTargets = steeringTargets;
    model.H = H;
    model.details = details;
end

function value = getFieldOrDefault(s, fieldName, defaultValue)
    if isfield(s, fieldName)
        value = s.(fieldName);
    else
        value = defaultValue;
    end
end

function params = defaultParameters()
%DEFAULTPARAMETERS 设置示例所需的系统参数
    params.N = 4;
    params.M = 8;
    params.K = 10;
    params.Dx = 10;
    params.d = 4;
    params.randomSeed = 20260319;

    params.userRegionX = [0, params.Dx];
    params.userRegionY = [6, 20];
    params.paYOffsetRange = [1, 10];
    params.paPlacementMode = 'uniform';

    params.lambda = 0.01;
    params.nEff = 1.6;
    params.alphaW = 0.08;

    params.alphaL = 0.5;
    params.nMedium = 1.0;
    params.a = 0.45;
    params.b = 0.30;
    params.v = 1.1;
end

function userXY = generateRandomUsers(params)
%GENERATERANDOMUSERS 根据用户数 K 和区域范围随机生成用户位置
    rng(params.randomSeed);
    x = params.userRegionX(1) + diff(params.userRegionX) * rand(params.K, 1);
    y = params.userRegionY(1) + diff(params.userRegionY) * rand(params.K, 1);
    userXY = [x, y];
end

function X = generatePALayout(params)
%GENERATEPALAYOUT 根据参数自动生成每条波导上的 PA 纵向坐标
    yMin = params.paYOffsetRange(1);
    yMax = params.paYOffsetRange(2);

    switch lower(params.paPlacementMode)
        case 'uniform'
            baseLine = linspace(yMin, yMax, params.M).';
            X = repmat(baseLine, 1, params.N);
        case 'random'
            rng(params.randomSeed + 1);
            X = yMin + (yMax - yMin) * rand(params.M, params.N);
            X = sort(X, 1, 'ascend');
        otherwise
            error('未知的 paPlacementMode: %s', params.paPlacementMode);
    end
end

function [users, feedPoints, paPositions] = buildGeometry(userXY, X, params)
%BUILDGEOMETRY 构造用户、波导馈电点和 PA 三维位置
    validateattributes(userXY, {'numeric'}, {'size', [params.K, 2]});
    validateattributes(X, {'numeric'}, {'size', [params.M, params.N]});

    users = [userXY, zeros(params.K, 1)];
    feedPoints = zeros(params.N, 3);
    paPositions = zeros(3, params.M, params.N);

    for n = 1:params.N
        xW = ((2 * n - 1) / (2 * params.N)) * params.Dx;
        feedPoints(n, :) = [xW, 0, params.d];

        for m = 1:params.M
            paPositions(:, m, n) = [xW; X(m, n); params.d];
        end
    end
end

function [theta, phi, targetUserIdx] = generateBeamAngles(paPositions, users)
%GENERATEBEAMANGLES 根据几何关系自动生成每个 PA 的俯仰角和方位角
    [~, M, N] = size(paPositions);
    K = size(users, 1);

    theta = zeros(M, N);
    phi = zeros(M, N);
    targetUserIdx = zeros(M, N);

    for n = 1:N
        for m = 1:M
            pnm = paPositions(:, m, n);
            distances = zeros(K, 1);

            for k = 1:K
                distances(k) = norm(users(k, :).'- pnm);
            end

            [~, kStar] = min(distances);
            targetUserIdx(m, n) = kStar;

            direction = users(kStar, :).'- pnm;
            direction = direction / norm(direction);

            theta(m, n) = acos(direction(3));
            phi(m, n) = atan2(direction(2), direction(1));
        end
    end
end

function g = guidedWaveChannel(xn, params)
%GUIDEDWAVECHANNEL 论文中的导波信道 g_n(x_n)
    lambdaG = params.lambda / params.nEff;
    xn = xn(:);

    amplitude = sqrt(1 / params.M) .* exp(-(params.alphaW / 2) .* xn);
    phase = exp(-1j * (2 * pi / lambdaG) .* xn);
    g = amplitude .* phase;
end

function [hFree, localCoords] = freeSpaceChannelForWaveguide(userPos, paPosWaveguide, thetaN, phiN, params)
%FREESPACECHANNELFORWAVEGUIDE 计算用户 k 与第 n 条波导之间的自由空间信道
    M = size(paPosWaveguide, 2);
    eta = params.lambda^2 / (4 * pi);
    hFree = zeros(M, 1);
    localCoords = zeros(3, M);

    for m = 1:M
        pnm = paPosWaveguide(:, m);
        localCoords(:, m) = globalToPALocal(userPos, pnm, thetaN(m), phiN(m));
        distance = norm(userPos - pnm);
        pattern = radiationPattern(localCoords(:, m), params);
        hFree(m) = sqrt(eta) * (params.alphaL ^ distance) * pattern;
    end
end

function [H, details] = compositeChannel(users, X, theta, phi, params)
%COMPOSITECHANNEL 计算所有用户的总复合信道 h_k(X)
    [~, feedPoints, paPositions] = buildGeometry(users(:, 1:2), X, params); %#ok<ASGLU>

    K = size(users, 1);
    H = zeros(params.N * params.M, K);
    details = struct();
    details.feedPoints = feedPoints;
    details.paPositions = paPositions;
    details.g = cell(params.N, 1);
    details.hFree = cell(K, params.N);
    details.hComposite = cell(K, 1);
    details.localCoords = cell(K, params.N);

    for n = 1:params.N
        details.g{n} = guidedWaveChannel(X(:, n), params);
    end

    for k = 1:K
        hk = zeros(params.N * params.M, 1);

        for n = 1:params.N
            idx = (n - 1) * params.M + (1:params.M);
            [hFreeKN, localCoordsKN] = freeSpaceChannelForWaveguide( ...
                users(k, :).', squeeze(paPositions(:, :, n)), ...
                theta(:, n), phi(:, n), params);

            hk(idx) = details.g{n} .* hFreeKN;
            details.hFree{k, n} = hFreeKN;
            details.localCoords{k, n} = localCoordsKN;
        end

        H(:, k) = hk;
        details.hComposite{k} = hk;
    end
end

function localCoord = globalToPALocal(userPos, paPos, theta, phi)
%GLOBALTOPALOCAL 全局坐标系到 PA 局部坐标系旋转
    delta = userPos - paPos;
    Rx = rotationMatrixX(theta - pi / 2);
    Rz = rotationMatrixZ(pi / 2 - phi);
    localCoord = Rx * Rz * delta;
end

function R = rotationMatrixX(theta)
%ROTATIONMATRIXX 绕 x 轴旋转矩阵
    R = [
        1,           0,            0;
        0,  cos(theta), -sin(theta);
        0,  sin(theta),  cos(theta)
    ];
end

function R = rotationMatrixZ(phi)
%ROTATIONMATRIXZ 绕 z 轴旋转矩阵
    R = [
        cos(phi), -sin(phi), 0;
        sin(phi),  cos(phi), 0;
              0,         0, 1
    ];
end

function value = radiationPattern(localCoord, params)
%RADIATIONPATTERN 计算 PA 局部坐标系下的高斯波束辐射模式
    x = localCoord(1);
    y = localCoord(2);
    z = localCoord(3);

    if y <= 0
        value = 0;
        return;
    end

    lambda = params.lambda;
    n = params.nMedium;
    k = 2 * pi / lambda;
    v = params.v;
    a = params.a;
    b = params.b;

    W1 = lambda * y / (pi * n * (v * a * lambda));
    W2 = lambda * y / (pi * n * (v * b * lambda));
    R1 = y;
    R2 = y;

    Theta1 = atan(lambda * y / (pi * n * (v * a * lambda)));
    Theta2 = atan(lambda * y / (pi * n * (v * b * lambda)));

    B = sqrt(2 / (pi * (v * a * lambda) * (v * b * lambda)));
    amplitude = sqrt((v * a * lambda) * (v * b * lambda) / (W1 * W2)) * B;
    quadraticEnvelope = exp(-((x ^ 2) / (W1 ^ 2) + (z ^ 2) / (W2 ^ 2)));
    phase = exp(-1j * k * n * ((x ^ 2) / (2 * R1) + (z ^ 2) / (2 * R2) + y) ...
        + 1j * (Theta1 + Theta2) / 2);

    value = amplitude * quadraticEnvelope * phase;
end
