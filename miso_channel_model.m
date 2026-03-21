function [channelState, params] = miso_channel_model()
%MISO_CHANNEL_MODEL 多用户 MISO 波导-自由空间复合信道建模
% 输出：
%   channelState.users               -> 用户三维坐标 [K,3]
%   channelState.X                   -> PA 纵向位置 [M,N]
%   channelState.waveguideFeedPoints -> 波导馈电点 [N,3]
%   channelState.paPositions         -> PA 三维位置 [3,M,N]
%   channelState.theta, phi          -> 每个 PA 的波束角 [M,N]
%   channelState.H                   -> 总复合信道矩阵 [NM,K]
%   channelState.details             -> 中间变量结构体
%   params                           -> 系统参数结构体

    params = defaultParameters();

    % 用户位置与 PA 位置
    userXY = generateRandomUsers(params);
    X = generatePALayout(params);
    [users, waveguideFeedPoints, paPositions] = buildGeometry(userXY, X, params);
    [theta, phi, steeringTargets] = generateBeamAngles(paPositions, users);

    % 示例：单条波导导波信道与自由空间信道
    g1 = guidedWaveChannel(X(:, 1), params);
    userIdx = 1;
    wgIdx = 1;
    [hFreeExample, localCoordsExample] = freeSpaceChannelForWaveguide( ...
        users(userIdx, :).', squeeze(paPositions(:, :, wgIdx)), ...
        theta(:, wgIdx), phi(:, wgIdx), params);

    % 所有用户的总复合信道
    [H, details] = compositeChannel(users, X, theta, phi, params);

    % 单个用户/单个 PA 的局部坐标示例
    k = 1; n = 1; m = 1;
    localCoordExample = globalToPALocal(users(k, :).', paPositions(:, m, n), theta(m, n), phi(m, n));

    channelState = struct();
    channelState.userXY = userXY;
    channelState.users = users;
    channelState.X = X;
    channelState.waveguideFeedPoints = waveguideFeedPoints;
    channelState.paPositions = paPositions;
    channelState.theta = theta;
    channelState.phi = phi;
    channelState.steeringTargets = steeringTargets;
    channelState.g1 = g1;
    channelState.hFreeExample = hFreeExample;
    channelState.localCoordsExample = localCoordsExample;
    channelState.localCoordExample = localCoordExample;
    channelState.H = H;
    channelState.details = details;
end

function params = defaultParameters()
%DEFAULTPARAMETERS 设置示例所需的系统参数
    params.N = 2;                 % 波导数量（可改）
    params.M = 3;                 % 每条波导的 PA 数量（可改）
    params.K = 3;                 % 用户数量（可改）
    params.NRF = 2;               % 基站 RF 链数量，支持至多 N_RF 条数据流
    params.Dx = 10;               % 波导在 x 方向覆盖宽度
    params.d = 4;                 % 波导部署高度
    params.randomSeed = 20260319; % 随机种子，保证结果可复现
    params.symbolSeed = 20260320; % 符号随机种子

    params.userRegionX = [0, params.Dx];
    params.userRegionY = [6, 20];
    params.paYOffsetRange = [1, 10];
    params.paPlacementMode = 'uniform'; % 可选：'uniform' / 'random'
    params.schedulingMode = 'first';    % 可选：'first' / 'random'
    params.beamformingMethod = 'MRT';   % 可选：'MRT' / 'ZF'
    params.totalTransmitPower = 1;      % 总发射功率约束 sum_k ||w_k||^2 <= P
    params.sigma2 = 1e-4;               % AWGN 噪声功率

    params.lambda = 0.01;         % 载波波长（m）
    params.nEff = 1.6;            % 波导有效折射率 n_eff
    params.alphaW = 0.08;         % 波导衰减系数 alpha_W（Np/m）

    params.alphaL = 0.96;         % LoS 存在系数/米
    params.nMedium = 1.0;         % 自由空间/介质折射率 n
    params.a = 0.45;              % 模式尺寸参数 a
    params.b = 0.30;              % 模式尺寸参数 b
    params.v = 1.1;               % 经验修正系数 v
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
% 输出 X 大小为 [M, N]，第 n 列对应第 n 条波导上的 y_{n,m}

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
                users(k, :).', squeeze(paPositions(:, :, n)), theta(:, n), phi(:, n), params);
            hk(idx) = details.g{n} .* hFreeKN;
            details.hFree{k, n} = hFreeKN;
            details.localCoords{k, n} = localCoordsKN;
        end
        H(:, k) = hk;
        details.hComposite{k} = hk;
    end
end

function localCoord = globalToPALocal(userPos, paPos, theta, phi)
%GLOBALTOPALOCAL 将全局坐标转换到 PA 局部坐标系
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
%RADIATIONPATTERN 计算局部坐标下的辐射模式值
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
