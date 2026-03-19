%% 多用户 MISO 波导-自由空间复合信道建模示例
% 本脚本对应论文中的三个部分：
% 1) 用户位置、波导位置和 PA 位置建模；
% 2) 波导内导波信道 g_n(x_n)；
% 3) 自由空间信道 \tilde{h}_{k,n}(x_n)；
% 4) 总复合信道 h_k(X)；
% 5) 全局坐标系到 PA 局部坐标系的旋转变换。
%
% 说明：
% - 脚本可直接运行，MATLAB R2016b+ 支持脚本末尾局部函数。
% - 若使用 Octave，建议版本支持脚本局部函数；否则可将每个函数拆分为独立 .m 文件。

clear; clc;

%% =========================
%  1. 系统参数配置
%  =========================
params = defaultParameters();

% 用户位置（每一行为一个用户 [x_k, y_k, 0]）
userXY = [
    2.0,  8.0;
    5.5, 11.5;
    8.0, 16.0
];

% X(:, n) = [y_{n,1}, ..., y_{n,M}]^T，表示第 n 条波导上 M 个 PA 的纵向坐标
X = [
    1.0,  1.5;
    3.0,  4.5;
    6.0,  8.0
];

% 每个 PA 的波束朝向（按 [M, N] 存储）
% theta \in (pi/2, pi]：俯仰角；phi \in (-pi, pi]：方位角
theta = deg2rad([
    135, 140;
    130, 135;
    125, 130
]);

phi = deg2rad([
     10, -15;
      5, -10;
      0,  -5
]);

%% =========================
%  2. 几何建模
%  =========================
[users, waveguideFeedPoints, paPositions] = buildGeometry(userXY, X, params);

disp('用户全局坐标 q_k = [x_k, y_k, 0]^T：');
disp(users);

disp('波导馈电点 [x_n^(W), 0, d]^T：');
disp(waveguideFeedPoints);

disp('第 1 条波导上所有 PA 的位置 p_{1,m}：');
disp(squeeze(paPositions(:, :, 1)).');

%% =========================
%  3. 单条波导导波信道 g_n(x_n)
%  =========================
g1 = guidedWaveChannel(X(:, 1), params);
disp('第 1 条波导的导波信道 g_1(x_1)：');
disp(g1);

%% =========================
%  4. 用户-PA 坐标变换与自由空间信道
%  =========================
userIdx = 1;
wgIdx = 1;

[hFreeExample, localCoordsExample] = freeSpaceChannelForWaveguide( ...
    users(userIdx, :).', squeeze(paPositions(:, :, wgIdx)), ...
    theta(:, wgIdx), phi(:, wgIdx), params);

disp('用户 1 相对于第 1 条波导上各 PA 的局部坐标 [x~, y~, z~]：');
disp(localCoordsExample.');

disp('用户 1 与第 1 条波导之间的自由空间信道 \\tilde{h}_{1,1}(x_1)：');
disp(hFreeExample);

%% =========================
%  5. 所有用户的总复合信道 h_k(X)
%  =========================
[H, details] = compositeChannel(users, X, theta, phi, params);

disp('所有用户的总复合信道矩阵 H（每列对应一个用户 h_k）：');
disp(H);

% details 结构体中保留了中间变量，便于后续优化、调试或画图：
% details.g{n}            -> 第 n 条波导的导波信道
% details.hFree{k, n}     -> 用户 k 与第 n 条波导的自由空间信道
% details.hComposite{k}   -> 用户 k 的总复合信道
% details.localCoords{k,n}-> 用户 k 在第 n 条波导各 PA 局部坐标系下的位置

%% =========================
%  6. 示例：单个用户、单个 PA 的局部坐标变换
%  =========================
k = 1; n = 1; m = 1;
qk = users(k, :).';
pnm = paPositions(:, m, n);
localCoord = globalToPALocal(qk, pnm, theta(m, n), phi(m, n));

fprintf('用户 %d 到波导 %d 的第 %d 个 PA 的局部坐标为:\n', k, n, m);
disp(localCoord.');

%% 局部函数定义

function params = defaultParameters()
%DEFAULTPARAMETERS 设置示例所需的系统参数
    params.N = 2;                 % 波导数量
    params.M = 3;                 % 每条波导的 PA 数量
    params.K = 3;                 % 用户数量
    params.Dx = 10;               % 波导在 x 方向覆盖宽度
    params.d = 4;                 % 波导部署高度

    params.lambda = 0.01;         % 载波波长（m）
    params.nEff = 1.6;            % 波导有效折射率 n_eff
    params.alphaW = 0.08;         % 波导衰减系数 alpha_W（Np/m）

    params.alphaL = 0.96;         % LoS 存在系数/米
    params.nMedium = 1.0;         % 自由空间/介质折射率 n
    params.a = 0.45;              % 模式尺寸参数 a
    params.b = 0.30;              % 模式尺寸参数 b
    params.v = 1.1;               % 经验修正系数 v
end

function [users, feedPoints, paPositions] = buildGeometry(userXY, X, params)
%BUILDGEOMETRY 构造用户、波导馈电点和 PA 三维位置
% 输入：
%   userXY : [K, 2]，每行为 [x_k, y_k]
%   X      : [M, N]，第 n 列是第 n 条波导上的 y_{n,m}
% 输出：
%   users      : [K, 3]
%   feedPoints : [N, 3]
%   paPositions: [3, M, N]

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

function g = guidedWaveChannel(xn, params)
%GUIDEDWAVECHANNEL 论文中的导波信道 g_n(x_n)
% g_m = sqrt(1/M) * exp(-alpha_W/2 * y_{n,m}) * exp(-j * 2pi/lambda_g * y_{n,m})

    lambdaG = params.lambda / params.nEff;
    xn = xn(:);

    amplitude = sqrt(1 / params.M) .* exp(-(params.alphaW / 2) .* xn);
    phase = exp(-1j * (2 * pi / lambdaG) .* xn);
    g = amplitude .* phase;
end

function [hFree, localCoords] = freeSpaceChannelForWaveguide(userPos, paPosWaveguide, thetaN, phiN, params)
%FREESPACECHANNELFORWAVEGUIDE 计算用户 k 与第 n 条波导之间的自由空间信道
% 输入：
%   userPos        : [3,1]
%   paPosWaveguide : [3,M]
%   thetaN, phiN   : [M,1]
% 输出：
%   hFree          : [M,1]
%   localCoords    : [3,M]

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
% 输出：
%   H : [N*M, K]，第 k 列对应用户 k 的 h_k(X)

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
%GLOBALTOPALOCAL 实现论文公式：
% [x~, y~, z~]^T = R_x(theta - pi/2) * R_z(pi/2 - phi) * (q_k - p_{n,m})

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
%RADIATIONPATTERN 计算 PA 局部坐标系下的高斯波束辐射模式 Upsilon(x,y,z)
% 若用户位于 PA 后向半空间（即 y~ <= 0），则认为该方向不辐射，返回 0。

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

    w1 = params.v * params.a * lambda;
    w2 = params.v * params.b * lambda;

    W1 = lambda * y / (pi * n * w1);
    W2 = lambda * y / (pi * n * w2);
    R1 = y;
    R2 = y;
    Theta1 = atan(lambda * y / (pi * n * w1));
    Theta2 = atan(lambda * y / (pi * n * w2));

    % 归一化常数 B：满足 \int |Upsilon|^2 dx dz = 1
    B = sqrt(2 / (pi * w1 * w2));

    amplitude = sqrt((w1 * w2) / (W1 * W2)) * B;
    quadraticEnvelope = exp(-((x ^ 2) / (W1 ^ 2) + (z ^ 2) / (W2 ^ 2)));
    phase = exp(-1j * k * n * ((x ^ 2) / (2 * R1) + (z ^ 2) / (2 * R2) + y) ...
        + 1j * (Theta1 + Theta2) / 2);

    value = amplitude * quadraticEnvelope * phase;
end
