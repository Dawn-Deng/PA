function varargout = Channel_model(varargin)
%CHANNEL_MODEL 五文件版本中的系统模型/几何/信道统一入口。
% 用法：
%   [params, state, info] = Channel_model();
%   [params, state, info] = Channel_model(paramsOverrideStruct);
%   state = Channel_model('update_state', state, params);
%   yRef  = Channel_model('reference_positions', params);
%   xProj = Channel_model('project_waveguide_positions', xRaw, params);

    if nargin == 0
        [params, state, info] = initializeSystem();
        varargout = {params, state, info};
        return;
    end

    if nargin == 1 && isstruct(varargin{1})
        [params, state, info] = initializeSystem(varargin{1});
        varargout = {params, state, info};
        return;
    end

    action = lower(string(varargin{1}));
    switch action
        case "initialize"
            override = struct();
            if nargin >= 2 && isstruct(varargin{2})
                override = varargin{2};
            end
            [params, state, info] = initializeSystem(override);
            varargout = {params, state, info};
        case "update_state"
            varargout = {updateState(varargin{2}, varargin{3})};
        case "reference_positions"
            varargout = {referencePositions(varargin{2})};
        case "project_waveguide_positions"
            varargout = {projectWaveguidePositions(varargin{2}, varargin{3})};
        otherwise
            error('Channel_model: 未知操作 %s', action);
    end
end

function [params, state, info] = initializeSystem(paramsOverride)
    if nargin < 1
        paramsOverride = struct();
    end

    params = defaultParameters();
    params = applyOverrides(params, paramsOverride);
    validateParameters(params); % 工程化增强：统一参数校验

    state = struct();
    state.users = generateUsers(params);
    state.X = repmat(referencePositions(params), 1, params.N);
    state.theta = pi * ones(params.M, params.N);
    state.phi = zeros(params.M, params.N);
    state.S = [];
    state.W = [];
    state.candidatePool = [];
    state.Emax = [];
    state.iteration = 0;
    state = updateState(state, params);

    info = struct();
    info.referencePositions = referencePositions(params);
    info.description = 'Baseline geometry before initialization.';
    info.userGenerationMode = params.userGeneration.mode;
end

function params = defaultParameters()
    params.N = 4;
    params.M = 8;
    params.K = 64;
    params.NRF = 4;
    params.Kmax = 4;
    params.KServ = min(params.NRF, params.Kmax);

    params.Dx = 10;
    params.Dy = 10;
    params.d = 4;
    params.deltaMin = 1.0;

    params.lambda = 0.01;
    params.nEff = 1.6;
    params.alphaW = 0.08;
    params.alphaL = 0.96;
    params.nMedium = 1.0;
    params.a = 0.45;
    params.b = 0.30;
    params.v = 1.1;

    params.Pmax = 1;
    params.sigma2 = 1e-4;
    params.lambdaMov = 1e-4;

    params.userRegionX = [0, params.Dx];
    params.userRegionY = [6, 20];
    params.randomSeed = 20260321;
    params.symbolSeed = 20260322;

    % 工程化增强：场景预设，默认保持 uniform 兼容旧流程
    params.userGeneration = struct();
    params.userGeneration.mode = 'uniform';
    params.userGeneration.hotspotCenters = [params.Dx*0.3, 9; params.Dx*0.7, 17];
    params.userGeneration.hotspotStd = [1.2, 1.5];

    params.IW = 40;
    params.wmmseTol = 1e-6;
    params.muBisectionTol = 1e-8;
    params.muBisectionMaxIter = 60;

    params.Tmax = 15;
    params.epsilonOuter = 1e-5;

    params.angleStepThetaInit = 0.20;
    params.angleStepPhiInit = 0.20;
    params.angleStepThetaMin = 1e-3;
    params.angleStepPhiMin = 1e-3;
    params.betaTheta = 0.5;
    params.betaPhi = 0.5;
    params.epsilonTheta = 1e-7;
    params.angleAnchorUsers = 4;
    params.angleReanchorFailureThreshold = 2;

    params.positionFiniteDiff = 1e-3;
    params.positionLineSearchInit = 1.0;
    params.positionLineSearchBeta = 0.5;
    params.positionLineSearchMin = 1e-4;
    params.positionMemory = 7;
    params.epsilonX = 1e-7;

    params.TS = 2;
    params.Lin = 2;
    params.Lout = 4;
    params.maxSwapPerUpdate = 4;
    params.epsilonS = 1e-7;

    % 工程化增强：日志与结果导出配置
    params.verbosity = 1;
    params.saveResults = false;
    params.savePath = 'results_latest.mat';
end

function params = applyOverrides(params, override)
    if isempty(override)
        return;
    end
    f = fieldnames(override);
    for i = 1:numel(f)
        params.(f{i}) = override.(f{i});
    end
end

function users = generateUsers(params)
% 工程化增强：支持 uniform / clustered / corridor 与外部用户坐标输入
    if isfield(params, 'userPositions') && ~isempty(params.userPositions)
        users = params.userPositions;
        return;
    end

    rng(params.randomSeed);
    modeName = lower(string(params.userGeneration.mode));
    switch modeName
        case "uniform"
            x = params.userRegionX(1) + diff(params.userRegionX) * rand(params.K, 1);
            y = params.userRegionY(1) + diff(params.userRegionY) * rand(params.K, 1);
        case "clustered"
            centers = params.userGeneration.hotspotCenters;
            stdVec = params.userGeneration.hotspotStd;
            numCenters = size(centers, 1);
            assign = randi(numCenters, params.K, 1);
            x = zeros(params.K, 1);
            y = zeros(params.K, 1);
            for k = 1:params.K
                c = assign(k);
                x(k) = centers(c, 1) + stdVec(1) * randn();
                y(k) = centers(c, 2) + stdVec(min(2, numel(stdVec))) * randn();
            end
            x = min(max(x, params.userRegionX(1)), params.userRegionX(2));
            y = min(max(y, params.userRegionY(1)), params.userRegionY(2));
        case "corridor"
            x = params.userRegionX(1) + diff(params.userRegionX) * rand(params.K, 1);
            yCenter = mean(params.userRegionY);
            y = yCenter + 0.15 * diff(params.userRegionY) * randn(params.K, 1);
            y = min(max(y, params.userRegionY(1)), params.userRegionY(2));
        otherwise
            error('Unknown userGeneration.mode: %s', modeName);
    end
    users = [x, y, zeros(params.K, 1)];
end

function yRef = referencePositions(params)
    if params.M == 1
        yRef = 0;
        return;
    end

    yRef = zeros(params.M, 1);
    for m = 1:params.M
        yRef(m) = (m - 1) * params.deltaMin + (m - 1) / (params.M - 1) * ...
            (params.Dy - (params.M - 1) * params.deltaMin);
    end
end

function state = updateState(state, params)
    [feedPoints, paPositions] = buildGeometry(state.X, params);
    [H, details] = buildCompositeChannel(state.users, state.X, state.theta, state.phi, paPositions, params);
    state.waveguideFeedPoints = feedPoints;
    state.paPositions = paPositions;
    state.H = H;
    state.channelMatrix = H';
    state.channelDetails = details;
    sanityCheckState(state, params); % 工程化增强：统一状态数值健康检查
end

function [feedPoints, paPositions] = buildGeometry(X, params)
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

function [H, details] = buildCompositeChannel(users, X, theta, phi, paPositions, params)
    K = size(users, 1);
    H = zeros(params.N * params.M, K);
    details = struct();
    details.g = cell(params.N, 1);
    details.hFree = cell(K, params.N);
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
    end
end

function g = guidedWaveChannel(xn, params)
    lambdaG = params.lambda / params.nEff;
    xn = xn(:);
    g = sqrt(1 / params.M) .* exp(-(params.alphaW / 2) .* xn) .* exp(-1j * (2 * pi / lambdaG) .* xn);
end

function [hFree, localCoords] = freeSpaceChannelForWaveguide(userPos, paPosWaveguide, thetaN, phiN, params)
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

function localCoord = globalToPALocal(userPos, paPos, theta, phi)
    delta = userPos - paPos;
    Rx = [1,0,0; 0,cos(theta - pi/2),-sin(theta - pi/2); 0,sin(theta - pi/2),cos(theta - pi/2)];
    Rz = [cos(pi/2 - phi),-sin(pi/2 - phi),0; sin(pi/2 - phi),cos(pi/2 - phi),0; 0,0,1];
    localCoord = Rx * Rz * delta;
end

function value = radiationPattern(localCoord, params)
    x = localCoord(1); y = localCoord(2); z = localCoord(3);
    if y <= 0
        value = 0;
        return;
    end
    lambda = params.lambda; n = params.nMedium; k = 2 * pi / lambda;
    v = params.v; a = params.a; b = params.b;
    W1 = lambda * y / (pi * n * (v * a * lambda));
    W2 = lambda * y / (pi * n * (v * b * lambda));
    Theta1 = atan(lambda * y / (pi * n * (v * a * lambda)));
    Theta2 = atan(lambda * y / (pi * n * (v * b * lambda)));
    B = sqrt(2 / (pi * (v * a * lambda) * (v * b * lambda)));
    amplitude = sqrt((v * a * lambda) * (v * b * lambda) / (W1 * W2)) * B;
    envelope = exp(-((x ^ 2) / (W1 ^ 2) + (z ^ 2) / (W2 ^ 2)));
    phase = exp(-1j * k * n * ((x ^ 2) / (2 * y) + (z ^ 2) / (2 * y) + y) + 1j * (Theta1 + Theta2) / 2);
    value = amplitude * envelope * phase;
end

function xProj = projectWaveguidePositions(xRaw, params)
    M = numel(xRaw);
    z = sort(xRaw(:) - (0:M-1)' * params.deltaMin);
    z = min(max(z, 0), params.Dy - (M - 1) * params.deltaMin);
    z = isotonicProjection(z);
    z = min(max(z, 0), params.Dy - (M - 1) * params.deltaMin);
    xProj = z + (0:M-1)' * params.deltaMin;
end

function zProj = isotonicProjection(z)
    n = numel(z);
    level = z(:).';
    weight = ones(1, n);
    idx = 1;
    while idx < numel(level)
        if level(idx) > level(idx + 1)
            newLevel = (weight(idx) * level(idx) + weight(idx + 1) * level(idx + 1)) / ...
                (weight(idx) + weight(idx + 1));
            level(idx) = newLevel;
            weight(idx) = weight(idx) + weight(idx + 1);
            level(idx + 1) = [];
            weight(idx + 1) = [];
            if idx > 1
                idx = idx - 1;
            end
        else
            idx = idx + 1;
        end
    end

    zProj = zeros(n, 1);
    startIdx = 1;
    for block = 1:numel(level)
        blockLength = weight(block);
        zProj(startIdx:startIdx + blockLength - 1) = level(block);
        startIdx = startIdx + blockLength;
    end
end
