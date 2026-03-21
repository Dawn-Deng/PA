function signalState = miso_signal_model(channelState, params)
%MISO_SIGNAL_MODEL 多用户 MISO 下行信号建模
% 输入：
%   channelState.H -> 复合信道矩阵 [NM,K]，第 k 列为 h_k(X)
%   params         -> 系统参数结构体
% 输出：
%   signalState    -> 发射、接收、SINR、速率等结果

    channelMatrix = channelState.H';
    scheduledUsers = selectScheduledUsers(params);
    symbols = generateNormalizedSymbols(numel(scheduledUsers), params.symbolSeed);
    [W, beamformerInfo] = designDigitalBeamformers(channelMatrix, scheduledUsers, params);
    slotResult = simulateTransmissionSlot(channelMatrix, scheduledUsers, W, symbols, params.sigma2);

    signalState = struct();
    signalState.channelMatrix = channelMatrix;
    signalState.scheduledUsers = scheduledUsers;
    signalState.symbols = symbols;
    signalState.beamformerInfo = beamformerInfo;
    signalState.W = W;
    signalState.slotResult = slotResult;
end

function scheduledUsers = selectScheduledUsers(params)
%SELECTSCHEDULEDUSERS 生成时隙 t 的调度用户集合 S^(t)
    candidatePool = params.candidateUserPool(:).';
    KServ = min([params.NRF, params.Kmax, numel(candidatePool)]);

    switch lower(params.schedulingMode)
        case 'initialization'
            if ~isempty(params.initialServiceSet)
                scheduledUsers = unique(params.initialServiceSet(:).', 'stable');
                if numel(scheduledUsers) < KServ
                    fallbackUsers = setdiff(candidatePool, scheduledUsers, 'stable');
                    scheduledUsers = [scheduledUsers, fallbackUsers(1:KServ - numel(scheduledUsers))];
                end
            else
                scheduledUsers = candidatePool(1:KServ);
            end
        case 'first'
            scheduledUsers = candidatePool(1:KServ);
        case 'random'
            rng(params.randomSeed + 2);
            perm = randperm(numel(candidatePool), KServ);
            scheduledUsers = sort(candidatePool(perm));
        otherwise
            error('未知的 schedulingMode: %s', params.schedulingMode);
    end
end

function symbols = generateNormalizedSymbols(numStreams, seed)
%GENERATENORMALIZEDSYMBOLS 生成满足 E{|s_k|^2}=1 的复高斯归一化数据符号
    rng(seed);
    symbols = (randn(numStreams, 1) + 1j * randn(numStreams, 1)) / sqrt(2);
end

function [W, info] = designDigitalBeamformers(channelMatrix, scheduledUsers, params)
%DESIGNDIGITALBEAMFORMERS 为调度用户设计数字波束形成矩阵
    Hs = channelMatrix(scheduledUsers, :);   % [Ns, NM]，每行是 h_k^H
    effectiveH = Hs';                        % [NM, Ns]，每列是 h_k
    numStreams = numel(scheduledUsers);

    switch upper(params.beamformingMethod)
        case 'MRT'
            W = effectiveH;
            info.method = 'MRT';
        case 'ZF'
            gramMatrix = Hs * Hs';
            W = Hs' * pinv(gramMatrix);
            info.method = 'ZF';
        otherwise
            error('未知的 beamformingMethod: %s', params.beamformingMethod);
    end

    for i = 1:numStreams
        columnNorm = norm(W(:, i));
        if columnNorm > 0
            W(:, i) = W(:, i) / columnNorm;
        end
    end

    W = sqrt(params.totalTransmitPower / numStreams) * W;
    info.totalPower = norm(W, 'fro')^2;
end

function result = simulateTransmissionSlot(channelMatrix, scheduledUsers, W, symbols, sigma2)
%SIMULATETRANSMISSIONSLOT 依据论文公式计算 x_rad、y、SINR 和速率
    numStreams = numel(scheduledUsers);
    numAntennas = size(channelMatrix, 2);

    if size(W, 2) ~= numStreams
        error('W 的列数必须等于调度用户数。');
    end
    if numel(symbols) ~= numStreams
        error('symbols 的长度必须等于调度用户数。');
    end

    symbols = symbols(:);
    xRad = reshape(W * symbols, [numAntennas, 1]);

    y = zeros(numStreams, 1);
    sinr = zeros(numStreams, 1);
    rate = zeros(numStreams, 1);
    noise = sqrt(sigma2 / 2) * (randn(numStreams, 1) + 1j * randn(numStreams, 1));
    userMetrics(numStreams, 1) = struct( ...
        'userIndex', [], 'desiredSignal', [], 'interference', [], ...
        'noise', [], 'receivedSignal', [], 'sinr', [], 'rate', []);

    for i = 1:numStreams
        userIdx = scheduledUsers(i);
        hkH = channelMatrix(userIdx, :);
        coupling = hkH * W;
        desiredSignal = coupling(i) * symbols(i);
        interference = coupling * symbols - desiredSignal;

        y(i) = desiredSignal + interference + noise(i);
        signalPower = abs(coupling(i))^2;
        interferencePower = sum(abs(coupling).^2) - signalPower;
        sinr(i) = signalPower / (interferencePower + sigma2);
        rate(i) = log2(1 + sinr(i));

        userMetrics(i).userIndex = userIdx;
        userMetrics(i).desiredSignal = desiredSignal;
        userMetrics(i).interference = interference;
        userMetrics(i).noise = noise(i);
        userMetrics(i).receivedSignal = y(i);
        userMetrics(i).sinr = sinr(i);
        userMetrics(i).rate = rate(i);
    end

    result = struct();
    result.xRad = xRad;
    result.y = y;
    result.sinr = sinr;
    result.rate = rate;
    result.noise = noise;
    result.scheduledUsers = scheduledUsers;
    result.symbols = symbols;
    result.beamformers = W;
    result.userMetrics = userMetrics;
end
