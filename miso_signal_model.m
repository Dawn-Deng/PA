function signalState = miso_signal_model(H, params)
%MISO_SIGNAL_MODEL 根据信论文公式建模多用户发射、接收、SINR 与频谱效率
% 1) 发射信号：x_rad^(t) = sum_k w_k^(t) s_k^(t)
% 2) 接收信号：y_k^(t) = h_k^H w_k s_k + sum_{j!=k} h_k^H w_j s_j + n_k
% 3) SINR 与频谱效率：gamma_k^(t), R_k^(t) = log2(1 + gamma_k^(t))

    validateattributes(H, {'numeric'}, {'2d', 'nrows', params.N * params.M, 'ncols', params.K});
    validateattributes(params, {'struct'}, {'nonempty'});

    scheduledUsers = scheduleUsers(H, params);
    scheduledChannels = H(:, scheduledUsers);
    numStreams = numel(scheduledUsers);

    symbols = generateNormalizedSymbols(numStreams, params);
    W = buildDigitalBeamformer(scheduledChannels, params);
    xRad = W * symbols;
    noiseSamples = generateNoiseSamples(numStreams, params);
    userMetrics = evaluateScheduledUsers(H, scheduledUsers, W, symbols, noiseSamples, params);

    signalState = struct();
    signalState.slotIndex = params.slotIndex;
    signalState.scheduledUsers = scheduledUsers;
    signalState.symbols = symbols;
    signalState.W = W;
    signalState.xRad = xRad;
    signalState.noiseSamples = noiseSamples;
    signalState.userMetrics = userMetrics;
    signalState.userMetricsTable = buildUserMetricsTable(userMetrics);
end

function scheduledUsers = scheduleUsers(H, params)
%SCHEDULEUSERS 生成每个传输时隙的调度用户集合 S^(t)

    numScheduled = min(params.NRF, params.K);

    switch lower(params.schedulingMode)
        case 'strongest'
            channelStrength = vecnorm(H, 2, 1);
            [~, order] = sort(channelStrength, 'descend');
            scheduledUsers = sort(order(1:numScheduled));

        case 'first'
            scheduledUsers = 1:numScheduled;

        otherwise
            error('未知的 schedulingMode: %s', params.schedulingMode);
    end

    scheduledUsers = scheduledUsers(:);
end

function symbols = generateNormalizedSymbols(numStreams, params)
%GENERATENORMALIZEDSYMBOLS 生成满足 E[|s_k^(t)|^2] = 1 的归一化数据符号

    rng(params.randomSeed + params.slotIndex);
    rawSymbols = sign(randn(numStreams, 1)) + 1j * sign(randn(numStreams, 1));
    symbols = rawSymbols / sqrt(2);
end

function W = buildDigitalBeamformer(HScheduled, params)
%BUILDDIGITALBEAMFORMER 构造数字波束赋形向量 w_k^(t)

    [numAntennas, numStreams] = size(HScheduled);
    W = zeros(numAntennas, numStreams);

    switch lower(params.precoderMode)
        case 'mrt'
            powerPerStream = params.totalTransmitPower / max(numStreams, 1);

            for idx = 1:numStreams
                hk = HScheduled(:, idx);
                hkNorm = norm(hk);

                if hkNorm > 0
                    W(:, idx) = sqrt(powerPerStream) * hk / hkNorm;
                end
            end

        otherwise
            error('未知的 precoderMode: %s', params.precoderMode);
    end
end

function noiseSamples = generateNoiseSamples(numStreams, params)
%GENERATENOISESAMPLES 生成加性白高斯噪声样本 n_k^(t)

    rng(params.randomSeed + 10 * params.slotIndex);
    noiseSamples = sqrt(params.sigma2 / 2) * ...
        (randn(numStreams, 1) + 1j * randn(numStreams, 1));
end

function userMetrics = evaluateScheduledUsers(H, scheduledUsers, W, symbols, noiseSamples, params)
%EVALUATESCHEDULEDUSERS 计算接收信号、SINR 与频谱效率

    numStreams = numel(scheduledUsers);
    userMetrics = repmat(struct( ...
        'userIndex', 0, ...
        'desiredSignal', 0, ...
        'interference', 0, ...
        'noiseSample', 0, ...
        'receivedSignal', 0, ...
        'sinr', 0, ...
        'rate', 0), numStreams, 1);

    for idx = 1:numStreams
        userIndex = scheduledUsers(idx);
        hk = H(:, userIndex);
        desiredSignal = (hk' * W(:, idx)) * symbols(idx);

        interference = 0;
        interferencePower = 0;
        for jdx = 1:numStreams
            if jdx == idx
                continue;
            end

            interference = interference + (hk' * W(:, jdx)) * symbols(jdx);
            interferencePower = interferencePower + abs(hk' * W(:, jdx)) ^ 2;
        end

        receivedSignal = desiredSignal + interference + noiseSamples(idx);
        desiredPower = abs(hk' * W(:, idx)) ^ 2;
        sinr = desiredPower / (interferencePower + params.sigma2);
        rate = log2(1 + sinr);

        userMetrics(idx).userIndex = userIndex;
        userMetrics(idx).desiredSignal = desiredSignal;
        userMetrics(idx).interference = interference;
        userMetrics(idx).noiseSample = noiseSamples(idx);
        userMetrics(idx).receivedSignal = receivedSignal;
        userMetrics(idx).sinr = sinr;
        userMetrics(idx).rate = rate;
    end
end

function metricsTable = buildUserMetricsTable(userMetrics)
%BUILDUSERMETRICSTABLE 将结构体结果整理为便于查看的表格

    userIndex = transpose([userMetrics.userIndex]);
    desiredSignal = transpose([userMetrics.desiredSignal]);
    interference = transpose([userMetrics.interference]);
    noiseSample = transpose([userMetrics.noiseSample]);
    receivedSignal = transpose([userMetrics.receivedSignal]);
    sinr = transpose([userMetrics.sinr]);
    rate = transpose([userMetrics.rate]);

    metricsTable = table( ...
        userIndex, desiredSignal, interference, noiseSample, receivedSignal, sinr, rate, ...
        'VariableNames', {'userIndex', 'desiredSignal', 'interference', 'noiseSample', ...
        'receivedSignal', 'sinr', 'rate'});
end
