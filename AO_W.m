function [state, info] = AO_W(state, params)
%AO_W 交替优化中的独立 W-block：仅负责通过 WMMSE 更新 W。

    Hs = state.channelMatrix(state.S, :);
    H = Hs';
    numStreams = numel(state.S);

    if isempty(state.W) || size(state.W, 2) ~= numStreams
        W = initializeW(H, params);
    else
        W = state.W;
    end

    sumRateHistory = zeros(params.IW + 1, 1);
    metrics = evaluateWBlockMetrics(state, params, W, state.S);
    sumRateHistory(1) = metrics.sumRate;
    mseHistory = zeros(params.IW, numStreams);
    converged = false;

    for iter = 1:params.IW
        u = updateReceivers(Hs, W, params.sigma2);
        e = computeMSE(Hs, W, u, params.sigma2);
        v = 1 ./ max(real(e), eps);
        Wnew = updateW(H, u, v, params);
        newMetrics = evaluateWBlockMetrics(state, params, Wnew, state.S);
        mseHistory(iter, :) = real(e).';
        sumRateHistory(iter + 1) = newMetrics.sumRate;
        W = Wnew;

        if sumRateHistory(iter + 1) - sumRateHistory(iter) <= params.wmmseTol
            converged = true;
            break;
        end
    end

    finalMetrics = evaluateWBlockMetrics(state, params, W, state.S);
    state.W = W;
    state.sinr = finalMetrics.sinr;
    state.rate = finalMetrics.rate;
    state.sumRate = finalMetrics.sumRate;

    info = struct();
    info.receivers = updateReceivers(Hs, state.W, params.sigma2);
    info.mse = computeMSE(Hs, state.W, info.receivers, params.sigma2);
    info.weights = 1 ./ max(real(info.mse), eps);
    info.sumRateHistory = sumRateHistory(1:find(sumRateHistory ~= 0, 1, 'last'));
    info.mseHistory = mseHistory;
    info.converged = converged;
    info.iterations = numel(info.sumRateHistory) - 1;
end

function W0 = initializeW(H, params)
    A = H * H';
    B = H;
    W0 = solvePowerConstrained(A, B, params);
    if real(trace(W0 * W0')) < eps
        W0 = (1 / sqrt(size(H, 2))) * eye(size(H, 1), size(H, 2));
    end
end

function u = updateReceivers(Hs, W, sigma2)
    numStreams = size(Hs, 1);
    u = zeros(numStreams, 1);
    for k = 1:numStreams
        coupling = Hs(k, :) * W;
        u(k) = coupling(k) / (sum(abs(coupling).^2) + sigma2);
    end
end

function e = computeMSE(Hs, W, u, sigma2)
    numStreams = size(Hs, 1);
    e = zeros(numStreams, 1);
    for k = 1:numStreams
        coupling = Hs(k, :) * W;
        totalPower = sum(abs(coupling).^2) + sigma2;
        e(k) = abs(u(k))^2 * totalPower - 2 * real(u(k) * coupling(k)) + 1;
    end
end

function W = updateW(H, u, v, params)
    numAnt = size(H, 1);
    numStreams = size(H, 2);
    A = zeros(numAnt, numAnt);
    B = zeros(numAnt, numStreams);
    for k = 1:numStreams
        hk = H(:, k);
        A = A + v(k) * abs(u(k))^2 * (hk * hk');
        B(:, k) = v(k) * conj(u(k)) * hk;
    end
    W = solvePowerConstrained(A, B, params);
end

function W = solvePowerConstrained(A, B, params)
    Wfree = (A + 1e-12 * eye(size(A, 1))) \ B;
    if real(trace(Wfree * Wfree')) <= params.Pmax
        W = Wfree;
        return;
    end
    muLow = 0;
    muHigh = 1;
    while real(trace(computeW(A, B, muHigh) * computeW(A, B, muHigh)')) > params.Pmax
        muHigh = 2 * muHigh;
    end
    for iter = 1:params.muBisectionMaxIter
        muMid = 0.5 * (muLow + muHigh);
        Wmid = computeW(A, B, muMid);
        powerMid = real(trace(Wmid * Wmid'));
        if abs(powerMid - params.Pmax) <= params.muBisectionTol
            W = Wmid;
            return;
        end
        if powerMid > params.Pmax
            muLow = muMid;
        else
            muHigh = muMid;
        end
    end
    W = computeW(A, B, muHigh);
end

function W = computeW(A, B, mu)
    W = (A + mu * eye(size(A, 1))) \ B;
end

function metrics = evaluateWBlockMetrics(state, params, W, S)
    Hs = state.channelMatrix(S, :);
    numStreams = numel(S);
    sinr = zeros(numStreams, 1);
    rate = zeros(numStreams, 1);

    for k = 1:numStreams
        coupling = Hs(k, :) * W;
        signalPower = abs(coupling(k))^2;
        interferencePower = sum(abs(coupling).^2) - signalPower;
        sinr(k) = signalPower / (interferencePower + params.sigma2);
        rate(k) = log2(1 + sinr(k));
    end

    metrics = struct('sinr', sinr, 'rate', rate, 'sumRate', sum(rate), ...
        'scheduledUsers', S(:).', 'W', W);
end
