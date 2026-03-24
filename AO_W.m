function [state, info] = AO_W(state, params)
%AO_W W-block update in the alternating-optimization loop via WMMSE.
% With (S^(t), X^(t), theta^(t), phi^(t)) fixed, this function solves
%
%   W^(t+1) = arg max_W R_sum(S^(t), X^(t), W, theta^(t), phi^(t))
%             s.t. tr(WW^H) <= Pmax
%
% by the standard WMMSE equivalence.

    scheduledUsers = state.S(:).';
    Hs = state.channelMatrix(scheduledUsers, :);
    H = Hs';
    numStreams = numel(scheduledUsers);

    if isempty(state.W) || size(state.W, 2) ~= numStreams
        % Initialization stage does not optimize W. First AO W-block only
        % creates feasible warm-start, then enters standard WMMSE updates.
        W = initializeWFromChannel(Hs, params.Pmax);
        warmStartReason = 'channelAwareWarmStart';
    else
        W = state.W;
        warmStartReason = 'reusePreviousW';
    end
    warmStartPower = real(trace(W * W'));
    warmStartColumnNorms = sqrt(sum(abs(W).^2, 1));

    initialMetrics = evaluateWBlockMetrics(Hs, W, params.sigma2, scheduledUsers);
    sumRateHistory = zeros(params.IW + 1, 1);
    sumRateHistory(1) = initialMetrics.sumRate;
    deltaSumRateHistory = zeros(params.IW, 1);       % inner trace
    transmitPowerHistory = zeros(params.IW, 1);      % inner trace
    muHistory = nan(params.IW, 1);                   % inner trace
    mseHistory = zeros(params.IW, numStreams);
    receiverHistory = zeros(numStreams, params.IW);
    weightHistory = zeros(numStreams, params.IW);
    converged = false;
    numIterations = params.IW;
    stopReason = 'IW';
    numericalGuardTriggered = false;

    for iter = 1:params.IW
        u = updateReceivers(Hs, W, params.sigma2);
        e = computeMSE(Hs, W, u, params.sigma2);
        v = 1 ./ max(real(e), eps);

        [W, wUpdateTrace] = updateW(H, u, v, params);
        currentMetrics = evaluateWBlockMetrics(Hs, W, params.sigma2, scheduledUsers);

        if ~isfinite(currentMetrics.sumRate)
            numericalGuardTriggered = true;
            currentMetrics = initialMetrics;
            stopReason = 'numericalGuardTriggered';
            numIterations = iter;
            break;
        end

        receiverHistory(:, iter) = u;
        weightHistory(:, iter) = v;
        mseHistory(iter, :) = real(e).';
        sumRateHistory(iter + 1) = currentMetrics.sumRate;
        deltaSumRateHistory(iter) = sumRateHistory(iter + 1) - sumRateHistory(iter);
        transmitPowerHistory(iter) = real(trace(W * W'));
        muHistory(iter) = wUpdateTrace.muFinal;

        if deltaSumRateHistory(iter) <= params.wmmseTol
            converged = true;
            numIterations = iter;
            stopReason = 'tol'; % reject/break reason
            break;
        end
    end

    finalMetrics = evaluateWBlockMetrics(Hs, W, params.sigma2, scheduledUsers);
    state.W = W;
    state.sinr = finalMetrics.sinr;
    state.rate = finalMetrics.rate;
    state.sumRate = finalMetrics.sumRate;

    info = struct();
    info.receivers = updateReceivers(Hs, state.W, params.sigma2);
    info.mse = computeMSE(Hs, state.W, info.receivers, params.sigma2);
    info.weights = 1 ./ max(real(info.mse), eps);
    info.sumRateHistory = sumRateHistory(1:numIterations + 1);
    info.deltaSumRateHistory = deltaSumRateHistory(1:numIterations); % inner trace
    info.transmitPowerHistory = transmitPowerHistory(1:numIterations); % inner trace
    info.muHistory = muHistory(1:numIterations); % inner trace
    info.mseHistory = mseHistory(1:numIterations, :);
    info.receiverHistory = receiverHistory(:, 1:numIterations);
    info.weightHistory = weightHistory(:, 1:numIterations);
    info.converged = converged;
    info.iterations = numIterations;
    info.scheduledUsers = scheduledUsers;
    info.stopReason = stopReason; % reject/break reason
    info.warmStartReason = warmStartReason;
    info.warmStartPower = warmStartPower;
    info.warmStartColumnNorms = warmStartColumnNorms;
    info.finalTransmitPower = real(trace(state.W * state.W'));
    info.numericalGuardTriggered = numericalGuardTriggered; % diagnostic info
end

function W0 = initializeWFromChannel(Hs, Pmax)
    numStreams = size(Hs, 1);
    numAnt = size(Hs, 2);
    W0 = zeros(numAnt, numStreams);
    tinyNorm = 1e-12;
    for k = 1:numStreams
        hk = Hs(k, :).';
        hkNorm = norm(hk);
        if hkNorm > tinyNorm
            W0(:, k) = conj(hk) / hkNorm;
        else
            W0(:, k) = zeros(numAnt, 1);
        end
    end
    currentPower = real(trace(W0 * W0'));
    if isfinite(currentPower) && currentPower > 0
        W0 = sqrt(Pmax / currentPower) * W0;
    end
end

function u = updateReceivers(Hs, W, sigma2)
    numStreams = size(Hs, 1);
    u = zeros(numStreams, 1);
    for k = 1:numStreams
        hkH = Hs(k, :);
        coupling = hkH * W;
        desiredTerm = coupling(k);
        u(k) = desiredTerm / (sum(abs(coupling).^2) + sigma2);
    end
end

function e = computeMSE(Hs, W, u, sigma2)
    numStreams = size(Hs, 1);
    e = zeros(numStreams, 1);
    for k = 1:numStreams
        hkH = Hs(k, :);
        coupling = hkH * W;
        desiredTerm = coupling(k);
        totalPower = sum(abs(coupling).^2) + sigma2;
        e(k) = abs(u(k))^2 * totalPower - 2 * real(u(k) * desiredTerm) + 1;
    end
end

function [W, traceInfo] = updateW(H, u, v, params)
    numAnt = size(H, 1);
    numStreams = size(H, 2);
    A = zeros(numAnt, numAnt);
    B = zeros(numAnt, numStreams);
    for k = 1:numStreams
        hk = H(:, k);
        A = A + v(k) * abs(u(k))^2 * (hk * hk');
        B(:, k) = v(k) * conj(u(k)) * hk;
    end
    [W, muFinal, usedBisection] = solvePowerConstrained(A, B, params);
    traceInfo = struct('muFinal', muFinal, 'usedBisection', usedBisection);
end

function [W, muFinal, usedBisection] = solvePowerConstrained(A, B, params)
    Wfree = computeW(A, B, 0);
    if real(trace(Wfree * Wfree')) <= params.Pmax
        W = Wfree;
        muFinal = 0;
        usedBisection = false;
        return;
    end

    usedBisection = true;
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
            muFinal = muMid;
            return;
        end
        if powerMid > params.Pmax
            muLow = muMid;
        else
            muHigh = muMid;
        end
    end

    W = computeW(A, B, muHigh);
    muFinal = muHigh;
end

function W = computeW(A, B, mu)
    W = (A + (mu + 1e-12) * eye(size(A, 1))) \ B;
end

function metrics = evaluateWBlockMetrics(Hs, W, sigma2, scheduledUsers)
    numStreams = numel(scheduledUsers);
    sinr = zeros(numStreams, 1);
    rate = zeros(numStreams, 1);

    for k = 1:numStreams
        hkH = Hs(k, :);
        coupling = hkH * W;
        desiredTerm = coupling(k);
        signalPower = abs(desiredTerm)^2;
        interferencePower = sum(abs(coupling).^2) - signalPower;
        sinr(k) = signalPower / (interferencePower + sigma2);
        rate(k) = log2(1 + sinr(k));
    end

    metrics = struct( ...
        'sinr', sinr, ...
        'rate', rate, ...
        'sumRate', sum(rate), ...
        'scheduledUsers', scheduledUsers, ...
        'W', W);
end
