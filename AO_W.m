function [state, info] = AO_W(state, params)
%AO_W W-block update in the alternating-optimization loop via WMMSE.
% With (S^(t), X^(t), theta^(t), phi^(t)) fixed, this function solves
%
%   W^(t+1) = arg max_W R_sum(S^(t), X^(t), W, theta^(t), phi^(t))
%             s.t. tr(WW^H) <= Pmax
%
% by the standard WMMSE equivalence:
%   1) update MMSE receivers u_k
%   2) update WMMSE weights v_k = 1 / e_k
%   3) update W under the total power constraint via bisection on mu
%
% The inner loop terminates when either IW iterations are reached or the
% improvement of the sum rate is no larger than params.wmmseTol.

    scheduledUsers = state.S(:).';
    Hs = state.channelMatrix(scheduledUsers, :);
    H = Hs';
    numStreams = numel(scheduledUsers);

    if isempty(state.W) || size(state.W, 2) ~= numStreams
        % Initialization stage does not optimize W. When AO enters the
        % W-block for the first time, we only build a feasible non-optimized
        % starting point; the first effective W update is the updateW(...)
        % call inside the inner WMMSE loop below.
        W = initializeWPlaceholder(size(H, 1), numStreams, params.Pmax);
    else
        W = state.W;
    end

    initialMetrics = evaluateWBlockMetrics(Hs, W, params.sigma2, scheduledUsers);
    sumRateHistory = zeros(params.IW + 1, 1);
    sumRateHistory(1) = initialMetrics.sumRate;
    mseHistory = zeros(params.IW, numStreams);
    receiverHistory = zeros(numStreams, params.IW);
    weightHistory = zeros(numStreams, params.IW);
    converged = false;
    numIterations = params.IW;

    for iter = 1:params.IW
        u = updateReceivers(Hs, W, params.sigma2);
        e = computeMSE(Hs, W, u, params.sigma2);
        v = 1 ./ max(real(e), eps);

        W = updateW(H, u, v, params);
        currentMetrics = evaluateWBlockMetrics(Hs, W, params.sigma2, scheduledUsers);

        receiverHistory(:, iter) = u;
        weightHistory(:, iter) = v;
        mseHistory(iter, :) = real(e).';
        sumRateHistory(iter + 1) = currentMetrics.sumRate;

        if sumRateHistory(iter + 1) - sumRateHistory(iter) <= params.wmmseTol
            converged = true;
            numIterations = iter;
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
    info.mseHistory = mseHistory(1:numIterations, :);
    info.receiverHistory = receiverHistory(:, 1:numIterations);
    info.weightHistory = weightHistory(:, 1:numIterations);
    info.converged = converged;
    info.iterations = numIterations;
    info.scheduledUsers = scheduledUsers;
end

function W0 = initializeWPlaceholder(numAnt, numStreams, Pmax)
% Builds a feasible placeholder only. This is not interpreted as the first
% optimized W-update in the paper; the first effective update happens in
% the AO W-block through updateW(...) after u_k and v_k are formed.
    W0 = zeros(numAnt, numStreams);
    activeDiag = min(numAnt, numStreams);
    if activeDiag > 0
        W0(1:activeDiag, 1:activeDiag) = eye(activeDiag);
        currentPower = real(trace(W0 * W0'));
        if currentPower > 0
            W0 = sqrt(Pmax / currentPower) * W0;
        end
    end
end

function u = updateReceivers(Hs, W, sigma2)
% Implements the MMSE receiver update
%   u_k = (h_k^H w_k) / (sum_j |h_k^H w_j|^2 + sigma^2).
    numStreams = size(Hs, 1);
    u = zeros(numStreams, 1);
    for k = 1:numStreams
        hkH = Hs(k, :);          % hkH = h_k^H
        coupling = hkH * W;      % entries are h_k^H w_j
        desiredTerm = coupling(k); % desiredTerm = h_k^H w_k
        u(k) = desiredTerm / (sum(abs(coupling).^2) + sigma2);
    end
end

function e = computeMSE(Hs, W, u, sigma2)
% Implements
%   e_k = |u_k|^2 (sum_j |h_k^H w_j|^2 + sigma^2)
%         - 2 Re{u_k h_k^H w_k} + 1.
    numStreams = size(Hs, 1);
    e = zeros(numStreams, 1);
    for k = 1:numStreams
        hkH = Hs(k, :);          % hkH = h_k^H
        coupling = hkH * W;      % entries are h_k^H w_j
        desiredTerm = coupling(k); % desiredTerm = h_k^H w_k
        totalPower = sum(abs(coupling).^2) + sigma2;
        e(k) = abs(u(k))^2 * totalPower - 2 * real(u(k) * desiredTerm) + 1;
    end
end

function W = updateW(H, u, v, params)
% For fixed {u_k} and {v_k}, updates W by
%   w_k = (sum_j v_j |u_j|^2 h_j h_j^H + mu I)^(-1) v_k u_k^* h_k.
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
    Wfree = computeW(A, B, 0);
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
    W = (A + (mu + 1e-12) * eye(size(A, 1))) \ B;
end

function metrics = evaluateWBlockMetrics(Hs, W, sigma2, scheduledUsers)
    numStreams = numel(scheduledUsers);
    sinr = zeros(numStreams, 1);
    rate = zeros(numStreams, 1);

    for k = 1:numStreams
        hkH = Hs(k, :);          % hkH = h_k^H
        coupling = hkH * W;      % entries are h_k^H w_j
        desiredTerm = coupling(k); % desiredTerm = h_k^H w_k
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
