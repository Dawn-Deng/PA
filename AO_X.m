function [state, info, memory] = AO_X(state, params, memory)
%AO_X Position-block update with detailed line-search trace.

    if nargin < 3 || isempty(memory)
        memory = initializePositionMemory(params);
    end

    acceptedCount = 0;
    waveguideInfo = repmat(struct( ...
        'waveguideIndex', [], ...
        'accepted', false, ...
        'rejectReason', '', ... % reject reason
        'lineSearchSteps', 0, ...
        'alphaAccepted', [], ...
        'sumRateBefore', [], ...
        'sumRateAfter', [], ...
        'gradCurrent', [], ...
        'direction', [], ...
        'bestCoarseIndex', [], ...
        'bestCoarseRate', [], ...
        'bestFinalRate', [], ...
        'refinedUsed', false, ...
        'shortRefineIterations', 0, ...
        'numericalGuardTriggered', false, ...
        'projectionCollapsedCount', 0, ...
        'lineSearchTrace', []), params.N, 1); % inner trace

    for n = 1:params.N
        xCurrent = state.X(:, n);
        fCurrent = state.sumRate;
        gradCurrent = numericalGradient(state, params, n, xCurrent);
        direction = -lbfgsDirection(gradCurrent, memory{n});
        if norm(direction) <= 1e-12
            direction = -gradCurrent;
        end

        alpha = params.positionLineSearchInit;
        accepted = false;
        lineSearchSteps = 0;
        candidateStateBest = state;
        acceptedAlpha = [];
        rejectReason = 'noImprovement';
        projectionCollapsedCount = 0;
        lsTrace = repmat(struct( ...
            'alphaTrial', [], ...
            'xCandidate', [], ...
            'candidateRate', [], ...
            'candidateRateCoarse', [], ...
            'projectionDistance', [], ...
            'isBestCoarseCandidate', false, ...
            'refinedUsed', false, ...
            'candidateRateFinal', [], ...
            'accepted', false, ...
            'rejectReason', '', ...
            'rejectReasonFinal', ''), 0, 1);

        while alpha >= params.positionLineSearchMin
            lineSearchSteps = lineSearchSteps + 1;
            xBar = xCurrent + alpha * direction;
            xCandidate = Channel_model('project_waveguide_positions', xBar, params);
            projectionDistance = norm(xCandidate - xCurrent);
            if projectionDistance <= 1e-12
                projectionCollapsedCount = projectionCollapsedCount + 1;
            end

            candidateState = state;
            candidateState.X(:, n) = xCandidate;
            candidateState = Channel_model('update_state', candidateState, params);
            candidateMetricsCoarse = Signal_model('evaluate', candidateState, params, state.W, state.S);
            coarseRate = candidateMetricsCoarse.sumRate;

            thisRejectReason = '';
            if ~isfinite(coarseRate)
                thisRejectReason = 'invalidCoarseMetrics';
            elseif projectionDistance <= 1e-12
                thisRejectReason = 'projectionCollapsedMove'; % reject reason
            elseif coarseRate < fCurrent + params.epsilonX
                thisRejectReason = 'belowTolerance'; % reject reason
            end

            lsTrace(end + 1) = struct( ...
                'alphaTrial', alpha, ...
                'xCandidate', xCandidate, ...
                'candidateRate', coarseRate, ...
                'candidateRateCoarse', coarseRate, ...
                'projectionDistance', projectionDistance, ...
                'isBestCoarseCandidate', false, ...
                'refinedUsed', false, ...
                'candidateRateFinal', [], ...
                'accepted', false, ...
                'rejectReason', thisRejectReason, ...
                'rejectReasonFinal', ''); %#ok<AGROW>

            alpha = alpha * params.positionLineSearchBeta;
        end

        coarseRates = [lsTrace.candidateRateCoarse];
        projectionDistances = [lsTrace.projectionDistance];
        validMask = isfinite(coarseRates) & (projectionDistances > 1e-12);
        bestCoarseIndex = [];
        bestCoarseRate = -inf;
        bestFinalRate = -inf;
        refinedUsed = false;
        shortRefineIterations = 0;
        refineGuardTriggered = false;

        if any(validMask)
            validIdx = find(validMask);
            [bestCoarseRate, relIdx] = max(coarseRates(validIdx));
            bestCoarseIndex = validIdx(relIdx);
            lsTrace(bestCoarseIndex).isBestCoarseCandidate = true;

            candidateStateRefine = state;
            candidateStateRefine.X(:, n) = lsTrace(bestCoarseIndex).xCandidate;
            candidateStateRefine = Channel_model('update_state', candidateStateRefine, params);
            maxIterShort = 3;
            [Wrefined, refinedMetrics, refineInfo] = refinePositionCandidateWShort( ...
                candidateStateRefine, params, state.W, maxIterShort);

            refinedUsed = true;
            shortRefineIterations = refineInfo.shortIterations;
            refineGuardTriggered = refineInfo.numericalGuardTriggered;
            bestFinalRate = refinedMetrics.sumRate;

            lsTrace(bestCoarseIndex).refinedUsed = true;
            lsTrace(bestCoarseIndex).candidateRateFinal = refinedMetrics.sumRate;

            if refineInfo.numericalGuardTriggered
                rejectReason = 'numericalGuardTriggered';
                lsTrace(bestCoarseIndex).rejectReasonFinal = 'numericalGuardTriggered';
            elseif refinedMetrics.sumRate >= fCurrent + params.epsilonX
                candidateStateRefine.W = Wrefined;
                candidateStateRefine.sinr = refinedMetrics.sinr;
                candidateStateRefine.rate = refinedMetrics.rate;
                candidateStateRefine.sumRate = refinedMetrics.sumRate;
                candidateStateBest = candidateStateRefine;
                accepted = true;
                acceptedAlpha = lsTrace(bestCoarseIndex).alphaTrial;
                rejectReason = '';
                lsTrace(bestCoarseIndex).accepted = true;
                lsTrace(bestCoarseIndex).rejectReasonFinal = '';
            else
                rejectReason = 'refinedBelowTolerance';
                lsTrace(bestCoarseIndex).rejectReasonFinal = 'refinedBelowTolerance';
            end
        else
            if ~isempty(lsTrace) && all(projectionDistances <= 1e-12)
                rejectReason = 'projectionCollapsedMove';
            elseif alpha < params.positionLineSearchMin
                rejectReason = 'stepTooSmall'; % reject reason
            else
                rejectReason = 'stepTooSmall';
            end
        end

        if accepted
            xNew = candidateStateBest.X(:, n);
            gradNew = numericalGradient(candidateStateBest, params, n, xNew);
            memory{n} = updateLbfgsMemory(memory{n}, xNew - xCurrent, gradNew - gradCurrent, params.positionMemory);
            state = candidateStateBest;
            acceptedCount = acceptedCount + 1;
        end

        waveguideInfo(n).waveguideIndex = n;
        waveguideInfo(n).accepted = accepted;
        waveguideInfo(n).rejectReason = rejectReason;
        waveguideInfo(n).lineSearchSteps = lineSearchSteps;
        waveguideInfo(n).alphaAccepted = acceptedAlpha;
        waveguideInfo(n).sumRateBefore = fCurrent;
        waveguideInfo(n).sumRateAfter = state.sumRate;
        waveguideInfo(n).gradCurrent = gradCurrent;
        waveguideInfo(n).direction = direction;
        waveguideInfo(n).bestCoarseIndex = bestCoarseIndex;
        waveguideInfo(n).bestCoarseRate = bestCoarseRate;
        waveguideInfo(n).bestFinalRate = bestFinalRate;
        waveguideInfo(n).refinedUsed = refinedUsed;
        waveguideInfo(n).shortRefineIterations = shortRefineIterations;
        waveguideInfo(n).numericalGuardTriggered = refineGuardTriggered;
        waveguideInfo(n).projectionCollapsedCount = projectionCollapsedCount;
        waveguideInfo(n).lineSearchTrace = lsTrace;
    end

    info = struct( ...
        'acceptedCount', acceptedCount, ...
        'sumRate', state.sumRate, ...
        'waveguideInfo', waveguideInfo);
end

function memory = initializePositionMemory(params)
    memory = cell(params.N, 1);
    for n = 1:params.N
        memory{n}.S = [];
        memory{n}.Y = [];
    end
end

function grad = numericalGradient(state, params, n, x)
    grad = zeros(params.M, 1);
    delta = params.positionFiniteDiff;

    for m = 1:params.M
        xp = x;
        xm = x;
        xp(m) = xp(m) + delta;
        xm(m) = xm(m) - delta;

        stateP = state;
        stateM = state;
        stateP.X(:, n) = Channel_model('project_waveguide_positions', xp, params);
        stateM.X(:, n) = Channel_model('project_waveguide_positions', xm, params);
        stateP = Channel_model('update_state', stateP, params);
        stateM = Channel_model('update_state', stateM, params);

        fp = Signal_model('evaluate', stateP, params, state.W, state.S).sumRate;
        fm = Signal_model('evaluate', stateM, params, state.W, state.S).sumRate;
        grad(m) = (fp - fm) / (2 * delta);
    end
end

function direction = lbfgsDirection(grad, memory)
    if isempty(memory.S)
        direction = grad;
        return;
    end

    S = memory.S;
    Y = memory.Y;
    historySize = size(S, 2);
    q = grad;
    alpha = zeros(historySize, 1);
    rho = zeros(historySize, 1);

    for idx = historySize:-1:1
        rho(idx) = 1 / max(Y(:, idx)' * S(:, idx), eps);
        alpha(idx) = rho(idx) * (S(:, idx)' * q);
        q = q - alpha(idx) * Y(:, idx);
    end

    gamma = (S(:, end)' * Y(:, end)) / max(Y(:, end)' * Y(:, end), eps);
    r = gamma * q;

    for idx = 1:historySize
        beta = rho(idx) * (Y(:, idx)' * r);
        r = r + S(:, idx) * (alpha(idx) - beta);
    end

    direction = r;
end

function memory = updateLbfgsMemory(memory, s, y, maxMemory)
    if s' * y <= 1e-12
        return;
    end

    if isempty(memory.S)
        memory.S = s;
        memory.Y = y;
        return;
    end

    memory.S = [memory.S, s];
    memory.Y = [memory.Y, y];

    if size(memory.S, 2) > maxMemory
        memory.S(:, 1) = [];
        memory.Y(:, 1) = [];
    end
end

function [Wrefined, refinedMetrics, refineInfo] = refinePositionCandidateWShort(candidateState, params, Winit, maxIterShort)
    scheduledUsers = candidateState.S(:).';
    Hs = candidateState.channelMatrix(scheduledUsers, :);
    H = Hs';
    numStreams = numel(scheduledUsers);

    useWinit = ~isempty(Winit) && size(Winit, 2) == numStreams && all(isfinite(Winit(:)));
    if useWinit
        W = Winit;
    else
        W = initializeWFromChannelShort(Hs, params.Pmax);
    end

    numericalGuardTriggered = false;
    shortIterations = 0;

    for iter = 1:maxIterShort
        u = updateReceiversShort(Hs, W, params.sigma2);
        e = computeMSEShort(Hs, W, u, params.sigma2);
        v = 1 ./ max(real(e), eps);
        [W, ~] = updateWShort(H, u, v, params);

        shortIterations = iter;
        currentPower = real(trace(W * W'));
        nearZeroW = norm(W, 'fro') <= 1e-12;
        if any(~isfinite(W(:))) || ~isfinite(currentPower) || currentPower > params.Pmax * (1 + 1e-6) || nearZeroW
            numericalGuardTriggered = true;
            break;
        end
    end

    if numericalGuardTriggered
        Wrefined = W;
        refinedMetrics = struct('sinr', [], 'rate', [], 'sumRate', -inf);
    else
        refinedMetrics = Signal_model('evaluate', candidateState, params, W, candidateState.S);
        if ~isfinite(refinedMetrics.sumRate)
            numericalGuardTriggered = true;
            refinedMetrics = struct('sinr', [], 'rate', [], 'sumRate', -inf);
        end
        Wrefined = W;
    end

    refineInfo = struct( ...
        'shortIterations', shortIterations, ...
        'numericalGuardTriggered', numericalGuardTriggered, ...
        'finalPower', real(trace(W * W')));
end

function W0 = initializeWFromChannelShort(Hs, Pmax)
    numStreams = size(Hs, 1);
    numAnt = size(Hs, 2);
    W0 = zeros(numAnt, numStreams);
    tinyNorm = 1e-12;
    for k = 1:numStreams
        hk = Hs(k, :).';
        hkNorm = norm(hk);
        if hkNorm > tinyNorm
            W0(:, k) = conj(hk) / hkNorm;
        end
    end
    currentPower = real(trace(W0 * W0'));
    if isfinite(currentPower) && currentPower > 0
        W0 = sqrt(Pmax / currentPower) * W0;
    end
end

function u = updateReceiversShort(Hs, W, sigma2)
    numStreams = size(Hs, 1);
    u = zeros(numStreams, 1);
    for k = 1:numStreams
        hkH = Hs(k, :);
        coupling = hkH * W;
        desiredTerm = coupling(k);
        u(k) = desiredTerm / (sum(abs(coupling).^2) + sigma2);
    end
end

function e = computeMSEShort(Hs, W, u, sigma2)
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

function [W, traceInfo] = updateWShort(H, u, v, params)
    numAnt = size(H, 1);
    numStreams = size(H, 2);
    A = zeros(numAnt, numAnt);
    B = zeros(numAnt, numStreams);
    for k = 1:numStreams
        hk = H(:, k);
        A = A + v(k) * abs(u(k))^2 * (hk * hk');
        B(:, k) = v(k) * conj(u(k)) * hk;
    end
    [W, muFinal, usedBisection] = solvePowerConstrainedShort(A, B, params);
    traceInfo = struct('muFinal', muFinal, 'usedBisection', usedBisection);
end

function [W, muFinal, usedBisection] = solvePowerConstrainedShort(A, B, params)
    Wfree = computeWShort(A, B, 0);
    if real(trace(Wfree * Wfree')) <= params.Pmax
        W = Wfree;
        muFinal = 0;
        usedBisection = false;
        return;
    end

    usedBisection = true;
    muLow = 0;
    muHigh = 1;
    while real(trace(computeWShort(A, B, muHigh) * computeWShort(A, B, muHigh)')) > params.Pmax
        muHigh = 2 * muHigh;
    end

    for iter = 1:params.muBisectionMaxIter
        muMid = 0.5 * (muLow + muHigh);
        Wmid = computeWShort(A, B, muMid);
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

    W = computeWShort(A, B, muHigh);
    muFinal = muHigh;
end

function W = computeWShort(A, B, mu)
    W = (A + (mu + 1e-12) * eye(size(A, 1))) \ B;
end
