function [state, info] = AO_S(state, params, t)
%AO_S User-set block update with full swapTrace.

    info = struct( ...
        'triggered', mod(t, params.TS) == 0, ...
        'acceptedSwaps', 0, ...
        'bestDelta', 0, ...
        'swapTrace', [], ... % inner trace
        'breakReason', 'notTriggered');

    if ~info.triggered
        return;
    end

    swapTrace = repmat(struct( ...
        'swapIter', [], ...
        'weakUserSet', [], ...
        'strongUserSet', [], ...
        'evaluatedPairs', [], ...
        'bestDelta', [], ...
        'bestDeltaCoarse', [], ...
        'bestDeltaFinal', [], ...
        'numCandidatesEvaluated', [], ...
        'numCandidatesRefined', [], ...
        'bestPair', [], ...
        'acceptedUsingRefine', false, ...
        'accepted', false, ...
        'sumRateBefore', [], ...
        'sumRateAfter', [], ...
        'breakReason', ''), params.maxSwapPerUpdate, 1);

    for swapIter = 1:params.maxSwapPerUpdate
        sumRateBefore = state.sumRate;
        [weakUserPositions, weakUsers] = selectWeakUsers(state, params);
        strongExternal = selectStrongExternalUsers(state, params);

        swapTrace(swapIter).swapIter = swapIter;
        swapTrace(swapIter).weakUserSet = weakUsers;
        swapTrace(swapIter).strongUserSet = strongExternal;
        swapTrace(swapIter).sumRateBefore = sumRateBefore;

        if isempty(weakUserPositions) || isempty(strongExternal)
            swapTrace(swapIter).bestDelta = -inf;
            swapTrace(swapIter).bestDeltaCoarse = -inf;
            swapTrace(swapIter).bestDeltaFinal = -inf;
            swapTrace(swapIter).numCandidatesEvaluated = 0;
            swapTrace(swapIter).numCandidatesRefined = 0;
            swapTrace(swapIter).acceptedUsingRefine = false;
            swapTrace(swapIter).accepted = false;
            swapTrace(swapIter).sumRateAfter = state.sumRate;
            swapTrace(swapIter).breakReason = 'noWeakOrStrongCandidates';
            info.breakReason = 'noWeakOrStrongCandidates';
            break;
        end

        [bestState, bestDelta, evaluatedPairs, bestPair, evalSummary] = evaluateRestrictedSwapNeighborhood(state, params, weakUserPositions, weakUsers, strongExternal);
        swapTrace(swapIter).evaluatedPairs = evaluatedPairs;
        swapTrace(swapIter).bestDelta = bestDelta;
        swapTrace(swapIter).bestDeltaCoarse = evalSummary.bestDeltaCoarse;
        swapTrace(swapIter).bestDeltaFinal = evalSummary.bestDeltaFinal;
        swapTrace(swapIter).numCandidatesEvaluated = evalSummary.numCandidatesEvaluated;
        swapTrace(swapIter).numCandidatesRefined = evalSummary.numCandidatesRefined;
        swapTrace(swapIter).bestPair = bestPair;
        swapTrace(swapIter).acceptedUsingRefine = false;

        if bestDelta >= params.epsilonS
            state = bestState;
            info.acceptedSwaps = info.acceptedSwaps + 1;
            swapTrace(swapIter).accepted = true;
            swapTrace(swapIter).acceptedUsingRefine = evalSummary.bestCandidateUsedRefine;
            swapTrace(swapIter).sumRateAfter = state.sumRate;
            swapTrace(swapIter).breakReason = '';
            info.breakReason = '';
        else
            swapTrace(swapIter).accepted = false;
            swapTrace(swapIter).sumRateAfter = state.sumRate;
            swapTrace(swapIter).breakReason = 'bestDeltaBelowThreshold';
            info.breakReason = 'bestDeltaBelowThreshold';
            break;
        end

        if swapIter == params.maxSwapPerUpdate
            info.breakReason = 'reachedMaxSwapPerUpdate';
        end
    end

    validMask = arrayfun(@(x) ~isempty(x.swapIter), swapTrace);
    info.swapTrace = swapTrace(validMask);
    if isempty(info.swapTrace)
        info.swapTrace = repmat(swapTrace(1), 0, 1);
    end
    if ~isempty(info.swapTrace)
        info.bestDelta = info.swapTrace(end).bestDelta;
    end
end

function [weakUserPositions, weakUsers] = selectWeakUsers(state, params)
    [~, orderAsc] = sort(state.rate, 'ascend');
    weakUserPositions = orderAsc(1:min(params.Lin, numel(state.S)));
    weakUsers = state.S(weakUserPositions);
end

function strongExternal = selectStrongExternalUsers(state, params)
    external = setdiff(state.candidatePool, state.S, 'stable');
    if isempty(external)
        strongExternal = [];
        return;
    end

    [~, extOrder] = sort(state.Emax(external), 'descend');
    strongExternal = external(extOrder(1:min(params.Lout, numel(external))));
end

function [bestState, bestDelta, evaluatedPairs, bestPair, evalSummary] = evaluateRestrictedSwapNeighborhood(state, params, weakUserPositions, weakUsers, strongExternal)
    bestDelta = -inf;
    bestDeltaCoarse = -inf;
    bestState = state;
    bestPair = struct('weakUser', [], 'strongUser', [], 'position', []);
    bestCandidateUsedRefine = false;
    pairCount = numel(weakUserPositions) * numel(strongExternal);
    evaluatedPairs = repmat(struct( ...
        'weakUser', [], ...
        'strongUser', [], ...
        'candidateSet', [], ...
        'delta', [], ...
        'deltaCoarse', [], ...
        'deltaFinal', [], ...
        'coarseSumRate', [], ...
        'finalSumRate', [], ...
        'refinedUsed', false, ...
        'candidateWPower', [], ...
        'shortRefineIterations', 0, ...
        'numericalGuardTriggered', false), pairCount, 1);
    pairIndex = 0;
    shortMaxIter = 3;
    numRefineCandidates = min(2, pairCount);
    pairData = repmat(struct( ...
        'weakUser', [], ...
        'strongUser', [], ...
        'position', [], ...
        'candidateSet', [], ...
        'Wcoarse', [], ...
        'coarseMetrics', [], ...
        'deltaCoarse', -inf, ...
        'candidateWPower', NaN, ...
        'refinedUsed', false, ...
        'Wfinal', [], ...
        'finalMetrics', [], ...
        'deltaFinal', -inf, ...
        'shortRefineIterations', 0, ...
        'numericalGuardTriggered', false), pairCount, 1);

    for weakIdx = 1:numel(weakUserPositions)
        pos = weakUserPositions(weakIdx);
        weakUser = weakUsers(weakIdx);
        for strongIdx = 1:numel(strongExternal)
            strongUser = strongExternal(strongIdx);
            candidateUsers = state.S;
            candidateUsers(pos) = strongUser;

            Wcandidate = buildCandidateW(state, params, candidateUsers);
            [coarseGuardTriggered, candidateWPower] = isInvalidCandidateW(Wcandidate, params.Pmax);
            if coarseGuardTriggered
                coarseMetrics = struct('sinr', state.sinr, 'rate', state.rate, 'sumRate', -inf);
                deltaCoarse = -inf;
            else
                coarseMetrics = Signal_model('evaluate', state, params, Wcandidate, candidateUsers);
                if ~isfinite(coarseMetrics.sumRate)
                    coarseGuardTriggered = true;
                    deltaCoarse = -inf;
                else
                    deltaCoarse = coarseMetrics.sumRate - state.sumRate;
                end
            end

            pairIndex = pairIndex + 1;
            pairData(pairIndex).weakUser = weakUser;
            pairData(pairIndex).strongUser = strongUser;
            pairData(pairIndex).position = pos;
            pairData(pairIndex).candidateSet = candidateUsers;
            pairData(pairIndex).Wcoarse = Wcandidate;
            pairData(pairIndex).coarseMetrics = coarseMetrics;
            pairData(pairIndex).deltaCoarse = deltaCoarse;
            pairData(pairIndex).candidateWPower = candidateWPower;
            pairData(pairIndex).Wfinal = Wcandidate;
            pairData(pairIndex).finalMetrics = coarseMetrics;
            pairData(pairIndex).deltaFinal = deltaCoarse;
            pairData(pairIndex).numericalGuardTriggered = coarseGuardTriggered;
            if deltaCoarse > bestDeltaCoarse
                bestDeltaCoarse = deltaCoarse;
            end
        end
    end

    coarseDeltas = [pairData.deltaCoarse];
    [~, order] = sort(coarseDeltas, 'descend');
    refineIndices = order(1:numRefineCandidates);

    for ridx = 1:numel(refineIndices)
        idx = refineIndices(ridx);
        if ~isfinite(pairData(idx).deltaCoarse)
            continue;
        end
        [Wrefined, refinedMetrics, refineInfo] = refineCandidateWShort(state, params, pairData(idx).candidateSet, pairData(idx).Wcoarse, shortMaxIter);
        pairData(idx).refinedUsed = true;
        pairData(idx).shortRefineIterations = refineInfo.shortIterations;
        guardTriggered = pairData(idx).numericalGuardTriggered || refineInfo.numericalGuardTriggered;
        [refinedInvalid, refinedPower] = isInvalidCandidateW(Wrefined, params.Pmax);
        pairData(idx).candidateWPower = refinedPower;
        guardTriggered = guardTriggered || refinedInvalid;
        if guardTriggered || ~isfinite(refinedMetrics.sumRate)
            pairData(idx).numericalGuardTriggered = true;
            pairData(idx).finalMetrics = struct('sinr', state.sinr, 'rate', state.rate, 'sumRate', -inf);
            pairData(idx).deltaFinal = -inf;
            pairData(idx).Wfinal = Wrefined;
        else
            pairData(idx).numericalGuardTriggered = false;
            pairData(idx).finalMetrics = refinedMetrics;
            pairData(idx).deltaFinal = refinedMetrics.sumRate - state.sumRate;
            pairData(idx).Wfinal = Wrefined;
        end
    end

    for i = 1:pairCount
        evaluatedPairs(i).weakUser = pairData(i).weakUser;
        evaluatedPairs(i).strongUser = pairData(i).strongUser;
        evaluatedPairs(i).candidateSet = pairData(i).candidateSet;
        evaluatedPairs(i).delta = pairData(i).deltaFinal;
        evaluatedPairs(i).deltaCoarse = pairData(i).deltaCoarse;
        evaluatedPairs(i).deltaFinal = pairData(i).deltaFinal;
        evaluatedPairs(i).coarseSumRate = pairData(i).coarseMetrics.sumRate;
        evaluatedPairs(i).finalSumRate = pairData(i).finalMetrics.sumRate;
        evaluatedPairs(i).refinedUsed = pairData(i).refinedUsed;
        evaluatedPairs(i).candidateWPower = pairData(i).candidateWPower;
        evaluatedPairs(i).shortRefineIterations = pairData(i).shortRefineIterations;
        evaluatedPairs(i).numericalGuardTriggered = pairData(i).numericalGuardTriggered;

        if pairData(i).deltaFinal > bestDelta
            bestDelta = pairData(i).deltaFinal;
            bestState = state;
            bestState.S = pairData(i).candidateSet;
            bestState.sinr = pairData(i).finalMetrics.sinr;
            bestState.rate = pairData(i).finalMetrics.rate;
            bestState.sumRate = pairData(i).finalMetrics.sumRate;
            bestState.W = pairData(i).Wfinal;
            bestPair = struct('weakUser', pairData(i).weakUser, 'strongUser', pairData(i).strongUser, 'position', pairData(i).position);
            bestCandidateUsedRefine = pairData(i).refinedUsed;
        end
    end

    evalSummary = struct();
    evalSummary.numCandidatesEvaluated = pairCount;
    evalSummary.numCandidatesRefined = sum([pairData.refinedUsed]);
    evalSummary.bestDeltaCoarse = bestDeltaCoarse;
    evalSummary.bestDeltaFinal = bestDelta;
    evalSummary.bestCandidateUsedRefine = bestCandidateUsedRefine;
end

function Wcandidate = buildCandidateW(state, params, candidateUsers)
    HsCandidate = state.channelMatrix(candidateUsers, :);
    numStreams = size(HsCandidate, 1);
    numAnt = size(HsCandidate, 2);
    Wcandidate = zeros(numAnt, numStreams);
    tinyNorm = 1e-12;
    for k = 1:numStreams
        hk = HsCandidate(k, :).';
        hkNorm = norm(hk);
        if hkNorm > tinyNorm
            Wcandidate(:, k) = conj(hk) / hkNorm;
        else
            Wcandidate(:, k) = zeros(numAnt, 1);
        end
    end
    totalPower = real(trace(Wcandidate * Wcandidate'));
    if isfinite(totalPower) && totalPower > 0
        Wcandidate = sqrt(params.Pmax / totalPower) * Wcandidate;
    end
end

function [Wrefined, refinedMetrics, refineInfo] = refineCandidateWShort(state, params, candidateUsers, Winit, maxIterShort)
    Hs = state.channelMatrix(candidateUsers, :);
    H = Hs';
    Wrefined = Winit;
    refinedMetrics = Signal_model('evaluate', state, params, Wrefined, candidateUsers);
    numericalGuardTriggered = false;
    performedIter = 0;
    for iter = 1:maxIterShort
        [invalidW, ~] = isInvalidCandidateW(Wrefined, params.Pmax);
        if invalidW
            numericalGuardTriggered = true;
            break;
        end
        u = wmmseUpdateReceivers(Hs, Wrefined, params.sigma2);
        e = wmmseComputeMSE(Hs, Wrefined, u, params.sigma2);
        if any(~isfinite(e))
            numericalGuardTriggered = true;
            break;
        end
        v = 1 ./ max(real(e), eps);
        Wnext = wmmseUpdateW(H, u, v, params);
        [invalidNext, ~] = isInvalidCandidateW(Wnext, params.Pmax);
        if invalidNext
            numericalGuardTriggered = true;
            break;
        end
        metricsNext = Signal_model('evaluate', state, params, Wnext, candidateUsers);
        if ~isfinite(metricsNext.sumRate)
            numericalGuardTriggered = true;
            break;
        end
        Wrefined = Wnext;
        refinedMetrics = metricsNext;
        performedIter = iter;
    end
    refineInfo = struct('shortIterations', performedIter, 'numericalGuardTriggered', numericalGuardTriggered);
end

function u = wmmseUpdateReceivers(Hs, W, sigma2)
    numStreams = size(Hs, 1);
    u = zeros(numStreams, 1);
    for k = 1:numStreams
        hkH = Hs(k, :);
        coupling = hkH * W;
        desiredTerm = coupling(k);
        u(k) = desiredTerm / (sum(abs(coupling).^2) + sigma2);
    end
end

function e = wmmseComputeMSE(Hs, W, u, sigma2)
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

function W = wmmseUpdateW(H, u, v, params)
    numAnt = size(H, 1);
    numStreams = size(H, 2);
    A = zeros(numAnt, numAnt);
    B = zeros(numAnt, numStreams);
    for k = 1:numStreams
        hk = H(:, k);
        A = A + v(k) * abs(u(k))^2 * (hk * hk');
        B(:, k) = v(k) * conj(u(k)) * hk;
    end
    W = wmmseSolvePowerConstrained(A, B, params);
end

function W = wmmseSolvePowerConstrained(A, B, params)
    Wfree = wmmseComputeW(A, B, 0);
    if real(trace(Wfree * Wfree')) <= params.Pmax
        W = Wfree;
        return;
    end
    muLow = 0;
    muHigh = 1;
    while real(trace(wmmseComputeW(A, B, muHigh) * wmmseComputeW(A, B, muHigh)')) > params.Pmax
        muHigh = 2 * muHigh;
    end
    for iter = 1:params.muBisectionMaxIter
        muMid = 0.5 * (muLow + muHigh);
        Wmid = wmmseComputeW(A, B, muMid);
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
    W = wmmseComputeW(A, B, muHigh);
end

function W = wmmseComputeW(A, B, mu)
    W = (A + (mu + 1e-12) * eye(size(A, 1))) \ B;
end

function [invalid, powerVal] = isInvalidCandidateW(W, pmax)
    powerVal = real(trace(W * W'));
    invalid = any(~isfinite(W(:))) || ~isfinite(powerVal) || powerVal > pmax * (1 + 1e-6) || norm(W, 'fro') <= 1e-12;
end
