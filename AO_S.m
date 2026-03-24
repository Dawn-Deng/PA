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
        'bestPair', [], ...
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
            swapTrace(swapIter).accepted = false;
            swapTrace(swapIter).sumRateAfter = state.sumRate;
            swapTrace(swapIter).breakReason = 'noWeakOrStrongCandidates';
            info.breakReason = 'noWeakOrStrongCandidates';
            break;
        end

        [bestState, bestDelta, evaluatedPairs, bestPair] = evaluateRestrictedSwapNeighborhood(state, params, weakUserPositions, weakUsers, strongExternal);
        swapTrace(swapIter).evaluatedPairs = evaluatedPairs;
        swapTrace(swapIter).bestDelta = bestDelta;
        swapTrace(swapIter).bestPair = bestPair;

        if bestDelta >= params.epsilonS
            state = bestState;
            info.acceptedSwaps = info.acceptedSwaps + 1;
            swapTrace(swapIter).accepted = true;
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

function [bestState, bestDelta, evaluatedPairs, bestPair] = evaluateRestrictedSwapNeighborhood(state, params, weakUserPositions, weakUsers, strongExternal)
    bestDelta = -inf;
    bestState = state;
    bestPair = struct('weakUser', [], 'strongUser', [], 'position', []);
    pairCount = numel(weakUserPositions) * numel(strongExternal);
    evaluatedPairs = repmat(struct( ...
        'weakUser', [], ...
        'strongUser', [], ...
        'candidateSet', [], ...
        'delta', [], ...
        'candidateWPower', [], ...
        'numericalGuardTriggered', false), pairCount, 1);
    pairIndex = 0;

    for weakIdx = 1:numel(weakUserPositions)
        pos = weakUserPositions(weakIdx);
        weakUser = weakUsers(weakIdx);
        for strongIdx = 1:numel(strongExternal)
            strongUser = strongExternal(strongIdx);
            candidateUsers = state.S;
            candidateUsers(pos) = strongUser;

            Wcandidate = buildCandidateW(state, params, candidateUsers);
            candidateWPower = real(trace(Wcandidate * Wcandidate'));
            candidateWNorm = norm(Wcandidate, 'fro');
            numericalGuardTriggered = any(~isfinite(Wcandidate(:))) || ~isfinite(candidateWPower) || ...
                candidateWPower > params.Pmax * (1 + 1e-6) || candidateWNorm <= 1e-12;
            if numericalGuardTriggered
                delta = -inf;
                metrics = struct('sinr', state.sinr, 'rate', state.rate, 'sumRate', -inf);
            else
                metrics = Signal_model('evaluate', state, params, Wcandidate, candidateUsers);
                if ~isfinite(metrics.sumRate)
                    numericalGuardTriggered = true;
                    delta = -inf;
                else
                    delta = metrics.sumRate - state.sumRate;
                end
            end

            pairIndex = pairIndex + 1;
            evaluatedPairs(pairIndex).weakUser = weakUser;
            evaluatedPairs(pairIndex).strongUser = strongUser;
            evaluatedPairs(pairIndex).candidateSet = candidateUsers;
            evaluatedPairs(pairIndex).delta = delta;
            evaluatedPairs(pairIndex).candidateWPower = candidateWPower;
            evaluatedPairs(pairIndex).numericalGuardTriggered = numericalGuardTriggered;

            if delta > bestDelta
                bestDelta = delta;
                bestState = state;
                bestState.S = candidateUsers;
                bestState.sinr = metrics.sinr;
                bestState.rate = metrics.rate;
                bestState.sumRate = metrics.sumRate;
                bestState.W = Wcandidate;
                bestPair = struct('weakUser', weakUser, 'strongUser', strongUser, 'position', pos);
            end
        end
    end
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
