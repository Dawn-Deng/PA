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
        'delta', []), pairCount, 1);
    pairIndex = 0;

    for weakIdx = 1:numel(weakUserPositions)
        pos = weakUserPositions(weakIdx);
        weakUser = weakUsers(weakIdx);
        for strongIdx = 1:numel(strongExternal)
            strongUser = strongExternal(strongIdx);
            candidateUsers = state.S;
            candidateUsers(pos) = strongUser;

            metrics = Signal_model('evaluate', state, params, state.W, candidateUsers);
            delta = metrics.sumRate - state.sumRate;

            pairIndex = pairIndex + 1;
            evaluatedPairs(pairIndex).weakUser = weakUser;
            evaluatedPairs(pairIndex).strongUser = strongUser;
            evaluatedPairs(pairIndex).candidateSet = candidateUsers;
            evaluatedPairs(pairIndex).delta = delta;

            if delta > bestDelta
                bestDelta = delta;
                bestState = state;
                bestState.S = candidateUsers;
                bestState.sinr = metrics.sinr;
                bestState.rate = metrics.rate;
                bestState.sumRate = metrics.sumRate;
                bestPair = struct('weakUser', weakUser, 'strongUser', strongUser, 'position', pos);
            end
        end
    end
end
