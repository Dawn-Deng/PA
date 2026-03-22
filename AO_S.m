function [state, info] = AO_S(state, params, t)
%AO_S User-set block update in the alternating-optimization loop.
% A periodically triggered restricted single-swap best-improvement search is
% used. Between set-update iterations, the user set remains unchanged. When
% triggered, swap candidates are scored using the true sum rate with the
% current W fixed, without re-running a full WMMSE update for each swap.

    info = struct( ...
        'triggered', mod(t, params.TS) == 0, ...
        'acceptedSwaps', 0, ...
        'bestDelta', 0, ...
        'weakUserSet', [], ...
        'strongUserSet', [], ...
        'evaluatedPairs', []);

    if ~info.triggered
        return;
    end

    for swapIter = 1:params.maxSwapPerUpdate
        [weakUserPositions, weakUsers] = selectWeakUsers(state, params);
        strongExternal = selectStrongExternalUsers(state, params);
        info.weakUserSet = weakUsers;
        info.strongUserSet = strongExternal;

        if isempty(weakUserPositions) || isempty(strongExternal)
            info.bestDelta = -inf;
            break;
        end

        [bestState, bestDelta, evaluatedPairs] = evaluateRestrictedSwapNeighborhood(state, params, weakUserPositions, weakUsers, strongExternal);
        info.bestDelta = bestDelta;
        info.evaluatedPairs = evaluatedPairs;

        if bestDelta >= params.epsilonS
            state = bestState;
            info.acceptedSwaps = info.acceptedSwaps + 1;
        else
            break;
        end
    end
end

function [weakUserPositions, weakUsers] = selectWeakUsers(state, params)
% Internal weak-user set I_weak^(t): the Lin served users with the lowest
% individual rates under the current (X, W, theta, phi).
    [~, orderAsc] = sort(state.rate, 'ascend');
    weakUserPositions = orderAsc(1:min(params.Lin, numel(state.S)));
    weakUsers = state.S(weakUserPositions);
end

function strongExternal = selectStrongExternalUsers(state, params)
% External strong-user set J_strong^(t): the Lout users in C \ S sorted by
% descending initialization-stage Emax, used only for coarse screening.
    external = setdiff(state.candidatePool, state.S, 'stable');
    if isempty(external)
        strongExternal = [];
        return;
    end

    [~, extOrder] = sort(state.Emax(external), 'descend');
    strongExternal = external(extOrder(1:min(params.Lout, numel(external))));
end

function [bestState, bestDelta, evaluatedPairs] = evaluateRestrictedSwapNeighborhood(state, params, weakUserPositions, weakUsers, strongExternal)
% Enumerate all restricted single-swap pairs in
%   I_weak^(t) x J_strong^(t),
% and evaluate each candidate by the true sum rate with the current W held
% fixed. The best-improvement pair is selected.
    bestDelta = -inf;
    bestState = state;
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
            end
        end
    end
end
