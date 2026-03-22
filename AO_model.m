function [state, aoInfo] = AO_model(state, params)
%AO_MODEL 交替优化主过程，信号计算统一通过 Signal_model 接口完成。

    blockHistory = repmat(struct('W', [], 'angle', [], 'position', [], 'userSet', []), params.Tmax, 1);
    sumRateHistory = zeros(params.Tmax + 1, 1);
    positionMemory = [];
    angleMemory = [];

    state.sumRate = 0;
    state.sinr = zeros(numel(state.S), 1);
    state.rate = zeros(numel(state.S), 1);
    sumRateHistory(1) = state.sumRate;

    for t = 1:params.Tmax
        previousRate = state.sumRate;
        [state, blockHistory(t).W] = AO_W(state, params);
        [state, blockHistory(t).angle, angleMemory] = AO_angle(state, params, angleMemory);
        [state, blockHistory(t).position, positionMemory] = AO_X(state, params, positionMemory);
        [state, blockHistory(t).userSet] = updateUserSet(state, params, t);
        sumRateHistory(t + 1) = state.sumRate;
        state.iteration = t;
        if abs(sumRateHistory(t + 1) - previousRate) < params.epsilonOuter
            aoInfo = finalizeAO(sumRateHistory(1:t + 1), blockHistory(1:t), true, t);
            return;
        end
    end

    aoInfo = finalizeAO(sumRateHistory, blockHistory, false, params.Tmax);
end

function aoInfo = finalizeAO(sumRateHistory, blockHistory, converged, iterations)
    aoInfo = struct();
    aoInfo.sumRateHistory = sumRateHistory;
    aoInfo.blockHistory = blockHistory;
    aoInfo.converged = converged;
    aoInfo.iterations = iterations;
end

function [state, info] = updateUserSet(state, params, t)
    info = struct('triggered', mod(t, params.TS) == 0, 'acceptedSwaps', 0, 'bestDelta', 0);
    if ~info.triggered
        return;
    end

    for swapIter = 1:params.maxSwapPerUpdate
        [~, orderAsc] = sort(state.rate, 'ascend');
        weakPositions = orderAsc(1:min(params.Lin, numel(state.S)));
        external = setdiff(state.candidatePool, state.S, 'stable');
        if isempty(external)
            break;
        end
        [~, extOrder] = sort(state.Emax(external), 'descend');
        strongExternal = external(extOrder(1:min(params.Lout, numel(external))));
        bestDelta = -inf;
        bestState = state;

        for w = 1:numel(weakPositions)
            pos = weakPositions(w);
            for j = 1:numel(strongExternal)
                candidateUsers = state.S;
                candidateUsers(pos) = strongExternal(j);
                metrics = Signal_model('evaluate', state, params, state.W, candidateUsers);
                delta = metrics.sumRate - state.sumRate;
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

        info.bestDelta = bestDelta;
        if bestDelta >= params.epsilonS
            state = bestState;
            info.acceptedSwaps = info.acceptedSwaps + 1;
        else
            break;
        end
    end
end
