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
        [state, blockHistory(t).userSet] = AO_S(state, params, t);
        sumRateHistory(t + 1) = state.sumRate;
        state.iteration = t;
        deltaSumRate = sumRateHistory(t + 1) - previousRate;
        % Outer stopping criterion after all four blocks:
        % stop when |R_sum^(t+1) - R_sum^(t)| < epsilonOuter, or when Tmax is reached.
        if abs(deltaSumRate) < params.epsilonOuter
            aoInfo = finalizeAO(sumRateHistory(1:t + 1), blockHistory(1:t), true, t, 'epsilonOuter');
            return;
        end
    end

    aoInfo = finalizeAO(sumRateHistory, blockHistory, false, params.Tmax, 'Tmax');
end

function aoInfo = finalizeAO(sumRateHistory, blockHistory, converged, iterations, terminationReason)
    aoInfo = struct();
    aoInfo.sumRateHistory = sumRateHistory;
    aoInfo.blockHistory = blockHistory;
    aoInfo.converged = converged;
    aoInfo.iterations = iterations;
    aoInfo.terminationReason = terminationReason;
end
