function [state, aoInfo] = AO_model(state, params)
%AO_MODEL 交替优化主过程，信号计算统一通过 Signal_model 接口完成。

    blockHistory = repmat(struct('W', [], 'angle', [], 'position', [], 'userSet', []), params.Tmax, 1);
    sumRateHistory = zeros(params.Tmax + 1, 1);
    positionMemory = [];
    angleMemory = [];

    runtimeStats = struct();
    runtimeStats.iteration = repmat(struct( ...
        'WTime', 0, 'angleTime', 0, 'XTime', 0, 'STime', 0, ...
        'WInnerIterations', 0, 'angleAcceptedCount', 0, 'angleReanchorCount', 0, ...
        'XAcceptedCount', 0, 'XLineSearchSteps', 0, ...
        'SAcceptedSwaps', 0, 'SEvaluatedPairs', 0), params.Tmax, 1);

    state.sumRate = 0;
    state.sinr = zeros(numel(state.S), 1);
    state.rate = zeros(numel(state.S), 1);
    sumRateHistory(1) = state.sumRate;

    for t = 1:params.Tmax
        previousRate = state.sumRate;

        tic;
        [state, blockHistory(t).W] = AO_W(state, params);
        runtimeStats.iteration(t).WTime = toc;

        tic;
        [state, blockHistory(t).angle, angleMemory] = AO_angle(state, params, angleMemory);
        runtimeStats.iteration(t).angleTime = toc;

        tic;
        [state, blockHistory(t).position, positionMemory] = AO_X(state, params, positionMemory);
        runtimeStats.iteration(t).XTime = toc;

        tic;
        [state, blockHistory(t).userSet] = AO_S(state, params, t);
        runtimeStats.iteration(t).STime = toc;

        runtimeStats.iteration(t).WInnerIterations = blockHistory(t).W.iterations;
        runtimeStats.iteration(t).angleAcceptedCount = blockHistory(t).angle.acceptedCount;
        runtimeStats.iteration(t).angleReanchorCount = blockHistory(t).angle.reanchorCount;
        runtimeStats.iteration(t).XAcceptedCount = blockHistory(t).position.acceptedCount;
        runtimeStats.iteration(t).XLineSearchSteps = sum([blockHistory(t).position.waveguideInfo.lineSearchSteps]);
        runtimeStats.iteration(t).SAcceptedSwaps = blockHistory(t).userSet.acceptedSwaps;
        runtimeStats.iteration(t).SEvaluatedPairs = countEvaluatedPairs(blockHistory(t).userSet.evaluatedPairs);

        sumRateHistory(t + 1) = state.sumRate;
        state.iteration = t;
        deltaSumRate = sumRateHistory(t + 1) - previousRate;
        % Outer stopping criterion after all four blocks:
        % stop when |R_sum^(t+1) - R_sum^(t)| < epsilonOuter, or when Tmax is reached.
        if abs(deltaSumRate) < params.epsilonOuter
            runtimeStats = finalizeRuntimeStats(runtimeStats.iteration(1:t));
            aoInfo = finalizeAO(sumRateHistory(1:t + 1), blockHistory(1:t), true, t, 'epsilonOuter', runtimeStats);
            return;
        end
    end

    runtimeStats = finalizeRuntimeStats(runtimeStats.iteration);
    aoInfo = finalizeAO(sumRateHistory, blockHistory, false, params.Tmax, 'Tmax', runtimeStats);
end

function count = countEvaluatedPairs(evaluatedPairs)
    if isempty(evaluatedPairs)
        count = 0;
        return;
    end
    count = sum(arrayfun(@(x) ~isempty(x.candidateSet), evaluatedPairs));
end

function runtimeStats = finalizeRuntimeStats(iterStats)
% 工程化增强：统一 runtime/profiling 汇总
    runtimeStats = struct();
    runtimeStats.iteration = iterStats;
    runtimeStats.totalWTime = sum([iterStats.WTime]);
    runtimeStats.totalAngleTime = sum([iterStats.angleTime]);
    runtimeStats.totalXTime = sum([iterStats.XTime]);
    runtimeStats.totalSTime = sum([iterStats.STime]);
    runtimeStats.totalTime = runtimeStats.totalWTime + runtimeStats.totalAngleTime + runtimeStats.totalXTime + runtimeStats.totalSTime;

    runtimeStats.totalWInnerIterations = sum([iterStats.WInnerIterations]);
    runtimeStats.totalAngleAccepted = sum([iterStats.angleAcceptedCount]);
    runtimeStats.totalReanchors = sum([iterStats.angleReanchorCount]);
    runtimeStats.totalXAccepted = sum([iterStats.XAcceptedCount]);
    runtimeStats.totalXLineSearchSteps = sum([iterStats.XLineSearchSteps]);
    runtimeStats.totalSAcceptedSwaps = sum([iterStats.SAcceptedSwaps]);
    runtimeStats.totalSEvaluatedPairs = sum([iterStats.SEvaluatedPairs]);
end

function aoInfo = finalizeAO(sumRateHistory, blockHistory, converged, iterations, terminationReason, runtimeStats)
    aoInfo = struct();
    aoInfo.sumRateHistory = sumRateHistory;
    aoInfo.blockHistory = blockHistory;
    aoInfo.converged = converged;
    aoInfo.iterations = iterations;
    aoInfo.terminationReason = terminationReason;
    aoInfo.runtimeStats = runtimeStats; % 工程化增强：统一统计输出
end
