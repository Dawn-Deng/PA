function [state, aoInfo] = AO_model(state, params)
%AO_MODEL 交替优化主过程 + 完整 iteration trace。

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

    iterationTrace = repmat(struct( ...
        'iteration', [], ...
        'sumRateBefore', [], ...
        'sumRateAfterW', [], ...
        'sumRateAfterAngle', [], ...
        'sumRateAfterX', [], ...
        'sumRateAfterS', [], ...
        'deltaW', [], 'deltaAngle', [], 'deltaX', [], 'deltaS', [], ...
        'snapshotS', [], 'snapshotX', [], 'snapshotTheta', [], 'snapshotPhi', ...
        'terminationCheck', []), params.Tmax, 1); % iteration trace

    state.sumRate = 0;
    state.sinr = zeros(numel(state.S), 1);
    state.rate = zeros(numel(state.S), 1);
    sumRateHistory(1) = state.sumRate;

    for t = 1:params.Tmax
        traceT = struct();
        traceT.iteration = t;
        traceT.sumRateBefore = state.sumRate;
        previousRate = state.sumRate;

        tic;
        [state, blockHistory(t).W] = AO_W(state, params);
        runtimeStats.iteration(t).WTime = toc;
        traceT.sumRateAfterW = state.sumRate;
        traceT.deltaW = traceT.sumRateAfterW - traceT.sumRateBefore;

        tic;
        [state, blockHistory(t).angle, angleMemory] = AO_angle(state, params, angleMemory);
        runtimeStats.iteration(t).angleTime = toc;
        traceT.sumRateAfterAngle = state.sumRate;
        traceT.deltaAngle = traceT.sumRateAfterAngle - traceT.sumRateAfterW;

        tic;
        [state, blockHistory(t).position, positionMemory] = AO_X(state, params, positionMemory);
        runtimeStats.iteration(t).XTime = toc;
        traceT.sumRateAfterX = state.sumRate;
        traceT.deltaX = traceT.sumRateAfterX - traceT.sumRateAfterAngle;

        tic;
        [state, blockHistory(t).userSet] = AO_S(state, params, t);
        runtimeStats.iteration(t).STime = toc;
        traceT.sumRateAfterS = state.sumRate;
        traceT.deltaS = traceT.sumRateAfterS - traceT.sumRateAfterX;

        runtimeStats.iteration(t).WInnerIterations = blockHistory(t).W.iterations;
        runtimeStats.iteration(t).angleAcceptedCount = blockHistory(t).angle.acceptedCount;
        runtimeStats.iteration(t).angleReanchorCount = blockHistory(t).angle.reanchorCount;
        runtimeStats.iteration(t).XAcceptedCount = blockHistory(t).position.acceptedCount;
        runtimeStats.iteration(t).XLineSearchSteps = sum([blockHistory(t).position.waveguideInfo.lineSearchSteps]);
        runtimeStats.iteration(t).SAcceptedSwaps = blockHistory(t).userSet.acceptedSwaps;
        runtimeStats.iteration(t).SEvaluatedPairs = countEvaluatedPairsFromSwapTrace(blockHistory(t).userSet.swapTrace);

        sumRateHistory(t + 1) = state.sumRate;
        state.iteration = t;
        deltaSumRate = sumRateHistory(t + 1) - previousRate;

        traceT.snapshotS = state.S;
        traceT.snapshotX = state.X;
        traceT.snapshotTheta = state.theta;
        traceT.snapshotPhi = state.phi;
        traceT.terminationCheck = struct( ...
            'deltaOuter', deltaSumRate, ...
            'epsilonOuter', params.epsilonOuter, ...
            'stop', abs(deltaSumRate) < params.epsilonOuter, ...
            'reason', ternary(abs(deltaSumRate) < params.epsilonOuter, 'epsilonOuter', 'continue'));
        iterationTrace(t) = traceT;

        if abs(deltaSumRate) < params.epsilonOuter
            runtimeStats = finalizeRuntimeStats(runtimeStats.iteration(1:t));
            aoInfo = finalizeAO(sumRateHistory(1:t + 1), blockHistory(1:t), true, t, 'epsilonOuter', runtimeStats, iterationTrace(1:t));
            return;
        end
    end

    runtimeStats = finalizeRuntimeStats(runtimeStats.iteration);
    aoInfo = finalizeAO(sumRateHistory, blockHistory, false, params.Tmax, 'Tmax', runtimeStats, iterationTrace);
end

function out = ternary(cond, a, b)
    if cond
        out = a;
    else
        out = b;
    end
end

function count = countEvaluatedPairsFromSwapTrace(swapTrace)
    count = 0;
    if isempty(swapTrace)
        return;
    end
    for i = 1:numel(swapTrace)
        if isfield(swapTrace(i), 'evaluatedPairs') && ~isempty(swapTrace(i).evaluatedPairs)
            count = count + numel(swapTrace(i).evaluatedPairs);
        end
    end
end

function runtimeStats = finalizeRuntimeStats(iterStats)
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

function aoInfo = finalizeAO(sumRateHistory, blockHistory, converged, iterations, terminationReason, runtimeStats, iterationTrace)
    aoInfo = struct();
    aoInfo.sumRateHistory = sumRateHistory;
    aoInfo.blockHistory = blockHistory;
    aoInfo.iterationTrace = iterationTrace; % iteration trace
    aoInfo.converged = converged;
    aoInfo.iterations = iterations;
    aoInfo.terminationReason = terminationReason;
    aoInfo.runtimeStats = runtimeStats;
end
