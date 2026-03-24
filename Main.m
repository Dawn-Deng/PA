%% Main.m
% 工程化 demo 入口：单次运行 + 全量过程追踪导出

clear; clc;

% 工程化增强：可在此处覆写默认参数（留空即使用默认）
paramsOverride = struct();
% paramsOverride.verbosity = 2;
% paramsOverride.saveResults = true;
% paramsOverride.savePath = 'results_demo.mat';

[params, state, channelInfo] = Channel_model(paramsOverride);
validateParameters(params);
sanityCheckState(state, params);

initialSnapshot = struct('X', state.X, 'theta', state.theta, 'phi', state.phi, 'users', state.users, 'S', state.S, 'W', state.W);

[state, initInfo] = Initialization(state, params);
sanityCheckState(state, params);

afterInitSnapshot = struct('X', state.X, 'theta', state.theta, 'phi', state.phi, 'S', state.S, 'W', state.W);

[state, aoInfo] = AO_model(state, params);
sanityCheckState(state, params);

[state, signalInfo] = Signal_model(state, params);
problemInfo = Problem_formulation(state, params, aoInfo);

% 统一结果结构：summary / trace / snapshots / runtime / diagnostics
results = struct();
results.summary = struct( ...
    'finalSumRate', signalInfo.sumRate, ...
    'converged', aoInfo.converged, ...
    'iterations', aoInfo.iterations, ...
    'terminationReason', aoInfo.terminationReason, ...
    'constraints', problemInfo.constraints);
results.trace = struct( ...
    'outerSumRateHistory', aoInfo.sumRateHistory, ...
    'iterationTrace', aoInfo.iterationTrace, ...
    'blockHistory', aoInfo.blockHistory); % 完整过程
results.snapshots = struct( ...
    'initial', initialSnapshot, ...
    'afterInit', afterInitSnapshot, ...
    'final', struct('X', state.X, 'theta', state.theta, 'phi', state.phi, 'S', state.S, 'W', state.W, 'sumRate', state.sumRate));
results.runtime = aoInfo.runtimeStats;
results.diagnostics = struct( ...
    'initDiagnostics', initInfo.diagnostics, ...
    'wStopReason', getLastWStopReason(aoInfo), ...
    'sBreakReasonLastIter', getLastSBreakReason(aoInfo), ...
    'numericalGuards', collectNumericalGuardInfo(aoInfo));

% 保留原有关键信息，兼容旧调用
results.params = params;
results.channelInfo = channelInfo;
results.initInfo = initInfo;
results.aoInfo = aoInfo;
results.signalInfo = signalInfo;
results.problemInfo = problemInfo;

if params.saveResults
    save(params.savePath, 'results');
end

printSummaryAndTrace(results, params.verbosity);


function reason = getLastWStopReason(aoInfo)
    if isempty(aoInfo.blockHistory)
        reason = 'noIteration';
        return;
    end
    reason = aoInfo.blockHistory(end).W.stopReason;
end
function reason = getLastSBreakReason(aoInfo)
    if isempty(aoInfo.blockHistory)
        reason = 'noIteration';
        return;
    end
    reason = aoInfo.blockHistory(end).userSet.breakReason;
end

function ng = collectNumericalGuardInfo(aoInfo)
    if isempty(aoInfo.blockHistory)
        ng = struct('wNumericalGuardTriggered', false);
        return;
    end
    ng = struct('wNumericalGuardTriggered', aoInfo.blockHistory(end).W.numericalGuardTriggered);
end

function printSummaryAndTrace(results, verbosity)
    if verbosity <= 0
        fprintf('sumRate=%.6f, iterations=%d, converged=%d\n', ...
            results.summary.finalSumRate, results.summary.iterations, results.summary.converged);
        return;
    end

    fprintf('===== Summary =====\n');
    fprintf('sumRate=%.6f, iterations=%d, converged=%d, reason=%s\n', ...
        results.summary.finalSumRate, results.summary.iterations, results.summary.converged, results.summary.terminationReason);
    fprintf('runtime total=%.4fs (W=%.4f, angle=%.4f, X=%.4f, S=%.4f)\n', ...
        results.runtime.totalTime, results.runtime.totalWTime, results.runtime.totalAngleTime, ...
        results.runtime.totalXTime, results.runtime.totalSTime);

    if verbosity >= 1
        fprintf('===== Iteration Block Trace =====\n');
        for t = 1:numel(results.trace.iterationTrace)
            tr = results.trace.iterationTrace(t);
            fprintf('t=%d | before=%.4f | W:+%.4e angle:+%.4e X:+%.4e S:+%.4e | afterS=%.4f | stop=%d\n', ...
                tr.iteration, tr.sumRateBefore, tr.deltaW, tr.deltaAngle, tr.deltaX, tr.deltaS, ...
                tr.sumRateAfterS, tr.terminationCheck.stop);
        end
    end

    if verbosity >= 2
        disp('Outer sum-rate history:');
        disp(results.trace.outerSumRateHistory(:).');
        disp('Runtime iteration table:');
        disp(struct2table(results.runtime.iteration));
        disp('Trace fields available: iterationTrace / blockHistory / snapshots / diagnostics');
    end
end
