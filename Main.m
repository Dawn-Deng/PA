%% Main.m
% 工程化 demo 入口：单次运行 + 可导出 results

clear; clc;

% 工程化增强：可在此处覆写默认参数（留空即使用默认）
paramsOverride = struct();
% paramsOverride.verbosity = 2;
% paramsOverride.saveResults = true;
% paramsOverride.savePath = 'results_demo.mat';

[params, state, channelInfo] = Channel_model(paramsOverride);
validateParameters(params); % 工程化增强：统一参数校验
sanityCheckState(state, params);

initialSnapshot = struct('X', state.X, 'theta', state.theta, 'phi', state.phi, 'users', state.users, 'S', state.S);

[state, initInfo] = Initialization(state, params);
sanityCheckState(state, params);

afterInitSnapshot = struct('X', state.X, 'theta', state.theta, 'phi', state.phi, 'S', state.S);

[state, aoInfo] = AO_model(state, params);
sanityCheckState(state, params);

[state, signalInfo] = Signal_model(state, params);
problemInfo = Problem_formulation(state, params, aoInfo);

results = struct();
results.params = params;
results.channelInfo = channelInfo;
results.initInfo = initInfo;
results.aoInfo = aoInfo;
results.signalInfo = signalInfo;
results.problemInfo = problemInfo;
results.runtimeStats = aoInfo.runtimeStats; % 结果导出/实验支持
results.snapshots = struct( ...
    'initial', initialSnapshot, ...
    'afterInit', afterInitSnapshot, ...
    'final', struct('X', state.X, 'theta', state.theta, 'phi', state.phi, 'S', state.S, 'sumRate', state.sumRate));
results.plotData = struct( ...
    'sumRateHistory', aoInfo.sumRateHistory, ...
    'finalSINR', signalInfo.sinr, ...
    'finalRate', signalInfo.rate, ...
    'blockAccepted', struct( ...
        'angle', [aoInfo.runtimeStats.iteration.angleAcceptedCount], ...
        'x', [aoInfo.runtimeStats.iteration.XAcceptedCount], ...
        's', [aoInfo.runtimeStats.iteration.SAcceptedSwaps]));

if params.saveResults
    save(params.savePath, 'results'); % 结果导出/实验支持
end

printSummary(results, params.verbosity);

function printSummary(results, verbosity)
    if verbosity <= 0
        fprintf('sumRate=%.6f, iterations=%d, converged=%d\n', ...
            results.signalInfo.sumRate, results.aoInfo.iterations, results.aoInfo.converged);
        return;
    end

    fprintf('===== Summary =====\n');
    fprintf('sumRate=%.6f, iterations=%d, converged=%d, reason=%s\n', ...
        results.signalInfo.sumRate, results.aoInfo.iterations, results.aoInfo.converged, results.aoInfo.terminationReason);
    fprintf('candidatePoolSize=%d, init utility mean=%.4g\n', ...
        results.initInfo.diagnostics.candidatePoolSize, results.initInfo.diagnostics.utility.mean);
    fprintf('runtime total=%.4fs (W=%.4f, angle=%.4f, X=%.4f, S=%.4f)\n', ...
        results.runtimeStats.totalTime, results.runtimeStats.totalWTime, results.runtimeStats.totalAngleTime, ...
        results.runtimeStats.totalXTime, results.runtimeStats.totalSTime);

    if verbosity >= 2
        disp('AO sum-rate history:');
        disp(results.aoInfo.sumRateHistory(:).');
        disp('Per-iteration block stats:');
        disp(struct2table(results.runtimeStats.iteration));
        disp('Final user metrics:');
        if ~isempty(results.signalInfo.userMetrics)
            disp(struct2table(results.signalInfo.userMetrics));
        else
            disp('No scheduled users / empty metrics.');
        end
    end
end
