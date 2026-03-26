%% Main.m
% 工程化 demo 入口：单次运行 + 全量过程追踪导出

clear; clc;

% 工程化增强：可在此处覆写默认参数（留空即使用默认）
paramsOverride = struct();
paramsOverride.verbosity = 1;

% 更友好的场景强度
paramsOverride.Pmax = 10;
paramsOverride.sigma2 = 1e-6;
paramsOverride.alphaW = 0.03;
paramsOverride.alphaL = 0.985;

% 静态候选池快速调试配置（中等规模）：NRF = N, KServ = NRF = N
paramsOverride.N = 8;
paramsOverride.M = 6;
paramsOverride.K = 64;
paramsOverride.NRF = paramsOverride.N;
paramsOverride.Kmax = paramsOverride.N;
paramsOverride.KServ = paramsOverride.N;

% 保持 clustered 用户分布
paramsOverride.userGeneration = struct();
paramsOverride.userGeneration.mode = 'clustered';
paramsOverride.userGeneration.hotspotCenters = [3, 8; 7, 10];
paramsOverride.userGeneration.hotspotStd = [0.8, 1.0];
paramsOverride.userRegionY = [6, 12];

% 给 AO 和 WMMSE 多一点迭代空间
paramsOverride.IW = 60;
paramsOverride.Tmax = 20;
paramsOverride.epsilonOuter = 1e-6;
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
            bh = getBlockHistoryEntry(results.trace.blockHistory, t);
            sDiag = summarizeSBlock(bh, results.params);
            fprintf('  S: trig=%d eval=%d refine=%d posC=%d posF=%d aboveEps=%d bestDc=%.2e bestDf=%.2e margin=%.2e pair=(%s->%s) break=%s\n', ...
                sDiag.triggered, sDiag.evaluatedPairs, sDiag.refinedPairs, sDiag.positiveCoarse, ...
                sDiag.positiveFinal, sDiag.aboveEpsilonS, sDiag.bestDeltaCoarse, sDiag.bestDeltaFinal, ...
                sDiag.bestMargin, sDiag.bestWeakUser, sDiag.bestStrongUser, sDiag.breakReason);
            fprintf('  Sdet: weak=%s strong=%s pool=%s scores=%s\n', ...
                sDiag.weakUserSetStr, sDiag.strongUserSetStr, sDiag.dynamicCandidatePoolHeadStr, sDiag.strongExternalScoresStr);
            xDiag = summarizeXBlock(bh);
            fprintf('  X: acc=%d ls=%d projCollapse=%d actB=%d actGap=%d bestWG=%d bestAlpha=%.2e bestCoarse=%.2e bestFinal=%.2e rejectMain=%s\n', ...
                xDiag.acceptedCount, xDiag.totalLineSearchSteps, xDiag.projectionCollapsedCount, ...
                xDiag.activeBoundaryCount, xDiag.activeSpacingCount, xDiag.bestWaveguide, xDiag.bestAlpha, ...
                xDiag.bestCoarseImprove, xDiag.bestFinalImprove, xDiag.dominantRejectReason);
        end
    end

    if verbosity >= 2
        disp('Outer sum-rate history:');
        disp(results.trace.outerSumRateHistory(:).');
        disp('Runtime iteration table:');
        disp(struct2table(results.runtime.iteration));
        disp('S diagnostics table:');
        disp(buildSTraceTable(results));
        disp('X diagnostics table:');
        disp(buildXTraceTable(results));
        disp('Trace fields available: iterationTrace / blockHistory / snapshots / diagnostics');
    end
end

function bh = getBlockHistoryEntry(blockHistory, idx)
    if isempty(blockHistory) || idx > numel(blockHistory)
        bh = struct();
    else
        bh = blockHistory(idx);
    end
end

function sDiag = summarizeSBlock(blockHistoryEntry, params)
    sDiag = struct( ...
        'triggered', 0, ...
        'evaluatedPairs', 0, ...
        'refinedPairs', 0, ...
        'positiveCoarse', 0, ...
        'positiveFinal', 0, ...
        'aboveEpsilonS', 0, ...
        'bestDeltaCoarse', NaN, ...
        'bestDeltaFinal', NaN, ...
        'bestMargin', NaN, ...
        'bestWeakUser', 'NaN', ...
        'bestStrongUser', 'NaN', ...
        'weakUserSetStr', '[]', ...
        'strongUserSetStr', '[]', ...
        'dynamicCandidatePoolHeadStr', '[]', ...
        'strongExternalScoresStr', '[]', ...
        'acceptedSwaps', 0, ...
        'breakReason', 'NA');

    if ~isstruct(blockHistoryEntry) || ~isfield(blockHistoryEntry, 'userSet') || isempty(blockHistoryEntry.userSet)
        return;
    end
    userSet = blockHistoryEntry.userSet;
    if isfield(userSet, 'triggered')
        sDiag.triggered = double(logical(userSet.triggered));
    elseif isfield(userSet, 'swapTrace')
        sDiag.triggered = double(~isempty(userSet.swapTrace));
    end
    if isfield(userSet, 'acceptedSwaps') && ~isempty(userSet.acceptedSwaps)
        sDiag.acceptedSwaps = userSet.acceptedSwaps;
    end
    if isfield(userSet, 'breakReason') && ~isempty(userSet.breakReason)
        sDiag.breakReason = char(string(userSet.breakReason));
    end
    if ~isfield(userSet, 'swapTrace') || isempty(userSet.swapTrace)
        return;
    end

    swapTrace = userSet.swapTrace;
    evaluated = [swapTrace.numCandidatesEvaluated];
    refined = [swapTrace.numCandidatesRefined];
    posC = [swapTrace.numCandidatesWithPositiveCoarseDelta];
    posF = [swapTrace.numCandidatesWithPositiveFinalDelta];
    above = [swapTrace.numCandidatesAboveEpsilonS];
    sDiag.evaluatedPairs = sum(evaluated(isfinite(evaluated)));
    sDiag.refinedPairs = sum(refined(isfinite(refined)));
    sDiag.positiveCoarse = sum(posC(isfinite(posC)));
    sDiag.positiveFinal = sum(posF(isfinite(posF)));
    sDiag.aboveEpsilonS = sum(above(isfinite(above)));

    bestDc = [swapTrace.bestDeltaCoarse];
    bestDf = [swapTrace.bestDeltaFinal];
    bestMargin = [swapTrace.bestDeltaMinusEpsilonS];
    if any(isfinite(bestDc))
        sDiag.bestDeltaCoarse = max(bestDc(isfinite(bestDc)));
    end
    if any(isfinite(bestDf))
        [sDiag.bestDeltaFinal, idxBest] = max(bestDf(isfinite(bestDf)));
        finiteIdx = find(isfinite(bestDf));
        idxBest = finiteIdx(idxBest);
        if isfield(swapTrace(idxBest), 'bestPair') && ~isempty(swapTrace(idxBest).bestPair)
            bp = swapTrace(idxBest).bestPair;
            if isfield(bp, 'weakUser') && ~isempty(bp.weakUser)
                sDiag.bestWeakUser = num2str(bp.weakUser);
            end
            if isfield(bp, 'strongUser') && ~isempty(bp.strongUser)
                sDiag.bestStrongUser = num2str(bp.strongUser);
            end
        end
    end
    if any(isfinite(bestMargin))
        sDiag.bestMargin = max(bestMargin(isfinite(bestMargin)));
    elseif isfinite(sDiag.bestDeltaFinal) && isfield(params, 'epsilonS')
        sDiag.bestMargin = sDiag.bestDeltaFinal - params.epsilonS;
    end

    idxDetail = numel(swapTrace);
    if any(isfinite(bestDf))
        [~, idxBestFinite] = max(bestDf(isfinite(bestDf)));
        finiteIdx = find(isfinite(bestDf));
        idxDetail = finiteIdx(idxBestFinite);
    end
    detailEntry = swapTrace(idxDetail);
    if isfield(detailEntry, 'weakUserSet') && ~isempty(detailEntry.weakUserSet)
        sDiag.weakUserSetStr = formatNumericVector(detailEntry.weakUserSet, '%d');
    end
    if isfield(detailEntry, 'strongUserSet') && ~isempty(detailEntry.strongUserSet)
        sDiag.strongUserSetStr = formatNumericVector(detailEntry.strongUserSet, '%d');
    end
    if isfield(detailEntry, 'dynamicCandidatePoolHead') && ~isempty(detailEntry.dynamicCandidatePoolHead)
        sDiag.dynamicCandidatePoolHeadStr = formatNumericVector(detailEntry.dynamicCandidatePoolHead, '%d');
    end
    if isfield(detailEntry, 'strongExternalScores') && ~isempty(detailEntry.strongExternalScores)
        sDiag.strongExternalScoresStr = formatNumericVector(detailEntry.strongExternalScores, '%.2e');
    end
end

function xDiag = summarizeXBlock(blockHistoryEntry)
    xDiag = struct( ...
        'acceptedCount', 0, ...
        'totalLineSearchSteps', 0, ...
        'projectionCollapsedCount', 0, ...
        'activeBoundaryCount', 0, ...
        'activeSpacingCount', 0, ...
        'bestWaveguide', NaN, ...
        'bestAlpha', NaN, ...
        'bestCoarseImprove', NaN, ...
        'bestFinalImprove', NaN, ...
        'dominantRejectReason', 'NA');

    if ~isstruct(blockHistoryEntry) || ~isfield(blockHistoryEntry, 'position') || isempty(blockHistoryEntry.position)
        return;
    end
    posInfo = blockHistoryEntry.position;
    if isfield(posInfo, 'acceptedCount') && ~isempty(posInfo.acceptedCount)
        xDiag.acceptedCount = posInfo.acceptedCount;
    end
    if ~isfield(posInfo, 'waveguideInfo') || isempty(posInfo.waveguideInfo)
        return;
    end
    wg = posInfo.waveguideInfo;
    lsSteps = [wg.lineSearchSteps];
    projCollapse = [wg.projectionCollapsedCount];
    xDiag.totalLineSearchSteps = sum(lsSteps(isfinite(lsSteps)));
    xDiag.projectionCollapsedCount = sum(projCollapse(isfinite(projCollapse)));
    if isfield(wg, 'activeBoundaryConstraint')
        xDiag.activeBoundaryCount = sum([wg.activeBoundaryConstraint]);
    end
    if isfield(wg, 'activeSpacingConstraint')
        xDiag.activeSpacingCount = sum([wg.activeSpacingConstraint]);
    end

    coarseImprove = [wg.bestCoarseRate] - [wg.sumRateBefore];
    finalImprove = [wg.bestFinalRate] - [wg.sumRateBefore];
    if any(isfinite(finalImprove))
        [xDiag.bestFinalImprove, idxBest] = max(finalImprove(isfinite(finalImprove)));
        finiteIdx = find(isfinite(finalImprove));
        idxBest = finiteIdx(idxBest);
        if isfield(wg(idxBest), 'waveguideIndex') && ~isempty(wg(idxBest).waveguideIndex)
            xDiag.bestWaveguide = wg(idxBest).waveguideIndex;
        else
            xDiag.bestWaveguide = idxBest;
        end
        if isfield(wg(idxBest), 'alphaAccepted') && ~isempty(wg(idxBest).alphaAccepted)
            xDiag.bestAlpha = wg(idxBest).alphaAccepted;
        end
    end
    if any(isfinite(coarseImprove))
        xDiag.bestCoarseImprove = max(coarseImprove(isfinite(coarseImprove)));
    end

    reasons = {};
    for i = 1:numel(wg)
        if isfield(wg(i), 'accepted') && wg(i).accepted
            continue;
        end
        if isfield(wg(i), 'rejectReason') && ~isempty(wg(i).rejectReason)
            reasons{end + 1} = char(string(wg(i).rejectReason)); %#ok<AGROW>
        end
    end
    if ~isempty(reasons)
        [uniqReasons, ~, ic] = unique(reasons);
        counts = accumarray(ic(:), 1);
        [~, imax] = max(counts);
        xDiag.dominantRejectReason = uniqReasons{imax};
    end
end

function sTable = buildSTraceTable(results)
    nIter = numel(results.trace.iterationTrace);
    iter = (1:nIter).';
    triggered = zeros(nIter, 1);
    evaluatedPairs = zeros(nIter, 1);
    refinedPairs = zeros(nIter, 1);
    positiveCoarse = zeros(nIter, 1);
    positiveFinal = zeros(nIter, 1);
    aboveEpsilonS = zeros(nIter, 1);
    bestDeltaCoarse = nan(nIter, 1);
    bestDeltaFinal = nan(nIter, 1);
    bestMargin = nan(nIter, 1);
    acceptedSwaps = zeros(nIter, 1);
    breakReason = strings(nIter, 1);
    weakUserSet = strings(nIter, 1);
    strongUserSet = strings(nIter, 1);
    dynamicCandidatePoolHead = strings(nIter, 1);
    strongExternalScores = strings(nIter, 1);
    for t = 1:nIter
        bh = getBlockHistoryEntry(results.trace.blockHistory, t);
        sDiag = summarizeSBlock(bh, results.params);
        triggered(t) = sDiag.triggered;
        evaluatedPairs(t) = sDiag.evaluatedPairs;
        refinedPairs(t) = sDiag.refinedPairs;
        positiveCoarse(t) = sDiag.positiveCoarse;
        positiveFinal(t) = sDiag.positiveFinal;
        aboveEpsilonS(t) = sDiag.aboveEpsilonS;
        bestDeltaCoarse(t) = sDiag.bestDeltaCoarse;
        bestDeltaFinal(t) = sDiag.bestDeltaFinal;
        bestMargin(t) = sDiag.bestMargin;
        acceptedSwaps(t) = sDiag.acceptedSwaps;
        breakReason(t) = string(sDiag.breakReason);
        weakUserSet(t) = string(sDiag.weakUserSetStr);
        strongUserSet(t) = string(sDiag.strongUserSetStr);
        dynamicCandidatePoolHead(t) = string(sDiag.dynamicCandidatePoolHeadStr);
        strongExternalScores(t) = string(sDiag.strongExternalScoresStr);
    end
    sTable = table(iter, triggered, evaluatedPairs, refinedPairs, positiveCoarse, ...
        positiveFinal, aboveEpsilonS, bestDeltaCoarse, bestDeltaFinal, bestMargin, ...
        acceptedSwaps, breakReason, weakUserSet, strongUserSet, ...
        dynamicCandidatePoolHead, strongExternalScores);
end

function xTable = buildXTraceTable(results)
    nIter = numel(results.trace.iterationTrace);
    iter = (1:nIter).';
    acceptedCount = zeros(nIter, 1);
    totalLineSearchSteps = zeros(nIter, 1);
    projectionCollapsedCount = zeros(nIter, 1);
    activeBoundaryCount = zeros(nIter, 1);
    activeSpacingCount = zeros(nIter, 1);
    bestWaveguide = nan(nIter, 1);
    bestAlpha = nan(nIter, 1);
    bestCoarseImprove = nan(nIter, 1);
    bestFinalImprove = nan(nIter, 1);
    dominantRejectReason = strings(nIter, 1);
    for t = 1:nIter
        bh = getBlockHistoryEntry(results.trace.blockHistory, t);
        xDiag = summarizeXBlock(bh);
        acceptedCount(t) = xDiag.acceptedCount;
        totalLineSearchSteps(t) = xDiag.totalLineSearchSteps;
        projectionCollapsedCount(t) = xDiag.projectionCollapsedCount;
        activeBoundaryCount(t) = xDiag.activeBoundaryCount;
        activeSpacingCount(t) = xDiag.activeSpacingCount;
        bestWaveguide(t) = xDiag.bestWaveguide;
        bestAlpha(t) = xDiag.bestAlpha;
        bestCoarseImprove(t) = xDiag.bestCoarseImprove;
        bestFinalImprove(t) = xDiag.bestFinalImprove;
        dominantRejectReason(t) = string(xDiag.dominantRejectReason);
    end
    xTable = table(iter, acceptedCount, totalLineSearchSteps, projectionCollapsedCount, ...
        activeBoundaryCount, activeSpacingCount, bestWaveguide, bestAlpha, ...
        bestCoarseImprove, bestFinalImprove, dominantRejectReason);
end

function out = formatNumericVector(vec, fmt)
    if nargin < 2 || isempty(fmt)
        fmt = '%g';
    end
    if isempty(vec)
        out = '[]';
        return;
    end
    vec = vec(:).';
    parts = arrayfun(@(x) sprintf(fmt, x), vec, 'UniformOutput', false);
    out = ['[', strjoin(parts, ' '), ']'];
end
