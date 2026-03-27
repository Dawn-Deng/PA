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
% 可选：结果可视化开关（默认关闭，不影响原流程）
paramsOverride.enablePlotResults = false;
paramsOverride.plotSaveDir = '';
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

% 可选调用：集中式结果可视化（默认关闭，保持最小侵入）
if isfield(params, 'enablePlotResults') && params.enablePlotResults
    try
        if isfield(params, 'plotSaveDir') && ~isempty(params.plotSaveDir)
            plot_all_results(results, params.plotSaveDir);
        else
            plot_all_results(results);
        end
    catch ME
        warning('Main:PlotAllResultsFailed', 'plot_all_results failed: %s', ME.message);
    end
end


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
            fprintf('  Sdet: weak=%s strong=%s pool=%s base=%s dyn=%s mode=%s stag=%d intens=%d twoSwap=%d/%d postRef=%d(+%.2e)\n', ...
                sDiag.weakUserSetStr, sDiag.strongUserSetStr, sDiag.dynamicCandidatePoolHeadStr, ...
                sDiag.strongExternalScoresStr, sDiag.dynamicScoreHeadStr, sDiag.scoreMode, ...
                sDiag.stagnationLevel, sDiag.intensificationTriggered, sDiag.limitedTwoSwapTried, ...
                sDiag.limitedTwoSwapAccepted, sDiag.postAcceptShortRefineTriggered, sDiag.postAcceptGain);
            xDiag = summarizeXBlock(bh);
            fprintf('  X: acc=%d ls=%d projCollapse=%d actB=%d actGap=%d bestWG=%d bestAlpha=%.2e bestCoarse=%.2e bestFinal=%.2e rejectMain=%s\n', ...
                xDiag.acceptedCount, xDiag.totalLineSearchSteps, xDiag.projectionCollapsedCount, ...
                xDiag.activeBoundaryCount, xDiag.activeSpacingCount, xDiag.bestWaveguide, xDiag.bestAlpha, ...
                xDiag.bestCoarseImprove, xDiag.bestFinalImprove, xDiag.dominantRejectReason);
            fprintf('  Xdet: gN=%.2e gDotD=%.2e rawDN=%.2e dirN=%.2e bClip=%d sClip=%d warm=%d aInit=%.2e aEma=%.2e projCorr=%.2e\n', ...
                xDiag.meanGradNorm, xDiag.meanGradDotDirection, xDiag.meanRawDirectionNorm, xDiag.meanDirectionNorm, ...
                xDiag.boundaryCleanupCount, xDiag.spacingCleanupCount, xDiag.warmStartUsedCount, ...
                xDiag.meanAlphaInit, xDiag.meanAlphaEma, xDiag.meanProjectionCorrection);
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
        'scoreMode', 'NA', ...
        'dynamicScoreHeadStr', '[]', ...
        'topRefineCoarseRanksStr', '[]', ...
        'topRefineFinalRanksStr', '[]', ...
        'stagnationLevel', 0, ...
        'intensificationTriggered', 0, ...
        'limitedTwoSwapTried', 0, ...
        'limitedTwoSwapAccepted', 0, ...
        'postAcceptShortRefineTriggered', 0, ...
        'postAcceptGain', 0, ...
        'tabuHeadPairsStr', '[]', ...
        'penaltyHeadPairsStr', '[]', ...
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
    if isfield(detailEntry, 'scoreMode') && ~isempty(detailEntry.scoreMode)
        sDiag.scoreMode = char(string(detailEntry.scoreMode));
    end
    if isfield(detailEntry, 'dynamicScoreHead') && ~isempty(detailEntry.dynamicScoreHead)
        sDiag.dynamicScoreHeadStr = formatNumericVector(detailEntry.dynamicScoreHead, '%.2e');
    end
    if isfield(detailEntry, 'topRefineCoarseRanks') && ~isempty(detailEntry.topRefineCoarseRanks)
        sDiag.topRefineCoarseRanksStr = formatNumericVector(detailEntry.topRefineCoarseRanks, '%d');
    end
    if isfield(detailEntry, 'topRefineFinalRanks') && ~isempty(detailEntry.topRefineFinalRanks)
        sDiag.topRefineFinalRanksStr = formatNumericVector(detailEntry.topRefineFinalRanks, '%d');
    end
    if isfield(detailEntry, 'stagnationLevel') && ~isempty(detailEntry.stagnationLevel)
        sDiag.stagnationLevel = detailEntry.stagnationLevel;
    end
    if isfield(detailEntry, 'intensificationTriggered') && ~isempty(detailEntry.intensificationTriggered)
        sDiag.intensificationTriggered = double(logical(detailEntry.intensificationTriggered));
    end
    if isfield(detailEntry, 'limitedTwoSwapTried') && ~isempty(detailEntry.limitedTwoSwapTried)
        sDiag.limitedTwoSwapTried = double(logical(detailEntry.limitedTwoSwapTried));
    end
    if isfield(detailEntry, 'limitedTwoSwapAccepted') && ~isempty(detailEntry.limitedTwoSwapAccepted)
        sDiag.limitedTwoSwapAccepted = double(logical(detailEntry.limitedTwoSwapAccepted));
    end
    if isfield(detailEntry, 'postAcceptShortRefineTriggered') && ~isempty(detailEntry.postAcceptShortRefineTriggered)
        sDiag.postAcceptShortRefineTriggered = double(logical(detailEntry.postAcceptShortRefineTriggered));
    end
    if isfield(detailEntry, 'postAcceptGain') && ~isempty(detailEntry.postAcceptGain)
        sDiag.postAcceptGain = detailEntry.postAcceptGain;
    end
    if isfield(detailEntry, 'tabuHeadPairs') && ~isempty(detailEntry.tabuHeadPairs)
        sDiag.tabuHeadPairsStr = mat2str(detailEntry.tabuHeadPairs);
    end
    if isfield(detailEntry, 'penaltyHeadPairs') && ~isempty(detailEntry.penaltyHeadPairs)
        sDiag.penaltyHeadPairsStr = formatNumericVector(detailEntry.penaltyHeadPairs, '%.2e');
    end
end

function xDiag = summarizeXBlock(blockHistoryEntry)
    xDiag = struct( ...
        'acceptedCount', 0, ...
        'totalLineSearchSteps', 0, ...
        'projectionCollapsedCount', 0, ...
        'activeBoundaryCount', 0, ...
        'activeSpacingCount', 0, ...
        'meanGradNorm', NaN, ...
        'meanGradDotDirection', NaN, ...
        'meanRawDirectionNorm', NaN, ...
        'meanDirectionNorm', NaN, ...
        'boundaryCleanupCount', 0, ...
        'spacingCleanupCount', 0, ...
        'warmStartUsedCount', 0, ...
        'meanAlphaInit', NaN, ...
        'meanAlphaEma', NaN, ...
        'meanProjectionCorrection', NaN, ...
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
    if isfield(wg, 'gradNorm')
        vals = [wg.gradNorm];
        xDiag.meanGradNorm = mean(vals(isfinite(vals)));
    end
    if isfield(wg, 'gradDotDirection')
        vals = [wg.gradDotDirection];
        xDiag.meanGradDotDirection = mean(vals(isfinite(vals)));
    end
    if isfield(wg, 'rawDirectionNorm')
        vals = [wg.rawDirectionNorm];
        xDiag.meanRawDirectionNorm = mean(vals(isfinite(vals)));
    end
    if isfield(wg, 'directionNorm')
        vals = [wg.directionNorm];
        xDiag.meanDirectionNorm = mean(vals(isfinite(vals)));
    end
    if isfield(wg, 'boundaryOutwardComponentsClipped')
        xDiag.boundaryCleanupCount = sum([wg.boundaryOutwardComponentsClipped]);
    end
    if isfield(wg, 'spacingConflictPairsClipped')
        xDiag.spacingCleanupCount = sum([wg.spacingConflictPairsClipped]);
    end
    if isfield(wg, 'alphaWarmStartUsed')
        xDiag.warmStartUsedCount = sum([wg.alphaWarmStartUsed]);
    end
    if isfield(wg, 'alphaInit')
        vals = [wg.alphaInit];
        xDiag.meanAlphaInit = mean(vals(isfinite(vals)));
    end
    if isfield(wg, 'alphaHistoryEma')
        vals = [wg.alphaHistoryEma];
        xDiag.meanAlphaEma = mean(vals(isfinite(vals)));
    end
    projCorr = [];
    for i = 1:numel(wg)
        if isfield(wg(i), 'lineSearchTrace') && ~isempty(wg(i).lineSearchTrace) && isfield(wg(i).lineSearchTrace, 'projectionCorrection')
            projCorr = [projCorr, [wg(i).lineSearchTrace.projectionCorrection]]; %#ok<AGROW>
        end
    end
    if any(isfinite(projCorr))
        xDiag.meanProjectionCorrection = mean(projCorr(isfinite(projCorr)));
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
    scoreMode = strings(nIter, 1);
    dynamicScoreHead = strings(nIter, 1);
    stagnationLevel = zeros(nIter, 1);
    intensificationTriggered = zeros(nIter, 1);
    limitedTwoSwapTried = zeros(nIter, 1);
    limitedTwoSwapAccepted = zeros(nIter, 1);
    postAcceptShortRefineTriggered = zeros(nIter, 1);
    postAcceptGain = zeros(nIter, 1);
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
        scoreMode(t) = string(sDiag.scoreMode);
        dynamicScoreHead(t) = string(sDiag.dynamicScoreHeadStr);
        stagnationLevel(t) = sDiag.stagnationLevel;
        intensificationTriggered(t) = sDiag.intensificationTriggered;
        limitedTwoSwapTried(t) = sDiag.limitedTwoSwapTried;
        limitedTwoSwapAccepted(t) = sDiag.limitedTwoSwapAccepted;
        postAcceptShortRefineTriggered(t) = sDiag.postAcceptShortRefineTriggered;
        postAcceptGain(t) = sDiag.postAcceptGain;
    end
    sTable = table(iter, triggered, evaluatedPairs, refinedPairs, positiveCoarse, ...
        positiveFinal, aboveEpsilonS, bestDeltaCoarse, bestDeltaFinal, bestMargin, ...
        acceptedSwaps, breakReason, weakUserSet, strongUserSet, ...
        dynamicCandidatePoolHead, strongExternalScores, scoreMode, dynamicScoreHead, ...
        stagnationLevel, intensificationTriggered, limitedTwoSwapTried, ...
        limitedTwoSwapAccepted, postAcceptShortRefineTriggered, postAcceptGain);
end

function xTable = buildXTraceTable(results)
    nIter = numel(results.trace.iterationTrace);
    iter = (1:nIter).';
    acceptedCount = zeros(nIter, 1);
    totalLineSearchSteps = zeros(nIter, 1);
    projectionCollapsedCount = zeros(nIter, 1);
    activeBoundaryCount = zeros(nIter, 1);
    activeSpacingCount = zeros(nIter, 1);
    meanGradNorm = nan(nIter, 1);
    meanGradDotDirection = nan(nIter, 1);
    meanRawDirectionNorm = nan(nIter, 1);
    meanDirectionNorm = nan(nIter, 1);
    boundaryCleanupCount = zeros(nIter, 1);
    spacingCleanupCount = zeros(nIter, 1);
    warmStartUsedCount = zeros(nIter, 1);
    meanAlphaInit = nan(nIter, 1);
    meanAlphaEma = nan(nIter, 1);
    meanProjectionCorrection = nan(nIter, 1);
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
        meanGradNorm(t) = xDiag.meanGradNorm;
        meanGradDotDirection(t) = xDiag.meanGradDotDirection;
        meanRawDirectionNorm(t) = xDiag.meanRawDirectionNorm;
        meanDirectionNorm(t) = xDiag.meanDirectionNorm;
        boundaryCleanupCount(t) = xDiag.boundaryCleanupCount;
        spacingCleanupCount(t) = xDiag.spacingCleanupCount;
        warmStartUsedCount(t) = xDiag.warmStartUsedCount;
        meanAlphaInit(t) = xDiag.meanAlphaInit;
        meanAlphaEma(t) = xDiag.meanAlphaEma;
        meanProjectionCorrection(t) = xDiag.meanProjectionCorrection;
        bestWaveguide(t) = xDiag.bestWaveguide;
        bestAlpha(t) = xDiag.bestAlpha;
        bestCoarseImprove(t) = xDiag.bestCoarseImprove;
        bestFinalImprove(t) = xDiag.bestFinalImprove;
        dominantRejectReason(t) = string(xDiag.dominantRejectReason);
    end
    xTable = table(iter, acceptedCount, totalLineSearchSteps, projectionCollapsedCount, ...
        activeBoundaryCount, activeSpacingCount, meanGradNorm, meanGradDotDirection, ...
        meanRawDirectionNorm, meanDirectionNorm, boundaryCleanupCount, spacingCleanupCount, ...
        warmStartUsedCount, meanAlphaInit, meanAlphaEma, meanProjectionCorrection, ...
        bestWaveguide, bestAlpha, ...
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
