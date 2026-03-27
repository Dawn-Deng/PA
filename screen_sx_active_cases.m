function [rankingTable, allCases] = screen_sx_active_cases(config)
%SCREEN_SX_ACTIVE_CASES 批量筛选“更能激活 S-block / X-block”的 case。
%
% 方式1：默认筛选
% [rankingTable, allCases] = screen_sx_active_cases();
%
% 方式2：自定义 seeds
% cfg = struct();
% cfg.seeds = 20260321:20260328;
% [rankingTable, allCases] = screen_sx_active_cases(cfg);
%
% 方式3：保存结果
% cfg = struct();
% cfg.saveResults = true;
% cfg.saveDir = 'screen_outputs';
% [rankingTable, allCases] = screen_sx_active_cases(cfg);

    if nargin < 1
        config = struct();
    end
    cfg = default_screen_config(config);

    % 说明：run_experiments.m 当前输出结构和 Main.m 最新 results 不一致，
    % 这里采用本文件内的 current pipeline 封装，保证字段与 Main.m 对齐。
    caseList = build_case_list(cfg);
    nCases = numel(caseList);

    allCases = repmat(struct( ...
        'caseId', '', ...
        'meta', struct(), ...
        'override', struct(), ...
        'success', false, ...
        'errorMessage', '', ...
        'results', [], ...
        'features', struct(), ...
        'screenScore', -inf), nCases, 1);

    fprintf('===== screen_sx_active_cases: start =====\n');
    fprintf('Total cases to run: %d\n', nCases);

    for i = 1:nCases
        one = caseList(i);
        allCases(i).caseId = one.caseId;
        allCases(i).meta = one.meta;
        allCases(i).override = one.override;

        try
            results = run_single_case_current_pipeline(one.override);
            feat = extract_case_features(results);
            allCases(i).results = results;
            allCases(i).features = feat;
            allCases(i).success = true;
            fprintf('[%3d/%3d] OK  | %s | finalSumRate=%.4f | Sswap=%d | Xacc=%d\n', ...
                i, nCases, one.caseId, feat.finalSumRate, feat.S_acceptedSwaps_total, feat.X_acceptedCount_total);
        catch ME
            allCases(i).success = false;
            allCases(i).errorMessage = ME.message;
            warning('screen_sx_active_cases:CaseFailed', '[%d/%d] Case failed (%s): %s', i, nCases, one.caseId, ME.message);
        end
    end

    rankingTable = build_ranking_table(allCases, cfg);
    allCases = attach_scores_to_cases(allCases, rankingTable);

    print_screen_summary(rankingTable, allCases, cfg);

    if cfg.saveResults
        save_screen_outputs(rankingTable, allCases, cfg);
    end
end

function cfg = default_screen_config(userCfg)
    cfg = struct();

    % 默认基准参数：与 Main.m 当前 demo 对齐
    cfg.baseParams = struct();
    cfg.baseParams.verbosity = 0;
    cfg.baseParams.enablePlotResults = false;
    cfg.baseParams.Pmax = 10;
    cfg.baseParams.sigma2 = 1e-6;
    cfg.baseParams.alphaW = 0.03;
    cfg.baseParams.alphaL = 0.985;
    cfg.baseParams.N = 8;
    cfg.baseParams.M = 6;
    cfg.baseParams.K = 64;
    cfg.baseParams.NRF = 8;
    cfg.baseParams.Kmax = 8;
    cfg.baseParams.KServ = 8;
    cfg.baseParams.IW = 60;
    cfg.baseParams.Tmax = 20;
    cfg.baseParams.epsilonOuter = 1e-6;
    cfg.baseParams.userRegionY = [6, 12];

    cfg.seeds = 20260321:20260328;
    cfg.casePresets = default_case_presets();
    cfg.sweep = struct('field', {}, 'values', {});

    cfg.topKPrint = 10;
    cfg.saveResults = false;
    cfg.saveDir = 'screen_outputs';
    cfg.saveTopKCases = 3;

    % 评分权重：S/X 活跃度优先，sumRate 次要。
    cfg.scoreWeights = struct( ...
        'S_acceptedSwaps_total', 0.25, ...
        'S_bestDeltaFinal_max_pos', 0.16, ...
        'X_acceptedCount_total', 0.25, ...
        'X_bestFinalImprove_max_pos', 0.16, ...
        'S_positiveFinal_total', 0.07, ...
        'X_totalLineSearchSteps', 0.07, ...
        'finalSumRate', 0.04);

    % 约束/收敛/异常惩罚
    cfg.penaltyConstraintViolation = 0.40; % 显著降权
    cfg.penaltyNotConverged = 0.08;        % 轻度降权
    cfg.invalidCaseScore = -1e9;           % 直接排后

    cfg = merge_struct(cfg, userCfg);
end

function presets = default_case_presets()
    presets = repmat(struct('name', '', 'override', struct()), 5, 1);

    presets(1).name = 'baseline_clustered';
    presets(1).override = struct('userGeneration', struct( ...
        'mode', 'clustered', 'hotspotCenters', [3, 8; 7, 10], 'hotspotStd', [0.8, 1.0]), ...
        'userRegionY', [6, 12]);

    presets(2).name = 'dual_clustered_separated';
    presets(2).override = struct('userGeneration', struct( ...
        'mode', 'clustered', 'hotspotCenters', [2.2, 7.0; 7.8, 11.0], 'hotspotStd', [0.5, 0.7]), ...
        'userRegionY', [6, 12]);

    presets(3).name = 'triple_clustered';
    presets(3).override = struct('userGeneration', struct( ...
        'mode', 'clustered', 'hotspotCenters', [2.2, 7.0; 5.0, 9.0; 7.8, 11.0], 'hotspotStd', [0.5, 0.7]), ...
        'userRegionY', [6, 12]);

    presets(4).name = 'uniform_control';
    presets(4).override = struct('userGeneration', struct('mode', 'uniform'), 'userRegionY', [6, 12]);

    presets(5).name = 'corridor_control';
    presets(5).override = struct('userGeneration', struct('mode', 'corridor'), 'userRegionY', [6, 12]);
end

function caseList = build_case_list(cfg)
    sweepCases = expand_sweep(cfg.sweep);
    nPresets = numel(cfg.casePresets);
    nSeeds = numel(cfg.seeds);
    nSweep = numel(sweepCases);
    total = nPresets * nSeeds * nSweep;

    caseList = repmat(struct('caseId', '', 'meta', struct(), 'override', struct()), total, 1);
    idx = 0;
    for p = 1:nPresets
        preset = cfg.casePresets(p);
        for s = 1:nSeeds
            seed = cfg.seeds(s);
            for c = 1:nSweep
                idx = idx + 1;
                ov = cfg.baseParams;
                ov.randomSeed = seed;
                ov.symbolSeed = seed + 101;
                ov = merge_struct(ov, preset.override);
                ov = merge_struct(ov, sweepCases(c));

                userMode = nested_get(ov, {'userGeneration', 'mode'}, 'NA');
                caseList(idx).caseId = sprintf('%s_seed%d_sw%d', preset.name, seed, c);
                caseList(idx).meta = struct( ...
                    'seed', seed, ...
                    'symbolSeed', seed + 101, ...
                    'presetName', preset.name, ...
                    'userMode', userMode, ...
                    'sweepIndex', c);
                caseList(idx).override = ov;
            end
        end
    end
end

function sweepCases = expand_sweep(sweep)
    if isempty(sweep)
        sweepCases = struct();
        return;
    end
    sweepCases = struct();
    for i = 1:numel(sweep)
        f = sweep(i).field;
        vals = sweep(i).values;
        if isempty(vals)
            continue;
        end
        if i == 1 || isempty(fieldnames(sweepCases))
            sweepCases = repmat(struct(), numel(vals), 1);
            for j = 1:numel(vals)
                sweepCases(j).(f) = vals(j);
            end
        else
            prev = sweepCases;
            sweepCases = repmat(struct(), numel(prev) * numel(vals), 1);
            k = 0;
            for a = 1:numel(prev)
                for b = 1:numel(vals)
                    k = k + 1;
                    sweepCases(k) = prev(a);
                    sweepCases(k).(f) = vals(b);
                end
            end
        end
    end
    if isempty(sweepCases)
        sweepCases = struct();
    end
end

function results = run_single_case_current_pipeline(paramsOverride)
    [params, state, channelInfo] = Channel_model(paramsOverride);
    validateParameters(params);
    sanityCheckState(state, params);

    initialSnapshot = struct('X', state.X, 'theta', state.theta, 'phi', state.phi, ...
        'users', state.users, 'S', state.S, 'W', state.W);

    [state, initInfo] = Initialization(state, params);
    sanityCheckState(state, params);

    afterInitSnapshot = struct('X', state.X, 'theta', state.theta, 'phi', state.phi, 'S', state.S, 'W', state.W);

    [state, aoInfo] = AO_model(state, params);
    sanityCheckState(state, params);

    [state, signalInfo] = Signal_model(state, params);
    problemInfo = Problem_formulation(state, params, aoInfo);

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
        'blockHistory', aoInfo.blockHistory);
    results.snapshots = struct( ...
        'initial', initialSnapshot, ...
        'afterInit', afterInitSnapshot, ...
        'final', struct('X', state.X, 'theta', state.theta, 'phi', state.phi, 'S', state.S, 'W', state.W, 'sumRate', state.sumRate));
    results.runtime = aoInfo.runtimeStats;
    results.diagnostics = struct( ...
        'initDiagnostics', initInfo.diagnostics, ...
        'wStopReason', get_last_w_stop_reason(aoInfo), ...
        'sBreakReasonLastIter', get_last_s_break_reason(aoInfo), ...
        'numericalGuards', collect_numerical_guard_info(aoInfo));
    results.params = params;
    results.channelInfo = channelInfo;
    results.initInfo = initInfo;
    results.aoInfo = aoInfo;
    results.signalInfo = signalInfo;
    results.problemInfo = problemInfo;
end

function feat = extract_case_features(results)
    feat = struct();

    sAgg = aggregate_s_metrics(results);
    xAgg = aggregate_x_metrics(results);

    feat = merge_struct(feat, sAgg);
    feat = merge_struct(feat, xAgg);

    feat.finalSumRate = nested_get(results, {'summary', 'finalSumRate'}, NaN);
    feat.converged = logical(nested_get(results, {'summary', 'converged'}, false));
    feat.iterations = nested_get(results, {'summary', 'iterations'}, NaN);
    feat.terminationReason = string(nested_get(results, {'summary', 'terminationReason'}, 'NA'));
    feat.constraintsAllSatisfied = logical(nested_get(results, {'summary', 'constraints', 'allSatisfied'}, false));
    feat.runtimeTotal = nested_get(results, {'runtime', 'totalTime'}, NaN);
end

function sAgg = aggregate_s_metrics(results)
    sAgg = struct( ...
        'S_triggered_total', 0, ...
        'S_acceptedSwaps_total', 0, ...
        'S_bestDeltaFinal_max', NaN, ...
        'S_bestDeltaFinal_mean', NaN, ...
        'S_positiveFinal_total', 0, ...
        'S_aboveEpsilonS_total', 0, ...
        'S_postAcceptGain_sum', 0, ...
        'S_intensificationTriggered_total', 0);

    bh = nested_get(results, {'trace', 'blockHistory'}, []);
    if isempty(bh)
        warning('screen_sx_active_cases:NoBlockHistoryS', 'Missing blockHistory when aggregating S metrics.');
        return;
    end

    bestFinalVals = [];
    for t = 1:numel(bh)
        us = nested_get(bh(t), {'userSet'}, []);
        if isempty(us)
            continue;
        end
        sAgg.S_triggered_total = sAgg.S_triggered_total + double(logical(nested_get(us, {'triggered'}, false)));
        sAgg.S_acceptedSwaps_total = sAgg.S_acceptedSwaps_total + safe_num(nested_get(us, {'acceptedSwaps'}, 0));
        sAgg.S_intensificationTriggered_total = sAgg.S_intensificationTriggered_total + double(logical(nested_get(us, {'intensificationTriggered'}, false)));

        swapTrace = nested_get(us, {'swapTrace'}, []);
        if isempty(swapTrace)
            continue;
        end
        for i = 1:numel(swapTrace)
            st = swapTrace(i);
            sAgg.S_positiveFinal_total = sAgg.S_positiveFinal_total + safe_num(nested_get(st, {'numCandidatesWithPositiveFinalDelta'}, 0));
            sAgg.S_aboveEpsilonS_total = sAgg.S_aboveEpsilonS_total + safe_num(nested_get(st, {'numCandidatesAboveEpsilonS'}, 0));
            sAgg.S_postAcceptGain_sum = sAgg.S_postAcceptGain_sum + safe_num(nested_get(st, {'postAcceptGain'}, 0));
            bdf = nested_get(st, {'bestDeltaFinal'}, NaN);
            if isfinite(bdf)
                bestFinalVals(end + 1, 1) = bdf; %#ok<AGROW>
            end
        end
    end

    if ~isempty(bestFinalVals)
        sAgg.S_bestDeltaFinal_max = max(bestFinalVals);
        sAgg.S_bestDeltaFinal_mean = mean(bestFinalVals);
    end
end

function xAgg = aggregate_x_metrics(results)
    xAgg = struct( ...
        'X_acceptedCount_total', 0, ...
        'X_bestFinalImprove_max', NaN, ...
        'X_bestFinalImprove_mean', NaN, ...
        'X_totalLineSearchSteps', 0, ...
        'X_bestCoarseImprove_max', NaN, ...
        'X_projectionCollapsedCount_total', 0, ...
        'X_activeBoundaryCount_total', 0, ...
        'X_activeSpacingCount_total', 0);

    bh = nested_get(results, {'trace', 'blockHistory'}, []);
    if isempty(bh)
        warning('screen_sx_active_cases:NoBlockHistoryX', 'Missing blockHistory when aggregating X metrics.');
        return;
    end

    bestFinalVals = [];
    bestCoarseVals = [];
    for t = 1:numel(bh)
        pos = nested_get(bh(t), {'position'}, []);
        if isempty(pos)
            continue;
        end
        xAgg.X_acceptedCount_total = xAgg.X_acceptedCount_total + safe_num(nested_get(pos, {'acceptedCount'}, 0));

        wg = nested_get(pos, {'waveguideInfo'}, []);
        if isempty(wg)
            continue;
        end
        for n = 1:numel(wg)
            w = wg(n);
            xAgg.X_totalLineSearchSteps = xAgg.X_totalLineSearchSteps + safe_num(nested_get(w, {'lineSearchSteps'}, 0));
            xAgg.X_projectionCollapsedCount_total = xAgg.X_projectionCollapsedCount_total + safe_num(nested_get(w, {'projectionCollapsedCount'}, 0));
            xAgg.X_activeBoundaryCount_total = xAgg.X_activeBoundaryCount_total + double(logical(nested_get(w, {'activeBoundaryConstraint'}, false)));
            xAgg.X_activeSpacingCount_total = xAgg.X_activeSpacingCount_total + double(logical(nested_get(w, {'activeSpacingConstraint'}, false)));

            bfi = nested_get(w, {'bestFinalRate'}, NaN) - nested_get(w, {'sumRateBefore'}, NaN);
            if isfinite(bfi)
                bestFinalVals(end + 1, 1) = bfi; %#ok<AGROW>
            end
            bci = nested_get(w, {'bestCoarseRate'}, NaN) - nested_get(w, {'sumRateBefore'}, NaN);
            if isfinite(bci)
                bestCoarseVals(end + 1, 1) = bci; %#ok<AGROW>
            end

            bci2 = nested_get(w, {'bestCoarseImprove'}, NaN);
            if isfinite(bci2)
                bestCoarseVals(end + 1, 1) = bci2; %#ok<AGROW>
            end
            bfi2 = nested_get(w, {'bestFinalImprove'}, NaN);
            if isfinite(bfi2)
                bestFinalVals(end + 1, 1) = bfi2; %#ok<AGROW>
            end
        end
    end

    if ~isempty(bestFinalVals)
        xAgg.X_bestFinalImprove_max = max(bestFinalVals);
        xAgg.X_bestFinalImprove_mean = mean(bestFinalVals);
    end
    if ~isempty(bestCoarseVals)
        xAgg.X_bestCoarseImprove_max = max(bestCoarseVals);
    end
end

function rankingTable = build_ranking_table(allCases, cfg)
    n = numel(allCases);
    caseId = strings(n,1);
    seed = nan(n,1);
    symbolSeed = nan(n,1);
    presetName = strings(n,1);
    userMode = strings(n,1);
    finalSumRate = nan(n,1);
    converged = false(n,1);
    iterations = nan(n,1);
    terminationReason = strings(n,1);
    constraintsAllSatisfied = false(n,1);
    runtimeTotal = nan(n,1);
    S_triggered_total = zeros(n,1);
    S_acceptedSwaps_total = zeros(n,1);
    S_bestDeltaFinal_max = nan(n,1);
    S_bestDeltaFinal_mean = nan(n,1);
    S_positiveFinal_total = zeros(n,1);
    S_aboveEpsilonS_total = zeros(n,1);
    S_postAcceptGain_sum = zeros(n,1);
    S_intensificationTriggered_total = zeros(n,1);
    X_acceptedCount_total = zeros(n,1);
    X_bestFinalImprove_max = nan(n,1);
    X_bestFinalImprove_mean = nan(n,1);
    X_totalLineSearchSteps = zeros(n,1);
    X_bestCoarseImprove_max = nan(n,1);
    success = false(n,1);
    errorMessage = strings(n,1);

    for i = 1:n
        c = allCases(i);
        caseId(i) = string(c.caseId);
        seed(i) = safe_num(nested_get(c.meta, {'seed'}, NaN));
        symbolSeed(i) = safe_num(nested_get(c.meta, {'symbolSeed'}, NaN));
        presetName(i) = string(nested_get(c.meta, {'presetName'}, 'NA'));
        userMode(i) = string(nested_get(c.meta, {'userMode'}, 'NA'));
        success(i) = logical(c.success);
        errorMessage(i) = string(c.errorMessage);

        if c.success
            f = c.features;
            finalSumRate(i) = safe_num(nested_get(f, {'finalSumRate'}, NaN));
            converged(i) = logical(nested_get(f, {'converged'}, false));
            iterations(i) = safe_num(nested_get(f, {'iterations'}, NaN));
            terminationReason(i) = string(nested_get(f, {'terminationReason'}, 'NA'));
            constraintsAllSatisfied(i) = logical(nested_get(f, {'constraintsAllSatisfied'}, false));
            runtimeTotal(i) = safe_num(nested_get(f, {'runtimeTotal'}, NaN));

            S_triggered_total(i) = safe_num(nested_get(f, {'S_triggered_total'}, 0));
            S_acceptedSwaps_total(i) = safe_num(nested_get(f, {'S_acceptedSwaps_total'}, 0));
            S_bestDeltaFinal_max(i) = safe_num(nested_get(f, {'S_bestDeltaFinal_max'}, NaN));
            S_bestDeltaFinal_mean(i) = safe_num(nested_get(f, {'S_bestDeltaFinal_mean'}, NaN));
            S_positiveFinal_total(i) = safe_num(nested_get(f, {'S_positiveFinal_total'}, 0));
            S_aboveEpsilonS_total(i) = safe_num(nested_get(f, {'S_aboveEpsilonS_total'}, 0));
            S_postAcceptGain_sum(i) = safe_num(nested_get(f, {'S_postAcceptGain_sum'}, 0));
            S_intensificationTriggered_total(i) = safe_num(nested_get(f, {'S_intensificationTriggered_total'}, 0));

            X_acceptedCount_total(i) = safe_num(nested_get(f, {'X_acceptedCount_total'}, 0));
            X_bestFinalImprove_max(i) = safe_num(nested_get(f, {'X_bestFinalImprove_max'}, NaN));
            X_bestFinalImprove_mean(i) = safe_num(nested_get(f, {'X_bestFinalImprove_mean'}, NaN));
            X_totalLineSearchSteps(i) = safe_num(nested_get(f, {'X_totalLineSearchSteps'}, 0));
            X_bestCoarseImprove_max(i) = safe_num(nested_get(f, {'X_bestCoarseImprove_max'}, NaN));
        else
            terminationReason(i) = "FAILED";
        end
    end

    T = table(caseId, seed, symbolSeed, presetName, userMode, ...
        finalSumRate, converged, iterations, terminationReason, constraintsAllSatisfied, runtimeTotal, ...
        S_triggered_total, S_acceptedSwaps_total, S_bestDeltaFinal_max, S_bestDeltaFinal_mean, ...
        S_positiveFinal_total, S_aboveEpsilonS_total, S_postAcceptGain_sum, S_intensificationTriggered_total, ...
        X_acceptedCount_total, X_bestFinalImprove_max, X_bestFinalImprove_mean, X_totalLineSearchSteps, X_bestCoarseImprove_max, ...
        success, errorMessage);

    T.screenScore = compute_screen_score(T, cfg);

    [~, ord] = sort(T.screenScore, 'descend');
    rankingTable = T(ord, :);
    rankingTable.rank = (1:height(rankingTable)).';
    rankingTable = movevars(rankingTable, 'rank', 'Before', 1);

    % 保留问题要求的核心列优先展示
    keyCols = {'rank','caseId','seed','symbolSeed','presetName','userMode', ...
        'finalSumRate','converged','iterations','constraintsAllSatisfied', ...
        'S_acceptedSwaps_total','S_bestDeltaFinal_max','S_positiveFinal_total', ...
        'X_acceptedCount_total','X_bestFinalImprove_max','X_totalLineSearchSteps','screenScore'};
    otherCols = setdiff(rankingTable.Properties.VariableNames, keyCols, 'stable');
    rankingTable = rankingTable(:, [keyCols, otherCols]);
end

function score = compute_screen_score(T, cfg)
    n = height(T);
    score = zeros(n, 1);

    % 正向特征（min-max 归一化）
    zSswap = normalize_feature_column(T.S_acceptedSwaps_total);
    zSbest = normalize_feature_column(max(0, T.S_bestDeltaFinal_max));
    zXacc  = normalize_feature_column(T.X_acceptedCount_total);
    zXbest = normalize_feature_column(max(0, T.X_bestFinalImprove_max));
    zSpos  = normalize_feature_column(T.S_positiveFinal_total);
    zXls   = normalize_feature_column(T.X_totalLineSearchSteps);
    zRate  = normalize_feature_column(T.finalSumRate);

    w = cfg.scoreWeights;
    score = score + w.S_acceptedSwaps_total   * zSswap;
    score = score + w.S_bestDeltaFinal_max_pos * zSbest;
    score = score + w.X_acceptedCount_total   * zXacc;
    score = score + w.X_bestFinalImprove_max_pos * zXbest;
    score = score + w.S_positiveFinal_total   * zSpos;
    score = score + w.X_totalLineSearchSteps  * zXls;
    score = score + w.finalSumRate            * zRate;

    % 约束/收敛惩罚
    badConstraint = ~T.constraintsAllSatisfied;
    score(badConstraint) = score(badConstraint) - cfg.penaltyConstraintViolation;

    notConv = ~T.converged;
    score(notConv) = score(notConv) - cfg.penaltyNotConverged;

    % 失败或异常 case：直接排后
    invalid = ~T.success | ~isfinite(T.finalSumRate);
    score(invalid) = cfg.invalidCaseScore;
end

function z = normalize_feature_column(x)
    x = double(x);
    z = zeros(size(x));
    valid = isfinite(x);
    if ~any(valid)
        return;
    end
    xmin = min(x(valid));
    xmax = max(x(valid));
    if xmax <= xmin + eps
        z(valid) = 1;
        return;
    end
    z(valid) = (x(valid) - xmin) ./ (xmax - xmin);
end

function allCases = attach_scores_to_cases(allCases, rankingTable)
    scoreMap = containers.Map('KeyType', 'char', 'ValueType', 'double');
    for i = 1:height(rankingTable)
        scoreMap(char(rankingTable.caseId(i))) = rankingTable.screenScore(i);
    end
    for i = 1:numel(allCases)
        k = allCases(i).caseId;
        if isKey(scoreMap, k)
            allCases(i).screenScore = scoreMap(k);
        end
    end
end

function print_screen_summary(rankingTable, allCases, cfg)
    n = numel(allCases);
    ok = sum([allCases.success]);
    fail = n - ok;

    fprintf('\n===== screen_sx_active_cases: summary =====\n');
    fprintf('Total cases: %d | Success: %d | Failed: %d\n', n, ok, fail);

    if fail > 0
        reasons = strings(fail,1);
        idx = 0;
        for i = 1:numel(allCases)
            if ~allCases(i).success
                idx = idx + 1;
                reasons(idx) = string(allCases(i).errorMessage);
            end
        end
        [u,~,ic] = unique(reasons);
        cnt = accumarray(ic, 1);
        fprintf('Top failure reasons:\n');
        [~, ord] = sort(cnt, 'descend');
        for k = 1:min(5, numel(ord))
            fprintf('  - x%d: %s\n', cnt(ord(k)), u(ord(k)));
        end
    end

    if isempty(rankingTable)
        fprintf('No ranking output available.\n');
        return;
    end

    topK = min(cfg.topKPrint, height(rankingTable));
    disp(rankingTable(1:topK, {'rank','caseId','presetName','seed','finalSumRate', ...
        'S_acceptedSwaps_total','X_acceptedCount_total','S_bestDeltaFinal_max', ...
        'X_bestFinalImprove_max','screenScore'}));

    best = rankingTable(1,:);
    fprintf('Recommended case: %s | preset=%s | seed=%d | score=%.4f\n', ...
        best.caseId, best.presetName, best.seed, best.screenScore);
end

function save_screen_outputs(rankingTable, allCases, cfg)
    if ~exist(cfg.saveDir, 'dir')
        mkdir(cfg.saveDir);
    end
    save(fullfile(cfg.saveDir, 'screen_results.mat'), 'rankingTable', 'allCases', 'cfg');

    topK = min(cfg.saveTopKCases, height(rankingTable));
    for k = 1:topK
        cid = char(rankingTable.caseId(k));
        idx = find(strcmp({allCases.caseId}, cid), 1, 'first');
        if isempty(idx) || ~allCases(idx).success
            continue;
        end
        oneCase = allCases(idx); %#ok<NASGU>
        safeName = regexprep(cid, '[^a-zA-Z0-9_\-]', '_');
        save(fullfile(cfg.saveDir, sprintf('top_case_%02d_%s.mat', k, safeName)), 'oneCase');
    end
end

function reason = get_last_w_stop_reason(aoInfo)
    reason = 'noIteration';
    bh = nested_get(aoInfo, {'blockHistory'}, []);
    if ~isempty(bh)
        reason = string(nested_get(bh(end), {'W','stopReason'}, 'NA'));
    end
end

function reason = get_last_s_break_reason(aoInfo)
    reason = 'noIteration';
    bh = nested_get(aoInfo, {'blockHistory'}, []);
    if ~isempty(bh)
        reason = string(nested_get(bh(end), {'userSet','breakReason'}, 'NA'));
    end
end

function ng = collect_numerical_guard_info(aoInfo)
    ng = struct('wNumericalGuardTriggered', false);
    bh = nested_get(aoInfo, {'blockHistory'}, []);
    if ~isempty(bh)
        ng.wNumericalGuardTriggered = logical(nested_get(bh(end), {'W','numericalGuardTriggered'}, false));
    end
end

function out = merge_struct(base, override)
    out = base;
    if isempty(override) || ~isstruct(override)
        return;
    end
    fn = fieldnames(override);
    for i = 1:numel(fn)
        f = fn{i};
        if isfield(out, f) && isstruct(out.(f)) && isstruct(override.(f))
            out.(f) = merge_struct(out.(f), override.(f));
        else
            out.(f) = override.(f);
        end
    end
end

function v = nested_get(s, fields, defaultVal)
    v = defaultVal;
    try
        cur = s;
        for i = 1:numel(fields)
            f = fields{i};
            if isstruct(cur) && isfield(cur, f)
                cur = cur.(f);
            else
                return;
            end
        end
        v = cur;
    catch
        v = defaultVal;
    end
end

function x = safe_num(v)
    if isempty(v) || ~isnumeric(v) || ~isscalar(v) || ~isfinite(v)
        x = 0;
    else
        x = double(v);
    end
end
