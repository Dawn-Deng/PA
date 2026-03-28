function figs = plot_all_results(resultsOrMatFile, saveDir)
%PLOT_ALL_RESULTS 集中式结果可视化入口（不改动原算法流程，可选调用）。
%
% 用法1：工作区已有 results
% figs = plot_all_results(results);
%
% 用法2：从 mat 文件读取
% figs = plot_all_results('results_latest.mat');
%
% 用法3：保存到指定目录
% figs = plot_all_results(results, 'figures');
%
% 用法4：自动从 base workspace 读取 results
% figs = plot_all_results();

    if nargin < 1
        resultsOrMatFile = [];
    end
    if nargin < 2
        saveDir = '';
    end

    results = parse_results_input(resultsOrMatFile);
    doSave = ~isempty(saveDir);
    if doSave && ~exist(saveDir, 'dir')
        mkdir(saveDir);
    end

    figs = struct();
    apply_default_plot_style();

    figs.sum_rate_convergence = safe_plot(@() plot_convergence_figure(results), 'sum_rate_convergence');
    save_figure_if_needed(figs.sum_rate_convergence, saveDir, 'sum_rate_convergence', doSave);

    figs.block_rate_gain = safe_plot(@() plot_block_gain_figure(results), 'block_rate_gain');
    save_figure_if_needed(figs.block_rate_gain, saveDir, 'block_rate_gain', doSave);

    initFigs = safe_plot(@() plot_init_figures(results), 'init_figures');
    if isstruct(initFigs)
        figs.init_emax = get_struct_field(initFigs, 'init_emax', []);
        figs.init_utility_matrix = get_struct_field(initFigs, 'init_utility_matrix', []);
        figs.init_movement_cost = get_struct_field(initFigs, 'init_movement_cost', []);
        figs.init_x_compare = get_struct_field(initFigs, 'init_x_compare', []);
        save_figure_if_needed(figs.init_emax, saveDir, 'init_emax', doSave);
        save_figure_if_needed(figs.init_utility_matrix, saveDir, 'init_utility_matrix', doSave);
        save_figure_if_needed(figs.init_movement_cost, saveDir, 'init_movement_cost', doSave);
        save_figure_if_needed(figs.init_x_compare, saveDir, 'init_x_compare', doSave);
    end

    signalFigs = safe_plot(@() plot_signal_figures(results), 'signal_figures');
    if isstruct(signalFigs)
        figs.final_rate = get_struct_field(signalFigs, 'final_rate', []);
        figs.final_sinr = get_struct_field(signalFigs, 'final_sinr', []);
        save_figure_if_needed(figs.final_rate, saveDir, 'final_rate', doSave);
        save_figure_if_needed(figs.final_sinr, saveDir, 'final_sinr', doSave);
    end

    figs.geometry_snapshots = safe_plot(@() plot_geometry_figures(results), 'geometry_snapshots');
    save_figure_if_needed(figs.geometry_snapshots, saveDir, 'geometry_snapshots', doSave);

    runtimeFigs = safe_plot(@() plot_runtime_figures(results), 'runtime_figures');
    if isstruct(runtimeFigs)
        figs.runtime_totals = get_struct_field(runtimeFigs, 'runtime_totals', []);
        figs.runtime_per_iteration = get_struct_field(runtimeFigs, 'runtime_per_iteration', []);
        save_figure_if_needed(figs.runtime_totals, saveDir, 'runtime_totals', doSave);
        save_figure_if_needed(figs.runtime_per_iteration, saveDir, 'runtime_per_iteration', doSave);
    end

    figs.wmmse_inner_last_iter = safe_plot(@() plot_wmmse_inner_figures(results), 'wmmse_inner_last_iter');
    save_figure_if_needed(figs.wmmse_inner_last_iter, saveDir, 'wmmse_inner_last_iter', doSave);

    figs.user_set_evolution = safe_plot(@() plot_user_set_figures(results), 'user_set_evolution');
    save_figure_if_needed(figs.user_set_evolution, saveDir, 'user_set_evolution', doSave);

    figs.x_position_evolution = safe_plot(@() plot_position_update_figures(results), 'x_position_evolution');
    save_figure_if_needed(figs.x_position_evolution, saveDir, 'x_position_evolution', doSave);
end

function out = safe_plot(plotFcn, figTag)
    out = [];
    try
        out = plotFcn();
    catch ME
        warning('plot_all_results:%sFailed', figTag, 'Failed to generate %s: %s', figTag, ME.message);
    end
end

function results = parse_results_input(resultsOrMatFile)
    if nargin == 0 || isempty(resultsOrMatFile)
        if evalin('base', 'exist(''results'', ''var'')')
            results = evalin('base', 'results');
        else
            error('plot_all_results:NoInput', 'No input provided and base workspace variable ''results'' not found.');
        end
        return;
    end

    if isstruct(resultsOrMatFile)
        results = resultsOrMatFile;
        return;
    end

    if ischar(resultsOrMatFile) || isstring(resultsOrMatFile)
        matPath = char(resultsOrMatFile);
        if ~exist(matPath, 'file')
            error('plot_all_results:MatNotFound', 'MAT file not found: %s', matPath);
        end
        S = load(matPath);
        if isfield(S, 'results')
            results = S.results;
        else
            names = fieldnames(S);
            idx = find(structfun(@isstruct, S), 1, 'first');
            if isempty(idx)
                error('plot_all_results:MatNoResults', 'No struct variable found in MAT file: %s', matPath);
            end
            warning('plot_all_results:FallbackStruct', ...
                'Variable ''results'' not found. Using struct variable ''%s''.', names{idx});
            results = S.(names{idx});
        end
        return;
    end

    error('plot_all_results:InvalidInput', ...
        'Input must be a struct, a MAT file path, or empty.');
end

function apply_default_plot_style()
    set(groot, 'defaultAxesFontSize', 11);
    set(groot, 'defaultAxesLineWidth', 1.0);
    set(groot, 'defaultLineLineWidth', 1.6);
    set(groot, 'defaultLineMarkerSize', 6);
    set(groot, 'defaultFigureColor', 'w');
end

function fig = plot_convergence_figure(results)
    hist = get_nested(results, {'trace','outerSumRateHistory'}, []);
    if isempty(hist)
        warning('plot_all_results:NoOuterHistory', 'Missing results.trace.outerSumRateHistory. Skip convergence figure.');
        fig = [];
        return;
    end
    hist = hist(:);
    fig = figure('Name', 'Sum Rate Convergence');
    plot(0:numel(hist)-1, hist, '-o');
    title('Sum Rate Convergence');
    xlabel('Outer iteration');
    ylabel('Sum rate');
    grid on; box on;
end

function fig = plot_block_gain_figure(results)
    it = get_nested(results, {'trace','iterationTrace'}, []);
    if isempty(it)
        warning('plot_all_results:NoIterationTrace', 'Missing results.trace.iterationTrace. Skip block gain figure.');
        fig = [];
        return;
    end

    [dW, dA, dX, dS] = extract_block_deltas(it);
    if isempty(dW)
        warning('plot_all_results:NoBlockDelta', 'Missing deltaW/deltaAngle/deltaX/deltaS. Skip block gain figure.');
        fig = [];
        return;
    end

    Y = [dW(:), dA(:), dX(:), dS(:)];
    fig = figure('Name', 'Per-Block Rate Gain');
    b = bar(Y, 'stacked'); %#ok<NASGU>
    title('Per-Block Rate Gain');
    xlabel('Outer iteration');
    ylabel('Rate gain');
    legend({'W','Angle','X','S'}, 'Location', 'best');
    grid on; box on;
end

function figs = plot_init_figures(results)
    figs = struct('init_emax', [], 'init_utility_matrix', [], 'init_movement_cost', [], 'init_x_compare', []);
    initInfo = get_nested(results, {'initInfo'}, []);
    if isempty(initInfo)
        warning('plot_all_results:NoInitInfo', 'Missing results.initInfo. Skip initialization figures.');
        return;
    end

    % Emax 排序图
    Emax = get_nested(initInfo, {'Emax'}, []);
    if ~isempty(Emax)
        Emax = Emax(:);
        [sortedVal, ord] = sort(Emax, 'descend');
        figs.init_emax = figure('Name', 'Initialization Emax Ranking');
        plot(1:numel(sortedVal), sortedVal, '-o');
        hold on;
        candidatePool = get_nested(initInfo, {'candidatePool'}, []);
        if ~isempty(candidatePool)
            inPool = ismember(ord, candidatePool(:));
            scatter(find(inPool), sortedVal(inPool), 50, 'r', 'filled');
            legend({'Emax sorted', 'Candidate pool'}, 'Location', 'best');
        else
            legend({'Emax sorted'}, 'Location', 'best');
        end
        hold off;
        title('Initialization Emax Ranking');
        xlabel('Ranked user index');
        ylabel('Emax');
        grid on; box on;
    else
        warning('plot_all_results:NoEmax', 'Missing results.initInfo.Emax.');
    end

    % utilityMatrix 热力图
    U = get_nested(initInfo, {'utilityMatrix'}, []);
    if ~isempty(U)
        figs.init_utility_matrix = figure('Name', 'Initialization Utility Matrix');
        imagesc(U);
        axis tight;
        colorbar;
        title('Initialization Utility Matrix');
        xlabel('PA index (flattened n,m)');
        ylabel('Candidate user row');
        grid on; box on;
        hold on;
        selectedPairs = get_nested(initInfo, {'matching','selectedPairs'}, []);
        candidatePool = get_nested(initInfo, {'candidatePool'}, []);
        M = get_nested(results, {'params','M'}, []);
        if ~isempty(selectedPairs) && ~isempty(candidatePool) && ~isempty(M)
            for i = 1:numel(selectedPairs)
                sp = selectedPairs(i);
                if ~isfield(sp, 'userIndex') || ~isfield(sp, 'waveguideIndex') || ~isfield(sp, 'paIndex')
                    continue;
                end
                row = find(candidatePool(:) == sp.userIndex, 1, 'first');
                col = (sp.waveguideIndex - 1) * M + sp.paIndex;
                if ~isempty(row) && isfinite(col) && col >= 1 && col <= size(U,2)
                    plot(col, row, 'ws', 'LineWidth', 1.5, 'MarkerSize', 8);
                end
            end
        end
        hold off;
    else
        warning('plot_all_results:NoUtility', 'Missing results.initInfo.utilityMatrix.');
    end

    % movementCost 统计图（可选）
    moveCost = get_nested(initInfo, {'movementCost'}, []);
    if ~isempty(moveCost)
        figs.init_movement_cost = figure('Name', 'Initialization Movement Cost');
        histogram(moveCost(:), 20);
        title('Initialization Movement Cost Distribution');
        xlabel('Movement cost');
        ylabel('Count');
        grid on; box on;
    end

    % Xbar0 vs X0（可选）
    Xbar0 = get_nested(initInfo, {'Xbar0'}, []);
    X0 = get_nested(initInfo, {'X0'}, []);
    if ~isempty(Xbar0) && ~isempty(X0)
        figs.init_x_compare = figure('Name', 'Initialization Xbar0 vs X0');
        plot(Xbar0(:), '-o', 'DisplayName', 'Xbar0'); hold on;
        plot(X0(:), '-s', 'DisplayName', 'X0'); hold off;
        title('Initialization Xbar0 vs X0');
        xlabel('PA index (flattened)');
        ylabel('Position Y');
        legend('Location', 'best');
        grid on; box on;
    end
end

function figs = plot_signal_figures(results)
    figs = struct('final_rate', [], 'final_sinr', []);
    signalInfo = get_nested(results, {'signalInfo'}, []);
    if isempty(signalInfo)
        warning('plot_all_results:NoSignalInfo', 'Missing results.signalInfo. Skip signal figures.');
        return;
    end

    S = get_nested(signalInfo, {'S'}, []);
    rates = get_nested(signalInfo, {'rate'}, []);
    sinrVal = get_nested(signalInfo, {'sinr'}, []);

    % user axis 优先使用服务用户编号
    if ~isempty(S)
        xAxis = S(:);
    elseif ~isempty(rates)
        xAxis = (1:numel(rates)).';
    elseif ~isempty(sinrVal)
        xAxis = (1:numel(sinrVal)).';
    else
        xAxis = [];
    end

    if ~isempty(rates)
        rates = rates(:);
        if ~isempty(xAxis) && numel(xAxis) == numel(rates)
            figs.final_rate = figure('Name', 'Final User Rate');
            bar(xAxis, rates);
            xlabel('User index');
        else
            figs.final_rate = figure('Name', 'Final User Rate');
            bar(rates);
            xlabel('Scheduled user order');
        end
        ylabel('Rate (bit/s/Hz)');
        title('Final User Rate');
        grid on; box on;
    else
        warning('plot_all_results:NoRate', 'Missing results.signalInfo.rate.');
    end

    if ~isempty(sinrVal)
        sinrVal = sinrVal(:);
        if ~isempty(xAxis) && numel(xAxis) == numel(sinrVal)
            figs.final_sinr = figure('Name', 'Final User SINR');
            stem(xAxis, sinrVal, 'filled');
            xlabel('User index');
        else
            figs.final_sinr = figure('Name', 'Final User SINR');
            stem(sinrVal, 'filled');
            xlabel('Scheduled user order');
        end
        ylabel('SINR (linear)');
        title('Final User SINR');
        grid on; box on;
    else
        warning('plot_all_results:NoSINR', 'Missing results.signalInfo.sinr.');
    end
end

function fig = plot_geometry_figures(results)
    params = get_nested(results, {'params'}, []);
    snaps = get_nested(results, {'snapshots'}, []);
    if isempty(params) || isempty(snaps)
        warning('plot_all_results:NoGeometryInfo', 'Missing results.params or results.snapshots. Skip geometry figure.');
        fig = [];
        return;
    end

    N = get_nested(params, {'N'}, []);
    M = get_nested(params, {'M'}, []);
    Dx = get_nested(params, {'Dx'}, []);
    d = get_nested(params, {'d'}, []);
    if isempty(N) || isempty(M) || isempty(Dx) || isempty(d)
        warning('plot_all_results:NoGeometryParams', 'Need params.N/M/Dx/d for geometry restoration.');
        fig = [];
        return;
    end

    users = get_nested(snaps, {'initial','users'}, []);
    if isempty(users)
        users = get_nested(results, {'channelInfo','users'}, []);
    end

    xW = zeros(N, 1);
    for n = 1:N
        xW(n) = ((2*n - 1) / (2*N)) * Dx;
    end

    phaseNames = {'initial','afterInit','final'};
    phaseTitles = {'Initial', 'After Initialization', 'Final'};

    fig = figure('Name', 'Geometry Snapshots');
    tiledlayout(1,3, 'Padding', 'compact', 'TileSpacing', 'compact');
    for p = 1:3
        nexttile;
        snap = get_nested(snaps, {phaseNames{p}}, []);
        if isempty(snap)
            title(sprintf('%s (missing)', phaseTitles{p}));
            axis off;
            continue;
        end

        X = get_nested(snap, {'X'}, []);
        theta = get_nested(snap, {'theta'}, []);
        phi = get_nested(snap, {'phi'}, []);
        S = get_nested(snap, {'S'}, get_nested(results, {'signalInfo','S'}, []));

        if isempty(X)
            title(sprintf('%s (X missing)', phaseTitles{p}));
            axis off;
            continue;
        end

        % 用户散点
        if ~isempty(users) && size(users,2) >= 2
            scatter(users(:,1), users(:,2), 20, [0.65 0.65 0.65], 'filled', 'DisplayName', 'Users');
            hold on;
            if ~isempty(S)
                validS = S(:);
                validS = validS(validS>=1 & validS<=size(users,1));
                if ~isempty(validS)
                    scatter(users(validS,1), users(validS,2), 35, 'r', 'filled', 'DisplayName', 'Served users');
                end
            end
        else
            hold on;
        end

        % 波导位置
        scatter(xW, zeros(size(xW)), 36, 'k', '^', 'filled', 'DisplayName', 'Waveguides');

        % PA 位置 pa=[xW; X(m,n); d] -> 画平面 (xW, X)
        paX = repmat(xW, 1, M); % N x M
        paY = X.';              % N x M
        scatter(paX(:), paY(:), 22, [0 0.45 0.74], 'filled', 'DisplayName', 'PAs');

        % 朝向箭头（可选）
        if ~isempty(theta) && ~isempty(phi) && isequal(size(theta), size(X)) && isequal(size(phi), size(X))
            thetaNM = theta.'; % N x M
            phiNM = phi.';     % N x M
            dirX = sin(thetaNM) .* cos(phiNM);
            dirY = sin(thetaNM) .* sin(phiNM);

            px = paX(:);
            py = paY(:);
            dx = dirX(:);
            dy = dirY(:);
            projNorm = hypot(dx, dy);
            mask = projNorm > 1e-3;

            if any(mask)
                dxn = dx(mask) ./ projNorm(mask);
                dyn = dy(mask) ./ projNorm(mask);
                arrowScale = 0.18;
                quiver(px(mask), py(mask), arrowScale * dxn, arrowScale * dyn, 0, ...
                    'Color', [0.2 0.2 0.2], 'LineWidth', 0.9, 'MaxHeadSize', 0.45, ...
                    'AutoScale', 'off', 'DisplayName', 'Orientation');
            end
        end

        title(phaseTitles{p});
        xlabel('x');
        ylabel('y');
        grid on; box on;
        axis tight;
        hold off;
    end
end

function figs = plot_runtime_figures(results)
    figs = struct('runtime_totals', [], 'runtime_per_iteration', []);
    rt = get_nested(results, {'runtime'}, []);
    if isempty(rt)
        warning('plot_all_results:NoRuntime', 'Missing results.runtime. Skip runtime figures.');
        return;
    end

    keys = {'totalWTime','totalAngleTime','totalXTime','totalSTime'};
    vals = nan(1,4);
    for i = 1:4
        vals(i) = get_nested(rt, {keys{i}}, NaN);
    end
    if all(isfinite(vals))
        figs.runtime_totals = figure('Name', 'Runtime Breakdown');
        bar(vals);
        set(gca, 'XTickLabel', {'W','Angle','X','S'});
        ylabel('Time (s)');
        title('Runtime Breakdown');
        grid on; box on;
    else
        warning('plot_all_results:NoRuntimeTotals', 'Missing one or more runtime total fields.');
    end

    iterRt = get_nested(rt, {'iteration'}, []);
    if ~isempty(iterRt)
        Wt = extract_struct_vec(iterRt, 'WTime');
        At = extract_struct_vec(iterRt, 'angleTime');
        Xt = extract_struct_vec(iterRt, 'XTime');
        St = extract_struct_vec(iterRt, 'STime');
        L = min([numel(Wt), numel(At), numel(Xt), numel(St)]);
        if L > 0
            figs.runtime_per_iteration = figure('Name', 'Runtime Per Iteration');
            plot(1:L, Wt(1:L), '-o', 'DisplayName', 'WTime'); hold on;
            plot(1:L, At(1:L), '-s', 'DisplayName', 'angleTime');
            plot(1:L, Xt(1:L), '-d', 'DisplayName', 'XTime');
            plot(1:L, St(1:L), '-^', 'DisplayName', 'STime'); hold off;
            title('Runtime Per Iteration');
            xlabel('Outer iteration');
            ylabel('Time (s)');
            legend('Location', 'best');
            grid on; box on;
        end
    end
end

function fig = plot_wmmse_inner_figures(results)
    fig = [];
    blockHist = get_nested(results, {'trace','blockHistory'}, []);
    if isempty(blockHist)
        warning('plot_all_results:NoBlockHistory', 'Missing results.trace.blockHistory.');
        return;
    end

    idx = [];
    for i = numel(blockHist):-1:1
        if isfield(blockHist(i), 'W') && ~isempty(blockHist(i).W)
            idx = i;
            break;
        end
    end
    if isempty(idx)
        warning('plot_all_results:NoWBlock', 'No W block found in blockHistory.');
        return;
    end

    Wblk = blockHist(idx).W;
    s1 = get_nested(Wblk, {'sumRateHistory'}, []);
    s2 = get_nested(Wblk, {'deltaSumRateHistory'}, []);
    s3 = get_nested(Wblk, {'transmitPowerHistory'}, []);
    s4 = get_nested(Wblk, {'muHistory'}, []);
    if isempty(s1) && isempty(s2) && isempty(s3) && isempty(s4)
        warning('plot_all_results:NoWInnerTrace', 'No W-inner trace fields available.');
        return;
    end

    fig = figure('Name', 'WMMSE Inner Convergence (Last W-Block)');
    tiledlayout(2,2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile;
    if ~isempty(s1), plot(s1, '-o'); title('sumRateHistory'); else, text(0.1,0.5,'missing'); axis off; end
    grid on; box on; xlabel('Inner iter'); ylabel('Sum rate');

    nexttile;
    if ~isempty(s2), plot(s2, '-o'); title('deltaSumRateHistory'); else, text(0.1,0.5,'missing'); axis off; end
    grid on; box on; xlabel('Inner iter'); ylabel('\Delta sum rate');

    nexttile;
    if ~isempty(s3), plot(s3, '-o'); title('transmitPowerHistory'); else, text(0.1,0.5,'missing'); axis off; end
    grid on; box on; xlabel('Inner iter'); ylabel('Power');

    nexttile;
    if ~isempty(s4), plot(s4, '-o'); title('muHistory'); else, text(0.1,0.5,'missing'); axis off; end
    grid on; box on; xlabel('Inner iter'); ylabel('\mu');
end

function fig = plot_user_set_figures(results)
    fig = [];
    iterTrace = get_nested(results, {'trace','iterationTrace'}, []);
    if isempty(iterTrace)
        warning('plot_all_results:NoIterationTraceUserSet', 'Missing iterationTrace for user-set evolution.');
        return;
    end

    allS = cell(numel(iterTrace),1);
    for t = 1:numel(iterTrace)
        s = get_nested(iterTrace(t), {'snapshotS'}, []);
        if isempty(s)
            s = [];
        end
        allS{t} = s(:).';
    end

    K = get_nested(results, {'params','K'}, 0);
    if K <= 0
        K = max([cellfun(@(x) max([x,0]), allS)]);
    end
    if K <= 0
        warning('plot_all_results:NoUserSetData', 'Cannot infer user index range for user-set evolution.');
        return;
    end

    M = zeros(numel(allS), K);
    for t = 1:numel(allS)
        s = allS{t};
        s = s(s>=1 & s<=K);
        M(t, s) = 1;
    end

    fig = figure('Name', 'User Set Evolution');
    imagesc(1:K, 1:size(M,1), M);
    colormap(gca, [1 1 1; 0.1 0.55 0.1]);
    colorbar('Ticks', [0,1], 'TickLabels', {'Not served','Served'});
    title('User Set Evolution');
    xlabel('User index');
    ylabel('Outer iteration');
    grid on; box on;
end

function fig = plot_position_update_figures(results)
    fig = [];
    iterTrace = get_nested(results, {'trace','iterationTrace'}, []);
    if isempty(iterTrace)
        warning('plot_all_results:NoIterationTraceX', 'Missing iterationTrace for X evolution.');
        return;
    end

    xCells = {};
    for t = 1:numel(iterTrace)
        Xt = get_nested(iterTrace(t), {'snapshotX'}, []);
        if ~isempty(Xt)
            xCells{end+1} = Xt(:).'; %#ok<AGROW>
        end
    end
    if isempty(xCells)
        warning('plot_all_results:NoSnapshotX', 'No snapshotX found in iterationTrace.');
        return;
    end

    L = min(cellfun(@numel, xCells));
    Xmat = zeros(numel(xCells), L);
    for i = 1:numel(xCells)
        Xmat(i,:) = xCells{i}(1:L);
    end

    fig = figure('Name', 'X Position Evolution');
    imagesc(1:L, 1:size(Xmat,1), Xmat);
    colorbar;
    title('X Position Evolution');
    xlabel('PA index (flattened)');
    ylabel('Outer iteration');
    grid on; box on;
end

function save_figure_if_needed(fig, saveDir, fileStem, doSave)
    if ~doSave || isempty(fig) || ~ishandle(fig)
        return;
    end
    try
        pngPath = fullfile(saveDir, [fileStem, '.png']);
        figPath = fullfile(saveDir, [fileStem, '.fig']);
        saveas(fig, pngPath);
        savefig(fig, figPath);
    catch ME
        warning('plot_all_results:SaveFailed', 'Failed to save figure %s: %s', fileStem, ME.message);
    end
end

function v = get_nested(s, fields, defaultVal)
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

function [dW, dA, dX, dS] = extract_block_deltas(iterationTrace)
    dW = extract_struct_vec(iterationTrace, 'deltaW');
    dA = extract_struct_vec(iterationTrace, 'deltaAngle');
    dX = extract_struct_vec(iterationTrace, 'deltaX');
    dS = extract_struct_vec(iterationTrace, 'deltaS');
    L = min([numel(dW), numel(dA), numel(dX), numel(dS)]);
    if L == 0
        dW = [];
        dA = [];
        dX = [];
        dS = [];
        return;
    end
    dW = dW(1:L);
    dA = dA(1:L);
    dX = dX(1:L);
    dS = dS(1:L);
end

function vec = extract_struct_vec(st, fieldName)
    vec = [];
    if isempty(st)
        return;
    end
    vec = nan(numel(st), 1);
    for i = 1:numel(st)
        if isfield(st(i), fieldName)
            x = st(i).(fieldName);
            if isempty(x)
                vec(i) = NaN;
            elseif isnumeric(x) && isscalar(x)
                vec(i) = x;
            else
                vec(i) = NaN;
            end
        end
    end
end

function val = get_struct_field(s, fieldName, defaultVal)
    if isstruct(s) && isfield(s, fieldName)
        val = s.(fieldName);
    else
        val = defaultVal;
    end
end
