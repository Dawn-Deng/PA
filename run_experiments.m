function allResults = run_experiments(config)
%RUN_EXPERIMENTS 工程化增强：批量实验入口（seeds/scenarios/参数扫描）
%
% config fields (optional):
%   seeds      : random seed list, default [20260321 20260322 20260323]
%   scenarios  : {'uniform','clustered','corridor'}
%   sweep      : struct array, e.g. struct('field','sigma2','values',[1e-4 3e-4])
%   baseParams : struct overrides
%   saveEach   : true/false
%   saveDir    : output dir

    if nargin < 1
        config = struct();
    end
    if ~isfield(config, 'seeds'), config.seeds = [20260321, 20260322, 20260323]; end
    if ~isfield(config, 'scenarios'), config.scenarios = {'uniform', 'clustered', 'corridor'}; end
    if ~isfield(config, 'sweep'), config.sweep = struct('field', {}, 'values', {}); end
    if ~isfield(config, 'baseParams'), config.baseParams = struct(); end
    if ~isfield(config, 'saveEach'), config.saveEach = false; end
    if ~isfield(config, 'saveDir'), config.saveDir = 'batch_results'; end

    if config.saveEach && ~exist(config.saveDir, 'dir')
        mkdir(config.saveDir);
    end

    sweepCases = expandSweep(config.sweep);
    caseIdx = 0;
    allResults = repmat(struct('meta', [], 'results', []), 0, 1);

    for s = 1:numel(config.scenarios)
        scenarioName = config.scenarios{s};
        for seed = config.seeds
            for c = 1:numel(sweepCases)
                caseIdx = caseIdx + 1;
                override = config.baseParams;
                override.randomSeed = seed;
                override.symbolSeed = seed + 101;
                override.verbosity = 0;
                if ~isfield(override, 'userGeneration') || isempty(override.userGeneration)
                    override.userGeneration = struct();
                end
                override.userGeneration.mode = scenarioName;

                fields = fieldnames(sweepCases(c));
                for f = 1:numel(fields)
                    override.(fields{f}) = sweepCases(c).(fields{f});
                end

                single = run_single_experiment(override);
                allResults(caseIdx).meta = struct('seed', seed, 'scenario', scenarioName, 'sweep', sweepCases(c));
                allResults(caseIdx).results = single;

                if config.saveEach
                    save(fullfile(config.saveDir, sprintf('exp_%04d.mat', caseIdx)), 'single');
                end
            end
        end
    end
end

function results = run_single_experiment(paramsOverride)
    [params, state, channelInfo] = Channel_model(paramsOverride);
    [state, initInfo] = Initialization(state, params);
    [state, aoInfo] = AO_model(state, params);
    [state, signalInfo] = Signal_model(state, params);
    problemInfo = Problem_formulation(state, params, aoInfo);

    results = struct();
    results.params = params;
    results.channelInfo = channelInfo;
    results.initInfo = initInfo;
    results.aoInfo = aoInfo;
    results.signalInfo = signalInfo;
    results.problemInfo = problemInfo;
    results.runtimeStats = aoInfo.runtimeStats;
    results.snapshots = struct( ...
        'finalX', state.X, 'finalTheta', state.theta, 'finalPhi', state.phi, 'finalS', state.S, ...
        'finalUsers', state.users, 'sumRateHistory', aoInfo.sumRateHistory);
end

function sweepCases = expandSweep(sweep)
    if isempty(sweep)
        sweepCases = struct();
        return;
    end

    sweepCases = struct();
    for i = 1:numel(sweep)
        f = sweep(i).field;
        vals = sweep(i).values;
        if i == 1
            sweepCases = repmat(struct(), numel(vals), 1);
            for j = 1:numel(vals)
                sweepCases(j).(f) = vals(j);
            end
        else
            prev = sweepCases;
            sweepCases = repmat(struct(), numel(prev) * numel(vals), 1);
            idx = 0;
            for p = 1:numel(prev)
                for j = 1:numel(vals)
                    idx = idx + 1;
                    sweepCases(idx) = prev(p);
                    sweepCases(idx).(f) = vals(j);
                end
            end
        end
    end
end
