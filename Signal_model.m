function varargout = Signal_model(varargin)
%SIGNAL_MODEL 仅负责信号生成、SINR/速率计算与相关输出。
% 用法：
%   [state, signalInfo] = Signal_model(state, params)
%   metrics = Signal_model('evaluate', state, params)
%   metrics = Signal_model('evaluate', state, params, W, S)

    if nargin >= 1 && (ischar(varargin{1}) || isstring(varargin{1}))
        action = lower(string(varargin{1}));
        switch action
            case "evaluate"
                state = varargin{2};
                params = varargin{3};
                W = [];
                S = [];
                if nargin >= 4
                    W = varargin{4};
                end
                if nargin >= 5
                    S = varargin{5};
                end
                varargout = {evaluateSignalMetrics(state, params, W, S)};
            otherwise
                error('Signal_model: 未知操作 %s', action);
        end
        return;
    end

    state = varargin{1};
    params = varargin{2};
    [state, signalInfo] = generateSignalOutput(state, params);
    varargout = {state, signalInfo};
end

function [state, signalInfo] = generateSignalOutput(state, params)
    metrics = evaluateSignalMetrics(state, params, state.W, state.S);
    state.sinr = metrics.sinr;
    state.rate = metrics.rate;
    state.sumRate = metrics.sumRate;

    numStreams = numel(metrics.scheduledUsers);
    channelMatrix = state.channelMatrix(metrics.scheduledUsers, :);
    symbols = [];
    noise = [];
    xRad = [];
    y = [];
    userMetrics = struct('userIndex', {}, 'desiredSignal', {}, 'interference', {}, ...
        'noise', {}, 'receivedSignal', {}, 'sinr', {}, 'rate', {});

    if numStreams > 0 && ~isempty(state.W)
        rng(params.symbolSeed);
        symbols = (randn(numStreams, 1) + 1j * randn(numStreams, 1)) / sqrt(2);
        noise = sqrt(params.sigma2 / 2) * (randn(numStreams, 1) + 1j * randn(numStreams, 1));
        xRad = state.W * symbols;
        y = zeros(numStreams, 1);
        userMetrics(numStreams, 1) = struct('userIndex', [], 'desiredSignal', [], 'interference', [], ...
            'noise', [], 'receivedSignal', [], 'sinr', [], 'rate', []);

        for k = 1:numStreams
            coupling = channelMatrix(k, :) * state.W;
            desiredSignal = coupling(k) * symbols(k);
            interference = coupling * symbols - desiredSignal;
            y(k) = desiredSignal + interference + noise(k);
            userMetrics(k).userIndex = metrics.scheduledUsers(k);
            userMetrics(k).desiredSignal = desiredSignal;
            userMetrics(k).interference = interference;
            userMetrics(k).noise = noise(k);
            userMetrics(k).receivedSignal = y(k);
            userMetrics(k).sinr = metrics.sinr(k);
            userMetrics(k).rate = metrics.rate(k);
        end
    end

    signalInfo = struct();
    signalInfo.S = metrics.scheduledUsers(:).';
    signalInfo.W = state.W;
    signalInfo.channelMatrix = channelMatrix;
    signalInfo.symbols = symbols;
    signalInfo.noise = noise;
    signalInfo.xRad = xRad;
    signalInfo.y = y;
    signalInfo.sinr = metrics.sinr;
    signalInfo.rate = metrics.rate;
    signalInfo.sumRate = metrics.sumRate;
    signalInfo.userMetrics = userMetrics;
end

function metrics = evaluateSignalMetrics(state, params, W, S)
    if nargin < 4 || isempty(S)
        S = state.S;
    end
    if nargin < 3 || isempty(W)
        W = state.W;
    end
    if isempty(S) || isempty(W)
        metrics = struct('sinr', zeros(numel(S), 1), 'rate', zeros(numel(S), 1), ...
            'sumRate', 0, 'scheduledUsers', S(:).', 'W', W);
        return;
    end

    Hs = state.channelMatrix(S, :);
    if isempty(Hs) || any(~isfinite(real(Hs(:)))) || any(~isfinite(imag(Hs(:))))
        % 数值鲁棒性保护：异常信道时返回稳定零指标
        metrics = struct('sinr', zeros(numel(S), 1), 'rate', zeros(numel(S), 1), ...
            'sumRate', 0, 'scheduledUsers', S(:).', 'W', W);
        return;
    end
    numStreams = numel(S);
    sinr = zeros(numStreams, 1);
    rate = zeros(numStreams, 1);
    for k = 1:numStreams
        coupling = Hs(k, :) * W;
        signalPower = abs(coupling(k))^2;
        interferencePower = sum(abs(coupling).^2) - signalPower;
        interferencePower = max(interferencePower, 0); % 数值鲁棒性保护
        denom = interferencePower + params.sigma2;
        sinr(k) = signalPower / max(denom, eps);
        rate(k) = log2(1 + sinr(k));
    end
    sumRate = sum(rate);
    if any(~isfinite(sinr)) || any(~isfinite(rate)) || ~isfinite(sumRate)
        % 数值鲁棒性保护：检测 NaN/Inf 并回退到稳定值
        sinr = zeros(numStreams, 1);
        rate = zeros(numStreams, 1);
        sumRate = 0;
    end

    metrics = struct('sinr', sinr, 'rate', rate, 'sumRate', sumRate, ...
        'scheduledUsers', S(:).', 'W', W);
end
