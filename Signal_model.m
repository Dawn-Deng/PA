function [state, signalInfo, aoInfo] = Signal_model(state, params)
%SIGNAL_MODEL 五文件版本中的“信号模型 + AO 主算法”统一实现。
% 包含：
% 1) 广义 BCD 外层循环；
% 2) Precoding WMMSE；
% 3) Angle update；
% 4) Position update；
% 5) User-set update；
% 6) 最终信号生成与 SINR/rate 输出。

    [state, aoInfo] = runAO(state, params);
    signalInfo = buildSignalOutput(state, params);
end

function [state, aoInfo] = runAO(state, params)
    blockHistory = repmat(struct('W', [], 'angle', [], 'position', [], 'userSet', []), params.Tmax, 1);
    sumRateHistory = zeros(params.Tmax + 1, 1);
    positionMemory = [];

    wInfo0 = [];
    state.sumRate = 0;
    state.sinr = zeros(numel(state.S), 1);
    state.rate = zeros(numel(state.S), 1);
    sumRateHistory(1) = state.sumRate;

    for t = 1:params.Tmax
        previousRate = state.sumRate;
        [state, blockHistory(t).W] = updatePrecodingWMMSE(state, params);
        [state, blockHistory(t).angle] = updateAngles(state, params);
        [state, blockHistory(t).position, positionMemory] = updatePositions(state, params, positionMemory);
        [state, blockHistory(t).userSet] = updateUserSet(state, params, t);
        sumRateHistory(t + 1) = state.sumRate;
        state.iteration = t;
        if abs(sumRateHistory(t + 1) - previousRate) < params.epsilonOuter
            aoInfo = finalizeAO(sumRateHistory(1:t + 1), blockHistory(1:t), wInfo0, true, t);
            return;
        end
    end

    aoInfo = finalizeAO(sumRateHistory, blockHistory, wInfo0, false, params.Tmax);
end

function aoInfo = finalizeAO(sumRateHistory, blockHistory, wInfo0, converged, iterations)
    aoInfo = struct();
    aoInfo.initialWInfo = wInfo0;
    aoInfo.sumRateHistory = sumRateHistory;
    aoInfo.blockHistory = blockHistory;
    aoInfo.converged = converged;
    aoInfo.iterations = iterations;
end

function metrics = evaluateSumRate(state, params, W, S)
    if nargin < 4 || isempty(S)
        S = state.S;
    end
    if isempty(S) || isempty(W)
        metrics = struct('sinr', [], 'rate', [], 'sumRate', 0, 'scheduledUsers', S, 'W', W);
        return;
    end
    Hs = state.channelMatrix(S, :);
    numStreams = numel(S);
    sinr = zeros(numStreams, 1);
    rate = zeros(numStreams, 1);
    for k = 1:numStreams
        coupling = Hs(k, :) * W;
        signalPower = abs(coupling(k))^2;
        interferencePower = sum(abs(coupling).^2) - signalPower;
        sinr(k) = signalPower / (interferencePower + params.sigma2);
        rate(k) = log2(1 + sinr(k));
    end
    metrics = struct('sinr', sinr, 'rate', rate, 'sumRate', sum(rate), 'scheduledUsers', S(:).', 'W', W);
end

function [state, info] = updateAngles(state, params)
    directions = [1,0; -1,0; 0,1; 0,-1; 1,1; 1,-1; -1,1; -1,-1];
    currentRate = state.sumRate;
    acceptedCount = 0;

    for n = 1:params.N
        for m = 1:params.M
            bestLocalState = state;
            bestLocalRate = currentRate;
            anchorAngles = buildAnchorSet(state, n, m);

            for a = 1:size(anchorAngles, 1)
                candidateState = state;
                candidateState.theta(m, n) = anchorAngles(a, 1);
                candidateState.phi(m, n) = anchorAngles(a, 2);
                candidateState = Channel_model('update_state', candidateState, params);
                candidateMetrics = evaluateSumRate(candidateState, params, state.W, state.S);
                if candidateMetrics.sumRate >= bestLocalRate + params.epsilonTheta
                    candidateState.sinr = candidateMetrics.sinr;
                    candidateState.rate = candidateMetrics.rate;
                    candidateState.sumRate = candidateMetrics.sumRate;
                    bestLocalState = candidateState;
                    bestLocalRate = candidateMetrics.sumRate;
                end
            end

            state = bestLocalState;
            currentRate = bestLocalRate;
            localStepTheta = params.angleStepThetaInit;
            localStepPhi = params.angleStepPhiInit;

            while localStepTheta > params.angleStepThetaMin || localStepPhi > params.angleStepPhiMin
                improved = false;
                currentTheta = state.theta(m, n);
                currentPhi = state.phi(m, n);
                for d = 1:size(directions, 1)
                    candTheta = currentTheta + directions(d, 1) * localStepTheta;
                    candPhi = currentPhi + directions(d, 2) * localStepPhi;
                    [candTheta, candPhi] = projectAngles(candTheta, candPhi);
                    candidateState = state;
                    candidateState.theta(m, n) = candTheta;
                    candidateState.phi(m, n) = candPhi;
                    candidateState = Channel_model('update_state', candidateState, params);
                    candidateMetrics = evaluateSumRate(candidateState, params, state.W, state.S);
                    if candidateMetrics.sumRate >= currentRate + params.epsilonTheta
                        candidateState.sinr = candidateMetrics.sinr;
                        candidateState.rate = candidateMetrics.rate;
                        candidateState.sumRate = candidateMetrics.sumRate;
                        state = candidateState;
                        currentRate = candidateMetrics.sumRate;
                        improved = true;
                        acceptedCount = acceptedCount + 1;
                        break;
                    end
                end
                if ~improved
                    localStepTheta = params.betaTheta * localStepTheta;
                    localStepPhi = params.betaPhi * localStepPhi;
                end
            end
        end
    end

    info = struct('acceptedCount', acceptedCount, 'sumRate', state.sumRate);
end

function anchors = buildAnchorSet(state, n, m)
    anchors = [state.theta(m, n), state.phi(m, n); pi, 0];
    paPos = state.paPositions(:, m, n);
    for idx = 1:numel(state.S)
        k = state.S(idx);
        direction = state.users(k, :).'- paPos;
        direction = direction / norm(direction);
        anchors = [anchors; acos(direction(3)), atan2(direction(2), direction(1))]; %#ok<AGROW>
    end
    for i = 1:size(anchors, 1)
        [anchors(i, 1), anchors(i, 2)] = projectAngles(anchors(i, 1), anchors(i, 2));
    end
    anchors = unique(round(anchors, 12), 'rows');
end

function [theta, phi] = projectAngles(theta, phi)
    theta = min(max(theta, pi / 2), pi);
    phi = mod(phi + pi, 2 * pi) - pi;
    if phi <= -pi
        phi = phi + 2 * pi;
    end
end

function [state, info, memory] = updatePositions(state, params, memory)
    if nargin < 3 || isempty(memory)
        memory = cell(params.N, 1);
        for n = 1:params.N
            memory{n}.S = [];
            memory{n}.Y = [];
        end
    end
    acceptedCount = 0;
    for n = 1:params.N
        xCurrent = state.X(:, n);
        fCurrent = state.sumRate;
        gradCurrent = numericalGradient(state, params, n, xCurrent);
        direction = -lbfgsDirection(gradCurrent, memory{n});
        if norm(direction) <= 1e-12
            direction = -gradCurrent;
        end

        alpha = params.positionLineSearchInit;
        accepted = false;
        candidateStateBest = state;
        while alpha >= params.positionLineSearchMin
            xBar = xCurrent + alpha * direction;
            xCandidate = Channel_model('project_waveguide_positions', xBar, params);
            candidateState = state;
            candidateState.X(:, n) = xCandidate;
            candidateState = Channel_model('update_state', candidateState, params);
            metrics = evaluateSumRate(candidateState, params, state.W, state.S);
            if metrics.sumRate >= fCurrent + params.epsilonX
                candidateState.sinr = metrics.sinr;
                candidateState.rate = metrics.rate;
                candidateState.sumRate = metrics.sumRate;
                candidateStateBest = candidateState;
                accepted = true;
                break;
            end
            alpha = alpha * params.positionLineSearchBeta;
        end

        if accepted
            xNew = candidateStateBest.X(:, n);
            gradNew = numericalGradient(candidateStateBest, params, n, xNew);
            memory{n} = updateLbfgsMemory(memory{n}, xNew - xCurrent, gradNew - gradCurrent, params.positionMemory);
            state = candidateStateBest;
            acceptedCount = acceptedCount + 1;
        end
    end
    info = struct('acceptedCount', acceptedCount, 'sumRate', state.sumRate);
end

function grad = numericalGradient(state, params, n, x)
    grad = zeros(params.M, 1);
    delta = params.positionFiniteDiff;
    for m = 1:params.M
        xp = x; xm = x;
        xp(m) = xp(m) + delta;
        xm(m) = xm(m) - delta;
        stateP = state; stateM = state;
        stateP.X(:, n) = Channel_model('project_waveguide_positions', xp, params);
        stateM.X(:, n) = Channel_model('project_waveguide_positions', xm, params);
        stateP = Channel_model('update_state', stateP, params);
        stateM = Channel_model('update_state', stateM, params);
        fp = evaluateSumRate(stateP, params, state.W, state.S).sumRate;
        fm = evaluateSumRate(stateM, params, state.W, state.S).sumRate;
        grad(m) = (fp - fm) / (2 * delta);
    end
end

function direction = lbfgsDirection(grad, memory)
    if isempty(memory.S)
        direction = grad;
        return;
    end
    S = memory.S; Y = memory.Y; m = size(S, 2);
    q = grad; alpha = zeros(m, 1); rho = zeros(m, 1);
    for i = m:-1:1
        rho(i) = 1 / max(Y(:, i)' * S(:, i), eps);
        alpha(i) = rho(i) * (S(:, i)' * q);
        q = q - alpha(i) * Y(:, i);
    end
    gamma = (S(:, end)' * Y(:, end)) / max(Y(:, end)' * Y(:, end), eps);
    r = gamma * q;
    for i = 1:m
        beta = rho(i) * (Y(:, i)' * r);
        r = r + S(:, i) * (alpha(i) - beta);
    end
    direction = r;
end

function memory = updateLbfgsMemory(memory, s, y, maxMemory)
    if s' * y <= 1e-12
        return;
    end
    if isempty(memory.S)
        memory.S = s;
        memory.Y = y;
        return;
    end
    memory.S = [memory.S, s];
    memory.Y = [memory.Y, y];
    if size(memory.S, 2) > maxMemory
        memory.S(:, 1) = [];
        memory.Y(:, 1) = [];
    end
end

function [state, info] = updateUserSet(state, params, t)
    info = struct('triggered', mod(t, params.TS) == 0, 'acceptedSwaps', 0, 'bestDelta', 0);
    if ~info.triggered
        return;
    end

    for swapIter = 1:params.maxSwapPerUpdate
        [~, orderAsc] = sort(state.rate, 'ascend');
        weakPositions = orderAsc(1:min(params.Lin, numel(state.S)));
        external = setdiff(state.candidatePool, state.S, 'stable');
        if isempty(external)
            break;
        end
        [~, extOrder] = sort(state.Emax(external), 'descend');
        strongExternal = external(extOrder(1:min(params.Lout, numel(external))));
        bestDelta = -inf;
        bestState = state;

        for w = 1:numel(weakPositions)
            pos = weakPositions(w);
            for j = 1:numel(strongExternal)
                candidateUsers = state.S;
                candidateUsers(pos) = strongExternal(j);
                metrics = evaluateSumRate(state, params, state.W, candidateUsers);
                delta = metrics.sumRate - state.sumRate;
                if delta > bestDelta
                    bestDelta = delta;
                    bestState = state;
                    bestState.S = candidateUsers;
                    bestState.sinr = metrics.sinr;
                    bestState.rate = metrics.rate;
                    bestState.sumRate = metrics.sumRate;
                end
            end
        end

        info.bestDelta = bestDelta;
        if bestDelta >= params.epsilonS
            state = bestState;
            info.acceptedSwaps = info.acceptedSwaps + 1;
        else
            break;
        end
    end
end

function signalInfo = buildSignalOutput(state, params)
    rng(params.symbolSeed);
    numStreams = numel(state.S);
    symbols = (randn(numStreams, 1) + 1j * randn(numStreams, 1)) / sqrt(2);
    noise = sqrt(params.sigma2 / 2) * (randn(numStreams, 1) + 1j * randn(numStreams, 1));
    xRad = state.W * symbols;
    Hs = state.channelMatrix(state.S, :);
    metrics = evaluateSumRate(state, params, state.W, state.S);
    y = zeros(numStreams, 1);
    userMetrics(numStreams, 1) = struct('userIndex', [], 'desiredSignal', [], 'interference', [], ...
        'noise', [], 'receivedSignal', [], 'sinr', [], 'rate', []);

    for k = 1:numStreams
        coupling = Hs(k, :) * state.W;
        desiredSignal = coupling(k) * symbols(k);
        interference = coupling * symbols - desiredSignal;
        y(k) = desiredSignal + interference + noise(k);
        userMetrics(k).userIndex = state.S(k);
        userMetrics(k).desiredSignal = desiredSignal;
        userMetrics(k).interference = interference;
        userMetrics(k).noise = noise(k);
        userMetrics(k).receivedSignal = y(k);
        userMetrics(k).sinr = metrics.sinr(k);
        userMetrics(k).rate = metrics.rate(k);
    end

    signalInfo = struct();
    signalInfo.S = state.S;
    signalInfo.W = state.W;
    signalInfo.channelMatrix = state.channelMatrix(state.S, :);
    signalInfo.symbols = symbols;
    signalInfo.noise = noise;
    signalInfo.xRad = xRad;
    signalInfo.y = y;
    signalInfo.sinr = metrics.sinr;
    signalInfo.rate = metrics.rate;
    signalInfo.sumRate = metrics.sumRate;
    signalInfo.userMetrics = userMetrics;
end
