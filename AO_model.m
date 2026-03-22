function [state, aoInfo] = AO_model(state, params)
%AO_MODEL 交替优化主过程，信号计算统一通过 Signal_model 接口完成。

    blockHistory = repmat(struct('W', [], 'angle', [], 'position', [], 'userSet', []), params.Tmax, 1);
    sumRateHistory = zeros(params.Tmax + 1, 1);
    positionMemory = [];
    angleMemory = [];

    state.sumRate = 0;
    state.sinr = zeros(numel(state.S), 1);
    state.rate = zeros(numel(state.S), 1);
    sumRateHistory(1) = state.sumRate;

    for t = 1:params.Tmax
        previousRate = state.sumRate;
        [state, blockHistory(t).W] = AO_W(state, params);
        [state, blockHistory(t).angle, angleMemory] = AO_angle(state, params, angleMemory);
        [state, blockHistory(t).position, positionMemory] = updatePositions(state, params, positionMemory);
        [state, blockHistory(t).userSet] = updateUserSet(state, params, t);
        sumRateHistory(t + 1) = state.sumRate;
        state.iteration = t;
        if abs(sumRateHistory(t + 1) - previousRate) < params.epsilonOuter
            aoInfo = finalizeAO(sumRateHistory(1:t + 1), blockHistory(1:t), true, t);
            return;
        end
    end

    aoInfo = finalizeAO(sumRateHistory, blockHistory, false, params.Tmax);
end

function aoInfo = finalizeAO(sumRateHistory, blockHistory, converged, iterations)
    aoInfo = struct();
    aoInfo.sumRateHistory = sumRateHistory;
    aoInfo.blockHistory = blockHistory;
    aoInfo.converged = converged;
    aoInfo.iterations = iterations;
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
            metrics = Signal_model('evaluate', candidateState, params, state.W, state.S);
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
        xp = x;
        xm = x;
        xp(m) = xp(m) + delta;
        xm(m) = xm(m) - delta;
        stateP = state;
        stateM = state;
        stateP.X(:, n) = Channel_model('project_waveguide_positions', xp, params);
        stateM.X(:, n) = Channel_model('project_waveguide_positions', xm, params);
        stateP = Channel_model('update_state', stateP, params);
        stateM = Channel_model('update_state', stateM, params);
        fp = Signal_model('evaluate', stateP, params, state.W, state.S).sumRate;
        fm = Signal_model('evaluate', stateM, params, state.W, state.S).sumRate;
        grad(m) = (fp - fm) / (2 * delta);
    end
end

function direction = lbfgsDirection(grad, memory)
    if isempty(memory.S)
        direction = grad;
        return;
    end
    S = memory.S;
    Y = memory.Y;
    m = size(S, 2);
    q = grad;
    alpha = zeros(m, 1);
    rho = zeros(m, 1);
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
                metrics = Signal_model('evaluate', state, params, state.W, candidateUsers);
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
