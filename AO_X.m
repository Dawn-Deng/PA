function [state, info, memory] = AO_X(state, params, memory)
%AO_X Position-block update with detailed line-search trace.

    if nargin < 3 || isempty(memory)
        memory = initializePositionMemory(params);
    end

    acceptedCount = 0;
    waveguideInfo = repmat(struct( ...
        'waveguideIndex', [], ...
        'accepted', false, ...
        'rejectReason', '', ... % reject reason
        'lineSearchSteps', 0, ...
        'alphaAccepted', [], ...
        'sumRateBefore', [], ...
        'sumRateAfter', [], ...
        'gradCurrent', [], ...
        'direction', [], ...
        'lineSearchTrace', []), params.N, 1); % inner trace

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
        lineSearchSteps = 0;
        candidateStateBest = state;
        acceptedAlpha = [];
        rejectReason = 'noImprovement';
        lsTrace = repmat(struct('alphaTrial', [], 'xCandidate', [], 'candidateRate', [], ...
            'accepted', false, 'rejectReason', ''), 0, 1);

        while alpha >= params.positionLineSearchMin
            lineSearchSteps = lineSearchSteps + 1;
            xBar = xCurrent + alpha * direction;
            xCandidate = Channel_model('project_waveguide_positions', xBar, params);

            candidateState = state;
            candidateState.X(:, n) = xCandidate;
            candidateState = Channel_model('update_state', candidateState, params);
            candidateMetrics = Signal_model('evaluate', candidateState, params, state.W, state.S);

            thisAccepted = candidateMetrics.sumRate >= fCurrent + params.epsilonX;
            thisRejectReason = '';
            if ~thisAccepted
                if norm(xCandidate - xCurrent) <= 1e-12
                    thisRejectReason = 'projectionCollapsedMove'; % reject reason
                else
                    thisRejectReason = 'belowTolerance'; % reject reason
                end
            end

            lsTrace(end + 1) = struct( ...
                'alphaTrial', alpha, ...
                'xCandidate', xCandidate, ...
                'candidateRate', candidateMetrics.sumRate, ...
                'accepted', thisAccepted, ...
                'rejectReason', thisRejectReason); %#ok<AGROW>

            if thisAccepted
                candidateState.sinr = candidateMetrics.sinr;
                candidateState.rate = candidateMetrics.rate;
                candidateState.sumRate = candidateMetrics.sumRate;
                candidateStateBest = candidateState;
                accepted = true;
                acceptedAlpha = alpha;
                rejectReason = '';
                break;
            end

            alpha = alpha * params.positionLineSearchBeta;
        end

        if ~accepted && alpha < params.positionLineSearchMin
            rejectReason = 'stepTooSmall'; % reject reason
        end

        if accepted
            xNew = candidateStateBest.X(:, n);
            gradNew = numericalGradient(candidateStateBest, params, n, xNew);
            memory{n} = updateLbfgsMemory(memory{n}, xNew - xCurrent, gradNew - gradCurrent, params.positionMemory);
            state = candidateStateBest;
            acceptedCount = acceptedCount + 1;
        end

        waveguideInfo(n).waveguideIndex = n;
        waveguideInfo(n).accepted = accepted;
        waveguideInfo(n).rejectReason = rejectReason;
        waveguideInfo(n).lineSearchSteps = lineSearchSteps;
        waveguideInfo(n).alphaAccepted = acceptedAlpha;
        waveguideInfo(n).sumRateBefore = fCurrent;
        waveguideInfo(n).sumRateAfter = state.sumRate;
        waveguideInfo(n).gradCurrent = gradCurrent;
        waveguideInfo(n).direction = direction;
        waveguideInfo(n).lineSearchTrace = lsTrace;
    end

    info = struct( ...
        'acceptedCount', acceptedCount, ...
        'sumRate', state.sumRate, ...
        'waveguideInfo', waveguideInfo);
end

function memory = initializePositionMemory(params)
    memory = cell(params.N, 1);
    for n = 1:params.N
        memory{n}.S = [];
        memory{n}.Y = [];
    end
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
    historySize = size(S, 2);
    q = grad;
    alpha = zeros(historySize, 1);
    rho = zeros(historySize, 1);

    for idx = historySize:-1:1
        rho(idx) = 1 / max(Y(:, idx)' * S(:, idx), eps);
        alpha(idx) = rho(idx) * (S(:, idx)' * q);
        q = q - alpha(idx) * Y(:, idx);
    end

    gamma = (S(:, end)' * Y(:, end)) / max(Y(:, end)' * Y(:, end), eps);
    r = gamma * q;

    for idx = 1:historySize
        beta = rho(idx) * (Y(:, idx)' * r);
        r = r + S(:, idx) * (alpha(idx) - beta);
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
