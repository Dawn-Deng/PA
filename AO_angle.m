function [state, info, memory] = AO_angle(state, params, memory)
%AO_ANGLE Angle-block update with detailed per-PA trace.

    if nargin < 3 || isempty(memory)
        memory = initializeAngleMemory(params);
    end

    directions = [1,0; -1,0; 0,1; 0,-1; 1,1; 1,-1; -1,1; -1,-1];
    userSetChanged = isempty(memory.previousServiceSet) || ~isequal(memory.previousServiceSet, state.S);
    maxRounds = inferAngleRoundLimit(params);
    reanchorFailureThreshold = getReanchorFailureThreshold(params);
    skipFailureThreshold = getSkipFailureThreshold(params);

    acceptedCount = 0;
    reanchorCount = 0;
    perPA = repmat(struct( ...
        'waveguideIndex', [], ...
        'paIndex', [], ...
        'reanchored', false, ...
        'acceptedMoves', 0, ...
        'rounds', 0, ...
        'finalTheta', [], ...
        'finalPhi', [], ...
        'finalSumRate', [], ...
        'failureStreakBefore', [], ...
        'failureStreakAfter', [], ...
        'trace', []), params.M, params.N); % inner trace

    for n = 1:params.N
        for m = 1:params.M
            paAcceptedMoves = 0;
            paRounds = 0;
            failureBefore = memory.failureStreak(m, n);
            shouldSkip = (~userSetChanged) && (memory.failureStreak(m, n) >= skipFailureThreshold);
            shouldReanchor = userSetChanged || memory.failureStreak(m, n) >= reanchorFailureThreshold;

            paTrace = struct();
            paTrace.anchorSet = [];
            paTrace.anchorRate = [];
            paTrace.anchorSelected = [state.theta(m, n), state.phi(m, n)];
            paTrace.anchorSelectedRate = state.sumRate;
            paTrace.roundTrace = repmat(struct( ...
                'roundIndex', [], ...
                'stepThetaBefore', [], ...
                'stepPhiBefore', [], ...
                'candidateAngles', [], ...
                'candidateRates', [], ...
                'bestRoundRate', [], ...
                'accepted', false, ...
                'rejectReason', ''), 0, 1);

            if shouldSkip
                memory.failureStreak(m, n) = memory.failureStreak(m, n) + 1;
                perPA(m, n).waveguideIndex = n;
                perPA(m, n).paIndex = m;
                perPA(m, n).reanchored = false;
                perPA(m, n).acceptedMoves = 0;
                perPA(m, n).rounds = 0;
                perPA(m, n).finalTheta = state.theta(m, n);
                perPA(m, n).finalPhi = state.phi(m, n);
                perPA(m, n).finalSumRate = state.sumRate;
                perPA(m, n).failureStreakBefore = failureBefore;
                perPA(m, n).failureStreakAfter = memory.failureStreak(m, n);
                perPA(m, n).trace = paTrace;
                continue;
            end

            if shouldReanchor
                [anchorSet, anchorRate] = evaluateAnchorSet(state, params, m, n, buildAnchorSet(state, params, m, n));
                [state, bestAnchor] = selectBestAnchorFromEvaluation(state, params, m, n, anchorSet, anchorRate);
                reanchorCount = reanchorCount + 1;
                paTrace.anchorSet = anchorSet;
                paTrace.anchorRate = anchorRate;
                paTrace.anchorSelected = bestAnchor;
                paTrace.anchorSelectedRate = state.sumRate;
            end

            centerTheta = state.theta(m, n);
            centerPhi = state.phi(m, n);
            centerRate = state.sumRate;
            stepTheta = params.angleStepThetaInit;
            stepPhi = params.angleStepPhiInit;

            while paRounds < maxRounds && (stepTheta > params.angleStepThetaMin || stepPhi > params.angleStepPhiMin)
                paRounds = paRounds + 1;

                [bestRoundState, bestRoundRate, foundBetter, candidateAngles, candidateRates] = ...
                    evaluatePollingSet(state, params, m, n, centerTheta, centerPhi, stepTheta, stepPhi, directions);

                stepThetaBefore = stepTheta;
                stepPhiBefore = stepPhi;
                acceptedThisRound = false;
                rejectReason = '';
                if foundBetter && bestRoundRate >= centerRate + params.epsilonTheta
                    state = bestRoundState;
                    centerTheta = state.theta(m, n);
                    centerPhi = state.phi(m, n);
                    centerRate = state.sumRate;
                    paAcceptedMoves = paAcceptedMoves + 1;
                    acceptedCount = acceptedCount + 1;
                    acceptedThisRound = true;
                else
                    if ~foundBetter
                        rejectReason = 'noBetterCandidate'; % reject reason
                    else
                        rejectReason = 'belowEpsilonTheta'; % reject reason
                    end
                    stepTheta = params.betaTheta * stepTheta;
                    stepPhi = params.betaPhi * stepPhi;
                end

                thisRound = struct( ...
                    'roundIndex', paRounds, ...
                    'stepThetaBefore', stepThetaBefore, ...
                    'stepPhiBefore', stepPhiBefore, ...
                    'candidateAngles', candidateAngles, ...
                    'candidateRates', candidateRates, ...
                    'bestRoundRate', bestRoundRate, ...
                    'accepted', acceptedThisRound, ...
                    'rejectReason', rejectReason);

                if isempty(paTrace.roundTrace)
                    paTrace.roundTrace = thisRound;
                else
                    paTrace.roundTrace(end + 1) = thisRound; %#ok<AGROW>
                end
            end

            if paAcceptedMoves > 0
                memory.failureStreak(m, n) = 0;
            else
                memory.failureStreak(m, n) = memory.failureStreak(m, n) + 1;
            end

            perPA(m, n).waveguideIndex = n;
            perPA(m, n).paIndex = m;
            perPA(m, n).reanchored = shouldReanchor;
            perPA(m, n).acceptedMoves = paAcceptedMoves;
            perPA(m, n).rounds = paRounds;
            perPA(m, n).finalTheta = state.theta(m, n);
            perPA(m, n).finalPhi = state.phi(m, n);
            perPA(m, n).finalSumRate = state.sumRate;
            perPA(m, n).failureStreakBefore = failureBefore;
            perPA(m, n).failureStreakAfter = memory.failureStreak(m, n);
            perPA(m, n).trace = paTrace;
        end
    end

    memory.previousServiceSet = state.S;
    info = struct( ...
        'acceptedCount', acceptedCount, ...
        'reanchorCount', reanchorCount, ...
        'userSetChangedTrigger', userSetChanged, ...
        'sumRate', state.sumRate, ...
        'perPA', perPA);
end

function memory = initializeAngleMemory(params)
    memory = struct();
    memory.previousServiceSet = [];
    memory.failureStreak = zeros(params.M, params.N);
end

function maxRounds = inferAngleRoundLimit(params)
    if isfield(params, 'angleMaxRounds') && ~isempty(params.angleMaxRounds)
        maxRounds = max(1, round(params.angleMaxRounds));
        return;
    end
    if isfield(params, 'ITheta')
        maxRounds = max(1, round(params.ITheta));
        return;
    end

    maxRoundsTheta = estimateContractionRounds(params.angleStepThetaInit, params.angleStepThetaMin, params.betaTheta);
    maxRoundsPhi = estimateContractionRounds(params.angleStepPhiInit, params.angleStepPhiMin, params.betaPhi);
    maxRounds = max(maxRoundsTheta, maxRoundsPhi) + 1;
end

function rounds = estimateContractionRounds(stepInit, stepMin, beta)
    if stepInit <= stepMin || beta <= 0 || beta >= 1
        rounds = 1;
        return;
    end
    rounds = ceil(log(stepMin / stepInit) / log(beta));
    rounds = max(rounds, 1);
end

function threshold = getReanchorFailureThreshold(params)
    if isfield(params, 'angleReanchorFailureThreshold')
        threshold = params.angleReanchorFailureThreshold;
    else
        threshold = 4;
    end
end

function threshold = getSkipFailureThreshold(params)
    if isfield(params, 'angleSkipFailureThreshold')
        threshold = params.angleSkipFailureThreshold;
    else
        threshold = 3;
    end
end

function [anchorSet, anchorRate] = evaluateAnchorSet(state, params, m, n, anchorSet)
    anchorRate = -inf(size(anchorSet, 1), 1);
    for idx = 1:size(anchorSet, 1)
        candidateState = setSinglePAAngle(state, params, m, n, anchorSet(idx, 1), anchorSet(idx, 2));
        candidateMetrics = Signal_model('evaluate', candidateState, params, state.W, state.S);
        anchorRate(idx) = candidateMetrics.sumRate;
    end
end

function [bestState, bestAnchor] = selectBestAnchorFromEvaluation(state, params, m, n, anchorSet, anchorRate)
    [~, idx] = max(anchorRate);
    bestAnchor = anchorSet(idx, :);
    bestState = setSinglePAAngle(state, params, m, n, bestAnchor(1), bestAnchor(2));
    metrics = Signal_model('evaluate', bestState, params, state.W, state.S);
    bestState.sinr = metrics.sinr;
    bestState.rate = metrics.rate;
    bestState.sumRate = metrics.sumRate;
end

function [bestState, bestRate, foundBetter, candidateAngles, candidateRates] = evaluatePollingSet(state, params, m, n, centerTheta, centerPhi, stepTheta, stepPhi, directions)
    bestState = state;
    bestRate = state.sumRate;
    foundBetter = false;

    numDirections = size(directions, 1);
    rawCandidates = zeros(numDirections, 2);
    for d = 1:numDirections
        candidateTheta = centerTheta + directions(d, 1) * stepTheta;
        candidatePhi = centerPhi + directions(d, 2) * stepPhi;
        [candidateTheta, candidatePhi] = projectAngles(candidateTheta, candidatePhi);
        rawCandidates(d, :) = [candidateTheta, candidatePhi];
    end

    dedupKey = round(rawCandidates, 12);
    [~, uniqueIdx] = unique(dedupKey, 'rows', 'stable');
    uniqueCandidates = rawCandidates(uniqueIdx, :);

    candidateAngles = uniqueCandidates;
    candidateRates = zeros(size(uniqueCandidates, 1), 1);

    for d = 1:size(uniqueCandidates, 1)
        candidateTheta = uniqueCandidates(d, 1);
        candidatePhi = uniqueCandidates(d, 2);

        candidateState = setSinglePAAngle(state, params, m, n, candidateTheta, candidatePhi);
        candidateMetrics = Signal_model('evaluate', candidateState, params, state.W, state.S);
        candidateRates(d) = candidateMetrics.sumRate;
        if candidateMetrics.sumRate > bestRate
            candidateState.sinr = candidateMetrics.sinr;
            candidateState.rate = candidateMetrics.rate;
            candidateState.sumRate = candidateMetrics.sumRate;
            bestState = candidateState;
            bestRate = candidateMetrics.sumRate;
            foundBetter = true;
        end
    end
end

function anchorSet = buildAnchorSet(state, params, m, n)
    anchorSet = [state.theta(m, n), state.phi(m, n); pi, 0];
    paPos = state.paPositions(:, m, n);
    representativeUsers = selectRepresentativeUsers(state, paPos, params);

    for idx = 1:numel(representativeUsers)
        userIndex = representativeUsers(idx);
        direction = state.users(userIndex, :).'- paPos;
        direction = direction / max(norm(direction), eps);
        thetaAnchor = acos(direction(3));
        phiAnchor = atan2(direction(2), direction(1));
        anchorSet = [anchorSet; thetaAnchor, phiAnchor]; %#ok<AGROW>
    end

    for idx = 1:size(anchorSet, 1)
        [anchorSet(idx, 1), anchorSet(idx, 2)] = projectAngles(anchorSet(idx, 1), anchorSet(idx, 2));
    end

    anchorSet = unique(round(anchorSet, 12), 'rows');
end

function representativeUsers = selectRepresentativeUsers(state, paPos, params)
    distances = zeros(numel(state.S), 1);
    for idx = 1:numel(state.S)
        userIndex = state.S(idx);
        distances(idx) = norm(state.users(userIndex, :).'- paPos);
    end
    [~, order] = sort(distances, 'ascend');

    if isfield(params, 'angleAnchorUsers')
        maxUsers = params.angleAnchorUsers;
    else
        maxUsers = 2;
    end
    representativeUsers = state.S(order(1:min(maxUsers, numel(order))));
end

function candidateState = setSinglePAAngle(state, params, m, n, theta, phi)
    candidateState = state;
    candidateState.theta(m, n) = theta;
    candidateState.phi(m, n) = phi;
    candidateState = Channel_model('update_state', candidateState, params);
end

function [theta, phi] = projectAngles(theta, phi)
    theta = min(max(theta, pi / 2), pi);
    phi = mod(phi + pi, 2 * pi) - pi;
    if phi <= -pi
        phi = phi + 2 * pi;
    end
end
