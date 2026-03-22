function [state, info, memory] = AO_angle(state, params, memory)
%AO_ANGLE Angle-block update in the alternating-optimization loop.
% With (S^(t), X^(t), W^(t+1)) fixed, this function updates the PA angles
% by a two-dimensional derivative-free pattern search together with a
% re-anchoring mechanism, matching the paper's angle-subproblem modeling.

    if nargin < 3 || isempty(memory)
        memory = initializeAngleMemory(params);
    end

    directions = [1,0; -1,0; 0,1; 0,-1; 1,1; 1,-1; -1,1; -1,-1];
    userSetChanged = isempty(memory.previousServiceSet) || ~isequal(memory.previousServiceSet, state.S);
    maxRounds = inferAngleRoundLimit(params);
    reanchorFailureThreshold = getReanchorFailureThreshold(params);

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
        'finalSumRate', []), params.M, params.N);

    for n = 1:params.N
        for m = 1:params.M
            paAcceptedMoves = 0;
            paRounds = 0;
            shouldReanchor = userSetChanged || memory.failureStreak(m, n) >= reanchorFailureThreshold;

            if shouldReanchor
                anchorSet = buildAnchorSet(state, params, m, n);
                [state, ~] = selectBestAnchor(state, params, m, n, anchorSet);
                reanchorCount = reanchorCount + 1;
            end

            centerTheta = state.theta(m, n);
            centerPhi = state.phi(m, n);
            centerRate = state.sumRate;
            stepTheta = params.angleStepThetaInit;
            stepPhi = params.angleStepPhiInit;

            while paRounds < maxRounds && (stepTheta > params.angleStepThetaMin || stepPhi > params.angleStepPhiMin)
                paRounds = paRounds + 1;
                [bestRoundState, bestRoundRate, foundBetter] = ...
                    evaluatePollingSet(state, params, m, n, centerTheta, centerPhi, stepTheta, stepPhi, directions);

                if foundBetter && bestRoundRate >= centerRate + params.epsilonTheta
                    state = bestRoundState;
                    centerTheta = state.theta(m, n);
                    centerPhi = state.phi(m, n);
                    centerRate = state.sumRate;
                    paAcceptedMoves = paAcceptedMoves + 1;
                    acceptedCount = acceptedCount + 1;
                else
                    stepTheta = params.betaTheta * stepTheta;
                    stepPhi = params.betaPhi * stepPhi;
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
    if isfield(params, 'ITheta')
        maxRounds = params.ITheta;
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
        threshold = 2;
    end
end

function [bestState, bestRate] = selectBestAnchor(state, params, m, n, anchorSet)
% Re-anchoring: choose the best search center from the anchor set using the
% true sum rate. The current angle pair is included in anchorSet, so the
% selected anchor is non-decreasing relative to the current center.
    bestState = state;
    bestRate = state.sumRate;
    for idx = 1:size(anchorSet, 1)
        candidateState = setSinglePAAngle(state, params, m, n, anchorSet(idx, 1), anchorSet(idx, 2));
        candidateMetrics = Signal_model('evaluate', candidateState, params, state.W, state.S);
        if candidateMetrics.sumRate >= bestRate
            candidateState.sinr = candidateMetrics.sinr;
            candidateState.rate = candidateMetrics.rate;
            candidateState.sumRate = candidateMetrics.sumRate;
            bestState = candidateState;
            bestRate = candidateMetrics.sumRate;
        end
    end
end

function [bestState, bestRate, foundBetter] = evaluatePollingSet(state, params, m, n, centerTheta, centerPhi, stepTheta, stepPhi, directions)
% Polling set:
%   P_p^(t,r) = { Pi_Omega(center + direction .* step) : direction in D }.
    bestState = state;
    bestRate = state.sumRate;
    foundBetter = false;

    for d = 1:size(directions, 1)
        candidateTheta = centerTheta + directions(d, 1) * stepTheta;
        candidatePhi = centerPhi + directions(d, 2) * stepPhi;
        [candidateTheta, candidatePhi] = projectAngles(candidateTheta, candidatePhi);

        candidateState = setSinglePAAngle(state, params, m, n, candidateTheta, candidatePhi);
        candidateMetrics = Signal_model('evaluate', candidateState, params, state.W, state.S);
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
% Anchor set B_p^(t):
%   { current angle, representative-user directions, (pi, 0) }.
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
        maxUsers = 3;
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
% Projection onto Omega:
%   pi/2 <= theta <= pi,  -pi < phi <= pi.
    theta = min(max(theta, pi / 2), pi);
    phi = mod(phi + pi, 2 * pi) - pi;
    if phi <= -pi
        phi = phi + 2 * pi;
    end
end
