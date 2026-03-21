function problemState = miso_problem_formulation(channelState, signalState, params)
%MISO_PROBLEM_FORMULATION 系统总优化目标与约束建模
% 对应论文中的联合优化问题：
%   max_{S, X, W, theta, phi} R_sum(S, X, W, theta, phi)
%
% 输出：
%   problemState.objective      -> 当前配置下的总频谱效率 R_sum
%   problemState.variables      -> 优化变量快照 S, X, W, theta, phi
%   problemState.constraints    -> 各约束的检查结果
%   problemState.description    -> 问题文本描述

    candidatePool = params.candidateUserPool(:).';
    KServ = min([params.NRF, params.Kmax, numel(candidatePool)]);

    S = signalState.scheduledUsers(:).';
    X = channelState.X;
    W = signalState.W;
    theta = channelState.theta;
    phi = channelState.phi;
    sinr = signalState.slotResult.sinr(:);
    rate = signalState.slotResult.rate(:);
    Rsum = sum(rate);

    constraints = evaluateConstraints(S, candidatePool, KServ, X, theta, phi, W, params);

    variables = struct();
    variables.S = S;
    variables.X = X;
    variables.W = W;
    variables.theta = theta;
    variables.phi = phi;
    variables.KServ = KServ;
    variables.candidatePool = candidatePool;

    description = struct();
    description.objective = 'maximize R_sum(S, X, W, theta, phi) = sum_{k in S} log2(1 + gamma_k)';
    description.problemType = 'mixed discrete-continuous strongly non-convex optimization';
    description.notes = {
        'S is selected from the candidate user pool C with fixed cardinality K_serv.'
        'X collects all PA positions on all waveguides.'
        'W is the digital precoding matrix under the total power constraint.'
        'theta and phi are the hardware-adjustable elevation/azimuth angles.'
    };

    problemState = struct();
    problemState.objective = Rsum;
    problemState.sinr = sinr;
    problemState.rate = rate;
    problemState.variables = variables;
    problemState.constraints = constraints;
    problemState.description = description;
end

function constraints = evaluateConstraints(S, candidatePool, KServ, X, theta, phi, W, params)
%EVALUATECONSTRAINTS 检查论文 (11b)-(11f) 约束是否满足
    spacingMatrix = diff(X, 1, 1);
    withinCandidatePool = all(ismember(S, candidatePool));
    fixedCardinality = numel(S) == KServ;
    minSpacingSatisfied = all(spacingMatrix(:) >= params.deltaMin - 1e-12);
    sortedPositions = all(spacingMatrix(:) >= -1e-12);
    withinDy = all(X(:) >= 0) && all(X(:) <= params.Dy);
    thetaRange = all(theta(:) >= pi / 2 - 1e-12) && all(theta(:) <= pi + 1e-12);
    phiRange = all(phi(:) > -pi - 1e-12) && all(phi(:) <= pi + 1e-12);
    transmitPower = real(trace(W * W'));
    powerSatisfied = transmitPower <= params.Pmax + 1e-12;

    constraints = struct();
    constraints.userSet = struct( ...
        'description', 'S subseteq C and |S| = K_serv', ...
        'withinCandidatePool', withinCandidatePool, ...
        'fixedCardinality', fixedCardinality, ...
        'satisfied', withinCandidatePool && fixedCardinality);
    constraints.spacing = struct( ...
        'description', 'y_{n,m} - y_{n,m-1} >= Delta', ...
        'minimumSpacing', min(spacingMatrix(:)), ...
        'requiredSpacing', params.deltaMin, ...
        'satisfied', minSpacingSatisfied);
    constraints.positionFeasibleSet = struct( ...
        'description', 'x_n in X_n with 0 <= y_{n,1} <= ... <= y_{n,M} <= D_y', ...
        'sortedPositions', sortedPositions, ...
        'withinDy', withinDy, ...
        'Dy', params.Dy, ...
        'satisfied', sortedPositions && withinDy);
    constraints.orientation = struct( ...
        'description', 'pi/2 <= theta <= pi and -pi < phi <= pi', ...
        'thetaSatisfied', thetaRange, ...
        'phiSatisfied', phiRange, ...
        'satisfied', thetaRange && phiRange);
    constraints.power = struct( ...
        'description', 'tr(W W^H) <= P_max', ...
        'transmitPower', transmitPower, ...
        'Pmax', params.Pmax, ...
        'satisfied', powerSatisfied);
    constraints.allSatisfied = constraints.userSet.satisfied && constraints.spacing.satisfied ...
        && constraints.positionFeasibleSet.satisfied && constraints.orientation.satisfied ...
        && constraints.power.satisfied;
end
