function [state, initInfo] = Initialization(state, params)
%INITIALIZATION 论文初始化阶段：G_pot, Emax, C, y_ref, U, Hungarian 匹配, X^(0), theta^(0), phi^(0)

    KServ = params.KServ;
    yRef = Channel_model('reference_positions', params);
    [Gpot, yStar, dStar, Emax] = computePotentialGain(state, params);
    candidatePool = selectCandidatePool(Emax, KServ);
    [utilityMatrix, movementCost] = buildUtility(candidatePool, yRef, yStar, Gpot, params);
    matching = solvePrimaryAssociationHungarian(utilityMatrix, candidatePool, params.N, params.M, KServ);
    [Xbar0, X0] = buildInitialPositions(matching, yRef, yStar, params);
    [theta0, phi0] = buildInitialAngles(state, matching, X0, params);

    state.candidatePool = candidatePool;
    state.Emax = Emax;
    state.S = matching.serviceSet;
    state.X = X0;
    state.theta = theta0;
    state.phi = phi0;
    state.W = [];
    state = Channel_model('update_state', state, params);

    initInfo = struct();
    initInfo.Gpot = Gpot;
    initInfo.yStar = yStar;
    initInfo.dStar = dStar;
    initInfo.Emax = Emax;
    initInfo.candidatePool = candidatePool;
    initInfo.referencePositions = yRef;
    initInfo.movementCost = movementCost;
    initInfo.utilityMatrix = utilityMatrix;
    initInfo.matching = matching;
    initInfo.Xbar0 = Xbar0;
    initInfo.X0 = X0;
    initInfo.theta0 = theta0;
    initInfo.phi0 = phi0;
end

function [Gpot, yStar, dStar, Emax] = computePotentialGain(state, params)
    K = size(state.users, 1);
    numPAs = params.N * params.M;
    Gpot = zeros(K, numPAs);
    yStar = zeros(K, params.N, params.M);
    dStar = zeros(K, params.N, params.M);
    waveguideX = state.waveguideFeedPoints(:, 1);

    denominator = (2 * log(params.alphaL))^2 - params.alphaW^2;
    if denominator <= 0
        error('Eq.(13) denominator is non-positive.');
    end

    for k = 1:K
        xk = state.users(k, 1); yk = state.users(k, 2); zk = state.users(k, 3);
        for n = 1:params.N
            xW = waveguideX(n);
            Akn = (xk - xW)^2 + (zk - params.d)^2;
            gammaStar = sqrt(Akn * params.alphaW^2 / denominator);
            yOpt = yk - gammaStar;
            dOpt = norm([xk; yk; zk] - [xW; yOpt; params.d]);
            gainOpt = sqrt(1 / params.M) * exp(-(params.alphaW / 2) * yOpt) * ...
                (params.alphaL ^ dOpt) * (params.lambda * params.nMedium * params.v * sqrt(2 * params.a * params.b) / (2 * dOpt));
            for m = 1:params.M
                idx = (n - 1) * params.M + m;
                Gpot(k, idx) = abs(gainOpt);
                yStar(k, n, m) = yOpt;
                dStar(k, n, m) = dOpt;
            end
        end
    end
    Emax = max(Gpot, [], 2);
end

function candidatePool = selectCandidatePool(Emax, KServ)
    [~, order] = sort(Emax, 'descend');
    candidatePool = order(1:min(numel(order), 2 * KServ)).';
end

function [utilityMatrix, movementCost] = buildUtility(candidatePool, yRef, yStar, Gpot, params)
    Nc = numel(candidatePool);
    P = params.N * params.M;
    utilityMatrix = zeros(Nc, P);
    movementCost = zeros(Nc, params.N, params.M);
    for c = 1:Nc
        k = candidatePool(c);
        for n = 1:params.N
            for m = 1:params.M
                idx = (n - 1) * params.M + m;
                dMov = abs(yRef(m) - yStar(k, n, m));
                movementCost(c, n, m) = dMov;
                utilityMatrix(c, idx) = Gpot(k, idx) - params.lambdaMov * dMov;
            end
        end
    end
end

function matching = solvePrimaryAssociationHungarian(utilityMatrix, candidatePool, N, M, KServ)
% (k,n,m) -> 二维 Hungarian 映射：
% 行 = 真实用户 Nc + dummy users (P-KServ)
% 列 = 真实 PA P + dummy PAs (Nc-KServ)
% 且 dummy-user -> dummy-PA 禁止，从而保证恰好 KServ 对真实匹配。
    Nc = numel(candidatePool);
    P = N * M;
    numDummyUsers = P - KServ;
    squareSize = Nc + numDummyUsers;
    prohibitCost = 1e9;

    cost = zeros(squareSize, squareSize);
    cost(1:Nc, 1:P) = -utilityMatrix;
    cost(1:Nc, P+1:end) = 0;
    cost(Nc+1:end, 1:P) = 0;
    cost(Nc+1:end, P+1:end) = prohibitCost;

    assignment = hungarianAssignment(cost);
    realPairs = [];
    for row = 1:Nc
        col = assignment(row);
        if col >= 1 && col <= P
            realPairs = [realPairs; row, col]; %#ok<AGROW>
        end
    end

    if size(realPairs, 1) ~= KServ
        error('Hungarian mapping failed to produce exactly K_serv real matches.');
    end

    assignmentTensor = zeros(Nc, N, M);
    selectedPairs = repmat(struct('userIndex', [], 'waveguideIndex', [], 'paIndex', [], 'utility', []), KServ, 1);
    serviceSet = zeros(1, KServ);

    for i = 1:KServ
        candRow = realPairs(i, 1);
        paFlat = realPairs(i, 2);
        n = ceil(paFlat / M);
        m = paFlat - (n - 1) * M;
        assignmentTensor(candRow, n, m) = 1;
        selectedPairs(i).userIndex = candidatePool(candRow);
        selectedPairs(i).waveguideIndex = n;
        selectedPairs(i).paIndex = m;
        selectedPairs(i).utility = utilityMatrix(candRow, paFlat);
        serviceSet(i) = candidatePool(candRow);
    end

    matching = struct();
    matching.assignmentTensor = assignmentTensor;
    matching.selectedPairs = selectedPairs;
    matching.serviceSet = serviceSet;
    matching.mappingExplanation = 'Rows: candidate users + dummy users; Columns: real PAs + dummy PAs; dummy-user to dummy-PA is forbidden.';
end

function assignment = hungarianAssignment(costMatrix)
    [nRows, nCols] = size(costMatrix);
    if nRows ~= nCols
        error('Hungarian requires a square matrix.');
    end
    n = nRows;
    u = zeros(n + 1, 1); v = zeros(n + 1, 1);
    p = zeros(n + 1, 1); way = zeros(n + 1, 1);

    for i = 1:n
        p(1) = i;
        j0 = 1;
        minv = inf(n + 1, 1);
        used = false(n + 1, 1);
        while true
            used(j0) = true;
            i0 = p(j0);
            delta = inf; j1 = 1;
            for j = 2:n + 1
                if ~used(j)
                    cur = costMatrix(i0, j - 1) - u(i0 + 1) - v(j);
                    if cur < minv(j)
                        minv(j) = cur;
                        way(j) = j0;
                    end
                    if minv(j) < delta
                        delta = minv(j);
                        j1 = j;
                    end
                end
            end
            for j = 1:n + 1
                if used(j)
                    u(p(j) + 1) = u(p(j) + 1) + delta;
                    v(j) = v(j) - delta;
                else
                    minv(j) = minv(j) - delta;
                end
            end
            j0 = j1;
            if p(j0) == 0
                break;
            end
        end
        while true
            j1 = way(j0);
            p(j0) = p(j1);
            j0 = j1;
            if j0 == 1
                break;
            end
        end
    end

    assignment = zeros(n, 1);
    for j = 2:n + 1
        if p(j) > 0
            assignment(p(j)) = j - 1;
        end
    end
end

function [Xbar0, X0] = buildInitialPositions(matching, yRef, yStar, params)
    Xbar0 = repmat(yRef, 1, params.N);
    for i = 1:numel(matching.selectedPairs)
        pair = matching.selectedPairs(i);
        Xbar0(pair.paIndex, pair.waveguideIndex) = yStar(pair.userIndex, pair.waveguideIndex, pair.paIndex);
    end
    X0 = zeros(size(Xbar0));
    for n = 1:params.N
        X0(:, n) = Channel_model('project_waveguide_positions', Xbar0(:, n), params);
    end
end

function [theta0, phi0] = buildInitialAngles(state, matching, X0, params)
    theta0 = pi * ones(params.M, params.N);
    phi0 = zeros(params.M, params.N);
    waveguideX = state.waveguideFeedPoints(:, 1);
    for i = 1:numel(matching.selectedPairs)
        pair = matching.selectedPairs(i);
        userPos = state.users(pair.userIndex, :).';
        paPos = [waveguideX(pair.waveguideIndex); X0(pair.paIndex, pair.waveguideIndex); params.d];
        direction = userPos - paPos;
        direction = direction / norm(direction);
        theta0(pair.paIndex, pair.waveguideIndex) = acos(direction(3));
        phi0(pair.paIndex, pair.waveguideIndex) = atan2(direction(2), direction(1));
    end
end
