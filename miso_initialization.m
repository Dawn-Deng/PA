function [initState, params] = miso_initialization(channelState, params)
%MISO_INITIALIZATION 根据论文初始化策略进行建模
% 包括：
% 1) 潜在增益矩阵 G_pot 与候选用户池 C；
% 2) 全局匹配得到初始服务集合 S^(0)；
% 3) 初始位置矩阵 X^(0)；
% 4) 初始姿态角 (theta^(0), phi^(0))。

    KServ = min([params.NRF, params.Kmax, params.K]);
    yRef = buildReferencePositions(params);
    [Gpot, yStar, dStar, Emax] = computePotentialGainMatrix(channelState, params);
    candidatePool = selectCandidateUserPool(Emax, KServ);
    [utilityMatrix, movementCost] = buildUtilityMatrix(candidatePool, yRef, yStar, Gpot, params);
    matchingResult = solvePrimaryAssociation(utilityMatrix, candidatePool, params.N, params.M, KServ);
    [Xbar0, X0] = buildInitialPositionMatrix(matchingResult, yRef, yStar, params);
    [theta0, phi0, paPositions0] = buildInitialOrientations(channelState, matchingResult, X0, params);

    params.candidateUserPool = candidatePool;
    params.initialServiceSet = matchingResult.serviceSet;

    initState = struct();
    initState.KServ = KServ;
    initState.Gpot = Gpot;
    initState.yStar = yStar;
    initState.dStar = dStar;
    initState.Emax = Emax;
    initState.candidatePool = candidatePool;
    initState.referencePositions = yRef;
    initState.movementCost = movementCost;
    initState.utilityMatrix = utilityMatrix;
    initState.matching = matchingResult;
    initState.initialServiceSet = matchingResult.serviceSet;
    initState.Xbar0 = Xbar0;
    initState.X0 = X0;
    initState.theta0 = theta0;
    initState.phi0 = phi0;
    initState.paPositions0 = paPositions0;
end

function yRef = buildReferencePositions(params)
%BUILDREFERENCEPOSITIONS 构造等间隔参考部署位置 y_ref
    if params.M == 1
        yRef = 0;
        return;
    end

    yRef = zeros(params.M, 1);
    for m = 1:params.M
        yRef(m) = (m - 1) * params.deltaMin + (m - 1) / (params.M - 1) * ...
            (params.Dy - (params.M - 1) * params.deltaMin);
    end
end

function [Gpot, yStar, dStar, Emax] = computePotentialGainMatrix(channelState, params)
%COMPUTEPOTENTIALGAINMATRIX 计算潜在增益矩阵 G_pot 与 E_max
    users = channelState.users;
    waveguideX = channelState.waveguideFeedPoints(:, 1);
    K = size(users, 1);
    N = params.N;
    M = params.M;

    Gpot = zeros(K, N * M);
    yStar = zeros(K, N, M);
    dStar = zeros(K, N, M);

    denominator = (2 * log(params.alphaL))^2 - params.alphaW^2;
    denominator = max(denominator, eps);

    for k = 1:K
        xk = users(k, 1);
        yk = users(k, 2);
        zk = users(k, 3);

        for n = 1:N
            xW = waveguideX(n);
            Akn = (xk - xW)^2 + (zk - params.d)^2;
            gammaStar = sqrt(Akn * params.alphaW^2 / denominator);
            yStarKN = yk - gammaStar;
            userPaDistance = norm([xk; yk; zk] - [xW; yStarKN; params.d]);
            potentialGain = sqrt(1 / M) * exp(-(params.alphaW / 2) * yStarKN) * ...
                (params.alphaL ^ userPaDistance) * ...
                (params.lambda * params.nMedium * params.v * sqrt(2 * params.a * params.b) / ...
                (2 * userPaDistance));

            for m = 1:M
                columnIdx = (n - 1) * M + m;
                Gpot(k, columnIdx) = abs(potentialGain);
                yStar(k, n, m) = yStarKN;
                dStar(k, n, m) = userPaDistance;
            end
        end
    end

    Emax = max(Gpot, [], 2);
end

function candidatePool = selectCandidateUserPool(Emax, KServ)
%SELECTCANDIDATEUSERPOOL 选取前 2K_serv 个用户作为候选集合 C
    [~, order] = sort(Emax, 'descend');
    numCandidates = min(numel(Emax), 2 * KServ);
    candidatePool = sort(order(1:numCandidates)).';
end

function [utilityMatrix, movementCost] = buildUtilityMatrix(candidatePool, yRef, yStar, Gpot, params)
%BUILDUTILITYMATRIX 构造匹配效用 U_{k,n,m}
    numCandidates = numel(candidatePool);
    numPAs = params.N * params.M;
    utilityMatrix = zeros(numCandidates, numPAs);
    movementCost = zeros(numCandidates, params.N, params.M);

    for c = 1:numCandidates
        k = candidatePool(c);
        for n = 1:params.N
            for m = 1:params.M
                paIdx = (n - 1) * params.M + m;
                moveCost = abs(yRef(m) - yStar(k, n, m));
                movementCost(c, n, m) = moveCost;
                utilityMatrix(c, paIdx) = Gpot(k, paIdx) - params.lambdaMov * moveCost;
            end
        end
    end
end

function matchingResult = solvePrimaryAssociation(utilityMatrix, candidatePool, N, M, KServ)
%SOLVEPRIMARYASSOCIATION 求解初始主关联的最大权匹配
% 这里采用深度优先搜索穷举 K_serv 个不冲突的匹配对，适合当前示例规模。
    numCandidates = numel(candidatePool);
    numPAs = N * M;

    bestScore = -inf;
    bestPairs = zeros(KServ, 2);
    usedUsers = false(numCandidates, 1);
    usedPAs = false(numPAs, 1);
    currentPairs = zeros(KServ, 2);

    search(1, 0);

    serviceSet = sort(candidatePool(bestPairs(:, 1))).';
    assignmentTensor = zeros(numCandidates, N, M);
    pairList = repmat(struct('userIndex', [], 'waveguideIndex', [], 'paIndex', [], 'utility', []), KServ, 1);

    for i = 1:KServ
        candIdx = bestPairs(i, 1);
        paFlat = bestPairs(i, 2);
        n = ceil(paFlat / M);
        m = paFlat - (n - 1) * M;
        assignmentTensor(candIdx, n, m) = 1;
        pairList(i).userIndex = candidatePool(candIdx);
        pairList(i).waveguideIndex = n;
        pairList(i).paIndex = m;
        pairList(i).utility = utilityMatrix(candIdx, paFlat);
    end

    matchingResult = struct();
    matchingResult.totalUtility = bestScore;
    matchingResult.assignmentTensor = assignmentTensor;
    matchingResult.selectedPairs = pairList;
    matchingResult.serviceSet = serviceSet;

    function search(depth, currentScore)
        if depth > KServ
            if currentScore > bestScore
                bestScore = currentScore;
                bestPairs = currentPairs;
            end
            return;
        end

        for candIdx = 1:numCandidates
            if usedUsers(candIdx)
                continue;
            end
            usedUsers(candIdx) = true;
            for paFlat = 1:numPAs
                if usedPAs(paFlat)
                    continue;
                end
                usedPAs(paFlat) = true;
                currentPairs(depth, :) = [candIdx, paFlat];
                search(depth + 1, currentScore + utilityMatrix(candIdx, paFlat));
                usedPAs(paFlat) = false;
            end
            usedUsers(candIdx) = false;
        end
    end
end

function [Xbar0, X0] = buildInitialPositionMatrix(matchingResult, yRef, yStar, params)
%BUILDINITIALPOSITIONMATRIX 构造名义位置矩阵 Xbar0 并投影到可行域得到 X0
    Xbar0 = repmat(yRef, 1, params.N);

    for i = 1:numel(matchingResult.selectedPairs)
        pair = matchingResult.selectedPairs(i);
        Xbar0(pair.paIndex, pair.waveguideIndex) = yStar(pair.userIndex, pair.waveguideIndex, pair.paIndex);
    end

    X0 = zeros(size(Xbar0));
    for n = 1:params.N
        X0(:, n) = projectWaveguidePositions(Xbar0(:, n), params);
    end
end

function projected = projectWaveguidePositions(yBar, params)
%PROJECTWAVEGUIDEPOSITIONS 将名义位置投影到 0<=y1<=...<=yM<=Dy 且间距>=Delta 的可行域
    M = numel(yBar);
    y = sort(yBar(:));

    lowerCap = (0:M-1)' * params.deltaMin;
    upperCap = params.Dy - (M-1:-1:0)' * params.deltaMin;

    y = min(max(y, lowerCap), upperCap);

    for m = 2:M
        y(m) = max(y(m), y(m-1) + params.deltaMin);
    end

    if y(end) > params.Dy
        y(end) = params.Dy;
        for m = M-1:-1:1
            y(m) = min(y(m), y(m+1) - params.deltaMin);
        end
    end

    y(1) = max(y(1), 0);
    for m = 2:M
        y(m) = max(y(m), y(m-1) + params.deltaMin);
    end

    projected = min(y, upperCap);
end

function [theta0, phi0, paPositions0] = buildInitialOrientations(channelState, matchingResult, X0, params)
%BUILDINITIALORIENTATIONS 为匹配 PA 指向对应用户，其余 PA 垂直向下
    users = channelState.users;
    waveguideX = channelState.waveguideFeedPoints(:, 1);

    theta0 = pi * ones(params.M, params.N);
    phi0 = zeros(params.M, params.N);
    paPositions0 = zeros(3, params.M, params.N);

    for n = 1:params.N
        for m = 1:params.M
            paPositions0(:, m, n) = [waveguideX(n); X0(m, n); params.d];
        end
    end

    for i = 1:numel(matchingResult.selectedPairs)
        pair = matchingResult.selectedPairs(i);
        userPos = users(pair.userIndex, :).';
        paPos = paPositions0(:, pair.paIndex, pair.waveguideIndex);
        direction = userPos - paPos;
        direction = direction / norm(direction);
        theta0(pair.paIndex, pair.waveguideIndex) = acos(direction(3));
        phi0(pair.paIndex, pair.waveguideIndex) = atan2(direction(2), direction(1));
    end
end
