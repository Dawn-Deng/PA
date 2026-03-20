function signal = miso_signal_model(H, signalConfig)
%MISO_SIGNAL_MODEL 根据论文已给出的信号公式计算单时隙信号部分。
% 输入：
%   H                    : [NM, K] 复合信道矩阵，第 k 列为 h_k(X)
%   signalConfig.N_RF    : RF 链路数
%   signalConfig.scheduledUsers : 已调度用户集合 S^(t)，用户索引向量
%   signalConfig.symbols : 与已调度用户对应的归一化符号向量
%   signalConfig.W       : [NM, |S|] 数字波束赋形矩阵，第 i 列对应 w_k^(t)
%   signalConfig.sigma2  : 噪声功率 sigma^2
%   signalConfig.noise   : [|S|,1] 可选，若给定则直接使用
%
% 输出 signal 结构体包含：
%   scheduledUsers, xRad, y, desired, interference, noise,
%   effectiveChannels, sinr, rate

    validateattributes(H, {'numeric'}, {'2d', 'nonempty'});

    requiredFields = {'N_RF', 'scheduledUsers', 'symbols', 'W', 'sigma2'};
    for idx = 1:numel(requiredFields)
        fieldName = requiredFields{idx};
        if ~isfield(signalConfig, fieldName)
            error('signalConfig 缺少必要字段: %s', fieldName);
        end
    end

    scheduledUsers = signalConfig.scheduledUsers(:);
    symbols = signalConfig.symbols(:);
    W = signalConfig.W;
    sigma2 = signalConfig.sigma2;
    numScheduled = numel(scheduledUsers);
    numUsers = size(H, 2);

    validateattributes(scheduledUsers, {'numeric'}, {'integer', 'positive', '<=', numUsers});
    validateattributes(symbols, {'numeric'}, {'column', 'numel', numScheduled});
    validateattributes(W, {'numeric'}, {'size', [size(H, 1), numScheduled]});
    validateattributes(sigma2, {'numeric'}, {'scalar', 'nonnegative'});

    if numScheduled > signalConfig.N_RF
        error('已调度用户数 |S^(t)|=%d 超过 N_RF=%d。', numScheduled, signalConfig.N_RF);
    end

    % 这里默认输入 symbols 已满足论文中的归一化条件 E[|s_k^(t)|^2]=1。
    xRad = W * symbols;

    if isfield(signalConfig, 'noise')
        noise = signalConfig.noise(:);
        validateattributes(noise, {'numeric'}, {'column', 'numel', numScheduled});
    else
        noiseSigma = sqrt(sigma2 / 2);
        noise = noiseSigma * (randn(numScheduled, 1) + 1j * randn(numScheduled, 1));
    end

    y = zeros(numScheduled, 1);
    desired = zeros(numScheduled, 1);
    interference = zeros(numScheduled, 1);
    sinr = zeros(numScheduled, 1);
    rate = zeros(numScheduled, 1);
    effectiveChannels = zeros(numScheduled, numScheduled);

    for idxUser = 1:numScheduled
        userId = scheduledUsers(idxUser);
        hk = H(:, userId);
        effectiveRow = hk' * W;
        effectiveChannels(idxUser, :) = effectiveRow;

        desired(idxUser) = effectiveRow(idxUser) * symbols(idxUser);
        interference(idxUser) = sum(effectiveRow .* symbols.') - desired(idxUser);
        y(idxUser) = desired(idxUser) + interference(idxUser) + noise(idxUser);

        desiredPower = abs(effectiveRow(idxUser)) ^ 2;
        interferencePower = sum(abs(effectiveRow) .^ 2) - desiredPower;
        sinr(idxUser) = desiredPower / (interferencePower + sigma2);
        rate(idxUser) = log2(1 + sinr(idxUser));
    end

    signal = struct();
    signal.scheduledUsers = scheduledUsers;
    signal.xRad = xRad;
    signal.y = y;
    signal.desired = desired;
    signal.interference = interference;
    signal.noise = noise;
    signal.effectiveChannels = effectiveChannels;
    signal.sinr = sinr;
    signal.rate = rate;
    signal.sumRate = sum(rate);
end
