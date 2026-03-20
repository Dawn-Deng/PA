function results = miso_system_main()
%MISO_SYSTEM_MAIN 主函数：统一调用信道建模与信号建模模块

    clc;

    params = defaultSystemParameters();

    channelState = miso_waveguide_channel_model(params);
    signalState = miso_signal_model(channelState.H, params);

    disp('=== 信道建模结果 ===');
    disp('用户全局坐标 q_k = [x_k, y_k, 0]^T：');
    disp(channelState.users);

    disp('波导馈电点 [x_n^(W), 0, d]^T：');
    disp(channelState.waveguideFeedPoints);

    disp('第 1 条波导上所有 PA 的位置 p_{1,m}：');
    disp(squeeze(channelState.paPositions(:, :, 1)).');

    disp('每个 PA 的指向目标用户索引：');
    disp(channelState.steeringTargets);

    disp('第 1 条波导的导波信道 g_1(x_1)：');
    disp(channelState.guidedWaveExample);

    disp('用户 1 相对于第 1 条波导上各 PA 的局部坐标 [x~, y~, z~]：');
    disp(channelState.localCoordsExample.');

    disp('用户 1 与第 1 条波导之间的自由空间信道 \\tilde{h}_{1,1}(x_1)：');
    disp(channelState.freeSpaceExample);

    disp('所有用户的总复合信道矩阵 H（每列对应一个用户 h_k）：');
    disp(channelState.H);

    disp('=== 信号建模结果 ===');
    fprintf('调度用户集合 S^(t) = [%s]\n', num2str(signalState.scheduledUsers.'));
    fprintf('发射信号 x_rad^(t) 的维度：%d x %d\n', size(signalState.xRad, 1), size(signalState.xRad, 2));
    disp('各已调度用户的归一化数据符号 s_k^(t)：');
    disp(signalState.symbols.');

    disp('数字波束赋形矩阵 W = [w_1, ..., w_|S|]：');
    disp(signalState.W);

    disp('每个已调度用户的接收信号、SINR 与频谱效率：');
    disp(signalState.userMetricsTable);

    results = struct();
    results.params = params;
    results.channelState = channelState;
    results.signalState = signalState;
end

function params = defaultSystemParameters()
%DEFAULTSYSTEMPARAMETERS 设置论文建模示例所需的默认参数

    params.N = 2;                 % 波导数量
    params.M = 3;                 % 每条波导的 PA 数量
    params.K = 4;                 % 用户总数
    params.NRF = 2;               % RF 链路数，对应最多并行数据流数
    params.Dx = 10;               % 波导在 x 方向覆盖宽度
    params.d = 4;                 % 波导部署高度
    params.randomSeed = 20260320; % 随机种子，保证结果可复现

    params.userRegionX = [0, params.Dx];
    params.userRegionY = [6, 20];
    params.paYOffsetRange = [1, 10];
    params.paPlacementMode = 'uniform'; % 可选：'uniform' / 'random'

    params.lambda = 0.01;         % 载波波长（m）
    params.nEff = 1.6;            % 波导有效折射率 n_eff
    params.alphaW = 0.08;         % 波导衰减系数 alpha_W（Np/m）
    params.alphaL = 0.96;         % LoS 存在系数/米
    params.nMedium = 1.0;         % 自由空间/介质折射率 n
    params.a = 0.45;              % 模式尺寸参数 a
    params.b = 0.30;              % 模式尺寸参数 b
    params.v = 1.1;               % 经验修正系数 v

    params.totalTransmitPower = 1.0; % 总发射功率
    params.sigma2 = 1e-3;            % 噪声功率
    params.slotIndex = 1;            % 传输时隙编号
    params.schedulingMode = 'strongest'; % 可选：'strongest' / 'first'
    params.precoderMode = 'mrt';        % 当前使用 MRT 数字预编码
end
