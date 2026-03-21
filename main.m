%% 主入口：先进行信道建模，再进行信号建模
clear; clc;

[channelState, params] = miso_channel_model();
signalState = miso_signal_model(channelState, params);

%% =========================
%  1. 信道部分结果展示
%  =========================
disp('===== 信道部分建模结果 =====');
disp('用户全局坐标 q_k = [x_k, y_k, 0]^T：');
disp(channelState.users);

disp('波导馈电点 [x_n^(W), 0, d]^T：');
disp(channelState.waveguideFeedPoints);

disp('第 1 条波导上所有 PA 的位置 p_{1,m}：');
disp(squeeze(channelState.paPositions(:, :, 1)).');

disp('每个 PA 的指向目标用户索引：');
disp(channelState.steeringTargets);

disp('第 1 条波导的导波信道 g_1(x_1)：');
disp(channelState.g1);

disp('用户 1 相对于第 1 条波导上各 PA 的局部坐标 [x~, y~, z~]：');
disp(channelState.localCoordsExample.');

disp('用户 1 与第 1 条波导之间的自由空间信道 \tilde{h}_{1,1}(x_1)：');
disp(channelState.hFreeExample);

disp('所有用户的总复合信道矩阵 H（每列对应一个用户 h_k）：');
disp(channelState.H);

fprintf('用户 1 到波导 1 的第 1 个 PA 的局部坐标为:\n');
disp(channelState.localCoordExample.');

%% =========================
%  2. 信号部分结果展示
%  =========================
disp('===== 信号部分建模结果 =====');
disp('用户行信道矩阵 channelMatrix，其中第 k 行为 h_k^H(X)：');
disp(signalState.channelMatrix);

disp('本时隙调度用户集合 S^(t)：');
disp(signalState.scheduledUsers);

disp('归一化发送符号向量 s^(t)：');
disp(signalState.symbols);

disp(['采用的数字波束形成方法: ', signalState.beamformerInfo.method]);
disp('数字波束形成矩阵 W：');
disp(signalState.W);

disp('发射信号 x_rad^(t)：');
disp(signalState.slotResult.xRad);

disp('各调度用户的接收信号 y_k^(t)：');
disp(signalState.slotResult.y);

disp('各调度用户的瞬时 SINR：');
disp(signalState.slotResult.sinr);

disp('各调度用户的瞬时频谱效率 R_k^(t) = log2(1 + SINR)：');
disp(signalState.slotResult.rate);

disp('每个调度用户的分解结果（期望信号 / 多用户干扰 / 噪声）：');
disp(struct2table(signalState.slotResult.userMetrics));
