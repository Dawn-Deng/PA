%% 多用户 MISO 波导系统：信号部分分模块建模示例
% 运行方式：直接执行本脚本。
% 你可以把论文中的每个小部分替换到对应子函数中，而主流程保持不变。

clear; clc;

config = struct();
config.signalParams = struct( ...
    'modulation', 'qpsk', ...
    'precoderType', 'mrt', ...
    'txPower', 1, ...
    'noisePower', 1e-3);

simulation = miso_signal_model_main(config);

fprintf('=== 信号建模结果概览 ===\n');
fprintf('用户数 K = %d, 波导数 N = %d, 每条波导 PA 数 M = %d\n', ...
    simulation.params.K, simulation.params.N, simulation.params.M);

fprintf('\n生成的用户符号 s:\n');
disp(simulation.symbols.');

fprintf('预编码矩阵 W 的尺寸: %d x %d\n', size(simulation.precoder, 1), size(simulation.precoder, 2));
fprintf('发射向量 x 的尺寸: %d x %d\n', size(simulation.txVector, 1), size(simulation.txVector, 2));

fprintf('\n各用户接收信号 y:\n');
disp(simulation.rxSignal.');

fprintf('\n各用户 SINR:\n');
disp(simulation.metrics.sinr.');

fprintf('系统总速率 sum-rate = %.4f bit/s/Hz\n', simulation.metrics.sumRate);
