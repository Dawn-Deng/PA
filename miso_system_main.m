function simulation = miso_system_main(config)
%MISO_SYSTEM_MAIN 主函数：并列调用“信道部分”和“信号部分”两个子模块。
% 输入：
%   config.channel : 传给 miso_waveguide_channel_model 的配置
%   config.signal  : 传给 miso_signal_model 的配置
%
% 输出 simulation 结构体包含：
%   simulation.channel : 信道子模块输出
%   simulation.signal  : 信号子模块输出

    if nargin < 1
        config = struct();
    end

    if ~isfield(config, 'channel')
        config.channel = struct();
    end
    if ~isfield(config, 'signal')
        error('主函数需要提供 config.signal，用于信号部分建模。');
    end

    channel = miso_waveguide_channel_model(config.channel);
    signal = miso_signal_model(channel.H, config.signal);

    simulation = struct();
    simulation.channel = channel;
    simulation.signal = signal;
end
