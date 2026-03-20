function simulation = miso_signal_model_main(config)
%MISO_SIGNAL_MODEL_MAIN 信号部分总控主函数。
% 建议把论文中的“信号部分”拆成以下可替换模块：
%   1) 信道/几何建模            -> miso_waveguide_channel_model
%   2) 多用户符号生成            -> signal_generate_symbols
%   3) 预编码/波束赋形           -> signal_build_precoder
%   4) 发射向量形成              -> signal_form_waveguide_input
%   5) 用户接收信号计算          -> signal_receive_users
%   6) 性能指标评估              -> signal_evaluate_metrics
%
% 这样你后续替换论文中的具体公式时，只需逐个修改子模块。

    if nargin < 1
        config = struct();
    end

    channelModel = miso_waveguide_channel_model(config);
    params = channelModel.params;

    params = applySignalDefaults(params, config);

    [symbols, bits] = signal_generate_symbols(params);
    [W, precoderInfo] = signal_build_precoder(channelModel.H, params);
    [xTx, txInfo] = signal_form_waveguide_input(W, symbols, params);
    [yRx, rxInfo] = signal_receive_users(channelModel.H, xTx, params);
    metrics = signal_evaluate_metrics(channelModel.H, W, symbols, txInfo.scale, params);

    simulation = struct();
    simulation.channel = channelModel;
    simulation.params = params;
    simulation.bits = bits;
    simulation.symbols = symbols;
    simulation.precoder = W;
    simulation.precoderInfo = precoderInfo;
    simulation.txVector = xTx;
    simulation.txInfo = txInfo;
    simulation.rxSignal = yRx;
    simulation.rxInfo = rxInfo;
    simulation.metrics = metrics;
end

function params = applySignalDefaults(params, config)
%APPLYSIGNALDEFAULTS 设置并覆盖信号处理部分参数
    if isfield(config, 'signalParams')
        override = config.signalParams;
        fields = fieldnames(override);
        for idx = 1:numel(fields)
            params.(fields{idx}) = override.(fields{idx});
        end
    end

    if ~isfield(params, 'modulation')
        params.modulation = 'qpsk';
    end
    if ~isfield(params, 'precoderType')
        params.precoderType = 'mrt';
    end
    if ~isfield(params, 'txPower')
        params.txPower = 1;
    end
    if ~isfield(params, 'noisePower')
        params.noisePower = 1e-3;
    end
end
