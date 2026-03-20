function [xTx, metadata] = signal_form_waveguide_input(W, symbols, params)
%SIGNAL_FORM_WAVEGUIDE_INPUT 形成馈入所有波导-PA 通道的发射向量。
% 模型：x = sqrt(Pt) * W * s，其中 E[||W s||^2] 由归一化控制。

    if ~isfield(params, 'txPower')
        params.txPower = 1;
    end

    xLinear = W * symbols;
    powerBeforeScaling = real(xLinear' * xLinear);

    if powerBeforeScaling <= 0
        scale = 0;
    else
        scale = sqrt(params.txPower / powerBeforeScaling);
    end

    xTx = scale * xLinear;

    metadata = struct();
    metadata.powerBeforeScaling = powerBeforeScaling;
    metadata.scale = scale;
    metadata.txPower = params.txPower;
end
