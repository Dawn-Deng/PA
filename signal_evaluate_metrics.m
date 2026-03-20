function metrics = signal_evaluate_metrics(H, W, symbols, txScale, params)
%SIGNAL_EVALUATE_METRICS 评估各用户的有用信号、干扰、SINR 与速率。

    if ~isfield(params, 'noisePower')
        params.noisePower = 1e-3;
    end

    effectiveChannel = txScale * (H' * W);
    K = size(effectiveChannel, 1);

    usefulPower = zeros(K, 1);
    interferencePower = zeros(K, 1);
    sinr = zeros(K, 1);
    rate = zeros(K, 1);

    symbolPower = abs(symbols) .^ 2;

    for k = 1:K
        rowPower = abs(effectiveChannel(k, :)) .^ 2;
        usefulPower(k) = rowPower(k) * symbolPower(k);
        interferencePower(k) = sum(rowPower .* symbolPower.') - usefulPower(k);
        sinr(k) = usefulPower(k) / (interferencePower(k) + params.noisePower);
        rate(k) = log2(1 + sinr(k));
    end

    metrics = struct();
    metrics.effectiveChannel = effectiveChannel;
    metrics.usefulPower = usefulPower;
    metrics.interferencePower = interferencePower;
    metrics.sinr = sinr;
    metrics.rate = rate;
    metrics.sumRate = sum(rate);
end
