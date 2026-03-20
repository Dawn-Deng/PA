function [yRx, metadata] = signal_receive_users(H, xTx, params)
%SIGNAL_RECEIVE_USERS 根据实际发射向量计算每个用户的接收信号。
% 对第 k 个用户：y_k = h_k^H x + n_k。

    if ~isfield(params, 'noisePower')
        params.noisePower = 1e-3;
    end

    noiseSigma = sqrt(params.noisePower / 2);
    noise = noiseSigma * (randn(size(H, 2), 1) + 1j * randn(size(H, 2), 1));

    yRx = H' * xTx + noise;

    metadata = struct();
    metadata.effectiveTxChannel = H';
    metadata.noise = noise;
    metadata.noisePower = params.noisePower;
end
