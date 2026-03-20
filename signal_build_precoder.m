function [W, metadata] = signal_build_precoder(H, params)
%SIGNAL_BUILD_PRECODER 根据信道矩阵构建预编码器。
% 输入 H 大小为 [NM, K]。
% 输出 W 大小为 [NM, K]，第 k 列对应用户 k 的预编码向量。

    if ~isfield(params, 'precoderType')
        params.precoderType = 'mrt';
    end

    switch lower(params.precoderType)
        case 'mrt'
            W = conj(H);

        case 'zf'
            gram = H' * H;
            regularizer = 1e-9 * eye(size(gram));
            W = H / (gram + regularizer);

        otherwise
            error('当前未实现预编码方式: %s', params.precoderType);
    end

    columnNorms = sqrt(sum(abs(W) .^ 2, 1));
    columnNorms(columnNorms == 0) = 1;
    W = W ./ columnNorms;

    metadata = struct();
    metadata.columnNorms = columnNorms;
    metadata.precoderType = params.precoderType;
end
