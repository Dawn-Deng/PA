function [symbols, bits] = signal_generate_symbols(params)
%SIGNAL_GENERATE_SYMBOLS 生成多用户下行基带符号。
% 输出：
%   symbols : [K,1] 复符号向量
%   bits    : [K,bitsPerSymbol] 对应比特

    if ~isfield(params, 'modulation')
        params.modulation = 'qpsk';
    end

    switch lower(params.modulation)
        case 'qpsk'
            bits = randi([0, 1], params.K, 2);
            iPart = 1 - 2 * bits(:, 1);
            qPart = 1 - 2 * bits(:, 2);
            symbols = (iPart + 1j * qPart) / sqrt(2);

        otherwise
            error('当前未实现调制方式: %s', params.modulation);
    end
end
