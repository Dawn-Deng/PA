function sanityCheckState(state, params)
%SANITYCHECKSTATE 工程化增强：统一状态数值健康检查

    checks = {};
    if isfield(state, 'X'), checks{end+1} = state.X; end %#ok<AGROW>
    if isfield(state, 'theta'), checks{end+1} = state.theta; end %#ok<AGROW>
    if isfield(state, 'phi'), checks{end+1} = state.phi; end %#ok<AGROW>
    if isfield(state, 'H'), checks{end+1} = state.H; end %#ok<AGROW>
    if isfield(state, 'channelMatrix'), checks{end+1} = state.channelMatrix; end %#ok<AGROW>
    if isfield(state, 'W') && ~isempty(state.W), checks{end+1} = state.W; end %#ok<AGROW>

    for i = 1:numel(checks)
        arr = checks{i};
        if any(~isfinite(real(arr(:)))) || any(~isfinite(imag(arr(:))))
            error('State contains NaN/Inf in numeric arrays.');
        end
    end

    if isfield(state, 'X')
        if any(state.X(:) < -1e-9) || any(state.X(:) > params.Dy + 1e-9)
            error('State.X is out of bounds [0, Dy].');
        end
    end
end
