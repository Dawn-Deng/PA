function validateParameters(params)
%VALIDATEPARAMETERS 工程化增强：统一参数合法性检查

    if params.Dy < (params.M - 1) * params.deltaMin
        error('Invalid geometry: Dy must satisfy Dy >= (M-1)*deltaMin.');
    end
    if params.KServ > min(params.NRF, params.Kmax)
        error('Invalid scheduling: KServ must be <= min(NRF, Kmax).');
    end
    if ~(params.Pmax > 0)
        error('Invalid power: Pmax must be > 0.');
    end
    if ~(params.sigma2 > 0)
        error('Invalid noise: sigma2 must be > 0.');
    end

    assertOpenUnitInterval(params.betaTheta, 'betaTheta');
    assertOpenUnitInterval(params.betaPhi, 'betaPhi');
    assertOpenUnitInterval(params.positionLineSearchBeta, 'positionLineSearchBeta');

    if params.positionMemory < 1 || mod(params.positionMemory,1) ~= 0
        error('positionMemory must be an integer >= 1.');
    end

    critical = [params.N, params.M, params.K, params.NRF, params.Kmax, params.KServ, ...
        params.Dx, params.Dy, params.d, params.deltaMin, params.Pmax, params.sigma2];
    if any(~isfinite(critical))
        error('Parameters contain NaN/Inf in critical fields.');
    end
end

function assertOpenUnitInterval(x, name)
    if ~(isfinite(x) && x > 0 && x < 1)
        error('%s must be in (0,1).', name);
    end
end
