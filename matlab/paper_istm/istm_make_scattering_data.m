function data = istm_make_scattering_data(side, tau, alpha, rho, reflectionValues)
%ISTM_MAKE_SCATTERING_DATA Container for one side of scattering data.
side = validatestring(side, {'right', 'left'});
data = struct();
data.side = side;
data.tau = tau(:).';
data.alpha = alpha(:).';
data.rho = rho(:).';
data.reflectionValues = reflectionValues(:).';
if numel(data.tau) ~= numel(data.alpha)
    error('tau and alpha must have the same length.');
end
end
