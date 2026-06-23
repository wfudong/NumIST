function evolved = istm_evolve_scattering(data, t)
%ISTM_EVOLVE_SCATTERING Apply the KdV scattering-data evolution law.
evolved = data;
if strcmp(data.side, 'right')
    sgn = 1;
else
    sgn = -1;
end
evolved.alpha = data.alpha .* exp(sgn * 8 * data.tau.^3 * t);
evolved.reflectionValues = data.reflectionValues .* exp(sgn * 8i * data.rho.^3 * t);
end
