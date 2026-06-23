function [theta, z, rho, sPlus, sMinus] = istm_reflection_grid_from_direct(scattering, numTheta)
%ISTM_REFLECTION_GRID_FROM_DIRECT Evaluate s^+ and s^- on the 9.4 grid.
if nargin < 2
    numTheta = 10000;
end
dtheta = 2 * pi / numTheta;
theta = -pi + ((0:numTheta-1) + 0.5) * dtheta;
z = exp(1i * theta);
rho = real(1i * (1 - z) ./ (2 * (1 + z)));
sPlus = scattering.reflectionPlus(rho);
sMinus = scattering.reflectionMinus(rho);
end
