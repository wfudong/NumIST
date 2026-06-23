function result = istm_recover_potential(rightData, leftData, xgrid, varargin)
%ISTM_RECOVER_POTENTIAL Recover q from J+ and J- by (9.2)-(9.3).
%
% Uses the numerically stable branch on each half-line. With the scattering
% data orientation produced by istm_paper_direct_scattering, the left-data
% system is stable for x < switch and the right-data system for x >= switch.

parser = inputParser;
parser.addParameter('Ns', 5, @(x) isnumeric(x) && isscalar(x) && x >= 0);
parser.addParameter('Switch', 0, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('DerivativeMethod', 'spline', @(s) ischar(s) || isstring(s));
parser.addParameter('ReflectionMatrixFactor', 1, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('ReflectionRhsFactor', 1, @(x) isnumeric(x) && isscalar(x));
parser.parse(varargin{:});

xgrid = xgrid(:);
rightBranch = recover_branch(rightData, xgrid, parser.Results.Ns, parser.Results.DerivativeMethod, ...
    parser.Results.ReflectionMatrixFactor, parser.Results.ReflectionRhsFactor);
leftBranch = recover_branch(leftData, xgrid, parser.Results.Ns, parser.Results.DerivativeMethod, ...
    parser.Results.ReflectionMatrixFactor, parser.Results.ReflectionRhsFactor);

q = leftBranch.q;
useRight = xgrid >= parser.Results.Switch;
q(useRight) = rightBranch.q(useRight);

result = struct();
result.x = xgrid;
result.q = real(q);
result.qRight = real(rightBranch.q);
result.qLeft = real(leftBranch.q);
result.a0 = rightBranch.c0;
result.b0 = leftBranch.c0;
end

function branch = recover_branch(data, xgrid, Ns, derivativeMethod, reflectionMatrixFactor, reflectionRhsFactor)
c0 = zeros(size(xgrid));
for j = 1:numel(xgrid)
    coeffs = finite_section_coefficients(data, xgrid(j), Ns, reflectionMatrixFactor, reflectionRhsFactor);
    c0(j) = coeffs(1);
end

switch lower(string(derivativeMethod))
    case "spline"
        [dc0, d2c0] = spline_derivatives(xgrid, c0);
    case "finite-difference"
        [dc0, d2c0] = finite_difference_derivatives(xgrid, c0);
    otherwise
        error('Unknown derivative method "%s".', derivativeMethod);
end
if strcmp(data.side, 'right')
    q = (d2c0 - dc0) ./ (1 + c0);
else
    q = (d2c0 + dc0) ./ (1 + c0);
end

function [dy, d2y] = finite_difference_derivatives(x, y)
x = x(:);
y = y(:);
n = numel(x);
if n < 5
    error('Need at least five points for finite-difference derivatives.');
end
h = x(2) - x(1);
dy = zeros(n, 1);
d2y = zeros(n, 1);

dy(1) = (-3*y(1) + 4*y(2) - y(3)) / (2*h);
dy(n) = (3*y(n) - 4*y(n-1) + y(n-2)) / (2*h);
d2y(1) = (2*y(1) - 5*y(2) + 4*y(3) - y(4)) / h^2;
d2y(n) = (2*y(n) - 5*y(n-1) + 4*y(n-2) - y(n-3)) / h^2;

dy(2) = (y(3) - y(1)) / (2*h);
dy(n-1) = (y(n) - y(n-2)) / (2*h);
d2y(2) = (y(3) - 2*y(2) + y(1)) / h^2;
d2y(n-1) = (y(n) - 2*y(n-1) + y(n-2)) / h^2;

for j = 3:(n-2)
    dy(j) = (-y(j+2) + 8*y(j+1) - 8*y(j-1) + y(j-2)) / (12*h);
    d2y(j) = (-y(j+2) + 16*y(j+1) - 30*y(j) + 16*y(j-1) - y(j-2)) / (12*h^2);
end
end

branch = struct('c0', c0, 'dc0', dc0, 'd2c0', d2c0, 'q', q);
end

function coeffs = finite_section_coefficients(data, x, Ns, reflectionMatrixFactor, reflectionRhsFactor)
dim = Ns + 1;
matrix = eye(dim);
rhs = zeros(dim, 1);

reflMatrix = zeros(2 * Ns + 1, 1);
for p = 0:(2 * Ns)
    reflMatrix(p + 1) = reflection_integral(data, x, p, p + 2);
end
reflRhs = zeros(Ns + 1, 1);
for m = 0:Ns
    reflRhs(m + 1) = reflection_integral(data, x, m, m + 1);
end

for m = 0:Ns
    rhs(m + 1) = sign_power(m + 1) * (discrete_sum(data, x, m, 1) ...
        + reflectionRhsFactor * reflRhs(m + 1) / (2 * pi));
    for n = 0:Ns
        p = m + n;
        matrix(m + 1, n + 1) = matrix(m + 1, n + 1) ...
            + sign_power(p) * (discrete_sum(data, x, p, 2) ...
            + reflectionMatrixFactor * reflMatrix(p + 1) / (2 * pi));
    end
end

coeffs = matrix \ rhs;
end

function value = discrete_sum(data, x, power, denomShift)
if isempty(data.tau)
    value = 0;
    return;
end
if strcmp(data.side, 'right')
    expo = exp(-2 * data.tau * x);
else
    expo = exp(2 * data.tau * x);
end
value = sum(data.alpha .* expo .* (0.5 - data.tau).^power ./ (0.5 + data.tau).^(power + denomShift));
end

function value = reflection_integral(data, x, m, n)
% Implements (9.4) for I(x)=int s(rho) exp(2i rho x) ... d rho.
if isempty(data.reflectionValues)
    value = 0;
    return;
end
if strcmp(data.side, 'right')
    xeff = x;
else
    xeff = -x;
end
rho = data.rho;
z = (0.5 + 1i * rho) ./ (0.5 - 1i * rho);
numTheta = numel(rho);
dtheta = 2 * pi / numTheta;
integrand = data.reflectionValues ...
    .* exp((z - 1) ./ (z + 1) * xeff) ...
    .* z.^(m + 1) ...
    .* (z + 1).^(n - m - 2);
value = dtheta * sum(integrand);
end

function s = sign_power(k)
if mod(k, 2) == 0
    s = 1;
else
    s = -1;
end
end

function [dy, d2y] = spline_derivatives(x, y)
pp = spline(x, y);
dpp = pp_derivative(pp, 1);
d2pp = pp_derivative(pp, 2);
dy = ppval(dpp, x);
d2y = ppval(d2pp, x);
dy = dy(:);
d2y = d2y(:);
end

function dpp = pp_derivative(pp, order)
dpp = pp;
for k = 1:order
    coefs = dpp.coefs;
    degree = size(coefs, 2) - 1;
    if degree <= 0
        dpp.coefs = zeros(size(coefs, 1), 1);
        dpp.order = 1;
        continue;
    end
    powers = degree:-1:1;
    dpp.coefs = coefs(:, 1:end-1) .* powers;
    dpp.order = dpp.order - 1;
end
end
