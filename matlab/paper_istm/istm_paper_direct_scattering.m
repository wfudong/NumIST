function scattering = istm_paper_direct_scattering(qfun, varargin)
%ISTM_PAPER_DIRECT_SCATTERING Direct scattering by the method of 152.pdf.

parser = inputParser;
parser.addParameter('HalfLength', 40, @(x) isnumeric(x) && isscalar(x) && x > 0);
parser.addParameter('NumGrid', 4001, @(x) isnumeric(x) && isscalar(x) && x >= 101);
parser.addParameter('NumCoefficients', 80, @(x) isnumeric(x) && isscalar(x) && x >= 1);
parser.addParameter('SppsTerms', 40, @(x) isnumeric(x) && isscalar(x) && x >= 1);
parser.addParameter('SppsSubintervalLength', 2.0, @(x) isnumeric(x) && isscalar(x) && x > 0);
parser.addParameter('RootGridSize', 4000, @(x) isnumeric(x) && isscalar(x) && x >= 100);
parser.addParameter('RootTolerance', 1e-10, @(x) isnumeric(x) && isscalar(x) && x > 0);
parser.parse(varargin{:});

coeff = istm_paper_coefficients(qfun, ...
    'HalfLength', parser.Results.HalfLength, ...
    'NumGrid', parser.Results.NumGrid, ...
    'NumCoefficients', parser.Results.NumCoefficients, ...
    'SppsTerms', parser.Results.SppsTerms, ...
    'SppsSubintervalLength', parser.Results.SppsSubintervalLength);

[zRoots, tau, lambda] = find_discrete_roots(coeff, parser.Results.RootGridSize, parser.Results.RootTolerance);
[alphaPlus, alphaMinus, dRatio, aPrime] = norming_constants(coeff, zRoots);

scattering = struct();
scattering.coeff = coeff;
scattering.zRoots = zRoots(:);
scattering.tau = tau(:);
scattering.lambda = lambda(:);
scattering.alphaPlus = real(alphaPlus(:));
scattering.alphaMinus = real(alphaMinus(:));
scattering.dRatio = dRatio(:);
scattering.aPrime = aPrime(:);
scattering.options = parser.Results;
scattering.reflectionPlus = @(rho) reflection_coeff(coeff, rho, +1);
scattering.reflectionMinus = @(rho) reflection_coeff(coeff, rho, -1);
end

function [zRoots, tau, lambda] = find_discrete_roots(coeff, rootGridSize, tol)
epsEdge = 1e-6;
zgrid = linspace(-1 + epsEdge, 1 - epsEdge, rootGridSize);
phi = real(istm_eval_egEG(zgrid, coeff).Phi);

zRoots = [];
for j = 1:(numel(zgrid) - 1)
    if ~isfinite(phi(j)) || ~isfinite(phi(j + 1))
        continue;
    end
    if abs(phi(j)) < tol
        zRoots(end + 1) = zgrid(j); %#ok<AGROW>
    elseif phi(j) * phi(j + 1) < 0
        root = fzero(@(z) real(istm_eval_egEG(z, coeff).Phi), [zgrid(j), zgrid(j + 1)]);
        zRoots(end + 1) = root; %#ok<AGROW>
    end
end

if isempty(zRoots)
    tau = [];
    lambda = [];
    return;
end

zRoots = unique_roots(sort(zRoots), 1e-6);
tau = (1 - zRoots) ./ (2 * (1 + zRoots));
positive = tau > 0 & isfinite(tau);
zRoots = zRoots(positive);
tau = tau(positive);
lambda = -tau.^2;
end

function rootsOut = unique_roots(rootsIn, tolerance)
rootsOut = [];
for k = 1:numel(rootsIn)
    if isempty(rootsOut) || abs(rootsIn(k) - rootsOut(end)) > tolerance
        rootsOut(end + 1) = rootsIn(k); %#ok<AGROW>
    end
end
end

function [alphaPlus, alphaMinus, dRatio, aPrime] = norming_constants(coeff, zRoots)
alphaPlus = zeros(size(zRoots));
alphaMinus = zeros(size(zRoots));
dRatio = zeros(size(zRoots));
aPrime = zeros(size(zRoots));
for k = 1:numel(zRoots)
    z = zRoots(k);
    evals = istm_eval_egEG(z, coeff);
    dRatio(k) = evals.g / evals.e;
    aPrime(k) = 2i * ((z + 1) / (z - 1))^2 * evals.Phi ...
        - 1i * (z + 1)^3 / (z - 1) * evals.Phiz;
    denom = 1i * aPrime(k);
    alphaPlus(k) = dRatio(k) / denom;
    alphaMinus(k) = (1 / dRatio(k)) / denom;
end
end

function s = reflection_coeff(coeff, rho, side)
rhoShape = size(rho);
rho = rho(:).';
z = (0.5 + 1i * rho) ./ (0.5 - 1i * rho);
base = istm_eval_egEG(z, coeff);
zc = conj(z);
conjSide = istm_eval_egEG(zc, coeff);
den = base.e .* base.G - base.E .* base.g;
if side > 0
    num = conjSide.e .* base.G - conjSide.E .* base.g;
else
    num = base.e .* conjSide.G - base.E .* conjSide.g;
end
% The finite-section formulas (6.4)-(6.7) use the reflection coefficient
% with the opposite Wronskian orientation from the intermediate numerator
% above.
s = reshape(-num ./ den, rhoShape);
end
