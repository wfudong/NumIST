function coeff = istm_paper_coefficients(qfun, varargin)
%ISTM_PAPER_COEFFICIENTS Appendix A coefficients for the direct problem.
%
% Computes a_n(0), b_n(0), c_n(0), d_n(0) using the recurrent integration
% procedure in Appendix A of 152.pdf. The Jost values at rho=i/2 are
% obtained by the SPPS method from sppsmethod.pdf.

parser = inputParser;
parser.addParameter('HalfLength', 40, @(x) isnumeric(x) && isscalar(x) && x > 0);
parser.addParameter('NumGrid', 4001, @(x) isnumeric(x) && isscalar(x) && x >= 101);
parser.addParameter('NumCoefficients', 80, @(x) isnumeric(x) && isscalar(x) && x >= 1);
parser.addParameter('SppsTerms', 40, @(x) isnumeric(x) && isscalar(x) && x >= 1);
parser.addParameter('SppsSubintervalLength', 2.0, @(x) isnumeric(x) && isscalar(x) && x > 0);
parser.parse(varargin{:});

L = parser.Results.HalfLength;
numGrid = parser.Results.NumGrid;
if mod(numGrid, 2) == 0
    numGrid = numGrid + 1;
end
N = parser.Results.NumCoefficients;

x = linspace(-L, L, numGrid).';
zeroIndex = find(abs(x) == min(abs(x)), 1);
if abs(x(zeroIndex)) > 1e-12
    error('Internal grid does not contain x=0.');
end

q = qfun(x);
r = q + 0.25;

% Right Jost solution e(i/2,x), propagated from +L to -L.
xdesc = flipud(x);
rdesc = flipud(r);
eRightDesc = spps_propagate_second_order(xdesc, rdesc, 1.0, ...
    exp(-L / 2), -0.5 * exp(-L / 2), ...
    parser.Results.SppsTerms, parser.Results.SppsSubintervalLength);
e = flipud(eRightDesc.y);
ep = flipud(eRightDesc.dy);

% Left Jost solution g(i/2,x), propagated from -L to +L.
gLeft = spps_propagate_second_order(x, r, 1.0, ...
    exp(-L / 2), 0.5 * exp(-L / 2), ...
    parser.Results.SppsTerms, parser.Results.SppsSubintervalLength);
g = gLeft.y;
gp = gLeft.dy;

etaIntegral = istm_cumulative_from_zero(x, 1 ./ (e.^2));
eta = e .* etaIntegral;
etap = ep .* etaIntegral + 1 ./ e;

gCum = cumtrapz(x, 1 ./ (g.^2));
gCum0 = gCum(zeroIndex);
xiIntegral = gCum0 - gCum;
xi = g .* xiIntegral;
xip = gp .* xiIntegral - 1 ./ g;

qCum = cumtrapz(x, q);
qMinusAll = qCum;
qPlusAll = qCum(end) - qCum;
qMinus = qCum(zeroIndex);
qPlus = qCum(end) - qCum(zeroIndex);

a = zeros(numGrid, N + 1);
ap = zeros(numGrid, N + 1);
b = zeros(numGrid, N + 1);
bp = zeros(numGrid, N + 1);
d = zeros(numGrid, N + 1);
c = zeros(numGrid, N + 1);

a(:, 1) = e .* exp(x / 2) - 1;
ap(:, 1) = exp(x / 2) .* (ep + 0.5 * e);
b(:, 1) = g .* exp(-x / 2) - 1;
bp(:, 1) = exp(-x / 2) .* (gp - 0.5 * g);
d(:, 1) = ap(:, 1) - 0.5 * a(:, 1) + 0.5 * qPlusAll;
c(:, 1) = bp(:, 1) + 0.5 * b(:, 1) - 0.5 * qMinusAll;

J1 = zeros(numGrid, 1);
J2 = zeros(numGrid, 1);
I1 = zeros(numGrid, 1);
I2 = zeros(numGrid, 1);
J1p = zeros(numGrid, 1);
J2p = zeros(numGrid, 1);
I1p = zeros(numGrid, 1);
I2p = zeros(numGrid, 1);

pEeMinus = e .* exp(-x / 2);
pEeMinusPrime = exp(-x / 2) .* (ep - 0.5 * e);
pEtaMinus = eta .* exp(-x / 2);
pEtaMinusPrime = exp(-x / 2) .* (etap - 0.5 * eta);
pGPlus = g .* exp(x / 2);
pGPlusPrime = exp(x / 2) .* (gp + 0.5 * g);
pXiPlus = xi .* exp(x / 2);
pXiPlusPrime = exp(x / 2) .* (xip + 0.5 * xi);

pEtaForA = exp(x / 2) .* eta;
pEtaForAPrime = exp(x / 2) .* (etap + 0.5 * eta);
pEForA = exp(x / 2) .* e;
pEForAPrime = exp(x / 2) .* (ep + 0.5 * e);
pXiForB = exp(-x / 2) .* xi;
pXiForBPrime = exp(-x / 2) .* (xip - 0.5 * xi);
pGForB = exp(-x / 2) .* g;
pGForBPrime = exp(-x / 2) .* (gp - 0.5 * g);

for n = 2:(N + 1)
    prev = n - 1;
    intJ1 = integral_to_right(x, pEeMinusPrime .* a(:, prev));
    intJ2 = integral_to_right(x, pEtaMinusPrime .* a(:, prev));
    intI1 = integral_to_left(x, pGPlusPrime .* b(:, prev));
    intI2 = integral_to_left(x, pXiPlusPrime .* b(:, prev));

    J1 = J1 - pEeMinus .* a(:, prev) - intJ1;
    J2 = J2 - pEtaMinus .* a(:, prev) - intJ2;
    I1 = I1 + pGPlus .* b(:, prev) - intI1;
    I2 = I2 + pXiPlus .* b(:, prev) - intI2;

    J1p = J1p - pEeMinus .* ap(:, prev);
    J2p = J2p - pEtaMinus .* ap(:, prev);
    I1p = I1p + pGPlus .* bp(:, prev);
    I2p = I2p + pXiPlus .* bp(:, prev);

    a(:, n) = a(:, 1) - 2 * pEtaForA .* J1 + 2 * pEForA .* J2;
    ap(:, n) = ap(:, 1) - 2 * pEtaForAPrime .* J1 - 2 * pEtaForA .* J1p ...
        + 2 * pEForAPrime .* J2 + 2 * pEForA .* J2p;

    b(:, n) = b(:, 1) + 2 * pXiForB .* I1 - 2 * pGForB .* I2;
    bp(:, n) = bp(:, 1) + 2 * pXiForBPrime .* I1 + 2 * pXiForB .* I1p ...
        - 2 * pGForBPrime .* I2 - 2 * pGForB .* I2p;

    d(:, n) = d(:, prev) + ap(:, n) - ap(:, prev) - 0.5 * (a(:, n) + a(:, prev));
    c(:, n) = c(:, prev) + bp(:, n) - bp(:, prev) + 0.5 * (b(:, n) + b(:, prev));
end

coeff = struct();
coeff.x = x;
coeff.q = q;
coeff.eHalf = e;
coeff.gHalf = g;
coeff.aAll = a;
coeff.bAll = b;
coeff.cAll = c;
coeff.dAll = d;
coeff.a0 = a(zeroIndex, :).';
coeff.b0 = b(zeroIndex, :).';
coeff.c0 = c(zeroIndex, :).';
coeff.d0 = d(zeroIndex, :).';
coeff.qPlus = qPlus;
coeff.qMinus = qMinus;
coeff.zeroIndex = zeroIndex;
coeff.options = parser.Results;
end

function out = integral_to_right(x, f)
cum = cumtrapz(x, f);
out = cum(end) - cum;
end

function out = integral_to_left(x, f)
out = cumtrapz(x, f);
end
