function sol = spps_ivp_second_order(x, r, lambda, y0, dy0, numTerms)
%SPPS_IVP_SECOND_ORDER Solve y'' = lambda*r(x)*y by SPPS.
%
% This is the p=1, q=0, u0=1 specialization of the SPPS construction in
% sppsmethod.pdf. x may be increasing or decreasing; cumulative trapezoidal
% integration gives the signed integrals from x(1).
%
% The returned fundamental solutions satisfy
%   u1(x(1)) = 1, u1'(x(1)) = 0,
%   u2(x(1)) = 0, u2'(x(1)) = 1.

x = x(:);
r = r(:);
if numel(x) ~= numel(r)
    error('x and r must have the same length.');
end
if numel(x) < 2
    error('x must contain at least two points.');
end
if numTerms < 1
    error('numTerms must be positive.');
end

n = numel(x);
maxPower = 2 * numTerms + 1;
xtilde = zeros(n, maxPower + 1);
xplain = zeros(n, maxPower + 1);
xtilde(:, 1) = 1;
xplain(:, 1) = 1;

for power = 1:maxPower
    prevTilde = xtilde(:, power);
    prevPlain = xplain(:, power);
    if mod(power, 2) == 1
        dTilde = prevTilde .* r;
        dPlain = prevPlain;
    else
        dTilde = prevTilde;
        dPlain = prevPlain .* r;
    end
    xtilde(:, power + 1) = cumtrapz(x, dTilde);
    xplain(:, power + 1) = cumtrapz(x, dPlain);
end

u1 = zeros(n, 1);
du1 = zeros(n, 1);
u2 = zeros(n, 1);
du2 = zeros(n, 1);
lambdaPower = 1;

for k = 0:numTerms
    idxEven = 2 * k + 1;
    idxOdd = 2 * k + 2;

    u1 = u1 + lambdaPower * xtilde(:, idxEven);
    if k > 0
        prevIdx = idxEven - 1;
        du1 = du1 + lambdaPower * derivative_formal_power(2 * k, xtilde(:, prevIdx), r, true);
    end

    u2 = u2 + lambdaPower * xplain(:, idxOdd);
    du2 = du2 + lambdaPower * derivative_formal_power(2 * k + 1, xplain(:, idxOdd - 1), r, false);

    lambdaPower = lambdaPower * lambda;
end

y = y0 * u1 + dy0 * u2;
dy = y0 * du1 + dy0 * du2;
sol = struct('x', x, 'y', y, 'dy', dy, 'u1', u1, 'du1', du1, 'u2', u2, 'du2', du2);
end

function dval = derivative_formal_power(power, previousPower, r, isTilde)
if isTilde
    if mod(power, 2) == 1
        dval = previousPower .* r;
    else
        dval = previousPower;
    end
else
    if mod(power, 2) == 1
        dval = previousPower;
    else
        dval = previousPower .* r;
    end
end
end
