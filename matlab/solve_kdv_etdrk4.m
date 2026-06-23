function sol = solve_kdv_etdrk4(q0, xgrid, times, dt, varargin)
%SOLVE_KDV_ETDRK4 Fourier ETDRK4 solver for u_t - 6*u*u_x + u_xxx = 0.
%
% q0 and xgrid must be vectors on one periodic cell. times must be sorted
% and begin with zero. The solver uses steps no larger than dt and takes a
% shorter final step when a requested time is not an exact multiple of dt.

parser = inputParser;
parser.addParameter('Dealias', true, @(x) islogical(x) && isscalar(x));
parser.addParameter('ContourPoints', 32, @(x) isnumeric(x) && isscalar(x) && x > 0);
parser.parse(varargin{:});

xgrid = xgrid(:);
q0 = q0(:);
times = times(:).';
n = numel(xgrid);
if numel(q0) ~= n
    error('q0 and xgrid must have the same length.');
end
if n < 8 || mod(n, 2) ~= 0
    error('The grid length must be even and at least 8.');
end
if any(diff(times) < 0) || abs(times(1)) > 100 * eps
    error('times must be sorted and start at zero.');
end

dx = xgrid(2) - xgrid(1);
period = n * dx;
k = fourier_wavenumbers(n, period);
L = 1i * k.^3;
m = parser.Results.ContourPoints;
[E, E2, Q, f1, f2, f3] = etdrk4_coefficients(L, dt, m);

if parser.Results.Dealias
    mask = abs(k) <= (2 / 3) * max(abs(k));
else
    mask = true(size(k));
end

v = fft(q0);
values = zeros(n, numel(times));
values(:, 1) = q0;
t = 0.0;

    function nv = nonlinear(vhat)
        u = real(ifft(vhat));
        nv = 3i * k .* fft(u.^2);
        nv(~mask) = 0;
    end

for outIdx = 2:numel(times)
    target = times(outIdx);
    while t < target - 100 * eps(max(1, target))
        h = min(dt, target - t);
        if abs(h - dt) <= 100 * eps(max(1, dt))
            Estep = E;
            E2step = E2;
            Qstep = Q;
            f1step = f1;
            f2step = f2;
            f3step = f3;
        else
            [Estep, E2step, Qstep, f1step, f2step, f3step] = etdrk4_coefficients(L, h, m);
        end
        Nv = nonlinear(v);
        a = E2step .* v + Qstep .* Nv;
        Na = nonlinear(a);
        b = E2step .* v + Qstep .* Na;
        Nb = nonlinear(b);
        c = E2step .* a + Qstep .* (2 * Nb - Nv);
        Nc = nonlinear(c);
        v = Estep .* v + f1step .* Nv + 2 * f2step .* (Na + Nb) + f3step .* Nc;
        v(~mask) = 0;
        t = t + h;
    end
    values(:, outIdx) = real(ifft(v));
end

sol = struct('xgrid', xgrid, 'times', times, 'values', values, 'dt', dt);
end

function [E, E2, Q, f1, f2, f3] = etdrk4_coefficients(L, h, m)
r = exp(1i * pi * (((1:m) - 0.5) / m));
LR = h * L + r;
E = exp(h * L);
E2 = exp(h * L / 2);
Q = h * mean((exp(LR / 2) - 1) ./ LR, 2);
f1 = h * mean((-4 - LR + exp(LR) .* (4 - 3 * LR + LR.^2)) ./ LR.^3, 2);
f2 = h * mean((2 + LR + exp(LR) .* (-2 + LR)) ./ LR.^3, 2);
f3 = h * mean((-4 - 3 * LR - LR.^2 + exp(LR) .* (4 - LR)) ./ LR.^3, 2);
end
