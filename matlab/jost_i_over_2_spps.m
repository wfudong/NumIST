function result = jost_i_over_2_spps(qfun, side, varargin)
%JOST_I_OVER_2_SPPS Approximate the Jost solution at xi = i/2.
%
% Section 9.1 of 152.pdf solves
%
%     y'' - q(x)*y - (1/4)*y = 0
%
% on a large finite interval, using the asymptotic Jost boundary data. This
% routine solves that IVP with the SPPS construction from sppsmethod.pdf,
% rewritten as y'' = r(x)y with r(x) = q(x)+1/4.
%
% side = 'right' returns e(i/2,x) on [0,L] using data at x=L.
% side = 'left'  returns g(i/2,x) on [-L,0] using data at x=-L.

parser = inputParser;
parser.addParameter('HalfLength', 20, @(x) isnumeric(x) && isscalar(x) && x > 0);
parser.addParameter('NumGrid', 801, @(x) isnumeric(x) && isscalar(x) && x >= 5);
parser.addParameter('NumTerms', 42, @(x) isnumeric(x) && isscalar(x) && x >= 1);
parser.addParameter('SubintervalLength', 2.0, @(x) isnumeric(x) && isscalar(x) && x > 0);
parser.parse(varargin{:});

L = parser.Results.HalfLength;
m = parser.Results.NumGrid;
numTerms = parser.Results.NumTerms;
subintervalLength = parser.Results.SubintervalLength;
side = validatestring(side, {'right', 'left'});

switch side
    case 'right'
        x = linspace(L, 0, m).';
        y0 = exp(-L / 2);
        dy0 = -0.5 * exp(-L / 2);
    case 'left'
        x = linspace(-L, 0, m).';
        y0 = exp(-L / 2);
        dy0 = 0.5 * exp(-L / 2);
end

q = qfun(x);
r = q + 0.25;
sol = propagate_spps_segments(x, r, 1.0, y0, dy0, numTerms, subintervalLength);
residual = equation_residual(x, sol.y, q);

result = sol;
result.side = side;
result.q = q;
result.maxResidual = max(abs(residual(3:end-2)));
result.yAtZero = sol.y(end);
result.dyAtZero = sol.dy(end);
result.residual = residual;
end

function sol = propagate_spps_segments(x, r, lambda, y0, dy0, numTerms, subintervalLength)
n = numel(x);
span = abs(x(end) - x(1));
numSegments = max(1, min(n - 1, ceil(span / subintervalLength)));
edges = unique(round(linspace(1, n, numSegments + 1)), 'stable');
if edges(end) ~= n
    edges(end + 1) = n;
end

y = nan(n, 1);
dy = nan(n, 1);
currentY = y0;
currentDy = dy0;

for k = 1:(numel(edges) - 1)
    first = edges(k);
    last = edges(k + 1);
    if last <= first
        continue;
    end
    idx = first:last;
    segment = spps_ivp_second_order(x(idx), r(idx), lambda, currentY, currentDy, numTerms);
    if k == 1
        y(idx) = segment.y;
        dy(idx) = segment.dy;
    else
        y(idx(2:end)) = segment.y(2:end);
        dy(idx(2:end)) = segment.dy(2:end);
    end
    currentY = segment.y(end);
    currentDy = segment.dy(end);
end

sol = struct('x', x, 'y', y, 'dy', dy);
end

function residual = equation_residual(x, y, q)
n = numel(x);
h = x(2) - x(1);
residual = nan(n, 1);
for j = 3:n-2
    yxx = (-y(j + 2) + 16 * y(j + 1) - 30 * y(j) + 16 * y(j - 1) - y(j - 2)) / (12 * h^2);
    residual(j) = yxx - (q(j) + 0.25) * y(j);
end
end
