function sol = spps_propagate_second_order(x, r, lambda, y0, dy0, numTerms, subintervalLength)
%SPPS_PROPAGATE_SECOND_ORDER Piecewise SPPS propagation for y'' = lambda*r(x)*y.
%
% This implements the subdivision strategy recommended in sppsmethod.pdf.
% x may be increasing or decreasing. The initial data are imposed at x(1).

x = x(:);
r = r(:);
if numel(x) ~= numel(r)
    error('x and r must have the same length.');
end
if nargin < 7 || isempty(subintervalLength)
    subintervalLength = 2.0;
end

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
