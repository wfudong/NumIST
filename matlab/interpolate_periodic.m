function yq = interpolate_periodic(xgrid, values, xq)
%INTERPOLATE_PERIODIC Linear interpolation on a uniform periodic grid.
xgrid = xgrid(:);
values = values(:);
shape = size(xq);
xq = xq(:);
n = numel(xgrid);
if numel(values) ~= n
    error('xgrid and values must have the same length.');
end
dx = xgrid(2) - xgrid(1);
period = n * dx;
xi = mod(xq - xgrid(1), period) / dx;
j0 = floor(xi);
theta = xi - j0;
j1 = mod(j0, n) + 1;
j2 = mod(j0 + 1, n) + 1;
yq = (1 - theta) .* values(j1) + theta .* values(j2);
yq = reshape(yq, shape);
end
