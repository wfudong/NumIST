function u = paper_example92_exact(x, t, c)
%PAPER_EXAMPLE92_EXACT Exact negative KdV one-soliton from Example 9.2.
if nargin < 3
    c = pi;
end
arg = sqrt(c) .* (x - c * t) / 2;
u = -c / 2 ./ cosh(arg).^2;
end
