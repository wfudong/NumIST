function q = paper_example93_initial(x)
%PAPER_EXAMPLE93_INITIAL Piecewise initial profile from Example 9.3.
q = zeros(size(x));
left = x < 0;
q(left) = exp(x(left)) .* cos(4 * x(left));
q(~left) = exp(-x(~left)) .* besselj(0, 2 * x(~left));
end
