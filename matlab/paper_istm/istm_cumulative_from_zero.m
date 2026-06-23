function value = istm_cumulative_from_zero(x, f)
%ISTM_CUMULATIVE_FROM_ZERO Signed cumulative trapezoidal integral from 0 to x.
x = x(:);
f = f(:);
cum = cumtrapz(x, f);
zeroValue = interp1(x, cum, 0, 'linear', 'extrap');
value = cum - zeroValue;
end
