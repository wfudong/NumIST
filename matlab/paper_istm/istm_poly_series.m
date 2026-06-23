function [s, ds] = istm_poly_series(coeffs, z)
%ISTM_POLY_SERIES Sum s(z)=sum_n (-1)^n coeffs(n+1) z^n and derivative.
coeffs = coeffs(:);
zShape = size(z);
z = z(:).';
nmax = numel(coeffs) - 1;
s = zeros(size(z));
ds = zeros(size(z));
for n = nmax:-1:0
    s = s .* z + (-1)^n * coeffs(n + 1);
end
if nmax > 0
    for n = nmax:-1:1
        ds = ds .* z + (-1)^n * n * coeffs(n + 1);
    end
end
s = reshape(s, zShape);
ds = reshape(ds, zShape);
end
