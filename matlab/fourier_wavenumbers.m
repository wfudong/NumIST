function k = fourier_wavenumbers(n, period)
%FOURIER_WAVENUMBERS Return FFT-ordered angular wavenumbers.
if mod(n, 2) ~= 0
    error('n must be even.');
end
k = (2 * pi / period) * [0:(n / 2), (-(n / 2) + 1):-1].';
end
