function out = istm_eval_egEG(z, data)
%ISTM_EVAL_EGEG Evaluate equations (5.1)-(5.4) and z-derivatives.
[sa, dsa] = istm_poly_series(data.a0, z);
[sb, dsb] = istm_poly_series(data.b0, z);
[sd, dsd] = istm_poly_series(data.d0, z);
[sc, dsc] = istm_poly_series(data.c0, z);

out.e = 1 + (z + 1) .* sa;
out.g = 1 + (z + 1) .* sb;
out.E = (z - 1) ./ (2 * (z + 1)) - 0.5 * data.qPlus + (z + 1) .* sd;
out.G = -(z - 1) ./ (2 * (z + 1)) + 0.5 * data.qMinus + (z + 1) .* sc;

out.ez = sa + (z + 1) .* dsa;
out.gz = sb + (z + 1) .* dsb;
out.Ez = 1 ./ (z + 1).^2 + sd + (z + 1) .* dsd;
out.Gz = -1 ./ (z + 1).^2 + sc + (z + 1) .* dsc;
out.Phi = out.e .* out.G - out.E .* out.g;
out.Phiz = out.ez .* out.G + out.e .* out.Gz - out.Ez .* out.g - out.E .* out.gz;
end
