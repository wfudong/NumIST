function summary = run_section9_istm_paper_method(varargin)
%RUN_SECTION9_ISTM_PAPER_METHOD Reproduce Section 9 using the ISTM algorithm.
%
% This driver follows 152.pdf:
% 1. Direct scattering via the Appendix A recurrent coefficients, with the
%    Jost values at rho=i/2 computed by SPPS.
% 2. KdV scattering-data evolution.
% 3. Inverse scattering by finite sections (9.2)-(9.3), with the
%    unit-circle quadrature (9.4).
%
% The older run_section9_examples.m script is a pseudospectral PDE
% reproduction; this file is the paper-method implementation.

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, 'paper_istm'));

parser = inputParser;
parser.addParameter('OutputDir', fullfile(thisDir, 'figures_section9_istm'), ...
    @(s) ischar(s) || isstring(s));
parser.addParameter('NumTheta', 10000, @(x) isnumeric(x) && isscalar(x) && x >= 100);
parser.addParameter('SavePdf', true, @(x) islogical(x) && isscalar(x));
parser.parse(varargin{:});

outdir = char(parser.Results.OutputDir);
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

summary = struct();

fprintf('Example 9.1: direct scattering...\n');
sc91 = istm_paper_direct_scattering(@paper_example91_initial, ...
    'HalfLength', 10, 'NumGrid', 10001, 'NumCoefficients', 90, ...
    'SppsTerms', 70, 'RootGridSize', 12000);
[~, ~, rho91, sp91, sm91] = istm_reflection_grid_from_direct(sc91, parser.Results.NumTheta);
right91 = istm_make_scattering_data('right', sc91.tau, sc91.alphaPlus, rho91, sp91);
left91 = istm_make_scattering_data('left', sc91.tau, sc91.alphaMinus, rho91, sm91);
times91 = [0, 0.08, 0.16];
x91 = linspace(-5, 7, 2001).';
u91 = recover_times(right91, left91, x91, times91, 5);
summary.example91 = struct('scattering', sc91, 'times', times91, 'x', x91, 'u', u91, ...
    't0Error', max(abs(u91(:, 1) - paper_example91_initial(x91))));
plot_stack(x91, u91, times91, [-5, 7], [-1, 1], ...
    fullfile(outdir, 'figure1_example91_istm'), parser.Results.SavePdf);

fprintf('Example 9.2: direct-scattering diagnostic and reflectionless inverse...\n');
sc92 = istm_paper_direct_scattering(@(x) paper_example92_exact(x, 0, pi), ...
    'HalfLength', 8, 'NumGrid', 5001, 'NumCoefficients', 60, ...
    'SppsTerms', 70, 'RootGridSize', 8000);
tau92 = sqrt(pi) / 2;
alpha92 = sqrt(pi);
right92 = istm_make_scattering_data('right', tau92, alpha92, [], []);
left92 = istm_make_scattering_data('left', tau92, alpha92, [], []);
times92 = [0, 0.5, 1];
x92 = linspace(-5, 7, 12001).';
u92 = recover_times(right92, left92, x92, times92, 5);
summary.example92 = struct('scattering', sc92, 'times', times92, 'x', x92, 'u', u92, ...
    't0Error', max(abs(u92(:, 1) - paper_example92_exact(x92, 0, pi))), ...
    't1Error', max(abs(u92(:, 3) - paper_example92_exact(x92, 1, pi))));
plot_stack(x92, u92, times92, [-5, 7], [-2, 0], ...
    fullfile(outdir, 'figure2_example92_istm'), parser.Results.SavePdf);

fprintf('Example 9.3: direct scattering...\n');
sc93 = istm_paper_direct_scattering(@paper_example93_initial, ...
    'HalfLength', 10, 'NumGrid', 5001, 'NumCoefficients', 90, ...
    'SppsTerms', 70, 'RootGridSize', 8000);
[~, ~, rho93, sp93, sm93] = istm_reflection_grid_from_direct(sc93, parser.Results.NumTheta);
right93 = istm_make_scattering_data('right', sc93.tau, sc93.alphaPlus, rho93, sp93);
left93 = istm_make_scattering_data('left', sc93.tau, sc93.alphaMinus, rho93, sm93);

x93initial = linspace(-8, 8, 1200).';
q93initial = paper_example93_initial(x93initial);
plot_initial_profile(x93initial, q93initial, fullfile(outdir, 'figure3_example93_initial_istm'), parser.Results.SavePdf);

times93 = [0, 0.015, 0.03];
x93 = linspace(-7, 7, 2001).';
u93 = recover_times(right93, left93, x93, times93, 9);
summary.example93 = struct('scattering', sc93, 'times', times93, 'x', x93, 'u', u93, ...
    't0Error', max(abs(u93(:, 1) - paper_example93_initial(x93))));
plot_stack(x93, u93, times93, [-7, 7], [-1, 1], ...
    fullfile(outdir, 'figure4_example93_istm'), parser.Results.SavePdf);

summary.summaryFile = write_summary(outdir, summary);
fprintf('Wrote ISTM figures to:\n  %s\n', outdir);
fprintf('Validation summary:\n  %s\n', summary.summaryFile);
end

function values = recover_times(rightData, leftData, x, times, Ns)
values = zeros(numel(x), numel(times));
for j = 1:numel(times)
    fprintf('  inverse scattering at t=%g, Ns=%d...\n', times(j), Ns);
    rt = istm_evolve_scattering(rightData, times(j));
    lt = istm_evolve_scattering(leftData, times(j));
    rec = istm_recover_potential(rt, lt, x, 'Ns', Ns);
    values(:, j) = rec.q;
end
end

function plot_stack(x, values, times, xlims, ylims, basename, savePdf)
fig = figure('Color', 'w', 'Position', [100, 100, 560, 650]);
tiledlayout(fig, numel(times), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for j = 1:numel(times)
    ax = nexttile;
    plot(ax, x, values(:, j), 'k-', 'LineWidth', 1.6);
    xlim(ax, xlims);
    ylim(ax, ylims);
    xlabel(ax, 'x');
    set(ax, 'Box', 'on', 'XTick', ceil(xlims(1) / 2) * 2:2:floor(xlims(2) / 2) * 2, ...
        'YTick', ceil(ylims(1)):floor(ylims(2)), 'FontSize', 9);
    grid(ax, 'off');
    text(ax, xlims(1) + 0.13 * diff(xlims), ylims(2) - 0.15 * diff(ylims), ...
        sprintf('t=%g', times(j)), 'FontSize', 10, 'Color', 'k');
    if j == 1
        text(ax, xlims(1) + 0.62 * diff(xlims), ylims(2) - 0.14 * diff(ylims), ...
            'Initial Data', 'FontSize', 10, 'Color', 'k');
    end
end
save_figure(fig, basename, savePdf);
close(fig);
end

function plot_initial_profile(x, y, basename, savePdf)
fig = figure('Color', 'w', 'Position', [120, 120, 560, 280]);
ax = axes(fig);
plot(ax, x, y, 'k-', 'LineWidth', 1.6);
xlim(ax, [-8, 8]);
ylim(ax, [-0.5, 1.05]);
xlabel(ax, 'x');
ylabel(ax, 'q(x)');
set(ax, 'Box', 'on', 'XTick', -8:2:8, 'YTick', [-0.5, 0, 0.5, 1], 'FontSize', 9);
grid(ax, 'off');
save_figure(fig, basename, savePdf);
close(fig);
end

function save_figure(fig, basename, savePdf)
try
    exportgraphics(fig, [basename, '.png'], 'Resolution', 220);
catch
    print(fig, [basename, '.png'], '-dpng', '-r220');
end
if savePdf
    try
        exportgraphics(fig, [basename, '.pdf'], 'ContentType', 'vector');
    catch
        print(fig, [basename, '.pdf'], '-dpdf', '-painters');
    end
end
end

function path = write_summary(outdir, summary)
path = fullfile(outdir, 'validation_summary_istm.txt');
fid = fopen(path, 'w');
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, 'Validation summary for Section 9 ISTM paper-method MATLAB implementation\n\n');

sc91 = summary.example91.scattering;
fprintf(fid, 'Example 9.1\n');
fprintf(fid, '  Paper: lambda1 ~= -0.0138384593995, alpha- ~= 0.2055954681199, alpha+ ~= 0.0416040800785\n');
fprintf(fid, '  Computed lambda1 = %.16f\n', sc91.lambda(1));
fprintf(fid, '  Computed alpha-  = %.16f\n', sc91.alphaMinus(1));
fprintf(fid, '  Computed alpha+  = %.16f\n', sc91.alphaPlus(1));
fprintf(fid, '  t=0 inverse reconstruction max error on [-5,7]: %.6e\n\n', summary.example91.t0Error);

sc92 = summary.example92.scattering;
fprintf(fid, 'Example 9.2\n');
fprintf(fid, '  Paper: exact lambda1 = -pi/4 = %.16f, alpha+- = sqrt(pi) = %.16f\n', -pi / 4, sqrt(pi));
if ~isempty(sc92.lambda)
    fprintf(fid, '  Direct-scattering diagnostic lambda1 = %.16f\n', sc92.lambda(1));
    fprintf(fid, '  Direct-scattering diagnostic alpha-  = %.16f\n', sc92.alphaMinus(1));
    fprintf(fid, '  Direct-scattering diagnostic alpha+  = %.16f\n', sc92.alphaPlus(1));
end
fprintf(fid, '  t=0 inverse reconstruction max error on [-5,7]: %.6e\n', summary.example92.t0Error);
fprintf(fid, '  t=1 inverse reconstruction max error on [-5,7]: %.6e\n\n', summary.example92.t1Error);

sc93 = summary.example93.scattering;
fprintf(fid, 'Example 9.3\n');
fprintf(fid, '  Computed number of discrete eigenvalues: %d\n', numel(sc93.lambda));
fprintf(fid, '  Paper reports t=0 inverse reconstruction max error about 6e-3 on [-7,7].\n');
fprintf(fid, '  Computed t=0 inverse reconstruction max error on [-7,7]: %.6e\n', summary.example93.t0Error);
end
