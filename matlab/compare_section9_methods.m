function summary = compare_section9_methods(varargin)
%COMPARE_SECTION9_METHODS Overlay paper ISTM and ETDRK results for Section 9.
%
% This script compares two independent numerical routes:
%
%   1. Paper ISTM method from 152.pdf:
%      direct scattering, scattering-data evolution, finite-section inverse.
%   2. Fourier pseudospectral ETDRK4 time stepping for the KdV PDE.
%
% It writes overlay plots and absolute-difference plots for Examples 9.1,
% 9.2, and 9.3.

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, 'paper_istm'));

parser = inputParser;
parser.addParameter('OutputDir', fullfile(thisDir, 'figures_method_comparison'), ...
    @(s) ischar(s) || isstring(s));
parser.addParameter('NumTheta', 10000, @(x) isnumeric(x) && isscalar(x) && x >= 100);
parser.addParameter('SavePdf', true, @(x) islogical(x) && isscalar(x));
parser.parse(varargin{:});

outdir = char(parser.Results.OutputDir);
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

summary = struct();

summary.example91 = compare_example91(parser.Results.NumTheta, outdir, parser.Results.SavePdf);
summary.example92 = compare_example92(outdir, parser.Results.SavePdf);
summary.example93 = compare_example93(parser.Results.NumTheta, outdir, parser.Results.SavePdf);
summary.summaryFile = write_comparison_summary(outdir, summary);

fprintf('Wrote method-comparison figures to:\n  %s\n', outdir);
fprintf('Comparison summary:\n  %s\n', summary.summaryFile);
end

function result = compare_example91(numTheta, outdir, savePdf)
cfg = struct();
cfg.name = 'Example 9.1';
cfg.fileBase = 'example91';
cfg.qfun = @paper_example91_initial;
cfg.times = [0, 0.08, 0.16];
cfg.xplot = linspace(-5, 7, 1601).';
cfg.xlims = [-5, 7];
cfg.ylims = [-1, 1];
cfg.etdrkBox = [-40, 40];
cfg.etdrkN = 4096;
cfg.etdrkDt = 2e-4;

fprintf('%s: computing ISTM data...\n', cfg.name);
sc = istm_paper_direct_scattering(cfg.qfun, ...
    'HalfLength', 10, 'NumGrid', 10001, 'NumCoefficients', 90, ...
    'SppsTerms', 70, 'RootGridSize', 12000);
[~, ~, rho, sPlus, sMinus] = istm_reflection_grid_from_direct(sc, numTheta);
rightData = istm_make_scattering_data('right', sc.tau, sc.alphaPlus, rho, sPlus);
leftData = istm_make_scattering_data('left', sc.tau, sc.alphaMinus, rho, sMinus);
uIstm = recover_istm_times(rightData, leftData, cfg.xplot, cfg.times, 5);

fprintf('%s: computing ETDRK data...\n', cfg.name);
uEtdrk = recover_etdrk_times(cfg.qfun, cfg.etdrkBox, cfg.etdrkN, cfg.times, cfg.etdrkDt, cfg.xplot);

result = finish_comparison(cfg, uIstm, uEtdrk, outdir, savePdf);
result.scattering = sc;
end

function result = compare_example92(outdir, savePdf)
cfg = struct();
cfg.name = 'Example 9.2';
cfg.fileBase = 'example92';
cfg.qfun = @(x) paper_example92_exact(x, 0, pi);
cfg.times = [0, 0.5, 1];
cfg.xplot = linspace(-5, 7, 3001).';
cfg.xlims = [-5, 7];
cfg.ylims = [-2, 0.2];
cfg.etdrkBox = [-40, 40];
cfg.etdrkN = 4096;
cfg.etdrkDt = 2e-4;

fprintf('%s: computing reflectionless ISTM data...\n', cfg.name);
tau = sqrt(pi) / 2;
alpha = sqrt(pi);
rightData = istm_make_scattering_data('right', tau, alpha, [], []);
leftData = istm_make_scattering_data('left', tau, alpha, [], []);
uIstm = recover_istm_times(rightData, leftData, cfg.xplot, cfg.times, 5);

fprintf('%s: computing ETDRK data...\n', cfg.name);
uEtdrk = recover_etdrk_times(cfg.qfun, cfg.etdrkBox, cfg.etdrkN, cfg.times, cfg.etdrkDt, cfg.xplot);

result = finish_comparison(cfg, uIstm, uEtdrk, outdir, savePdf);
end

function result = compare_example93(numTheta, outdir, savePdf)
cfg = struct();
cfg.name = 'Example 9.3';
cfg.fileBase = 'example93';
cfg.qfun = @paper_example93_initial;
cfg.times = [0, 0.015, 0.03];
cfg.xplot = linspace(-7, 7, 1601).';
cfg.xlims = [-7, 7];
cfg.ylims = [-1, 1];
cfg.etdrkBox = [-48, 48];
cfg.etdrkN = 4096;
cfg.etdrkDt = 5e-5;

fprintf('%s: computing ISTM data...\n', cfg.name);
sc = istm_paper_direct_scattering(cfg.qfun, ...
    'HalfLength', 10, 'NumGrid', 5001, 'NumCoefficients', 90, ...
    'SppsTerms', 70, 'RootGridSize', 8000);
[~, ~, rho, sPlus, sMinus] = istm_reflection_grid_from_direct(sc, numTheta);
rightData = istm_make_scattering_data('right', sc.tau, sc.alphaPlus, rho, sPlus);
leftData = istm_make_scattering_data('left', sc.tau, sc.alphaMinus, rho, sMinus);
uIstm = recover_istm_times(rightData, leftData, cfg.xplot, cfg.times, 9);

fprintf('%s: computing ETDRK data...\n', cfg.name);
uEtdrk = recover_etdrk_times(cfg.qfun, cfg.etdrkBox, cfg.etdrkN, cfg.times, cfg.etdrkDt, cfg.xplot);

result = finish_comparison(cfg, uIstm, uEtdrk, outdir, savePdf);
result.scattering = sc;
end

function values = recover_istm_times(rightData, leftData, xplot, times, Ns)
values = zeros(numel(xplot), numel(times));
for j = 1:numel(times)
    fprintf('  ISTM inverse at t=%g, Ns=%d...\n', times(j), Ns);
    rt = istm_evolve_scattering(rightData, times(j));
    lt = istm_evolve_scattering(leftData, times(j));
    rec = istm_recover_potential(rt, lt, xplot, 'Ns', Ns);
    values(:, j) = rec.q;
end
end

function values = recover_etdrk_times(qfun, box, n, times, dt, xplot)
xgrid = linspace(box(1), box(2), n + 1).';
xgrid(end) = [];
q0 = qfun(xgrid);
sol = solve_kdv_etdrk4(q0, xgrid, times, dt);
values = zeros(numel(xplot), numel(times));
for j = 1:numel(times)
    values(:, j) = interpolate_periodic(xgrid, sol.values(:, j), xplot);
end
end

function result = finish_comparison(cfg, uIstm, uEtdrk, outdir, savePdf)
diffValues = uIstm - uEtdrk;
maxAbs = max(abs(diffValues), [], 1);
rmsDiff = sqrt(mean(diffValues.^2, 1));

plot_overlay(cfg, uIstm, uEtdrk, fullfile(outdir, cfg.fileBase + "_overlay"), savePdf);
plot_difference(cfg, abs(diffValues), fullfile(outdir, cfg.fileBase + "_difference"), savePdf);

result = struct();
result.name = cfg.name;
result.times = cfg.times;
result.x = cfg.xplot;
result.uIstm = uIstm;
result.uEtdrk = uEtdrk;
result.maxAbsDifference = maxAbs;
result.rmsDifference = rmsDiff;
end

function plot_overlay(cfg, uIstm, uEtdrk, basename, savePdf)
fig = figure('Color', 'w', 'Position', [100, 100, 720, 720]);
tiledlayout(fig, numel(cfg.times), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for j = 1:numel(cfg.times)
    ax = nexttile;
    plot(ax, cfg.xplot, uIstm(:, j), 'k-', 'LineWidth', 1.7, 'DisplayName', 'ISTM finite section');
    hold(ax, 'on');
    plot(ax, cfg.xplot, uEtdrk(:, j), 'r--', 'LineWidth', 1.15, 'DisplayName', 'ETDRK4');
    hold(ax, 'off');
    xlim(ax, cfg.xlims);
    ylim(ax, cfg.ylims);
    xlabel(ax, 'x');
    ylabel(ax, 'u');
    title(ax, sprintf('%s, t=%g', cfg.name, cfg.times(j)), 'Interpreter', 'none');
    set(ax, 'Box', 'on', 'FontSize', 9);
    grid(ax, 'on');
    if j == 1
        legend(ax, 'Location', 'best');
    end
end
save_comparison_figure(fig, basename, savePdf);
close(fig);
end

function plot_difference(cfg, absDiff, basename, savePdf)
fig = figure('Color', 'w', 'Position', [120, 120, 720, 720]);
tiledlayout(fig, numel(cfg.times), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for j = 1:numel(cfg.times)
    ax = nexttile;
    plot(ax, cfg.xplot, absDiff(:, j), 'b-', 'LineWidth', 1.35);
    xlim(ax, cfg.xlims);
    xlabel(ax, 'x');
    ylabel(ax, '|difference|');
    title(ax, sprintf('%s, |ISTM - ETDRK4|, t=%g', cfg.name, cfg.times(j)), ...
        'Interpreter', 'none');
    set(ax, 'Box', 'on', 'FontSize', 9);
    grid(ax, 'on');
end
save_comparison_figure(fig, basename, savePdf);
close(fig);
end

function save_comparison_figure(fig, basename, savePdf)
try
    exportgraphics(fig, [char(basename), '.png'], 'Resolution', 220);
catch
    print(fig, [char(basename), '.png'], '-dpng', '-r220');
end
if savePdf
    try
        exportgraphics(fig, [char(basename), '.pdf'], 'ContentType', 'vector');
    catch
        print(fig, [char(basename), '.pdf'], '-dpdf', '-painters');
    end
end
end

function path = write_comparison_summary(outdir, summary)
path = fullfile(outdir, 'method_comparison_summary.txt');
fid = fopen(path, 'w');
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, 'Section 9 comparison: paper ISTM finite-section method vs ETDRK4\n\n');
names = {'example91', 'example92', 'example93'};
for k = 1:numel(names)
    item = summary.(names{k});
    fprintf(fid, '%s\n', item.name);
    for j = 1:numel(item.times)
        fprintf(fid, '  t=%0.6g  max|ISTM-ETDRK4|=% .6e  rms=% .6e\n', ...
            item.times(j), item.maxAbsDifference(j), item.rmsDifference(j));
    end
    fprintf(fid, '\n');
end
end
