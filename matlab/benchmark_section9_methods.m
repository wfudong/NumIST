function summary = benchmark_section9_methods(varargin)
%BENCHMARK_SECTION9_METHODS Time paper ISTM versus ETDRK4 for Section 9.
%
% The timings use the same numerical settings as compare_section9_methods.m,
% but omit figure generation.  The ISTM time is split into direct scattering,
% reflection-grid sampling, and finite-section inverse recovery.
%
% Example:
%   summary = benchmark_section9_methods('NumTheta', 10000, 'Repeats', 1);

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, 'paper_istm'));

parser = inputParser;
parser.addParameter('OutputDir', fullfile(thisDir, 'figures_speed_comparison'), ...
    @(s) ischar(s) || isstring(s));
parser.addParameter('NumTheta', 10000, @(x) isnumeric(x) && isscalar(x) && x >= 100);
parser.addParameter('Repeats', 1, @(x) isnumeric(x) && isscalar(x) && x >= 1);
parser.addParameter('NumX', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x >= 101));
parser.addParameter('MakePlot', true, @(x) islogical(x) && isscalar(x));
parser.addParameter('SaveMat', true, @(x) islogical(x) && isscalar(x));
parser.addParameter('IncludeExample92DirectDiagnostic', true, @(x) islogical(x) && isscalar(x));
parser.parse(varargin{:});

outdir = char(parser.Results.OutputDir);
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

numTheta = parser.Results.NumTheta;
repeats = parser.Results.Repeats;
numX = parser.Results.NumX;

fprintf('Benchmarking Section 9 methods with NumTheta=%d, Repeats=%d.\n', numTheta, repeats);
if ~isempty(numX)
    fprintf('Using NumX=%d points for each displayed x interval.\n', numX);
end
fprintf('Plot/PDF writing is excluded from the measured times.\n\n');

summary = struct();
summary.generatedAt = datetime('now');
summary.numTheta = numTheta;
summary.repeats = repeats;
summary.numXOverride = numX;
summary.notes = ['Example 9.2 uses exact reflectionless scattering data for ', ...
    'the plotted ISTM inverse; its direct scattering run is timed as a diagnostic.'];

summary.example91 = benchmark_example91(numTheta, repeats, numX);
summary.example92 = benchmark_example92(repeats, parser.Results.IncludeExample92DirectDiagnostic, numX);
summary.example93 = benchmark_example93(numTheta, repeats, numX);

summary.summaryFile = write_speed_summary(outdir, summary);
if parser.Results.MakePlot
    summary.plotFile = plot_speed_summary(outdir, summary);
end
if parser.Results.SaveMat
    save(fullfile(outdir, 'section9_speed_summary.mat'), 'summary');
end

fprintf('\nWrote speed summary to:\n  %s\n', summary.summaryFile);
if isfield(summary, 'plotFile')
    fprintf('Wrote speed plot to:\n  %s\n', summary.plotFile);
end
end

function result = benchmark_example91(numTheta, repeats, numX)
cfg = struct();
cfg.name = 'Example 9.1';
cfg.qfun = @paper_example91_initial;
cfg.times = [0, 0.08, 0.16];
cfg.xplot = linspace(-5, 7, choose_num_x(numX, 1601)).';
cfg.etdrkBox = [-40, 40];
cfg.etdrkN = 4096;
cfg.etdrkDt = 2e-4;
cfg.istmNs = 5;

result = run_reflection_benchmark(cfg, repeats, numTheta, ...
    @() istm_paper_direct_scattering(cfg.qfun, ...
        'HalfLength', 10, 'NumGrid', 10001, 'NumCoefficients', 90, ...
        'SppsTerms', 70, 'RootGridSize', 12000));
end

function result = benchmark_example92(repeats, includeDirectDiagnostic, numX)
cfg = struct();
cfg.name = 'Example 9.2';
cfg.qfun = @(x) paper_example92_exact(x, 0, pi);
cfg.times = [0, 0.5, 1];
cfg.xplot = linspace(-5, 7, choose_num_x(numX, 3001)).';
cfg.etdrkBox = [-40, 40];
cfg.etdrkN = 4096;
cfg.etdrkDt = 2e-4;
cfg.istmNs = 5;

result = run_reflectionless_benchmark(cfg, repeats);

if includeDirectDiagnostic
    directSamples = zeros(repeats, 1);
    for r = 1:repeats
        fprintf('%s: direct-scattering diagnostic repeat %d/%d...\n', cfg.name, r, repeats);
        timer = tic;
        sc = istm_paper_direct_scattering(cfg.qfun, ...
            'HalfLength', 8, 'NumGrid', 5001, 'NumCoefficients', 60, ...
            'SppsTerms', 70, 'RootGridSize', 8000);
        directSamples(r) = toc(timer);
    end
    result.istmDirectDiagnosticSeconds = median(directSamples);
    result.istmDirectDiagnosticSamples = directSamples;
    result.diagnosticLambda = sc.lambda;
end
end

function result = benchmark_example93(numTheta, repeats, numX)
cfg = struct();
cfg.name = 'Example 9.3';
cfg.qfun = @paper_example93_initial;
cfg.times = [0, 0.015, 0.03];
cfg.xplot = linspace(-7, 7, choose_num_x(numX, 1601)).';
cfg.etdrkBox = [-48, 48];
cfg.etdrkN = 4096;
cfg.etdrkDt = 5e-5;
cfg.istmNs = 9;

result = run_reflection_benchmark(cfg, repeats, numTheta, ...
    @() istm_paper_direct_scattering(cfg.qfun, ...
        'HalfLength', 10, 'NumGrid', 5001, 'NumCoefficients', 90, ...
        'SppsTerms', 70, 'RootGridSize', 8000));
end

function n = choose_num_x(numX, defaultValue)
if isempty(numX)
    n = defaultValue;
else
    n = round(numX);
end
end

function result = run_reflection_benchmark(cfg, repeats, numTheta, directFun)
samples = empty_samples(repeats);
for r = 1:repeats
    fprintf('%s: ISTM direct scattering repeat %d/%d...\n', cfg.name, r, repeats);
    timer = tic;
    sc = directFun();
    samples.istmDirect(r) = toc(timer);

    fprintf('%s: ISTM reflection grid repeat %d/%d...\n', cfg.name, r, repeats);
    timer = tic;
    [~, ~, rho, sPlus, sMinus] = istm_reflection_grid_from_direct(sc, numTheta);
    rightData = istm_make_scattering_data('right', sc.tau, sc.alphaPlus, rho, sPlus);
    leftData = istm_make_scattering_data('left', sc.tau, sc.alphaMinus, rho, sMinus);
    samples.istmReflectionGrid(r) = toc(timer);

    fprintf('%s: ISTM inverse repeat %d/%d...\n', cfg.name, r, repeats);
    timer = tic;
    uIstm = recover_istm_times(rightData, leftData, cfg.xplot, cfg.times, cfg.istmNs);
    samples.istmInverse(r) = toc(timer);

    fprintf('%s: ETDRK4 repeat %d/%d...\n', cfg.name, r, repeats);
    timer = tic;
    uEtdrk = recover_etdrk_times(cfg.qfun, cfg.etdrkBox, cfg.etdrkN, ...
        cfg.times, cfg.etdrkDt, cfg.xplot);
    samples.etdrk(r) = toc(timer);
end

result = finish_timing_result(cfg, samples, uIstm, uEtdrk);
result.numTheta = numTheta;
result.numEigenvalues = numel(sc.lambda);
end

function result = run_reflectionless_benchmark(cfg, repeats)
samples = empty_samples(repeats);
for r = 1:repeats
    fprintf('%s: reflectionless ISTM inverse repeat %d/%d...\n', cfg.name, r, repeats);
    tau = sqrt(pi) / 2;
    alpha = sqrt(pi);
    rightData = istm_make_scattering_data('right', tau, alpha, [], []);
    leftData = istm_make_scattering_data('left', tau, alpha, [], []);

    timer = tic;
    uIstm = recover_istm_times(rightData, leftData, cfg.xplot, cfg.times, cfg.istmNs);
    samples.istmInverse(r) = toc(timer);

    fprintf('%s: ETDRK4 repeat %d/%d...\n', cfg.name, r, repeats);
    timer = tic;
    uEtdrk = recover_etdrk_times(cfg.qfun, cfg.etdrkBox, cfg.etdrkN, ...
        cfg.times, cfg.etdrkDt, cfg.xplot);
    samples.etdrk(r) = toc(timer);
end

result = finish_timing_result(cfg, samples, uIstm, uEtdrk);
result.numTheta = 0;
result.numEigenvalues = 1;
end

function samples = empty_samples(repeats)
samples = struct();
samples.istmDirect = zeros(repeats, 1);
samples.istmReflectionGrid = zeros(repeats, 1);
samples.istmInverse = zeros(repeats, 1);
samples.etdrk = zeros(repeats, 1);
end

function result = finish_timing_result(cfg, samples, uIstm, uEtdrk)
istmDirect = median(samples.istmDirect);
istmReflectionGrid = median(samples.istmReflectionGrid);
istmInverse = median(samples.istmInverse);
etdrk = median(samples.etdrk);
istmTotal = istmDirect + istmReflectionGrid + istmInverse;

diffValues = uIstm - uEtdrk;
result = struct();
result.name = cfg.name;
result.times = cfg.times;
result.numX = numel(cfg.xplot);
result.istmNs = cfg.istmNs;
result.etdrkN = cfg.etdrkN;
result.etdrkDt = cfg.etdrkDt;
result.samples = samples;
result.istmDirectSeconds = istmDirect;
result.istmReflectionGridSeconds = istmReflectionGrid;
result.istmInverseSeconds = istmInverse;
result.istmTotalSeconds = istmTotal;
result.etdrkSeconds = etdrk;
result.speedRatioIstmOverEtdrk = istmTotal / etdrk;
result.maxAbsDifference = max(abs(diffValues), [], 1);
result.rmsDifference = sqrt(mean(diffValues.^2, 1));
end

function values = recover_istm_times(rightData, leftData, xplot, times, Ns)
values = zeros(numel(xplot), numel(times));
for j = 1:numel(times)
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

function summaryFile = write_speed_summary(outdir, summary)
summaryFile = fullfile(outdir, 'section9_speed_summary.txt');
fid = fopen(summaryFile, 'w');
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, 'Section 9 speed comparison: paper ISTM finite-section method vs ETDRK4\n');
fprintf(fid, 'Generated: %s\n', char(summary.generatedAt));
fprintf(fid, 'NumTheta: %d\n', summary.numTheta);
fprintf(fid, 'Repeats: %d (reported values are medians)\n', summary.repeats);
if ~isempty(summary.numXOverride)
    fprintf(fid, 'NumX override: %d points on each displayed x interval\n', summary.numXOverride);
end
fprintf(fid, 'Note: plot/PDF writing is excluded.\n');
fprintf(fid, 'Note: %s\n\n', summary.notes);

fprintf(fid, '%-12s %10s %10s %10s %10s %10s %10s\n', ...
    'Example', 'ISTM-dir', 'ISTM-rho', 'ISTM-inv', 'ISTM-tot', 'ETDRK4', 'ratio');
fprintf(fid, '%-12s %10s %10s %10s %10s %10s %10s\n', ...
    '', 'seconds', 'seconds', 'seconds', 'seconds', 'seconds', 'I/E');

names = {'example91', 'example92', 'example93'};
for k = 1:numel(names)
    item = summary.(names{k});
    fprintf(fid, '%-12s %10.4f %10.4f %10.4f %10.4f %10.4f %10.2f\n', ...
        item.name, item.istmDirectSeconds, item.istmReflectionGridSeconds, ...
        item.istmInverseSeconds, item.istmTotalSeconds, item.etdrkSeconds, ...
        item.speedRatioIstmOverEtdrk);
end

if isfield(summary.example92, 'istmDirectDiagnosticSeconds')
    fprintf(fid, '\nExample 9.2 direct-scattering diagnostic: %.4f seconds\n', ...
        summary.example92.istmDirectDiagnosticSeconds);
end

fprintf(fid, '\nAccuracy diagnostic from the same run:\n');
for k = 1:numel(names)
    item = summary.(names{k});
    fprintf(fid, '%s\n', item.name);
    for j = 1:numel(item.times)
        fprintf(fid, '  t=%0.6g  max|ISTM-ETDRK4|=% .6e  rms=% .6e\n', ...
            item.times(j), item.maxAbsDifference(j), item.rmsDifference(j));
    end
end
end

function plotFile = plot_speed_summary(outdir, summary)
names = {'example91', 'example92', 'example93'};
labels = strings(1, numel(names));
timings = zeros(numel(names), 2);
for k = 1:numel(names)
    item = summary.(names{k});
    labels(k) = string(item.name);
    timings(k, :) = [item.istmTotalSeconds, item.etdrkSeconds];
end

fig = figure('Color', 'w', 'Visible', 'off', 'Position', [100, 100, 760, 420]);
ax = axes(fig);
bar(ax, timings);
set(ax, 'YScale', 'log', 'XTickLabel', labels, 'Box', 'on', 'FontSize', 10);
ylabel(ax, 'seconds, log scale');
legend(ax, {'ISTM total', 'ETDRK4'}, 'Location', 'northwest');
title(ax, 'Section 9 runtime comparison');
grid(ax, 'on');

plotFile = fullfile(outdir, 'section9_speed_bar.png');
try
    exportgraphics(fig, plotFile, 'Resolution', 220);
catch
    print(fig, plotFile, '-dpng', '-r220');
end
close(fig);
end
