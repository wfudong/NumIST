function summary = run_section9_examples(varargin)
%RUN_SECTION9_EXAMPLES Reproduce the Section 9 examples from 152.pdf.
%
% The examples use the KdV convention in the paper,
%
%     u_t - 6*u*u_x + u_xxx = 0.
%
% Examples 9.1 and 9.3 are advanced with a Fourier pseudospectral ETDRK4
% method on a large periodic box. Example 9.2 is generated from the exact
% reflectionless solitary-wave formula reported in the paper. The auxiliary
% SPPS solver used by the paper for the Jost equation (9.1) is provided in
% jost_i_over_2_spps.m and is exercised by the optional SPPS check below.

parser = inputParser;
parser.addParameter('OutputDir', fullfile(fileparts(mfilename('fullpath')), 'figures_section9'), ...
    @(s) ischar(s) || isstring(s));
parser.addParameter('SavePdf', true, @(x) islogical(x) && isscalar(x));
parser.addParameter('RunSppsCheck', true, @(x) islogical(x) && isscalar(x));
parser.parse(varargin{:});

outdir = char(parser.Results.OutputDir);
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

summary = struct();
summary.example91 = run_example91(outdir, parser.Results.SavePdf);
summary.example92 = run_example92(outdir, parser.Results.SavePdf);
summary.example93 = run_example93(outdir, parser.Results.SavePdf);

if parser.Results.RunSppsCheck
    summary.spps = run_spps_checks();
else
    summary.spps = struct();
end

summary.summaryFile = write_validation_summary(outdir, summary);
fprintf('Wrote Section 9 MATLAB figures to:\n  %s\n', outdir);
fprintf('Validation summary:\n  %s\n', summary.summaryFile);
end

function result = run_example91(outdir, savePdf)
xmin = -40.0;
xmax = 40.0;
n = 4096;
xgrid = linspace(xmin, xmax, n + 1).';
xgrid(end) = [];
q0 = paper_example91_initial(xgrid);
times = [0.0, 0.08, 0.16];
sol = solve_kdv_etdrk4(q0, xgrid, times, 2e-4);

xplot = linspace(-5.0, 7.0, 900).';
labels = {'t=0', 't=0.08', 't=0.16'};
fig = figure('Color', 'w', 'Position', [100, 100, 560, 650]);
tiledlayout(fig, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for j = 1:numel(times)
    ax = nexttile;
    y = interpolate_periodic(xgrid, sol.values(:, j), xplot);
    paper_subplot(ax, xplot, y, [-5, 7], [-1, 1], labels{j}, j == 1);
end
save_figure(fig, fullfile(outdir, 'figure1_example91'), savePdf);
close(fig);

result = struct();
result.xgrid = xgrid;
result.times = times;
result.values = sol.values;
result.invariants = invariants(xgrid, sol.values, times);
result.description = 'q(x) = x exp(-x^2)';
end

function result = run_example92(outdir, savePdf)
c = pi;
times = [0.0, 0.5, 1.0];
xplot = linspace(-5.0, 7.0, 900).';
labels = {'t=0', 't=0.5', 't=1'};

fig = figure('Color', 'w', 'Position', [100, 100, 560, 650]);
tiledlayout(fig, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
values = zeros(numel(xplot), numel(times));
for j = 1:numel(times)
    ax = nexttile;
    values(:, j) = paper_example92_exact(xplot, times(j), c);
    paper_subplot(ax, xplot, values(:, j), [-5, 7], [-2, 0], labels{j}, j == 1);
end
save_figure(fig, fullfile(outdir, 'figure2_example92'), savePdf);
close(fig);

result = struct();
result.xplot = xplot;
result.times = times;
result.values = values;
result.c = c;
result.lambda1 = -c / 4;
result.alpha = sqrt(c);
result.description = '-c/2 sech^2(sqrt(c) x / 2), c = pi';
end

function result = run_example93(outdir, savePdf)
xinitial = linspace(-8.0, 8.0, 1200).';
yinitial = paper_example93_initial(xinitial);

fig3 = figure('Color', 'w', 'Position', [120, 120, 560, 280]);
ax = axes(fig3);
plot(ax, xinitial, yinitial, 'k-', 'LineWidth', 1.6);
xlim(ax, [-8, 8]);
ylim(ax, [-0.5, 1.05]);
xlabel(ax, 'x');
ylabel(ax, 'q(x)');
set(ax, 'Box', 'on', 'XTick', -8:2:8, 'YTick', [-0.5, 0, 0.5, 1.0], ...
    'FontSize', 9);
grid(ax, 'off');
save_figure(fig3, fullfile(outdir, 'figure3_example93_initial'), savePdf);
close(fig3);

xmin = -48.0;
xmax = 48.0;
n = 4096;
xgrid = linspace(xmin, xmax, n + 1).';
xgrid(end) = [];
q0 = paper_example93_initial(xgrid);
times = [0.0, 0.015, 0.03];
sol = solve_kdv_etdrk4(q0, xgrid, times, 5e-5);

xplot = linspace(-7.0, 7.0, 900).';
labels = {'t=0', 't=0.015', 't=0.03'};
fig4 = figure('Color', 'w', 'Position', [100, 100, 560, 650]);
tiledlayout(fig4, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for j = 1:numel(times)
    ax = nexttile;
    y = interpolate_periodic(xgrid, sol.values(:, j), xplot);
    paper_subplot(ax, xplot, y, [-7, 7], [-1, 1], labels{j}, j == 1);
end
save_figure(fig4, fullfile(outdir, 'figure4_example93'), savePdf);
close(fig4);

result = struct();
result.xinitial = xinitial;
result.yinitial = yinitial;
result.xgrid = xgrid;
result.times = times;
result.values = sol.values;
result.invariants = invariants(xgrid, sol.values, times);
result.description = 'exp(x) cos(4x), x<0; exp(-x) J0(2x), x>=0';
end

function checks = run_spps_checks()
checks = struct();
opts = {'HalfLength', 20, 'NumGrid', 801, 'NumTerms', 42};
checks.example91Right = jost_i_over_2_spps(@paper_example91_initial, 'right', opts{:});
checks.example91Left = jost_i_over_2_spps(@paper_example91_initial, 'left', opts{:});
checks.example93Right = jost_i_over_2_spps(@paper_example93_initial, 'right', opts{:});
checks.example93Left = jost_i_over_2_spps(@paper_example93_initial, 'left', opts{:});
end

function invtab = invariants(xgrid, values, times)
dx = xgrid(2) - xgrid(1);
mass = dx * sum(values, 1);
l2 = dx * sum(values.^2, 1);
invtab = table(times(:), mass(:), abs(mass(:) - mass(1)), l2(:), abs(l2(:) - l2(1)), ...
    'VariableNames', {'time', 'mass', 'massDrift', 'l2', 'l2Drift'});
end

function paper_subplot(ax, x, y, xlims, ylims, labeltext, initialText)
plot(ax, x, y, 'k-', 'LineWidth', 1.6);
xlim(ax, xlims);
ylim(ax, ylims);
xlabel(ax, 'x');
set(ax, 'Box', 'on', 'XTick', ceil(xlims(1) / 2) * 2:2:floor(xlims(2) / 2) * 2, ...
    'YTick', ceil(ylims(1)):floor(ylims(2)), 'FontSize', 9);
grid(ax, 'off');
text(ax, xlims(1) + 0.13 * diff(xlims), ylims(2) - 0.15 * diff(ylims), ...
    labeltext, 'FontSize', 10, 'Color', 'k');
if initialText
    if abs(ylims(1) + 2) < 100 * eps && abs(ylims(2)) < 100 * eps
        initialX = xlims(1) + 0.58 * diff(xlims);
        initialY = ylims(1) + 0.18 * diff(ylims);
    elseif abs(xlims(1) + 7) < 100 * eps && abs(xlims(2) - 7) < 100 * eps
        initialX = xlims(1) + 0.68 * diff(xlims);
        initialY = ylims(2) - 0.14 * diff(ylims);
    else
        initialX = xlims(1) + 0.48 * diff(xlims);
        initialY = ylims(2) - 0.10 * diff(ylims);
    end
    text(ax, initialX, initialY, 'Initial Data', 'FontSize', 10, 'Color', 'k');
end
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
    catch pdfError
        try
            print(fig, [basename, '.pdf'], '-dpdf', '-painters');
        catch
            warning('run_section9_examples:PdfExportSkipped', ...
                'Skipped PDF export for %s: %s', basename, pdfError.message);
        end
    end
end
end

function path = write_validation_summary(outdir, summary)
path = fullfile(outdir, 'validation_summary_matlab.txt');
fid = fopen(path, 'w');
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, 'Validation summary for MATLAB realization of 152.pdf Section 9 examples\n\n');
fprintf(fid, 'Example 9.1 / Figure 1\n');
fprintf(fid, '  Initial data: %s\n', summary.example91.description);
fprintf(fid, '  Paper time slices: t = 0, 0.08, 0.16\n');
fprintf(fid, '  Paper reports lambda1 ~= -0.0138384593995, alpha- ~= 0.2055954681199, alpha+ ~= 0.0416040800785.\n');
write_invariant_rows(fid, summary.example91.invariants);

fprintf(fid, '\nExample 9.2 / Figure 2\n');
fprintf(fid, '  Initial data: %s\n', summary.example92.description);
fprintf(fid, '  Exact lambda1 = -pi/4 = %.14f\n', summary.example92.lambda1);
fprintf(fid, '  Exact alpha+ = alpha- = sqrt(pi) = %.14f\n', summary.example92.alpha);
fprintf(fid, '  Paper time slices: t = 0, 0.5, 1\n');

fprintf(fid, '\nExample 9.3 / Figures 3 and 4\n');
fprintf(fid, '  Initial data: %s\n', summary.example93.description);
fprintf(fid, '  Paper time slices for Figure 4: t = 0, 0.015, 0.03\n');
write_invariant_rows(fid, summary.example93.invariants);

if isfield(summary, 'spps') && ~isempty(fieldnames(summary.spps))
    fprintf(fid, '\nSPPS checks for the auxiliary Jost equation y'''' - (q(x)+1/4)y = 0\n');
    names = fieldnames(summary.spps);
    for k = 1:numel(names)
        item = summary.spps.(names{k});
        fprintf(fid, '  %-16s side=%-5s max residual ~= %.3e, value at 0 ~= %.12e, derivative at 0 ~= %.12e\n', ...
            names{k}, item.side, item.maxResidual, item.yAtZero, item.dyAtZero);
    end
end
end

function write_invariant_rows(fid, invtab)
for j = 1:height(invtab)
    fprintf(fid, '  t=%0.5f mass=% .12e mass_drift=% .3e L2=% .12e L2_drift=% .3e\n', ...
        invtab.time(j), invtab.mass(j), invtab.massDrift(j), invtab.l2(j), invtab.l2Drift(j));
end
end
