%% MAIN_KDV_PLAYGROUND
% Edit this first block, then press Run.
%
% Equation:
%     u_t - 6*u*u_x + u_xxx = 0
%
% Notes:
% - The Fourier solver is periodic on [xmin, xmax). Choose a large box when
%   you want to mimic whole-line decay.
% - Smaller dt is safer for rough or large-amplitude data.

clear; close all; clc;

%% Things to play with
initialData = "custom";     % "example91", "soliton", "example93", "custom"
targetTimes = [0, 0.08, 0.16]; % Put any nonnegative times here, e.g. [0, 0.137]

xmin = -70;
xmax = 70;
numGrid = 4096;                % Must be even.
dt = 2e-4;                     % Maximum time step.

xPlotMin = -6;
xPlotMax = 8;
numPlotPoints = 1000;
yLimits = [];                  % [] for automatic limits, or e.g. [-1, 1].

plotMode = "stacked";          % "stacked" or "overlay"
saveFigure = true;
outputDir = fullfile(fileparts(mfilename('fullpath')), 'figures_playground');

% Used only when initialData = "custom".
% Try, for example:
customQ = @(x) -12 ./ cosh(x).^2;
% customQ = @(x) 0.8 * exp(-x.^2) .* cos(3*x);
% customQ = @(x) x .* exp(-x.^2);

%% Choose initial data
switch initialData
    case "example91"
        qfun = @paper_example91_initial;
        titleText = 'q(x) = x exp(-x^2)';
    case "soliton"
        c = pi;
        qfun = @(x) paper_example92_exact(x, 0, c);
        titleText = 'q(x) = -pi/2 sech^2(sqrt(pi) x / 2)';
        xmin = -40; xmax = 40; dt = min(dt, 2e-4);
    case "example93"
        qfun = @paper_example93_initial;
        titleText = 'piecewise exp/Bessel initial data';
        xmin = -48; xmax = 48; dt = min(dt, 5e-5);
        xPlotMin = -7; xPlotMax = 7;
    case "custom"
        qfun = customQ;
        titleText = 'custom initial data';
    otherwise
        error('Unknown initialData value: %s', initialData);
end

%% Solve
targetTimes = unique([0, targetTimes(:).'], 'stable');
if any(targetTimes < 0) || any(diff(targetTimes) < 0)
    error('targetTimes must be sorted and nonnegative.');
end
if mod(numGrid, 2) ~= 0
    error('numGrid must be even.');
end

xgrid = linspace(xmin, xmax, numGrid + 1).';
xgrid(end) = [];
q0 = qfun(xgrid);

fprintf('Solving %s on [%g, %g) with N=%d, dt<=%g...\n', ...
    titleText, xmin, xmax, numGrid, dt);
sol = solve_kdv_etdrk4(q0, xgrid, targetTimes, dt);

%% Plot
xplot = linspace(xPlotMin, xPlotMax, numPlotPoints).';
fig = figure('Color', 'w', 'Position', [100, 100, 760, 560]);

switch plotMode
    case "overlay"
        ax = axes(fig);
        hold(ax, 'on');
        for j = 1:numel(sol.times)
            y = interpolate_periodic(xgrid, sol.values(:, j), xplot);
            plot(ax, xplot, y, 'LineWidth', 1.5, 'DisplayName', time_label(sol.times(j)));
        end
        hold(ax, 'off');
        xlabel(ax, 'x');
        ylabel(ax, 'u(x,t)');
        title(ax, titleText, 'Interpreter', 'none');
        legend(ax, 'Location', 'best');
        grid(ax, 'on');
        box(ax, 'on');
        xlim(ax, [xPlotMin, xPlotMax]);
        if ~isempty(yLimits)
            ylim(ax, yLimits);
        end
    case "stacked"
        tiledlayout(fig, numel(sol.times), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
        for j = 1:numel(sol.times)
            ax = nexttile;
            y = interpolate_periodic(xgrid, sol.values(:, j), xplot);
            plot(ax, xplot, y, 'k-', 'LineWidth', 1.5);
            xlabel(ax, 'x');
            ylabel(ax, 'u');
            title(ax, sprintf('%s, %s', titleText, time_label(sol.times(j))), 'Interpreter', 'none');
            grid(ax, 'on');
            box(ax, 'on');
            xlim(ax, [xPlotMin, xPlotMax]);
            if ~isempty(yLimits)
                ylim(ax, yLimits);
            end
        end
    otherwise
        error('plotMode must be "stacked" or "overlay".');
end

%% Save optional output
if saveFigure
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    stamp = datestr(now, 'yyyymmdd_HHMMSS');
    outPng = fullfile(outputDir, ['kdv_playground_', stamp, '.png']);
    try
        exportgraphics(fig, outPng, 'Resolution', 220);
    catch
        print(fig, outPng, '-dpng', '-r220');
    end
    fprintf('Saved plot:\n  %s\n', outPng);
end

%% Quick invariant diagnostics
dx = xgrid(2) - xgrid(1);
mass = dx * sum(sol.values, 1);
l2 = dx * sum(sol.values.^2, 1);
diagnostics = table(sol.times(:), mass(:), abs(mass(:) - mass(1)), ...
    l2(:), abs(l2(:) - l2(1)), ...
    'VariableNames', {'time', 'mass', 'massDrift', 'l2', 'l2Drift'});
disp(diagnostics);

function label = time_label(t)
label = sprintf('t = %.12g', t);
end
