function summary = animate_section9_examples(varargin)
%ANIMATE_SECTION9_EXAMPLES Animate the three Section 9 examples from 152.pdf.
%
% Examples:
%
%   animate_section9_examples
%   animate_section9_examples('Example', 'example93')
%   animate_section9_examples('Example', 'all', 'NumFrames', 160)
%   animate_section9_examples('SaveGif', false)
%
% The KdV convention is the one used in the paper:
%
%     u_t - 6*u*u_x + u_xxx = 0.

parser = inputParser;
parser.addParameter('Example', 'all', @(s) ischar(s) || isstring(s));
parser.addParameter('NumFrames', 121, @(x) isnumeric(x) && isscalar(x) && x >= 2);
parser.addParameter('SaveGif', true, @(x) islogical(x) && isscalar(x));
parser.addParameter('OutputDir', fullfile(fileparts(mfilename('fullpath')), 'animations_section9'), ...
    @(s) ischar(s) || isstring(s));
parser.addParameter('FrameDelay', 0.05, @(x) isnumeric(x) && isscalar(x) && x >= 0);
parser.addParameter('PauseTime', 0.01, @(x) isnumeric(x) && isscalar(x) && x >= 0);
parser.addParameter('CloseFigure', false, @(x) islogical(x) && isscalar(x));
parser.parse(varargin{:});

requested = lower(string(parser.Results.Example));
if requested == "all"
    examples = ["example91", "example92", "example93"];
else
    examples = requested;
end

if parser.Results.SaveGif && ~exist(parser.Results.OutputDir, 'dir')
    mkdir(parser.Results.OutputDir);
end

summary = struct();
for k = 1:numel(examples)
    cfg = example_config(examples(k), parser.Results.NumFrames);
    fprintf('Preparing animation for %s, t in [0, %.12g]...\n', cfg.name, cfg.tEnd);
    item = solve_animation_data(cfg);
    item.gifPath = "";
    if parser.Results.SaveGif
        item.gifPath = fullfile(parser.Results.OutputDir, cfg.fileBase + ".gif");
    end

    animate_solution(item, ...
        'SaveGif', parser.Results.SaveGif, ...
        'FrameDelay', parser.Results.FrameDelay, ...
        'PauseTime', parser.Results.PauseTime, ...
        'CloseFigure', parser.Results.CloseFigure);

    summary.(cfg.fileBase) = item;
end
end

function cfg = example_config(exampleName, numFrames)
switch lower(string(exampleName))
    case "example91"
        cfg = struct();
        cfg.name = "Example 9.1";
        cfg.fileBase = "animation_example91";
        cfg.title = "Example 9.1: q(x) = x exp(-x^2)";
        cfg.qfun = @paper_example91_initial;
        cfg.xmin = -40;
        cfg.xmax = 40;
        cfg.numGrid = 4096;
        cfg.dt = 2e-4;
        cfg.tEnd = 0.16;
        cfg.times = linspace(0, cfg.tEnd, numFrames);
        cfg.xplot = linspace(-5, 7, 1000).';
        cfg.ylims = [-1, 1];
        cfg.exact = [];
    case "example92"
        cfg = struct();
        cfg.name = "Example 9.2";
        cfg.fileBase = "animation_example92";
        cfg.title = "Example 9.2: one-soliton, c = pi";
        cfg.qfun = @(x) paper_example92_exact(x, 0, pi);
        cfg.xmin = -40;
        cfg.xmax = 40;
        cfg.numGrid = 4096;
        cfg.dt = 2e-4;
        cfg.tEnd = 1.0;
        cfg.times = linspace(0, cfg.tEnd, numFrames);
        cfg.xplot = linspace(-5, 7, 1000).';
        cfg.ylims = [-2, 0.2];
        cfg.exact = @(x, t) paper_example92_exact(x, t, pi);
    case "example93"
        cfg = struct();
        cfg.name = "Example 9.3";
        cfg.fileBase = "animation_example93";
        cfg.title = "Example 9.3: piecewise exp/Bessel data";
        cfg.qfun = @paper_example93_initial;
        cfg.xmin = -48;
        cfg.xmax = 48;
        cfg.numGrid = 4096;
        cfg.dt = 5e-5;
        cfg.tEnd = 0.03;
        cfg.times = linspace(0, cfg.tEnd, numFrames);
        cfg.xplot = linspace(-7, 7, 1000).';
        cfg.ylims = [-1, 1];
        cfg.exact = [];
    otherwise
        error('Unknown example "%s". Use "example91", "example92", "example93", or "all".', exampleName);
end
end

function item = solve_animation_data(cfg)
xgrid = linspace(cfg.xmin, cfg.xmax, cfg.numGrid + 1).';
xgrid(end) = [];
q0 = cfg.qfun(xgrid);

valuesPlot = zeros(numel(cfg.xplot), numel(cfg.times));
if isempty(cfg.exact)
    sol = solve_kdv_etdrk4(q0, xgrid, cfg.times, cfg.dt);
    for j = 1:numel(cfg.times)
        valuesPlot(:, j) = interpolate_periodic(xgrid, sol.values(:, j), cfg.xplot);
    end
else
    sol = struct('xgrid', xgrid, 'times', cfg.times, 'values', []);
    for j = 1:numel(cfg.times)
        valuesPlot(:, j) = cfg.exact(cfg.xplot, cfg.times(j));
    end
end

item = struct();
item.cfg = cfg;
item.sol = sol;
item.initialPlot = cfg.qfun(cfg.xplot);
item.valuesPlot = valuesPlot;
item.gifPath = "";
end

function animate_solution(item, varargin)
parser = inputParser;
parser.addParameter('SaveGif', true, @(x) islogical(x) && isscalar(x));
parser.addParameter('FrameDelay', 0.05, @(x) isnumeric(x) && isscalar(x) && x >= 0);
parser.addParameter('PauseTime', 0.01, @(x) isnumeric(x) && isscalar(x) && x >= 0);
parser.addParameter('CloseFigure', false, @(x) islogical(x) && isscalar(x));
parser.parse(varargin{:});

cfg = item.cfg;
fig = figure('Color', 'w', 'Position', [100, 100, 760, 460], ...
    'Name', char(cfg.title));
ax = axes(fig);
hold(ax, 'on');
plot(ax, cfg.xplot, item.initialPlot, '-', 'Color', 0.75 * [1, 1, 1], ...
    'LineWidth', 1.1, 'DisplayName', 'initial data');
currentLine = plot(ax, cfg.xplot, item.valuesPlot(:, 1), 'k-', ...
    'LineWidth', 2.0, 'DisplayName', 'u(x,t)');
hold(ax, 'off');

xlabel(ax, 'x');
ylabel(ax, 'u(x,t)');
xlim(ax, [cfg.xplot(1), cfg.xplot(end)]);
ylim(ax, cfg.ylims);
grid(ax, 'on');
box(ax, 'on');
legend(ax, 'Location', 'best');

for j = 1:numel(cfg.times)
    set(currentLine, 'YData', item.valuesPlot(:, j));
    title(ax, sprintf('%s, t = %.6g', cfg.title, cfg.times(j)), ...
        'Interpreter', 'none');
    drawnow;

    if parser.Results.SaveGif
        append_gif_frame(fig, item.gifPath, j, parser.Results.FrameDelay);
    end
    if parser.Results.PauseTime > 0
        pause(parser.Results.PauseTime);
    end
end

if parser.Results.SaveGif
    fprintf('Saved GIF:\n  %s\n', item.gifPath);
end
if parser.Results.CloseFigure
    close(fig);
end
end

function append_gif_frame(fig, gifPath, frameIndex, delayTime)
frame = getframe(fig);
[rgbImage, map] = rgb2ind(frame2im(frame), 256);
if frameIndex == 1
    imwrite(rgbImage, map, gifPath, 'gif', 'LoopCount', inf, 'DelayTime', delayTime);
else
    imwrite(rgbImage, map, gifPath, 'gif', 'WriteMode', 'append', 'DelayTime', delayTime);
end
end
