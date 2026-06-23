function summary = animate_sech110_istm(varargin)
%ANIMATE_SECH110_ISTM Animate q(x)=-110 sech^2(x) by reflectionless ISTM.
%
% The initial data
%
%   q(x) = -110 sech^2(x) = -10*11 sech^2(x)
%
% is the reflectionless Poschl-Teller potential with ten bound states.
% For the KdV equation u_t - 6*u*u_x + u_xxx = 0, its scattering data are
%   tau_j = j, lambda_j = -j^2, j = 1,...,10,
% with zero reflection.  The animation below evolves these data by the ISTM
% exponential law and reconstructs u(x,t) from the reflectionless Marchenko
% determinant / Hirota tau function.
%
% Examples:
%   animate_sech110_istm
%   animate_sech110_istm('TFinal', 0.6, 'XRange', [-30, 270])

thisDir = fileparts(mfilename('fullpath'));

parser = inputParser;
parser.addParameter('TFinal', 0.50, @(x) isnumeric(x) && isscalar(x) && x > 0);
parser.addParameter('NumFrames', 600, @(x) isnumeric(x) && isscalar(x) && x >= 2);
parser.addParameter('XRange', [-25, 230], @(x) isnumeric(x) && numel(x) == 2 && x(1) < x(2));
parser.addParameter('NumX', 12000, @(x) isnumeric(x) && isscalar(x) && x >= 200);
parser.addParameter('OutputDir', fullfile(thisDir, 'animations_sech110_istm'), ...
    @(s) ischar(s) || isstring(s));
parser.addParameter('SaveGif', false, @(x) islogical(x) && isscalar(x));
parser.addParameter('SaveMp4', true, @(x) islogical(x) && isscalar(x));
parser.addParameter('SaveMat', true, @(x) islogical(x) && isscalar(x));
parser.addParameter('MakeWaterfall', true, @(x) islogical(x) && isscalar(x));
parser.addParameter('FrameDelay', 0.035, @(x) isnumeric(x) && isscalar(x) && x > 0);
parser.addParameter('FrameRate', 60, @(x) isnumeric(x) && isscalar(x) && x > 0);
parser.addParameter('VideoQuality', 95, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 100);
parser.addParameter('FigurePosition', [60, 60, 1600, 900], ...
    @(x) isnumeric(x) && numel(x) == 4 && all(x(3:4) > 0));
parser.addParameter('ShowFigure', true, @(x) islogical(x) && isscalar(x));
parser.parse(varargin{:});

opts = parser.Results;
outdir = char(opts.OutputDir);
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

N = 10;
tau = 1:N;
alpha = poschl_teller_norming_constants(tau);
data = reflectionless_precompute(tau, alpha);

x = linspace(opts.XRange(1), opts.XRange(2), round(opts.NumX)).';
times = linspace(0, opts.TFinal, round(opts.NumFrames));
storeSolution = opts.SaveMat || opts.MakeWaterfall || nargout > 0;
if storeSolution
    u = zeros(numel(x), numel(times));
else
    u = [];
end

fprintf('Animating ISTM evolution of -110 sech^2(x).\n');
fprintf('Discrete data: tau_j = 1,...,10, lambda_j = -tau_j^2, reflection = 0.\n');
fprintf('Computing %d frames on x in [%g,%g].\n', numel(times), x(1), x(end));

if storeSolution
    for j = 1:numel(times)
        u(:, j) = reflectionless_solution(x, times(j), data);
    end
end

% Validate the t=0 reconstruction on a focused grid where the initial data
% is not visually compressed by the large animation window.
xCheck = linspace(-8, 8, 2001).';
uCheck = reflectionless_solution(xCheck, 0, data);
qExact = -110 ./ cosh(xCheck).^2;
t0MaxError = max(abs(uCheck - qExact));

gifFile = fullfile(outdir, 'sech110_istm_10soliton.gif');
mp4File = fullfile(outdir, 'sech110_istm_10soliton.mp4');
waterfallFile = fullfile(outdir, 'sech110_istm_waterfall.png');
matFile = fullfile(outdir, 'sech110_istm_data.mat');
summaryFile = fullfile(outdir, 'sech110_istm_summary.txt');

if opts.SaveGif || opts.SaveMp4
    if storeSolution
        write_animation_files(x, times, u, opts, gifFile, mp4File);
    else
        write_animation_files_streamed(x, times, data, opts, gifFile, mp4File);
    end
end
if opts.MakeWaterfall
    write_waterfall_plot(x, times, u, waterfallFile);
end
if opts.SaveMat
    save(matFile, 'x', 'times', 'u', 'tau', 'alpha', 't0MaxError', 'opts');
end
write_summary(summaryFile, tau, alpha, opts, t0MaxError, gifFile, mp4File, waterfallFile, matFile);

summary = struct();
summary.x = x;
summary.times = times;
summary.u = u;
summary.tau = tau;
summary.alpha = alpha;
summary.t0MaxError = t0MaxError;
summary.gifFile = gifFile;
summary.mp4File = mp4File;
summary.waterfallFile = waterfallFile;
summary.matFile = matFile;
summary.summaryFile = summaryFile;

fprintf('t=0 ISTM reconstruction max error on [-8,8]: %.3e\n', t0MaxError);
fprintf('Wrote summary to:\n  %s\n', summaryFile);
if opts.SaveGif
    fprintf('Wrote GIF to:\n  %s\n', gifFile);
end
if opts.SaveMp4
    fprintf('Wrote MP4 to:\n  %s\n', mp4File);
end
if opts.MakeWaterfall
    fprintf('Wrote waterfall plot to:\n  %s\n', waterfallFile);
end
end

function alpha = poschl_teller_norming_constants(tau)
% Norming constants for -N(N+1)sech^2(x) in the convention used here.
tau = tau(:).';
alpha = zeros(size(tau));
for j = 1:numel(tau)
    others = tau;
    others(j) = [];
    alpha(j) = 2 * tau(j) * prod(abs((tau(j) + others) ./ (tau(j) - others)));
end
end

function data = reflectionless_precompute(tau, alpha)
tau = tau(:).';
alpha = alpha(:).';
N = numel(tau);
numTerms = 2^N;

masks = false(numTerms, N);
for bit = 1:N
    block = 2^(bit - 1);
    pattern = repmat([false(block, 1); true(block, 1)], 2^(N - bit), 1);
    masks(:, bit) = pattern;
end

logBeta = log(alpha ./ (2 * tau));
base = masks * logBeta(:);
for i = 1:N
    for j = (i + 1):N
        pairLog = 2 * log(abs((tau(i) - tau(j)) / (tau(i) + tau(j))));
        base = base + double(masks(:, i) & masks(:, j)) * pairLog;
    end
end

p = -2 * masks * tau(:);
r = 8 * masks * (tau(:).^3);

data = struct();
data.tau = tau;
data.alpha = alpha;
data.base = base;
data.p = p;
data.p2 = p.^2;
data.r = r;
end

function u = reflectionless_solution(x, t, data)
% Stable log-sum-exp evaluation of u=-2 d_x^2 log(tau_function).
x = x(:);
exponents = x * data.p.' + t * data.r.' + data.base.';
shift = max(exponents, [], 2);
weights = exp(exponents - shift);
normalizer = sum(weights, 2);
meanP = (weights * data.p) ./ normalizer;
meanP2 = (weights * data.p2) ./ normalizer;
u = -2 * (meanP2 - meanP.^2);
u = real(u);
end

function write_animation_files(x, times, u, opts, gifFile, mp4File)
if opts.ShowFigure
    fig = figure('Color', 'w', 'Position', opts.FigurePosition);
else
    fig = figure('Color', 'w', 'Visible', 'off', 'Position', opts.FigurePosition);
end
set(fig, 'Renderer', 'opengl');
ax = axes(fig);
mainLine = plot(ax, x, u(:, 1), 'k-', 'LineWidth', 1.7);
hold(ax, 'on');
zeroLine = plot(ax, [x(1), x(end)], [0, 0], 'Color', [0.68, 0.68, 0.68], 'LineWidth', 0.8);
initialLine = plot(ax, x, u(:, 1), 'Color', [0.55, 0.55, 0.55], 'LineStyle', ':', 'LineWidth', 1.0);
hold(ax, 'off');

xlim(ax, [x(1), x(end)]);
ylim(ax, [-225, 12]);
xlabel(ax, 'x');
ylabel(ax, 'u(x,t)');
title(ax, 'ISTM evolution of -110 sech^2(x)');
set(ax, 'Box', 'on', 'FontSize', 11);
grid(ax, 'on');

label = text(ax, x(1) + 0.02 * (x(end) - x(1)), -18, '', ...
    'FontSize', 11, 'Color', 'k', 'VerticalAlignment', 'top');
legend(ax, [mainLine, initialLine, zeroLine], ...
    {'ISTM solution', 'initial profile', 'zero level'}, 'Location', 'southeast');

video = [];
if opts.SaveMp4
    try
        video = VideoWriter(mp4File, 'MPEG-4');
        video.FrameRate = opts.FrameRate;
        video.Quality = opts.VideoQuality;
        open(video);
    catch err
        warning('Could not open MP4 writer "%s": %s', mp4File, err.message);
        video = [];
    end
end

for j = 1:numel(times)
    set(mainLine, 'YData', u(:, j));
    set(label, 'String', sprintf('t = %.4g,  ten solitons: speeds 4,16,...,400', times(j)));
    title(ax, sprintf('ISTM evolution of -110 sech^2(x), frame %d/%d', j, numel(times)));
    drawnow;

    frame = getframe(fig);
    if ~isempty(video)
        writeVideo(video, frame);
    end
    if opts.SaveGif
        [im, map] = rgb2ind(frame2im(frame), 256);
        if j == 1
            imwrite(im, map, gifFile, 'gif', 'LoopCount', Inf, 'DelayTime', opts.FrameDelay);
        else
            imwrite(im, map, gifFile, 'gif', 'WriteMode', 'append', 'DelayTime', opts.FrameDelay);
        end
    end
end

if ~isempty(video)
    close(video);
end
close(fig);
end

function write_animation_files_streamed(x, times, data, opts, gifFile, mp4File)
if opts.ShowFigure
    fig = figure('Color', 'w', 'Position', opts.FigurePosition);
else
    fig = figure('Color', 'w', 'Visible', 'off', 'Position', opts.FigurePosition);
end
set(fig, 'Renderer', 'opengl');

u0 = reflectionless_solution(x, times(1), data);
ax = axes(fig);
mainLine = plot(ax, x, u0, 'k-', 'LineWidth', 1.9);
hold(ax, 'on');
zeroLine = plot(ax, [x(1), x(end)], [0, 0], 'Color', [0.68, 0.68, 0.68], 'LineWidth', 0.8);
initialLine = plot(ax, x, u0, 'Color', [0.55, 0.55, 0.55], 'LineStyle', ':', 'LineWidth', 1.0);
hold(ax, 'off');

xlim(ax, [x(1), x(end)]);
ylim(ax, [-225, 12]);
xlabel(ax, 'x');
ylabel(ax, 'u(x,t)');
title(ax, 'ISTM evolution of -110 sech^2(x)');
set(ax, 'Box', 'on', 'FontSize', 12);
grid(ax, 'on');

label = text(ax, x(1) + 0.02 * (x(end) - x(1)), -18, '', ...
    'FontSize', 12, 'Color', 'k', 'VerticalAlignment', 'top');
legend(ax, [mainLine, initialLine, zeroLine], ...
    {'ISTM solution', 'initial profile', 'zero level'}, 'Location', 'southeast');

video = [];
if opts.SaveMp4
    try
        video = VideoWriter(mp4File, 'MPEG-4');
        video.FrameRate = opts.FrameRate;
        video.Quality = opts.VideoQuality;
        open(video);
    catch err
        warning('Could not open MP4 writer "%s": %s', mp4File, err.message);
        video = [];
    end
end

for j = 1:numel(times)
    if j == 1
        uj = u0;
    else
        uj = reflectionless_solution(x, times(j), data);
    end
    set(mainLine, 'YData', uj);
    set(label, 'String', sprintf('t = %.4g,  ten solitons: speeds 4,16,...,400', times(j)));
    title(ax, sprintf('ISTM evolution of -110 sech^2(x), frame %d/%d', j, numel(times)));
    drawnow;

    frame = getframe(fig);
    if ~isempty(video)
        writeVideo(video, frame);
    end
    if opts.SaveGif
        [im, map] = rgb2ind(frame2im(frame), 256);
        if j == 1
            imwrite(im, map, gifFile, 'gif', 'LoopCount', Inf, 'DelayTime', opts.FrameDelay);
        else
            imwrite(im, map, gifFile, 'gif', 'WriteMode', 'append', 'DelayTime', opts.FrameDelay);
        end
    end
end

if ~isempty(video)
    close(video);
end
close(fig);
end

function write_waterfall_plot(x, times, u, waterfallFile)
fig = figure('Color', 'w', 'Visible', 'off', 'Position', [100, 100, 980, 620]);
ax = axes(fig);
imagesc(ax, x, times, u.');
set(ax, 'YDir', 'normal', 'FontSize', 11, 'Box', 'on');
xlabel(ax, 'x');
ylabel(ax, 't');
title(ax, 'ISTM ten-soliton evolution from -110 sech^2(x)');
colormap(ax, turbo(256));
caxis(ax, [-205, 5]);
cb = colorbar(ax);
ylabel(cb, 'u(x,t)');
try
    exportgraphics(fig, waterfallFile, 'Resolution', 220);
catch
    print(fig, waterfallFile, '-dpng', '-r220');
end
close(fig);
end

function write_summary(summaryFile, tau, alpha, opts, t0MaxError, gifFile, mp4File, waterfallFile, matFile)
fid = fopen(summaryFile, 'w');
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, 'ISTM animation for q(x) = -110 sech^2(x)\n\n');
fprintf(fid, 'Equation: u_t - 6*u*u_x + u_xxx = 0\n');
fprintf(fid, 'Scattering data are reflectionless.\n');
fprintf(fid, 'tau_j: ');
fprintf(fid, '%g ', tau);
fprintf(fid, '\n');
fprintf(fid, 'lambda_j: ');
fprintf(fid, '%g ', -tau.^2);
fprintf(fid, '\n');
fprintf(fid, 'alpha_j:\n');
fprintf(fid, '  %.16g\n', alpha);
fprintf(fid, '\n');
fprintf(fid, 'TFinal: %.16g\n', opts.TFinal);
fprintf(fid, 'NumFrames: %d\n', round(opts.NumFrames));
fprintf(fid, 'XRange: [%.16g, %.16g]\n', opts.XRange(1), opts.XRange(2));
fprintf(fid, 'NumX: %d\n', round(opts.NumX));
fprintf(fid, 'dx: %.16g\n', diff(opts.XRange) / (round(opts.NumX) - 1));
fprintf(fid, 'FrameRate: %.16g\n', opts.FrameRate);
fprintf(fid, 'FigurePosition: [%.16g, %.16g, %.16g, %.16g]\n', opts.FigurePosition);
fprintf(fid, 't=0 max reconstruction error on [-8,8]: %.6e\n\n', t0MaxError);
if opts.SaveGif
    fprintf(fid, 'GIF: %s\n', gifFile);
end
if opts.SaveMp4
    fprintf(fid, 'MP4: %s\n', mp4File);
end
if opts.MakeWaterfall
    fprintf(fid, 'Waterfall plot: %s\n', waterfallFile);
end
if opts.SaveMat
    fprintf(fid, 'MAT data: %s\n', matFile);
end
end
