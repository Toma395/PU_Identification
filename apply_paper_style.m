function apply_paper_style(fig, varargin)
% Apply consistent paper/presentation style to all axes in a figure.
% Call BEFORE saveas, AFTER all plot commands.
%
% Usage:
%   apply_paper_style(fig)
%   apply_paper_style(fig, 'font_size',14, 'title_size',16, 'line_width',2, 'marker_size',8)

p = inputParser;
addParameter(p, 'font_size',   14);
addParameter(p, 'title_size',  16);
addParameter(p, 'line_width',   2.0);
addParameter(p, 'marker_size',  8);
parse(p, varargin{:});
opt = p.Results;

fig.Color = 'w';

all_ax = findall(fig, 'Type', 'axes');
for i = 1:numel(all_ax)
    ax = all_ax(i);
    % Skip transparent overlay axes (Color='none')
    if isequal(ax.Color, 'none'); continue; end

    ax.Color          = 'w';
    ax.FontSize       = opt.font_size;
    ax.LineWidth      = 1.0;
    ax.GridColor      = [0.70 0.70 0.70];
    ax.GridAlpha      = 0.80;
    ax.MinorGridColor = [0.80 0.80 0.80];
    ax.MinorGridAlpha = 0.50;

    if isgraphics(ax.Title,'text') && ~isempty(ax.Title.String)
        ax.Title.FontSize = opt.title_size;
    end
    ax.XLabel.FontSize = opt.font_size;
    ax.YLabel.FontSize = opt.font_size;

    % Raise LineWidth and MarkerSize of all Line objects
    lines = findall(ax, 'Type', 'line');
    for j = 1:numel(lines)
        l = lines(j);
        if l.LineWidth < opt.line_width
            l.LineWidth = opt.line_width;
        end
        if ~strcmp(l.Marker,'none') && l.MarkerSize < opt.marker_size
            l.MarkerSize = opt.marker_size;
        end
    end
end

% Legend styling (all legends in the figure)
all_leg = findall(fig, 'Type', 'Legend');
for i = 1:numel(all_leg)
    all_leg(i).Color   = 'w';
    all_leg(i).Box     = 'on';
    all_leg(i).FontSize = max(all_leg(i).FontSize, opt.font_size - 2);
end

% Force all text elements to near-black (prevents gray-on-white from dark theme reversal)
dark = [0.15 0.15 0.15];
for i = 1:numel(all_ax)
    ax = all_ax(i);
    if isequal(ax.Color, 'none'); continue; end
    ax.XColor = dark;
    ax.YColor = dark;
    ax.Title.Color  = dark;
    ax.XLabel.Color = dark;
    ax.YLabel.Color = dark;
end
for i = 1:numel(all_leg)
    all_leg(i).TextColor = dark;
end
all_cb = findall(fig, 'Type', 'ColorBar');
for i = 1:numel(all_cb)
    all_cb(i).Color = dark;
end
end
