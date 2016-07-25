function [fh, rasterh, meanh, stimh, labelh] = makeRasterIFRplot(title)

fh = figure('Position',[1 1 900 800], 'Color', [1 1 1]);
labelh = axes('Parent', fh,  'Position', [0 0 1 1], 'Visible', 'off'); hold on;
meanh = axes('Parent', fh, 'Position', [.1 .15 .8 .25], 'Color','w', 'TickDir', 'out'); hold on;
rasterh = axes('Parent', fh, 'Position', [.1 .4 .8 .4], 'Ydir', 'reverse', 'Visible', 'off'); hold on;
stimh = axes('Parent', fh, 'Position', [.1 .8 .8 .1], 'Visible', 'off'); hold on;
if(~isempty(title))
    text('parent', labelh, 'Position', [.1 .90], 'String', title, 'FontSize', 18);
end