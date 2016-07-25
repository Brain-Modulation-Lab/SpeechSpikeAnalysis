function plotSpikeWaveform(x, y, varargin)

hold on;
if nargin>2 
    ah = varargin{1};
    axes(ah);
else
    ah = gca;
end

mean_y = mean(y,1);
std_y = std(y,0,1);

lh = plot(x,mean_y,'-', 'LineWidth', 2);
plot(x,mean_y + std_y, 'LineStyle', '--', 'Color', lh.Color);
plot(x,mean_y - std_y, 'LineStyle', '--', 'Color', lh.Color);

