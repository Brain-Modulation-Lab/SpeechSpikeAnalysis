function [ fh ] = plot_rasterSpeechTaskBurstsOrder( figdir, figtitle, Stim, trange, BB, thresh, BurstColor, order, markers, MarkerColor )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

%respWind = round(respInterval(1)*spkSampRate):round(respInterval(2)*spkSampRate);
plott = linspace(Stim.respInterval(1),Stim.respInterval(2),size(Stim.DD,2));
[fh, rasterh, meanh, ~, labelh] = makeRasterIFRplot(figtitle);

if isempty(BurstColor)
    colormap parula;
    cm = colormap+(1-colormap)*0.5;
    ncolor = size(cm,1);
    mbb = max([BB.S]);
    BurstColor = arrayfun(@(x) cm(find(x<=linspace(0,mbb,ncolor),1,'first'),:), ...
        [BB(:).S],...
        'uniformoutput', false);
else
    BurstColor = arrayfun(@(x) BurstColor, 1:size(Stim.DD,1), ...
        'uniformoutput', false);
end

k=0;
for i = order'%1:size(Stim.DD,1)
    k=k+1;
    if ~isempty(BB(i).t) && BB(i).S > thresh
        rectangle('Parent',rasterh,'Position', ...
            [BB(i).t(1), ...
            k, ...
            BB(i).t(end)-BB(i).t(1), ...
            0.9], ...
            'FaceColor',BurstColor{i},'EdgeColor',BurstColor{i});
% max FR within burst
%         line('Parent',rasterh,'XData',BB(i).tmax*[1 1],...
%               'YData',[k (k+0.9)], 'LineWidth', 3, 'color', 'r');
    end
    
    if ~isempty(markers)
        line('Parent',rasterh,'XData',markers(i)*[1 1],...
              'YData',[k (k+0.9)], 'LineWidth', 3, 'color', MarkerColor);
    end
end

plotRastersandIFRs(plott, meanh, Stim.meanifr, Stim.stdifr, rasterh, Stim.DD(order,:), 0, [0 0 0]);

plot(rasterh, [0 0], [0 length(trange)], 'Color', [.4 .4 .4], 'LineWidth', 2);
plot(meanh, [0 0], meanh.YLim, 'Color', [.4 .4 .4], 'LineWidth', 2);
rasterh.XLim = [plott(1) plott(end)];
meanh.XLim = [plott(1) plott(end)];
plot(meanh, meanh.XLim, Stim.baseIFRCI(1)*[1 1], 'r');
plot(meanh, meanh.XLim, Stim.baseIFRCI(2)*[1 1], 'r');
plot(meanh, plott, meanh.YLim(2)*Stim.sig_excit, 'r', 'LineWidth', 3);
plot(meanh, plott, meanh.YLim(2)*Stim.sig_inhib, 'b', 'LineWidth', 3);

% saveas(fh,[figdir,'/',figtitle,'.pdf'], 'pdf');
% close(fh);
end

