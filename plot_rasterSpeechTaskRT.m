function [ fh ] = plot_rasterSpeechTaskRT( figdir, figtitle, Stim, RT, SpDur, trange )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

%respWind = round(respInterval(1)*spkSampRate):round(respInterval(2)*spkSampRate);
plott = linspace(Stim.respInterval(1),Stim.respInterval(2),size(Stim.DD,2));
[fh, rasterh, meanh, ~, labelh] = makeRasterIFRplot(figtitle);

plotRastersandIFRs(plott, meanh, Stim.meanifr, Stim.stdifr, rasterh, Stim.DD, 0, [0 0 0]);
% plotRastersandIFRs(plott, meanh, mean(Stim.ISITSdata,1), std(Stim.ISITSdata,[],1)/sqrt(size(Stim.ISITSdata,1)), rasterh, Stim.DD, 0, [0 0 0]);
plot(rasterh, [0 0], [0 length(trange)], 'Color', [.4 .4 .4], 'LineWidth', 2);
plot(meanh, [0 0], meanh.YLim, 'Color', [.4 .4 .4], 'LineWidth', 2);
rasterh.XLim = [plott(1) plott(end)];
meanh.XLim = [plott(1) plott(end)];
plot(meanh, meanh.XLim, Stim.baseIFRCI(1)*[1 1], 'r');
plot(meanh, meanh.XLim, Stim.baseIFRCI(2)*[1 1], 'r');
% plot(meanh, meanh.XLim, Stim.baseISITSCI(1)*[1 1], 'r');
% plot(meanh, meanh.XLim, Stim.baseISITSCI(2)*[1 1], 'r');
plot(meanh, plott, meanh.YLim(2)*Stim.sig_excit, 'r', 'LineWidth', 3);
plot(meanh, plott, meanh.YLim(2)*Stim.sig_inhib, 'b', 'LineWidth', 3);

plot(meanh, -mean(RT)*[1 1], meanh.YLim, 'k')
plot(meanh, -(mean(RT)+std(RT))*[1 1], meanh.YLim, 'k:')
plot(meanh, -(mean(RT)-std(RT))*[1 1], meanh.YLim, 'k:')

plot(meanh, mean(SpDur)*[1 1], meanh.YLim, 'k')
plot(meanh, (mean(SpDur)+std(SpDur))*[1 1], meanh.YLim, 'k:')
plot(meanh, (mean(SpDur)-std(SpDur))*[1 1], meanh.YLim, 'k:')

saveas(fh,[figdir,'/',figtitle,'.pdf'], 'pdf');
close(fh);
end

