figure; plot(mean(cell2mat(ISITSdata),2))

%figure; plot(mean(cell2mat(base),2))

hold on; plot(xlim, (mean(mean(cell2mat(ISITSbase),2),1)+icdf('Normal',0.05/38,0,1)*std(mean(cell2mat(ISITSbase),2),[],1))*[1 1], 'r')
hold on; plot(xlim, (mean(mean(cell2mat(ISITSbase),2),1)-icdf('Normal',0.05/38,0,1)*std(mean(cell2mat(ISITSbase),2),[],1))*[1 1], 'r')
hold on; plot(xlim, mean(mean(cell2mat(ISITSbase),2),1)*[1 1], 'k')


respWind = round(respInterval(1)*spkSampRate):round(respInterval(2)*spkSampRate);
plott = respWind./spkSampRate;
[fh, rasterh, meanh, ~, labelh] = makeRasterIFRplot('title');
meanifr = mean(cell2mat(IFRdata),2);
stdifr = std(cell2mat(IFRbase),[],2)./sqrt(length(trange));
plotRastersandIFRs(plott, meanh, meanifr, stdifr, rasterh, DD, 0, [0 0 0]);
plot(rasterh, [0 0], [0 length(trange)], 'Color', [.4 .4 .4], 'LineWidth', 2);
plot(meanh, [0 0], meanh.YLim, 'Color', [.4 .4 .4], 'LineWidth', 2);
rasterh.XLim = [plott(1) plott(end)];
meanh.XLim = [plott(1) plott(end)];
plot(meanh, meanh.XLim, (mean(mean(cell2mat(IFRbase),2),1)+icdf('Normal',0.05/38,0,1)*std(mean(cell2mat(IFRbase),2),[],1))*[1 1], 'r');
plot(meanh, meanh.XLim, (mean(mean(cell2mat(IFRbase),2),1)-icdf('Normal',0.05/38,0,1)*std(mean(cell2mat(IFRbase),2),[],1))*[1 1], 'r');
%plot(meanh, plott, sig, 'b', 'LineWidth', 3);

TrialData = cell2struct(TrialData, {'subject','session','depth','trials'}, 2);
SessionTable = cell2struct(SessionTable, {'subject','depth','session','wordlist'}, 2);
SpeechTrials = cell2struct(SpeechTrials, {'subject','fs','trial','AudioTrial','stimulus','target','CodingSession','response','C1error','Verror','C2error','SpeechOnset','SpeechOffset','StimulusRevealed','Comment'}, 2);