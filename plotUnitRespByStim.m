%Compare the responses of two sets of trials
respInterval = [-2 2.5];
colorlist = {[.3 .6 1], [1 0 0], [0 .8 0], [0 0 0],}; 

stimMat = L{wl};
block = MerData(rec_idx(wl));
% grouping by first consonant type, lips v tongue
consType = stimMat(:,4);
consType = cellfun(@min, consType);
trialFeature = consType;
features = unique(consType);
groups = {[1 3], [4 5]};
% grouping by 2nd consonant type, lips v tongue
% consType = stimMat(:,5);
% consType = cellfun(@min, consType);
% trialFeature = consType;
% features = unique(consType);
% groups = {[1 3], [4 5]};
%grouping by early v late consonant
% early = cell2mat(cellfun(@(x) ~isempty(x), stimMat(:,3), 'UniformOutput', false));
% late = ~early;
% trialFeature = late;
% features = [0 1];
% groups = {0, 1};
% grouping by word/nonword
% word = [stimMat{:,2}];
% nonword = ~word;
% trialFeature = word;
% features = [0 1];
% groups = {0, 1};

kk = 1;

[fh, rasterh, meanh, ~, labelh] = makeRasterIFRplot(sprintf('%s,Word List %d,%s electrode,Unit %d', subjectName, wl, electrodeList{trode}, n));
respWind = round(respInterval(1)*unit.spkSampRate):round(respInterval(2)*unit.spkSampRate);
spkm = NaN*zeros(unit.nTrials, length(respWind)); 
ifrm = NaN*zeros(size(spkm));

corrTrials = ~isnan(block.ResponseTimes); nGroupTrials = 0;
for ll = 1:length(groups) %loop through groups, plot them separately
    % Just putting together groups from the set of features 
    groupTrials = false(s,1); nGroupTrials(ll) = 0;
    for mm = 1:length(groups{ll}) % loop through features in a group
        featTrials = (trialFeature == groups{ll}(mm) & corrTrials(:));
        groupTrials = groupTrials(:) | featTrials(:);
    end
    
    trialList = find(groupTrials);
    unitTrials = unit.trialRange(1):unit.trialRange(2);
    [~,overlapTrials] = ismember(trialList, unitTrials);
    overlapTrials = overlapTrials(overlapTrials ~= 0);
    nGroupTrials(ll) = length(overlapTrials);
    for mm=1:nGroupTrials(ll)
        t = overlapTrials(mm);
        audioInd = unit.trial(t).EventInds(4) + round(Rec(rec_idx(wl)).AudioStart{unit.trial(t).stim}/sampRate*unit.spkSampRate);
        offset = unit.trial(t).EventInds(1);
        ifrm(kk,:) = unit.IFR(audioInd+respWind);
        spkm(kk,:) = unit.D(audioInd+respWind);
        kk = kk+1;
    end
    meanifr = nanmean(ifrm); stdifr = nanstd(ifrm)./sqrt(length(overlapTrials));
    dataRange(ll,:) = [(kk-nGroupTrials(ll)), (kk-1)];
    meanfunc = @(x)nanmean(x,1);
    semfunc = @(x) nanstd(x,0,1);
    if nGroupTrials(ll) ~= 0
        [ifrCI, bootstat] = bootci(1000, {meanfunc, ifrm(dataRange(ll,1):dataRange(ll,2),:)}, 'type', 'cper','Options', statset('UseParallel',true));
        %the SEM of the bootstrapped mean is the Standard Dev of the resampled means (Kass, Eden, Brown 2014; Ch 9)
        ifrsem = nanstd(bootstat,0,1); 
        plott = respWind./unit.spkSampRate;
        plotRastersandIFRs(plott, meanh, meanifr, ifrsem, rasterh, spkm((kk-nGroupTrials(ll)):sum(nGroupTrials), :), kk-nGroupTrials(ll), colorlist{ll});
        %plotRastersandIFRs(plott, meanh, meanifr, ifrCI(2,:)-meanifr, rasterh, spkm((kk-nGroupTrials(ll)):sum(nGroupTrials), :), kk-nGroupTrials(ll), colorlist{ll});
    end
end
% Annotate the plot
rasterh.XLim = [plott(1) plott(end)];
meanh.XLim = [plott(1) plott(end)];
plot(meanh, meanh.XLim, [unit.zBound(1) unit.zBound(1)], 'r');
plot(meanh, meanh.XLim, [unit.zBound(2) unit.zBound(2)], 'r');
xlabel(meanh, 'Time relative to speech onset (sec)', 'FontSize', 16);
ylabel(meanh, 'Average Firing Rate (Hz)', 'FontSize', 16);
plot(rasterh, [0 0], [0 unit.nTrials], 'Color', [.7 .7 .7], 'LineWidth', 2); 
plot(meanh, [0 0], meanh.YLim, 'Color', [.7 .7 .7], 'LineWidth', 2);

if ~sum(nGroupTrials == 0)
%     h = bootstrapMeanTest(ifrm(dataRange(1,1):dataRange(1,2),:), ifrm(dataRange(2,1):dataRange(2,2),:), .01);
%     sig = h*(meanh.YLim(2)-1); sig(sig==0)=NaN;
%     plot(meanh, plott, sig, 'k', 'LineWidth', 3);

    [h, p] = clusterPermuteTtest(ifrm(dataRange(1,1):dataRange(1,2),:), ifrm(dataRange(2,1):dataRange(2,2),:));
    sig = h*(meanh.YLim(2)-1); sig(sig==0)=NaN;
    plot(meanh, plott, sig, 'b', 'LineWidth', 3);
end