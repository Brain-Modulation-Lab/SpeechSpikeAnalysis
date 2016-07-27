function unit = UnitResponse(t, spkTimes, vSampRate, spkSampRate, filtSD, EventTimes, SkipEvents, nEventsPerTrial, trialRange)
% function unit = UnitResponse(t, spkTimes, vSampRate, spkSampRate, filtSD, Events, nEventsPerTrial)

dbg = 0;
down = round(vSampRate/spkSampRate);
spkSampRate = vSampRate/down;
tSpk = downsample(t, down);
EventInds = round(EventTimes * spkSampRate);

D = zeros(length(tSpk), 1);
D(round(spkTimes.*spkSampRate)) = 1; %binary spike vector

%sets up the gaussian filter for spike density estimation
%standard  deviation of the gaussian smoothing
stdg = filtSD; %40 samples is a gaussian with std 40ms
filtx = -4*stdg:1:4*stdg;
filty = normpdf(filtx,0,stdg);
filty = filty/sum(filty);

%20ms gaussian smoothing
IFR = conv2(D(:), filty(:), 'same');
%normalize to the number of spikes in the trial as in Gabbiani et al. 1999
IFR = IFR/(sum(IFR)/spkSampRate);
IFR = IFR .* sum(D);

if dbg
    colors = {'r','b','m','c'};
    figure; hold on;
end
%nTrials = floor((length(EventTimes) - SkipEvents)./nEventsPerTrial);
unit.trialRange = trialRange;
trialRange = trialRange(1):trialRange(2); %changing to be a vector
nTrials = trialRange(end) - trialRange(1) + 1;
for ii = 1:nTrials
    stim = trialRange(ii);
    eventI = (ii)*nEventsPerTrial + SkipEvents;
    trial(ii).EventInds = EventInds((eventI-3):(eventI+1));
    trial(ii).EventTimes = EventTimes((eventI-3):(eventI+1));
    trialWind = trial(ii).EventInds(1):trial(ii).EventInds(5);
    baseWind = (trial(ii).EventInds(4)-1*spkSampRate):(trial(ii).EventInds(4)); %defines the baseline period
    trial(ii).baselineIFR = nanmean(IFR(baseWind));
    trial(ii).basestdIFR = nanstd(IFR(baseWind));
    trial(ii).baselineSpkCount(ii) = sum(D(baseWind));
    trial(ii).baselineCountFR = trial(ii).baselineSpkCount.*(length(baseWind)/spkSampRate);
    trial(ii).IFRZ = (IFR(trialWind)-trial(ii).baselineIFR)./trial(ii).basestdIFR;
    trial(ii).baselineWindow = baseWind;
    trial(ii).stim = stim;
    if (sum(D(EventInds((eventI-3):(eventI+1)))) >=1)
        trial(ii).hasSpikes = true;
    else
        trial(ii).hasSpikes = false;
    end
    
    if dbg
        plot(trial(ii).EventTimes'*[1 1], [-10 0], 'Color', colors{mod(ii-1,4)+1});
        plot(tSpk(trial(ii).EventInds(4):trial(ii).EventInds(5)), IFR(trial(ii).EventInds(4):trial(ii).EventInds(5)), 'k');
    end
end

% Compute a resampled baseline firing rate statistic that gives z-scored
% bounds for the firing rate range outside of which a response should be
% considered as significant. Bounds + mean could be used to z-score trials.
iter = 100;
allBase = zeros(length(trial(1).baselineWindow), nTrials, iter);
for ii=1:nTrials
    allBase(:,ii,1) = IFR(trial(ii).baselineWindow);
    for jj=2:iter
        allBase(:,ii,jj) = circshift(IFR(trial(ii).baselineWindow), randi(length(trial(ii).baselineWindow)));
    end
end
allBase_mean = squeeze(mean(allBase, 2));
alpha = .05;
zBound  = [prctile(allBase_mean(:),alpha/2), prctile(allBase_mean(:),100*(1-alpha/2))];
zMean = mean(allBase_mean);
%allBase = [trial.baselineIFR];
% [f, x] = ecdf(allBase(:));
% zlowi = find(f < .05, 1, 'last');
% zhighi = find(f >= .95, 1, 'first');
% zBound = [x(zlowi) x(zhighi)];

% fill in the unit structure
unit.D = D;
unit.IFR = IFR;
unit.tSpk = tSpk;
unit.spkSampRate = spkSampRate;
unit.trial = trial;
unit.EventInds = EventInds;
unit.nTrials = nTrials;
unit.zBound = zBound;
unit.zMid = zMean;









