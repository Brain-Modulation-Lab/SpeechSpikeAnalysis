% iter = 100;
% allBase = zeros(length(trial(1).baselineWindow), nTrials, iter);
% for ii=1:nTrials
%     allBase(:,ii,1) = IFR(trial(ii).baselineWindow);
%     for jj=2:iter
%         allBase(:,ii,jj) = circshift(IFR(trial(ii).baselineWindow), randi(length(trial(ii).baselineWindow)));
%     end
% end
% 
% allBase_mean = squeeze(mean(allBase, 2));
% alpha = .05;
% zBound  = [prctile(allBase_mean(:),alpha/2), prctile(allBase_mean(:),100*(1-alpha/2))];

for ll=1:Unit.nTrials
    trialTimeEventTimes = (Unit.trial(ll).EventTimes(4) - Unit.trial(ll).EventTimes(1));
    trialTimeEventInds = (Unit.trial(ll).EventInds(4) - Unit.trial(ll).EventInds(1))/Unit.spkSampRate;
    recLen = length(Unit.trial(ll).IFRZ)/Unit.spkSampRate;
    baseLen = length(Unit.trial(ll).baselineWindow)/Unit.spkSampRate;
    postLen = (Unit.trial(ll).EventInds(5)-Unit.trial(ll).EventInds(4))/Unit.spkSampRate;
    trialIntervals = diff(Unit.trial(ll).EventTimes);
    %fprintf('EventTime: %d  EventInds: %d  recLen: %d baselineLen: %d postLen: %d\n',  trialTimeEventTimes, trialTimeEventInds, recLen, baseLen, postLen)
    fprintf('Trial Intervals: %f %f %f %f\n', trialIntervals);
end

% figure; hold on;
%     for j = 1:length(spikeInd)
%         y = V(electrodeNum, spikeWind+spikeInd(j)) - V(electrodeNum, spikeInd(j));
%         line(spikeWind./sampRate, y, 'Color', 'k');
%     end

% figure; hold on;
% plot(EventTimes1'*[1 1], [0 1], 'k');
% for ll=1:Unit.nTrials
%     plot(Unit.trial(ll).EventTimes(4)*[1,1], [1 2], 'r');
%     plot(Unit.tSpk(Unit.trial(ll).EventInds(4))*[1 1], [2 3], 'b');
% end

% for ll=1:119
%     delays(ll) = Unit.trial(ll).EventTimes(4) - Unit.trial(ll).EventTimes(1);
% end
% 
% hold on;
% for ll=1:119
%     plot(tAudioTrials(:,ll), AudioEnvTrials(:,ll)./max(AudioEnvTrials(:,ll))*10+40, 'r');
%     plot(Unit.tSpk(Unit.trial(ll).EventInds(1):Unit.trial(ll).EventInds(5)), Unit.IFR(Unit.trial(ll).EventInds(1):Unit.trial(ll).EventInds(5)), 'k');
% end
%plot(EventTimes1'*[1 1], ylim, 'b');