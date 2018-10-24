% % get cue locked responses from excit and inhib
% idx_allCueLocked = idxSigCorrSpOnset
% idx_allCueLocked = [idx_allCueLocked, idxSigCorrSpOnset]
% 
% % get speech locked responses from excit and inhib
% idx_allSpeechLocked = idxSigCorrCue
% idx_allSpeechLocked = [idx_allSpeechLocked, idxSigCorrCue]

temp1 = [];

ResponseType = 'PP';

for idx = idx_allCueLocked
   
    time_exists = find(arrayfun(@(x) ~isempty(x.t), DATA(idx).Trial.(ResponseType)));
    
    %temp1 = [temp1, mean(arrayfun(@(x) x.t(1), DATA(idx).Trial.(ResponseType)(time_exists)))];
    temp1 = [temp1, mean(arrayfun(@(x,y) x.t(1)-y, DATA(idx).Trial.(ResponseType)(time_exists),DATA(idx).RT(DATA(idx).trange(time_exists))'))];
end

temp2 = [];

ResponseType = 'BB';

for idx = idx_allSpeechLocked
   
    time_exists = find(arrayfun(@(x) ~isempty(x.t), DATA(idx).Trial.(ResponseType)));
    
    %temp2 = [temp2, mean(arrayfun(@(x) x.t(1), DATA(idx).Trial.(ResponseType)(time_exists)))];
    temp2 = [temp2, mean(arrayfun(@(x,y) x.t(1)-y, DATA(idx).Trial.(ResponseType)(time_exists),DATA(idx).RT(DATA(idx).trange(time_exists))'))];
    
end

mean(temp1)
std(temp1)/sqrt(length(temp1))
mean(temp2)
std(temp2)/sqrt(length(temp2))

[h, p] = ttest2(temp1, temp2)

[hist_lat1, t0] = hist(temp1, -1.25:.125:.25);
[hist_lat2, t0] = hist(temp2, -1.25:.125:.25);

figure; bar(t0, [hist_lat1; hist_lat2]')