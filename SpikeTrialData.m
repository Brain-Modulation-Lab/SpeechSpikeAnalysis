function [ Stim ] = SpikeTrialData( spkSampRate, IFR, ISITS, D, trange, basetimes, trialtimes, respInterval, respIntervalDisp, alpha, Nobs, minTsig )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

Stim.respInterval = respIntervalDisp;

%% Baseline data (hard-coded to take one second back from basetimes)
Stim.base_samples = round(spkSampRate*basetimes);
Stim.IFRbase = cell2mat(arrayfun(@(x) IFR((Stim.base_samples(x)-spkSampRate):Stim.base_samples(x)), ...
    trange, 'uniformoutput', false))';
Stim.ISITSbase = cell2mat(arrayfun(@(x) ISITS((Stim.base_samples(x)-spkSampRate):Stim.base_samples(x)), ...
    trange, 'uniformoutput', false))';

%% Trial data
Stim.trial_samples = round(spkSampRate*trialtimes);
Stim.DD = cell2mat(arrayfun(@(x) ...
    D((Stim.trial_samples(x)+round(spkSampRate*Stim.respInterval(1))): ...
    (Stim.trial_samples(x)+round(spkSampRate*Stim.respInterval(2)))), ...
    trange, 'uniformoutput', false))';
Stim.IFRdata = cell2mat(arrayfun(@(x) ...
    IFR((Stim.trial_samples(x)+round(spkSampRate*Stim.respInterval(1))): ...
    (Stim.trial_samples(x)+round(spkSampRate*Stim.respInterval(2)))), ...
    trange, 'uniformoutput', false))';
Stim.ISITSdata = cell2mat(arrayfun(@(x) ...
    ISITS((Stim.trial_samples(x)+round(spkSampRate*Stim.respInterval(1))): ...
    (Stim.trial_samples(x)+round(spkSampRate*Stim.respInterval(2)))), ...
    trange, 'uniformoutput', false))';

%% only used for statistical purposes (Nspike calc below)
DD = cell2mat(arrayfun(@(x) ...
    D((Stim.trial_samples(x)+round(spkSampRate*respInterval(1))): ...
    (Stim.trial_samples(x)+round(spkSampRate*respInterval(2)))), ...
    trange, 'uniformoutput', false))';
% IFRdata = cell2mat(arrayfun(@(x) ...
%     IFR((Stim.trial_samples(x)+round(spkSampRate*Stim.respIntervalDisp(1))): ...
%     (Stim.trial_samples(x)+round(spkSampRate*Stim.respInterval(2)))), ...
%     trange, 'uniformoutput', false))';
% ISITSdata = cell2mat(arrayfun(@(x) ...
%     ISITS((Stim.trial_samples(x)+round(spkSampRate*Stim.respIntervalDisp(1))): ...
%     (Stim.trial_samples(x)+round(spkSampRate*Stim.respInterval(2)))), ...
%     trange, 'uniformoutput', false))';

%% significance thresholds based on baseline mean and stdev corrected for the 
%% "numbner of observations" within the response window
%
%  WJL 7/23/2017 edit: change Nobs for ISITS to the mean number of ISIs
%  within the interval
%
Stim.basemeanIFR = mean(Stim.IFRbase(:));
Stim.baseIFRCI = Stim.basemeanIFR + ...
    [ icdf('Normal',alpha/Nobs,0,1)*std(mean(Stim.IFRbase,1),[],2) ...
    -icdf('Normal',alpha/Nobs,0,1)*std(mean(Stim.IFRbase,1),[],2) ];
Nisi = round(sum(mean(DD,1)))-1;
Stim.basemeanISITS = mean(Stim.ISITSbase(:));
Stim.baseISITSCI = Stim.basemeanISITS + ...
    [ icdf('Normal',alpha/Nisi,0,1)*std(mean(Stim.ISITSbase,1),[],2) ...
    -icdf('Normal',alpha/Nisi,0,1)*std(mean(Stim.ISITSbase,1),[],2) ];

%% mean and std error for display purposes
Stim.meanifr = mean(Stim.IFRdata,1);
Stim.stdifr = std(Stim.IFRdata,[],1)./sqrt(length(trange));

%% identify timepoints with significant increases and decreases in firing
Stim.sig_excit = round(Stim.meanifr>Stim.baseIFRCI(2));
Psig = bwconncomp(Stim.sig_excit);
for i = 1:Psig.NumObjects
    if length(Psig.PixelIdxList{i})<minTsig
        Stim.sig_excit(Psig.PixelIdxList{i})=0;
    end
end
Stim.sig_excit(Stim.sig_excit==0)=NaN;

Stim.sig_inhib = round(mean(Stim.ISITSdata,1)>Stim.baseISITSCI(2));
Psig = bwconncomp(Stim.sig_inhib);
for i = 1:Psig.NumObjects
    if length(Psig.PixelIdxList{i})<minTsig
        Stim.sig_inhib(Psig.PixelIdxList{i})=0;
    end
end
Stim.sig_inhib(Stim.sig_inhib==0)=NaN;

end

