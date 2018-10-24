function [IFR, SpikeBin] = InstantFR(cfg)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

spkinds = cfg.spkinds;
Nsamples = cfg.Nsamples;
stdg = cfg.filtSD; % this standard deviation is in samples
fs = cfg.fs; % sampling frequency
starts = cfg.starts;
ends = cfg.ends;

%% Generate IFR
%t = (1:round(spkSampRate*S.RecDuration))/spkSampRate; %timeseries
SpikeBin = zeros(Nsamples, 1);
spkinds(spkinds < 1) = 1;
SpikeBin(spkinds) = 1; %binary spike vector
%sets up the gaussian filter for spike density estimation
%standard  deviation of the gaussian smoothing
filtx = -4*stdg:1:4*stdg;
filty = normpdf(filtx,0,stdg);
filty = filty/sum(filty);
%20ms gaussian smoothing
IFR = conv2(SpikeBin(:), filty(:), 'same');
%normalize to the number of spikes in the trial as in Gabbiani et al. 1999
IFR = IFR/(sum(IFR)/fs);
IFR = IFR .* sum(SpikeBin);

IFR(1:(round(fs*starts)-1)) = NaN;
IFR((round(fs*ends)+1):end) = NaN;
SpikeBin(1:(round(fs*starts)-1)) = NaN;
SpikeBin((round(fs*ends)+1):end) = NaN;

end