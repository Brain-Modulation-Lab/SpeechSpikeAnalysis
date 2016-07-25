function [spkVects, closei] = removeDuplicateSpikes(spkTimes, tThresh, SU)
% function spkVects = removeDuplicateSpikes(spkTimes, tThresh, SU)
%
% Spikes are removed from MUs if they are within almost symultanteous w a
% SU spike, since it is likely a double trigger of the same waveform.
%
% spkTimes - nx2 vector of unit numbers and spkTime doublets
% tThresh - the smallest time difference that is considered ok. ISIs between
% units smaller than it will be removed
% SU - the single units of the recording. These are used as the template 

nSU = length(SU);
%tThresh = .0005;

if isempty(SU) %arbitary
    SU = 1;
end

isi = diff(spkTimes(:,2));

closei = find(isi < tThresh);
disp(sprintf('Found %d ISI violations in recording', length(closei)));
%removeCount = 0;
removei = [];
for ii = 1:length(closei)
    units = spkTimes(closei:(closei+1), 1);
    if units(1) == SU
        removei = cat(1, removei, closei+1);
    else 
        removei = cat(1, removei, closei);
    end
end
   
keep = true(size(spkTimes,1),1);
keep(removei) = 0;

spkVects = spkTimes(keep,:);
