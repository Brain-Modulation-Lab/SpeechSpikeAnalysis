function spikeTimes = readSpikeTimes(fn)
% function spikeTimes = readSpikeTimes(fn)
%
% fn is a filename in the current directory to read spike times from
% spikeTimes is a cell array containing a vector of spike times for each
% unit in the file.

fid = fopen(fn, 'r');

spkT = textscan(fid, '%f, %f');

if length(spkT) == 2 
    spikeTimes = cat(2, spkT{1}, spkT{2});
else
    spikeTimes = [];
end

