function spikeTimes = readSpikeTimes(fn)
% function spikeTimes = readSpikeTimes(fn)
%
% fn is a filename in the current directory to read spike times from
% spikeTimes is a cell array containing a vector of spike times for each
% unit in the file.

fid = fopen(fn, 'r');

% spkT = textscan(fid, '%f, %f');
% Accomadate larger .txt files
% Columns 1-4: Unit,Timestamps,PC1,PC2 
% Note: Timestamps must be in second column!!
delimiter = ',';
formatSpec = '%f%f%f%f%[^\n\r]';
spkT = textscan(fid, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);

% Output Timestamps only
if length(spkT) > 1 
    spikeTimes = cat(2, spkT{1}, spkT{2});
else
    spikeTimes = [];
end

