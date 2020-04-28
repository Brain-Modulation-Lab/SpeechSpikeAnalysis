function [ isits ] = ISIts(D, fs)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

I = isi(D);

isits = zeros(size(D));
k = find(D);

for i=1:length(I)
    isits(k(i):(k(i)+I(i))) = I(i);
end

isits = isits/fs;
isits(isnan(D)) = NaN;

end

