function [ color ] = getColor( f, frange, cm )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

bins = (frange(1):(frange(2)-frange(1))/size(cm,1):frange(2));
colorID = max(1, sum(f > bins(1:size(cm,1)))); 

color = cm(colorID, :);
end

