function [h, p] = clusterPermuteTtest(A, B)
% function [h, p] = clusterPermuteTtest(data1, data2)
% data matrices have time along the x, and trial along the y

nrep = 1000;
ntrials = [size(A, 1), size(B,1)];
%mintrials = min(ntrials);

[h, p, ~, stats] = ttest2(A, B); %operates on each column
clusts = bwlabel(h);
nclust = length(unique(clusts));
groupt = [];
for jj=1:nclust
        ts = sum([stats.tstat(clusts == jj)]);
        groupt = cat(1,groupt, ts);
end

C = cat(1, A, B); %make a single pool of all trials to permute
nrows = size(C,1);

tsums = [];
for ii=1:nrep
    order = randperm(nrows);
    A1 = C(order(1:ntrials(1)), :);
    bi = order((ntrials(1)+1):nrows);
    B1 = C(bi, :);
    
    [h1, ~, ~, stat1] = ttest2(A1,B1);
    clust1 = bwlabel(h1);
    nclust = length(unique(clust1));
    for jj=1:nclust
        ts = sum([stat1.tstat(clust1 == jj)]);
        tsums = cat(1,tsums, ts);
    end
end