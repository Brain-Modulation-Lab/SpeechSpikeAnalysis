function [h, p] = clusterPermuteTtest(A, B)
% function [h, p] = clusterPermuteTtest(data1, data2)
% data matrices have time along the x, and trial along the y
plotdistb = 1;

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

tsums = tsums(tsums ~= 0); %get rid of the zeros in this vector
sig_thresh = quantile(tsums, [.025 1-.025]);
hgroup = (groupt <= sig_thresh(1)) | (groupt >= sig_thresh(2)); %above or below permuted significance

for ii=1:length(hgroup)
   h(clusts == ii) = hgroup(ii);
   p(clusts == ii) = p(clusts == ii);
end

if plotdistb && ~isempty(tsums)
    figure;
    [f, stepx] = ecdf(tsums);
    stairs(stepx,f, 'k'); hold on;
    ah = gca;
    for ii=1:length(groupt)
        plot([groupt(ii), groupt(ii)], ah.YLim,'r');
    end
end

    

