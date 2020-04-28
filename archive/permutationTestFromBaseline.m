function [h, p] = permutationTestFromBaseline(zdata)
% function [h, p] = permutationTestFromBaseline(zdata)
%
% zdata - trials are arranged in rows, time along the columns
% 
nperm = 1000;
%nt = size(zdata,2);
signm = [-1; 1];
%means = zeros(nperm, nt);

[h, p, ~, stats] = ttest(zdata, 0, 'alpha', .050); %operates on each column
clusts = bwlabel(h);
nclust = length(unique(clusts));
groupt = [];
for jj=1:nclust
    ts = sum([stats.tstat(clusts == jj)]);
    groupt = cat(1,groupt, ts);
end
groupt = abs(groupt);

tsums = cell(nperm, 1);
parfor ii=1:nperm
    s = signm(randi([1, 2], size(zdata,1), 1));
    sm = repmat(s, 1, size(zdata,2));
    [h1, ~, ~, stats] = ttest(zdata.*sm, 0, 'alpha', .050); %operates on each column
    clusts1 = bwlabel(h1);
    nclust = length(unique(clusts1));
    for jj=1:nclust
        ts = sum([stats.tstat(clusts1 == jj)]);
        tsums{ii} = cat(1, tsums{ii}, ts);
    end   
end
tsums = cell2mat(tsums); tsums = tsums(:);
tsums = tsums(tsums ~= 0); %get rid of the zeros in this vector
tsums = abs(tsums);

sig_thresh = quantile(tsums, .95);
hgroup = groupt >= sig_thresh; %above or below permuted significance

for ii=1:length(hgroup)
   h(clusts == ii) = hgroup(ii);
   p(clusts == ii) = p(clusts == ii);
end


% means = sort(means);
% meanz = nanmean(zdata);
% p = zeros(nt,1); h = zeros(nt,1);
% parfor ii=1:nt
%     fi = find(meanz(ii) <= means(:,ii), 1, 'first');
%     if isempty(fi) fi=nperm; end
%     p(ii) = fi/nperm;
% end
% h = p > .975 | p < .025;
% i=0;


    