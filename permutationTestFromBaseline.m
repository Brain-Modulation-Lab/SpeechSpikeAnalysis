function [h, p] = permutationTestFromBaseline(zdata)
% function [h, p] = permutationTestFromBaseline(zdata)
%
% zdata - trials are arranged in rows, time along the columns
% 
nperm = 2000;
nt = size(zdata,2);
signm = [-1; 1];
means = zeros(nperm, nt);

parfor ii=1:nperm
    s = signm(randi([1, 2], size(zdata,1), 1));
    sm = repmat(s, 1, size(zdata,2));
    m = nanmean(zdata.*sm);
    means(ii,:) = m;
end
means = sort(means);
meanz = nanmean(zdata);
p = zeros(nt,1); h = zeros(nt,1);
parfor ii=1:nt
    fi = find(meanz(ii) <= means(:,ii), 1, 'first');
    if isempty(fi) fi=nperm; end
    p(ii) = fi/nperm;
end
h = p > .975 | p < .025;
i=0;


    