function [h] = bootstrapMeanTest(data1, data2, alpha)
% function [h] = bootstrapMeanTest(data1, data2, alpha)
%
% Test for a difference in means of the two data sets
% Data can be vectors, like firing rates, each row being a
% separate observation.  
%
% h is the result, has the same number of columns as data, 
% and indicates whether the 95% CI of the Y = mean(data1)-mean(data)
% spans 0. h = 0 means it includes 0, h = 1 means that it does not (that
% the means are different).

%meandiff = @(x, y)nanmean(x,1)-nanmean(y,1);
%diffCI = bootci(2000, {@bootMeanDiff, data1, data2}, 'type', 'cper','Options', statset('UseParallel',true), 'alpha', alpha);
nresamp = 5000;
[bootstat, ~] = bootstrp(nresamp, @nanmean, data1, 'Options', statset('UseParallel',true));
[bootstat2, ~] = bootstrp(nresamp, @nanmean, data2, 'Options', statset('UseParallel',true));
d = bootstat - bootstat2;
CI = quantile(d, [alpha/2 1-alpha/2], 1);
meand = nanmean(d);
h = CI(1,:) > 0 | CI(2,:) < 0;
meandiff = nanmean(data1) - nanmean(data2);

% Get p-values for FDR correction
% pval = zeros(1, size(data1,2));
% for ii = 1:size(data1,2)
%     ds = sort(d(:,ii));
%     pval(ii) = find(ds >= meandiff(ii), 1,'first')/nresamp;
% end
% fdr = mafdr(pval, 'BHFDR',1);    


dbg=0;
if dbg
    figure;
    t = 1:length(meand);
    plot(t, meand,'k');
    hold on;
    plot(t', CI', 'r');
end


