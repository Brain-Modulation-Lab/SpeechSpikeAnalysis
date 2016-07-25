function diff = bootMeanDiff(d1, d2)

diff = nanmean(d1, 1) - nanmean(d2, 1);

