RespType = cell(size(STATS));
RespType(idx_excit) = {'excit'};
RespType(idx_inhib) = {'inhib'};
RespType(idx_mix) = {'mix'};
RespType(idx_nr) = {'nr'};

LockingType = cell(size(STATS));
LockingType(union(idxPauseSigCorrCue, idxBurstSigCorrCue)) = {'speech'};
LockingType(union(idxPauseSigCorrSpOnset, idxBurstSigCorrSpOnset)) = {'cue'};
LockingType(union(idxPauseSigCorrBoth, idxBurstSigCorrBoth)) = {'both'};
LockingType(find(cellfun(@isempty, LockingType))) = {'nr'};

[~,~,stats] = anovan([STATS(:).depth],{{STATS(:).elec} RespType LockingType},'model','interaction',...
    'varnames',{'track','RespType','LockingType'});

idxSTN = setdiff(1:length(STATS), idxSNR);

[~,~,stats] = anovan([STATS(idxSTN).depth],{{STATS(idxSTN).elec} RespType(idxSTN)},'model','interaction',...
    'varnames',{'track','RespType'});

multcompare(stats, 'dimension', [1 2])

p = kruskalwallis([STATS(intersect(idxCent, idxSTN)).depth], RespType(intersect(idxCent, idxSTN)))
p = kruskalwallis([STATS(intersect(idxPost, idxSTN)).depth], RespType(intersect(idxPost, idxSTN)))
p = kruskalwallis([STATS(intersect(idxMed, idxSTN)).depth], RespType(intersect(idxMed, idxSTN)))

p = kruskalwallis([STATS(intersect(idxCent, idxSTN)).depth], LockingType(intersect(idxCent, idxSTN)))
p = kruskalwallis([STATS(intersect(idxPost, idxSTN)).depth], LockingType(intersect(idxPost, idxSTN)))
p = kruskalwallis([STATS(intersect(idxMed, idxSTN)).depth], LockingType(intersect(idxMed, idxSTN)))