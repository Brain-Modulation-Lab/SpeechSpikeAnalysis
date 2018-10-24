idxSNR = find([STATS(:).depth] < 0)

%% Percent depth -- Response type (excit, inhib, mix, nr) -- scatter
figure; hold on;
idxCent = find(arrayfun(@(x) strcmp(x.elec, 'Cent'), STATS));
idxCent = setdiff(idxCent, idxSNR);
idx = intersect(idxCent, idx_excit);
plot3(zeros(length(idx),1), zeros(length(idx),1), 100*[STATS(idx).depth], 'r.', 'markersize', 25)
idx = intersect(idxCent, idx_inhib);
plot3(zeros(length(idx),1), zeros(length(idx),1), 100*[STATS(idx).depth], 'bo', 'markersize', 12)
idx = intersect(idxCent, idx_mix);
plot3(zeros(length(idx),1), zeros(length(idx),1), 100*[STATS(idx).depth], 'gx', 'markersize', 12)
idx = intersect(idxCent, idx_nr);
plot3(zeros(length(idx),1), zeros(length(idx),1), 100*[STATS(idx).depth], 'k^', 'markersize', 12)
idxPost = find(arrayfun(@(x) strcmp(x.elec, 'Post'), STATS));
idxPost = setdiff(idxPost, idxSNR);
idx = intersect(idxPost, idx_excit);
plot3(zeros(length(idx),1), -ones(length(idx),1), 100*[STATS(idx).depth], 'r.', 'markersize', 25)
idx = intersect(idxPost, idx_inhib);
plot3(zeros(length(idx),1), -ones(length(idx),1), 100*[STATS(idx).depth], 'bo', 'markersize', 12)
idx = intersect(idxPost, idx_mix);
plot3(zeros(length(idx),1), -ones(length(idx),1), 100*[STATS(idx).depth], 'gx', 'markersize', 12)
idx = intersect(idxPost, idx_nr);
plot3(zeros(length(idx),1), -ones(length(idx),1), 100*[STATS(idx).depth], 'k^', 'markersize', 12)
idxMed = find(arrayfun(@(x) strcmp(x.elec, 'Med'), STATS));
idxMed = setdiff(idxMed, idxSNR);
idx = intersect(idxMed, idx_excit);
plot3(ones(length(idx),1), zeros(length(idx),1), 100*[STATS(idx).depth], 'r.', 'markersize', 25)
idx = intersect(idxMed, idx_inhib);
plot3(ones(length(idx),1), zeros(length(idx),1), 100*[STATS(idx).depth], 'bo', 'markersize', 12)
idx = intersect(idxMed, idx_mix);
plot3(ones(length(idx),1), zeros(length(idx),1), 100*[STATS(idx).depth], 'gx', 'markersize', 12)
idx = intersect(idxMed, idx_nr);
plot3(ones(length(idx),1), zeros(length(idx),1), 100*[STATS(idx).depth], 'k^', 'markersize', 12)

xlabel('M-L')
ylabel('A-P')
zlabel('D-V')
legend({'excit','inhib','mix','NR'})
xlim([-0.2 1])
ylim([-1 0])
thisAxis = gca;
thisAxis.XTick = [];
thisAxis.YTick = [];


%% Percent depth -- Response type (excit, inhib, mix, nr) -- median, IQR
figure; hold on;
idxCent = find(arrayfun(@(x) strcmp(x.elec, 'Cent'), STATS));
idxCent = setdiff(idxCent, idxSNR);
idx = intersect(idxCent, idx_excit);
plot3(0, 0, 100*median([STATS(idx).depth]), 'r.', 'markersize', 10)
plot3([0 0], [0 0], 100*quantile([STATS(idx).depth], [0.25 0.75]), 'r.', 'markersize', 10)
idx = intersect(idxCent, idx_inhib);
plot3(0, 0, 100*median([STATS(idx).depth]), 'b.', 'markersize', 10)
plot3([0 0], [0 0], 100*quantile([STATS(idx).depth], [0.25 0.75]), 'b.', 'markersize', 10)
idx = intersect(idxCent, idx_mix);
plot3(0, 0, 100*median([STATS(idx).depth]), 'g.', 'markersize', 10)
plot3([0 0], [0 0], 100*quantile([STATS(idx).depth], [0.25 0.75]), 'g.', 'markersize', 10)
idx = intersect(idxCent, idx_nr);
plot3(0, 0, 100*median([STATS(idx).depth]), 'k.', 'markersize', 10)
plot3([0 0], [0 0], 100*quantile([STATS(idx).depth], [0.25 0.75]), 'k.', 'markersize', 10)
idxPost = find(arrayfun(@(x) strcmp(x.elec, 'Post'), STATS));
idxPost = setdiff(idxPost, idxSNR);
idx = intersect(idxPost, idx_excit);
plot3(0, -1, 100*median([STATS(idx).depth]), 'r.', 'markersize', 10)
plot3([0 0], [-1 -1], 100*quantile([STATS(idx).depth], [0.25 0.75]), 'r.', 'markersize', 10)
idx = intersect(idxPost, idx_inhib);
plot3(0, -1, 100*median([STATS(idx).depth]), 'b.', 'markersize', 10)
plot3([0 0], [-1 -1], 100*quantile([STATS(idx).depth], [0.25 0.75]), 'b.', 'markersize', 10)
idx = intersect(idxPost, idx_mix);
plot3(0, -1, 100*median([STATS(idx).depth]), 'g.', 'markersize', 10)
plot3([0 0], [-1 -1], 100*quantile([STATS(idx).depth], [0.25 0.75]), 'g.', 'markersize', 10)
idx = intersect(idxPost, idx_nr);
plot3(0, -1, 100*median([STATS(idx).depth]), 'k.', 'markersize', 10)
plot3([0 0], [-1 -1], 100*quantile([STATS(idx).depth], [0.25 0.75]), 'k.', 'markersize', 10)
idxMed = find(arrayfun(@(x) strcmp(x.elec, 'Med'), STATS));
idxMed = setdiff(idxMed, idxSNR);
idx = intersect(idxMed, idx_excit);
plot3(1, 0, 100*median([STATS(idx).depth]), 'r.', 'markersize', 10)
plot3([1 1], [0 0], 100*quantile([STATS(idx).depth], [0.25 0.75]), 'r.', 'markersize', 10)
idx = intersect(idxMed, idx_inhib);
plot3(1, 0, 100*median([STATS(idx).depth]), 'b.', 'markersize', 10)
plot3([1 1], [0 0], 100*quantile([STATS(idx).depth], [0.25 0.75]), 'b.', 'markersize',10)
idx = intersect(idxMed, idx_mix);
plot3(1, 0, 100*median([STATS(idx).depth]), 'g.', 'markersize', 10)
plot3([1 1], [0 0], 100*quantile([STATS(idx).depth], [0.25 0.75]), 'g.', 'markersize', 10)
idx = intersect(idxMed, idx_nr);
plot3(1, 0, 100*median([STATS(idx).depth]), 'k.', 'markersize', 10)
plot3([1 1], [0 0], 100*quantile([STATS(idx).depth], [0.25 0.75]), 'k.', 'markersize', 10)

xlabel('M-L')
ylabel('A-P')
zlabel('D-V')
legend({'excit','inhib','mix','NR'})
xlim([-0.2 1])
ylim([-1 0])
thisAxis = gca;
thisAxis.XTick = [];
thisAxis.YTick = [];


DispCamPos.cp = campos;
DispCamPos.cva = camva;
DispCamPos.ct = camtarget;
DispCamPos.uv = camup;

set(gca,'CameraPosition',DispCamPos.cp,...
'CameraTarget',DispCamPos.ct,...
'CameraViewAngle',DispCamPos.cva,...
'CameraUpVector',DispCamPos.uv);


%% Percent depth -- Locking type (cue, speech onset, nr) -- scatter
% Cent

alpha=0.05;
SortCat = union(union(Asorts,Bsorts),Csorts);

idxPauseSigCorrSpOnset = find(arrayfun(@(x) ...
    x.PauseSpOnset.(CorrType).pval<=alpha&&x.CuePause.(CorrType).pval>alpha, DATA));
idxPauseSigCorrCue = find(arrayfun(@(x) ...
    x.CuePause.(CorrType).pval<=alpha&&x.PauseSpOnset.(CorrType).pval>alpha, DATA));
idxPauseSigCorrBoth = find(arrayfun(@(x) ...
    x.CuePause.(CorrType).pval<=alpha&&x.PauseSpOnset.(CorrType).pval<=alpha, DATA));
idxPauseSigCorrSpOnset = intersect(idxPauseSigCorrSpOnset, union(idx_inhib_Cue, idx_mix_Cue));
idxPauseSigCorrSpOnset = intersect(idxPauseSigCorrSpOnset, SortCat);
idxPauseSigCorrCue = intersect(idxPauseSigCorrCue, union(idx_inhib_SpOnset, idx_mix_SpOnset));
idxPauseSigCorrCue = intersect(idxPauseSigCorrCue, SortCat);
idxPauseSigCorrBoth = intersect(idxPauseSigCorrBoth, union(union(idx_inhib_Cue, idx_inhib_SpOnset), union(idx_mix_Cue, idx_mix_SpOnset)));
idxPauseSigCorrBoth = intersect(idxPauseSigCorrBoth, SortCat);

idxBurstSigCorrSpOnset = find(arrayfun(@(x) ...
    x.BurstSpOnset.(CorrType).pval<=alpha&&x.CueBurst.(CorrType).pval>alpha, DATA));
idxBurstSigCorrCue = find(arrayfun(@(x) ...
    x.CueBurst.(CorrType).pval<=alpha&&x.BurstSpOnset.(CorrType).pval>alpha, DATA));
idxBurstSigCorrBoth = find(arrayfun(@(x) ...
    x.CueBurst.(CorrType).pval<=alpha&&x.BurstSpOnset.(CorrType).pval<=alpha, DATA));
idxBurstSigCorrSpOnset = intersect(idxBurstSigCorrSpOnset, union(idx_excit_Cue, idx_mix_Cue));
idxBurstSigCorrSpOnset = intersect(idxBurstSigCorrSpOnset, SortCat);
idxBurstSigCorrCue = intersect(idxBurstSigCorrCue, union(idx_excit_SpOnset, idx_mix_SpOnset));
idxBurstSigCorrCue = intersect(idxBurstSigCorrCue, SortCat);
idxBurstSigCorrBoth = intersect(idxBurstSigCorrBoth, union(union(idx_excit_Cue, idx_excit_SpOnset), union(idx_mix_Cue, idx_mix_SpOnset)));
idxBurstSigCorrBoth = intersect(idxBurstSigCorrBoth, SortCat);

figure; hold on;
idxCent = find(arrayfun(@(x) strcmp(x.elec, 'Cent'), STATS));
idxCent = setdiff(idxCent, idxSNR);
idx = intersect(idxCent, union(idxPauseSigCorrCue, idxBurstSigCorrCue));
plot3(zeros(length(idx),1), zeros(length(idx),1), 100*[STATS(idx).depth], 'r.', 'markersize', 25)
fprintf('Central track: Speech onset-locked n = %d\n', length(idx));
idx = intersect(idxCent, union(idxPauseSigCorrSpOnset, idxBurstSigCorrSpOnset));
plot3(zeros(length(idx),1), zeros(length(idx),1), 100*[STATS(idx).depth], 'bo', 'markersize', 12)
fprintf('Central track: Cue-locked n = %d\n', length(idx));
idx = intersect(idxCent, union(idxPauseSigCorrBoth, idxBurstSigCorrBoth));
plot3(zeros(length(idx),1), zeros(length(idx),1), 100*[STATS(idx).depth], 'gx', 'markersize', 12)
fprintf('Central track: Speech onset- and Cue-locked n = %d\n', length(idx));
idx = setdiff(idxCent, ...
    intersect(idxCent, [idxPauseSigCorrCue, idxBurstSigCorrCue, idxPauseSigCorrSpOnset, idxBurstSigCorrSpOnset, idxPauseSigCorrBoth, idxBurstSigCorrBoth]));
plot3(zeros(length(idx),1), zeros(length(idx),1), 100*[STATS(idx).depth], 'k^', 'markersize', 12)
fprintf('Central track: No significant locking n = %d\n', length(idx));
% Post
idxPost = find(arrayfun(@(x) strcmp(x.elec, 'Post'), STATS));
idxPost = setdiff(idxPost, idxSNR);
idx = intersect(idxPost, union(idxPauseSigCorrCue, idxBurstSigCorrCue));
plot3(zeros(length(idx),1), -ones(length(idx),1), 100*[STATS(idx).depth], 'r.', 'markersize', 25)
fprintf('Post track: Speech onset-locked n = %d\n', length(idx));
idx = intersect(idxPost, union(idxPauseSigCorrSpOnset, idxBurstSigCorrSpOnset));
plot3(zeros(length(idx),1), -ones(length(idx),1), 100*[STATS(idx).depth], 'bo', 'markersize', 12)
fprintf('Post track: Cue-locked n = %d\n', length(idx));
idx = intersect(idxPost, union(idxPauseSigCorrBoth, idxBurstSigCorrBoth));
plot3(zeros(length(idx),1), -ones(length(idx),1), 100*[STATS(idx).depth], 'gx', 'markersize', 12)
fprintf('Post track: Speech onset- and Cue-locked n = %d\n', length(idx));
idx = setdiff(idxPost, ...
    intersect(idxPost, [idxPauseSigCorrCue, idxBurstSigCorrCue, idxPauseSigCorrSpOnset, idxBurstSigCorrSpOnset, idxPauseSigCorrBoth, idxBurstSigCorrBoth]));
plot3(zeros(length(idx),1), -ones(length(idx),1), 100*[STATS(idx).depth], 'k^', 'markersize', 12)
fprintf('Post track: No significant locking n = %d\n', length(idx));
% Med
idxMed = find(arrayfun(@(x) strcmp(x.elec, 'Med'), STATS));
idxMed = setdiff(idxMed, idxSNR);
idx = intersect(idxMed, union(idxPauseSigCorrCue, idxBurstSigCorrCue));
plot3(ones(length(idx),1), zeros(length(idx),1), 100*[STATS(idx).depth], 'r.', 'markersize', 25)
fprintf('Medial track: Speech onset-locked n = %d\n', length(idx));
idx = intersect(idxMed, union(idxPauseSigCorrSpOnset, idxBurstSigCorrSpOnset));
plot3(ones(length(idx),1), zeros(length(idx),1), 100*[STATS(idx).depth], 'bo', 'markersize', 12)
fprintf('Medial track: Cue-locked n = %d\n', length(idx));
idx = intersect(idxMed, union(idxPauseSigCorrBoth, idxBurstSigCorrBoth));
plot3(ones(length(idx),1), zeros(length(idx),1), 100*[STATS(idx).depth], 'gx', 'markersize', 12)
fprintf('Medial track: Speech onset- and Cue-locked n = %d\n', length(idx));
idx = setdiff(idxMed, ...
    intersect(idxMed, [idxPauseSigCorrCue, idxBurstSigCorrCue, idxPauseSigCorrSpOnset, idxBurstSigCorrSpOnset, idxPauseSigCorrBoth, idxBurstSigCorrBoth]));
plot3(ones(length(idx),1), zeros(length(idx),1), 100*[STATS(idx).depth], 'k^', 'markersize', 12)
fprintf('Medial track: No significant locking n = %d\n', length(idx));

fprintf('SNR: Speech onset-locked n = %d\n', length(intersect(idxSNR, union(idxPauseSigCorrCue, idxBurstSigCorrCue))));
fprintf('SNR: Cue-locked n = %d\n', length(intersect(idxSNR, union(idxPauseSigCorrSpOnset, idxBurstSigCorrSpOnset))));
fprintf('SNR: Speech onset- and Cue-locked n = %d\n', length(intersect(idxSNR, union(idxPauseSigCorrBoth, idxBurstSigCorrBoth))));

xlabel('M-L')
ylabel('A-P')
zlabel('D-V')
legend({'speech onset','cue','both','NR'})
xlim([-0.2 1])
ylim([-1 0])
thisAxis = gca;
thisAxis.XTick = [];
thisAxis.YTick = [];

%% Percent depth -- Locking type (cue, speech onset, nr) --  median, IQR
% Cent
figure; hold on;
idxCent = find(arrayfun(@(x) strcmp(x.elec, 'Cent'), STATS));
idxCent = setdiff(idxCent, idxSNR);
idx = intersect(idxCent, union(idxPauseSigCorrCue, idxBurstSigCorrCue));
plot3(0, 0, 100*median([STATS(idx).depth]), 'r.', 'markersize', 10)
plot3([0 0], [0 0], 100*quantile([STATS(idx).depth], [0.25 0.75]), 'r.', 'markersize', 10)
idx = intersect(idxCent, union(idxPauseSigCorrSpOnset, idxBurstSigCorrSpOnset));
plot3(0, 0, 100*median([STATS(idx).depth]), 'b.', 'markersize', 10)
plot3([0 0], [0 0], 100*quantile([STATS(idx).depth], [0.25 0.75]), 'b.', 'markersize', 10)
idx = intersect(idxCent, union(idxPauseSigCorrBoth, idxBurstSigCorrBoth));
plot3(0, 0, 100*median([STATS(idx).depth]), 'g.', 'markersize', 10)
plot3([0 0], [0 0], 100*quantile([STATS(idx).depth], [0.25 0.75]), 'g.', 'markersize', 10)
idx = setdiff(idxCent, ...
    intersect(idxCent, [idxPauseSigCorrCue, idxBurstSigCorrCue, idxPauseSigCorrSpOnset, idxBurstSigCorrSpOnset, idxPauseSigCorrBoth, idxBurstSigCorrBoth]));
plot3(0, 0, 100*median([STATS(idx).depth]), 'k.', 'markersize', 10)
plot3([0 0], [0 0], 100*quantile([STATS(idx).depth], [0.25 0.75]), 'k.', 'markersize', 10)
% Post
idxPost = find(arrayfun(@(x) strcmp(x.elec, 'Post'), STATS));
idxPost = setdiff(idxPost, idxSNR);
idx = intersect(idxPost, union(idxPauseSigCorrCue, idxBurstSigCorrCue));
plot3(0, -1, 100*median([STATS(idx).depth]), 'r.', 'markersize', 10)
plot3([0 0], [-1 -1], 100*quantile([STATS(idx).depth], [0.25 0.75]), 'r.', 'markersize', 10)
idx = intersect(idxPost, union(idxPauseSigCorrSpOnset, idxBurstSigCorrSpOnset));
plot3(0, -1, 100*median([STATS(idx).depth]), 'b.', 'markersize', 10)
plot3([0 0], [-1 -1], 100*quantile([STATS(idx).depth], [0.25 0.75]), 'b.', 'markersize', 10)
idx = intersect(idxPost, union(idxPauseSigCorrBoth, idxBurstSigCorrBoth));
plot3(0, -1, 100*median([STATS(idx).depth]), 'g.', 'markersize', 10)
plot3([0 0], [-1 -1], 100*quantile([STATS(idx).depth], [0.25 0.75]), 'g.', 'markersize', 10)
idx = setdiff(idxPost, ...
    intersect(idxPost, [idxPauseSigCorrCue, idxBurstSigCorrCue, idxPauseSigCorrSpOnset, idxBurstSigCorrSpOnset, idxPauseSigCorrBoth, idxBurstSigCorrBoth]));
plot3(0, -1, 100*median([STATS(idx).depth]), 'k.', 'markersize', 10)
plot3([0 0], [-1 -1], 100*quantile([STATS(idx).depth], [0.25 0.75]), 'k.', 'markersize', 10)
% Med
idxMed = find(arrayfun(@(x) strcmp(x.elec, 'Med'), STATS));
idxMed = setdiff(idxMed, idxSNR);
idx = intersect(idxMed, union(idxPauseSigCorrCue, idxBurstSigCorrCue));
plot3(1, 0, 100*median([STATS(idx).depth]), 'r.', 'markersize', 10)
plot3([1 1], [0 0], 100*quantile([STATS(idx).depth], [0.25 0.75]), 'r.', 'markersize', 10)
idx = intersect(idxMed, union(idxPauseSigCorrSpOnset, idxBurstSigCorrSpOnset));
plot3(1, 0, 100*median([STATS(idx).depth]), 'b.', 'markersize', 10)
plot3([1 1], [0 0], 100*quantile([STATS(idx).depth], [0.25 0.75]), 'b.', 'markersize', 10)
idx = intersect(idxMed, union(idxPauseSigCorrBoth, idxBurstSigCorrBoth));
plot3(1, 0, 100*median([STATS(idx).depth]), 'g.', 'markersize', 10)
plot3([1 1], [0 0], 100*quantile([STATS(idx).depth], [0.25 0.75]), 'g.', 'markersize', 10)
idx = setdiff(idxMed, ...
    intersect(idxMed, [idxPauseSigCorrCue, idxBurstSigCorrCue, idxPauseSigCorrSpOnset, idxBurstSigCorrSpOnset, idxPauseSigCorrBoth, idxBurstSigCorrBoth]));
plot3(1, 0, 100*median([STATS(idx).depth]), 'k.', 'markersize', 10)
plot3([1 1], [0 0], 100*quantile([STATS(idx).depth], [0.25 0.75]), 'k.', 'markersize', 10)

xlabel('M-L')
ylabel('A-P')
zlabel('D-V')
legend({'speech onset','cue','both','NR'})
xlim([-0.2 1])
ylim([-1 0])
thisAxis = gca;
thisAxis.XTick = [];
thisAxis.YTick = [];



