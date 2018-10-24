

figure;  plot3(atlases.XYZ{1,1}.mm(:,1), atlases.XYZ{1,1}.mm(:,2), atlases.XYZ{1,1}.mm(:,3), 'b.')
hold on; plot3(atlases.XYZ{2,1}.mm(:,1), atlases.XYZ{2,1}.mm(:,2), atlases.XYZ{2,1}.mm(:,3), 'r.')
hold on; plot3(atlases.XYZ{3,1}.mm(:,1), atlases.XYZ{3,1}.mm(:,2), atlases.XYZ{3,1}.mm(:,3), 'g.')
hold on; plot3(atlases.XYZ{1,2}.mm(:,1), atlases.XYZ{1,2}.mm(:,2), atlases.XYZ{1,2}.mm(:,3), 'b.')
hold on; plot3(atlases.XYZ{2,2}.mm(:,1), atlases.XYZ{2,2}.mm(:,2), atlases.XYZ{2,2}.mm(:,3), 'r.')
hold on; plot3(atlases.XYZ{3,2}.mm(:,1), atlases.XYZ{3,2}.mm(:,2), atlases.XYZ{3,2}.mm(:,3), 'g.')


T1 = delaunay(atlases.XYZ{1,1}.mm(:,1), atlases.XYZ{1,1}.mm(:,2), atlases.XYZ{1,1}.mm(:,3));
T2 = delaunay(atlases.XYZ{2,1}.mm(:,1), atlases.XYZ{2,1}.mm(:,2), atlases.XYZ{2,1}.mm(:,3));
T3 = delaunay(atlases.XYZ{3,1}.mm(:,1), atlases.XYZ{3,1}.mm(:,2), atlases.XYZ{3,1}.mm(:,3));
T4 = delaunay(atlases.XYZ{1,2}.mm(:,1), atlases.XYZ{1,2}.mm(:,2), atlases.XYZ{1,2}.mm(:,3));
T5 = delaunay(atlases.XYZ{2,2}.mm(:,1), atlases.XYZ{2,2}.mm(:,2), atlases.XYZ{2,2}.mm(:,3));
T6 = delaunay(atlases.XYZ{3,2}.mm(:,1), atlases.XYZ{3,2}.mm(:,2), atlases.XYZ{3,2}.mm(:,3));
figure;
clear alpha
%% right STN
hold on; trisurf(T1,atlases.XYZ{1,1}.mm(:,1), atlases.XYZ{1,1}.mm(:,2), atlases.XYZ{1,1}.mm(:,3),'FaceColor','b','EdgeColor','none');
alpha 0.1;
hold on; trisurf(T2,atlases.XYZ{2,1}.mm(:,1), atlases.XYZ{2,1}.mm(:,2), atlases.XYZ{2,1}.mm(:,3),'FaceColor','r','EdgeColor','none');
alpha 0.1;
hold on; trisurf(T3,atlases.XYZ{3,1}.mm(:,1), atlases.XYZ{3,1}.mm(:,2), atlases.XYZ{3,1}.mm(:,3),'FaceColor','g','EdgeColor','none');
alpha 0.1;
%% left STN
hold on; trisurf(T4,atlases.XYZ{1,2}.mm(:,1), atlases.XYZ{1,2}.mm(:,2), atlases.XYZ{1,2}.mm(:,3),'FaceColor','b','EdgeColor','none');
alpha 0.1;
hold on; trisurf(T5,atlases.XYZ{2,2}.mm(:,1), atlases.XYZ{2,2}.mm(:,2), atlases.XYZ{2,2}.mm(:,3),'FaceColor','r','EdgeColor','none');
alpha 0.1;
hold on; trisurf(T6,atlases.XYZ{3,2}.mm(:,1), atlases.XYZ{3,2}.mm(:,2), atlases.XYZ{3,2}.mm(:,3),'FaceColor','g','EdgeColor','none');
alpha 0.1;

axis equal
xlabel('M-L')
ylabel('A-P')
zlabel('D-V')


idxSNR = find([STATS(:).depth] < 0)
idxSTN = setdiff(1:length(STATS), idxSNR)


figure;
temp=cell2mat(arrayfun(@(x) x.MNIcoords', STATS(intersect(idx_excit, idxSTN)), 'uniformoutput', false))';
hold on; plot3(temp(:,1), temp(:,2), temp(:,3), 'r.', 'markersize', 25)
temp=cell2mat(arrayfun(@(x) x.MNIcoords', STATS(intersect(idx_inhib, idxSTN)), 'uniformoutput', false))';
hold on; plot3(temp(:,1), temp(:,2), temp(:,3), 'bo', 'markersize', 12)
temp=cell2mat(arrayfun(@(x) x.MNIcoords', STATS(intersect(idx_mix, idxSTN)), 'uniformoutput', false))';
hold on; plot3(temp(:,1), temp(:,2), temp(:,3), 'gx', 'markersize', 12)
temp=cell2mat(arrayfun(@(x) x.MNIcoords', STATS(intersect(idx_nr, idxSTN)), 'uniformoutput', false))';
hold on; plot3(temp(:,1), temp(:,2), temp(:,3), 'k^', 'markersize', 12)
axis equal
xlabel('M-L')
ylabel('A-P')
zlabel('D-V')
legend({'excit','inhib','mix','NR'})



figure;
idx = intersect(idxSTN, union(idxPauseSigCorrCue, idxBurstSigCorrCue));
temp=cell2mat(arrayfun(@(x) x.MNIcoords', STATS(idx), 'uniformoutput', false))';
hold on; plot3(temp(:,1), temp(:,2), temp(:,3), 'r.', 'markersize', 25)
idx = intersect(idxSTN, union(idxPauseSigCorrSpOnset, idxBurstSigCorrSpOnset));
temp=cell2mat(arrayfun(@(x) x.MNIcoords', STATS(idx), 'uniformoutput', false))';
hold on; plot3(temp(:,1), temp(:,2), temp(:,3), 'bo', 'markersize', 12)
idx = intersect(idxSTN, union(idxPauseSigCorrBoth, idxBurstSigCorrBoth)); 
temp=cell2mat(arrayfun(@(x) x.MNIcoords', STATS(idx), 'uniformoutput', false))';
hold on; plot3(temp(:,1), temp(:,2), temp(:,3), 'gx', 'markersize', 12)
idx = setdiff(idxSTN, ...
    intersect(idxSTN, [idxPauseSigCorrCue, idxBurstSigCorrCue, idxPauseSigCorrSpOnset, idxBurstSigCorrSpOnset, idxPauseSigCorrBoth, idxBurstSigCorrBoth]));
temp=cell2mat(arrayfun(@(x) x.MNIcoords', STATS(idx), 'uniformoutput', false))';
hold on; plot3(temp(:,1), temp(:,2), temp(:,3), 'k^', 'markersize', 12)
axis equal
xlabel('M-L')
ylabel('A-P')
zlabel('D-V')
legend({'speech-locked','cue-locked','both','neither'})



temp = cell2mat(arrayfun(@(x) x.MNIcoords', STATS, 'uniformoutput', false))';

class1 = {};
class1(idx_excit) = {'excit'};
class1(idx_inhib) = {'inhib'};
class1(idx_mix) = {'mix'};
class1(idx_nr) = {'nr'};
class1=class1';
class1(find(arrayfun(@(x) isempty(x.MNIcoords), STATS))) = [];

lda = fitcdiscr(temp,class1);
ldaClass = resubPredict(lda);

ldaResubErr = resubLoss(lda)

[ldaResubCM,grpOrder] = confusionmat(class1,ldaClass)

