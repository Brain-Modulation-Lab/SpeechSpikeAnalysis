Asorts = find(strcmp({DATA(:).grade}, 'A'));

Bsorts = find(strcmp({DATA(:).grade}, 'B'));

Csorts = find(strcmp({DATA(:).grade}, 'C'));

BsortsSU = intersect(Bsorts,...
    find(strcmp({DATA(:).RecType}, 'SU')));

BsortsMU = intersect(Bsorts,...
    find(strcmp({DATA(:).RecType}, 'MU')));

CsortsSU = intersect(Csorts,...
    find(strcmp({DATA(:).RecType}, 'SU')));

CsortsMU = intersect(Csorts,...
    find(strcmp({DATA(:).RecType}, 'MU')));


STATS = [];

%SortCat = Asorts;
SortCat = union(union(Asorts,Bsorts),Csorts);
%SortCat = union(Asorts,Bsorts);
%ElecOrder = {{'Cent'}, {'Post'}, {'Med', 'Lat'}};

idx_excit = [];
idx_inhib = [];
idx_mix = [];
idx_nr = [];

alignment = 'SpeechOnsetBaseline';

for i=1:length(SortCat)
    
    STATS(i).SubjectID = DATA(SortCat(i)).SubjectID;
    %STATS(i).depth = DATA(SortCat(i)).depth;
    STATS(i).session = DATA(SortCat(i)).session;
    STATS(i).RecType = DATA(SortCat(i)).RecType;
    STATS(i).grade = DATA(SortCat(i)).grade;
    
    %eln = find(cellfun(@(x) max(strcmp(DATA(SortCat(i)).elec, x)), ElecOrder));
    %STATS(i).depth = TrialData(DATA(SortCat(i)).idx_trials).MicroNormDepth{eln};
    
%     idxLocs = find(strcmp(DATA(SortCat(i)).SubjectID,{Locs(:).newID}));
%     
%     if ~isempty(idxLocs)
%         
%         if ~isempty(intersect(...
%             find(cellfun(@(x) strcmp(x, 'Top border'), {Locs(idxLocs).mermarkers(:).markertype})),...
%             find(cellfun(@(x) ~isempty(strfind(x,lower(DATA(SortCat(i)).elec))), {Locs(idxLocs).mermarkers(:).tract}))...
%             ))
%             topBorder = Locs(idxLocs).mermarkers(intersect(...
%                 find(cellfun(@(x) strcmp(x, 'Top border'), {Locs(idxLocs).mermarkers(:).markertype})),...
%                 find(cellfun(@(x) ~isempty(strfind(x,lower(DATA(SortCat(i)).elec))), {Locs(idxLocs).mermarkers(:).tract}))...
%                 )).depth;
%         else
%             fprintf('Top Border not found: %s %s\n', DATA(SortCat(i)).SubjectID, DATA(SortCat(i)).elec);
%         end
%         if ~isempty(intersect(...
%             find(cellfun(@(x) strcmp(x, 'Bottom border'), {Locs(idxLocs).mermarkers(:).markertype})),...
%             find(cellfun(@(x) ~isempty(strfind(x,lower(DATA(SortCat(i)).elec))), {Locs(idxLocs).mermarkers(:).tract}))...
%             ))
%             bottomBorder = Locs(idxLocs).mermarkers(intersect(...
%                 find(cellfun(@(x) strcmp(x, 'Bottom border'), {Locs(idxLocs).mermarkers(:).markertype})),...
%                 find(cellfun(@(x) ~isempty(strfind(x,lower(DATA(SortCat(i)).elec))), {Locs(idxLocs).mermarkers(:).tract}))...
%                 )).depth;
%         else
%             fprintf('Bottom Border not found: %s %s\n', DATA(SortCat(i)).SubjectID, DATA(SortCat(i)).elec);
%         end
%         if ~isempty(find([Locs(idxLocs).mermarkers(:).depth]==DATA(SortCat(i)).depth))
%             MNIcoords = Locs(idxLocs).mermarkers(intersect(...
%                 find([Locs(idxLocs).mermarkers(:).depth]==DATA(SortCat(i)).depth),...
%                 find(cellfun(@(x) ~isempty(strfind(x,lower(DATA(SortCat(i)).elec))), {Locs(idxLocs).mermarkers(:).tract}))...
%                 )).coords_mm;
%             STATS(i).depth = 1-(topBorder-DATA(SortCat(i)).depth)/(topBorder-bottomBorder);
%             STATS(i).MNIcoords = MNIcoords;
%         else
%             fprintf('Depth not found: %s %s %f\n', DATA(SortCat(i)).SubjectID, DATA(SortCat(i)).elec, DATA(SortCat(i)).depth);
%             STATS(i).depth = [];
%             STATS(i).MNIcoords = [];
%         end
%         
%     else
%         fprintf('lead-dbs for %s not found.\n', DATA(SortCat(i)).SubjectID);
%         STATS(i).depth = [];
%         STATS(i).MNIcoords = [];
%     end
    STATS(i).elec = DATA(SortCat(i)).elec;
    STATS(i).basemeanFR = DATA(SortCat(i)).(alignment).basemeanIFR;
    
    excit = find(DATA(SortCat(i)).(alignment).sig_excit==1);
    if ~isempty(excit)
        STATS(i).Excit = excit(1)/spkSampRate + ...
            DATA(SortCat(i)).(alignment).respInterval(1);
    else
        STATS(i).Excit = [];%'-';
    end
    
    inhib = find(DATA(SortCat(i)).(alignment).sig_inhib==1);
    if ~isempty(inhib)
        STATS(i).Inhib = inhib(1)/spkSampRate + ...
            DATA(SortCat(i)).(alignment).respInterval(1);
    else
        STATS(i).Inhib = [];%'-';
    end
    
    if ~isempty(excit) && ~isempty(inhib)
        idx_mix(end+1) = SortCat(i);
    elseif ~isempty(excit)
        idx_excit(end+1) = SortCat(i);
    elseif ~isempty(inhib)
        idx_inhib(end+1) = SortCat(i);
    else
        idx_nr(end+1) = SortCat(i);
    end
    
end


for i=1:length(SortCat)
STATS(i).CueBurstRho = DATA(SortCat(i)).CueBurst.Spearman.rho;
STATS(i).CueBurstPval = DATA(SortCat(i)).CueBurst.Spearman.pval;

STATS(i).SpOnsetBurstRho = DATA(SortCat(i)).BurstSpOnset.Spearman.rho;
STATS(i).SpOnsetBurstPval = DATA(SortCat(i)).BurstSpOnset.Spearman.pval;

STATS(i).CuePauseRho = DATA(SortCat(i)).CuePause.Spearman.rho;
STATS(i).CuePausePval = DATA(SortCat(i)).CuePause.Spearman.pval;

STATS(i).SpOnsetPauseRho = DATA(SortCat(i)).PauseSpOnset.Spearman.rho;
STATS(i).SpOnsetPausePval = DATA(SortCat(i)).PauseSpOnset.Spearman.pval;
end



%% determine mean RT and SpDur

SubjectsCat = unique({STATS(:).SubjectID});
RT = [];
SpDur = [];
RejectedTrials = [];
for s=1:length(SubjectsCat)
    idx = find(strcmp(SubjectsCat{s}, {DATA(:).SubjectID}));
    sessions = unique([DATA(idx).session]);
    thisRT = [];
    thisSpDur = [];
    for session = sessions
        idx2 = intersect(idx, find([DATA(:).session]==session));
        RejectedTrials(s, session) = numel(DATA(idx2(1)).trials.ResponseReject.all);
        RT(s,session) = cat(1, thisRT, nanmean(DATA(idx2(1)).RT(setdiff(1:60, DATA(idx2(1)).trials.ResponseReject.all))));
        SpDur(s,session) = cat(1, thisSpDur, nanmean(DATA(idx2(1)).SpDur(setdiff(1:60, DATA(idx2(1)).trials.ResponseReject.all))));
    end
end
RT(RT==0) = NaN;
SpDur(SpDur==0) = NaN;


NormPre = round(1000*(meanRT+0.5));
NormPost = round(1000*(meanSpDur+0.5));


%% group raster for excit units
metaIFRZexcit = zeros(length(idx_excit),NormPre+NormPost);
for i = 1:length(idx_excit)  
    
        z = (DATA(idx_excit(i)).(alignment).meanifr-DATA(idx_excit(i)).(alignment).basemeanIFR)/...
            std(mean(DATA(idx_excit(i)).(alignment).IFRbase,1),[],2);
        
        pre_idx = 1:(-round(1000*(DATA(idx_excit(i)).(alignment).respInterval(1))));
        sig_pre_norm = zeros(1,NormPre);
        temp=bwconncomp(DATA(idx_excit(i)).(alignment).sig_excit(pre_idx)==1);
        for j=1:temp.NumObjects
            idx1=max([1,round(NormPre*min(temp.PixelIdxList{j})/length(pre_idx))]);
            idx2=min([NormPre,round(NormPre*max(temp.PixelIdxList{j})/length(pre_idx))]);
            sig_pre_norm(idx1:idx2) = 1;
        end
        z_pre_norm = resample(z(pre_idx), NormPre, length(pre_idx)).*sig_pre_norm;
        
        post_idx = setdiff(1:length(z), pre_idx);
        sig_post_norm = zeros(1,NormPost);
        temp=bwconncomp(DATA(idx_excit(i)).(alignment).sig_excit(post_idx)==1);
        for j=1:temp.NumObjects
            idx1=max([1,round(NormPost*min(temp.PixelIdxList{j})/length(post_idx))]);
            idx2=min([NormPost,round(NormPost*max(temp.PixelIdxList{j})/length(post_idx))]);
            sig_post_norm(idx1:idx2) = 1;
        end
        z_post_norm = resample(z(post_idx), NormPost, length(post_idx)).*sig_post_norm;
        
        sig = cat(2, z_pre_norm, z_post_norm);       
        metaIFRZexcit(i,:) = sig';

end
idx_sig = {};
for i=1:size(metaIFRZexcit,1)
    [m, idx_sig{i}] = max(metaIFRZexcit(i,:));
    if m==0
        idx_sig{i} = size(metaIFRZexcit,2);
    end
end
[~, order] = sort(cell2mat(idx_sig));
figure; 
subplot('Position', [0 0 1 0.25])
imagesc((1:size(metaIFRZexcit,2))/1000-NormPre/1000, ...
    1:size(metaIFRZexcit,1), metaIFRZexcit(order,:)) 
colormap jet;
hold on; plot([0 0], ylim.*[1 1], 'w')
%CLIM = get(gca, 'clim')
%CLIM = max(abs(CLIM))*[-1 1];
CLIM = set(gca, 'clim', [-20 33]);
cm = colormap;
cm(1:24,:) = cm(24:-1:1,:);
cm(25:64,:) = cm(64:-1:25,:);
cm(25,:)=[0 0 0];
colormap(cm)

%% group raster for inhib units
idx_inhib_stn = intersect(idx_inhib,idxSTN)
metaIFRZinhib = zeros(length(idx_inhib_stn),NormPre+NormPost);

for i = 1:length(idx_inhib_stn)  

        meanisits  = mean(DATA(idx_inhib_stn(i)).(alignment).ISITSdata,1);
        z = -(meanisits-DATA(idx_inhib_stn(i)).(alignment).basemeanISITS)/...
            std(mean(DATA(idx_inhib_stn(i)).(alignment).ISITSbase,1),[],2);
        
        pre_idx = 1:(-round(1000*(DATA(idx_inhib_stn(i)).(alignment).respInterval(1))));
        sig_pre_norm = zeros(1,NormPre);
        temp=bwconncomp(DATA(idx_inhib_stn(i)).(alignment).sig_inhib(pre_idx)==1);
        for j=1:temp.NumObjects
            idx1=max([1,round(NormPre*min(temp.PixelIdxList{j})/length(pre_idx))]);
            idx2=min([NormPre,round(NormPre*max(temp.PixelIdxList{j})/length(pre_idx))]);
            sig_pre_norm(idx1:idx2) = 1;
        end
        z_pre_norm = resample(z(pre_idx), NormPre, length(pre_idx)).*sig_pre_norm;
        
        post_idx = setdiff(1:length(z), pre_idx);
        sig_post_norm = zeros(1,NormPost);
        temp=bwconncomp(DATA(idx_inhib_stn(i)).(alignment).sig_inhib(post_idx)==1);
        for j=1:temp.NumObjects
            idx1=max([1,round(NormPost*min(temp.PixelIdxList{j})/length(post_idx))]);
            idx2=min([NormPost,round(NormPost*max(temp.PixelIdxList{j})/length(post_idx))]);
            sig_post_norm(idx1:idx2) = 1;
        end
        z_post_norm = resample(z(post_idx), NormPost, length(post_idx)).*sig_post_norm;
        
        sig = cat(2, z_pre_norm, z_post_norm);       
        metaIFRZinhib(i,:) = sig';
end
idx_sig = {};
for i=1:size(metaIFRZinhib,1)
    [m, idx_sig{i}] = min(metaIFRZinhib(i,:));
    if m==0
        idx_sig{i} = size(metaIFRZinhib,2);
    end
end
[~, order] = sort(cell2mat(idx_sig));
figure; 
subplot('Position', [0 0 1 0.25])
imagesc((1:size(metaIFRZinhib,2))/1000-NormPre/1000, ...
    1:size(metaIFRZinhib,1), metaIFRZinhib(order,:)) 
colormap jet;
hold on; plot([0 0], ylim.*[1 1], 'w')
%CLIM = get(gca, 'clim');
% CLIM = max(abs(CLIM))*[-1 1];
% CLIM = set(gca, 'clim', CLIM);
CLIM = set(gca, 'clim', [-20 33]);
%colormap(cm3)
cm = colormap;
cm(1:24,:) = cm(24:-1:1,:);
cm(25:64,:) = cm(64:-1:25,:);
cm(25,:)=[0 0 0];
colormap(cm)


%% group raster for mix units
%idx_mix_stn = intersect(idx_mix,idxSTN);
metaIFRZmix = zeros(length(idx_mix),NormPre+NormPost);
for i = 1:length(idx_mix)  
        
        % excit pre
        z = (DATA(idx_mix(i)).(alignment).meanifr-DATA(idx_mix(i)).(alignment).basemeanIFR)/...
            std(mean(DATA(idx_mix(i)).(alignment).IFRbase,1),[],2);
        
        pre_idx = 1:(-round(1000*(DATA(idx_mix(i)).(alignment).respInterval(1))));
        sig_pre_norm = zeros(1,NormPre);

        temp=bwconncomp(DATA(idx_mix(i)).(alignment).sig_excit(pre_idx)==1);
        for j=1:temp.NumObjects
            idx1=max([1,round(NormPre*min(temp.PixelIdxList{j})/length(pre_idx))]);
            idx2=min([NormPre,round(NormPre*max(temp.PixelIdxList{j})/length(pre_idx))]);
            sig_pre_norm(idx1:idx2) = 1;
        end
        z_pre_norm1 = resample(z(pre_idx), NormPre, length(pre_idx)).*sig_pre_norm;
        
        % excit post
        post_idx = setdiff(1:length(z), pre_idx);
        sig_post_norm = zeros(1,NormPost);

        temp=bwconncomp(DATA(idx_mix(i)).(alignment).sig_excit(post_idx)==1);
        for j=1:temp.NumObjects
            idx1=max([1,round(NormPost*min(temp.PixelIdxList{j})/length(post_idx))]);
            idx2=min([NormPost,round(NormPost*max(temp.PixelIdxList{j})/length(post_idx))]);
            sig_post_norm(idx1:idx2) = 1;
        end
        z_post_norm1 = resample(z(post_idx), NormPost, length(post_idx)).*sig_post_norm;
        
        % inhib pre
        sig_pre_norm = zeros(1,NormPre);
        meanisits  = mean(DATA(idx_mix(i)).(alignment).ISITSdata,1);
        z = -(meanisits-DATA(idx_mix(i)).(alignment).basemeanISITS)/...
            std(mean(DATA(idx_mix(i)).(alignment).ISITSbase,1),[],2);
        
        temp=bwconncomp(DATA(idx_mix(i)).(alignment).sig_inhib(pre_idx)==1);
        for j=1:temp.NumObjects
            idx1=max([1,round(NormPre*min(temp.PixelIdxList{j})/length(pre_idx))]);
            idx2=min([NormPre,round(NormPre*max(temp.PixelIdxList{j})/length(pre_idx))]);
            sig_pre_norm(idx1:idx2) = 1;
        end
        z_pre_norm2 = resample(z(pre_idx), NormPre, length(pre_idx)).*sig_pre_norm;
        
        %inhib post
        sig_post_norm = zeros(1,NormPost);
        meanisits  = mean(DATA(idx_mix(i)).(alignment).ISITSdata,1);
        z = -(meanisits-DATA(idx_mix(i)).(alignment).basemeanISITS)/...
            std(mean(DATA(idx_mix(i)).(alignment).ISITSbase,1),[],2);
        
        temp=bwconncomp(DATA(idx_mix(i)).(alignment).sig_inhib(post_idx)==1);
        for j=1:temp.NumObjects
            idx1=max([1,round(NormPost*min(temp.PixelIdxList{j})/length(post_idx))]);
            idx2=min([NormPost,round(NormPost*max(temp.PixelIdxList{j})/length(post_idx))]);
            sig_post_norm(idx1:idx2) = 1;
        end
        z_post_norm2 = resample(z(post_idx), NormPost, length(post_idx)).*sig_post_norm;      
        
        z_pre_norm = z_pre_norm1 + z_pre_norm2;
        z_post_norm = z_post_norm1 + z_post_norm2;

        sig = cat(2, z_pre_norm, z_post_norm);     

        
        metaIFRZmix(i,:) = sig';

end

figure; 
subplot('Position', [0 0 1 0.25])
imagesc((1:size(metaIFRZmix,2))/1000-NormPre/1000, ...
    1:size(metaIFRZmix,1), metaIFRZmix) 
colormap jet;
hold on; plot([0 0], ylim.*[1 1], 'w')
%CLIM = get(gca, 'clim');
% CLIM = max(abs(CLIM))*[-1 1];
% CLIM = set(gca, 'clim', CLIM);
CLIM = set(gca, 'clim', [-20 33]);
colormap(cm3)
cm = colormap;
cm(1:24,:) = cm(24:-1:1,:);
cm(25:64,:) = cm(64:-1:25,:);
cm(25,:)=[0 0 0];
colormap(cm)

for idx=union(idx_excit,idx_mix);


    %% Testing for locking to RT using Cue-locked data
    delta = {};
    
    Psig = bwconncomp(~isnan(DATA(idx).sig_excit));
    %figure; plot(DATA(idx).sig_excit)
    Tsig = min([length(Psig.PixelIdxList{1}), 500]);
    
    for i=1:length(DATA_Cue(idx).trange)
        for j=1:(4*spkSampRate)
            t = DATA_Cue(idx).trial_samples(i);
            delta{i}(j) = ...
                sum(DATA_Cue(idx).IFR((t:(t+Tsig-1))+j)) ...
                / sum(DATA_Cue(idx).IFR(((t-Tsig):(t-1))+j));
        end
    end
    
    SigLat = [];
    for i=1:length(delta)
        [~, SigLat(i)] = max(delta{i});
    end
    SigLat = (SigLat-1*spkSampRate)'/spkSampRate;
    
    % i=18;
    % t = DATA_Cue(idx).trial_samples(i);
    % figure; plot((1:1001)/spkSampRate, DATA_Cue(idx).IFR(t:(t+1000)))
    % hold on; plot((1:length(delta{i}))/spkSampRate, delta{i})
    % hold on; plot(SigLat(i)*[1 1], ylim)
    
    [RHO,PVAL] = corr(SigLat, DATA_Cue(idx).RT)
    
    STATS(find(SortCat==idx)).CuelockRHO = RHO;
    STATS(find(SortCat==idx)).CuelockPVAL = PVAL;
    
    [b] = regress(DATA_Cue(idx).RT, [ones(size(SigLat)), SigLat]);
    figure; plot(SigLat, DATA_Cue(idx).RT, 'bo')
    hold on; plot(xlim, b(2)*xlim + b(1), 'b')
    
    
    %% Testing for locking to Cue using ITI-locked data
    delta = {};
    
    Psig = bwconncomp(~isnan(DATA_Cue(idx).sig_excit));
    %figure; plot(DATA_Cue(idx).sig_excit)
    Tsig = min([length(Psig.PixelIdxList{1}), 500]);
    
    for i=1:length(DATA_ITI(idx).trange)
        for j=1:(3.5*spkSampRate)
            t = DATA_ITI(idx).trial_samples(i);
            delta{i}(j) = ...
                sum(DATA_ITI(idx).IFR((t:(t+Tsig-1))+j)) ...
                / sum(DATA_ITI(idx).IFR(((t-Tsig):(t-1))+j));
        end
    end
    
    SigLat = [];
    for i=1:length(delta)
        [~, SigLat(i)] = max(delta{i});
    end
    %SigLat = SigLat'/spkSampRate;
    SigLat = (SigLat-0.5*spkSampRate)'/spkSampRate;
    
    % i=5;
    % t = DATA_ITI_Cue(idx).trial_samples(i);
    % figure; plot(DATA_ITI(idx).IFR(t:(t+1000)))
    % hold on; plot(delta{i})
    % hold on; plot(SigLat(5)*[1 1], ylim)
    
    
    [RHO,PVAL] = corr(SigLat, DATA_ITI(idx).CueT)
    
    STATS(find(SortCat==idx)).RTlockRHO = RHO;
    STATS(find(SortCat==idx)).RTlockPVAL = PVAL;
    
    [b] = regress(DATA_ITI(idx).CueT, [ones(size(SigLat)), SigLat]);
    hold on; plot(SigLat, DATA_ITI(idx).CueT, 'ro')
    hold on; plot(xlim, b(2)*xlim + b(1), 'r')
    
    title(['Rec ',num2str(idx),' ',DATA(idx).SubjectID,' ',DATA(idx).session,' ',DATA(idx).session]);
    
end