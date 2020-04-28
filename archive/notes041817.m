Asorts = intersect(find(strcmp({UnitList(:).grade}, 'A')),...
    find(cellfun(@(x) length(x)>10, {DATA(:).trange})));

Bsorts = intersect(find(strcmp({UnitList(:).grade}, 'B')),...
    find(cellfun(@(x) length(x)>10, {DATA(:).trange})));

BsortsSU = intersect(Bsorts,...
    find(strcmp({UnitList(:).RecType}, 'SU')));

BsortsMU = intersect(Bsorts,...
    find(strcmp({UnitList(:).RecType}, 'MU')));

Csorts = intersect(find(strcmp({UnitList(:).grade}, 'C')),...
    find(cellfun(@(x) length(x)>10, {DATA(:).trange})));

CsortsSU = intersect(Csorts,...
    find(strcmp({UnitList(:).RecType}, 'SU')));

CsortsMU = intersect(Csorts,...
    find(strcmp({UnitList(:).RecType}, 'MU')));


STATS = [];

%SortCat = Asorts;
%SortCat = union(union(Asorts,Bsorts),Csorts);
SortCat = union(Asorts,BsortsSU);
ElecOrder = {{'Cent'}, {'Post'}, {'Med', 'Lat'}};

idx_excit = [];
idx_inhib = [];
idx_mix = [];
idx_nr = [];

for i=1:length(SortCat)
    
    STATS(i).SubjectID = DATA(SortCat(i)).SubjectID;
    %STATS(i).depth = UnitList(SortCat(i)).depth;
    
    eln = find(cellfun(@(x) max(strcmp(DATA(SortCat(i)).elec, x)), ElecOrder));
    
    %STATS(i).depth = TrialData(DATA(SortCat(i)).idx_trials).MicroNormDepth{eln};
    
    STATS(i).elec = DATA(SortCat(i)).elec;
    STATS(i).basemeanFR = DATA(SortCat(i)).basemeanIFR;
    
    excit = find(DATA(SortCat(i)).sig_excit==1);
    if ~isempty(excit)
        STATS(i).Excit = excit(1)/spkSampRate + respInterval(1);
    else
        STATS(i).Excit = [];%'-';
    end
    
    inhib = find(DATA(SortCat(i)).sig_inhib==1);
    if ~isempty(inhib)
        STATS(i).Inhib = inhib(1)/spkSampRate + respInterval(1);
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