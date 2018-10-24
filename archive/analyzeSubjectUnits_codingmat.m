coding_dir = '/Users/Witek/Dropbox/Richardson Lab/Speech Project/Decoding/coding_files';
spkDir = 'spkData';

% Analyze Subject's Units
eval([SubjectID,'subjectInfo;']);
%load(recFN);
respInterval = [-2 2.5];
MerData = struct([]);
dbg= 1;

for ii = 1:length(rec_idx)
    
    fprintf('Rec %d.\n', rec_idx(ii));
    
    %% Load EventTimes from coding files
    load([coding_dir,'/',SubjectID,'_Session',num2str(taskNum(ii))]);

    nTrials = floor((length(EventTimes)-SkipEvents)/4);
    EventTimes1 = [EventTimes EventTimes(end)+2]; % add extra event time in case session ended in the middle of trial
    
    %% Load ResponseTimes from SpeechTrials 
    AudioStart = SpeechTrials( ...
        find(strcmp([subjectName,'_Session',num2str(taskNum(ii)),'.mat'], ...
        SpeechTrials(:,1))),12);
    
    ResponseTimes=[];
    
    for trial=1:nTrials
        StimulusEvent1 = SkipEvents + 4*trial;
        if trial<=length(AudioStart) && ~isempty(AudioStart{trial})
            ResponseTimes(trial) = EventTimes1(StimulusEvent1) + AudioStart{trial};
        else 
            ResponseTimes(trial) = NaN;
        end
    end
    
    %% Loop through the electrodes, assemble units from each
    for jj=1:length(electrodeList)
        fprintf('   %s.\n', electrodeList{jj});
        spkFile = dir([spkDir,'/',subjectName,'/', subjectName, '_', tasks{ii}, '_', electrodeList{jj}, '*']);
        if ~isempty(spkFile) && ~isempty(unitList.(electrodeList{jj}){ii}.units)
            spikeFN = [spkDir,'/',subjectName,'/', spkFile(1).name];
            %Read in the associated spike time data
            spikeTimes = readSpikeTimes(spikeFN);
            %nUnits = max(unique(spikeTimes(:,1)));
            nUnits = length(unitList.(electrodeList{jj}){ii}.units);
            if isfield(Rec(rec_idx(ii)),'Raw')
                V = Rec(rec_idx(ii)).Raw; % Raw waveforms
                timeFull = (0:(size(V.ts,2)-1))/Vraw_sampRate; %timeseries
            elseif isfield(Rec(rec_idx(ii)),'Vraw')
                if isfloat(Rec(rec_idx(ii)).Vraw)
                    V = Rec(rec_idx(ii)).Vraw; % Raw waveforms
                    timeFull = (0:(size(V,2)-1))/Vraw_sampRate; %timeseries
                else
                    V = Rec(rec_idx(ii)).Vraw; % Raw waveforms
                    timeFull = (0:(size(V.ts,2)-1))/Vraw_sampRate; %timeseries
                end
            else
                disp('NO Raw or Vraw in Rec!!!');
            end
            
            %timeFull = (0:(size(V,2)-1))/Vraw_sampRate; %timeseries
            if length(timeFull)<2
                disp('timeFull is short');
            end
            %figure; ah = axes();
            fprintf('%s contains %d units.\n', spikeFN, nUnits);
            for n=1:nUnits

                tags = unitList.(electrodeList{jj}){ii}.units{n};
                ts = [];
                for tt=1:length(tags)    
                    ts = cat(1, ts, spikeTimes(spikeTimes(:,1) == tags(tt), 2));
                end
                ts = sort(ts);
                trialRange = unitList.(electrodeList{jj}){ii}.trials{n};
                if isempty(trialRange)
                    trialRange = [1 60];
                elseif trialRange(1) < 60 && trialRange(2) > 60
                    trialRange(2) = 60;
                elseif trialRange(1) >= 60
                    trialRange = [];
                end
                
                if ~isempty(trialRange)
                    unit = UnitResponse_noempty( ...
                        timeFull,  ...      % time vector from NeurOmega
                        ts,  ...            % spike times (secs)
                        Vraw_sampRate,  ... % NeurOmega spike samp rate (44000Hz)
                        spkSampRate,  ...   % spike IFR samp rate (500Hz)
                        filtSD,  ...        % gaussian width for IFR convolution, in samples (20)
                        EventTimes1,  ...
                        SkipEvents,  ...
                        4,  ...             % nEventsPerTrial
                        trialRange);
                    
                    if ~isempty(unit)
                        if unit.nTrials < 5
                            unit = [];
                        end
                    end
                else
                    unit = [];
                end
                                
                if ~isempty(unit)
                    
                    nRec = nRec+1;
                    
                    AdjustRange{ii} = [];
                    
                    respWind = round(respInterval(1)*unit.spkSampRate):round(respInterval(2)*unit.spkSampRate);
                    ifrzm = [];%NaN*zeros(unit.nTrials, length(respWind));
                    spkm = [];%NaN*zeros(size(ifrzm));
                    ifrm = [];%NaN*zeros(size(ifrzm));
                    zi = find(respWind == 0);
                                   
                    k_ind=0;
                    for kk = 1:unit.nTrials
                        %for kk = 1:unit.nTrials
                        tnum = unit.trial(kk).stim;
                        if ~isempty(AudioStart{tnum})               %%%% !!!
                            
                            k_ind = k_ind+1;
                            
                            audioInd = unit.trial(kk).EventInds(4) + round(AudioStart{tnum}*unit.spkSampRate);
                            offset = unit.trial(kk).EventInds(1);
                            trial_inds = audioInd+respWind-offset;
                            if(trial_inds(1) < 1)
                                fprintf('%f time between trial start and cue\n', (unit.trial(kk).EventInds(4)-unit.trial(kk).EventInds(1))/unit.spkSampRate);
                                fprintf('%f Audio latency\n', AudioStart{tnum});
                                disp('Trial indices are negative');
                            end
                            ifrm(k_ind,:) = unit.IFR(audioInd+respWind);
                            spkm(k_ind,:) = unit.D(audioInd+respWind);
                            % if indexing by trial, need to protect for overrun from the recorded data
                            %zInterval = [audioInd-offsetUnit.
                            %if trial_inds(1) < 1 %if any indices are below 1
                            %    trial_inds = trial_inds + trial_inds(1) + 1;
                            %end
                            %if trial_inds(end) > length(Unit.trial(jj).IFRZ)
                            %    zinds = zinds(zinds <= length(Unit.trial(jj).IFRZ));
                            %end
                            
                            
                            %ifrzm(jj,zinds) = Unit.trial(jj).IFRZ(zinds);
                        end
                    end
                    % Compute means, stds, z-scores, etc
                    basemean = unit.zMid; basestd = unit.zStd;
                    ifrz = (ifrm-basemean)/basestd;
                    meanifrz = nanmean(ifrz); stdifrz = nanstd(ifrz)./sqrt(unit.nTrials);
                    meanifr = nanmean(ifrm); stdifr = nanstd(ifrm)./sqrt(unit.nTrials);
                    
                    fname=sprintf('Rec %d %s, %5.3f, %s, Unit %d', ...
                        nRec, subjectName, Rec(rec_idx(ii)).Depth, electrodeList{jj}, n);
                    [fh, rasterh, meanh, ~, labelh] = makeRasterIFRplot(fname);
                    plott = respWind./unit.spkSampRate;
                    plotRastersandIFRs(plott, meanh, meanifr, stdifr, rasterh, spkm, 0, [0 0 0]);
                    %plotRastersandIFRs(plott, meanh, meanifrz, stdifrz, rasterh, spkm, 0, [0 0 0]);
                    % Test if the response is significantly different from baseline
                    [h, p] = permutationTestFromBaseline(ifrz);
                    sig = h*(meanh.YLim(2)-1); sig(sig==0)=NaN;
                    plot(meanh, plott, sig, 'b', 'LineWidth', 3);
                    
                    % Mark up
                    plot(rasterh, [0 0], [0 unit.nTrials], 'Color', [.4 .4 .4], 'LineWidth', 2);
                    plot(meanh, [0 0], meanh.YLim, 'Color', [.4 .4 .4], 'LineWidth', 2);
                    rasterh.XLim = [plott(1) plott(end)];
                    meanh.XLim = [plott(1) plott(end)];
                    plot(meanh, meanh.XLim, [unit.zBound(1) unit.zBound(1)], 'r');
                    plot(meanh, meanh.XLim, [unit.zBound(2) unit.zBound(2)], 'r');
                    %plot(meanh, meanh.XLim, [-2 -2], 'r');
                    %plot(meanh, meanh.XLim, [2 2], 'r');
                    xlabel(meanh, 'Time relative to speech onset (sec)', 'FontSize', 16);
                    ylabel(meanh, 'Average Firing Rate (Hz)', 'FontSize', 16);
                    
                    saveas(fh,['figures/',fname,'.pdf'], 'pdf');
                    close(fh);
      
                    %getAudioResponse
                    
                else
                    disp('   NO SPIKE TRIALS');
                    unit.trialRange = [];
                    unit.D = [];
                    unit.IFR = [];
                    unit.tSpk = [];
                    unit.spkSampRate = [];
                    unit.trial = [];
                    unit.EventInds = [];
                    unit.nTrials = [];
                    unit.zBound = [];
                    unit.zMid = [];
                    unit.zStd = [];
                end
                
                unit.sig = sig;
                unit.ifrz = ifrz;
                unit.spkm = spkm;
                unit.ifrm = ifrm;
                
                MerData(rec_idx(ii)).(electrodeList{jj}).Units(n) = unit;
            end
        end
    end
    MerData(rec_idx(ii)).EventTimes = EventTimes1;
    MerData(rec_idx(ii)).ResponseTimes = ResponseTimes;
    MerData(rec_idx(ii)).Audio = [];%Audio;
    MerData(rec_idx(ii)).AudioEnv = [];%AudioEnv;
    %%
end

save(['MerData_60/',subjectName 'MerData'], 'MerData', '-v7.3');
