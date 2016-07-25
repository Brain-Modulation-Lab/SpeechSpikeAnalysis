% Analyze Subject's Units
subjectInfo;
%load(recFN);
respInterval = [-2 2.5];
MerData = struct([]);
dbg= 1;

for ii = 1:length(rec_idx)
    
    % Figure out task session event timing
    EventTimes = sort([Rec(rec_idx(ii)).DigUp Rec(rec_idx(ii)).DigDown]);
    if ii==1
        [EventTimesTrellis] = GetEventData([dataFN '.nev']);
    end
    Event0 = EventTimesTrellis(find(EventTimesTrellis>cutoff_times(rec_idx(ii)),1,'first'));
    tstart= Event0 - EventTimes(1); % This is now in trellis time
    tend = tstart + length(Rec(rec_idx(ii)).Vraw.ts)/Vraw_sampRate;
    SkipEvents = Rec(rec_idx(ii)).SkipEvents;
    nTrials = 120;
    EventTimes1 = [EventTimes EventTimes(end)+2]; % Why add an extra event here?
    
    for trial=1:nTrials
        StimulusEvent1 = SkipEvents + 4*trial;
        if ~isempty(Rec(rec_idx(ii)).AudioStart{trial})
            % This is a neuroOmega time
            ResponseTimes(trial) = EventTimes1(StimulusEvent1) + Rec(rec_idx(ii)).AudioStart{trial}/sampRate;
        else 
            ResponseTimes(trial) = NaN;
        end
    end
    
    % What do have/plot the audio signal too
%     if ii ==1 && exist([dataFN '.ns5'])
%         [~, AudioFull, AnalogElectrodeIDs] = GetAnalogData([dataFN '.ns5'], sampRate, 10269, [], []);
%     elseif ii==1
%             AudioFull = NaN*zeros(round((EventTimesTrellis(end)-EventTimesTrellis(1))*sampRate),1);
%     end
    Audio = AudioFull(round(tstart*sampRate):round(tend*sampRate));
    AudioEnv = abs(hilbert(highpassfilter(double(Audio),sampRate,100)));
    AudioEnv = smooth(AudioEnv,1500); %50 ms
%    tAudio = (0:(length(Audio)-1))/sampRate;
    %tAudioFull = (0:(length(AudioFull)-1))/sampRate - tstart;
%     clear('AudioEnvTrials', 'tAudioTrials');
%     for trial=1:length(ResponseTimes)
%         if ~isnan(ResponseTimes(trial))
%             audioWind = (round(sampRate*ResponseTimes(trial))-30000*respInterval(1)):(round(sampRate*ResponseTimes(trial))+sampRate*respInterval(2));
%             AudioEnvTrials(:,trial) = AudioEnv(audioWind);
%             tAudioTrials(:,trial) = tAudio(audioWind);
%         end
%     end
    
    % Loop through the electrodes, assemble units from each 
    for jj=1:length(electrodeList)
        spkFile = dir([subjectName '_' tasks{ii} '_' electrodeList{jj}, '*']);
        if ~isempty(spkFile)
            spikeFN = spkFile(1).name;
            %% Read in the associated spike time data
            spikeTimes = readSpikeTimes(spikeFN);
            nUnits = max(unique(spikeTimes(:,1)));
            V = Rec(rec_idx(ii)).Vraw.ts; % Raw waveforms
            timeFull = (0:(size(V,2)-1))/Vraw_sampRate; %timeseries
            %figure; ah = axes();
            fprintf('%s contains %d units.\n', spikeFN, nUnits);
            for n=1:nUnits
                ts = spikeTimes(spikeTimes(:,1) == n, 2);
%                 spikeInd = round(ts*Vraw_sampRate);
%                 spikeWind = -20:80;
%                 y = zeros(length(spikeInd), length(spikeWind));
%                 for k = 1:length(spikeInd)
%                     y(k,:) = V(j, spikeWind+spikeInd(k)) - V(j, spikeInd(j));
%                     %line(spikeWind./Vraw_sampRate, y, 'Color', 'k');
%                 end
                %plotSpikeWaveform(spikeWind./Vraw_sampRate, y, ah);
                
                unit = UnitResponse(timeFull, ts, Vraw_sampRate, [], 20, EventTimes1, SkipEvents, 4);
                
                respWind = round(respInterval(1)*unit.spkSampRate):round(respInterval(2)*unit.spkSampRate);
                ifrzm = NaN*zeros(unit.nTrials, length(respWind)); 
                spkm = NaN*zeros(size(ifrzm));
                ifrm = NaN*zeros(size(ifrzm));
                zi = find(respWind == 0);
                for kk = 1:unit.nTrials
                    if ~isempty(Rec(rec_idx(ii)).AudioStart{kk})
                        audioInd = unit.trial(kk).EventInds(4) + round(Rec(rec_idx(ii)).AudioStart{kk}/sampRate*unit.spkSampRate);
                        offset = unit.trial(kk).EventInds(1);
                        trial_inds = audioInd+respWind-offset;
                        if(trial_inds(1) < 1)
                            fprintf('%f time between trial start and cue\n', (unit.trial(kk).EventInds(4)-unit.trial(kk).EventInds(1))/unit.spkSampRate);
                            fprintf('%f Audio latency\n', Rec(rec_idx(ii)).AudioStart{kk}/sampRate);
                            disp('Trial indices are negative');             
                        end
                        ifrm(kk,:) = unit.IFR(audioInd+respWind);
                        spkm(kk,:) = unit.D(audioInd+respWind); 
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
                
                %meanifrz = nanmean(ifrzm); stdifrz = nanstd(ifrzm)./sqrt(Unit.nTrials);
                meanifr = nanmean(ifrm); stdifr = nanstd(ifrm)./sqrt(unit.nTrials);
                [fh, rasterh, meanh, ~, labelh] = makeRasterIFRplot(sprintf('%s, Unit %d', spikeFN, n));
                plott = respWind./unit.spkSampRate;
                plotRastersandIFRs(plott, meanh, meanifr, stdifr, rasterh, spkm, 0, [0 0 0]);
                %plotRastersandIFRs(plott, meanh, meanifrz, stdifrz, rasterh, spkm, 0, [0 0 0]);
                % Mark up
                plot(rasterh, [0 0], [0 unit.nTrials], 'Color', [.4 .4 .4], 'LineWidth', 2); 
                plot(meanh, [0 0], meanh.YLim, 'Color', [.4 .4 .4], 'LineWidth', 2);
                rasterh.XLim = [plott(1) plott(end)];
                meanh.XLim = [plott(1) plott(end)];
                plot(meanh, meanh.XLim, [unit.zBound(1) unit.zBound(1)], 'r');
                plot(meanh, meanh.XLim, [unit.zBound(2) unit.zBound(2)], 'r');
                xlabel(meanh, 'Time relative to speech onset (sec)', 'FontSize', 16);
                ylabel(meanh, 'Average Firing Rate (Hz)', 'FontSize', 16);
                
                MerData(rec_idx(ii)).(electrodeList{jj}).Units(n) = unit;
                %getAudioResponse
            end
        end
    end
    MerData(rec_idx(ii)).EventTimes = EventTimes1;
    MerData(rec_idx(ii)).ResponseTimes = ResponseTimes;
    MerData(rec_idx(ii)).Audio = Audio;
    MerData(rec_idx(ii)).AudioEnv = AudioEnv;
end


%% 
for wl = 1:length(rec_idx)
    for jj=1:length(electrodeList)
        electrode = MerData(rec_idx(wl)).(electrodeList{jj});
        if isfield(electrode, 'Units')
            nUnits = length(MerData(rec_idx(wl)).(electrodeList{jj}).Units);
        else
            nUnits = 0;
        end
        for n=1:nUnits
            unit = MerData(rec_idx(wl)).(electrodeList{jj}).Units(n);
            plotUnitRespByStim;
        end
    end
end