%% Loop through the electrodes, assemble units from each 
    for jj=1:length(electrodeList)
        spkFile = dir([subjectName '_' tasks{ii} '_' electrodeList{jj}, '*']);
        if ~isempty(spkFile) && ~isempty(unitList.(electrodeList{jj}){ii}.units)
            spikeFN = spkFile(1).name;
            %Read in the associated spike time data
            spikeTimes = readSpikeTimes(spikeFN);
            %nUnits = max(unique(spikeTimes(:,1)));
            nUnits = length(unitList.(electrodeList{jj}){ii}.units);
            V = Rec(rec_idx(ii)).Vraw; % Raw waveforms
            timeFull = (0:(size(V,2)-1))/Vraw_sampRate; %timeseries
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
                if isempty(trialRange); trialRange = [1 120]; end

%                 spikeInd = round(ts*Vraw_sampRate);
%                 spikeWind = -20:80;
%                 y = zeros(length(spikeInd), length(spikeWind));
%                 for k = 1:length(spikeInd)
%                     y(k,:) = V(j, spikeWind+spikeInd(k)) - V(j, spikeInd(j));
%                     %line(spikeWind./Vraw_sampRate, y, 'Color', 'k');
%                 end
                %plotSpikeWaveform(spikeWind./Vraw_sampRate, y, ah);
                
                unit = UnitResponse(timeFull, ts, Vraw_sampRate, spkSampRate, filtSD, EventTimes1, SkipEvents, 4, trialRange);
                
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
                % Compute means, stds, z-scores, etc
                basemean = unit.zMid; basestd = unit.zStd;
                ifrz = (ifrm-basemean)/basestd;
                meanifrz = nanmean(ifrz); stdifrz = nanstd(ifrz)./sqrt(unit.nTrials);
                meanifr = nanmean(ifrm); stdifr = nanstd(ifrm)./sqrt(unit.nTrials);
                
                [fh, rasterh, meanh, ~, labelh] = makeRasterIFRplot(sprintf('%s, Unit %d', spikeFN, n));
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
                
                MerData(rec_idx(ii)).(electrodeList{jj}).Units(n) = unit;
                %getAudioResponse
            end
        end
    end
    MerData(rec_idx(ii)).EventTimes = EventTimes1;
    MerData(rec_idx(ii)).ResponseTimes = ResponseTimes;
    MerData(rec_idx(ii)).Audio = Audio;
    MerData(rec_idx(ii)).AudioEnv = AudioEnv;
    %%