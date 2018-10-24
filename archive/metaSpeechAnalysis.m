respInterval = [-1 2];

spkSampRate = 500;
filtSD = 20;
minTsig = 125; % minimum significant FR change (samples)

alphabet = 'abcdefghijklmopqrstuvwxyz';

DATA = [];

for idx = 1:length(UnitList)
    
    if ~strcmp(UnitList(idx).grade, 'F')
        
        SubjectID = UnitList(idx).subject;
        session = UnitList(idx).session;
        elec = UnitList(idx).elec;
        
        fprintf('Rec %d: %s Session%d %s...\n', idx, SubjectID, session, elec);
        
        %     load(['Rec/',SubjectID,'Rec']);
        
        eval([SubjectID,'subjectInfo;']);
        
        idx_trials = intersect(find(strcmp(SubjectID,{TrialData(:).subject})), ...
            find([TrialData(:).session]==session));
        trials = TrialData(idx_trials).trials;
        
        % Find the data corresponding to the current unit
        SpikeDataSession  = cellfun(@(x) cellfun(@str2num, regexp(x, '[1-9]', 'match')), {SpikeData(:).session});
        if max(~cellfun(@isempty, regexp(tasks, 'CS')))
            wordlist = SessionTable( ...
                find(cellfun(@(x,y) strcmp(SubjectID, x) && session==y, ...
                {SessionTable(:).subject}, {SessionTable(:).session}))).wordlist;
            idx_spike = intersect(intersect(find(strcmp(SubjectID,{SpikeData(:).subject})), ...
                find(SpikeDataSession==wordlist)),find(strcmp(elec,{SpikeData(:).elec})));
        elseif max(~cellfun(@isempty, regexp(tasks, 'Session')))
            idx_spike = intersect(intersect(find(strcmp(SubjectID,{SpikeData(:).subject})), ...
                find(SpikeDataSession==session)),find(strcmp(elec,{SpikeData(:).elec})));
        end
        
        if length(idx_spike)>1
            [s,isort] = max([SpikeData(idx_spike).sortv]);
            idx_spike = idx_spike(isort);
            fprintf('   Multiple sort versions.  Using latest sort (v%d).\n', s);
        end
        
        % Get spike timestamps
        ts = [];
        for tag = UnitList(idx).tags
            ts = cat(1, ts, SpikeData(idx_spike).PLXdata.(['sig001',alphabet(tag)]).ts);
        end
        ts = sort(ts);
        vSampRate = SpikeData(idx_spike).PLXdata.cont_fs;
        
        % Generate IFR
        t = (1:round(vSampRate*SpikeData(idx_spike).PLXdata.Duration))/vSampRate; %timeseries
        down = round(vSampRate/spkSampRate);
        t = downsample(t, down);
        D = zeros(length(t), 1);
        spkinds = round(ts.*spkSampRate);
        spkinds(spkinds < 1) = 1;
        D(spkinds) = 1; %binary spike vector
        %sets up the gaussian filter for spike density estimation
        %standard  deviation of the gaussian smoothing
        stdg = filtSD; % this standard deviation is in samples
        filtx = -4*stdg:1:4*stdg;
        filty = normpdf(filtx,0,stdg);
        filty = filty/sum(filty);
        %20ms gaussian smoothing
        IFR = conv2(D(:), filty(:), 'same');
        %normalize to the number of spikes in the trial as in Gabbiani et al. 1999
        IFR = IFR/(sum(IFR)/spkSampRate);
        IFR = IFR .* sum(D);
        
        ISITS = ISIts(D,spkSampRate);
        
        trange = 1:60;
        trange = setdiff(trange, find(isnan(trials.SpOnset)));
        trange = setdiff(trange, union(trials.BaseRejectNoise, trials.BaseRejectSpk));
        
        trange = setdiff(trange, find((trials.BaseBack-1)<UnitList(idx).tstart));
        if strcmp(UnitList(idx).tend,'end')
            tend = spkSampRate*length(D);
        else
            tend = UnitList(idx).tend;
        end
        trange = setdiff(trange, find((trials.SpOnset+respInterval(2))>tend));
        
        if length(trange)>10
            
            base_samples = round(spkSampRate*trials.BaseBack);
            IFRbase = cell2mat(arrayfun(@(x) IFR((base_samples(x)-spkSampRate):base_samples(x)), ...
                trange, 'uniformoutput', false))';
            ISITSbase = cell2mat(arrayfun(@(x) ISITS((base_samples(x)-spkSampRate):base_samples(x)), ...
                trange, 'uniformoutput', false))';
            
            trial_samples = round(spkSampRate*trials.SpOnset);
            DD = cell2mat(arrayfun(@(x) ...
                D((trial_samples(x)+spkSampRate*respInterval(1)): ...
                (trial_samples(x)+spkSampRate*respInterval(2))), ...
                trange, 'uniformoutput', false))';
            IFRdata = cell2mat(arrayfun(@(x) ...
                IFR((trial_samples(x)+spkSampRate*respInterval(1)): ...
                (trial_samples(x)+spkSampRate*respInterval(2))), ...
                trange, 'uniformoutput', false))';
            ISITSdata = cell2mat(arrayfun(@(x) ...
                ISITS((trial_samples(x)+spkSampRate*respInterval(1)): ...
                (trial_samples(x)+spkSampRate*respInterval(2))), ...
                trange, 'uniformoutput', false))';
            
            basemeanIFR = mean(IFRbase(:));
            baseIFRCI = basemeanIFR + ...
                [ icdf('Normal',0.05/38,0,1)*std(mean(IFRbase,1),[],2) ...
                -icdf('Normal',0.05/38,0,1)*std(mean(IFRbase,1),[],2) ];
            
            basemeanISITS = mean(ISITSbase(:));
            baseISITSCI = basemeanISITS + ...
                [ icdf('Normal',0.05/38,0,1)*std(mean(ISITSbase,1),[],2) ...
                -icdf('Normal',0.05/38,0,1)*std(mean(ISITSbase,1),[],2) ];
            
            meanifr = mean(IFRdata,1);
            stdifr = std(IFRdata,[],1)./sqrt(length(trange));
            
            sig_excit = round(meanifr>baseIFRCI(2));
            Psig = bwconncomp(sig_excit);
            for i = 1:Psig.NumObjects
                if length(Psig.PixelIdxList{i})<minTsig
                    sig_excit(Psig.PixelIdxList{i})=0;
                end
            end
            sig_excit(sig_excit==0)=NaN;
            
            sig_inhib = round(mean(ISITSdata,1)>baseISITSCI(2));
            Psig = bwconncomp(sig_inhib);
            for i = 1:Psig.NumObjects
                if length(Psig.PixelIdxList{i})<minTsig
                    sig_inhib(Psig.PixelIdxList{i})=0;
                end
            end
            sig_inhib(sig_inhib==0)=NaN;
            
            DATA(idx).SubjectID = SubjectID;
            DATA(idx).session = session;
            DATA(idx).elec = elec;
            DATA(idx).idx_trials = idx_trials;
            DATA(idx).idx_spike = idx_spike;
            DATA(idx).DD = DD;
            DATA(idx).IFR = IFR;
            DATA(idx).ISITS = ISITS;
            DATA(idx).trange = trange;
            DATA(idx).base_samples = base_samples;
            DATA(idx).trial_samples = trial_samples;
            DATA(idx).IFRbase = IFRbase;
            DATA(idx).ISITSbase = ISITSbase;
            DATA(idx).IFRdata = IFRdata;
            DATA(idx).ISITSdata = ISITSdata;
            DATA(idx).basemeanIFR = basemeanIFR;
            DATA(idx).baseIFRCI = baseIFRCI;
            DATA(idx).baseISITS = basemeanISITS;
            DATA(idx).baseISITSCI = baseISITSCI;
            DATA(idx).sig_excit = sig_excit;
            DATA(idx).sig_inhib = sig_inhib;
            
%             figtitle = sprintf('Rec %d %s, %5.3f, %s, Unit %d', ...
%                 idx, SubjectID, UnitList(idx).depth, elec, UnitList(idx).unit);
%             
%             respWind = round(respInterval(1)*spkSampRate):round(respInterval(2)*spkSampRate);
%             plott = respWind./spkSampRate;
%             [fh, rasterh, meanh, ~, labelh] = makeRasterIFRplot(figtitle);
%             plotRastersandIFRs(plott, meanh, meanifr, stdifr, rasterh, DD, 0, [0 0 0]);
%             plot(rasterh, [0 0], [0 length(trange)], 'Color', [.4 .4 .4], 'LineWidth', 2);
%             plot(meanh, [0 0], meanh.YLim, 'Color', [.4 .4 .4], 'LineWidth', 2);
%             rasterh.XLim = [plott(1) plott(end)];
%             meanh.XLim = [plott(1) plott(end)];
%             plot(meanh, meanh.XLim, baseIFRCI(1)*[1 1], 'r');
%             plot(meanh, meanh.XLim, baseIFRCI(2)*[1 1], 'r');
%             plot(meanh, plott, meanh.YLim(2)*sig_excit, 'r', 'LineWidth', 3);
%             plot(meanh, plott, meanh.YLim(2)*sig_inhib, 'b', 'LineWidth', 3);
%             
%             saveas(fh,['figures/',figtitle,'.pdf'], 'pdf');
%             close(fh);
        end
        
    end
    
end