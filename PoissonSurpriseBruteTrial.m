
%%%% combination method
figdir = '/Users/brainmodulationlab/Dropbox (Brain Modulation Lab)/Speech Project/MerData/figures_DBS2000';


for idx=1:length(DATA)
    Stim = DATA(idx).Trial;
    clear BB PP BBbase PPbase;
    N0=4;
    fs = 1000;
    
    DDbase = arrayfun(@(x) DATA(idx).D((x-fs+1):x), DATA(idx).SpeechOnset.base_samples(DATA(idx).trange), ...
        'uniformoutput', false);
    
    %r = fs*nnz(cell2mat(DDbase))/length(cell2mat(DDbase));
    
    for trial=1:length(Stim.DD)
        %r = Stim.basemeanIFR;
        r = fs*nnz(Stim.DD{trial})/length(Stim.DD{trial});
        
        clear B P
        
        %% check for long trials and truncate to SpeechEnd + 1s
        D = Stim.DD{trial};

        if length(D)>(Stim.idxEventTimes{trial}(end)+fs)
            fprintf('Trial segment shortered from %4.2f secs to %4.2f secs.\n', length(D)/fs, ...
                (Stim.idxEventTimes{trial}(end)+fs)/fs)
            D = D(1:(Stim.idxEventTimes{trial}(end)+fs));
        end
        
        D = D(Stim.idxEventTimes{trial}(2):end);  %% Start detecting at command stim (Cue)
        
        if nnz(DDbase{trial}) > N0;
            
            fprintf('Unit %d, trial %d baseline: %d spikes...', idx, trial, nnz(DDbase{trial}));
            [Sb, Sp, comb, ~] = PoissonBrute(DDbase{trial}, fs, r, N0);
            fprintf(' done.\n');
            
            [Smax, Sind] = max(Sb);
            Sind = Sind(1);
            ii = comb{Sind};
            %B.ii = setdiff(ii, ii(end))';
            B.ii = ii';
            B.S = double(Smax);
            
            [Smax, Sind] = max(Sp);
            Sind = Sind(1);
            ii = comb{Sind};
            %P.ii = setdiff(ii, ii(end))';
            P.ii = ii';
            P.S = double(Smax);
            
        else
            B.ii = [];
            B.S = 0;
            P.ii = [];
            P.S = 0;
        end
        BBbase(trial) = B;
        PPbase(trial) = P;
        
        if nnz(D) > N0;
            
            fprintf('Unit %d, trial %d: %d spikes...', idx, trial, nnz(D));
            [Sb, Sp, comb, t] = PoissonBrute(D, fs, r, N0);
            fprintf(' done.\n');
            
            [Smax, Sind] = max(Sb);
            Sind = Sind(1);
            ii = comb{Sind};
            %B.ii = setdiff(ii, ii(end))';
            B.ii = ii';
            B.t = t(B.ii);
            [~,tmax] = max(DATA(idx).IFR(Stim.start_samples(DATA(idx).trange(trial)) + Stim.idxEventTimes{trial}(2) +...
            (round(fs*B.t(1)):round(fs*B.t(end)))));
            B.tmax = B.t(1) + tmax/fs;
            B.S = double(Smax);
            
            [Smax, Sind] = max(Sp);
            Sind = Sind(1);
            ii = comb{Sind};
            %P.ii = setdiff(ii, ii(end))';
            P.ii = ii';
            P.t = t(P.ii);
            m = max(DATA(idx).ISITS(Stim.start_samples(DATA(idx).trange(trial)) + Stim.idxEventTimes{trial}(2) +...
                (round(fs*P.t(1)):round(fs*P.t(end)))));
            Ptmax = bwconncomp(DATA(idx).ISITS(Stim.start_samples(DATA(idx).trange(trial)) + Stim.idxEventTimes{trial}(2) +...
                (round(fs*P.t(1)):round(fs*P.t(end))))==m);
            [~,Pidx] = max(cellfun(@length, Ptmax.PixelIdxList));
            tmax = floor(mean(Ptmax.PixelIdxList{Pidx}));
            P.tmax = P.t(1) + tmax/fs;
            P.S = double(Smax);
            
        else
            B.ii = [];
            B.t = [];
            B.tmax = [];
            B.S = 0;
            P.ii = [];
            P.t = [];
            P.tmax = [];
            P.S = 0;
        end
        BB(trial) = B;
        PP(trial) = P;
    end
    
%     figtitle = sprintf('Cue2: Rec %d %s, %5.3f, %s, Unit x, (%s)', ...
%         idx, DATA(idx).SubjectID, DATA(idx).depth, DATA(idx).elec, DATA(idx).grade);
%     
%     RT = DATA(idx).RT(DATA(idx).trange);
%     [~,order] = sort(RT);
%     
%     plot_rasterSpeechTaskBurstsOrder( [figdir,'/TrialBaselinedBurstMarked'], figtitle, ...
%         DATA(idx).Cue2, DATA(idx).trange, BB, quantile([BBbase(:).S],0.75), [], order, RT, 'b');
% 
%     figtitle = sprintf('Cue2: Rec %d %s, %5.3f, %s, Unit x, (%s)', ...
%         idx, DATA(idx).SubjectID, DATA(idx).depth, DATA(idx).elec, DATA(idx).grade);
%     
%     RT = DATA(idx).RT(DATA(idx).trange);
%     [~,order] = sort(RT);
%     
%     plot_rasterSpeechTaskBurstsOrder( [figdir,'/TrialBaselinedPauseMarked'], figtitle, ...
%         DATA(idx).Cue2, DATA(idx).trange, PP, quantile([PPbase(:).S],0.75), [], order, RT, 'b');
    
    DATA(idx).Trial.BBbase = BBbase;
    DATA(idx).Trial.PPbase = PPbase;
    DATA(idx).Trial.BB = BB;
    DATA(idx).Trial.PP = PP;
end