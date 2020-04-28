figdir = '/Users/brainmodulationlab/Dropbox (Brain Modulation Lab)/Speech Project/MerData/figures_DBS2000';

alignment = 'Trial';

CorrType = 'Spearman';

for idx=1:length(DATA)

    RT = DATA(idx).RT(DATA(idx).trange); % + DATA(idx).SpDur(DATA(idx).trange);

    PPthresh = quantile([DATA(idx).(alignment).PPbase(:).S],0.75);
    
    clear CuePause PauseSpOnset;
    for trial=1:length(DATA(idx).trange)
        
%         timeCue = DATA(idx).(alignment).timeCue{trial};
%         timeSpOnset = DATA(idx).(alignment).timeSpOnset{trial};
%         
%         indDD = find(DATA(idx).(alignment).DD{trial});
%         ii = min(DATA(idx).(alignment).PP(trial).ii);
%         if ~isempty(ii) && DATA(idx).(alignment).PP(trial).S > PPthresh
        if DATA(idx).(alignment).PP(trial).S > PPthresh
%             CuePause.T(trial) = timeCue(indDD(ii));
%             PauseSpOnset.T(trial) = -timeSpOnset(indDD(ii));
            %% Start of burst
            CuePause.T(trial) = DATA(idx).(alignment).PP(trial).t(1);
            PauseSpOnset.T(trial) = RT(trial)-DATA(idx).(alignment).PP(trial).t(1);
            %% Max of burst
%             CuePause.T(trial) = DATA(idx).(alignment).PP(trial).tmax;
%             PauseSpOnset.T(trial) = RT(trial)-DATA(idx).(alignment).PP(trial).tmax;
            %% End of burst
%             CuePause.T(trial) = DATA(idx).(alignment).PP(trial).t(end);
%             PauseSpOnset.T(trial) = RT(trial)-DATA(idx).(alignment).PP(trial).t(end);
        else
            CuePause.T(trial) = NaN;
            PauseSpOnset.T(trial) = NaN;
            RT(trial) = NaN;
        end
    end
    
    % get rid of trials with no Pause detected
    RT(isnan(RT)) = [];
    PauseSpOnset.T(isnan(PauseSpOnset.T)) = [];
    CuePause.T(isnan(CuePause.T)) = [];
    
    PauseSpOnset.T =PauseSpOnset.T';
    CuePause.T = CuePause.T';
    
    %PauseSpOnset.T = RT - CuePause.T;

    if length(RT) > 3
        [PauseSpOnset.Pearson.rho,PauseSpOnset.Pearson.pval] = corr(RT, PauseSpOnset.T, 'type', 'Pearson');
        [PauseSpOnset.Spearman.rho,PauseSpOnset.Spearman.pval] = corr(RT, PauseSpOnset.T, 'type', 'Spearman');
        
        [PauseSpOnset.b] = regress(PauseSpOnset.T, [ones(size(RT)), RT]);
%         fh = figure; plot(RT, PauseSpOnset.T, 'bo')
%         hold on; plot(xlim, PauseSpOnset.b(2)*xlim + PauseSpOnset.b(1), 'b')
        
        [CuePause.Pearson.rho,CuePause.Pearson.pval] = corr(RT, CuePause.T, 'type', 'Pearson');
        [CuePause.Spearman.rho,CuePause.Spearman.pval] = corr(RT, CuePause.T, 'type', 'Spearman');
        
        [CuePause.b] = regress(CuePause.T, [ones(size(RT)), RT]);
%         hold on; plot(RT, CuePause.T, 'ro')
%         hold on; plot(xlim, CuePause.b(2)*xlim + CuePause.b(1), 'r')
%         XLIM=xlim;
%         YLIM=ylim;
%         hold on; text(XLIM(1)+0.1*XLIM*[-1; 1], YLIM(1)+0.9*YLIM*[-1; 1], ...
%             ['RT vs Pause-to-SpeechOnset: rho= ', ...
%             num2str(PauseSpOnset.(CorrType).rho, '%4.2f'), ...
%             ' pval= ', num2str(PauseSpOnset.(CorrType).pval, '%2.1e'), ...
%             sprintf('\n'), ...
%             'RT vs Cue-to-Pause: rho= ', ...
%             num2str(CuePause.(CorrType).rho, '%4.2f'), ...
%             ' pval= ', num2str(CuePause.(CorrType).pval, '%2.1e')]);
%         
%         figtitle = sprintf(' Rec %d %s, %5.3f, %s, Unit x, (%s)', ...
%             idx, DATA(idx).SubjectID, DATA(idx).depth, DATA(idx).elec, DATA(idx).grade);
%         
%         xlabel('Reaction Time (s)');
%         ylabel('Interval duration (s)')
%         title(figtitle);
%         
%         saveas(fh,[figdir,'/TimingCorrelationsTrialBaselinedPauseStart/',figtitle,'.pdf'], 'pdf');
%         close(fh);
    else
        CuePause.(CorrType).rho=0;
        CuePause.(CorrType).pval=1;
        PauseSpOnset.(CorrType).rho = 0;
        PauseSpOnset.(CorrType).pval = 1;
    end
%     DATA(idx).CuePause = CuePause;
%     DATA(idx).PauseSpOnset = PauseSpOnset;
end

alpha=0.05;
%SortCat = Csorts;
%SortCat = union(Asorts,Bsorts);
SortCat = union(union(Asorts,Bsorts),Csorts);

SortCat = intersect(SortCat, idxSTN);

idxSigCorrSpOnset = find(arrayfun(@(x) ...
    x.PauseSpOnset.(CorrType).pval<=alpha&&x.CuePause.(CorrType).pval>alpha, DATA));
idxSigCorrCue = find(arrayfun(@(x) ...
    x.CuePause.(CorrType).pval<=alpha&&x.PauseSpOnset.(CorrType).pval>alpha, DATA));
idxSigCorrBoth = find(arrayfun(@(x) ...
    x.CuePause.(CorrType).pval<=alpha&&x.PauseSpOnset.(CorrType).pval<=alpha, DATA));

idxSigCorrSpOnset = intersect(idxSigCorrSpOnset, union(idx_inhib_Cue, idx_mix_Cue));
idxSigCorrSpOnset = intersect(idxSigCorrSpOnset, SortCat);

idxSigCorrCue = intersect(idxSigCorrCue, union(idx_inhib_SpOnset, idx_mix_SpOnset));
idxSigCorrCue = intersect(idxSigCorrCue, SortCat);

idxSigCorrBoth = intersect(idxSigCorrBoth, union(union(idx_inhib_Cue, idx_inhib_SpOnset), union(idx_mix_Cue, idx_mix_SpOnset)));
idxSigCorrBoth = intersect(idxSigCorrBoth, SortCat);


figure; plot(arrayfun(@(x) double(x.PauseSpOnset.(CorrType).rho), DATA(idxSigCorrSpOnset)), ...
    arrayfun(@(x) double(x.CuePause.(CorrType).rho), DATA(idxSigCorrSpOnset)) ...
    , 'bo');
%     , 'b.', 'markersize', 32    );


hold on; plot(arrayfun(@(x) double(x.PauseSpOnset.(CorrType).rho), DATA(idxSigCorrCue)), ...
    arrayfun(@(x) double(x.CuePause.(CorrType).rho), DATA(idxSigCorrCue)) ...
    , 'ro');
%     , 'r.', 'markersize', 32);


hold on; plot(arrayfun(@(x) double(x.PauseSpOnset.(CorrType).rho), DATA(idxSigCorrBoth)), ...
    arrayfun(@(x) double(x.CuePause.(CorrType).rho), DATA(idxSigCorrBoth)) ...
    , 'ko');
%     , 'k.', 'markersize', 32);

hold on; plot([-1 1], [-1 1])

xlabel('Correlation Coefficient (RT vs Pause-to-SpeechOnset)')

ylabel('Correlation Coefficient (RT vs Cue-to-Pause)')

title('Pause correlations');