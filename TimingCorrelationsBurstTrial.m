figdir = '/Users/brainmodulationlab/Dropbox (Brain Modulation Lab)/Speech Project/MerData/figures_DBS2000';

alignment = 'Trial';

CorrType = 'Spearman';

for idx=1:length(DATA)

    RT = DATA(idx).RT(DATA(idx).trange); % + DATA(idx).SpDur(DATA(idx).trange);
    
    BBthresh = quantile([DATA(idx).(alignment).BBbase(:).S],0.75);

    clear CueBurst BurstSpOnset;
    for trial=1:length(DATA(idx).trange)
        
%         timeCue = DATA(idx).(alignment).timeCue{trial} + DATA(idx).(alignment).idxEventTimes{trial}(2)*fs;
%         timeSpOnset = DATA(idx).(alignment).timeSpOnset{trial} + DATA(idx).(alignment).idxEventTimes{trial}(2)*fs;
        
%         indDD = find(DATA(idx).(alignment).DD{trial});
%         ii = min(DATA(idx).(alignment).BB(trial).ii);
%         if ~isempty(ii) && DATA(idx).(alignment).BB(trial).S > BBthresh
        if DATA(idx).(alignment).BB(trial).S > BBthresh
%             CueBurst.T(trial) = timeCue(indDD(ii));
%             BurstSpOnset.T(trial) = -timeSpOnset(indDD(ii));
            %% Start of burst
            CueBurst.T(trial) = DATA(idx).(alignment).BB(trial).t(1);
            BurstSpOnset.T(trial) = RT(trial)-DATA(idx).(alignment).BB(trial).t(1);
%             %% Max of burst
%             CueBurst.T(trial) = DATA(idx).(alignment).BB(trial).tmax;
%             BurstSpOnset.T(trial) = RT(trial)-DATA(idx).(alignment).BB(trial).tmax;
%             %% End of burst
%             CueBurst.T(trial) = DATA(idx).(alignment).BB(trial).t(end);
%             BurstSpOnset.T(trial) = RT(trial)-DATA(idx).(alignment).BB(trial).t(end);
        else
            CueBurst.T(trial) = NaN;
            BurstSpOnset.T(trial) = NaN;
            RT(trial) = NaN;
        end
    end
    
    % get rid of trials with no burst detected
    RT(isnan(RT)) = [];
    BurstSpOnset.T(isnan(BurstSpOnset.T)) = [];
    CueBurst.T(isnan(CueBurst.T)) = [];
    
    BurstSpOnset.T =BurstSpOnset.T';
    CueBurst.T = CueBurst.T';
    
    %BurstSpOnset.T = RT - CueBurst.T;

    if length(RT) > 3
        [BurstSpOnset.Pearson.rho,BurstSpOnset.Pearson.pval] = corr(RT, BurstSpOnset.T, 'type', 'Pearson');
        [BurstSpOnset.Spearman.rho,BurstSpOnset.Spearman.pval] = corr(RT, BurstSpOnset.T, 'type', 'Spearman');
        
        [BurstSpOnset.b] = regress(BurstSpOnset.T, [ones(size(RT)), RT]);
%         fh = figure; plot(RT, BurstSpOnset.T, 'bo')
%         hold on; plot(xlim, BurstSpOnset.b(2)*xlim + BurstSpOnset.b(1), 'b')
        
        [CueBurst.Pearson.rho,CueBurst.Pearson.pval] = corr(RT, CueBurst.T, 'type', 'Pearson');
        [CueBurst.Spearman.rho,CueBurst.Spearman.pval] = corr(RT, CueBurst.T, 'type', 'Spearman');
        
        [CueBurst.b] = regress(CueBurst.T, [ones(size(RT)), RT]);
%         hold on; plot(RT, CueBurst.T, 'ro')
%         hold on; plot(xlim, CueBurst.b(2)*xlim + CueBurst.b(1), 'r')
%         XLIM=xlim;
%         YLIM=ylim;
%         hold on; text(XLIM(1)+0.1*XLIM*[-1; 1], YLIM(1)+0.9*YLIM*[-1; 1], ...
%             ['RT vs Burst-to-SpeechOnset: rho= ', ...
%             num2str(BurstSpOnset.(CorrType).rho, '%4.2f'), ...
%             ' pval= ', num2str(BurstSpOnset.(CorrType).pval, '%2.1e'), ...
%             sprintf('\n'), ...
%             'RT vs Cue-to-Burst: rho= ', ...
%             num2str(CueBurst.(CorrType).rho, '%4.2f'), ...
%             ' pval= ', num2str(CueBurst.(CorrType).pval, '%2.1e')]);
%         
%         figtitle = sprintf('Rec %d %s, %5.3f, %s, Unit x, (%s)', ...
%             idx, DATA(idx).SubjectID, DATA(idx).depth, DATA(idx).elec, DATA(idx).grade);
%         
%         xlabel('Reaction Time (s)');
%         ylabel('Interval duration (s)')
%         title(figtitle);
%         
%         saveas(fh,[figdir,'/TimingCorrelationsTrialBaselinedBurst/',figtitle,'.pdf'], 'pdf');
%         close(fh);
    else
        CueBurst.Pearson.rho=0;
        CueBurst.Pearson.pval=1;
        BurstSpOnset.Pearson.rho = 0;
        BurstSpOnset.Pearson.pval = 1;
        CueBurst.Spearman.rho=0;
        CueBurst.Spearman.pval=1;
        BurstSpOnset.Spearman.rho = 0;
        BurstSpOnset.Spearman.pval = 1;
    end
%     DATA(idx).CueBurst = CueBurst;
%     DATA(idx).BurstSpOnset = BurstSpOnset;
end

alpha=0.05;

% SortCat = Csorts;
% SortCat = union(Asorts,Bsorts);
SortCat = union(union(Asorts,Bsorts),Csorts);

SortCat = intersect(SortCat, idxSTN);

idxSigCorrSpOnset = find(arrayfun(@(x) ...
    x.BurstSpOnset.(CorrType).pval<=alpha&&x.CueBurst.(CorrType).pval>alpha, DATA));
idxSigCorrCue = find(arrayfun(@(x) ...
    x.CueBurst.(CorrType).pval<=alpha&&x.BurstSpOnset.(CorrType).pval>alpha, DATA));
idxSigCorrBoth = find(arrayfun(@(x) ...
    x.CueBurst.(CorrType).pval<=alpha&&x.BurstSpOnset.(CorrType).pval<=alpha, DATA));

idxSigCorrSpOnset = intersect(idxSigCorrSpOnset, union(idx_excit_Cue, idx_mix_Cue));
idxSigCorrSpOnset = intersect(idxSigCorrSpOnset, SortCat);

idxSigCorrCue = intersect(idxSigCorrCue, union(idx_excit_SpOnset, idx_mix_SpOnset));
idxSigCorrCue = intersect(idxSigCorrCue, SortCat);

idxSigCorrBoth = intersect(idxSigCorrBoth, union(union(idx_excit_Cue, idx_excit_SpOnset), union(idx_mix_Cue, idx_mix_SpOnset)));
idxSigCorrBoth = intersect(idxSigCorrBoth, SortCat);

figure; plot(arrayfun(@(x) double(x.BurstSpOnset.(CorrType).rho), DATA(idxSigCorrSpOnset)), ...
    arrayfun(@(x) double(x.CueBurst.(CorrType).rho), DATA(idxSigCorrSpOnset)) ...
    , 'bo');
%     , 'b.', 'markersize', 32);
    
hold on; plot(arrayfun(@(x) double(x.BurstSpOnset.(CorrType).rho), DATA(idxSigCorrCue)), ...
    arrayfun(@(x) double(x.CueBurst.(CorrType).rho), DATA(idxSigCorrCue)) ...
    , 'ro');
%     , 'r.', 'markersize', 32);

hold on; plot(arrayfun(@(x) double(x.BurstSpOnset.(CorrType).rho), DATA(idxSigCorrBoth)), ...
    arrayfun(@(x) double(x.CueBurst.(CorrType).rho), DATA(idxSigCorrBoth)) ...
    , 'ko');
%     , 'k.', 'markersize', 32);
hold on; plot([-1 1], [-1 1])

xlabel('Correlation Coefficient (RT vs Burst-to-SpeechOnset)')

ylabel('Correlation Coefficient (RT vs Cue-to-Burst)')

title('Burst correlations');
