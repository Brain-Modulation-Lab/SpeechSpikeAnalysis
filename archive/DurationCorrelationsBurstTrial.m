figdir = '/Users/brainmodulationlab/Dropbox (Brain Modulation Lab)/Speech Project/MerData/figures_DBS2000';

alignment = 'Trial';

CorrType = 'Spearman';

for idx=1:length(DATA)
    clear SpDur BurstDur;
    SpDur = DATA(idx).SpDur(DATA(idx).trange);
    
    BBthresh = quantile([DATA(idx).(alignment).BBbase(:).S],0.75);

    for trial=1:length(DATA(idx).trange)
        

        if ~isempty(ii) && DATA(idx).(alignment).BB(trial).S > BBthresh

            BurstDur.T(trial) = DATA(idx).(alignment).BB(trial).t(end) - DATA(idx).(alignment).BB(trial).t(1);

        else
            BurstDur.T(trial) = NaN;
            SpDur(trial) = NaN;
        end
    end
    
    % get rid of trials with no burst detected
    SpDur(isnan(SpDur)) = [];
    BurstDur.T(isnan(BurstDur.T)) = [];
    
    BurstDur.T = BurstDur.T';
    

    if length(SpDur) > 3
        
        [BurstDur.Pearson.rho,BurstDur.Pearson.pval] = corr(SpDur, BurstDur.T, 'type', 'Pearson');
        [BurstDur.Spearman.rho,BurstDur.Spearman.pval] = corr(SpDur, BurstDur.T, 'type', 'Spearman');
        [BurstDur.b] = regress(BurstDur.T, [ones(size(SpDur)), SpDur]);
        fh = figure; plot(SpDur, BurstDur.T, 'bo')
        hold on; plot(xlim, BurstDur.b(2)*xlim + BurstDur.b(1), 'b')
        
        XLIM=xlim;
        YLIM=ylim;
        hold on; text(XLIM(1)+0.1*XLIM*[-1; 1], YLIM(1)+0.9*YLIM*[-1; 1], ...
            ['SpDur vs Burst Buration: rho= ', ...
            num2str(BurstDur.(CorrType).rho, '%4.2f'), ...
            ' pval= ', num2str(BurstDur.(CorrType).pval, '%2.1e')]);
        
        figtitle = sprintf('Rec %d %s, %5.3f, %s, Unit x, (%s)', ...
            idx, DATA(idx).SubjectID, DATA(idx).depth, DATA(idx).elec, DATA(idx).grade);
        
        xlabel('Speech Duration (s)');
        ylabel('Burst Duration (s)')
        title(figtitle);
        
%         saveas(fh,[figdir,'/DurationCorrelationsTrialBaselinedBurst/',figtitle,'.pdf'], 'pdf');
        close(fh);
    else
        BurstDur.Pearson.rho=0;
        BurstDur.Pearson.pval=1;
        BurstDur.Spearman.rho=0;
        BurstDur.Spearman.pval=1;

    end
    DATA(idx).BurstDur = BurstDur;
end

alpha=0.05;

%SortCat = union(Asorts,Bsorts);
SortCat = union(union(Asorts,Bsorts),Csorts);

idxSigCorrSpDur = find(arrayfun(@(x) ...
    x.BurstDur.(CorrType).pval<alpha, DATA));


idxSigCorrSpDur = intersect(idxSigCorrSpDur, union(idx_excit_Cue, idx_mix_Cue));
idxSigCorrSpDur = intersect(idxSigCorrSpDur, SortCat);
idxSigCorrSpDur = intersect(idxSigCorrSpDur, idxSTN);

