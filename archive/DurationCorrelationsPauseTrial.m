figdir = '/Users/brainmodulationlab/Dropbox (Brain Modulation Lab)/Speech Project/MerData/figures_DBS2000';

alignment = 'Trial';

CorrType = 'Spearman';

for idx=1:length(DATA)
    clear SpDur PauseDur;
    SpDur = DATA(idx).SpDur(DATA(idx).trange);
    
    PPthresh = quantile([DATA(idx).(alignment).PPbase(:).S],0.75);

    for trial=1:length(DATA(idx).trange)
        

        if ~isempty(ii) && DATA(idx).(alignment).PP(trial).S > PPthresh

            PauseDur.T(trial) = DATA(idx).(alignment).PP(trial).t(end) - DATA(idx).(alignment).PP(trial).t(1);

        else
            PauseDur.T(trial) = NaN;
            SpDur(trial) = NaN;
        end
    end
    
    % get rid of trials with no burst detected
    SpDur(isnan(SpDur)) = [];
    PauseDur.T(isnan(PauseDur.T)) = [];
    
    PauseDur.T = PauseDur.T';
    

    if length(SpDur) > 3
        
        [PauseDur.Pearson.rho,PauseDur.Pearson.pval] = corr(SpDur, PauseDur.T, 'type', 'Pearson');
        [PauseDur.Spearman.rho,PauseDur.Spearman.pval] = corr(SpDur, PauseDur.T, 'type', 'Spearman');
        [PauseDur.b] = regress(PauseDur.T, [ones(size(SpDur)), SpDur]);
        fh = figure; plot(SpDur, PauseDur.T, 'bo')
        hold on; plot(xlim, PauseDur.b(2)*xlim + PauseDur.b(1), 'b')
        
        XLIM=xlim;
        YLIM=ylim;
        hold on; text(XLIM(1)+0.1*XLIM*[-1; 1], YLIM(1)+0.9*YLIM*[-1; 1], ...
            ['SpDur vs Pause Buration: rho= ', ...
            num2str(PauseDur.(CorrType).rho, '%4.2f'), ...
            ' pval= ', num2str(PauseDur.(CorrType).pval, '%2.1e')]);
        
        figtitle = sprintf('Rec %d %s, %5.3f, %s, Unit x, (%s)', ...
            idx, DATA(idx).SubjectID, DATA(idx).depth, DATA(idx).elec, DATA(idx).grade);
        
        xlabel('Speech Duration (s)');
        ylabel('Pause Duration (s)')
        title(figtitle);
        
%         saveas(fh,[figdir,'/DurationCorrelationsTrialBaselinedPause/',figtitle,'.pdf'], 'pdf');
        close(fh);
    else
        PauseDur.Pearson.rho=0;
        PauseDur.Pearson.pval=1;
        PauseDur.Spearman.rho=0;
        PauseDur.Spearman.pval=1;

    end
    DATA(idx).PauseDur = PauseDur;
end

alpha=0.05;

%SortCat = union(Asorts,Bsorts);
SortCat = union(union(Asorts,Bsorts),Csorts);

idxSigCorrSpDur = find(arrayfun(@(x) ...
    x.PauseDur.(CorrType).pval<alpha, DATA));


idxSigCorrSpDur = intersect(idxSigCorrSpDur, union(idx_inhib_Cue, idx_mix_Cue));
idxSigCorrSpDur = intersect(idxSigCorrSpDur, SortCat);
idxSigCorrSpDur = intersect(idxSigCorrSpDur, idxSTN);
