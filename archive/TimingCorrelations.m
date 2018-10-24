figdir = '/Users/brainmodulationlab/Dropbox (Brain Modulation Lab)/Speech Project/MerData/figures_DBS2000';

alignment = 'Cue2';

for idx=1:length(DATA)
    
    time = linspace(DATA(idx).(alignment).respInterval(1), ...
        DATA(idx).(alignment).respInterval(2), ...
        size(DATA(idx).(alignment).DD,2));
    
    RT = DATA(idx).RT(DATA(idx).trange);

    indDD = arrayfun(@(x) find(DATA(idx).(alignment).DD(x,:)), ...
        1:size(DATA(idx).(alignment).DD,1), 'uniformoutput', false);
    
    CueBurst.T = cellfun( @(x,y) time(indDD{y}(min(x))), ...
        {DATA(idx).(alignment).BB(:).ii}, num2cell(1:length(indDD)), ...
        'uniformoutput', false)';
    % convert timing data for trials with no burst detected to NaN
    idxNoBurst = find(cellfun(@isempty, CueBurst.T));  
    CueBurst.T(idxNoBurst) = {NaN};
    CueBurst.T = cellfun(@double, CueBurst.T);
    
    BurstSpOnset.T = RT - CueBurst.T;
    
    % ...or just get rid of them. 
    RT(idxNoBurst) = [];
    BurstSpOnset.T(idxNoBurst) = [];
    CueBurst.T(idxNoBurst) = [];
    
    
    if length(RT) > 3
        [BurstSpOnset.rho,BurstSpOnset.pval] = corr(RT, BurstSpOnset.T);
        
        [BurstSpOnset.b] = regress(BurstSpOnset.T, [ones(size(RT)), RT]);
        fh = figure; plot(RT, BurstSpOnset.T, 'bo')
        hold on; plot(xlim, BurstSpOnset.b(2)*xlim + BurstSpOnset.b(1), 'b')
        
        [CueBurst.rho,CueBurst.pval] = corr(RT, CueBurst.T);
        
        [CueBurst.b] = regress(CueBurst.T, [ones(size(RT)), RT]);
        hold on; plot(RT, CueBurst.T, 'ro')
        hold on; plot(xlim, CueBurst.b(2)*xlim + CueBurst.b(1), 'r')
        XLIM=xlim;
        YLIM=ylim;
        hold on; text(XLIM(1)+0.1*XLIM*[-1; 1], YLIM(1)+0.9*YLIM*[-1; 1], ...
            ['RT vs Burst-to-SpeechOnset: rho= ', ...
            num2str(BurstSpOnset.rho, '%4.2f'), ...
            ' pval= ', num2str(BurstSpOnset.pval, '%2.1e'), ...
            sprintf('\n'), ...
            'RT vs Cue-to-Burst: rho= ', ...
            num2str(CueBurst.rho, '%4.2f'), ...
            ' pval= ', num2str(CueBurst.pval, '%2.1e')]);
        
        figtitle = sprintf('SpeechOnset PSmarked: Rec %d %s, %5.3f, %s, Unit x, (%s)', ...
            idx, DATA(idx).SubjectID, DATA(idx).depth, DATA(idx).elec, DATA(idx).grade);
        
        xlabel('Reaction Time (s)');
        ylabel('Interval duration (s)')
        title(figtitle);
        
        saveas(fh,[figdir,'/TimingCorrelationsCue2Burst/',figtitle,'.pdf'], 'pdf');
        close(fh);
    else
        CueBurst.rho=0;
        CueBurst.pval=1;
        BurstSpOnset.rho = 0;
        BurstSpOnset.pval = 1;
    end
    DATA(idx).CueBurst = CueBurst;
    DATA(idx).BurstSpOnset = BurstSpOnset;
end

alpha=0.05;

idxSigCorrSpOnset = find(arrayfun(@(x) ...
    x.BurstSpOnset.pval<=alpha&&x.CueBurst.pval>alpha, DATA));
idxSigCorrCue = find(arrayfun(@(x) ...
    x.CueBurst.pval<=alpha&&x.BurstSpOnset.pval>alpha, DATA));
idxSigCorrBoth = find(arrayfun(@(x) ...
    x.CueBurst.pval<=alpha&&x.BurstSpOnset.pval<=alpha, DATA));

figure; plot(arrayfun(@(x) double(x.BurstSpOnset.rho), DATA(idxSigCorrSpOnset)), ...
    arrayfun(@(x) double(x.CueBurst.rho), DATA(idxSigCorrSpOnset)), 'bo');

hold on; plot(arrayfun(@(x) double(x.BurstSpOnset.rho), DATA(idxSigCorrCue)), ...
    arrayfun(@(x) double(x.CueBurst.rho), DATA(idxSigCorrCue)), 'ro');

hold on; plot(arrayfun(@(x) double(x.BurstSpOnset.rho), DATA(idxSigCorrBoth)), ...
    arrayfun(@(x) double(x.CueBurst.rho), DATA(idxSigCorrBoth)), 'ko');
hold on; plot([-1 1], [-1 1])

xlabel('Correlation Coefficient (RT vs Burst-to-SpeechOnset)')

ylabel('Correlation Coefficient (RT vs Cue-to-Burst)')