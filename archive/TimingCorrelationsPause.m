figdir = '/Users/Witek/Dropbox (Brain Modulation Lab)/Speech Project/MerData/figures_DBS2000';

for idx=1:length(DATA)
    
    time = linspace(DATA(idx).SpeechOnset.respInterval(1), ...
        DATA(idx).SpeechOnset.respInterval(2), ...
        size(DATA(idx).SpeechOnset.DD,2));
    
    RT = DATA(idx).RT(DATA(idx).trange);

    indDD = arrayfun(@(x) find(DATA(idx).SpeechOnset.DD(x,:)), ...
        1:size(DATA(idx).SpeechOnset.DD,1), 'uniformoutput', false);
    
    PauseSpOnset.T = cellfun( @(x,y) time(indDD{y}(min(x))), ...
        {DATA(idx).SpeechOnset.PP(:).ii}, num2cell(1:length(indDD)), ...
        'uniformoutput', false)';
    % convert timing data for trials with no burst detected to NaN
    idxNoPause = find(cellfun(@isempty, PauseSpOnset.T));  
    PauseSpOnset.T(idxNoPause) = {NaN};
    PauseSpOnset.T = cellfun(@double, PauseSpOnset.T);
    CuePause.T = RT - PauseSpOnset.T;
    
    % ...or just get rid of them. 
    RT(idxNoPause) = [];
    PauseSpOnset.T(idxNoPause) = [];
    CuePause.T(idxNoPause) = [];
    
    
    if length(RT) > 3
        [PauseSpOnset.rho,PauseSpOnset.pval] = corr(RT, PauseSpOnset.T);
        
        [PauseSpOnset.b] = regress(PauseSpOnset.T, [ones(size(RT)), RT]);
        fh = figure; plot(RT, PauseSpOnset.T, 'bo')
        hold on; plot(xlim, PauseSpOnset.b(2)*xlim + PauseSpOnset.b(1), 'b')
        
        [CuePause.rho,CuePause.pval] = corr(RT, CuePause.T);
        
        [CuePause.b] = regress(CuePause.T, [ones(size(RT)), RT]);
        hold on; plot(RT, CuePause.T, 'ro')
        hold on; plot(xlim, CuePause.b(2)*xlim + CuePause.b(1), 'r')
        XLIM=xlim;
        YLIM=ylim;
        hold on; text(XLIM(1)+0.1*XLIM*[-1; 1], YLIM(1)+0.9*YLIM*[-1; 1], ...
            ['RT vs Pause-to-SpeechOnset: rho= ', ...
            num2str(PauseSpOnset.rho, '%4.2f'), ...
            ' pval= ', num2str(PauseSpOnset.pval, '%2.1e'), ...
            sprintf('\n'), ...
            'RT vs Cue-to-Pause: rho= ', ...
            num2str(CuePause.rho, '%4.2f'), ...
            ' pval= ', num2str(CuePause.pval, '%2.1e')]);
        
        figtitle = sprintf('SpeechOnset PSmarked: Rec %d %s, %5.3f, %s, Unit x, (%s)', ...
            idx, DATA(idx).SubjectID, DATA(idx).depth, DATA(idx).elec, DATA(idx).grade);
        
        xlabel('Reaction Time (s)');
        ylabel('Interval duration (s)')
        title(figtitle);
        
        saveas(fh,[figdir,'/TimingCorrelationsPause/',figtitle,'.pdf'], 'pdf');
        close(fh);
    else
        CuePause.rho=0;
        CuePause.pval=1;
        PauseSpOnset.rho = 0;
        PauseSpOnset.pval = 1;
    end
    DATA(idx).CuePause = CuePause;
    DATA(idx).PauseSpOnset = PauseSpOnset;
end

idxSigCorrSpOnset = find(arrayfun(@(x) ...
    x.PauseSpOnset.pval<=0.05&&x.CuePause.pval>0.05, DATA));
idxSigCorrCue = find(arrayfun(@(x) ...
    x.CuePause.pval<=0.05&&x.PauseSpOnset.pval>0.05, DATA));
idxSigCorrBoth = find(arrayfun(@(x) ...
    x.CuePause.pval<=0.05&&x.PauseSpOnset.pval<=0.05, DATA));

figure; plot(arrayfun(@(x) double(x.PauseSpOnset.rho), DATA(idxSigCorrSpOnset)), ...
    arrayfun(@(x) double(x.CuePause.rho), DATA(idxSigCorrSpOnset)), 'bo');

hold on; plot(arrayfun(@(x) double(x.PauseSpOnset.rho), DATA(idxSigCorrCue)), ...
    arrayfun(@(x) double(x.CuePause.rho), DATA(idxSigCorrCue)), 'ro');

hold on; plot(arrayfun(@(x) double(x.PauseSpOnset.rho), DATA(idxSigCorrBoth)), ...
    arrayfun(@(x) double(x.CuePause.rho), DATA(idxSigCorrBoth)), 'ko');
hold on; plot([-1 1], [-1 1])

xlabel('Correlation Coefficient (RT vs Pause-to-SpeechOnset)')

ylabel('Correlation Coefficient (RT vs Cue-to-Pause)')