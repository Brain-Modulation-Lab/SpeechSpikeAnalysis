
%%%% combination method
figdir = '/Users/Witek/Dropbox (Brain Modulation Lab)/Speech Project/MerData/figures_DBS2000';


for idx=1:length(DATA)
    Stim = DATA(idx).Cue2;
    clear PP;
    N0=4;
    fs = 1000;
    
    for trial=1:size(Stim.DD,1)
        %r = Stim.basemeanIFR;
        if nnz(Stim.DD(trial,:)) > N0;
            Tavg = length(Stim.DD(trial,:))/(fs*nnz(Stim.DD(trial,:)));
            r = 1/Tavg;
            indD = find(Stim.DD(trial,:));
            indDmax = length(indD);
            Ncomb = sum(1:(indDmax-N0+1));
            comb = cell([1 Ncomb]);
            k=0;
            for i=1:(indDmax-N0+1)
                for j=(i+N0-1):indDmax
                    k=k+1;
                    comb{k} = i:j;
                end
            end
            
            n = cellfun(@(x) length(x)-1, comb);
            T = cellfun(@(x) (indD(x(end))-indD(x(1)))/fs, comb);
            S = arrayfun(@(x,y) -log10(poisscdf(x,r*y)), n, T);
            [Smax, Sind] = max(S);
            Sind = Sind(1);
            ii = comb{Sind};
            P.ii = setdiff(ii, ii(end))';
            P.S = Smax;
            
        else
            P.ii = [];
            P.S = 0;
        end
        PP(trial) = P;
    end
    
    figtitle = sprintf('SpeechOnset PSmarked: Rec %d %s, %5.3f, %s, Unit x, (%s)', ...
        idx, DATA(idx).SubjectID, DATA(idx).depth, DATA(idx).elec, DATA(idx).grade);
    
    plot_rasterSpeechTaskBursts( [figdir,'/SpeechOnsetPSPausemarked'], figtitle, Stim, DATA(idx).trange, PP, [] );
    
    DATA(idx).SpeechOnset.PP = PP;
end