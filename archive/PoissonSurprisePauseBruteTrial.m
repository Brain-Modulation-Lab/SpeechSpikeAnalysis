
%%%% combination method
figdir = '/Users/brainmodulationlab/Dropbox (Brain Modulation Lab)/Speech Project/MerData/figures_DBS2000';


for idx=1:length(DATA)
    Stim = DATA(idx).Trial;
    clear PP;
    N0=4;
    fs = 1000;
    
    for trial=1:length(Stim.DD)
        %r = Stim.basemeanIFR;
        
        %% check for long trials and truncate to SpeechEnd + 1s
        D = Stim.DD{trial};
        if length(D)>(Stim.idxEventTimes{trial}(end)+fs)
            fprintf('Trial segment shortered from %4.2f secs to %4.2f secs.\n', length(D)/fs, (Stim.idxEventTimes{trial}(end)+fs)/fs)
            D = D(1:(Stim.idxEventTimes{trial}(end)+fs));
        end
        
        if nnz(D) > N0;
            Tavg = length(D)/(fs*nnz(D));
            r = single(1/Tavg);
            indD = find(D);
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
            
            n = cellfun(@(x) single(length(x)-1), comb);
            T = cellfun(@(x) single((indD(x(end))-indD(x(1)))/fs), comb);
            
            %S = arrayfun(@(x,y) -log10(poisscdf(x,r*y,'upper')), n, T);
            
            S = zeros(Ncomb,1);
            
            fprintf('Unit %d, trial %d: Ncomb=%d...', idx, trial, Ncomb);
            parfor i=1:Ncomb
                S(i) = -log10(poisscdf(n(i),r*T(i)));
            end
            fprintf(' done.\n');
            
            [Smax, Sind] = max(S);
            Sind = Sind(1);
            ii = comb{Sind};
            P.ii = setdiff(ii, ii(end))';
            P.S = double(Smax);
            
        else
            P.ii = [];
            P.S = 0;
        end
        PP(trial) = P;
    end
    
%     figtitle = sprintf('SpeechOnset: Rec %d %s, %5.3f, %s, Unit x, (%s)', ...
%         idx, DATA(idx).SubjectID, DATA(idx).depth, DATA(idx).elec, DATA(idx).grade);
%     
%     RT = DATA(idx).RT(DATA(idx).trange);
%     [~,order] = sort(RT);
%     
%     plot_rasterSpeechTaskBurstsOrder( [figdir,'/Cue2PSmarked'], figtitle, Stim, DATA(idx).trange, PP, [], order, RT, 'b');
    
    DATA(idx).Trial.PP = PP;
end