
respInterval = [-2 2.5];


for ii = 1:length(DATA)
    
    
    fprintf('Rec %s, depth %5.3f, %s, Unit %d.\n', DATA{ii,1}, DATA{ii,2}, DATA{ii,3}, DATA{ii,4});
    
    unit = DATA{ii,5};
    
    AdjustRange{ii} = [];
    
    %  Adjust for trials with no spikes
    if nnz(~[unit.trial(:).hasSpikes])>0
        disp('   empty trial');
    end
    unit.nTrials = unit.nTrials-nnz(~[unit.trial(:).hasSpikes]);
    spike_trials = find([unit.trial(:).hasSpikes]);
    
    respWind = round(respInterval(1)*unit.spkSampRate):round(respInterval(2)*unit.spkSampRate);

    ifrzm = [];%NaN*zeros(unit.nTrials, length(respWind));
    spkm = [];%NaN*zeros(size(ifrzm));
    ifrm = [];%NaN*zeros(size(ifrzm));
    zi = find(respWind == 0);
    
    k_ind=0;
    for kk = 1:length(spike_trials)   %1:unit.nTrials
        
        if kk == 115
            disp('db')
        end
        
        RespInd = round(DATA{ii,6}(...
            find(DATA{ii,6} > unit.trial(spike_trials(kk)).EventInds(4)/unit.spkSampRate,1,'first')) ...
            *unit.spkSampRate);
        
        if RespInd > unit.trial(spike_trials(kk)).EventInds(5)
            fprintf('Skip trial %d, no response!\n', spike_trials(kk))
            RespInd = NaN;
        end
        
        if ~isempty(RespInd)
            if ~isnan(RespInd)&&~isnan(RespInd)
                audioInd = RespInd;
                k_ind=k_ind+1;
                ifrm(k_ind,:) = unit.IFR(audioInd+respWind);
                spkm(k_ind,:) = unit.D(audioInd+respWind);
            end
        end
        
        if nansum(unit.D(audioInd+respWind))==0
            fprintf('Empty trial %d, no spikes!\n', spike_trials(kk))
            AdjustRange{ii} = cat(1, AdjustRange{ii}, spike_trials(kk));
        end
        
    end
    
    % Compute means, stds, z-scores, etc
    basemean = unit.zMid; basestd = unit.zStd;
    ifrz = (ifrm-basemean)/basestd;
    meanifrz = nanmean(ifrz); stdifrz = nanstd(ifrz)./sqrt(unit.nTrials);
    meanifr = nanmean(ifrm); stdifr = nanstd(ifrm)./sqrt(unit.nTrials);
    
    fname=sprintf('Rec %d %s, %5.3f, %s, Unit %d', ...
        ii, DATA{ii,1}, DATA{ii,2}, DATA{ii,3}, DATA{ii,4});
%     [fh, rasterh, meanh, ~, labelh] = makeRasterIFRplot(fname);
    
    
    
%     plott = respWind./unit.spkSampRate;
%     
%     
%     plotRastersandIFRs(plott, meanh, meanifr, stdifr, rasterh, spkm, 0, [0 0 0]);
%     
    
    
    %plotRastersandIFRs(plott, meanh, meanifrz, stdifrz, rasterh, spkm, 0, [0 0 0]);
    % Test if the response is significantly different from baseline
%     [h, p] = permutationTestFromBaseline(ifrz);
%     sig = h*(meanh.YLim(2)-1); sig(sig==0)=NaN;
%     plot(meanh, plott, sig, 'b', 'LineWidth', 3);
%     
%     DATA{ii,7} = sig;
    
    DATA{ii,8} = ifrz;
    DATA{ii,9} = spkm;
    DATA{ii,10} = ifrm;
    
%     % Mark up
%     plot(rasterh, [0 0], [0 unit.nTrials], 'Color', [.4 .4 .4], 'LineWidth', 2);
%     plot(meanh, [0 0], meanh.YLim, 'Color', [.4 .4 .4], 'LineWidth', 2);
%     rasterh.XLim = [plott(1) plott(end)];
%     meanh.XLim = [plott(1) plott(end)];
%     plot(meanh, meanh.XLim, [unit.zBound(1) unit.zBound(1)], 'r');
%     plot(meanh, meanh.XLim, [unit.zBound(2) unit.zBound(2)], 'r');
%     %plot(meanh, meanh.XLim, [-2 -2], 'r');
%     %plot(meanh, meanh.XLim, [2 2], 'r');
%     xlabel(meanh, 'Time relative to speech onset (sec)', 'FontSize', 16);
%     ylabel(meanh, 'Average Firing Rate (Hz)', 'FontSize', 16);
%     
%     saveas(fh,['figures/',fname,'.pdf'], 'pdf');
%     close(fh);
    
    %%
end

