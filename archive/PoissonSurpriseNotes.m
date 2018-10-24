% % 
% % for LIM=1:1000
% % 
% %     P(LIM) = PoissonSurprise(10, 1, 10, LIM+100);
% %     
% % end
% 
%  idx=58;
% 
% Stim = DATA(idx).SpeechOnset;
% clear BB;
% trial = 1;
% LIM=100;
% fs = 1000;
% %r = Stim.basemeanIFR;
% Tavg = length(Stim.DD(trial,:))/(fs*nnz(Stim.DD(trial,:)));
% r = 1/Tavg;
% indD = find(Stim.DD(trial,:));
% indDmax = length(indD);
% 
% % figure; plot(Stim.IFRdata(trial,:))
% % hold on; plot(indD'*[1 1], ylim, 'k')
% % hold on; plot(xlim, r*[1 1])
% 
% %% first pass thresholding based on ISI
% N1 = 5;
% T1 =[];
% for i=1:(indDmax-N1-1)
%     T1(i) = length(Stim.DD(trial,indD(i):(indD(i+2)-1)))/(fs*(N1-1));
% end
% P1 = bwconncomp(T1<Tavg/2);
% P1idx = cellfun(@(x) min(x), P1.PixelIdxList);
% 
% % figure;
% % for i=1:P1.NumObjects
% %     hold on; area(P1.PixelIdxList{i}, ...
% %         ones(size(P1.PixelIdxList{i})), ...
% %         'facecolor', 'y')
% % end
% % hold on; plot(T1)
% % hold on; plot(P1idx'*[1 1], ylim, 'r')
% 
% 
% %% Calculate the Poisson Surprise for each burst
% S = [];
% for i=1:P1.NumObjects
%     S(i) = -log10(PoissonSurprise(r, (indD(max(P1.PixelIdxList{i})+1) - ...
%         indD(min(P1.PixelIdxList{i})))/fs, length(P1.PixelIdxList{i}), LIM));
% end
% 
% % for i=1:P1.NumObjects
% %     hold on; text(P1idx(i), 0.5, num2str(S(i), '%3.1f'));
% % end
% 
% %%%%
% %%  Second pass: add spikes to the end of each burst
% B=[];
% for i=1:P1.NumObjects
%     
%     B(i).ii = P1.PixelIdxList{i};
%     B(i).S = S(i);
%     
%     if max(B(i).ii)+1 < indDmax
%         
%         temp = -log10(PoissonSurprise(r, (indD(max(B(i).ii)+2) - ...
%             indD(min(B(i).ii)))/fs, length(B(i).ii)+1, LIM));
%         
%         while temp > B(i).S && (max(B(i).ii)+1) <= length(indD)
%             
%             fprintf('Trailing spike added...\n');
%             
%             B(i).ii = cat(1, B(i).ii, max(B(i).ii)+1);
%             B(i).S = temp;
%             
%             if max(B(i).ii)+1 < indDmax
%                 temp = -log10(PoissonSurprise(r, (indD(max(B(i).ii)+2) - ...
%                     indD(min(B(i).ii)))/fs, length(B(i).ii)+1, LIM));
%             end
%         end
%         
%     end
%     
% end
% 
% % for i=1:P1.NumObjects
% %     hold on; area(B(i).ii, ...
% %         ones(size(B(i).ii))/4, ...
% %         'facecolor', 'g')
% % end
% % 
% % for i=1:length(B)
% %     hold on; text(B(i).ii(1), 0.3, num2str(B(i).S, '%3.1f'));
% % end
% 
% %% Third pass: subract spikes from the beginning of each burst
% for i=1:length(B)
% 
%     temp = -log10(PoissonSurprise(r, (indD(max(B(i).ii)+1) - ...
%     (indD(min(B(i).ii)+1)))/fs, length(B(i).ii)-1, LIM));
%     
%     while temp > B(i).S && length(B(i).ii)>3
%         
%         fprintf('Leading spike subtracted...\n');
%         
%         B(i).ii = setdiff(B(i).ii, min(B(i).ii));
%         B(i).S = temp;
%         
%         temp = -log10(PoissonSurprise(r, (indD(max(B(i).ii)+1) - ...
%             (indD(min(B(i).ii))+1))/fs, length(B(i).ii)-1, LIM));
%     end
%     
% end
% 
% % for i=1:P1.NumObjects
% %     hold on; area(B(i).ii, ...
% %         ones(size(B(i).ii))/5, ...
% %         'facecolor', 'c')
% % end
% % 
% % for i=1:length(B)
% %     hold on; text(B(i).ii(1), 0.15, num2str(B(i).S, '%3.1f'));
% % end
% 
% %% apply threshold to S value
% %Bthresh = 2;
% 
% %% ...or take the largest S value
% [Smax, i] = max([B(:).S]);
% B = B(i);
% fprintf('Trial %d: S = %4.2f\n', trial, Smax)
% 
% BB(trial) = B;
% 
% % final result
% figure;
% for i=1:length(B)
%     hold on; area(indD(cat(1, B(i).ii, max(B(i).ii)+1)), ...
%         max(Stim.IFRdata(trial,:))*ones(size(B(i).ii)+[1 0]), ...
%         'facecolor', 'y')
% end
% 
% for i=1:length(B)
%     hold on; text(indD(B(i).ii(1)), max(Stim.IFRdata(trial,:))/2, ...
%         num2str(B(i).S, '%3.1f'));
% end
% hold on; plot(Stim.IFRdata(trial,:))
% hold on; plot(indD'*[1 1], ylim, 'k')
% hold on; plot(xlim, r*[1 1])


%%
%%%% combination method
figdir = '/Users/Witek/Dropbox (Brain Modulation Lab)/Speech Project/MerData/figures_DBS2000';


for idx=1:length(DATA)
    Stim = DATA(idx).SpeechOnset;
    clear BB;
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
            S = arrayfun(@(x,y) -log10(poisscdf(x,r*sum(y),'upper')), n, T);
            [Smax, Sind] = max(S);
            Sind = Sind(1);
            ii = comb{Sind};
            B.ii = setdiff(ii, ii(end))';
            B.S = Smax;
            
        else
            B.ii = [];
            B.S = 0;
        end
        BB(trial) = B;
    end
    
    figtitle = sprintf('SpeechOnset PSmarked: Rec %d %s, %5.3f, %s, Unit x, (%s)', ...
        idx, DATA(idx).SubjectID, DATA(idx).depth, DATA(idx).elec, DATA(idx).grade);
    
    plot_rasterSpeechTaskBursts( [figdir,'/SpeechOnsetPSmarked'], figtitle, Stim, DATA(idx).trange, BB, [] );
    
    DATA(idx).SpeechOnset.BB = BB;
end