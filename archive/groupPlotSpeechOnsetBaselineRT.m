figdir = '/Users/brainmodulationlab/Dropbox (Brain Modulation Lab)/Speech Project/MerData/figures_DBS2000';

spkSampRate = 1000;

for idx=1:27%length(DATA)

    fs = 1000;

    SpeechOnsetBaseline =  DATA(idx).SpeechOnset;
    
    %% Trial data
    SpeechOnsetBaseline.respInterval = round(1000*[-nanmean(DATA(idx).RT)-0.5 nanmean(DATA(idx).SpDur)+0.5])/1000;
    SpeechOnsetBaseline.DD = cell2mat(arrayfun(@(x) ...
        DATA(idx).D((SpeechOnsetBaseline.trial_samples(x)+round(spkSampRate*SpeechOnsetBaseline.respInterval(1))): ...
        (SpeechOnsetBaseline.trial_samples(x)+round(spkSampRate*SpeechOnsetBaseline.respInterval(2)))), ...
        DATA(idx).trange, 'uniformoutput', false))';
    SpeechOnsetBaseline.IFRdata = cell2mat(arrayfun(@(x) ...
        DATA(idx).IFR((SpeechOnsetBaseline.trial_samples(x)+round(spkSampRate*SpeechOnsetBaseline.respInterval(1))): ...
        (SpeechOnsetBaseline.trial_samples(x)+round(spkSampRate*SpeechOnsetBaseline.respInterval(2)))), ...
        DATA(idx).trange, 'uniformoutput', false))';
    SpeechOnsetBaseline.ISITSdata = cell2mat(arrayfun(@(x) ...
        DATA(idx).ISITS((SpeechOnsetBaseline.trial_samples(x)+round(spkSampRate*SpeechOnsetBaseline.respInterval(1))): ...
        (SpeechOnsetBaseline.trial_samples(x)+round(spkSampRate*SpeechOnsetBaseline.respInterval(2)))), ...
        DATA(idx).trange, 'uniformoutput', false))';
    
    SpeechOnsetBaseline.meanifr =  mean(SpeechOnsetBaseline.IFRdata,1);
    SpeechOnsetBaseline.stdifr = std(SpeechOnsetBaseline.IFRdata,[],1)/sqrt(size(SpeechOnsetBaseline.IFRdata,1));
    
    SpeechOnsetBaseline.sig_excit = cat(2, nan(1,0.5*spkSampRate), SpeechOnsetBaseline.sig_excit);
    SpeechOnsetBaseline.sig_inhib = cat(2, nan(1,0.5*spkSampRate), SpeechOnsetBaseline.sig_inhib);
    
    DATA(idx).SpeechOnsetBaseline = SpeechOnsetBaseline;
    
    figtitle = sprintf('SpeechOnset: Rec %d %s, %5.3f, %s, Unit x, (%s)', ...
        idx, DATA(idx).SubjectID, DATA(idx).depth, DATA(idx).elec, DATA(idx).grade);
    
%     RT = DATA(idx).RT(DATA(idx).trange);
%     [~,order] = sort(RT);
%     
%     plot_rasterSpeechTaskTrialBurstsOrder( [figdir,'/TrialBurstMarked'], figtitle, ...
%         DATA(idx).Cue2, DATA(idx).Trial, DATA(idx).trange, [], order, RT, 'b');

    RT = DATA(idx).RT(DATA(idx).trange);
    %RT = DATA(idx).RT;
    SpDur =  DATA(idx).SpDur(DATA(idx).trange);
    
    figtitle = sprintf('SpeechOnset: Rec %d %s, %5.3f, %s, Unit %d, (%s)', ...
        idx,  DATA(idx).SubjectID,  DATA(idx).depth,  DATA(idx).elec, DATA(idx).unit,  DATA(idx).grade);
    plot_rasterSpeechTaskRT( [figdir,'/SpeechOnsetIFR'], figtitle,  DATA(idx).SpeechOnsetBaseline,  RT, SpDur, DATA(idx).trange );

end