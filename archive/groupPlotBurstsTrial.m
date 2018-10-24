
figdir = '/Users/brainmodulationlab/Dropbox (Brain Modulation Lab)/Speech Project/MerData/figures_DBS2000';


for idx=1:length(DATA)
    clear BB;

    fs = 1000;

    figtitle = sprintf('Cue2 Trial: Rec %d %s, %5.3f, %s, Unit x, (%s)', ...
        idx, DATA(idx).SubjectID, DATA(idx).depth, DATA(idx).elec, DATA(idx).grade);
    
%     RT = DATA(idx).RT(DATA(idx).trange);
%     [~,order] = sort(RT);
%     
%     plot_rasterSpeechTaskTrialBurstsOrder( [figdir,'/TrialBurstMarked'], figtitle, ...
%         DATA(idx).Cue2, DATA(idx).Trial, DATA(idx).trange, [], order, RT, 'b');

    figtitle = sprintf('SpeechOnset: Rec %d %s, %5.3f, %s, Unit %d, (%s)', ...
        idx,  DATA(idx).SubjectID,  DATA(idx).depth,  DATA(idx).elec, DATA(idx).unit,  DATA(idx).grade);
    plot_rasterSpeechTask( [figdir,'/SpeechOnsetISITS'], figtitle,  DATA(idx).SpeechOnset,  DATA(idx).trange );

end