
figdir = '/Users/brainmodulationlab/Dropbox (Brain Modulation Lab)/Speech Project/MerData/figures_DBS2000';


for idx=1:length(DATA)
    clear PP;

    fs = 1000;

    figtitle = sprintf('Cue2 Trial: Rec %d %s, %5.3f, %s, Unit x, (%s)', ...
        idx, DATA(idx).SubjectID, DATA(idx).depth, DATA(idx).elec, DATA(idx).grade);
    
    RT = DATA(idx).RT(DATA(idx).trange);
    [~,order] = sort(RT);
    
    plot_rasterSpeechTaskTrialPauseOrder( [figdir,'/TrialPauseMarked'], figtitle, ...
        DATA(idx).Cue2, DATA(idx).Trial, DATA(idx).trange, [], order, RT, 'b');

end