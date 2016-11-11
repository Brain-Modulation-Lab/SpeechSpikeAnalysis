%AutoAudioOnsetMarking

%% Mark Response Times
%=======================================
TrialNum = 120;

YLIM = [100 15000];

EventTimes1 = [EventTimes EventTimes(end)+2];

f1 = figure;
for trial=1:TrialNum
    figure(f1); colormap jet;
%     cm = colormap;
%     cm = cm(end:-1:1,:);
%     colormap(cm);
    StimulusEvent1 = SkipEvents + 4*trial;
    StimulusEvent2 = SkipEvents + 4*trial + 1;
    
    AudioTrials{trial} = Audio(round(30000*EventTimes1(StimulusEvent1)): ...
        round(30000*EventTimes1(StimulusEvent2)));
    
    AudioTrials{trial} = 2*(AudioTrials{trial} - mean(AudioTrials{trial}))/ ...
        (max(AudioTrials{trial}) - min(AudioTrials{trial}));
    
    thisAudioenv = AudioEnv(round(30000*EventTimes1(StimulusEvent1)): ...
        round(30000*EventTimes1(StimulusEvent2)));
        
    [s,f,t,p] = spectrogram(AudioTrials{trial}, 512,7/8*512, [], 30000, 'power', 'yaxis'); %compute the spectrogram of the response
    p0 = mean(reshape(p(1:100, :), 1, [])); %baseline signal
    Lp = 20*log10(p./p0); %decibel calculation
    
    hold off
    pcolor(t, f, Lp); shading interp; %set(gca,'YScale','log');
    hold on;
%     plot(AudioTrials{trial})
%     plot(thisAudioenv/100)
    %hold on; text(length(AudioTrials{trial})/10, 0.75*max(AudioTrials{trial}), [num2str(trial), ': ', L{trial}], 'FontSize', 20);
    AudioStart{trial} = 1;
    button=1;
    %display the envelope/waveform 
    while button == 1
        player = audioplayer(AudioTrials{trial}(AudioStart{trial}:end), 30000);
        hold on; h = plot([1 1], YLIM, 'k', 'linewidth', 1);
        %% setup the timer for the audioplayer object
        player.TimerFcn = {@plotMarker, player, h, 30000}; % timer callback function
        player.TimerPeriod = 0.01; % period of the timer in seconds
        play(player);
        
        
        [x,~,button] = ginput(1); %wait for mouse input to mark when analysis should begin
        t_click = round(x);
        earliest_ind = find(t>=t_click, 1,'first');
        
        onset = detectSpeechOnset(Lp, t, f, 1:60, 0); %give the speech calculation the powers in dB
        if onset > 0
            onset_ind = t(onset);
            hold on; plot([onset_ind onset_ind], YLIM, 'k', 'linewidth', 1);
        end
        
        %% setup the timer for the audioplayer object
        %     player.TimerFcn = {@plotMarker, player, h, fs}; % timer callback function
        %     player.TimerPeriod = 0.01; % period of the timer in seconds
        play(player);
    end
    %input('trial done');
end


k=0;
for trial=1:TrialNum
    StimulusEvent1 = SkipEvents + 4*trial;
    if ~isempty(AudioStart{trial})
        k=k+1;
        ResponseTimes(k) = EventTimes1(StimulusEvent1) + AudioStart{trial}/30000;
    end
end

clear AudioEnvTrials;
for trial=1:length(ResponseTimes)
    AudioEnvTrials(:,trial) = AudioEnv((round(30000*ResponseTimes(trial))-30000*2.5):(round(30000*ResponseTimes(trial))+30000*2.5));
end