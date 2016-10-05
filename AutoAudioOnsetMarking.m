%AutoAudioOnsetMarking

%% Mark Response Times
%=======================================
TrialNum = 120;

EventTimes1 = [EventTimes EventTimes(end)+2];

f1 = figure;
for trial=1:TrialNum
    figure(f1);
    StimulusEvent1 = SkipEvents + 4*trial;
    StimulusEvent2 = SkipEvents + 4*trial + 1;
    
    AudioTrials{trial} = Audio(round(30000*EventTimes1(StimulusEvent1)): ...
        round(30000*EventTimes1(StimulusEvent2)));
    
    AudioTrials{trial} = 2*(AudioTrials{trial} - mean(AudioTrials{trial}))/ ...
        (max(AudioTrials{trial}) - min(AudioTrials{trial}));
    
    thisAudioenv = AudioEnv(round(30000*EventTimes1(StimulusEvent1)): ...
        round(30000*EventTimes1(StimulusEvent2)));
    
    hold off
    t_aud = ((1:length(AudioTrials{trial}))-1)/30000;
    plot(t_aud, AudioTrials{trial})
    hold on; plot(t_aud, thisAudioenv/100)
    hold on; text(t_aud(end)/10, 0.75*max(AudioTrials{trial}), [num2str(trial), ': ', L{trial}], 'FontSize', 20);
    AudioStart{trial} = 1;
    button=1;
    %display the envelope/waveform 
    player = audioplayer(AudioTrials{trial}(AudioStart{trial}:end), 30000); play(player);
    [x,~,button] = ginput(1); %wait for mouse input to mark when analysis should begin
    t_click = round(x);
    [s,f,t,p] = spectrogram(AudioTrials{trial}, 512,7/8*512, [], 30000, 'power', 'yaxis'); %compute the spectrogram of the response
    earliest_ind = find(t>=t_click, 1,'first');
    p0 = mean(reshape(p(1:100, :), 1, [])); %baseline signal
    Lp = 20*log10(p./p0); %decibel calculation
    onset = detectSpeechOnset(Lp, t, f, 1:60, earliest_ind); %give the speech calculation the powers in dB
    onset_ind = t(onset)*30000;
    hold on; plot([onset_ind onset_ind], [-1 3], 'k');
    figure(f1); plot([t(onset) t(onset)], [-1 3], 'k');
    disp('trial done');
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