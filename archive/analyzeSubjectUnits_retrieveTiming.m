coding_dir = '/Users/Witek/Dropbox/Richardson Lab/Speech Project/Decoding/coding_files';
spkDir = 'spkData';

% Analyze Subject's Units
eval([SubjectID,'subjectInfo;']);
%load(recFN);
respInterval = [-2 2.5];
%MerData = struct([]);
dbg= 1;

for ii = 1:length(rec_idx)
    
    fprintf('Rec %d.\n', rec_idx(ii));
    
    %% Load EventTimes from coding files
    load([coding_dir,'/',SubjectID,'_Session',num2str(taskNum(ii))]);

    nTrials = floor((length(EventTimes)-SkipEvents)/4);
    EventTimes1 = [EventTimes EventTimes(end)+2]; % add extra event time in case session ended in the middle of trial
    
    %% Load ResponseTimes from SpeechTrials 
    AudioStart = SpeechTrials( ...
        find(strcmp([subjectName,'_Session',num2str(taskNum(ii)),'.mat'], ...
        SpeechTrials(:,1))),12);
    
    ResponseTimes=[];
    
    for trial=1:nTrials
        StimulusEvent1 = SkipEvents + 4*trial;
        if trial<=length(AudioStart) && ~isempty(AudioStart{trial})
            ResponseTimes(trial) = EventTimes1(StimulusEvent1) + AudioStart{trial};
        else 
            ResponseTimes(trial) = NaN;
        end
    end
    
    %% Load OffsetTimes from SpeechTrials 
    AudioStop = SpeechTrials( ...
        find(strcmp([subjectName,'_Session',num2str(taskNum(ii)),'.mat'], ...
        SpeechTrials(:,1))),13);
    
    OffsetTimes=[];
    
    for trial=1:nTrials
        StimulusEvent1 = SkipEvents + 4*trial;
        if trial<=length(AudioStop) && ~isempty(AudioStop{trial})
            OffsetTimes(trial) = EventTimes1(StimulusEvent1) + AudioStop{trial};
        else 
            OffsetTimes(trial) = NaN;
        end
    end
    
        MerData(rec_idx(ii)).OffsetTimes = OffsetTimes;
    
end