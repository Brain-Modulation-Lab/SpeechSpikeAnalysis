figdir = '/Users/brainmodulationlab/Dropbox (Brain Modulation Lab)/Speech Project/MerData/figures_DBS2000';

% respInterval = [-1 2];

spkSampRate = 1000;
filtSD = 25; % in SAMPLES
minTsig = 100; % minimum significant FR change (SAMPLES)
alpha = 0.05; % significance threshold for FR changes
%Nobs = round((respInterval(2)-respInterval(1))/(2*filtSD/spkSampRate));


alphabet = 'abcdefghijklmopqrstuvwxyz';

% DATA = [];

Rec = 0;

for s=1:length(subjects)
    
    
    dat_files = dir([subjects{s},'/Preprocessed Data/',subjects{s},'*.mat']);
    
    for f=1:length(dat_files)
        fprintf('Processing %s...\n', dat_files(f).name);
        
        S = load([subjects{s},'/Preprocessed Data/',dat_files(f).name], ...
            'spikes', 'trials', 'track_labels', 'skips', 'RecDuration', 'depth');
        
        tokens = strsplit(dat_files(f).name,'_');
        session = cellfun(@str2num, regexp(tokens{2}, '[1-9]', 'match'));

%         % add RecDuration variable if one does not already exist
%         if isfield(S, 'macro') && ~isfield(S, 'RecDuration')
%             RecDuration = length(S.macro)/S.nfs;
%             save([subjects{s},'/Preprocessed Data/',dat_files(f).name], 'RecDuration', '-append');
%         else
%             fprintf('  no macro.\n');
%         end
        
        if isfield(S, 'spikes')
            
            for elec=1:length(S.spikes)
                fprintf('   %s electrode.\n', S.track_labels{elec});
                if ~isempty(S.spikes(elec).unit)
                    
                    for u=1:length(S.spikes(elec).unit)
                        
                        unit = S.spikes(elec).unit(u);
                        
                        fprintf('      Unit %d: %s (%s)\n', u, unit.RecType, unit.grade);
                        
                        %% Generate IFR
                        t = (1:round(spkSampRate*S.RecDuration))/spkSampRate; %timeseries
                        D = zeros(length(t), 1);
                        spkinds = round(unit.ts.*spkSampRate);
                        spkinds(spkinds < 1) = 1;
                        D(spkinds) = 1; %binary spike vector
                        %sets up the gaussian filter for spike density estimation
                        %standard  deviation of the gaussian smoothing
                        stdg = filtSD; % this standard deviation is in samples
                        filtx = -4*stdg:1:4*stdg;
                        filty = normpdf(filtx,0,stdg);
                        filty = filty/sum(filty);
                        %20ms gaussian smoothing
                        IFR = conv2(D(:), filty(:), 'same');
                        %normalize to the number of spikes in the trial as in Gabbiani et al. 1999
                        IFR = IFR/(sum(IFR)/spkSampRate);
                        IFR = IFR .* sum(D);
                        
                        %% Generate ISITS
                        ISITS = ISIts(D,spkSampRate);
                        
                        %% Determine trial range
                        trange = 1:60;
                        trange = setdiff(trange, find(isnan(S.trials.SpOnset)));
                        trange = setdiff(trange, union(S.trials.BaseRejectNoise, S.trials.BaseRejectSpk));
                        trange = setdiff(trange, S.trials.ResponseReject.all);
                        
                        % Calculate the response interval based on RT and
                        % Speech duration
                        RT = arrayfun(@minus, S.trials.SpOnset, S.trials.CommandStim(1:min([60 length(S.trials.CommandStim)]))');
                        SpDur = arrayfun(@minus, S.trials.SpEnd, S.trials.SpOnset);
                        
                        trange = setdiff(trange, find((S.trials.BaseBack-1)<unit.tstart));
                        if strcmp(unit.tend,'end')
                            tend = spkSampRate*length(D);
                        else
                            tend = unit.tend;
                        end
                        trange = setdiff(trange, find((S.trials.SpOnset+2)>tend));
                        
                        %% Compute per-trial data
                        if length(trange)>10
                            
                            Rec = Rec+1;
                            
                            respInterval = round(1000*[-nanmean(RT) nanmean(SpDur)+0.5])/1000;
                            Nobs = round((respInterval(2)-respInterval(1))/(2*filtSD/spkSampRate));
                            SpeechOnset = SpikeTrialData(spkSampRate, IFR, ISITS, D, trange, S.trials.BaseBack, S.trials.SpOnset, respInterval, respInterval, ...
                                alpha, Nobs, minTsig);
                            
                            respInterval = round(1000*[-0.5 nanmean(RT)])/1000;
                            respIntervalDisp = round(1000*[-0.5 nanmean(RT)])/1000;
                            Nobs = round((respInterval(2)-respInterval(1))/(2*filtSD/spkSampRate));
                            Cue = SpikeTrialData(spkSampRate, IFR, ISITS, D, trange, S.trials.BaseBack, S.trials.CommandStim, respInterval, respIntervalDisp, ...
                                alpha, Nobs, minTsig);
                            
                            respInterval = round(1000*[0 nanmean(RT)+nanmean(SpDur)+0.5])/1000;
                            respIntervalDisp = round(1000*[-0.5 nanmean(RT)+nanmean(SpDur)+0.5])/1000;
                            Nobs = round((respInterval(2)-respInterval(1))/(2*filtSD/spkSampRate));
                            Cue2 = SpikeTrialData(spkSampRate, IFR, ISITS, D, trange, S.trials.BaseBack, S.trials.CommandStim, respInterval, respIntervalDisp, ...
                                alpha, Nobs, minTsig);
                            
                            %% Whole Trial data
                            clear Trial;
                            
                            Trial.start_samples = round(spkSampRate*S.trials.Cue1Stim);
                            Trial.end_samples = round(spkSampRate*S.trials.ITIStim);
                            Trial.DD = arrayfun(@(x) ...
                                D(Trial.start_samples(x):Trial.end_samples(x)), ...
                                trange, 'uniformoutput', false);
                            Trial.timeCue = arrayfun(@(x,y) ...
                                linspace(S.trials.Cue1Stim(x)-S.trials.CommandStim(x),S.trials.ITIStim(x)-S.trials.CommandStim(x),length(Trial.DD{y})), ...
                                trange, 1:length(trange), 'uniformoutput', false);
                            Trial.timeSpOnset = arrayfun(@(x,y) ...
                                linspace(S.trials.Cue1Stim(x)-S.trials.SpOnset(x),S.trials.ITIStim(x)-S.trials.SpOnset(x),length(Trial.DD{y})), ...
                                trange, 1:length(trange), 'uniformoutput', false);
                            Trial.idxEventTimes = arrayfun(@(x) ...
                                round(spkSampRate*([S.trials.Cue2Stim(x) S.trials.CommandStim(x) S.trials.SpOnset(x) S.trials.SpEnd(x)] - ...
                                S.trials.Cue1Stim(x)))+1, ...
                                trange, 'uniformoutput', false);
                            
                            DATA(Rec).SubjectID = subjects{s};
                            DATA(Rec).session = session;
                            DATA(Rec).depth = S.depth;
                            DATA(Rec).elec = S.track_labels{elec};
                            DATA(Rec).unit = u;
                            DATA(Rec).grade = unit.grade;
                            DATA(Rec).RecType = unit.RecType;
                            DATA(Rec).trange = trange;
                            DATA(Rec).D = D;
                            DATA(Rec).IFR = IFR;
                            DATA(Rec).ISITS = ISITS;
                            DATA(Rec).PreCue = arrayfun(@minus, S.trials.CommandStim(1:min([60 length(S.trials.CommandStim)]))', ...
                                S.trials.Cue2Stim(1:min([60 length(S.trials.CommandStim)]))');
                            DATA(Rec).RT = RT;
                            DATA(Rec).SpDur = SpDur;
                            
                            DATA(Rec).Trial = Trial;
                            DATA(Rec).SpeechOnset = SpeechOnset;
                            DATA(Rec).Cue = Cue;
                            DATA(Rec).Cue2 = Cue2;

                            DATA(Rec).trials = S.trials;


                            figtitle = sprintf('SpeechOnset: Rec %d %s, %5.3f, %s, Unit %d, (%s)', ...
                                Rec, subjects{s}, S.depth, S.track_labels{elec}, u, unit.grade);
                            plot_rasterSpeechTaskRT( [figdir,'/SpeechOnset'], figtitle, SpeechOnset, RT, SpDur, trange );
                            
                            figtitle = sprintf('Cue: Rec %d %s, %5.3f, %s, Unit %d, (%s)', ...
                                Rec, subjects{s}, S.depth, S.track_labels{elec}, u, unit.grade);
                            plot_rasterSpeechTaskRT( [figdir,'/Cue'], figtitle, Cue, RT, SpDur, trange );
                            
                            figtitle = sprintf('Cue2: Rec %d %s, %5.3f, %s, Unit %d, (%s)', ...
                                Rec, subjects{s}, S.depth, S.track_labels{elec}, u, unit.grade);
                            plot_rasterSpeechTaskRT( [figdir,'/Cue2'], figtitle, Cue2, RT, SpDur, trange );

                        end
                                    
                        
                        clear unit t D spkinds stdg filtx filty IFR ISITS ...
                            tend trange  SpeechOnset Cue ITI RT SpDur;
                        
                    end
                    
                else
                    fprintf('      no units.\n');
                end
            end
        else
            fprintf('  no spikes.\n');
        end
        
    end
    
end
        

