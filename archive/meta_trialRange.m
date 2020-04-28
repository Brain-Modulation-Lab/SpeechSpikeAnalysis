elecs = {'Cent', 'Post', 'Med', 'Lat'};

DATA={};
k=0;
subjects={'JA120215','CS010616','RC021616','MS041916', ...
            'PR042116','JH061616','GY063016','BA080416',...
            'WK082316','TU081116','RK102516','PH101816',...
            'KS102716'};
%subjects={'MS041916'};

nRec = 0;

for s=1:length(subjects)
    
    SubjectID = subjects{s};
    
    load(['Rec/',SubjectID,'Rec']);
    
    load(['MerData_60/',SubjectID,'MerData']);
    
    analyzeSubjectUnits_retrieveTiming;
    
    for idx=length(MerData):-1:1
        for e = 1:length(elecs)
            if isfield(MerData(idx),elecs{e})
                if ~isempty(MerData(idx).(elecs{e}))
                    for u=1:length(MerData(idx).(elecs{e}).Units)
                        if ~isempty(MerData(idx).(elecs{e}).Units(u).trial)
                            fprintf('%s (%f) %s: unit %d\n', SubjectID, Rec(idx).Depth, elecs{e}, u);
                            k=k+1;
                            DATA(k,:) = {SubjectID, Rec(idx).Depth, elecs{e}, u, ...
                                MerData(idx).(elecs{e}).Units(u), ...
                                MerData(idx).ResponseTimes, ...
                                MerData(idx).OffsetTimes, ...
                                MerData(idx).(elecs{e}).Units(u).sig, ...
                                MerData(idx).(elecs{e}).Units(u).ifrz, ...
                                MerData(idx).(elecs{e}).Units(u).spkm, ...
                                MerData(idx).(elecs{e}).Units(u).ifrm};
                        end
                    end
                end
            end
        end
    end
    
end

% for i=1:length(DATA)
% temp = zeros(1, 120);
% temp(DATA{i,5}(1):DATA{i,5}(2)) = DATA{i,6};
% DATA{i,8} = temp;
% end
% 
% for i=1:length(DATA)
% temp = zeros(1, 120);
% temp(DATA{i,5}(1):DATA{i,5}(2)) = DATA{i,7};
% DATA{i,9} = temp;
% end
% 
% for i=1:length(DATA)
% temp = zeros(1, 120);
% temp(DATA{i,5}(1):DATA{i,5}(2)) = ones(size(DATA{i,7}));
% DATA{i,10} = temp;
% end
