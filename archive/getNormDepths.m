
for i = 1:length(NormDepth)
    if i~= 24 && i~=25
    temp = strsplit(NormDepth(i).session, '_');
    
    subject = temp{1};
    wordlist = cellfun(@str2num, regexp(temp{2}, '[1-9]', 'match'));
    
    session = SessionTable( ...
                find(cellfun(@(x,y) strcmp(subject, x) && wordlist==y, ...
                {SessionTable(:).subject}, {SessionTable(:).wordlist}))).session;
            
    idx = intersect(find(arrayfun(@(x) strcmp(subject, x.subject), TrialData)), find(arrayfun(@(x) x.session==session, TrialData)));
    
    TrialData(idx).MicroNormDepth = NormDepth(i).MicroNormDepth;
    end
end