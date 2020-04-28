
respInterval = [-2 2.5];


for ii = 1:length(DATA)
    

    fprintf('Rec %s, depth %5.3f, %s, Unit %d.\n', DATA{ii,1}, DATA{ii,4}, DATA{ii,5}, DATA{ii,6});
    
    SessionID = [DATA{ii,1},'_Session',num2str(DATA{ii,2}),'.mat'];
    
    idx = find(strcmp(SessionID, SpeechTrials(:,1)));
    
    
   
    DATA{ii,7} = unit;
    DATA{ii,8} = ResponseTimes;
    DATA{ii,9} = sig;
    DATA{ii,10} = ifrz;
    DATA{ii,11} = spkm;
    DATA{ii,12} = ifrm;
    
    
    %%
end

