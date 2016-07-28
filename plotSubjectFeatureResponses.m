%% Section for plotting conditional response
for wl = 1:length(rec_idx)
    for jj=1:length(electrodeList)
        electrode = MerData(rec_idx(wl)).(electrodeList{jj});
        if isfield(electrode, 'Units')
            nUnits = length(MerData(rec_idx(wl)).(electrodeList{jj}).Units);
        else
            nUnits = 0;
        end
        for n=1:nUnits
            unit = MerData(rec_idx(wl)).(electrodeList{jj}).Units(n);
            trode = jj;
            plotUnitRespByStim;
        end
    end
end