function tonset = detectSpeechOnset(spect, t, f, baseWind)

base = spect(:,baseWind);
basestd = std(base,0,2);
basemean = mean(base,2);
thresh = basemean + 3*basestd;
thresh = repmat(thresh, 1, size(spect,2));
loud = spect > thresh;
howloud = sum(loud);

isloud = howloud > 20;
loudlabel = bwlabel(isloud);
tthresh = 5;
onseti = -1;
for ii = 1:(length(unique(loudlabel))-1)
    sect = find(loudlabel == ii);
    if length(sect) >= tthresh && onseti==-1
        onseti = sect(1);
    end
end


if onseti ~= -1
    tonset = onseti;

    figure;
    pcolor(t, f, spect); shading flat;
    ah = gca;
    hold on; plot([t(tonset) t(tonset)], ah.YLim, 'k');
else
    disp('Nothing detected');
    tonset = -1;
end


