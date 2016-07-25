function waves = assembleSpikeWaves(inds, V, wind)
% function waves = assembleSpikeWaves(inds, V, wind)

waves = zeros(length(inds), length(wind));
for ii = 1:length(inds)
    ty = V(inds(ii)+wind);
    waves(ii,:) = ty;
end

    