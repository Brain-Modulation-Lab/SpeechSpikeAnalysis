GLM = {};

for rec=1:length(DATA)

    fprintf('Rec %d\n', rec);
    GLManalysis;
    GLM{rec,1} = B;
    GLM{rec,2} = psig;
    
end