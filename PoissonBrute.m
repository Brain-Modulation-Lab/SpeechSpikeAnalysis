function [ Sb, Sp, comb, t ] = PoissonBrute( D, fs, r, N0 )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%Tavg = length(D)/(fs*nnz(D));
%r = single(1/Tavg);
r = single(r);
indD = find(D);
indDmax = length(indD);
Ncomb = sum(1:(indDmax-N0+1));
comb = cell([1 Ncomb]);
k=0;
for i=1:(indDmax-N0+1)
    for j=(i+N0-1):indDmax
        k=k+1;
        comb{k} = i:j;
    end
end

n = cellfun(@(x) single(length(x)-1), comb);
T = cellfun(@(x) single((indD(x(end))-indD(x(1)))/fs), comb);

%S = arrayfun(@(x,y) -log10(poisscdf(x,r*y,'upper')), n, T);

Sb = zeros(Ncomb,1);
Sp = zeros(Ncomb,1);

parfor i=1:Ncomb
    Sb(i) = -log10(poisscdf(n(i),r*T(i),'upper'));
    Sp(i) = -log10(poisscdf(n(i),r*T(i)));
end

t = indD/fs;

end

