function P = PoissonSurprise( r, T, n, LIM )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

P = exp(-r*T)*sum(arrayfun(@(i) (r*T)^i/factorial(i), n:LIM));

%P = poisscdf(n, r*T, 'upper');

end

