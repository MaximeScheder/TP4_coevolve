function [delta_E, delta_F] = mutation_cost(a, X)
% The function here will compute the cost of mutation in our problem
%   paramters :
%       a : fields for our model with energy E = -a*X
%       X : A configuration of the system
% It will return a graph mapping the true fitness against the empirical
% cost.

if iscolumn(a)
    a = a';
end

if iscolumn(X)
    X = X';
end

delta_E = zeros(numel(X),1);
delta_F = zeros(numel(X),1);
dlmwrite('X.dat',X,'\t')

E = isingenergy(X',a');
FitnessCOOPPB(1,1,'X.dat','Dist_cons_X.dat')
Consensus = load('Dist_cons_X.dat');
F = -Consensus(:,2);
X_mutated = X;

for i=1:length(X)
    
    X_mutated(i)=1-X_mutated(i);
    dlmwrite('X.dat',X_mutated,'\t')
    delta_E(i) = isingenergy(X_mutated',a')-E;
    FitnessCOOPPB(1,1,'X.dat', 'Dist_cons_X.dat')
    Consensus = load('Dist_cons_X.dat');
    delta_F(i) = -Consensus(1,2)-F;
    X_mutated = X;
end

delta_E = abs(delta_E);
delta_F = abs(delta_F);

end