function [cons, Delta_E, Delta_F] = mutation_cost(Mend, filename)
% The function here will compute the cost of mutation in our problem
%   parameters :
%       a : fields for our model with energy E = -a*X
%       X : A configuration of the system
%       filename : path of the file where are the fields
% It will return the cost of single mutations in
%the mechanical network Delta_F and inferred model Delta_E.
% [cons, DE,~]=mutation_cost(15000)

%% Load
%In a you need to load the field inferred via the *independent* model with
%non-linearity
addpath('C:\Users\msche\OneDrive\Documents\EPFL\TP4\PlosNonLinFonction\NewFreshStart')
parameters = load(filename);
a = parameters.out.a;
X = load('C:\Users\msche\OneDrive\Documents\EPFL\TP4\DCA_ex\shear05_402_PB.dat');
listpos = importdata('C:\Users\msche\OneDrive\Documents\EPFL\TP4\PlosNonLinFonction\NewFreshStart\mutation_costs\listpos.dat'); % list of link positions in decreasing order of importance for mutation costs, to identify neutral positions


%% Initialize

delta_E = zeros(numel(X(1,:)),1);
delta_F = zeros(numel(X(1,:)),1);
Delta_E = zeros(numel(X(1,:)),1);
Delta_F = zeros(numel(X(1,:)),1);

%% Definition of the non-linearity
V = @(E) monotone(parameters.out.x', parameters.out.bin_centers, E);

%% Mutation costs in the mechanical network

% for m = 1:Mend
%     
%     fitness_0 = FitnessCOOPPB(1,1,X(m,:));
%     X_mutated = X(m,:);
% 
%     for i=1:length(X(1,:))
%         X_mutated(i)=1-X_mutated(i);
%         fitness_i = FitnessCOOPPB(1,1,X_mutated);
%         delta_F(i) = - (fitness_0 - fitness_i);
%         X_mutated = X(m,:);
%     end
%     
%     Delta_F = Delta_F + delta_F;
% end

%% Mutatation costs in the inferred model

for m = 1:Mend
    
    grid = X(m,:); % Need to flip configurations for the gauge of parameters
    gridp = grid(:);

    for r=1:length(X(1,:))
        
        % Performing mutations with respect to reference position of links
        % where mutations cost very little in order to perform single
        % mutations via swap
        for t=1:length(listpos)    
            if (gridp(listpos(t)) == 1 - gridp(r)) && (listpos(t) ~= r)
                if t > 262
                    disp('Index listpos > 262')
                end
                pos_ref1 = listpos(t);
                break 
            end
        end
  
        gridpp = gridp;
        gridpp(r) = gridp(pos_ref1);
        %gridpp(pos_ref1) = gridp(r);
        delta_E(r) = non_linear_energy(gridpp,a, V) - non_linear_energy(gridp,a, V); % effective single site cost in the inferred model
    end
    
    Delta_E = Delta_E + delta_E;
end

%% Resulting average over configurations
 Delta_E = Delta_E/Mend;
 Delta_F = Delta_F/Mend;
 cons = mean(X);        %measure the conservation of springs in the sequences

end