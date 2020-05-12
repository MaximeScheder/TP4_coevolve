function fit_nonlinear_indep(data)

%data=load('C:\Users\msche\OneDrive\Documents\EPFL\TP4\DCA_ex\shear05_402_PB.dat');
fields=load('C:\Users\msche\OneDrive\Documents\EPFL\TP4\PlosNonLinFonction\NewFreshStart\mutation_costs\field_conservation.dat');
data_train=data(1:1.0856e+05,:);
data_test=data(1.0856e+05+1:end,:);

n = size(data_train, 2);
in.a = -fields;
%This mess is needed because the parametrization of V depends on the range of energies.
E = -data_train*fields;
%bin_centers = linspace(min(E), max(E), 20);
bin_centers = linspace(-1200, 200, 1000);
x_linear = zeros(length(bin_centers)+2, 1);
x_linear(1) = bin_centers(1) - (bin_centers(2) - bin_centers(1))/2;
x_linear(2) = 1;
in.x = x_linear;
in.bin_centers = bin_centers;

%To play with the following parameters,
%especially the learning rates to find good convergence
options.learning_rate_a = 0.01;
options.learning_rate_x = 0.000005;
options.iter = 100000;
options.M_samples = 50000;
options.gibbs_steps = 50;
%Option to parallelise the code, in the comment with 15 workers
%parpool(15);
out = fit_semi_independent(data_train, in, options);

%Number of samples on which you want to compute the likelihood
N_samples=10;
L_train = likelihood_semi_independent(data_train, N_samples, out.a, out.x, out.bin_centers, linspace(0,1,10));
L_test = likelihood_semi_independent(data_test, N_samples, out.a, out.x, out.bin_centers, linspace(0,1,10));
out.a = -out.a;
%Choose which variables to save, here the output of fit_semi_independent
%(learned a, learned parameters of V and bin_centers, apparently the same
% as the one given in the input?) and the likelihoods on the train and test
%sets. All the other codes are commented by the authors!
save('Grad_basique_fix_n=50.mat', 'out', 'L_train', 'L_test');
end