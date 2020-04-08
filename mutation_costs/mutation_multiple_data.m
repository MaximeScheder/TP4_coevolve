function mutation_multiple_data(Mend, colorplo)
 
true_F = load('single_costs_mechanical.dat');

[~, model_E, ~] = mutation_cost(Mend);

% Compute correlation coefficient
c = corrcoef(model_E,true_F);

figure
hold on
plot(model_E, true_F,'ro')
x = linspace(0, max(max(model_E)), 100);
fmax = max(max(true_F));
m = fmax/max(max(model_E));
plot(x,m*x,'b:')
xlabel('\Delta E')
ylabel('\Delta F')
title('Mutation cost (correlation coefficient \rho=' + c(1,2) + ')')
grid on