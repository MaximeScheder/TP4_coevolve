% Cette fonction sert à ploter les résultats sous forme jolie

function plot_result(filename, data, Mend)
addpath('C:\Users\msche\OneDrive\Documents\EPFL\TP4\PlosNonLinFonction\NewFreshStart\mutation_costs')

out = load(filename);
h = out.out.a;
bins = out.out.bin_centers;
x = out.out.x;

E = -data*h;
E_min = min(E);
E_max = max(E);

if E_min < 0
    E_min = E_min*2;
else
    E_min = E_min * 0.5;
end

if E_max < 0
    E_max = E_max*0.5;
else
    E_max = E_max*2;
end

E = linspace(E_min, E_max, 1000);
V = @(E) monotone(x', bins, E);

figure
title(['N_P = ', filename(16:end-4)])

true_F = load('single_costs_mechanical.dat');

[~, model_E, ~] = mutation_cost(Mend, data, filename);

% Compute correlation coefficient
c = corrcoef(model_E,true_F);

subplot(1,2,1)
hold on
scatter(model_E, true_F, 'ko', 'MarkerFaceColor', 'k')
x = linspace(0, max(max(model_E)), 100);
fmax = max(max(true_F));
m = fmax/max(max(model_E));
plot(x,m*x,'b:', 'LineWidth', 2)
xlabel('\Delta E')
ylabel('\Delta F')
title(['Mutation cost (correlation coefficient \rho=' num2str(c(1,2)) ')'])
grid on

subplot(1,2,2)
plot(E, V(E), 'k-', 'LineWidth', 2)
xlabel('E')
ylabel('V(E)')
grid on

figure
E_boundary = linspace(-1000, 1000, 2000);
plot(E_boundary, V(E_boundary), 'k-', 'LineWidth', 2)
xlabel('E')
ylabel('V(E)')
title('V(E) long range')
grid on

end