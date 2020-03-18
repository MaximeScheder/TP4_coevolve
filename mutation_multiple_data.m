function mutation_multiple_data(a, data, colorplo)
    
model_E = zeros(size(data));
true_F = zeros(size(data));

n = size(data,1);

for i=1:n
    [model_E(i,:), true_F(i,:)]= mutation_cost(a, data(i,:));
end

figure
hold on
plot(model_E, true_F,'ro')
x = linspace(0, max(max(model_E)), 100);
fmax = max(max(true_F));
m = fmax/max(max(model_E));
plot(x,m*x,'b:')
xlabel('\Delta E')
ylabel('\Delta F')
title('Mutation Cost')
grid on