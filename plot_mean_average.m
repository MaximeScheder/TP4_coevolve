function plot_mean_average()
data = load('mu_std.mat');
mean = data.average_occupation;
std_ = data.std_occupation;

std_ = [0, std_];
mean = [0.8806, mean];
x = 0:1:(length(mean)-1);

figure
hold on
grid on
l = plot(x, 360/408*ones(size(x)), 'r-', 'LineWidth', 2);
m = plot(x, mean, 'b:', 'LineWidth', 2);
xlabel('Number of steps')
ylabel('average of occupation')
inbetween = [mean-std_, fliplr(mean+std_)];
x2 = [x, fliplr(x)];
f = fill(x2, inbetween, 'k');
alpha(f, 0.2)
legend(l, '360/408 springs')

end