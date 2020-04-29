function pars = grad_descent(grad, pars0, learning_rate, ...
                             samples_batch, iter)
%grad_descent Gradient descent. Passes "samples_batch" from 
%iteration T to the gradient function which generates (and utilizes) 
%"samples_batch" at iteration T+1.
%
% Syntax: pars = grad_descent(grad, pars0, learning_rate, ...
%                             samples_batch, iter)
%
% Inputs:
%   grad: A callable [g, samples_batch] = grad(pars, samples_batch),
%         where g is the gradient of the optimized function at "pars".
%   pars0: Initial guess.
%   learning_rate: Learning rate.
%   samples_batch: An array of samples which are being updated by "grad".
%   iter: Number of iterations.
%
% Outputs:
%   pars: Whathever gradient descent converged to in "iter" iterations.

figure
ButtonHandle = uicontrol('Style', 'PushButton', ...
                         'String', 'Stop loop', ...
                         'Callback', 'delete(gcbf)');
% figure
% ButtonReduce = uicontrol('Style', 'PushButton', ...
%                          'String', 'Reduce', ...
%                          'Callback', 'delete(gcbf)'); 
                                   
pars = pars0;
[g, samples_batch] = grad(pars0, samples_batch);
iteration = 0;

% List of average number of sequences
average_occupation = [];
std_occupation = [];

% parameters for the rate adaptation
beta = 0.5;
alpha = 1.2;
l_max = 1;
l_min = 1e-20;
g_old = g;

% reducution rate to make the gradient to converge
r = 0.9;

while iteration <= iter
    
      if ~ishandle(ButtonHandle)
        save('mu_std.mat', 'average_occupation', 'std_occupation', 'samples_batch')
        disp('Loop stopped by user');
        break;
      end
      
%       if ~ishandle(ButtonReduce)
%           learning_rate = r.*learning_rate;
%       end
    pause(0.01)
    
    % methode for adaptation of learning_rates
    pars = pars - learning_rate.*g;
    learning_rate = adapt_rate(learning_rate, g, g_old, alpha, beta, l_max, l_min);
    g_old = g;
    
    % methode without adaptation
    %pars = pars - learning_rate.*g;
    
    [g, samples_batch] = grad(pars, samples_batch);
    [mean, var] = compute_mean_batch(samples_batch);
    
    %mean and std of the occupation mean
    average_occupation = [average_occupation, mean];
    std_occupation = [std_occupation, var];
    
    iteration = iteration + 1;
    disp(['Iteration: ', num2str(iteration), ...
          ', ||Gradient||: ', num2str(norm(g))]);
end
end