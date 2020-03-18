function pars = grad_descent(grad, pars0, learning_rate, ...
                             samples_batch, iter, bins)
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

pars = pars0;
[g, samples_batch] = grad(pars0, samples_batch);
iteration = 0;

figure
ButtonHandle = uicontrol('Style', 'PushButton', ...
                         'String', 'Stop loop', ...
                         'Callback', 'delete(gcbf)');                
beta = 0.5;
alpha = 1.2
l_max = 0.05;
l_min = 1e-9;

g_old = g

while iteration <= iter
    
      if ~ishandle(ButtonHandle)
        disp('Loop stopped by user');
        break;
      end
    pause(0.01)
    
    has_move = true;
    old_lik = likelihood_semi_independent(samples_batch, 100, pars(1:402), pars(403:end), bins, linspace(0,1,10));
    while ~ has_move
            pars_candidate = pars - learning_rate.*g;
            target = likelihood_semi_independent(samples_batch, 100, pars_candidate(1:402), pars_candidate(403:end), bins, linspace(0,1,100));
        
            if target < old_lik
                pars = pars_candidate;
                has_move = true;
            else
                learning_rate = beta*learning_rate;
        end
    end
    
    pars = pars - learning_rate.*g;
    
    [g, samples_batch] = grad(pars, samples_batch);
    learning_rate = adapt_rate(learning_rate, g, g_old, alpha, beta, l_max, l_min);
    g_old = g;

    iteration = iteration + 1;
    
    disp(['Iteration: ', num2str(iteration), ...
          ', ||Gradient||: ', num2str(max(abs(g))), ', cost: ', num2str(old_lik)]);
end
end