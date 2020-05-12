function samples = sample_semi_independent(samples, a, V, n_steps)
%sample_semi_independent Samples semiparametric independent model.
%Applies "n_steps" of Gibbs sampling steps to every row in samples.
%
% Syntax: samples = sample_semi_independent(samples, a, V, n_steps)
%
% Inputs:
%   samples: Initial samples for Gibbs sampling, 
%            size is number of samples x number of neurons.
%   a: The 'biases' of the semiparametric independent model.
%   V: A function V: R -> R of the semiparametric independent model.
%   n_steps: Number of Gibbs sampling steps to apply to every sample.
%
% Outputs:
%   samples: (Approximate) samples from the semiparametric 
%            independent model.

% Initialize.
[M,n] = size(samples);
neuron_id = 1;
% Perform n_steps of Gibbs sampling.
for j = 1:n_steps
    % Calculate energy when current neuron = 0/1.
    %choose a random spring
    neuron_id = random('Discrete Uniform', 402);
    neuron_swap = random('Discrete Uniform', 402);
    
    not_ok = (samples(:, neuron_swap)~=1-samples(:, neuron_id));
    neuron_swap = neuron_swap*ones(M, 1);
    right = neuron_swap;
    left = neuron_swap;
    row = (1:M)';
    
    neuron_id = sub2ind([M, n], row, neuron_id*ones(M,1));
    
    while ~prod(~not_ok)
        right(not_ok) = right(not_ok) + 1;
        left(not_ok) = left(not_ok) - 1;
        
        right(right==403) = 1;
        left(left==0) = 402;
        
        indx_left = sub2ind([M, n], row, left);
        indx_right = sub2ind([M, n], row, right);
        
        ok_left = samples(indx_left) == 1-samples(neuron_id);
        neuron_swap(ok_left) = left(ok_left);
        not_ok(ok_left) = 0;
        ok_right = samples(indx_right) == 1-samples(neuron_id);
        neuron_swap(ok_right) = right(ok_right);
        not_ok(ok_right) = 0;
        
        
        
    end
    
    new_samples = samples;
    new_samples(neuron_id) = 1-samples(neuron_id);
    indx = sub2ind([M, n], row, neuron_swap);
    new_samples(indx) = 1-samples(indx);
    
    E0 = non_linear_energy(samples, -a, V);
    E1 = non_linear_energy(new_samples, -a, V);
    deltaE = E1-E0;
    % Calculate current neurons' spiking probability, flip neurons, and
    %update the energy of samples.
    
    p_spike = 1./(1+exp(deltaE));
    p_spike(E1 <= E0) = 1;
    accept = rand(M,1) < p_spike;
    idx = sub2ind([M, n], row(accept), neuron_swap(accept));
    samples(idx) = new_samples(idx);
    
    samples(neuron_id(accept)) = new_samples(neuron_id(accept));

end

end
