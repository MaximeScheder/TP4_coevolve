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
    
    no_pair = (samples(:, neuron_swap)==1-samples(neuron_id));
    neuron_swap = neuron_swap*ones(402, 1);
    right = neuron_swap;
    left = neuron_swap;
    
    while prod(no_pair)
        right(no_pair) = right + 1;
        left(no_pair) = left - 1;
        
        right(right==403) = 1;
        left(left==0) = 402;
        
        neuron_swap(no_pair) = comp_left(no_pair);
        comp_left = samples(:, left) == 1-samples(:, neuron_id);
        no_pair(no_pair.*comp_left) = 1;
        neuron_swap(no_pair)
        comp_right = samples(:, right) == 1-samples(:, neuron_id);
        no_pair(no_pair.*comp_right) = 1;
        
        
        
    end
    
    new_samples = samples;
    new_samples(:, neuron_id) = 1-samples(:, neuron_id);
    new_samples(:, neuron_swap) = 1-samples(:, neuron_swap)
    
    deltaE = non_linear_energy();
    E0 = non_linear_energy(samples, -a, V);
    E1 = non_linear_energy(new_samples, -a, V);
    deltaE = E1-E0;
    % Calculate current neurons' spiking probability, flip neurons, and
    %update the energy of samples.
    
    p_spike = 1./(1+exp(V(E1) - V(E0)));
    p_spike(E1 <= E0) = 1;
    accept = rand(M,1) < p_spike;
    samples(accept,:) = new_samples(accept, :);

end

end

