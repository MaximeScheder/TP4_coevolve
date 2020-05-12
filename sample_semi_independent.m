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
E_samples = samples*a;
% Perform n_steps of Gibbs sampling.
for j = 1:n_steps
    % Calculate energy when current neuron = 0/1.
    %choose a random spring
    neuron_id = random('Discrete Uniform', 402);
    neuron_swap = random('Discrete Uniform', 402);
    if neuron_swap == neuron_id
        neuron_swap = neuron_swap + 1;
    end
    if neuron_swap == 403
        neuron_swap = 1;
    end
    
    new_seq = samples;
    new_seq(:, neuron_id) = 1-samples(:, neuron_id);
    new_seq(:, neuron_swap) = 1-samples(:,neuron_swap);
    deltaE = a(neuron_id);
    E0 = E_samples - samples(:, neuron_id).*deltaE + (1 - samples(:, neuron_swap)).*a(neuron_swap);
    E1 = E_samples + (1 - samples(:, neuron_id)).*deltaE - samples(:, neuron_swap).*a(neuron_swap);
    % Calculate current neurons' spiking probability, flip neurons, and
    %update the energy of samples.
    V_E1 = V(E1);
    V_E0 = V(E0);
    p_spike = 1./(1+exp(V_E1 - V_E0));
    %p_spike(V_E1 < V_E0) = 1;
    flipped = rand(M,1) < p_spike;
    to_keep = samples(:, neuron_id)==samples(:, neuron_swap);
    temp = samples(:, neuron_id);
    flipped(to_keep) = temp(to_keep);
    temp_swap = samples(:, neuron_swap);
    temp_swap(~to_keep) = 1-flipped(~to_keep);
    E_samples = E_samples + (flipped - samples(:,neuron_id)).*deltaE + (temp_swap-samples(:, neuron_swap))*a(neuron_swap);
    samples(:,neuron_id) = flipped;
    samples(:, neuron_swap) = temp_swap;
    % Move on to the next neuron.
    neuron_id = neuron_id + 1;
end

end

