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

% Perform n_steps of Gibbs sampling.
    % Go through all samples to perform a change that conserve the number
    % of springs
    for s = 1:size(samples, 1)
        % Calculate energy when current neuron = 0/1.
        %choose a random spring
        sequence = samples(s, :);
        %select two random position
        neuron_id = random('Discrete Uniform', 402);
        neuron_swap = random('Discrete Uniform', 402);
        no_pair = true;
        
        %ensure here that the swap conserve the number of springs
        if sequence(neuron_id) == (1-sequence(neuron_swap))
            no_pair = false;
        else
            left = neuron_swap;
            right = neuron_swap;
        end
        % Explore right and left from before to find a correct spring
        while no_pair
            left = left - 1;
            right = right + 1;
            
            if right == 402
                right = 1;
            end
            
            if left == 0
                left = 402;
            end
            
            if sequence(neuron_id) == (1-sequence(left))
                neuron_swap = left;
                no_pair = false;
            elseif sequence(neuron_id) == (1-sequence(right))
                neuron_swap = right;
                no_pair = false;
            end
        end
        
        new_seq = sequence;
        new_seq(neuron_id) = 1-new_seq(neuron_id);
        new_seq(neuron_swap) = 1-new_seq(neuron_swap);

        deltaE = non_linear_energy(new_seq, -a, V)-non_linear_energy(sequence, -a, V);

        if deltaE <= 0
            samples(s, :) = new_seq;
        else
            p_spike = 1./(1+exp(deltaE));
            flipped = rand() < p_spike;
            if flipped
                samples(s,:) = new_seq;
            end
        end
    end
end

