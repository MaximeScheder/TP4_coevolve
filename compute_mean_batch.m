function [average, var] = compute_mean_batch(sample_batch)
samples = [];
    if ndims(sample_batch) > 2
        samples = sample_batch(:,:,1);
        for i=2:size(sample_batch, 3)
            samples = [samples ; sample_batch(:,:,i)];
        end
    end

average = mean(mean(samples, 2));
var = std(mean(samples, 2));
end
