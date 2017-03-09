function [sample_mean, sample_cov, sample_cpm, sample_median] = get_mean_cov_ver2(sample)
  sample_sum    = sum(sample);
  sample_norm   = ones(length(sample), 1) * sample_sum / 1000000;
  sample_cpm    = sample ./ sample_norm;
  sample_mean   = mean(sample_cpm,2);
  sample_median = median(sample_cpm,2);
  sample_std    = std(sample_cpm, 0, 2);
  sample_cov    = sample_std ./ sample_mean;
end
