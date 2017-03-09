clear;

% B: total number of distinct barcode sequences
B = NUM_DISTINCT_BARCODES_ON_PLATE;

data = load('raw.all.pad');

%p_noise = 0.028945645;
p_noise = 0.01;
for cell = 1:size(data,1)
  x = data(cell,:);

  % calculate null model (M_0)
  % P(D|M_0) = (1/B)^(sum x_i) * permute(B, B_0)
  model_probabilities = [null_probability(x, B)];

  B_0 = sum(x > 0);
  for j = 1:B_0

    % calculate probabilities for multinomial
    p_model{j} = [];
    current = 1;
    for k = 1:j
      current = current * 0.9;
      p_model{j} = [p_model{j} current];
    end
    p_model{j} = p_model{j} / sum(p_model{j});

    % compute multinomial probability
    p_multinomial = mnpdf(x(1:j), p_model{j});
    
    % compute null probability for other datapoints
    p_remainder = null_probability( x(j+1:numel(x)), B-j);

    % final pval is product of two pvals
    pval = p_multinomial * p_remainder;
    model_probabilities = [model_probabilities pval];
  end

  % get maximum probability;
  % if this event is not 2x as likely as the null, use the null
  [min, min_index] = max(model_probabilities);
  if (min_index > 0 & ...
      min / model_probabilities(1) < 2)
    min_index = 1;
  end
  model_index(cell) = min_index - 1;
end

model_index = model_index';

save ('raw.all.pad.states.txt', 'model_index', '-ascii');
