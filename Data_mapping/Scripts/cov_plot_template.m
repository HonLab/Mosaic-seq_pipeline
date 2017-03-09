clear;

%
% original data
%
all_data = importdata('MATRIX_FILE', '\t', 1);
sc = all_data.data;
[num_genes, num_cells] = size(sc);

  sc_sum  = sum(sc);
  sc_norm = ones(length(sc), 1) * sc_sum / 1000000;
  sc_cpm  = sc ./ sc_norm;
  sc_mean = mean(sc_cpm,2);
  sc_std  = std(sc_cpm, 0, 2);
  sc_cov  = sc_std ./ sc_mean;

  x = log2(sc_mean+1);
  y = sc_cov;
  yy = smooth(x,y,0.1,'loess');
  
  ind = find(sc_mean > 100);
  cpm100_mean_cov = mean(sc_cov(ind));

  hold on;
  plot(x, y, 'r.');
  plot(x, yy, 'k.', 'MarkerSize', 1)
  plot([0 12], [1 1], 'g-')
  plot([log2(100) log2(100)], [0 100], 'k:');
  text(8, 16, sprintf('CoV (cpm > 100): %.2f', cpm100_mean_cov));
  hold off;

  xlim([0 12]);
  ylim([0 20]);

  title('NAME');

  ylabel('CoV')
  xlabel('log2(CPM + 1)')

print -dpdf OUT_FILE
