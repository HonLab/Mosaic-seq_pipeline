clear;

lib_index  = [1 2 3 4 5 6 7 8 9 10 13 14 15 16 17 18 19 20 21 22 23 ...
             32];

PIM1_index = 16285;
ACTB_index = 18365;

select_index = PIM1_index;

final_array       = [];
final_norm_array  = [];
final_group       = [];
final_norm_group  = [];

group_name = {};

for i = 1:length(lib_index)

  sprintf('i = %d', i);
  current_index       = lib_index(i);
  current_norm_file   = sprintf('./sorted/RX134-%d.dropseq_pipe/normalized.all.matrix.txt.no_header',current_index);
  current_file        = sprintf('./sorted/RX134-%d.dropseq_pipe/combined.all.uniq.vals.no_header',current_index);

  current_matrix      = load(current_file);
  matrix_norm_cpm = load(current_norm_file);

  [matrix_mean, matrix_cov, matrix_cpm, matrix_median] = get_mean_cov_ver2(current_matrix);
%  [matrix_norm_mean, matrix_norm_cov, matrix_norm_cpm, matrix_norm_median] = get_mean_cov_ver2(current_norm_matrix);
  
  PIM1_array          = matrix_cpm (select_index,:);
  PIM1_norm_array     = matrix_norm_cpm (select_index,:);

  final_array         = [final_array PIM1_array];
  final_norm_array    = [final_norm_array PIM1_norm_array];

  final_group         = [final_group i*ones(1,length(PIM1_array))];
  final_norm_group    = [final_norm_group i*ones(1,length(PIM1_norm_array))];
  
  name          = sprintf('RX134-%d',current_index);
  group_name(i) = {name};

end



save ('test_data.mat');

subplot(2,1,1);
boxplot (final_array, final_group);
title ('PIM1 Before Normalization');
set (gca, 'XTickLabel', group_name, 'XTickLabelRotation', 45);

subplot(2,1,2);
boxplot (final_norm_array, final_norm_group);
title ('PIM1 After Normalization');
set (gca, 'XTickLabel', group_name, 'XTickLabelRotation', 45);


print ('Normalization_vs_nonNorm_PIM1.pdf', '-dpdf');
