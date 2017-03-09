clear;

%specify file names
matrix_file = 'combined.all.uniq.vals.no_header';

%read the data
matrix      = load(matrix_file);

one_cell_cpm  = (sum(matrix,2) + 1) / sum(sum(matrix)) * 1000000;

[matrix_mean, matrix_cov, matrix_cpm, matrix_median] = get_mean_cov_ver2(matrix);

%median_array = sum(matrix,2) ./ sum(sum(matrix_cpm)) * 1000000;

[gene_num, cell_num] = size(matrix);

normalized_matrix = [];
for i = 1:cell_num
  normalized_matrix(:,i) = matrix_cpm(:,i) ./ one_cell_cpm;
end

%save ('normalized.all.matrix.txt.no_header', 'normalized_matrix', '-ascii');
dlmwrite('normalized.all.matrix.txt.no_header', normalized_matrix,...
         'delimiter', '\t', 'precision', '%.3f');
