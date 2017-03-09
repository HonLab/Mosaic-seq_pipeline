clear;

PIM1_ind  = 16285;
ACTB_ind  = 18366;
FADS1_ind = 27995;
SMYD3_ind = 5053;
%region_id = {'PIM1.SE2.HS2',...
%             'PIM1.SE2.HS1',...
%             'PIM1.SE2.HS10',...
%             'PIM1.SE2.HS11',...
%             'PIM1.SE2.HS12',...
%             'PIM1.SE2.HS3',...
%             'PIM1.SE2.HS4',...
%             'PIM1.SE2.HS5',...
%             'PIM1.SE2.HS6',...
%             'PIM1.SE2.HS9'};
region_id = {'SMYD3.SE2.HS10'};
experiments = {'Batch_A.rep1', 'Batch_A.rep2', 'Batch_A.rep3'};

%
% load all data
%
data_cpm = [];
cell_barcode = {};
for experiment = 1:length(experiments)
  experiment

  % read expression data
  current_cpm = load(sprintf('./combined.%s/normalized.all.matrix.txt.no_header',...
                             experiments{experiment}));
  data_cpm = [data_cpm current_cpm];
  

  % read cell barcodes for matrix
  cell_bc_file         = fopen(sprintf('./combined.%s/all_sgRNA.uniq.cell_barcodes',...
                                       experiments{experiment}));
  cell_bc_array        = textscan (cell_bc_file, '%s', 'Delimiter', '\t'); 
  cell_barcode = [cell_barcode cell_bc_array{1}'];
end
[num_genes, num_cells] = size(data_cpm);

%
% process sgRNAs
%
for i=1:length(region_id)
  % get all sgRNA barcodes
  sgRNA_pos_barcode = [];
  for experiment = 1:length(experiments)
    full_name            = sprintf('./combined.%s/%s.sgRNA1.cell_barcodes',...
                                   experiments{experiment},...
                                   region_id{i});
    sgRNA_pos_file       = fopen(full_name,'r');
    sgRNA_pos_cell_array = textscan (sgRNA_pos_file, '%s', 'Delimiter','\n');  
    sgRNA_pos_barcode    = [sgRNA_pos_barcode sgRNA_pos_cell_array{1}'];
  end

  index_sgRNA = get_cell_barcode_index(cell_barcode, sgRNA_pos_barcode);

  sgRNA_bool = zeros(num_cells, 1);
  sgRNA_bool(index_sgRNA) = 1;

  ind = SMYD3_ind;
  cpm = data_cpm(ind, :);

  % permute the data to reduce artifacts
  rng(0);
  perm_order = randperm( length(cpm) );
  cpm_perm = cpm(perm_order);
  sgRNA_bool_perm = sgRNA_bool(perm_order);

  % sort
  [cpm_sorted, sort_ind] = sort(cpm_perm);
  sgRNA_bool_sorted = sgRNA_bool_perm(sort_ind);
  
  marked_lt_median = sum(sgRNA_bool_sorted(1:floor(num_cells/2)));
  marked_total = sum(sgRNA_bool_sorted);
  p = hygecdf(marked_lt_median - 1,...
              num_cells,...
              marked_total,...
              floor(num_cells/2),...
              'upper');

  fprintf('region = %s, hygecdf(%d, %d, %d, %d) = %f, percent marked_lt_median = %f\n',...
          region_id{i},...
          marked_lt_median - 1,...
          num_cells,...
          marked_total,...
          floor(num_cells/2),...
          p,...
          marked_lt_median / marked_total * 100);

  % plot
  FigHandle = figure('Position', [100, 100, 500, 100]);
  hold on;
  x = find ( sgRNA_bool_sorted == 1 );
  for k = 1:length(x)
    plot([ x(k) x(k) ], [0 1], 'k-');
  end
  hold off;

  xlabel('cells sorted by expr');
  ylabel('sgRNA');
  
  xlim([0 length(cpm)]);

  % save file
%  print('-dpdf', sprintf('sgRNA_bool_sorted.combined.%s.cpm.ordered.pdf',...
%                         region_id{i}));
end
