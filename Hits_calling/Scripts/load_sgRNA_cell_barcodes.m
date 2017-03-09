function [all_sgRNA_cell_barcodes] = load_sgRNA_cell_barcodes(experiment, ...
                                                         region,...
                                                         sgRNA)

  prefix = '../step1.data/Normalized_data';

  % get all replicates
  dir_search_string = sprintf('%s/combined.%s.*',...
                              prefix,...
                              experiment);
  dir_files = dir(dir_search_string);

  all_sgRNA_cell_barcodes = [];
  for rep_ind = 1:length(dir_files)
    ind = findstr(dir_files(rep_ind).name, '.rep');
    rep = dir_files(rep_ind).name(ind+1 : ind+1+3);

    full_name = sprintf('%s/combined.%s.%s/%s.%s.cell_barcodes',...
                        prefix,...
                        experiment,...
                        rep,...
                        region,...
                        sgRNA);

    sgRNA_pos_file       = fopen(full_name,'r');
    sgRNA_pos_cell_array = textscan (sgRNA_pos_file, '%s', 'Delimiter','\n');  
    sgRNA_cell_barcode   = sgRNA_pos_cell_array{1};
    
    all_sgRNA_cell_barcodes = [all_sgRNA_cell_barcodes ; sgRNA_cell_barcode];
  end
end
