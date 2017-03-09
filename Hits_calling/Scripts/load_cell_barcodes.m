function [cell_barcode] = load_cell_barcodes(experiment)
  cell_bc_file         = fopen(sprintf('../step1.data/Normalized_data/combined.%s/all_sgRNA.uniq.cell_barcodes',...
                                       experiment));
  cell_bc_array        = textscan (cell_bc_file, '%s', 'Delimiter', '\t'); 
  cell_barcode         = cell_bc_array{1};
end
