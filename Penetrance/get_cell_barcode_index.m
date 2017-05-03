function [index_sgRNA] = get_cell_barcode_index(cell_barcode, sgRNA_pos_barcode)
  for i = 1 : length(sgRNA_pos_barcode)
    index_sgRNA(i) = find(strcmp(cell_barcode, sgRNA_pos_barcode(i)));
  end
end
