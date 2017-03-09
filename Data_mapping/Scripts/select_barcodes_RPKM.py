#!/usr/bin/env python3
import sys

# read selected barcodes as a list
selected_barcodes = [i.split('\t')[1] for i in open(sys.argv[1], 'r')]

selected_columns = list()
selected_barcodes_reordered = list()
with open(sys.argv[2], 'r') as f:
    for line in f:
        if line.startswith('#'):
            barcode_list = line.rstrip().split('\t')
            barcode_list[0] = barcode_list[0][1:]

            for j, i in enumerate(barcode_list):
                if i in selected_barcodes:
                    selected_columns.append(j)
                    selected_barcodes_reordered.append(i)
            print('#Geneid',
                  '\t'.join(selected_barcodes_reordered), sep = '\t')

        elif not line.startswith('#'):
            i = line.rstrip().split('\t')
            output_exprs = list()

            for j, ii in enumerate(i[1:]):
                if j in selected_columns:
                    output_exprs.append(ii)
            print(i[0], '\t'.join([str(iii) for iii in output_exprs]),
                  sep = '\t')
