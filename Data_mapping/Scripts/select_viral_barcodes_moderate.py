#!/usr/bin/env python3
import sys
import os

'''
Moderate:
* cell can have any number of known barcodes
* each barcode is supported by at least 3 reads
* cell can have an unknown barcode if the unknown barcode is common (i.e., from
an error in the barcode generation step)
'''

# Usage: select_viral_barcodes_moderate.py STAMPs_info known_barcode_list

read_count_threshold = 3
unknown_barcodes_threshold = 1

known_viral_barcodes = [line.rstrip() for line in open(sys.argv[2], 'r')]

known_viral_barcodes_dict = dict()
output_files = list()

for i in known_viral_barcodes:
    known_viral_barcodes_dict[i] = open(sys.argv[1].replace('.txt', '') +
                                        '.' + i + '.moderate' + '.txt', 'w')
    output_files.append(sys.argv[1].replace('.txt', '') +
                        '.' + i + '.moderate' + '.txt')

with open(sys.argv[1], 'r') as f:
    for line in f:
        i = line.rstrip().split('\t')

        observed_barcodes =  i[7].split(',')

        if observed_barcodes[0] != 'NA':
            read_count_known_barcode = [int(ii) for ii in i[8].split(',')]
            num_unknown_barcodes = int(i[11])

            if num_unknown_barcodes <= unknown_barcodes_threshold:
                for index, value in enumerate(read_count_known_barcode):
                    if value >= read_count_threshold:
                        observed_barcode_selected = observed_barcodes[index]
                        known_viral_barcodes_dict[
                            observed_barcode_selected].write(line.rstrip() +
                                                             '\n')

for i in known_viral_barcodes:
    known_viral_barcodes_dict[i].close()

for i in output_files:
    if os.path.exists(i) and os.stat(i).st_size == 0:
        os.remove(i)
