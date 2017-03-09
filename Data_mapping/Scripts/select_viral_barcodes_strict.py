#!/usr/bin/env python3
import sys
import os

'''
Strict:
* cell has exactly 1 known barcode
* each barcode is supported by at least 5 reads
* cell has no unknown barcodes
'''

# Usage: select_viral_barcodes_strict.py STAMPs_info known_barcode_list

known_barcodes_threshold = 1
read_count_threshold = 5
unknown_barcodes_threshold = 0

known_viral_barcodes = [line.rstrip() for line in open(sys.argv[2], 'r')]

known_viral_barcodes_dict = dict()
output_files = list()

for i in known_viral_barcodes:
    known_viral_barcodes_dict[i] = open(sys.argv[1].replace('.txt', '') +
                                        '.' + i + '.strict' + '.txt', 'w')
    output_files.append(sys.argv[1].replace('.txt', '') +
                        '.' + i + '.strict' + '.txt')

with open(sys.argv[1], 'r') as f:
    for line in f:
        i = line.rstrip().split('\t')

        observed_barcodes =  i[7].split(',')

        if len(observed_barcodes) <= known_barcodes_threshold and \
                observed_barcodes[0] != 'NA':
            read_count_known_barcode = int(i[8])
            num_unknown_barcodes = int(i[11])

            if read_count_known_barcode >= read_count_threshold and \
                    num_unknown_barcodes <= unknown_barcodes_threshold:
                # print(line.rstrip())

                observed_barcode_selected = observed_barcodes[0]
                known_viral_barcodes_dict[observed_barcode_selected].write(
                    line.rstrip() + '\n')

for i in known_viral_barcodes:
    known_viral_barcodes_dict[i].close()

for i in output_files:
    if os.path.exists(i) and os.stat(i).st_size == 0:
        os.remove(i)
