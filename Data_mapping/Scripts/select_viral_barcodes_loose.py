#!/usr/bin/env python3
import sys
import os

'''
Loose:
* cell can have any number of known barcodes
* each barcode is supported by at least 1 read
* cells can have any number of unknown barcodes
'''

# Usage: select_viral_barcodes_loose.py STAMPs_info known_barcode_list

read_count_threshold = 1

known_viral_barcodes = [line.rstrip() for line in open(sys.argv[2], 'r')]

known_viral_barcodes_dict = dict()
output_files = list()

for i in known_viral_barcodes:
    known_viral_barcodes_dict[i] = open(sys.argv[1].replace('.txt', '') +
                                        '.' + i + '.loose' + '.txt', 'w')
    output_files.append(sys.argv[1].replace('.txt', '') +
                        '.' + i + '.loose' + '.txt')

with open(sys.argv[1], 'r') as f:
    for line in f:
        i = line.rstrip().split('\t')

        observed_barcodes =  i[7].split(',')

        if observed_barcodes[0] != 'NA':

            observed_barcode_selected = observed_barcodes[0]
            known_viral_barcodes_dict[observed_barcode_selected].write(
                line.rstrip() + '\n')

for i in known_viral_barcodes:
    known_viral_barcodes_dict[i].close()

for i in output_files:
    if os.path.exists(i) and os.stat(i).st_size == 0:
        os.remove(i)
