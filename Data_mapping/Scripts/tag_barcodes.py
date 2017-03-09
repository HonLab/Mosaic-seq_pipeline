#!/usr/bin/env python3
import sys
import pysam

barcode_base_quality_threshold = 10
num_of_bases_below_quality_threshold = 1

barcode_counter = [int()] * 21

f_unaligned_bam = pysam.AlignmentFile(sys.argv[1], check_sq = False)
f_unaligned_tagged_bam = pysam.AlignmentFile(sys.argv[1][: -4] +
    '_CellMolecule_tagged.bam', 'wb', template = f_unaligned_bam)

for read in f_unaligned_bam.fetch(until_eof = True):
    if read.is_read1:
        cell_barcode = read.seq[0: 12]
        # cell_barcode_quality = read.qual[0:12]
        # cell_barcode_phred_score = [ord(i)-33 for i in cell_barcode_quality]
        molecule_barcode = read.seq[12: 20]

        barcode_quality_phred_score = [ord(i) - 33 for i in
            read.qual.decode()[: 20]]

        num_of_low_quality_bases = int()
        for i,j in enumerate(barcode_quality_phred_score):
            if j < barcode_base_quality_threshold:
                num_of_low_quality_bases += 1
                barcode_counter[num_of_low_quality_bases] += 1

        if not num_of_low_quality_bases:
            barcode_counter[0] += 1

    elif read.is_read2:
        # print(read.tags)

        if num_of_low_quality_bases <= num_of_bases_below_quality_threshold:
            read.flag = 4
            # print(read.tags)

            read.tags = [('RG', 'A'.encode('utf-8'), 'Z'.encode('utf-8')),
                ('XC', cell_barcode.encode('utf-8')),
                ('XM', molecule_barcode.encode('utf-8'))]
            f_unaligned_tagged_bam.write(read)

f_unaligned_bam.close()
f_unaligned_tagged_bam.close()

f_output = open(sys.argv[1][: -4] + '_CellMolecule_tagged_summary.txt','w')
for i in range(len(barcode_counter)):
    f_output.write(str(i) + '\t' + str(barcode_counter[i]) + '\n')
f_output.close()
