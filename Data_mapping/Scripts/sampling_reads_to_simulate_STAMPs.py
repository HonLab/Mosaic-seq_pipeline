#!/usr/bin/env python3
import sys
import os
import random
import pysam

# input bam file
bam_file = sys.argv[1]
# read the reads number of each real STAMP
STAMPs_info = [int(i.split('\t')[5]) for i in open(sys.argv[2], 'r')]
# number of simulated STAMPs
num_of_STAMPs = int(sys.argv[3])
# gtf annotation file for featureCounts
annotation_file = sys.argv[4]

def index_bam(bam_file):
    """to be filled"""
    f = pysam.AlignmentFile(bam_file)
    if not f.has_index():
        os.system('samtools index ' + bam_file)
    f.close()

index_bam(bam_file)
total_reads_num = sum([int(i.split('\t')[2]) for i in pysam.idxstats(bam_file)])

fractions = [i / total_reads_num for i in STAMPs_info[: num_of_STAMPs]]

bam_file_list = list()
for j, i in enumerate(fractions):
    sampling_fraction = i
    output_bam = 'sampled_reads_' + str(range(num_of_STAMPs)[j] + 1) + '.bam'
    bam_file_list.append(output_bam)

    seed = random.randrange(1, 1000000)
    samtools_cmd_line = ' '.join(['samtools view -bs',
                                  str(seed + sampling_fraction),
                                  bam_file, '-o', output_bam])
    os.system(samtools_cmd_line)

featureCounts_cmd_line = ' '.join(['featureCounts -a', annotation_file, '-o',
                                   'sampling_results',
                                   ' '.join(bam_file_list)])
os.system(featureCounts_cmd_line)
