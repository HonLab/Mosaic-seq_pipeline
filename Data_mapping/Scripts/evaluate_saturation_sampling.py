#!/usr/bin/env python3
import sys
import os
import re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

seed = 10
annotation_file = sys.argv[2]
bam_file = sys.argv[1]

def count_transcript_num(featurecounts_file):
    count_table = [int(line.rstrip().split('\t')[6])
                   for line in open(featurecounts_file, 'r')
                   if not line.startswith(('#', 'Geneid'))]

    transcript_num = len([i for i in count_table if i > 0])
    return transcript_num

def calc_RPKM(featurecounts_file):
    # count_table = [int(line.rstrip().split('\t')[6])
    #                for line in open(featurecounts_file, 'r')
    #                if not line.startswith(('#', 'Geneid'))]

    transcript_id = list()
    transcript_length = list()
    count_table = list()

    with open(featurecounts_file, 'r') as f:
        for line in f:
            if not line.startswith(('#','Geneid')):
                i = line.rstrip().split()
                transcript_id.append(i[0])
                count_table.append(int(i[6]))
                transcript_length.append(int(i[5]))
    total_read_num = sum(count_table)

    rpkm_table = [i * (10 ** 9) / (total_read_num * transcript_length[j])
                  for j, i in enumerate(count_table)]
    return transcript_id, transcript_length, rpkm_table

f_output = open(bam_file + '_saturation_transcript_number.txt','w')

rpkm_sampling_table = list()

for i in range(1,10,1):

    sampling_fraction = str(seed + i/10)
    # output_bam = sys.argv[1][: -4] + '_' + str(sampling_fraction) + '.bam'
    output_bam = re.sub('\.bam$', '', bam_file) + \
                 '_' + str(sampling_fraction) + '.bam'

    cmd_line = ' '.join(['samtools view -bs',
                         sampling_fraction, bam_file, '-o', output_bam])
    os.system(cmd_line)

    output_featurecounts = output_bam[:-4] + '.featureCounts'
    cmd_line = ' '.join(['featureCounts -a', annotation_file, '-o',
                         output_featurecounts, output_bam])
    os.system(cmd_line)

    transcript_num = count_transcript_num(output_featurecounts)
    # print(bam_file[:-4], float(sampling_fraction[2:]),
    #       transcript_num, sep = '\t')
    transcript_num_sampling_result = '\t'.join([str(i) for i in
                                                [re.sub('\.bam$', '', bam_file),
                                                float(sampling_fraction[2:]),
                                                transcript_num]])
    f_output.write(transcript_num_sampling_result + '\n')
    os.remove(output_bam)

    rpkm_sampling_table.append(calc_RPKM(output_featurecounts)[2])

output_featurecounts = bam_file[:-4] + '_10.10.featureCounts'
os.system('featureCounts -a ' + annotation_file  +  ' -o ' +
          output_featurecounts + ' ' + bam_file)

transcript_num = count_transcript_num(output_featurecounts)
# print(bam_file[:-4], float(1), transcript_num, sep = '\t')
transcript_num_sampling_result = '\t'.join([str(i) for i in
                                            [re.sub('\.bam$', '', bam_file),
                                            float(1),
                                            transcript_num]])
f_output.write(transcript_num_sampling_result + '\n')
f_output.close()

unsampling_results = calc_RPKM(output_featurecounts)

rpkm_sampling_table.append(unsampling_results[2])
transcript_id = unsampling_results[0]
transcript_length = unsampling_results[1]


f_output = open(bam_file + '_saturation_RPKM.txt','w')
for j, i in enumerate(transcript_id):

    output_line = list()
    output_line.append(i)
    output_line.append(transcript_length[j])

    for ii in rpkm_sampling_table:
        output_line.append(ii[j])
    output_line = '\t'.join([str(i) for i in output_line])
    f_output.write(output_line + '\n')
f_output.close()


saturation_transcript_num_plot_x = [float(line.split('\t')[1])
                                    for line in
                                    open(bam_file +
                                         '_saturation_transcript_number.txt')]
saturation_transcript_num_plot_y = [int(line.split('\t')[2])
                                    for line in
                                    open(bam_file +
                                         '_saturation_transcript_number.txt')]

with PdfPages(bam_file + '_saturation_transcript_number.pdf') as pdf:

    plt.plot(saturation_transcript_num_plot_x,
             saturation_transcript_num_plot_y,
             c='steelblue')

    # plt.axis([0, 1000, 0, 1])
    plt.xlabel('Percentage of reads')
    plt.ylabel('Transcript number')
    # plt.title()
    pdf.savefig()
    plt.close()


















# print(calc_RPKM(output_featurecounts))
