#!/usr/bin/env python3
import sys
import pysam
import editdistance

def calc_edit_distance(string_a, string_b):
    """to be filled"""
    if len(string_a) > len(string_b):
        string_a,string_b = string_b,string_a

    distances = range(len(string_a) + 1)

    for index2, char2 in enumerate(string_b):
        new_distances = [index2+1]
        for index1, char1 in enumerate(string_a):
            if char1 == char2:
                new_distances.append(distances[index1])
            else:
                new_distances.append(1 + min((distances[index1],
                                              distances[index1+1],
                                              new_distances[-1])))
        distances = new_distances
    return distances[-1]

# criteria
# the boundary of real and noisy STAMPs
real_STAMPs_threshold = int(sys.argv[3])
noisy_STAMPs_threshold = int(sys.argv[4])
edit_distance_threthold = 1
minimal_mapping_quality = 10

bam_file = sys.argv[2]

'''
with open(sys.argv[1], 'r') as f:
    selected_STAMPs_real = [next(f).split('\t')[1] for x in range(real_STAMPs_threshold)]
'''

with open(sys.argv[1], 'r') as f:
    selected_STAMPs = [next(f).split('\t')[1]
                       for x in range(noisy_STAMPs_threshold)]

selected_STAMPs_real = selected_STAMPs[: real_STAMPs_threshold]
selected_STAMPs_intermediate = selected_STAMPs[real_STAMPs_threshold :
                                               noisy_STAMPs_threshold]


f = pysam.AlignmentFile(bam_file, 'rb')
output_bam_noisy = bam_file.replace('.bam', '') + '_mapQ' + \
                   str(minimal_mapping_quality) + '_below' + \
                   str(noisy_STAMPs_threshold) + '.bam'
f_output_noisy = pysam.AlignmentFile(output_bam_noisy, 'wb', template = f)

output_bam_real = bam_file.replace('.bam', '') + '_mapQ' + \
                  str(minimal_mapping_quality) + '_above' + \
                  str(real_STAMPs_threshold) + '.bam'
f_output_real = pysam.AlignmentFile(output_bam_real, 'wb', template = f)


for read in f.fetch(until_eof = True):
    if read.mapping_quality >= minimal_mapping_quality:
        real_STAMPs_read_indicator = int()
        for i in read.tags:
            if i[0] == 'XC':
                cell_barcode = i[1]

        for i in selected_STAMPs_real:
            edit_distance = editdistance.eval(i, cell_barcode)

            if edit_distance <= edit_distance_threthold:
                real_STAMPs_read_indicator += 1
                break

        if real_STAMPs_read_indicator:
            f_output_real.write(read)

        else:
            noisy_STAMPs_read_indicator = int()
            for i in selected_STAMPs_intermediate:
                edit_distance = editdistance.eval(i, cell_barcode)

                if edit_distance <= edit_distance_threthold:
                    noisy_STAMPs_read_indicator += 1
                    break
            if not noisy_STAMPs_read_indicator:
                f_output_noisy.write(read)

f.close()
f_output_noisy.close()
f_output_real.close()
