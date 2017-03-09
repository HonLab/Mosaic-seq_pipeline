#!/usr/bin/env python3
import sys

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


def compare_observed_barcode_with_known_barcodes(observed_barcodes,
                                                 known_barcodes_dict,
                                                 edit_distance_threthold):
    """to be filled"""

    for i in known_barcodes_dict:
        if calc_edit_distance(i, observed_barcodes) <= edit_distance_threthold:
            return i


known_virus_barcodes = ['CGGGTCGAGTGT',
                        'GCGGGGAGTGCA',
                        'GGTGGGAATAAC',
                        'ATGAGTCTCTCA',
                        'CACGAGTGGAAA',
                        'TGAGTCTTCACT',
                        'CGCAAGGTGGGG',
                        'GCATTGGCGGCA',
                        'TCGGCGGTATTA',
                        'CGGGTCGAGTGT']

num_of_STAMPs = 500
edit_distance_threthold = 1

with open(sys.argv[1], 'r') as f:
    selected_STAMPs = [next(f) for x in range(num_of_STAMPs)]

selected_cell_barcodes = [i.rstrip().split('\t')[1] for i in selected_STAMPs]
# print(selected_STAMPs)
# print(selected_cell_barcodes)

for j, i in enumerate(selected_cell_barcodes):
    known_virus_barcodes_dict = {i:int() for i in known_virus_barcodes}
    observed_virus_barcodes = [ii.split('\t')[1]
                               for ii in open('split.bam/' + i +
                                              '_barcodes.txt', 'r')]

    # matched_known_virus_barcodes = list()
    matched_known_virus_barcodes_occurrence = int()
    noise_virus_barcodes = list()
    noise_virus_barcodes_occurrence = int()

    '''
    for iii in observed_virus_barcodes:
        if iii in known_virus_barcodes_dict:
            known_virus_barcodes_dict[iii] += 1
            matched_known_virus_barcodes.append(iii)
            matched_known_virus_barcodes_occurrence += 1
        else:
            noise_virus_barcodes.append(iii)
            noise_virus_barcodes_occurrence += 1
    '''

    for iii in observed_virus_barcodes:
        matched_known_virus_barcode = \
            compare_observed_barcode_with_known_barcodes(
                iii, known_virus_barcodes_dict,
                edit_distance_threthold)

        if matched_known_virus_barcode:
            matched_known_virus_barcodes_occurrence += 1
            known_virus_barcodes_dict[matched_known_virus_barcode] += 1
            # print(matched_known_virus_barcode)
        else:
            noise_virus_barcodes.append(iii)
            noise_virus_barcodes_occurrence += 1
            # print(noise_virus_barcodes)
    # print(known_virus_barcodes_dict)

    matched_known_virus_barcodes_num = len([i for i in known_virus_barcodes_dict
                                            if known_virus_barcodes_dict[i]])
    noise_virus_barcodes_num = len(set(noise_virus_barcodes))

    matched_known_barcodes_results = list()

    sorted_matched_barcodes = \
        sorted(known_virus_barcodes_dict.items(),
               key=lambda
               known_virus_barcodes_dict: known_virus_barcodes_dict[1],
               reverse = True)
    # print(sorted_matched_barcodes)

    if matched_known_virus_barcodes_num:
        for jj in sorted_matched_barcodes:
            jj = jj[0]

            if known_virus_barcodes_dict[jj]:
                matched_known_barcodes_results.append(
                    [jj, known_virus_barcodes_dict[jj]])
                matched_known_barcodes_sequence = \
                    [m[0] for m in matched_known_barcodes_results]
                matched_known_barcodes_count = \
                    [m[1] for m in matched_known_barcodes_results]

                matched_known_barcodes_output = \
                    ','.join(matched_known_barcodes_sequence) + '\t' + \
                    ','.join([str(i) for i in matched_known_barcodes_count])
    else:
        matched_known_barcodes_output = '\t'.join(['NA'] * 2)


    print(selected_STAMPs[j].rstrip(),
          matched_known_barcodes_output,
          matched_known_virus_barcodes_num,
          matched_known_virus_barcodes_occurrence,
          noise_virus_barcodes_num,
          noise_virus_barcodes_occurrence, sep = '\t')
