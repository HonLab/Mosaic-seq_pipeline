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


known_viral_barcodes =[line.rstrip() for line in open(sys.argv[2], 'r')]

num_of_STAMPs = 2000
edit_distance_threthold = 1

with open(sys.argv[1], 'r') as f:
    selected_STAMPs = [next(f) for x in range(num_of_STAMPs)]

selected_cell_barcodes = [i.rstrip().split('\t')[1] for i in selected_STAMPs]
# print(selected_STAMPs)
# print(selected_cell_barcodes)

for j, i in enumerate(selected_cell_barcodes):
    known_viral_barcodes_dict = {i:int() for i in known_viral_barcodes}
    observed_viral_barcodes = [ii.split('\t')[1]
                               for ii in open('bam/split.bam/' + i +
                                              '_barcodes.txt', 'r')]

    # matched_known_viral_barcodes = list()
    matched_known_viral_barcodes_occurrence = int()
    noise_viral_barcodes = list()
    noise_viral_barcodes_occurrence = int()
    noise_viral_barcodes_dict = dict()

    '''
    for iii in observed_viral_barcodes:
        if iii in known_viral_barcodes_dict:
            known_viral_barcodes_dict[iii] += 1
            matched_known_viral_barcodes.append(iii)
            matched_known_viral_barcodes_occurrence += 1
        else:
            noise_viral_barcodes.append(iii)
            noise_viral_barcodes_occurrence += 1
    '''

    for iii in observed_viral_barcodes:
        matched_known_viral_barcode = \
            compare_observed_barcode_with_known_barcodes(
                iii, known_viral_barcodes_dict,
                edit_distance_threthold)

        if matched_known_viral_barcode:
            matched_known_viral_barcodes_occurrence += 1
            known_viral_barcodes_dict[matched_known_viral_barcode] += 1
            # print(matched_known_viral_barcode)
        else:
            noise_viral_barcodes.append(iii)
            noise_viral_barcodes_occurrence += 1

            if iii in noise_viral_barcodes_dict:
                noise_viral_barcodes_dict[iii] += 1
            else:
                noise_viral_barcodes_dict[iii] = int()
                noise_viral_barcodes_dict[iii] += 1

    # print(noise_viral_barcodes)
    # print(noise_viral_barcodes_dict)

    noise_viral_barcodes_container = dict()

    if noise_viral_barcodes:
        uniq_noise_viral_barcodes = set(noise_viral_barcodes)

        for k in uniq_noise_viral_barcodes:
            selected_noise_barcode = k
            selected_noise_barcode_read_count = noise_viral_barcodes_dict[k]
            del noise_viral_barcodes_dict[k]

            indicator = int()
            for kk in noise_viral_barcodes_dict:
                edit_distance = calc_edit_distance(k, kk)

                if edit_distance <= edit_distance_threthold:
                    noise_viral_barcodes_dict[kk] += selected_noise_barcode_read_count
                    indicator = 1
                    noise_viral_barcodes_container[kk] = noise_viral_barcodes_dict[kk]
            if not indicator:
                noise_viral_barcodes_container[selected_noise_barcode] = selected_noise_barcode_read_count

    # print(noise_viral_barcodes_dict)
    # print(noise_viral_barcodes_container)

    matched_known_viral_barcodes_num = len([i for i in known_viral_barcodes_dict
                                            if known_viral_barcodes_dict[i]])
    noise_viral_barcodes_num = len(set(noise_viral_barcodes))
    matched_known_barcodes_results = list()

    sorted_matched_barcodes = \
        sorted(known_viral_barcodes_dict.items(),
               key=lambda
               known_viral_barcodes_dict: known_viral_barcodes_dict[1],
               reverse = True)
    # print(sorted_matched_barcodes)

    if matched_known_viral_barcodes_num:
        for jj in sorted_matched_barcodes:
            jj = jj[0]

            if known_viral_barcodes_dict[jj]:
                matched_known_barcodes_results.append(
                    [jj, known_viral_barcodes_dict[jj]])
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
          matched_known_viral_barcodes_num,
          matched_known_viral_barcodes_occurrence,
          len(noise_viral_barcodes_container),
          noise_viral_barcodes_occurrence, sep = '\t')


    '''
    print(selected_STAMPs[j].rstrip(),
          matched_known_barcodes_output,
          matched_known_viral_barcodes_num,
          matched_known_viral_barcodes_occurrence,
          noise_viral_barcodes_num,
          noise_viral_barcodes_occurrence, sep = '\t')
    '''
