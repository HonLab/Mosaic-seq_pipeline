#!/usr/bin/env python3
import sys
import re
# import statistics


# gene_length = {i.split('\t')[0]: int(i.rstrip().split('\t')[5])
#     for i in open(sys.argv[1], 'r') if not i.startswith(('#', 'Geneid'))}

# gene_list = [i.split('\t')[0] for i in open(sys.argv[1], 'r')
#     if not i.startswith(('#', 'Geneid'))]

gene_length = dict()
gene_list = list()

with open(sys.argv[1], 'r') as f:
    for line in f:
        if line.startswith('Geneid'):
            i = line.rstrip().split('\t')
            barcode_list = [re.sub('\.bam$', '', ii) for ii in i[6:]]
            barcodes = {i: [] for i in barcode_list}
            # break
        elif not line.startswith(('#', 'Geneid')):
            i = line.rstrip().split('\t')
            # print(i[6: ])

            gene_length[i[0]] = int(i[5])
            gene_list.append(i[0])

            for j, ii in enumerate(i[6: ]):
                barcodes[barcode_list[j]].append(int(ii))

total_reads_of_each_barcode = [sum(barcodes[i]) for i in barcode_list]
total_gene_length_of_each_barcode = sum(gene_length.values())

f_output = open(sys.argv[1].replace('.txt', '') + '_TPM.txt', 'w')
output_line = '#Geneid' + '\t' + '\t'.join(barcode_list)
f_output.write(output_line + '\n')
# print('#'+ '\t'.join(barcode_list))

total_RPKM_of_each_barcode = dict()
gene_length_list = [gene_length[i] for i in gene_list]

for j, i in enumerate(barcode_list):
    total_RPKM_of_each_barcode[i] = \
        sum([ji[1] * (10 ** 9) / \
        (total_reads_of_each_barcode[j] * gene_length_list[ji[0]])
        for ji in enumerate(barcodes[i])])

for j, i in enumerate(gene_list):
    expr_of_each_transcript_of_each_barcode  = list()

    for jj, ii in enumerate(barcode_list):
        RPKM = barcodes[ii][j] * (10 ** 9) / (total_reads_of_each_barcode[jj] *
               gene_length[i])
        # RPKM = barcodes[ii][j] / \
        #        (gene_length[i] / \
        #         1000 * total_reads_of_each_barcode[jj] / 1000000)

        TPM = RPKM / total_RPKM_of_each_barcode[ii] * 10 ** 6
        expr_of_each_transcript_of_each_barcode.append(TPM)

    output_line = [i] + [str(i)
                         for i in expr_of_each_transcript_of_each_barcode]
    f_output.write('\t'.join(output_line) + '\n')
    # print(i,
    #       '\t'.join(str(i) for i in expr_of_each_transcript_of_each_barcode),
    #       sep = '\t')

f_output.close()
