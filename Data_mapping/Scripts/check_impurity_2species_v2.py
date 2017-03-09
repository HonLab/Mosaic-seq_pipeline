#!/usr/bin/env python3
import os
import sys
import pysam
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def check_impurity(bam_file, minimal_mapping_quality):
    """to be filled"""
    STAMPs = dict()

    f = pysam.AlignmentFile(bam_file)

    for read in f.fetch(until_eof = True):
        if read.mapping_quality >= minimal_mapping_quality:
            mapped_ref_name = f.get_reference_name(read.reference_id)

            for i in read.tags:
                if i[0] == 'XC':
                    cell_barcode = i[1]
                    cell_barcode = cell_barcode.encode(encoding='UTF-8')
                    # print(cell_barcode)

            if cell_barcode in STAMPs:
                if mapped_ref_name.startswith('hg19'):
                    STAMPs[cell_barcode][0] += 1
                elif mapped_ref_name.startswith('mm9'):
                        STAMPs[cell_barcode][1] += 1
            else:
                STAMPs[cell_barcode] = [int(), int()]
                if mapped_ref_name.startswith('hg19'):
                    STAMPs[cell_barcode][0] += 1
                elif mapped_ref_name.startswith('mm9'):
                    STAMPs[cell_barcode][1] += 1
    f.close()


    STAMPs_less = {i: sum(STAMPs[i]) for i in STAMPs}
    STAMPs_sorted = sorted(STAMPs_less.items(),
        key=lambda STAMPs_less: STAMPs_less[1], reverse = True)

    total_reads = sum([sum(STAMPs[i]) for i in STAMPs])

    f_output = open(bam_file + '_mapQ' + str(minimal_mapping_quality) +
        '_impurity_info.txt','w')

    count = int()
    for i in STAMPs_sorted:
        count += 1

        j = i[0]
        STAMPs_info = '\t'.join([str(count), j.decode(), str(STAMPs[j][0]),
            str(STAMPs[j][1]), str(sum(STAMPs[j])), str(total_reads)])
        f_output.write(STAMPs_info + '\n')

    f_output.close()


    scatter_plot_x = [STAMPs[i][0] for i in STAMPs]
    scatter_plot_y = [STAMPs[i][1] for i in STAMPs]
    # colors = ["blue" for i in STAMPs]

    with PdfPages(bam_file + '_mapQ' + str(minimal_mapping_quality) +
        '_impurity_scatter_plot.pdf') as pdf:

        plt.scatter(scatter_plot_x, scatter_plot_y,
            c='steelblue', edgecolors = None) # edgecolors = 'face'
        plt.axis([0, max(scatter_plot_x), 0, max(scatter_plot_y)])
        plt.xlabel('No. of reads mapped on human reference')
        plt.ylabel('No. of reads mapped on mouse reference')
        # plt.title()
        pdf.savefig()
        plt.close()

    saturation_plot_y_before_cumulation = [sum(STAMPs[i[0]]) / total_reads for i in STAMPs_sorted]
    saturation_plot_y = list()
    saturation_plot_y_individual = float()

    for i in saturation_plot_y_before_cumulation:
        saturation_plot_y_individual += i
        saturation_plot_y.append(saturation_plot_y_individual)

    saturation_plot_x = range(len(saturation_plot_y))

    with PdfPages(bam_file + '_mapQ' + str(minimal_mapping_quality) +
        '_impurity_saturation_plot.pdf') as pdf:

        plt.plot(saturation_plot_x, saturation_plot_y,
            c='steelblue')
        plt.axis([0, 1000, 0, 1])
        plt.xlabel('STAMPs (ordered largest to smallest)')
        plt.ylabel('Cumulative fraction of reads')
        # plt.title()
        pdf.savefig()
        plt.close()


def is_bam_sorted_by_coordinate(bam_file):
    """to be filled"""
    f = pysam.AlignmentFile(bam_file)
    if f.header['HD']['SO'] == 'coordinate':
        return True
    f.close()

def index_bam(bam_file):
    """to be filled"""
    f = pysam.AlignmentFile(bam_file)
    if not f.has_index():
        os.system('samtools index ' + bam_file)
    f.close()


mapping_quality_threshold = 10

if __name__ == "__main__":
    if is_bam_sorted_by_coordinate(sys.argv[1]):
        index_bam(sys.argv[1])
        check_impurity(sys.argv[1], mapping_quality_threshold)
