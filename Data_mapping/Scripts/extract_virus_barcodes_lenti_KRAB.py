#!/usr/bin/env python3
import sys
import pysam
import os

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

def retrieve_aligned_reads(bam_file, mapping_quality_threshold,
                           target_region_reference,
                           target_region_start,
                           target_region_end):
    """to be filled"""
    f = pysam.AlignmentFile(bam_file)
    aligned_reads = f.fetch(target_region_reference,
                    target_region_start, target_region_end)

    retrieved_reads = list()

    for read in aligned_reads:
        if read.reference_start <= target_region_start and \
           read.reference_end >= target_region_end and \
           read.mapping_quality >= mapping_quality_threshold:

            read_left = target_region_start - read.reference_start + \
                        read.query_alignment_start
            read_right = read.query_alignment_end - \
                         (read.reference_end - target_region_end)

            retrieved_reads.append(
            [read.qname, read.query_alignment_sequence[read_left: read_right],
             read.query_alignment_qualities[read_left: read_right]])

    return retrieved_reads


#target_region_reference = 'lentiGuide-MS2-puro-barcode-empty'
target_region_reference = 'lenti_KRAB'
# 0 based, [)
target_region_start = 5172 - 1
target_region_end = 5183

mapping_quality_threshold = 10

if __name__ == "__main__":
    if is_bam_sorted_by_coordinate(sys.argv[1]):
        index_bam(sys.argv[1])
        aligned_reads = retrieve_aligned_reads(sys.argv[1],
                                               mapping_quality_threshold,
                                               target_region_reference,
                                               target_region_start,
                                               target_region_end)

        for i in aligned_reads:
            i[2] = ','.join(str(ii) for ii in i[2])

        for i in aligned_reads:
            print('\t'.join(i))
