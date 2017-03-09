#!/usr/bin/env python3
import sys

refFlat = '/project/GCRB/Hon_lab/s166631/01.data/reference/hg19_ucsc_lentiGuide-MS2-puro-barcode-empty_STAR_2.5.0a/hg19_chr1-Y_chrM_lentiGuide-MS2-puro-barcode-empty.refFlat'

gene_symbols = {i.split('\t')[1]: i.split('\t')[0] for i in open(refFlat, 'r')}

with open(sys.argv[1], 'r') as f:
    for line in f:
        if line.startswith('#'):
            print('#Genesymbol', line.rstrip()[1:], sep = '\t')
        else:
            i = line.rstrip().split('\t')
            print(gene_symbols[i[0]], line.rstrip(), sep = '\t')
