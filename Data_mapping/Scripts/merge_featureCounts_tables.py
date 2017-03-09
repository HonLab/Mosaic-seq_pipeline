#!/usr/bin/env python3
import sys

transcript_exprs = {i.split('\t')[0]: []
                    for i in open(sys.argv[1], 'r')
                    if not i.startswith(('#', 'Geneid'))}
transcripts = [i.split('\t')[0]
               for i in open(sys.argv[1], 'r')
               if not i.startswith(('#', 'Geneid'))]
# num_transcripts = len(transcript_exprs)
transcripts_info = [i.split('\t')[: 6]
                      for i in open(sys.argv[1], 'r')
                      if not i.startswith(('#', 'Geneid'))]

transcript_exprs_header = list()

for expr_table in sys.argv[1:]:
    with open(expr_table, 'r') as f:
        num_line = int()
        for line in f:
            if line.startswith('Geneid'):
                i = line.rstrip().split('\t')[6:]
                transcript_exprs_header.append(i)

            elif not line.startswith(('#', 'Geneid')):
                i = line.rstrip().split('\t')[6:]
                transcript_exprs[transcripts[num_line]].append(i)
                num_line += 1

merged_header = list()
for i in transcript_exprs_header:
    merged_header += i
print('Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length',
      '\t'.join(merged_header), sep='\t')

for index, transcript in enumerate(transcripts):
    merged_exprs = []

    for i in transcript_exprs[transcript]:
        merged_exprs += i
    print('\t'.join(transcripts_info[index]), '\t'.join(merged_exprs), sep='\t')
