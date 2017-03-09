#!/bin/tcsh

#SBATCH -J DropSeq_mapping_1
#SBATCH -p super
#SBATCH -N 1
#SBATCH -t 0-24:00:00
#SBATCH --mail-user=shiqi.xie@utsouthwestern.edu
#SBATCH --mail-type=end
#SBATCH -o preprocess.out.%j
#SBATCH -e preprocess.err.%j

# scripts and directories
set SHARED_DIR        = /project/GCRB/Hon_lab/shared
set SCRIPT_DIR        = /project/GCRB/Hon_lab/shared/pipelines/Mosaic-Seq/2016-10-22/scripts
set DROPSEQ_TOOLS_DIR = $SHARED_DIR/software/Drop-seq_tools-1.0
set STAR              = $SHARED_DIR/software/STAR-STAR_2.4.2a/source/STAR
set DROPSEQ_PIPE      = $DROPSEQ_TOOLS_DIR/Drop-seq_alignment.with_threads.sh
set DetectBeadSynthesisErrors = /project/GCRB/Hon_lab/shared/software/Drop-seq_tools-1.11/DetectBeadSynthesisErrors 
set FEATURE_COUNTS    = /project/GCRB/Hon_lab/shared/software/subread-1.4.6-p4-Linux-x86_64/bin/featureCounts

# python scripts
set evaluate_STAMPs_by_cumulative_reads = $SCRIPT_DIR/evaluate_STAMPs_by_cumulative_reads2.py
set extract_virus_barcodes              = $SCRIPT_DIR/extract_virus_barcodes.py
set mark_STAMPs_with_viral_barcodes     = $SCRIPT_DIR/mark_STAMPs_with_viral_barcodes_custom.py
set calc_CPM                            = $SCRIPT_DIR/calc_CPM.py
set calc_TPM                            = $SCRIPT_DIR/calc_TPM.py
set calc_RPKM                           = $SCRIPT_DIR/calc_RPKM.py
set extract_noisy_reads_dual            = $SCRIPT_DIR/extract_noisy_reads_dual.py
set perform_wilcoxon_test               = $SCRIPT_DIR/perform_wilcoxon_test.R

# genome parameters
set STAR_DIRECTORY    = /project/GCRB/Hon_lab/s166631/01.data/reference/hg19_ucsc_lentiGuide-MS2-puro-barcode-empty_STAR_2.5.0a
set GENOME            = $STAR_DIRECTORY/hg19_chr1-Y_chrM_lentiGuide-MS2-puro-barcode-empty.fasta
set GENOME_SIZE_FILE  = $STAR_DIRECTORY/hg19_chr1-Y_chrM_lentiGuide-MS2-puro-barcode-empty.fasta.fai
set GENOME_GTF        = /project/GCRB/Hon_lab/shared/pipelines/Drop-Seq/2016-01-29/db/hg19_chr1-Y_chrM_lentiGuide-MS2-puro-barcode-empty.gtf

# pipeline parameters
set SEQ_DIR                     = /project/GCRB/Hon_lab/shared/data/sequencing_data/2017/2017-01-23-Mcdermott-NextSeq.bcl2/Drop-Seq
set BARCODE_LIST                = $PWD/known_barcodes.list
set MIN_READS_FOR_BIG_BARCODES  = 1000
set MIN_READS_FOR_SUB_BARCODES  = 100
set EST_NUM_CELL_BCS            = 2000   # about 5 * (# expected cells)
set PRIMER_SEQ                  = AAGCAGTGGTATCAACGCAGAGTAC
#set NUM_STAMPS_SLOPE_CUTOFF     = 0.001
set HIGH_COVERAGE_CELL_FRACTION = 0.05

module load picard/1.117
module load java/oracle/jdk1.7.0_51
module load samtools/intel/1.1
#module add python/3.4.x-anaconda
setenv PATH /home2/s160875/.conda/envs/py35/bin:$PATH

foreach lib_plate (`cat run1.txt | grep -v \#`)
  set lib  = `echo $lib_plate | cut -f1 -d,`
  set plate = `echo $lib_plate | cut -f2 -d,`
  echo "$lib, $plate"

  set out = $lib.dropseq_pipe
  rm -rf $out
  mkdir $out $out/temp $out/barcode_counts $out/raw $out/reports $out/genic_counts $out/bam $out/STAMPs $out/genic_counts/500cells $out/genic_counts/STAMPs $out/genic_counts/STAMPs.high_cov $out/bam/depth

  foreach read (1 2)
    zcat $SEQ_DIR/$lib*L001*R$read*001.fastq.gz\
         $SEQ_DIR/$lib*L002*R$read*001.fastq.gz\
         $SEQ_DIR/$lib*L003*R$read*001.fastq.gz\
         $SEQ_DIR/$lib*L004*R$read*001.fastq.gz\
      | gzip\
      > $out/raw/$lib.combined.R$read.gz
  end

  # make BAM file
  java -Xmx120g -jar $PICARD_DIR/FastqToSam.jar\
    SAMPLE_NAME=$lib\
    FASTQ=$out/raw/$lib.combined.R1.gz\
    FASTQ2=$out/raw/$lib.combined.R2.gz\
    OUTPUT=$out/bam/$lib.unaligned.bam\
    SORT_ORDER=queryname

  # align
  $DROPSEQ_PIPE\
    -g $STAR_DIRECTORY\
    -r $GENOME\
    -d $DROPSEQ_TOOLS_DIR\
    -o $out\
    -t $out/temp\
    -s $STAR\
    -p\
    -T $SLURM_CPUS_ON_NODE\
    $out/bam/$lib.unaligned.bam

  # move bam output to bam folder
  mv $out/star_gene_exon_tagged.ba? $out/bam/
  mv $out/adapter_trimming_report.txt                $out/reports/
  mv $out/unaligned_tagged_Cellular.bam_summary.txt  $out/reports/
  mv $out/unaligned_tagged_Molecular.bam_summary.txt $out/reports/
  mv $out/polyA_trimming_report.txt $out/reports/

  # get histogram, original
  set out1 = $out/bam/star_gene_exon_tagged.bam
  $DROPSEQ_TOOLS_DIR/BAMTagHistogram -m 4g\
    I=$out1\
    O=$out/barcode_counts/$lib.barcode_counts\
    TAG=XC

  # remove bad quality mapping
  set out2 = $out/bam/star_gene_exon_tagged.q10.bam
  samtools view -b -q10 -o $out2 $out1

  # get histogram, after quality filtering
  $DROPSEQ_TOOLS_DIR/BAMTagHistogram -m 4g\
    I=$out2\
    O=$out/barcode_counts/$lib.barcode_counts.q10\
    TAG=XC

  # fix bead synthesis error
  set out3 = $out/bam/star_gene_exon_tagged.q10.cleaned.bam
  $DetectBeadSynthesisErrors -m 4g\
    I=$out2\
    O=$out3\
    OUTPUT_STATS=$out/reports/$lib.synthesis_stats.txt\
    SUMMARY=$out/reports/$lib.synthesis_stats.summary.txt\
    NUM_BARCODES=$EST_NUM_CELL_BCS\
    PRIMER_SEQUENCE=$PRIMER_SEQ

  # get histogram, after bead synthesis fix
  $DROPSEQ_TOOLS_DIR/BAMTagHistogram -m 4g\
    I=$out3\
    O=$out/barcode_counts/$lib.barcode_counts.q10.cleaned\
    TAG=XC

  # merge cell barcodes
  set out4 = $out/bam/star_gene_exon_tagged.q10.cleaned.cell_bc_merged.bam
  samtools view -H $out1 > $out1.header
  $SCRIPT_DIR/merge_cell_bcs3.pl $out3\
                       $MIN_READS_FOR_BIG_BARCODES\
                       $MIN_READS_FOR_SUB_BARCODES\
                       $out1.header\
    | samtools view -bS -t $GENOME_SIZE_FILE -o $out4 -

  # get histogram, after cell bc merging
  $DROPSEQ_TOOLS_DIR/BAMTagHistogram -m 4g\
    I=$out4\
    O=$out/barcode_counts/$lib.barcode_counts.q10.cleaned.cell_bc_merged\
    TAG=XC

  # remove PCR duplicates, edit distance = 1
  set out5 = $out/bam/star_gene_exon_tagged.q10.cleaned.cell_bc_merged.nodup.bam
  $SCRIPT_DIR/rmdup3.pl $out4\
                        $out1.header\
    | samtools view -bS -t $GENOME_SIZE_FILE -o $out5 -

  # get histogram, after PCR duplicate removal
  $DROPSEQ_TOOLS_DIR/BAMTagHistogram -m 4g\
    I=$out5\
    O=$out/barcode_counts/$lib.barcode_counts.q10.cleaned.cell_bc_merged.nodup\
    TAG=XC

  echo "# reads, $out1\t"`samtools view $out1 | wc -l` >  $out/reports/$lib.processing.summary
  echo "# reads, $out2\t"`samtools view $out2 | wc -l` >> $out/reports/$lib.processing.summary
  echo "# reads, $out3\t"`samtools view $out3 | wc -l` >> $out/reports/$lib.processing.summary
  echo "# reads, $out4\t"`samtools view $out4 | wc -l` >> $out/reports/$lib.processing.summary
  echo "# reads, $out5\t"`samtools view $out5 | wc -l` >> $out/reports/$lib.processing.summary

  # make 1 BAM file per cell
  set split_out = $out/bam/split.bam
  rm -rf $split_out
  mkdir $split_out

  set bc_counts = $out/barcode_counts/$lib.barcode_counts.q10.cleaned.cell_bc_merged.nodup
  $SCRIPT_DIR/split_bam_per_cell.pl $out5\
                                    $out1.header\
                                    $bc_counts\
                                    $split_out\
                                    $EST_NUM_CELL_BCS

  cut -f2 $bc_counts\
    | head -$EST_NUM_CELL_BCS\
    | perl -ne 'chomp; print("$_.bam ");'\
    > $out/bam/bam_list

  # count reads, 
  cd $split_out
  $FEATURE_COUNTS -O\
                  -T $SLURM_CPUS_ON_NODE\
                  -a $GENOME_GTF\
                  -o ../../genic_counts/500cells/$lib.q10.cleaned.cell_bc_merged.nodup.exons\
                  `cat ../bam_list`

  foreach bam (`ls *.bam`)
    $extract_virus_barcodes $bam > `echo $bam`_barcodes.txt;
    samtools view $bam\
      | wc -l\
      > $bam.counts
  end
  cd ../../..

  # get viral barcode information
  cd $split_out
  foreach bam (`ls *.bam | sed -e 's/.bam//'`)
    $extract_virus_barcodes $bam.bam > `echo $bam`_barcodes.txt;
  end
  cd ../../..

  # create list of STAMPs
 set BARCODES_FILE = STAMPs_info_mapQ10_top500_virus_barcodes_editdist1.txt
  $evaluate_STAMPs_by_cumulative_reads $out5\
                                       $out/STAMPs
  cd $out
  $mark_STAMPs_with_viral_barcodes STAMPs/STAMPs_info_mapQ10.txt\
                                   $BARCODE_LIST\
                                   $EST_NUM_CELL_BCS\
    > STAMPs/$BARCODES_FILE
  cd ..

  # get CPM, TPM, RPKM
  $calc_CPM  $out/genic_counts/500cells/$lib.q10.cleaned.cell_bc_merged.nodup.exons
  $calc_TPM  $out/genic_counts/500cells/$lib.q10.cleaned.cell_bc_merged.nodup.exons
  $calc_RPKM $out/genic_counts/500cells/$lib.q10.cleaned.cell_bc_merged.nodup.exons

  # get sequencing stats
  set seqx4 = `zcat $out/raw/$lib.combined.R1.gz | wc -l`
  echo $seqx4\
    | perl -ne 'print($_ / 4);'\
    > $out/reports/$lib.num_reads_seq
  tail -1 $out/reports/$lib.processing.summary\
    | cut -f2\
    > $out/reports/$lib.num_nodup
  head -2 $out/reports/$lib.processing.summary\
    | tail -1\
    | cut -f2\
    > $out/reports/$lib.num_mapped

  # make a saturation plot of STAMPs
  cat $SCRIPT_DIR/myplot_saturation.template_flex.m\
    | sed -e "s/CELL_BARCODE_COUNT_FILE/..\/STAMPs\/$BARCODES_FILE/"\
    | sed -e "s/NODUP_SIZE_FILE/$lib.num_nodup/"\
    | sed -e "s/NAME/$lib/"\
    | sed -e "s/OUT_FILE/$lib.saturation.pdf/"\
    | sed -e "s/EST_NUM_CELL_BCS/$EST_NUM_CELL_BCS/"\
    > $out/reports/saturation.m
  cd $out/reports
  matlab  < saturation.m > saturation.m.out
  cd ../..

  # get barcodes belonging to final STAMPs
  cat $out/STAMPs/$BARCODES_FILE\
    | $SCRIPT_DIR/get_STAMPs_flex.pl $EST_NUM_CELL_BCS\
                                     $out/reports/$lib.num_nodup\
    > $out/STAMPs/$BARCODES_FILE.STAMPs

  # get read counts corresponding to STAMPs
  $SCRIPT_DIR/get_STAMP_counts.pl $out/genic_counts/500cells/$lib.q10.cleaned.cell_bc_merged.nodup.exons\
                                  $out/STAMPs/$BARCODES_FILE.STAMPs\
    > $out/genic_counts/STAMPs/$lib.q10.cleaned.cell_bc_merged.nodup.exons.STAMPs

  # get CPM, TPM, RPKM of STAMPs
  $calc_CPM  $out/genic_counts/STAMPs/$lib.q10.cleaned.cell_bc_merged.nodup.exons.STAMPs
  $calc_TPM  $out/genic_counts/STAMPs/$lib.q10.cleaned.cell_bc_merged.nodup.exons.STAMPs
  $calc_RPKM $out/genic_counts/STAMPs/$lib.q10.cleaned.cell_bc_merged.nodup.exons.STAMPs

  cut -f7- $out/genic_counts/STAMPs/$lib.q10.cleaned.cell_bc_merged.nodup.exons.STAMPs\
    > $out/genic_counts/STAMPs/$lib.q10.cleaned.cell_bc_merged.nodup.exons.STAMPs.counts_only

  cat $out/genic_counts/STAMPs/$lib.q10.cleaned.cell_bc_merged.nodup.exons.STAMPs.counts_only\
    | $SCRIPT_DIR/get_stats.pl $lib\
                               $out/reports/$lib.num_reads_seq\
                               $out/reports/$lib.num_nodup\
    > $out/reports/$lib.stats

  # make CoV plot
  set TEMP_FILE = temp
  cat $SCRIPT_DIR/cov_plot_template.m\
    | sed -e "s/MATRIX_FILE/$TEMP_FILE/"\
    | sed -e "s/NAME/$lib/"\
    | sed -e "s/OUT_FILE/$lib.CoV_vs_mean.pdf/"\
    > $out/reports/CoV.STAMPs.m
  cd $out/reports
  cp ../genic_counts/STAMPs/$lib.q10.cleaned.cell_bc_merged.nodup.exons.STAMPs.counts_only $TEMP_FILE
  matlab < CoV.STAMPs.m > CoV.STAMPs.m.out
  rm $TEMP_FILE
  cd ../..

  # get high coverage cells
  cat $out/genic_counts/STAMPs/$lib.q10.cleaned.cell_bc_merged.nodup.exons.STAMPs\
    | $SCRIPT_DIR/get_cells_high_coverage.pl $HIGH_COVERAGE_CELL_FRACTION\
    > $out/genic_counts/STAMPs.high_cov/$lib.q10.cleaned.cell_bc_merged.nodup.exons.STAMPs.high_cov

  # get CPM, TPM, RPKM of high coverage STAMPs
  $calc_CPM  $out/genic_counts/STAMPs.high_cov/$lib.q10.cleaned.cell_bc_merged.nodup.exons.STAMPs.high_cov
  $calc_TPM  $out/genic_counts/STAMPs.high_cov/$lib.q10.cleaned.cell_bc_merged.nodup.exons.STAMPs.high_cov
  $calc_RPKM $out/genic_counts/STAMPs.high_cov/$lib.q10.cleaned.cell_bc_merged.nodup.exons.STAMPs.high_cov

  cut -f7- $out/genic_counts/STAMPs.high_cov/$lib.q10.cleaned.cell_bc_merged.nodup.exons.STAMPs.high_cov\
    > $out/genic_counts/STAMPs.high_cov/$lib.q10.cleaned.cell_bc_merged.nodup.exons.STAMPs.high_cov.counts_only

  # make CoV plot
  set TEMP_FILE = temp
  cat $SCRIPT_DIR/cov_plot_template.m\
    | sed -e "s/MATRIX_FILE/$TEMP_FILE/"\
    | sed -e "s/NAME/$lib/"\
    | sed -e "s/OUT_FILE/$lib.CoV_vs_mean.high_cov.pdf/"\
    > $out/reports/CoV.STAMPs.high_cov.m
  cd $out/reports
  cp ../genic_counts/STAMPs.high_cov/$lib.q10.cleaned.cell_bc_merged.nodup.exons.STAMPs.high_cov.counts_only $TEMP_FILE
  matlab < CoV.STAMPs.high_cov.m > CoV.STAMPs.high_cov.m.out
  rm $TEMP_FILE
  cd ../..

  # get read counts
  foreach kind (500cells STAMPs STAMPs.high_cov)
    $SCRIPT_DIR/count_reads.pl $out/genic_counts/$kind\
                               $out/bam/split.bam\
                               $out/bam/depth/$kind
  end

  # get stats
  cat $out/genic_counts/STAMPs/$lib.q10.cleaned.cell_bc_merged.nodup.exons.STAMPs.counts_only\
    | $SCRIPT_DIR/get_stats2.pl $lib\
                                $out/reports/$lib.num_reads_seq\
                                $out/reports/$lib.num_mapped\
                                $out/reports/$lib.num_nodup\
                                $out/bam/depth\
    > $out/reports/$lib.stats2


  ##################################################
  #
  # DECIDE WHICH CELLS HAVE WHICH BARCODES
  #
  ##################################################
  echo "\n\nDECIDING WHICH CELLS HAVE WHICH BARCODES\n\n"
  mkdir $out/split_by_barcode $out/split_by_barcode/region_to_stats
 cd $out/split_by_barcode

  foreach region_sgRNA_viral_bc (`cat ../../plate$plate.txt`)
    set region_sgRNA = `echo $region_sgRNA_viral_bc | cut -f1 -d\;`
    set region       = `echo $region_sgRNA          | perl -ne 's/.sgRNA.+//; print;'`
    set viral_bc     = `echo $region_sgRNA_viral_bc | cut -f2 -d\;`

    if (-e $region_sgRNA) then
    else
      mkdir region_to_stats/$region_sgRNA
    endif

    echo $lib, $plate, $region, $viral_bc

    cat ../STAMPs/STAMPs_info_mapQ10_top500_virus_barcodes_editdist1.txt.STAMPs\
      | grep -vP "\tNA\t"\
      | grep $viral_bc\
      | cut -f2,8,9\
      > region_to_stats/$region_sgRNA/$lib.$plate.$viral_bc.raw
  end

  #
  # use model to determine the number of states
  #
  cat region_to_stats/*/*.raw\
    | sort\
    | uniq\
    > raw.all

  cut -f3 raw.all\
    | sed -e 's/,/\t/g'\
    | $SCRIPT_DIR/pad.pl\
    > raw.all.pad

#  cat raw.all\
#    | perl -ne 'chomp; ($cell, $bcs, $counts) = split(/\t/, $_); @bcs = split(/,/, $bcs); print(join("\t", $cell, $bcs, $counts, scalar(@bcs)) . "\n");'\
#    > raw.all.pad.states.txt.int.merged

  # set up Matlab scripts
  set num_distinct_barcodes_on_plate = `cat ../../plate$plate.txt | wc -l`
  cp $SCRIPT_DIR/null_probability.m .
  cat $SCRIPT_DIR/my_model_template.m\
    | sed -e "s/NUM_DISTINCT_BARCODES_ON_PLATE/$num_distinct_barcodes_on_plate/"\
    > my_model.m

  # execute matlab
  matlab -nodesktop -nosplash < my_model.m > mymodel.m.out

  cat raw.all.pad.states.txt\
    | perl -ne 'split; $_[0] = int($_[0]); print("$_[0]\n");'\
    > raw.all.pad.states.txt.int

  paste raw.all\
        raw.all.pad.states.txt.int\
    > raw.all.pad.states.txt.int.merged

  #
  # post-processing
  #

  # combine file for each region
  foreach region_sgRNA_viral_bc (`cat ../../plate$plate.txt`)
    set region_sgRNA = `echo $region_sgRNA_viral_bc | cut -f1 -d\;`
    set region       = `echo $region_sgRNA          | perl -ne 's/.sgRNA.+//; print;'`
    set viral_bc     = `echo $region_sgRNA_viral_bc | cut -f2 -d\;`

    echo $lib, $plate, $region, $viral_bc

    cat region_to_stats/$region_sgRNA/$lib.$plate.$viral_bc.raw\
      | $SCRIPT_DIR/filter_hits.pl $viral_bc\
                                   raw.all.pad.states.txt.int.merged\
      > region_to_stats/$region_sgRNA/$lib.$plate.$viral_bc.raw.hits
  end

  # get order of cells
  head -1 ../genic_counts/STAMPs/*.q10.cleaned.cell_bc_merged.nodup.exons.STAMPs\
    | cut -f7-\
    > cells.order

  # combine if multiple samples within each folder
  foreach region_sgRNA (`ls region_to_stats`)
    cat region_to_stats/$region_sgRNA/*hits\
      | sort -k4,4n\
      > region_to_stats/$region_sgRNA/combined

    $SCRIPT_DIR/get_labels.pl region_to_stats/$region_sgRNA/combined\
                              cells.order\
      > region_to_stats/$region_sgRNA/combined.labels

    grep -v Geneid region_to_stats/$region_sgRNA/combined.labels\
      | cut -f2\
      > region_to_stats/$region_sgRNA/combined.labels.depth

    grep -v Geneid region_to_stats/$region_sgRNA/combined.labels\
      | cut -f3\
      > region_to_stats/$region_sgRNA/combined.labels.percent

    grep -v Geneid region_to_stats/$region_sgRNA/combined.labels\
      | cut -f4\
      > region_to_stats/$region_sgRNA/combined.labels.num_states
  end

  # combine all regions together
  ls region_to_stats\
    | grep sgRNA\
    > sgRNAs.order
  foreach kind (depth percent num_states)
    echo $kind
    paste region_to_stats/*.sgRNA*/combined.labels.$kind\
      > region_to_stats/merged.$kind
  end


  # SPLIT MATRIX BY SGRNA
  mkdir ../genic_counts/split_by_sgRNA

  foreach region_sgRNA_viral_bc (`cat ../../plate$plate.txt`)
    set region_sgRNA = `echo $region_sgRNA_viral_bc | cut -f1 -d\;`
    set region       = `echo $region_sgRNA          | perl -ne 's/.sgRNA.+//; print;'`
    set viral_bc     = `echo $region_sgRNA_viral_bc | cut -f2 -d\;`

    cut -f1 region_to_stats/$region_sgRNA/combined\
      | $SCRIPT_DIR/filter_matrix_by_cell.pl ../genic_counts/STAMPs/$lib.q10.cleaned.cell_bc_merged.nodup.exons.STAMPs\
      > ../genic_counts/split_by_sgRNA/$region_sgRNA
  end

  #########################
  # ensure that control cells are not found in non-control cells
  #########################
  cd ../genic_counts
  mkdir split_by_sgRNA.control_uniq

  # get list of control cells
  paste split_by_sgRNA/Control.sgRNA*\
    | head -1\
    | sed -e 's/\t/\n/g'\
    | grep bam\
    | sort\
    | uniq\
    > split_by_sgRNA.control_uniq/cell_barcodes.control

  # get list of non-control cells
  foreach region_sgRNA_viral_bc (`cat ../../plate$plate.txt | grep -v Control`)
    set region_sgRNA = `echo $region_sgRNA_viral_bc | cut -f1 -d\;`
    set region       = `echo $region_sgRNA          | perl -ne 's/.sgRNA.+//; print;'`
    set viral_bc     = `echo $region_sgRNA_viral_bc | cut -f2 -d\;`

    head -1 split_by_sgRNA/$region_sgRNA\
      | sed -e 's/\t/\n/g'\
      | grep bam\
      | sort\
      | uniq\
      >> split_by_sgRNA.control_uniq/cell_barcodes.non_control
  end

  cat split_by_sgRNA.control_uniq/cell_barcodes.non_control\
    | sort\
    | uniq\
    > split_by_sgRNA.control_uniq/cell_barcodes.non_control.uniq

  # get list of cell barcodes that are in common between control and non-control
  cat split_by_sgRNA.control_uniq/cell_barcodes.control\
      split_by_sgRNA.control_uniq/cell_barcodes.non_control.uniq\
    | sort\
    | uniq -c\
    | perl -ne 'split; if($_[0] != 1){print("$_[1]\n");}'\
    > split_by_sgRNA.control_uniq/cell_barcodes.overlap

  # print matrix without overlapped cells
  foreach region_sgRNA_viral_bc (`cat ../../plate$plate.txt`)
    set region_sgRNA = `echo $region_sgRNA_viral_bc | cut -f1 -d\;`
    set region       = `echo $region_sgRNA          | perl -ne 's/.sgRNA.+//; print;'`
    set viral_bc     = `echo $region_sgRNA_viral_bc | cut -f2 -d\;`

    cat split_by_sgRNA.control_uniq/cell_barcodes.overlap\
      | $SCRIPT_DIR/remove_overlapped_cells.pl split_by_sgRNA/$region_sgRNA\
      > split_by_sgRNA.control_uniq/$region_sgRNA

    cut -f1 split_by_sgRNA.control_uniq/$region_sgRNA\
      > split_by_sgRNA.control_uniq/ids

    # get only values, with no gene column
    cat split_by_sgRNA.control_uniq/$region_sgRNA\
      | perl -ne '@a = split; shift(@a); print(join("\t", @a) . "\n");'\
      > split_by_sgRNA.control_uniq/$region_sgRNA.vals
  end

  # cd out
  cd ../..
end
