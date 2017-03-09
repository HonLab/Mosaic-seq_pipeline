#!/bin/tcsh

#SBATCH -J barcode_split
#SBATCH -p super
#SBATCH -N 1
#SBATCH -t 0-12:00:00
#SBATCH --mail-user=Gary.Hon@utsouthwestern.edu
#SBATCH --mail-type=end
#SBATCH -o viral_split.out.%j
#SBATCH -e viral_split.err.%j

# scripts and directories
set SHARED_DIR        = /project/GCRB/Hon_lab/shared
set SCRIPT_DIR        = /project/GCRB/Hon_lab/shared/pipelines/Drop-Seq/2016-04-04/scripts
set BARCODE_LIST      = $PWD/known_barcodes.list

module add python/3.4.x-anaconda

foreach name (JD106-old GH108 GH109 GH110 GH111 GH112 JD106-A)# JD106-B RX134-1 RX134-1A RX134-1B RX134-1C RX135-1 RX136-1 RX137-1 RX138-1 RX139-2)
  echo $name

  set out = $name.dropseq_pipe

  set STAMPS_SELECTED = $out/STAMPs/STAMPs_info_mapQ10_top500_virus_barcodes_editdist1.txt.STAMPs
  rm -rf $out/STAMPs/split
  mkdir  $out/STAMPs/split

  cd $out/STAMPs/split
  cp ../STAMPs_info_mapQ10_top500_virus_barcodes_editdist1.txt.STAMPs $name.STAMPs.txt

  $SCRIPT_DIR/select_viral_barcodes_moderate2.py $name.STAMPs.txt\
                                                $BARCODE_LIST

  $SCRIPT_DIR/select_viral_barcodes_loose.py $name.STAMPs.txt\
                                             $BARCODE_LIST

  $SCRIPT_DIR/select_viral_barcodes_strict.py $name.STAMPs.txt\
                                              $BARCODE_LIST

  cd ../../..

  rm -rf $out/genic_counts/STAMPs.high_cov.viral_split
  mkdir  $out/genic_counts/STAMPs.high_cov.viral_split
  rm -rf $out/genic_counts/STAMPs.viral_split
  mkdir  $out/genic_counts/STAMPs.viral_split

  foreach stamp (STAMPs STAMPs.high_cov)
    foreach kind (loose moderate strict)
      foreach bc (`cat $BARCODE_LIST`)
        set bc_file = $out/STAMPs/split/$name.STAMPs.$bc.$kind.txt
        if (-e $bc_file) then
          $SCRIPT_DIR/select_barcodes_read_count.py $bc_file\
                                                    $out/genic_counts/$stamp/$name.q10.cleaned.cell_bc_merged.nodup.exons.$stamp\
            > $out/genic_counts/$stamp.viral_split/$name.STAMPs.$bc.$kind.readcount.txt

          cut -f7- $out/genic_counts/$stamp.viral_split/$name.STAMPs.$bc.$kind.readcount.txt\
            > $out/genic_counts/$stamp.viral_split/$name.STAMPs.$bc.$kind.readcount.txt.counts_only
        endif
      end

      # make a placeholder CPM file containing all IDs
      echo -n "#Geneid\t" > $out/genic_counts/$stamp.viral_split/$name.STAMPs.all.$kind.dummy_CPM.txt
      paste $out/genic_counts/$stamp.viral_split/$name.STAMPs.*.$kind.*txt\
        | head -1\
        | sed -e 's/\t/\n/g'\
        | grep bam\
        | sed -e 's/.bam//g'\
        | sort\
        | uniq\
        | perl -ne 'chomp; print("$_\t");'\
        >> $out/genic_counts/$stamp.viral_split/$name.STAMPs.all.$kind.dummy_CPM.txt
    end
  end

  # get read counts
  foreach kind (STAMPs.high_cov.viral_split STAMPs.viral_split)
    $SCRIPT_DIR/count_reads.pl $out/genic_counts/$kind\
                               $out/bam/split.bam\
                               $out/bam/depth/$kind
  end

  # get stats
  cat $out/genic_counts/STAMPs/$name.q10.cleaned.cell_bc_merged.nodup.exons.STAMPs.counts_only\
    | $SCRIPT_DIR/get_stats2.pl $name\
                                $out/reports/$name.num_reads_seq\
                                $out/reports/$name.num_mapped\
                                $out/reports/$name.num_nodup\
                                $out/bam/depth\
    > $out/reports/$name.stats2
end
