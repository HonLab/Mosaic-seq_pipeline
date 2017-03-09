#!/bin/tcsh

#SBATCH --job-name=matrix_norm_2                          # job name
#SBATCH --partition=super                                 # select partion from 128GB, 256GB, 384GB, GPU and super
#SBATCH --nodes=1                                         # number of nodes requested by user
#SBATCH --time=0-20:00:00                                 # run time, format: D-H:M:S (max wallclock time)
#SBATCH --output=serialJob.%j.out                         # standard output file name
#SBATCH --error=serialJob.%j.time                         # standard error output file name
#SBATCH --mail-user=Shiqi.xie1@utsouthwestern.edu         # specify an email address
#SBATCH --mail-type=end                                   # send email when job status change (start, end, abortion and etc.)

echo "hello world"

#create symbolic links to every library
set TARGET_DIR = /project/GCRB/Hon_lab/s159317/projects/Drop-Kick-Seq/data/mapping/Drop-Seq.pipe=2016-11-01

foreach DATASET (\
		    2016-04-21-CRI-NextSeq\
		    2016-10-05-McDermott-NextSeq\
#		    2016-04-22-McDermott-NextSeq.run1\
#		    2016-11-08-McDermott-NextSeq.KRAB.links\
#		    2016-04-22-McDermott-NextSeq.run2\
#		    2016-10-24-McDermott-NextSeq\
#		    2016-09-02-McDermott-NextSeq\
#		    2016-09-15-McDermott-NextSeq\
#		    2016-10-31-McDermott-NextSeq\
		)
    foreach LIB (`ls $TARGET_DIR/$DATASET/ | grep dropseq.pipe`)
	ln -s $TARGET_DIR/$DATASET/$LIB $LIB.norm
    end
end

###################################################################################################
#
#for each library, sort by gene position, get the data matrix, normalize by median of every gene
#
###################################################################################################
#rm -rf sorted control_cells
mkdir sorted

foreach lib (`ls -d *.norm`)
  echo "sorting $lib"

  mkdir sorted/$lib
  set dir = $lib/genic_counts/split_by_sgRNA.control_uniq
  #sort genes by chromosomal coordinate
  foreach file (`ls $dir | grep sgRNA | grep -v vals`)
    # header
    echo $file
    head -1 $dir/$file\
        | ./append_lib_name.pl $lib\
      > sorted/$lib/$file

    # sort by gene location
    tail -n+2 $dir/$file\
      | perl -ne 'chomp; @line = split(/\t/, $_); @col1 = split(/;/, $line[0]); $tss = ($col1[3] eq "+") ? $col1[4] : $col1[5]; print(join("\t", $col1[2], $tss, $_) . "\n");'\
      | sort -k1,1V -k2,2n\
      | cut -f3-\
      >> sorted/$lib/$file

    cat sorted/$lib/$file\
      | perl -ne '@a = split; shift(@a); print(join("\t", @a) . "\n");'\
      > sorted/$lib/$file.vals

  end

#    ls -1 sorted/$lib/*vals | split -l 1000 -d - ./sorted/$lib/lists
#    bash -c  "cd sorted/$lib; for list in lists*; do paste $(cat $list) > merge${list##lists}; cd ../..; done"
#    paste ./sorted/$lib/merge* > ./sorted/$lib/all_sgRNA.cell_barcodes

    # get control cells
    cat ./$lib/genic_counts/split_by_sgRNA.control_uniq/cell_barcodes.control\
      | sort\
      | uniq\
      > sorted/$lib/$lib.control.cell_barcode

    #get all unique sgRNA barcodes
    paste ./sorted/$lib/*vals\
    | head -1\
    | sed -e 's/\t/\n/g'\
    | sort\
    | uniq\
    | perl -ne 'if(/^\s*$/){next;} print;'\
    > ./sorted/$lib/all_sgRNA.uniq.cell_barcodes
    
    #use cell barcode file to filter the uniq cells
    paste ./sorted/$lib/*.vals\
	| ./get_uniq_cells.pl\
	  ./sorted/$lib/all_sgRNA.uniq.cell_barcodes\
	  | perl -ne '$_=~ s/^\s+|\s+$//g; print($_."\n");'\
	> ./sorted/$lib/combined.all.uniq.vals
	
    #get rid of the header
    tail -n+2 ./sorted/$lib/combined.all.uniq.vals\
	> ./sorted/$lib/combined.all.uniq.vals.no_header

    #normalize the data by the meidan of each gene
    cd ./sorted/$lib/
    cp ../../*.m .
    matlab -nodesktop -nosplash -nodisplay < normalization.m

    #paste the header to the normalized matrix
    head -1 combined.all.uniq.vals\
	| cat - normalized.all.matrix.txt.no_header\
	| perl -ne '$_=~ s/^\s+|\s+$//g; print($_."\n");'\
	| perl -ne '$_=~ s/\s+/\t/g; print($_."\n");'\
	> normalized.all.matrix.txt

    cd ../../
    
end

