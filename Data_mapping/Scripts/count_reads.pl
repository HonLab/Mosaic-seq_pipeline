#!/usr/bin/perl

use strict;

MAIN : {
    my ($in_dir, $split_bam_dir, $out_root) = @ARGV;
    if ((not defined $in_dir) ||
	(not defined $split_bam_dir) ||
	(not defined $out_root)) {
	die ("Usage: ./count_reads.pl <in dir> <split bam dir> <out dir>\n");
    }

    if (not -e $in_dir) {
	die ("in_dir ($in_dir) must exist\n");
    }
    if (not -e $split_bam_dir) {
	die ("split_bam_dir ($split_bam_dir) must exist\n");
    }

    # read all bam files
    my $cell_barcode_seq_to_depth;
    open(LS, "ls $split_bam_dir | grep counts |") || die ("could not ls into $split_bam_dir\n");
    while (my $bam = <LS>) {
	chomp $bam;

	my ($cell_barcode_seq) = split(/\./, $bam);
	open(FILE, "$split_bam_dir/$bam") || die ("could not open file ($split_bam_dir/$bam)\n");
	while (my $depth = <FILE>) {
	    chomp $depth;
	    $cell_barcode_seq_to_depth->{$cell_barcode_seq} = $depth;
	}
	close(FILE) || die ("could not close file ($split_bam_dir/$bam)\n");
    }
    close(LS) || die ("could not close $split_bam_dir\n");

    # get query cell barcodes
    my @depths;
    open(LS, "ls $in_dir |") || die ("could not ls ($in_dir)\n");
    my $initial_ls = <LS>;
    close(LS);# || die ("could not close ls\n");

#    print STDERR ("initial ls ($initial_ls), defined = " . (defined $initial_ls) . "\n");

    if (defined $initial_ls) {

	# if viral, only do moderate
	if ($in_dir =~ /viral/) {
	    open(LS, "ls $in_dir/*moderate*CPM.txt |") || die ("could not ls ($in_dir/*CPM.txt)\n");
	} else {
	    open(LS, "ls $in_dir/*CPM.txt |") || die ("could not ls ($in_dir/*CPM.txt)\n");
	}
	while (my $file = <LS>) {
	    chomp $file;

	    # get first line
	    open(FILE, "$file") || die ("could not open file ($file)\n");
	    while (my $line = <FILE>) {
		chomp $line;

		my @line = split(/\t/, $line);
		for(my $i=1 ; $i < scalar(@line) ; $i++) {
		    push(@depths, $cell_barcode_seq_to_depth->{$line[$i]});
		}

		# only do first line
		last;
	    }
	    close(FILE) || die ("could not close file ($file)\n");

	    # only do once
	    last;
	}
	close(LS) || die ("could not close ls ($in_dir/*CPM.txt)\n");
    }

    # calculate median
    my @sorted = sort {$a <=> $b} @depths;
    my $num_cells = scalar( @sorted );
    my $median_depth = 0;
    if (scalar(@sorted) > 0) {
	if (scalar(@sorted) % 2 == 0) {
	    # even, so take average of two middle points
	    $median_depth = ($sorted[int($num_cells/2)] + $sorted[int($num_cells/2)-1])/2;
	} else {
	    # odd, so take mid
	    $median_depth = $sorted[int($num_cells/2)]
	}
    }

    # print
    open(FILE, ">$out_root") || die ("could not open out_root ($out_root)\n");
    print FILE (join("\t", @depths) . "\n");
    close(FILE) || die ("could not close file ($out_root)\n");

    open(FILE, ">$out_root.median") || die ("could not open out_root ($out_root.median)\n");
    print FILE ("$median_depth\n");
    close(FILE) || die ("could not close file ($out_root.median)\n");
}
