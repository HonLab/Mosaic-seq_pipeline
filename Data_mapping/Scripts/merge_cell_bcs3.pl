#!/usr/bin/perl

use strict;

################################################################################
#
# MAIN
#
################################################################################
MAIN : {
    my ($bam_file, $min_reads_for_big_bc, $min_reads_for_sub_bc, $header_file) = @ARGV;
    if ((not defined $bam_file) ||
	(not defined $min_reads_for_big_bc) ||
	(not defined $min_reads_for_sub_bc) ||
	(not defined $header_file)) {
	die ("Usage: ./merge_cell_bcs.pl <bam file> <min reads for big bc> <min reads for sub bc>\n");
    }

    # print header
    open(FILE, $header_file) || die ("could not open file ($header_file)\n");
    while (my $line = <FILE>) {
	print($line);
    }
    close(FILE) || die ("could not close file ($header_file)\n");

    # get all the barcodes
    my $barcode_count;
    open(FILE, "samtools view $bam_file |") || die ("could not open file ($bam_file)\n");
    while (my $line = <FILE>) {
	my ($cell_bc_first, $cell_bc_last) = ($line =~ /\tXC:Z:(\S+)(\S)\t/);
	$barcode_count->{$cell_bc_first}++;
    }
    close(FILE) || die ("could not close file ($bam_file)\n");

    # print stats
    print STDERR ("number of unique barcodes = " . scalar(keys %$barcode_count) . "\n");
    my $num_gt_100  = 0;
    my $num_gt_1000  = 0;
    my $num_gt_10000 = 0;
    foreach my $cell_bc (keys %$barcode_count) {
	if ($barcode_count->{$cell_bc} >= 100) {
	    $num_gt_100++;
	}
	if ($barcode_count->{$cell_bc} >= 1000) {
	    $num_gt_1000++;
	}
	if ($barcode_count->{$cell_bc} >= 10000) {
	    $num_gt_10000++;
	}
    }
    print STDERR ("number of unique barcodes >100  = $num_gt_100\n");
    print STDERR ("number of unique barcodes >1k   = $num_gt_1000\n");
    print STDERR ("number of unique barcodes >10k  = $num_gt_10000\n");

    # find out which barcodes should be merged
    my @sorted_barcode_count = sort {$barcode_count->{$b} <=> $barcode_count->{$a}} keys (%$barcode_count);
    my $big_id_to_small_id;
    my $small_id_to_big_id;
    my %exclude;
    my %big;
    for(my $i=0 ; $i < scalar(@sorted_barcode_count) ; $i++) {
	my $big_id = $sorted_barcode_count[$i];
	if ($exclude{$big_id}) {
	    next;
	}

	if ($barcode_count->{$big_id} < $min_reads_for_big_bc) {
	    last;
	}

	$exclude{$big_id} = 1;
	#$small_id_to_big_id->{$big_id} = $big_id;
	$big{$big_id} = 1;

#	print STDERR ("big id     = $big_id (count = " . $barcode_count->{$big_id} . ")\n");

	my $extra = 0;
	for(my $j=$i+1 ; $j < scalar(@sorted_barcode_count) ; $j++) {
	    my $small_id = $sorted_barcode_count[$j];
	    if ($exclude{$small_id}) {
		next;
	    }

	    if ($barcode_count->{$small_id} < $min_reads_for_sub_bc) {
		last;
	    }

	    if (edit_dist($big_id, $small_id) <= 1) {
		$big_id_to_small_id->{$big_id}->{$small_id} = 1;
		$small_id_to_big_id->{$small_id} = $big_id;
		$exclude{$small_id} = 1;

		$extra += $barcode_count->{$small_id};
#		print STDERR ("  small id = $small_id (count = " . $barcode_count->{$small_id} . ")\n");
	    }
	}

	my $percent = sprintf("%.4f", $extra / $barcode_count->{$big_id} * 100);
#	print STDERR ("  extra added = $extra (+ $percent %)\n");
    }

    print STDERR ("\n\ntotal number of big ids = " . scalar(keys %$big_id_to_small_id) . "\n\n");

    # print sam file with updated cell barcodes
    open(FILE, "samtools view $bam_file |") || die ("could not open file ($bam_file)\n");
    my $counts;
    while (my $line = <FILE>) {
	my ($cell_bc_first, $cell_bc_last) = ($line =~ /\tXC:Z:(\S+)(\S)\t/);

	if ($big{$cell_bc_first}) {
	    print($line);
	    $counts->{"big_id"}++;
	    next;
	}

	if (not defined $small_id_to_big_id->{$cell_bc_first}) {
	    #die ("small_id_to_big_id { $cell_bc_first } must be defined\n");
	    print($line);
	    $counts->{"unchanged small"}++;
	    next;
	}

	my $new_cell_id = $small_id_to_big_id->{$cell_bc_first} . $cell_bc_last;

	my $new_line = $line;
	$new_line =~ s/\tXC:Z:\S+\t/\tXC:Z:$new_cell_id\t/;
#	print("oln: $line");
	print($new_line);

	$counts->{"small to big"}++;
    }
    close(FILE) || die ("could not close file ($bam_file)\n");

    my $sum = $counts->{"big_id"} + $counts->{"unchanged small"} + $counts->{"small to big"};
    foreach my $kind (sort keys %$counts) {
	print STDERR (join("\t", $kind, sprintf("%.4f", $counts->{$kind} / $sum * 100)) . "\n");
    }
}

################################################################################
#
# sub edit_dist
#
################################################################################
sub edit_dist {
    my ($bc1, $bc2) = @_;

    # check args
    if (not defined $bc1) {
	die ("bc1 must be defined\n");
    }
    if (not defined $bc2) {
	die ("bc2 must be defined\n");
    }

    if (length($bc1) != length($bc2)) {
	die ("lengths must be same ($bc1, $bc2)\n");
    }

    my $edit_dist = 0;
    for(my $i=0 ; $i < length($bc1) ; $i++) {
	if (substr($bc1, $i, 1) ne substr($bc2, $i, 1)) {
	    $edit_dist++;
	}
    }

    return $edit_dist;
}

################################################################################
#
# sub get_uniq_bcs
#
################################################################################
sub get_uniq_bcs {
    my ($distances) = @_;

    # check args
    if (not defined $distances) {
	die ("distances must be defined\n");
    }

    my %exclude = ();
    my @cells = keys %$distances;
    my $num_uniq_bcs = 0;

    while (scalar(keys %exclude) != scalar(@cells)) {

	my %cell_to_edit_dist_1;
	for(my $i=0 ; $i < scalar(@cells) ; $i++) {
	    if ($exclude{$i}) {
		next;
	    }

	    my $num_edit_dist_1 = 0;
	    for(my $j=0 ; $j < scalar(@cells) ; $j++) {
		if (($j == $i) ||
		    ($exclude{$j})) {
		    next;
		}

		if ($distances->{$i}->{$j} <= 1) {
		    $num_edit_dist_1++;
		}
	    }

	    $cell_to_edit_dist_1{$i} = $num_edit_dist_1;
	}

	# choose the cell with most edit distance of 1
	if (scalar(keys %cell_to_edit_dist_1) > 0) {

	    my $max_val = undef;
	    my $max_i = undef;
	    foreach my $i (keys %cell_to_edit_dist_1) {
		if ((not defined $max_i) ||
		    ($cell_to_edit_dist_1{$i} > $max_val)) {
		    $max_val = $cell_to_edit_dist_1{$i};
		    $max_i = $i;
		}
	    }

	    # found best one, so increment
	    $num_uniq_bcs++;

	    # exclude all other cells
	    $exclude{$max_i} = 1;
	    for(my $j=0 ; $j < scalar(@cells) ; $j++) {
		if ($distances->{$max_i}->{$j} == 1) {
		    $exclude{$j} = 1;
		}
	    }

	# all cells have edit dist greater than 1, so output them all
	} else {

	    $num_uniq_bcs += (scalar @cells) - scalar(keys %exclude);
	    last;
	}
    }

    return $num_uniq_bcs;
}
