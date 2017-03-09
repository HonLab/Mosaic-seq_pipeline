#!/usr/bin/perl

use strict;

################################################################################
#
# MAIN
#
################################################################################
MAIN : {
    my ($bam_file, $header_file) = @ARGV;
    if ((not defined $bam_file) ||
	(not defined $header_file)) {
	die ("Usage: ./rmdup2.pl <bam file> <header file>\n");
    }

    # print header
    open(FILE, $header_file) || die ("could not open file ($header_file)\n");
    while (my $line = <FILE>) {
	print($line);
    }
    close(FILE) || die ("could not close file ($header_file)\n");

    my $current_lines;
    my $current_chr = undef;
    my $current_loc = undef;
    my $qual_count;
#    while (my $line = <>) {
    open(FILE, "samtools view $bam_file |") || die ("could not open file ($bam_file)\n");
    while (my $line = <FILE>) {
	my (@line) = split(/\t/, $line);
	my $flag = $line[1];
	my $chr  = $line[2];
	my $loc  = $line[3];
	my $qual = $line[4];

	if ($qual < 10) {
	    $qual_count->{bad}++;
	    next;
	}
	$qual_count->{good}++;

	if (($current_chr ne $chr) ||
	    ($current_loc ne $loc)) {

	    if (defined $current_chr) {
		process_reads($current_lines);
	    }

	    $current_chr = $chr;
	    $current_loc = $loc;
	    $current_lines = undef;
	}

	my $strand = ($flag & 16) ? "-" : "+";
	push(@{ $current_lines->{$strand} }, $line);
    }
    close(FILE) || die ("could not close file ($bam_file)\n");

    if (defined $current_chr) {
	process_reads($current_lines);
    }
}

################################################################################
#
# sub process_reads
#
################################################################################
sub process_reads {
    my ($current_lines) = @_;

    # check args
    if (not defined $current_lines) {
	die ("current_lines must be defined\n");
    }

#    print ("process_reads called with:\n");
#    foreach my $strand ("+", "-") {
#	foreach my $line (@{ $current_lines->{$strand} }) {
#	    print("  $line\n");
#	}
#    }

    foreach my $strand ("+", "-") {

	# has by cell and mol ID
	my $cell_to_line;
	foreach my $line (@{ $current_lines->{$strand} }) {
	    my ($cell_id) = ($line =~ /XC:Z:(\S+)/);
	    my ($mol_id) = ($line =~ /XM:Z:(\S+)/);
	    my @line = split(/\t/, $line);
	    my $qual = $line[4];

	    if ((not defined $cell_to_line->{$cell_id}) ||
		(not defined $cell_to_line->{$cell_id}->{$mol_id}) ||
		($qual > $cell_to_line->{$cell_id}->{$mol_id})) {

		$cell_to_line->{$cell_id}->{$mol_id}->{qual} = $qual;
		$cell_to_line->{$cell_id}->{$mol_id}->{line} = $line;
	    }
	}

	# get distances between barcodes for given cell
	foreach my $cell_id (keys %$cell_to_line) {
	    my @molecular_bcs = keys %{ $cell_to_line->{$cell_id} };

#		print("molecular_bcs = @molecular_bcs\n");

	    # get distances
	    my $distances;
	    for(my $i=0 ; $i < scalar(@molecular_bcs) ; $i++) {
		for(my $j=$i ; $j < scalar(@molecular_bcs) ; $j++) {
		    my $dist = edit_dist($molecular_bcs[$i], $molecular_bcs[$j]);

#			print("i = $i, j = $j, dist = $dist\n");

		    $distances->{$i}->{$j} = $dist;
		    $distances->{$j}->{$i} = $dist;
		}
	    }

	    my $uniq_bcs = get_uniq_bcs($distances);

	    foreach my $i (@$uniq_bcs) {
		my $mol_id = $molecular_bcs[$i];
		print($cell_to_line->{$cell_id}->{$mol_id}->{line});
	    }
	}
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
    my @uniq_bcs;

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

	my $max_val = undef;
	my $max_i = undef;
	foreach my $i (keys %cell_to_edit_dist_1) {
	    if ((not defined $max_i) ||
		($cell_to_edit_dist_1{$i} > $max_val)) {
		$max_val = $cell_to_edit_dist_1{$i};
		$max_i = $i;
	    }
	}

	push(@uniq_bcs, $max_i);

	# found best one, so increment
	$num_uniq_bcs++;

	# exclude all other cells
	$exclude{$max_i} = 1;
	for(my $j=0 ; $j < scalar(@cells) ; $j++) {
	    if ($distances->{$max_i}->{$j} == 1) {
		$exclude{$j} = 1;
	    }
	}
    }

    return \@uniq_bcs;
}
