#!/usr/bin/perl

use strict;

#my $FLAG = 0;
#$| = 1;

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
    my $line_count = 0;
#    while (my $line = <>) {
    open(FILE, "samtools view $bam_file |") || die ("could not open file ($bam_file)\n");
#    open(FILE, $bam_file) || die ("could not open file ($bam_file)\n");
    while (my $line = <FILE>) {
	if (++$line_count % 10000 == 0) {
	    print STDERR (".");
	}
	if ($line_count % 100000 == 0) {
	    print STDERR (" $line_count\n");
	}

#	if ($line_count >= 315000) {
#	    $FLAG = 1;
#	}

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

#    print ("process_reads called:\n");
#    foreach my $strand ("+", "-") {
#	if (defined $current_lines->{$strand}) {
#	    print("$strand strand: " . scalar( @{ $current_lines->{$strand} }) . " reads\n");
#	    foreach my $line (@{ $current_lines->{$strand} }) {
#		print("  $line\n");
#	    }
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

#	        print ("cell_id = $cell_id\n");
#		print("molecular_bcs (" . scalar(@molecular_bcs) . ") = @molecular_bcs\n");

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

#	    if ($FLAG) {print STDERR ("calling get_uniq_bcs\n");}
#	    if (scalar(@molecular_bcs) >= 1000) {
#		print STDERR ("extra large size: " . scalar(@molecular_bcs) . ", cell = $cell_id\n");
#	    }
	    my $uniq_bcs = get_uniq_bcs($distances);
#	    if ($FLAG) {print STDERR ("done calling get_uniq_bcs\n");}

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

    my %cell_to_edit_dist_1;
    for(my $i=0 ; $i < scalar(@cells) ; $i++) {
	my $num_edit_dist_1 = 0;
	for(my $j=0 ; $j < scalar(@cells) ; $j++) {
	    if ($j == $i) {
		next;
	    }

	    if ($distances->{$i}->{$j} <= 1) {
		$num_edit_dist_1++;
	    }
	}

	$cell_to_edit_dist_1{$i} = $num_edit_dist_1;
    }
    my @sorted_ids = sort {$cell_to_edit_dist_1{$b} <=> $cell_to_edit_dist_1{$a}} keys %cell_to_edit_dist_1;

#    while (scalar(keys %exclude) != scalar(@cells)) {
#	print ("exclude size = " . scalar(keys %exclude) . ", final =  " . scalar(@cells) . "\n");
    foreach my $max_i (@sorted_ids) {
	if ($exclude{$max_i}) {
	    next;
	}

#	my $max_val = undef;
#	my $max_i = undef;
#	foreach my $i (keys %cell_to_edit_dist_1) {
#	    if ($exclude{$i}) {
#		next;
#	    }
#	    
#	    if ((not defined $max_i) ||
#		($cell_to_edit_dist_1{$i} > $max_val)) {
#		$max_val = $cell_to_edit_dist_1{$i};
#		$max_i = $i;
#	    }
#	}

	push(@uniq_bcs, $max_i);

	# found best one, so increment
	$num_uniq_bcs++;

	# exclude all other cells
	$exclude{$max_i} = 1;
	for(my $j=0 ; $j < scalar(@cells) ; $j++) {
	    if ($distances->{$max_i}->{$j} <= 1) {
		$exclude{$j} = 1;
	    }
	}
    }

    return \@uniq_bcs;
}
