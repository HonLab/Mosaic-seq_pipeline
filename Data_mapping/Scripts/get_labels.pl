#!/usr/bin/perl

use strict;

MAIN : {
    my ($subset_file, $ids_list_file) = @ARGV;
    if ((not defined $subset_file) ||
	(not defined $ids_list_file)) {
	die ("Usage: ./get_labels.pl <subset file> <ids_list_file>\n");
    }

    my $data;
    open(FILE, $subset_file) || die ("could not open file ($subset_file)\n");
    while (my $line = <FILE>) {
	chomp $line;
	my ($cell_bc, $viral_bcs, $viral_bc_counts, $bc_count, $num_states) = split(/\t/, $line);
	my $id = $cell_bc . ".bam";

	my @viral_bc_counts = split(/,/, $viral_bc_counts);
	my $sum = 0;
	foreach my $c (@viral_bc_counts) {
	    $sum += $c;
	}

	$data->{$id}->{bc_count   } = $bc_count;
	$data->{$id}->{bc_percent } = sprintf("%.4f", $bc_count / $sum * 100);
	$data->{$id}->{num_states } = $num_states;
    }
    close(FILE) || die ("could not close file ($subset_file)\n");

    open(FILE, $ids_list_file) || die ("could not open file ($ids_list_file)\n");
    my $line = <FILE>;
    chomp $line;
    my @ids = split(/\t/, $line);
    close(FILE) || die ("could not close file ($subset_file)\n");

    foreach my $id (@ids) {
	print(join("\t",
		   $id,
		   $data->{$id}->{bc_count  } || 0,
		   $data->{$id}->{bc_percent} || 0,
		   $data->{$id}->{num_states} || 0,
	      ) . "\n");
    }
}
