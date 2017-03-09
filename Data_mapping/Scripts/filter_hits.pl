#!/usr/bin/perl

use strict;

MAIN : {
    my ($query_bc, $counts_file) = @ARGV;
    if (not defined $query_bc) {
	die ("Usage: cat file | ./filter_hits.pl <query_bc>\n");
    }

    my $all_num_hits;
    open(FILE, $counts_file) || die ("could not open file ($counts_file)\n");
    while (my $line = <FILE>) {
	chomp $line;
	my ($cell_bc, $viral_bc, $counts, $num_hits) = split(/\t/, $line);

	my $id = join("\t", $cell_bc, $viral_bc, $counts);
	$all_num_hits->{ $id } = $num_hits;
    }
    close(FILE) || die ("could not close file ($counts_file)\n");

    while (my $line = <STDIN>) {
	chomp $line;
	my ($cell_bc, $all_sgRNA_bcs, $counts) = split(/\t/, $line);
	my $id = join("\t", $cell_bc, $all_sgRNA_bcs, $counts);

	my $num_hits = $all_num_hits->{$id};
	if (not defined $num_hits) {
	    die("num hits not defined for: $line\n");
	}

	my @counts = split(/,/, $counts);
	my @all_sgRNA_bcs = split(/,/, $all_sgRNA_bcs);
	for(my $i=0 ; $i < $num_hits ; $i++) {
	    if ($all_sgRNA_bcs[$i] eq $query_bc) {
		print(join("\t", $line, $counts[$i], $num_hits) . "\n");
		last;
	    }
	}
    }
}
