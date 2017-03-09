#!/usr/bin/perl

use strict;

MAIN : {
    my ($high_coverage_fraction) = @ARGV;
    if (not defined $high_coverage_fraction) {
	die ("Usage: cat file | ./get_stats.pl <high coverage fraction>\n");
    }

    # read counts
    my $counts;
    while (my $line = <STDIN>) {
	chomp $line;
	my @line = split(/\t/, $line);
	push(@{ $counts }, \@line);
    }

    # get: #STAMPs with >5% of genes having at least 10 counts
    my $num_genes_ge10_reads;
    for(my $i=1 ; $i < scalar( @{ $counts }) ; $i++) { # first line is header
	for(my $j=6 ; $j < scalar( @{ $counts->[$i] }) ; $j++) { # first 6 columns are headers
	    if ($counts->[$i]->[$j] >= 10) {
		$num_genes_ge10_reads->[$j]++;
	    }
	}
    }

    # get STAMPs with >x% of genes having at least 10 counts
    my $num_genes = scalar( @{ $counts } );
    my $keep;
    for(my $j=6 ; $j < scalar(@$num_genes_ge10_reads) ; $j++) { # first 6 columns are headers
	if ($num_genes_ge10_reads->[$j] >= $num_genes * $high_coverage_fraction) {
	    $keep->{$j}++;
	}
    }
    $keep->{0} = 1;
    $keep->{1} = 1;
    $keep->{2} = 1;
    $keep->{3} = 1;
    $keep->{4} = 1;
    $keep->{5} = 1;

    for(my $i=0 ; $i < scalar( @{ $counts }) ; $i++) { # first line is header
	my @line;
	for(my $j=0 ; $j < scalar( @{ $counts->[$i] }) ; $j++) {
	    if ($keep->{$j}) {
		push(@line, $counts->[$i]->[$j]);
	    }
	}

	print(join("\t", @line) . "\n");
    }
}
