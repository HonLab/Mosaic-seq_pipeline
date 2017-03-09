#!/usr/bin/perl

use strict;

MAIN : {
    my ($matrix_file) = @ARGV;
    if (not defined $matrix_file) {
	die ("Usage: cat ids | ./filter_matrix_by_cell.pl <matrix_file>\n");
    }

    # read counts
    my $counts;
    open(FILE, $matrix_file) || die ("could not open file ($matrix_file)\n");
    while (my $line = <FILE>) {
	chomp $line;
	my @line = split(/\t/, $line);
	push(@{ $counts }, \@line);
    }
    close(FILE) || die ("could not close file ($matrix_file)\n");

    # read ids
    my %exclude_ids;
    while (my $line = <STDIN>) {
	chomp $line;
	$exclude_ids{$line}++;
    }

    my $keep;
    for(my $j=1 ; $j < scalar( @{ $counts->[0] }) ; $j++) { # first line is header
	my $id = $counts->[0]->[$j];

	if ($exclude_ids{$id}) {
	    $keep->{$j} = 0;
	} else {
	    $keep->{$j} = 1;
	}
    }
    $keep->{0} = 1;

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
