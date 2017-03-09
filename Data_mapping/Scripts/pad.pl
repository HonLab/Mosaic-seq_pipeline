#!/usr/bin/perl

use strict;

MAIN : {
    my $max = 0;
    my @lines;
    while (my $line = <>) {
	chomp $line;
	push(@lines, $line);

	my @line = split(/\t/, $line);
	if (scalar(@line) > $max) {
	    $max = scalar(@line);
	}
    }

    foreach my $line (@lines) {
	my @line = split(/\t/, $line);
	my $orig = scalar(@line);

	while (scalar(@line) < $max) {
	    push(@line, 0);
	}

	print(join("\t", @line) . "\n");
    }
}
