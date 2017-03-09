#!/usr/bin/perl

use strict;

MAIN : {
    my ($lib) = @ARGV;
    if (not defined $lib) {
	die ("lib must be defined\n");
    }

    foreach my $line (<STDIN>) {
	chomp $line;
	my @line = split(/\t/, $line);

	for(my $i=0 ; $i < scalar(@line) ; $i++) {
	    $line[$i] = $lib . "." . $line[$i];
	}
	print(join("\t", @line) . "\n");
    }
}
