#!/usr/bin/perl

use strict;

MAIN : {
    my ($suffix) = @ARGV;
    foreach my $prefix (<STDIN>) {
	chomp $prefix;
	print($prefix . "/" . $suffix . "\n");
    }
}
