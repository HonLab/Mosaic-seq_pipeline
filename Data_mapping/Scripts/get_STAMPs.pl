#!/usr/bin/perl

use strict;

MAIN : {
    my ($num_stamps_slope_cutoff, $num_reads_file) = @ARGV;
    if (not defined $num_stamps_slope_cutoff) {
	die ("Usage: cat file | ./get_STAMPs.pl <num stamps slope cutoff>\n");
    }

    # get total num reads
    open(FILE, $num_reads_file) || die ("could not open file ($num_reads_file)\n");
    my $total_num_reads = <FILE>;
    chomp $total_num_reads;
    close(FILE) || die ("could not close file ($num_reads_file)\n");

    my $cumsum = 0;
    my $printed = 0;
    while (my $line = <STDIN>) {
	my ($cell_id, $cell_bc, $num_reads) = split(/\t/, $line);

	my $last_cumsum = $cumsum;
	$cumsum += $num_reads / $total_num_reads;

	my $diff = $cumsum - $last_cumsum;
	if ($diff >= $num_stamps_slope_cutoff) {
	    print($line);
	    $printed++;
	} else {
	    last;
	}
    }

    print STDERR ("num STAMPs = $printed\n");
}
