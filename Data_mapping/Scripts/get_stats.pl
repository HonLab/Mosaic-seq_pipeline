#!/usr/bin/perl

use strict;

MAIN : {
    my ($name, $total_reads_file, $nodup_reads_file) = @ARGV;
    if ((not defined $name) ||
	(not defined $total_reads_file) ||
	(not defined $nodup_reads_file)) {
	die ("Usage: cat file | ./get_stats.pl <name> <total reads file> <nodup reads file>\n");
    }

    # get total num reads sequenced
    open(FILE, $total_reads_file) || die ("could not open file ($total_reads_file)\n");
    my $num_reads_seq = <FILE>;
    chomp $num_reads_seq;
    close(FILE) || die ("could not close file ($total_reads_file)\n");

    # get nodup
    open(FILE, $nodup_reads_file) || die ("could not open file ($nodup_reads_file)\n");
    my $num_nodup = <FILE>;
    chomp $num_nodup;
    close(FILE) || die ("could not close file ($nodup_reads_file)\n");

    # percent mono
    my $percent_map_q_mono = sprintf("%.2f", $num_nodup / $num_reads_seq * 100);

    # read counts
    my $counts;
    while (my $line = <STDIN>) {
	chomp $line;

	# line with just numbers
	if ($line =~ /^[\d\s]+$/) {
	    my @line = split(/\t/, $line);
	    push(@{ $counts }, \@line);
	}
    }

    # get:
    # 1. total num reads in genic table
    # 2. num reads / cell
    # 3. # STAMPs with >5% of genes having at least 10 counts
    my $num_genic = 0;
    my $num_reads_per_cell;
    my $num_genes_ge10_reads;
    for(my $i=0 ; $i < scalar( @{ $counts }) ; $i++) {
	for(my $j=0 ; $j < scalar( @{ $counts->[$i] }) ; $j++) {
	    $num_genic += $counts->[$i]->[$j];
	    $num_reads_per_cell->[$j] += $counts->[$i]->[$j];

	    if ($counts->[$i]->[$j] >= 10) {
		$num_genes_ge10_reads->[$j]++;
	    }
	}
    }
    my $genic_nodup_ratio = sprintf("%.2f", $num_genic / $num_nodup);

    # get:
    # 1. num cells
    # 2. median num reads/cell
    my $num_cells = scalar( @$num_reads_per_cell );
    my @sorted    = sort {$a <=> $b} @$num_reads_per_cell;
    my $median_num_reads_per_cell;
    if (scalar(@sorted) % 2 == 0) {
	# even, so take average of two middle points
	$median_num_reads_per_cell = ($sorted[int($num_cells/2)] + $sorted[int($num_cells/2)-1])/2;
    } else {
	# odd, so take mid
	$median_num_reads_per_cell = $sorted[int($num_cells/2)]
    }

    # get num STAMPs with >5% of genes having at least 10 counts
    my $num_genes = scalar( @{ $counts } );
    my $num_good_STAMPs = 0;
    for(my $j=0 ; $j < scalar(@$num_genes_ge10_reads) ; $j++) {
	if ($num_genes_ge10_reads->[$j] >= $num_genes * 0.05) {
	    $num_good_STAMPs++;
	}
    }
    my $percent_good_STAMPs = sprintf("%.2f", $num_good_STAMPs / $num_cells * 100);

    print(join("\t",
	       "#lib",
	       "#reads sequenced",
	       "#reads map_q_mono",
	       "%map_q_mono",
	       "#genic reads (in STAMPs)",
	       "ratio genic/map_q_mono",
	       "# STAMPs",
	       "median #reads / STAMP",
	       "# cells with >5% of genes having at least 10 counts",
	       "% good STAMPs") . "\n");

    print(join("\t",
	       $name,
	       $num_reads_seq,
	       $num_nodup,
	       $percent_map_q_mono . "%",
	       $num_genic,
	       $genic_nodup_ratio,
	       $num_cells,
	       $median_num_reads_per_cell,
	       $num_good_STAMPs,
	       $percent_good_STAMPs . "%",
	  ) . "\n");
}
