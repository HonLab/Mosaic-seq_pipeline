#!/usr/bin/perl

use strict;

MAIN : {
    my ($name, $total_reads_file, $mapped_reads_file, $nodup_reads_file, $median_dir) = @ARGV;
    if ((not defined $name) ||
	(not defined $total_reads_file) ||
	(not defined $mapped_reads_file) ||
	(not defined $nodup_reads_file) ||
	(not defined $median_dir)) {
	die ("Usage: cat file | ./get_stats.pl <name> <total reads file> <mapped reads file> <nodup reads file> <median dir>\n");
    }

    # get total num reads sequenced
    open(FILE, $total_reads_file) || die ("could not open file ($total_reads_file)\n");
    my $num_reads_seq = <FILE>;
    chomp $num_reads_seq;
    close(FILE) || die ("could not close file ($total_reads_file)\n");

    # get nodup
    open(FILE, $mapped_reads_file) || die ("could not open file ($mapped_reads_file)\n");
    my $num_mapped = <FILE>;
    chomp $num_mapped;
    close(FILE) || die ("could not close file ($mapped_reads_file)\n");

    # get nodup
    open(FILE, $nodup_reads_file) || die ("could not open file ($nodup_reads_file)\n");
    my $num_nodup = <FILE>;
    chomp $num_nodup;
    close(FILE) || die ("could not close file ($nodup_reads_file)\n");

    # percentages
    my $percent_mapped       = sprintf("%.2f", $num_mapped / $num_reads_seq * 100);
    my $percent_mono         = sprintf("%.2f", $num_nodup  / $num_mapped    * 100);
    my $percent_mono_div_all = sprintf("%.2f", $num_nodup  / $num_reads_seq * 100);

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
#    # 2. num reads / cell
#    # 3. # STAMPs with >5% of genes having at least 10 counts
    my $num_genic = 0;
#    my $num_reads_per_cell;
#    my $num_genes_ge10_reads;
    for(my $i=0 ; $i < scalar( @{ $counts }) ; $i++) {
	for(my $j=0 ; $j < scalar( @{ $counts->[$i] }) ; $j++) {
	    $num_genic += $counts->[$i]->[$j];
#	    $num_reads_per_cell->[$j] += $counts->[$i]->[$j];
#
#	    if ($counts->[$i]->[$j] >= 10) {
#		$num_genes_ge10_reads->[$j]++;
#	    }
	}
    }
    my $genic_nodup_ratio = sprintf("%.2f", $num_genic / $num_nodup);

#    # get:
#    # 1. num cells
#    # 2. median num reads/cell
#    my $num_cells = scalar( @$num_reads_per_cell );
#    my @sorted    = sort {$a <=> $b} @$num_reads_per_cell;
#    my $median_num_reads_per_cell;
#    if (scalar(@sorted) % 2 == 0) {
#	# even, so take average of two middle points
#	$median_num_reads_per_cell = ($sorted[int($num_cells/2)] + $sorted[int($num_cells/2)-1])/2;
#    } else {
#	# odd, so take mid
#	$median_num_reads_per_cell = $sorted[int($num_cells/2)]
#    }
#
#    # get num STAMPs with >5% of genes having at least 10 counts
#    my $num_genes = scalar( @{ $counts } );
#    my $num_good_STAMPs = 0;
#    for(my $j=0 ; $j < scalar(@$num_genes_ge10_reads) ; $j++) {
#	if ($num_genes_ge10_reads->[$j] >= $num_genes * 0.05) {
#	    $num_good_STAMPs++;
#	}
#    }
#    my $percent_good_STAMPs = sprintf("%.2f", $num_good_STAMPs / $num_cells * 100);

    # get num STAMPs
    my $num_STAMPs;
    my $median_reads;
    foreach my $kind ("STAMPs", "STAMPs.high_cov", "STAMPs.high_cov.viral_split", "STAMPs.viral_split") {
	$num_STAMPs->{  $kind} = 0;
	$median_reads->{$kind} = 0;

	# get num STAMPs
	my $count_file = "$median_dir/$kind";
	if (-e $count_file) {
	    open(FILE, $count_file) || die ("could not open file ($count_file)\n");
	    my $line = <FILE>;

	    if ($line !~ /^\s*$/) {
		my @counts = split(/\t/, $line);
		$num_STAMPs->{ $kind } = scalar(@counts);
	    }
	    close(FILE) || die ("could not close file ($count_file)\n");
        }

	# get median reads
	my $median_file = "$median_dir/$kind.median";
	if (-e $count_file) {
	    open(FILE, $median_file) || die ("could not open file ($median_file)\n");
	    my $median = <FILE>;
	    chomp $median;
	    $median_reads->{ $kind } = $median;
	    close(FILE) || die ("could not close file ($median_file)\n");
        }
    }

    print(join("\t",
	       "#lib",
	       "#reads sequenced",
	       "#reads mapped",
	       "%mapped",
	       "#reads mono",
	       "%mono (/mapped)",
	       "%mono (/all)",
	       "#genic reads (in STAMPs, Refseq)",
	       "ratio genic/mono",

	       "# STAMPs",
	       "median #reads / STAMP",
	       "# cells with >5% of genes having at least 10 counts (high cov)",
	       "median # reads / high cov STAMP",
	       "# cells (high cov) viral split",
	       "median # reads / high cov viral split STAMP",
	       "# cells (any cov) viral split",
	       "median # reads / any cov viral split STAMP") . "\n");

    print(join("\t",
	       $name,
	       $num_reads_seq,
	       $num_mapped,
	       $percent_mapped . "%",
	       $num_nodup,
	       $percent_mono . "%",
	       $percent_mono_div_all . "%",
	       $num_genic,
	       $genic_nodup_ratio,

	       $num_STAMPs->{"STAMPs"},
	       $median_reads->{"STAMPs"},
	       $num_STAMPs->{"STAMPs.high_cov"},
	       $median_reads->{"STAMPs.high_cov"},
	       $num_STAMPs->{"STAMPs.high_cov.viral_split"},
	       $median_reads->{"STAMPs.high_cov.viral_split"},
	       $num_STAMPs->{"STAMPs.viral_split"},
	       $median_reads->{"STAMPs.viral_split"}
	  ) . "\n");
}
