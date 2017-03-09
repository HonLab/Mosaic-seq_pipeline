#!/usr/bin/perl

use strict;

MAIN : {
    my ($exon_counts_file, $stamps_file) = @ARGV;
    if ((not defined $exon_counts_file) ||
	(not defined $stamps_file)) {
	die ("Usage: ./get_STAMP_counts.pl <exon_counts_file> <stamps_file>\n");
    }

    my $keep_barcodes;
    open(FILE, $stamps_file) || die ("could not open file ($stamps_file)\n");
    while (my $line = <FILE>) {
	my ($d1, $barcode) = split(/\t/, $line);
	$keep_barcodes->{$barcode}++;
    }
    close(FILE) || die ("could not close file ($stamps_file)\n");

    my $header = 1;
    my $keep_index;
    my @header;
    open(FILE, $exon_counts_file) || die ("could not open file ($exon_counts_file)\n");
    while (my $line = <FILE>) {
	chomp $line;
	my @line = split(/\t/, $line);

	if ($header) {
	    if ($line =~ /featureCounts/) {
		next;
	    } else {
		# process the header
		for(my $i=0 ; $i < scalar(@line) ; $i++) {
		    if ($line[$i] =~ /^([ACGTN]+)(.bam)?$/) {
			if ($keep_barcodes->{$1}) {
			    $keep_index->[$i]++;
			    push(@header, $line[$i]);
			}

		    # keep all other columns
		    } else {
			$keep_index->[$i]++;
			push(@header, $line[$i]);
		    }
		}

		print(join("\t", @header) . "\n");

		$header = 0;
	    }

	# process data lines
	} else {
	    my @out;
	    for(my $i=0 ; $i < scalar(@line) ; $i++) {
		if ($keep_index->[$i]) {
		    push(@out, $line[$i]);
		}
	    }

	    print(join("\t", @out) . "\n");
	}
    }
    close(FILE) || die ("could not close file ($exon_counts_file)\n");

}
