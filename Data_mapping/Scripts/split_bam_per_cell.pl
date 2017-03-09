#!/usr/bin/perl

use strict;

MAIN : {
    my ($bam_file, $header_file, $bc_file, $out_dir, $est_num_cell_bcs) = @ARGV;
    if ((not defined $bam_file) ||
	(not defined $header_file) ||
	(not defined $bc_file) ||
	(not defined $out_dir) ||
	(not defined $est_num_cell_bcs)) {
	die ("Usage: ./split_to_sam.pl <bam file> <bc file> <out dir> <est # cell bcs>\n");
    }

    $est_num_cell_bcs += 100;

    # read header
    my @header_lines;
    print STDERR ("reading header\n");
    open(FILE, $header_file) || die ("could not open file ($header_file)\n");
    while (my $line = <FILE>) {
	push(@header_lines, $line);
    }
    close(FILE) || die ("could not close file ($header_file)\n");

    # get list of the top barcodes
    print STDERR ("getting top barcodes\n");
    my $bc_to_lines;
    open(FILE, $bc_file) || die ("could not open file ($bc_file)\n");
    while (my $line = <FILE>) {
	chomp $line;
	my ($count, $cell_bc) = split(/\t/, $line);

	if (scalar (keys %$bc_to_lines) < $est_num_cell_bcs) {
	    $bc_to_lines->{$cell_bc}->{keep} = 1;
	} else {
	    last;
	}
    }
    close(FILE) || die ("could not close file ($bc_file)\n");

    # read bam file
    print STDERR ("reading bam file\n");
    open(FILE, "samtools view $bam_file |") || die ("could not open file ($bam_file)\n");
    while (my $line = <FILE>) {
	my ($id, $flag, $chr, $loc, $map, $cigar, $d7, $d8, $d9, $d10, $d11, @rest) = split(/\t/, $line);
	foreach my $rest (@rest) {
	    my ($tag, $type, $val) = split(/:/, $rest);
	    if ($tag eq "XC") {
		if (defined $bc_to_lines->{$val}) {
		    push(@{ $bc_to_lines->{$val}->{lines} }, $line);
		}
		last;
	    }
	}
    }
    close(FILE) || die ("could not close file ($bam_file)\n");

    print STDERR ("printing\n");
    foreach my $cell_id (keys %{ $bc_to_lines }) {
	print STDERR ("printing $cell_id\n");
	my $out_file = "$out_dir/$cell_id.bam";

	open(SAMTOOLS, "| samtools view -bS -o $out_file -") || die ("could not open samtools\n");
	print SAMTOOLS (@header_lines);
	print SAMTOOLS (@{ $bc_to_lines->{$cell_id}->{lines} });
	close(SAMTOOLS) || print STDERR  ("could not close samtools on file $out_file\n");
    }
}
