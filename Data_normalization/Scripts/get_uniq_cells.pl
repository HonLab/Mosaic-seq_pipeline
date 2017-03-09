#!/usr/bin/perl

use strict;

 MAIN:{
     my ($uniq_bc_file) = @ARGV;
     my @uniq_bc_array;
     
     open (my $fh1, $uniq_bc_file) || die ("Could not open file $uniq_bc_file");
     while (my $line1 = <$fh1>){
	 chomp $line1;
	 push(@uniq_bc_array, $line1);
     }
     close ($fh1) || die ("Could not close file $fh1");

     my $input_line1 = <STDIN>;
     chomp $input_line1;
     my @id_array = split(/\t/, $input_line1);

     my @uniq_bc_counter = get_uniq_bc(\@uniq_bc_array, \@id_array);

     print_uniq_bc($input_line1, \@uniq_bc_counter);

     while (my $line2 = <STDIN>){
	 chomp $line2;
	 print_uniq_bc($line2, \@uniq_bc_counter);
     }
}


sub get_uniq_bc{
 
    my ($uniq_bc_ref, $id_ref) = @_;
    my @uniq_bc_array = @{$uniq_bc_ref};
    my @id_array = @{$id_ref};
    
    my @uniq_bc_counter;
    
    for (my $i = 0; $i < scalar(@uniq_bc_array); $i++) {
	my $current_query_id = @uniq_bc_array[$i];
	
	for (my $j = 0; $j < scalar(@id_array); $j++){
	    my $current_id = @id_array[$j];

	    if ($current_id eq $current_query_id){
		push(@uniq_bc_counter, $j);
		last;
	    }
	}
    }
    return(@uniq_bc_counter);
}

sub print_uniq_bc{
    my ($line, $uniq_bc_ref) = @_;

    my @uniq_bc_counter = @{$uniq_bc_ref};
    my @data_array = split (/\t/, $line);

    for (my $i = 0; $i < scalar(@uniq_bc_counter); $i++){
	my $output_id = @uniq_bc_counter[$i];
	print(@data_array[$output_id]."\t");
    }
    print("\n");
}
