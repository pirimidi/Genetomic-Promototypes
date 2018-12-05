#!/usr/bin/perl

use strict;
use DNA::Promoter;

if (@ARGV<1) {
  die "\n--> Error: not correct input!\n--> Usage: perl test.pl config_file\n\n";
}

my($config_file) = &DNA::IO::validate_cf($ARGV[0]);

my(@FIELD) = (DNA::IO::FIELD1, DNA::IO::FIELD2, "Finalized_map", "Affinity_file", "Out_dir"); 
my(@TYPE) = (DNA::IO::TYPE1, DNA::IO::TYPE2, "FILE", "FILE"); 

my(%valid_data) = &DNA::IO::evaluate($config_file, \@FIELD, \@TYPE);  # Evaluate input validity.

my($trial) = $valid_data{"File_name"};

if (!($trial =~ /^$|\s*/)) { print "\n***** DEFINED: $trial. *****\n\n"; }


my($prom_obj) = new DNA::Promoter();	
%{ $prom_obj->{parameters} } = %valid_data;	        # Assign input data to prom_obj.

$prom_obj->get_conditions();

my(%hash) = %{ $prom_obj->{tf_conds} };
foreach my $k (sort keys %hash) { print "$k => $hash{$k}\n"; }
