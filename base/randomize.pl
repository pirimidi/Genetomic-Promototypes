#-----------------------------------------------------------
# Author: Mirko Palla.
# Date: February 20, 2006.
# For: Task 3 at the Church Lab - Genetics Department,
# Harvard Medical School.
#
# Purpose: This program takes a gene specification and a
# configuration file, then it constructs a corresponding 
# output in textfile with a set of gene sequences generated 
# by random P-value proof kmer substitution using the 
# sliding window method in Perl.
#
# Input: gene_name, config_file.
#
# Output: a text file with a set of randomized gene 
# sequences.
#------------------------------------------------------------

#!/usr/bin/perl

use strict;

use DNA::Promoter;
use DNA::Mutant_set;

#------------------------------ User input handling -----------------------------------

my($gene) = $ARGV[0];

if (@ARGV<2) {
  die "\n--> Error: not correct input!\n--> Usage: perl randomize.pl gene config_file\n\n";
}

my($config_file) = &DNA::IO::validate_cf($ARGV[1]);

#----------------------------- Config. file handling ----------------------------------

my(@FIELD) = (DNA::IO::FIELD1, DNA::IO::FIELD2, DNA::IO::FIELD3, "Out_dir"); 
my(@TYPE) = (DNA::IO::TYPE1, DNA::IO::TYPE2, DNA::IO::TYPE3); 
              
my(%valid_data) = &DNA::IO::evaluate($config_file, \@FIELD, \@TYPE);  # Evaluate input validity.
              
#------------------------- Valid IO parameter assignment ------------------------------

print "\n--> Process started - randomize.pl";		     # Process start.

my($base_dir) = $valid_data{Base_dir};			     # Variable assignment.
my($out_dir) = $valid_data{Out_dir};  

my($conv_file) = $valid_data{Conv_file};

my($prom) = &DNA::Func::gene_to_promoter($gene, $conv_file, $base_dir); # Convert gene to promoter.

#----------------------------- Promoter object creation -------------------------------

my($prom_obj) = new DNA::Promoter($gene, $prom);	     # Create promoter object.
%{ $prom_obj->{parameters} } = %valid_data;                  # Assign input data to prom_obj.

#--------------------------- Promoter sequence retrieving -----------------------------

$prom_obj->promoter_search();   		             # Retrieve promoter sequence.

#---------------------------- Mutant_set object creation ------------------------------

my($mutants) = new DNA::Mutant_set($prom_obj);		     # Create Mutant_set object.
%{ $mutants->{parameters} } = %valid_data;	   	     # Assign input data to mutants.

#------------------------- Promoter sequence randomization -----------------------------

$mutants->sliding_window(1); 				     # Randomize gene sequences.

#--------------------------------- Output handling -------------------------------------

$mutants->output("randomized");			             # IO - output file creation.

print "\n--> Location - directory: $out_dir";
print "\n--> Process finished - randomize.pl\n\n";
