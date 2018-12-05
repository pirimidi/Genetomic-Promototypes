#-----------------------------------------------------------
# Author: Mirko Palla.
# Date: August 11, 2005.
# For: Task 1 at the Church Lab - Genetics Department,
# Harvard Medical School.
#
# Purpose: This program takes a gene specification and a 
# configuration file, then it constructs a corresponding 
# output in a textfile with the two sets (1. randomly 
# shuffled, 2. randomly generated) of gene sequence using 
# the sliding window method in Perl.
#
# Input: gene_name, config_file.
#
# Output: a text file with two sets of gene sequences.
#------------------------------------------------------------

#!/usr/bin/perl

use strict;

use DNA::Promoter;
use DNA::Mutant_set;

#------------------------------ User input handling -----------------------------------

my($gene) = $ARGV[0];

if (@ARGV<2) {
  die "\n--> Error: not correct input!\n--> Usage: perl mutagen.pl gene config_file\n\n";
}

my($config_file) = &DNA::IO::validate_cf($ARGV[1]);

#----------------------------- Config. file handling ----------------------------------

my(@FIELD) = (DNA::IO::FIELD1, DNA::IO::FIELD2, DNA::IO::FIELD3, "Out_dir"); 
my(@TYPE) = (DNA::IO::TYPE1, DNA::IO::TYPE2, DNA::IO::TYPE3); 
              
my(%valid_data) = &DNA::IO::evaluate($config_file, \@FIELD, \@TYPE);  # Evaluate input validity.
     
#------------------------- Valid IO parameter assignment ------------------------------
         
print "\n--> Process started - mutagen.pl";		     # Process start.

my($base_dir) = $valid_data{Base_dir};			     # Variable assignment.
my($out_dir) = $valid_data{Out_dir};  

my($conv_file) = $valid_data{Conv_file};

my($prom) = &DNA::Func::gene_to_promoter($gene, $conv_file, $base_dir); # Convert gene to promoter.

#------------------------------ Promoter object creation -------------------------------

my($prom_obj) = new DNA::Promoter($gene, $prom);	     # Create promoter object.
%{ $prom_obj->{parameters} } = %valid_data;                  # Assign input data to prom_obj.

#---------------------------- Promoter sequence retrieving -----------------------------

$prom_obj->promoter_search();			             # Retrieve promoter sequence.

#----------------------------- Mutant_set object creation --------------------------------

my($mutants) = new DNA::Mutant_set($prom_obj);	     	     # Create Mutant_set object.
%{ $mutants->{parameters} } = %valid_data;	    	     # Assign input data to mutants.

#--------------------------  Promoter sequence permutation -----------------------------

$mutants->sliding_window(0); 				     # Generate permuted gene sequences.

#-------------------------- Promoter sequence randomization ----------------------------

$mutants->sliding_window(1); 				     # Randomize gene sequences.
                                                                                        
#--------------------------------- Output handling ------------------------------------

$mutants->output("permuted");				     # IO - output file creation.
$mutants->output("randomized");	     

print "\n--> Location - directory: $out_dir";
print "\n--> Process finished - mutagen.pl\n\n";
