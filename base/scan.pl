#-----------------------------------------------------------
# Author: Mirko Palla.
# Date: February 23, 2006.
# For: Task 3 at the Church Lab - Genetics Department,
# Harvard Medical School
#
# Purpose: This program takes a gene specification and a 
# configuration file, then it scans the corresponding
# promoter sequence for kmer occurences, and finally constructs 
# an output in a textfile with the co-occuring permuted 
# kmer sequences using the sliding window method in Perl.
#
# Input: gene_name, config_file.
#
# Output: a text file with a set of gene sequences 
# containing co-occuring kmer mutations.
#------------------------------------------------------------

#!/usr/bin/perl

use strict;

use DNA::Promoter;					# Import promoter package.
use DNA::Mutant_set;					# Import Mutant_set package.

#------------------------------ User input handling -----------------------------------

my($gene) = $ARGV[0];			

if (@ARGV<2) {
  die "\n--> Error: not correct input!\n--> Usage: perl scan.pl gene config_file\n\n";
}

my($config_file) = &DNA::IO::validate_cf($ARGV[1]);

#----------------------------- Config. file handling ----------------------------------

my(@FIELD) = (DNA::IO::FIELD1, "Kmer_l", "BP_diff", "Out_dir"); 
my(@TYPE) = (DNA::IO::TYPE1, "KMER_L", "INT"); 
              
my(%valid_data) = &DNA::IO::evaluate($config_file, \@FIELD, \@TYPE);  # Evaluate input validity.

#------------------------- Valid IO parameter assignment ------------------------------

print "\n--> Process started - scan.pl";		     # Process start.

my($base_dir) = $valid_data{Base_dir};			     # Variable assignment.
my($out_dir) = $valid_data{Out_dir};  

my($conv_file) = $valid_data{Conv_file};

my($prom) = &DNA::Func::gene_to_promoter($gene, $conv_file, $base_dir); # Convert gene to promoter.

#------------------------------ Promoter object creation ------------------------------

my($prom_obj) = new DNA::Promoter($gene, $prom);	     # Create promoter object.
%{ $prom_obj->{parameters} } = %valid_data;    	             # Assign input data to prom_obj.

#---------------------------- Promoter sequence retrieving ----------------------------

$prom_obj->promoter_search();		             	     # Retrieve promoter sequence.

#------------------------------- Mutant object creation -------------------------------

my($mutants) = new DNA::Mutant_set($prom_obj);		     # Create Mutant_set object.
%{ $mutants->{parameters} } = %valid_data;		     # Assign input data to obj.

#--------------------------------- Scanning routine -----------------------------------

my(%seq_count) = $mutants->scan();		             # Document occurences in a hash.
$mutants->perm_scan(\%seq_count, 1);		             # Generate permuted gene sequence.

#--------------------------------- Output handling ------------------------------------

$mutants->output("scanned");			             # IO - output file creation.

print "\n--> Location - directory: $out_dir";
print "\n--> Process finished - scan.pl\n\n";
