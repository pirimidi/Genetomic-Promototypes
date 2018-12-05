#-----------------------------------------------------------
# Author: Mirko Palla.
# Date: June 21, 2006.
# For: Task 5 at the Church Lab - Genetics Department, 
# Harvard Medical School.
# 
# Purpose: This program takes a set of PSSM matrices and 
# constructs a randomized kmer library (hash) - containing 
# ten different kmers - of given matrix lengths in Perl.
#
# Input: PSSM file directory
# 
# Output: randomized kmer library of each PSSM length
#------------------------------------------------------------ 

#!/usr/bin/perl

use strict;

use DNA::Func;					             # Import Func package.
use DNA::Mutant_set;					     # Import Mutant_set package.

#------------------------------- User input handling ---------------------------------

if (@ARGV<1) {
  die "\n--> Error: not correct input!\n--> Usage: perl random_kmer.pl config_file\n\n";
}

my($config_file) = &DNA::IO::validate_cf($ARGV[0]);

#------------------------------ Config. file handling --------------------------------

my(@FIELD) = ("Base_dir", "PSSMS", "Matrix", "Out_dir"); 
my(@TYPE) = ("DIR", "DIR", "BINARY"); 
              
my(%valid_data) = &DNA::IO::evaluate($config_file, \@FIELD, \@TYPE);  # Evaluate input validity.

#-------------------------- Valid IO parameter assignment ----------------------------

print "\n--> Process started - random_kmer.pl";		     # Process start.

my($base_dir) = $valid_data{Base_dir};			     # Variable assignment.
my($out_dir) = $valid_data{Out_dir};	
my($pssm_dir) = $valid_data{PSSMS};

#---------------------------- Mutant_set object creation -----------------------------

my($mutants) = new DNA::Mutant_set();			     # Create Mutant_set object.
%{ $mutants->{parameters} } = %valid_data;		     # Assign input data to mutants.

#--------------------------- Random kmer library creation ----------------------------

$mutants->get_r_kmer_library();

#------------------------------- Output file creation --------------------------------

$mutants->random_kmer_output();

print "\n--> Location - directory: $out_dir";
print "\n--> Process finished - random_kmer.pl\n\n";
