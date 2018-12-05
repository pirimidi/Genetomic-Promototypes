#-----------------------------------------------------------
# Author: Mirko Palla.
# Date: February 9, 2006.
# For: Task 3 at the Church Lab - Genetics Department, 
# Harvard Medical School.
# 
# Purpose: This program takes a transcription factor map 
# directory via configuration file, constructs a database: 
# [TF -> promoter list] and records it in a text file in Perl.
#
# Input: TF map directory [cerevisiae-output] via config_file
# 
# Output: a text file with TF -> promoter set database
#------------------------------------------------------------ 

#!/usr/bin/perl

use strict;
use DNA::Stem;						 # Import Stem package.

#------------------------------ User input handling -----------------------------------

if (@ARGV<1) {
  die "\n--> Error: not correct input!\n--> Usage: perl motif_db_nir.pl config_file\n\n";
}

my($config_file) = &DNA::IO::validate_cf($ARGV[0]);

#----------------------------- Config. file handling ----------------------------------

my(@FIELD) = ("Base_dir", "TF_dir", "Out_dir");              
my(@TYPE) = ("DIR", "DIR"); 
              
my(%valid_data) = &DNA::IO::evaluate($config_file, \@FIELD, \@TYPE);  # Evaluate input validity.
    
#------------------------- Valid IO parameter assignment ------------------------------
          
print "\n--> Process started - motif_db_nir.pl";	     # Process start.

my($map_dir) = $valid_data{TF_dir};			     # Variable assignment.
my($out_dir) = $valid_data{Out_dir};			     

#------------------------------- Stem object creation ---------------------------------

my($stem) = new DNA::Stem();  			     	     # Create Stem object.
%{ $stem->{parameters} } = %valid_data;                      # Assign input data to stem.

#--------------------------------- Database creation ----------------------------------

my($TF_hash) = &DNA::Func::motif_db_nir($map_dir);	     # Create database.

#---------------------------------- Output handling -----------------------------------

$stem->db_output($TF_hash);		     		     # IO - output file creation.

print "\n--> Location - directory: $out_dir";
print "\n--> Process finished - motif_db_nir.pl\n\n";
