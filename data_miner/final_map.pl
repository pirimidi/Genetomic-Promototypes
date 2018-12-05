#-----------------------------------------------------------
# Author: Mirko Palla.
# Date: June 8, 2006.
# For: Task 5 at the Church Lab - Genetics Department, 
# Harvard Medical School.
# 
# Purpose: This program takes a configuration file specifying 
# a set of ".p_db" files [Maps], and constructs finalized 
# map formats in Perl.
#
# Input: config_file -> Standard: [Maps]
# 
# Output: a text files in finalized ".db" format
#------------------------------------------------------------ 

#!/usr/bin/perl

use strict;
use DNA::Fact_list;					 # Import Fact_list package.

#------------------------------ User input handling -----------------------------------

if (@ARGV<1) {
  die "\n--> Error: not correct input!\n--> Usage: perl final_map.pl config_file\n\n";
}

my($config_file) = &DNA::IO::validate_cf($ARGV[0]);

#----------------------------- Config. file handling ----------------------------------

my(@FIELD) = (DNA::IO::FIELD1, "Affinity_file", "Data_bases", "Out_dir"); 
my(@TYPE) = (DNA::IO::TYPE1, "FILE"); 
             
my(%valid_data) = &DNA::IO::evaluate($config_file, \@FIELD, \@TYPE);  # Evaluate input validity.
   
#------------------------- Valid IO parameter assignment ------------------------------
          
print "\n--> Process started - final_map.pl";	     # Process start.
          
my($out_dir) = $valid_data{Out_dir};			     # Variable assignment. 

#------------------------------ Promoter object creation ------------------------------

my($factors) = new DNA::Fact_list();		     # Create TF list object.
%{ $factors->{parameters} } = %valid_data;	     # Assign input data to Fact_list.

#----------------------------- Affinitiy matrix creation --------------------------------

my($tf_conds, $proms) = $factors->get_matrix();	     # Get binding affinity matrix of yeast.

#----------------------------------- Ernest's list ------------------------------------

$factors->final_map($tf_conds, $proms);	             # Create finalized gene to promoter map. 

print "\n--> Location - directory: $out_dir";
print "\n--> Process finished - final_map.pl\n\n";
