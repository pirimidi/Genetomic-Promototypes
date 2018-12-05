#-----------------------------------------------------------
# Author: Mirko Palla.
# Date: March 15, 2006.
# For: Task 3 at the Church Lab - Genetics Department, 
# Harvard Medical School.
# 
# Purpose: This program takes a transcription factor map
# and a configuration file, constructs a database:
# [TF -> promoter list] and records it in a text file in Perl.
#
# Input: TF map file [ernest], config_file
# 
# Output: a text file with TF -> promoter set database
#------------------------------------------------------------ 

#!/usr/bin/perl

use strict;

use DNA::Promoter;					# Import Promoter package.
use DNA::Fact_list;					# Import Fact_list package.

#------------------------------ User input handling -----------------------------------

if (@ARGV<2) {
  die "\n--> Error: not correct input!\n--> Usage: perl motif_db_ernest.pl map_file config_file\n\n";
}

my($map_file) = $ARGV[0];  
my($config_file) = &DNA::IO::validate_cf($ARGV[1]);

#----------------------------- Config. file handling ----------------------------------

my(@FIELD) = (DNA::IO::FIELD1, "Matrix", "PSSMS", "Out_dir"); 
my(@TYPE) = (DNA::IO::TYPE1, "BINARY", "DIR"); 

my(%valid_data) = &DNA::IO::evaluate($config_file, \@FIELD, \@TYPE);  # Evaluate input validity.

#------------------------- Valid IO parameter assignment ------------------------------

print "\n--> Process started - motif_db_ernest.pl";		              # Process start.

my($out_dir) = $valid_data{Out_dir};			     

#------------------------------ Promoter object creation ------------------------------

my($prom_obj) = new DNA::Promoter();			     # Create Promoter object.
%{ $prom_obj->{parameters} } = %valid_data;		     # Assign input data to prom_obj.

my($prom_list) = $prom_obj->promoter_list();  		     # Get promoter set of yeast. 

#------------------------------- TF-Promoter set creation -----------------------------

my($factors) = new DNA::Fact_list();			     # Create Fact_list object.
%{ $factors->{parameters} } = %valid_data;                   # Assign input data to factors.

my(@e_factors) = $factors->all_tf($map_file);	             # Get all TFs from Ernest's list. 

#----------------------------------- Database merging ---------------------------------

my($TF_hash) = $factors->motif_db_ernest($prom_list, @e_factors); # Create database. 

#----------------------------------- Output handling ----------------------------------

$factors->db_output($TF_hash, $map_file);   		     # IO - output file creation.

print "\n--> Location - directory: $out_dir";
print "\n--> Process finished - motif_db_ernest.pl\n\n";
