#-------------------------------------------------------------
# Author: Mirko Palla.
# Date: August 1, 2006.
# For: Task 6 at the Church Lab - Genetics Department, 
# Harvard Medical School.
# 
# Purpose: This program takes a configuration file specifying 
# a TF database file [Finalized_db], and constructs a finalized 
# map with motif shift correction based on PSSMs in Perl.
#
# Input: config_file -> Pre-finalized: [Finalized_db]
# 
# Output: a text file in post-finalized ".db" format
#-------------------------------------------------------------- 

#!/usr/bin/perl

use strict;

use DNA::Promoter;				   	 # Import Promoter package.
use DNA::Fact_list;					 # Import Fact_list package.

#------------------------------ User input handling -----------------------------------

if (@ARGV<1) {
  die "\n--> Error: not correct input!\n--> Usage: perl correct_motif.pl config_file\n\n";
}

my($config_file) = &DNA::IO::validate_cf($ARGV[0]);

#----------------------------- Config. file handling ----------------------------------

my(@FIELD) = (DNA::IO::FIELD1, DNA::IO::FIELD2, "Out_dir"); 
my(@TYPE) = (DNA::IO::TYPE1, DNA::IO::TYPE2); 
             
my(%valid_data) = &DNA::IO::evaluate($config_file, \@FIELD, \@TYPE);  # Evaluate input validity.
   
#-------------------------- Valid IO parameter assignment ----------------------------
          
print "\n--> Process started - correct_motif.pl";	 # Process start.
          
my($out_dir) = $valid_data{Out_dir};			 # Variable assignment. 

#----------------------------- Promoter object creation ------------------------------

my($prom_obj) = new DNA::Promoter();			 # Create Promoter object.
%{ $prom_obj->{parameters} } = %valid_data;		 # Assign input data to prom_obj.

my($prom_list) = $prom_obj->promoter_list();  		 # Get promoter set of yeast. 

#---------------------------- Factor list object creation -----------------------------

my($factors) = new DNA::Fact_list();		         # Create TF list object.
%{ $factors->{parameters} } = %valid_data;	         # Assign input data to Fact_list.

#------------------------------- Motif shift correction -------------------------------

$factors->correct_motif($prom_list);	                 # Create 'motif shift' corrected TF database. 

print "\n--> Location - directory: $out_dir";
print "\n--> Process finished - correct_motif.pl\n\n";
