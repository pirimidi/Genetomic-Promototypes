#-----------------------------------------------------------
# Author: Mirko Palla.
# Date: May 26, 2006.
# For: Task 4 at the Church Lab - Genetics Department, 
# Harvard Medical School.
# 
# Purpose: This program takes a set of transcription factor 
# maps via a configuration file, constructs a database:
# [promoter -> bound TF set] and records it in a text file 
# in Perl.
#
# Input: TF map name [maps] via config_file
# 
# Output: a text file with promoter -> bound TF set database
#------------------------------------------------------------ 

#!/usr/bin/perl

use strict;
use DNA::Fact_list;					     # Import Fact_list package.

#------------------------------ User input handling -----------------------------------

if (@ARGV<1) {
  die "\n--> Error: not correct input!\n--> Usage: perl promoter_db.pl config_file\n\n";
}

my($config_file) = &DNA::IO::validate_cf($ARGV[0]);

#----------------------------- Config. file handling ----------------------------------

my(@FIELD) = (DNA::IO::FIELD1, "Out_dir"); 
my(@TYPE) = (DNA::IO::TYPE1); 
             
my(%valid_data) = &DNA::IO::evaluate($config_file, \@FIELD, \@TYPE);  # Evaluate input validity.
   
#------------------------- Valid IO parameter assignment ------------------------------
          
print "\n--> Process started - promoter_db.pl";		     # Process start.

my($out_dir) = $valid_data{Out_dir};			     # Variable assignment.  

#----------------------------- TF list object creation ---------------------------------

my($factors) = new DNA::Fact_list();	                     # Create Fact_list object.
%{ $factors->{parameters} } = %valid_data;                   # Assign input data to factors.

#--------------------------------- Dataset creation -----------------------------------

$factors->promoter_db();				     # Create database.

print "\n--> Location - directory: $out_dir";
print "\n--> Process finished - promoter_db.pl\n\n";
