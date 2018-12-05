#-----------------------------------------------------------
# Author: Mirko Palla.
# Date: May 16, 2006.
# For: Task 4 at the Church Lab - Genetics Department, 
# Harvard Medical School.
# 
# Purpose: This program takes a configuration file specifying 
# a set of ".gff" files [New_maps], and constructs standarized 
# map formats [Old_maps] in Perl.
#
# Input: config_file -> Ernest: [New_maps]
# 
# Output: a text files in standarized ".gff" format
#------------------------------------------------------------ 

#!/usr/bin/perl

use strict;

use Cwd;
use DNA::Promoter;					 # Import Promoter package.

#------------------------------ User input handling -----------------------------------

if (@ARGV<1) {
  die "\n--> Error: not correct input!\n--> Usage: perl standard_map.pl config_file\n\n";
}

my($config_file) = &DNA::IO::validate_cf($ARGV[0]);

#----------------------------- Config. file handling ----------------------------------

my(@FIELD) = (DNA::IO::FIELD1); 
my(@TYPE) = (DNA::IO::TYPE1); 
             
my(%valid_data) = &DNA::IO::evaluate($config_file, \@FIELD, \@TYPE);  # Evaluate input validity.
   
#------------------------------ Promoter object creation ------------------------------

print "\n--> Process started - standard_map.pl";	     # Process start.

my($prom_obj) = new DNA::Promoter();			     # Create Promoter object.
%{ $prom_obj->{parameters} } = %valid_data;		     # Assign input data to obj.

#------------------------------- Promoter set creation --------------------------------

my($prom_list) = $prom_obj->promoter_list(0);	     	     # Get promoter set of yeast.

# ---------------------------------- Ernest's list ------------------------------------

$prom_obj->standard_map($prom_list);		             # Create gene to promoter map. 

my($out_dir) = cwd();
print "\n--> Location - directory: $out_dir";
print "\n--> Process finished - standard_map.pl\n\n";
