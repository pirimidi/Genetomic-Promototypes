#-----------------------------------------------------------
# Author: Mirko Palla.
# Date: August 4, 2006.
# For: Task 6 at the Church Lab - Genetics Department, 
# Harvard Medical School.
# 
# Purpose: This program takes a transcription factor map
# and a configuration file, creates a list of transcription
# factors (Fact_list) and promoter objects (Promoter) for each 
# promoter for S. cerevisiae. Then calculates the relative bp 
# distance, orientation and order for each TF pair in all pro-
# moter regions. Then it creates two ".txt" files containing 
# all distance and orientation values respectively. Finally,
# using these two text files (list of integers) it executes some 
# MATLAB functions (hist) to display the numeric distribution via
# a histogram with bin-size of 10 and orientation types (++, --, 
# -+, +-) respectively in Perl.
# 
# Input: TF map file [Finalized_db], config_file
# 
# Output: two histograms displaying relative distance and 
# orientation for all TF pairs in the yeast genome
#------------------------------------------------------------ 

#!/usr/bin/perl

use strict;

use DNA::Promoter_set;					# Import Promoter_set package.

#------------------------------ User input handling -----------------------------------

if (@ARGV<1) {
  die "\n--> Error: not correct input!\n--> Usage: perl relative_distance.pl config_file\n\n";
}

my($config_file) = &DNA::IO::validate_cf($ARGV[0]);

#----------------------------- Config. file handling ----------------------------------

my(@FIELD) = (DNA::IO::FIELD1, DNA::IO::FIELD7, "Finalized_db", "Out_dir"); 
my(@TYPE) = (DNA::IO::TYPE1, DNA::IO::TYPE7, "FILE"); 

my(%valid_data) = &DNA::IO::evaluate($config_file, \@FIELD, \@TYPE);  # Evaluate input validity.

#------------------------- Valid IO parameter assignment ------------------------------

print "\n--> Process started - relative_distance.pl";        # Process start.

my($out_dir) = $valid_data{Out_dir};			     

#------------------------------- Promoter set creation --------------------------------

my($prom_set) = new DNA::Promoter_set();		     # Create Promoter set object.
%{ $prom_set->{parameters} } = %valid_data;		     # Assign input data to the set. 

#------------------------------- TF-Promoter set creation -----------------------------

$prom_set->all_tf_list_promoter();			     # Collect all TF list (plus promoter object).
$prom_set->relative_distance();				     # Calculate relative distance (orientation) statistics.

#----------------------------------- Output handling ----------------------------------

$prom_set->statistics_output();				     # IO - histogram creation.

print "\n--> Location - directory: $out_dir";
print "\n--> Process finished - relative_distance.pl\n\n";
