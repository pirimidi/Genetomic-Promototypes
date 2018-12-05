#-----------------------------------------------------------
# Author: Mirko Palla.
# Date: March 21, 2006.
# For: Task 3 at the Church Lab - Genetics Department, 
# Harvard Medical School.
# 
# Purpose: This program takes a gene specification and a
# configuration file, and constructs a corresponding 
# output in a textfile with the appropriate motif overlap
# pair description holding member TF data (name, location 
# on chromosome, P-value and its motif) in Perl.
#
# Input: gene_name, config_file.
# 
# Output: a text file with motif overlap pair data.
#------------------------------------------------------------ 

#!/usr/bin/perl

use strict;

use DNA::Promoter;					# Import promoter package. 
use DNA::Fact_list;					# Import Fact_list.
use DNA::Overlap_set;					# Import Overlap_set package.

#------------------------------- User input handling ----------------------------------

my($gene) = $ARGV[0];

if (@ARGV<2) {
  die "\n--> Error: not correct input!\n--> Usage: perl overlap.pl gene config_file\n\n";
}

my($config_file) = &DNA::IO::validate_cf($ARGV[1]);

#------------------------------ Config. file handling ---------------------------------

my(@FIELD) = (DNA::IO::FIELD1, DNA::IO::FIELD2, DNA::IO::FIELD7, "Overlap_gap", "Out_dir"); 
my(@TYPE) = (DNA::IO::TYPE1, DNA::IO::TYPE2, DNA::IO::TYPE7, "INT"); 
              
my(%valid_data) = &DNA::IO::evaluate($config_file, \@FIELD, \@TYPE);  # Evaluate input validity.
              
#------------------------- Valid IO parameter assignment ------------------------------

print "\n--> Process started - overlap.pl";		     # Process start.

my($base_dir) = $valid_data{Base_dir};			     # Variable assignment.
my($out_dir) = $valid_data{Out_dir};  

my($conv_file) = $valid_data{Conv_file};

my($prom) = &DNA::Func::gene_to_promoter($gene, $conv_file, $base_dir); # Convert gene to promoter.

#------------------------------ Promoter object creation ------------------------------

my($prom_obj) = new DNA::Promoter($gene, $prom);	# Create promoter object.
%{ $prom_obj->{parameters} } = %valid_data;	        # Assign input data to prom_obj.

#--------------------------- Promoter sequence retrieving -----------------------------

$prom_obj->promoter_search();		                # Retrieve promoter sequence.

#-------------------------------- Data set handling -----------------------------------

my($factors) = new DNA::Fact_list($prom_obj);		# Create Fact_list object.
%{ $factors->{parameters} } = %valid_data;	        # Assign input data to factors.

$factors->TF_search();					# Retrieve all TFs bound to given promoter.

#--------------------- Motif overlap assignment / output file creation ----------------------

my($overlaps) = new DNA::Overlap_set($factors);         # Create Overlap_set object.
%{ $overlaps->{parameters} } = %valid_data;	        # Assign input data to overlaps.

$overlaps->overlap(); 
$overlaps->overlap_output();

print "\n--> Location - directory: $out_dir";
print "\n--> Process finished - overlap.pl\n\n";
