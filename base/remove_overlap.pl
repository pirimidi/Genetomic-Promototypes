#--------------------------------------------------------------
# Author: Mirko Palla.
# Date: March 27, 2006.
# For: Task 3 at the Church Lab - Genetics Department, 
# Harvard Medical School.
# 
# Purpose: This program takes a gene specification and a
# configuration file, and constructs a corresponding 
# output in a textfile with a set of mutagenized promoter 
# sequences, where for every motif overlap pair: motif A is 
# "knocked out" using PSSM values, while motif B is conserved 
# [as it was observed originally] and vice and versa in Perl.
#
# Input: gene_name, config_file.
# 
# Output: a text file pair with mutagenized promoter sequences.
#--------------------------------------------------------------- 

#!/usr/bin/perl

use strict;

use DNA::Promoter;					# Import Promoter package. 
use DNA::Fact_list; 					# Import Fact_list package.
use DNA::Overlap_set; 					# Import Overlap_set package.

#------------------------------ User input handling -----------------------------------

my($gene) = $ARGV[0];

if (@ARGV<2) {
  die "\n--> Error: not correct input!\n--> Usage: perl remove_overlap.pl gene config_file\n\n";
}

my($config_file) = &DNA::IO::validate_cf($ARGV[1]);

#----------------------------- Config. file handling ----------------------------------


my(@FIELD) = (DNA::IO::FIELD1, DNA::IO::FIELD2, DNA::IO::FIELD4, DNA::IO::FIELD7, "Out_dir"); 
my(@TYPE) = (DNA::IO::TYPE1, DNA::IO::TYPE2, DNA::IO::TYPE4, DNA::IO::TYPE7); 
              
my(%valid_data) = &DNA::IO::evaluate($config_file, \@FIELD, \@TYPE);  # Evaluate input validity.
              
#------------------------- Valid IO parameter assignment ------------------------------

print "\n--> Process started - remove_overlap.pl";	     # Process start.

my($base_dir) = $valid_data{Base_dir};			     # Variable assignment.
my($out_dir) = $valid_data{Out_dir};  

my($conv_file) = $valid_data{Conv_file};

my($prom) = &DNA::Func::gene_to_promoter($gene, $conv_file, $base_dir); # Convert gene to promoter.

#----------------------------- Promoter object creation -------------------------------

my($prom_obj) = new DNA::Promoter($gene, $prom);	# Create promoter object.
%{ $prom_obj->{parameters} } = %valid_data;	        # Assign input data to prom_obj.

#--------------------------- Promoter sequence retrieving -----------------------------
            
$prom_obj->promoter_search();		                # Retrieve promoter sequence.

#-------------------------------- Data set handling -----------------------------------

my($factors) = new DNA::Fact_list($prom_obj);		# Create Fact_list object.
%{ $factors->{parameters} } = %valid_data;	        # Assign input data to factors.

$factors->TF_search();					# Retrieve all TFs bound to given promoter.
$factors->TF_search_output();	
                                            
#----------------------------- Motif overlap assignment -------------------------------

my($overlaps) = new DNA::Overlap_set($factors);         # Create Overlap_set object.
%{ $overlaps->{parameters} } = %valid_data;             # Assign input data to overlaps.

$overlaps->overlap();
$overlaps->remove_overlap();			        # Remove motif overlap pairs.

#------------------------------- Output file creation ---------------------------------

$overlaps->remove_overlap_output();

print "\n--> Location - directory: $out_dir";
print "\n--> Process finished - remove_overlap.pl\n\n";
