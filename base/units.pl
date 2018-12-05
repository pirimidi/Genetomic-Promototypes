#-----------------------------------------------------------
# Author: Mirko Palla.
# Date: February 23, 2006.
# For: Task 3 at the Church Lab - Genetics Department, 
# Harvard Medical School.
# 
# Purpose: This program takes a gene specification and a
# configuration file, and constructs a corresponding 
# output in a textfile with the appropriate motif unit
# description holding member TF data (name, location on 
# chromosome, P-value and its motif) in Perl.
#
# Input: gene_name, config_file.
# 
# Output: a text file with motif unit data.
#------------------------------------------------------------ 

#!/usr/bin/perl

use strict;

use DNA::Promoter;					# Import promoter package. 
use DNA::Fact_list;					# Import Fact_list.
use DNA::Unit_set;					# Import Unit_set package.

#------------------------------- User input handling ----------------------------------

my($gene) = $ARGV[0];

if (@ARGV<2) {
  die "\n--> Error: not correct input!\n--> Usage: perl units.pl gene config_file\n\n";
}

my($config_file) = &DNA::IO::validate_cf($ARGV[1]);

#------------------------------ Config. file handling ---------------------------------

my(@FIELD) = (DNA::IO::FIELD1, DNA::IO::FIELD2, DNA::IO::FIELD7, DNA::IO::FIELD8, "Out_dir"); 
my(@TYPE) = (DNA::IO::TYPE1, DNA::IO::TYPE2, DNA::IO::TYPE7, DNA::IO::TYPE8); 
              
my(%valid_data) = &DNA::IO::evaluate($config_file, \@FIELD, \@TYPE);  # Evaluate input validity.
              
#------------------------- Valid IO parameter assignment ------------------------------

print "\n--> Process started - units.pl";		     # Process start.

my($base_dir) = $valid_data{Base_dir};			     # Variable assignment.
my($out_dir) = $valid_data{Out_dir};	

my($conv_file) = $valid_data{Conv_file};
my($mode) = $valid_data{Mode};

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

#------------------------ Unit assignment / output file creation ----------------------

my($units) = new DNA::Unit_set($prom_obj);		      # Create Unit_set object.
%{ $units->{parameters} } = %valid_data;	              # Assign input data to units.

# Mode 1: every copy of a motif is a unit.

if ($mode == 1) { 
  $units->mode_1($factors); 
  $units->mode1_output();
}
  
# Mode 2: all copies of a motif is a unit.
                        
if ($mode == 2) { 
  $units->mode_2($factors);				     # Units, a hash reference.
  $units->mode2_output();
}

# Mode 3: all copies of a motif in distance 'd' from eachother is a unit.
                        
if ($mode == 3) { 
  $units->mode_3($factors); 
  $units->mode3_output();
}

print "\n--> Location - directory: $out_dir";
print "\n--> Process finished - units.pl\n\n";
