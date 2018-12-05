#-----------------------------------------------------------
# Author: Mirko Palla.
# Date: June 12, 2006.
# For: Task 5 at the Church Lab - Genetics Department, 
# Harvard Medical School.
# 
# Purpose: This program takes a configuration file specifying 
# a binding affinity <-> condition matrix, and constructs a 
# text file documenting all TFs bound to given promoters with 
# given affinity p-value threshold in Perl.
#
# Input: config_file -> binding affinity matrix
# 
# Output: a text files summarizing TF binding 
#------------------------------------------------------------ 

#!/usr/bin/perl

use strict;
use DNA::Promoter;				     # Import Promoter package.

#------------------------------ User input handling -----------------------------------

if (@ARGV<1) {
  die "\n--> Error: not correct input!\n--> Usage: perl binding_summary.pl config_file\n\n";
}

my($config_file) = &DNA::IO::validate_cf($ARGV[0]);

#----------------------------- Config. file handling ----------------------------------

my(@FIELD) = (DNA::IO::FIELD1, "Affinity_file", "P_cutoff", "Out_dir"); 
my(@TYPE) = (DNA::IO::TYPE1, "FILE"); 
             
my(%valid_data) = &DNA::IO::evaluate($config_file, \@FIELD, \@TYPE);  # Evaluate input validity.
   
#------------------------- Valid IO parameter assignment ------------------------------
          
print "\n--> Process started - binding_summary.pl";  # Process start.
          
my($out_dir) = $valid_data{Out_dir};		     # Variable assignment. 

#------------------------------ Promoter object creation ------------------------------

my($prom_obj) = new DNA::Promoter();		     # Create a promoter object.
%{ $prom_obj->{parameters} } = %valid_data;	     # Assign input data to prom_obj.

#----------------------------- Affinitiy matrix creation ------------------------------

my($tf_conds, $proms) = $prom_obj->get_matrix();    # Get binding affinity matrix of yeast.

#----------------------------- Promoter list collection -------------------------------

my($TF_hash) = $prom_obj->binding_summary($tf_conds, $proms);  # Create file summarizing TF binding.

#--------------------------------- Output handling ------------------------------------

$prom_obj->binding_summary_output($TF_hash);	    # Create gene to promoter map.

print "\n--> Location - directory: $out_dir";
print "\n--> Process finished - binding_summary.pl\n\n";
