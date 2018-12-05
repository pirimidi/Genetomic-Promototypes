#-----------------------------------------------------------
# Author: Mirko Palla.
# Date: May 30, 2006.
# For: Task 4 at the Church Lab - Genetics Department, 
# Harvard Medical School.
# 
# Purpose: This program takes a gene specification and a
# configuration file, and constructs output in a textfile 
# with the list of promoters it is bound to in Perl.
#
# Input: gene_name, config_file -> [Maps: *.m_db]
# 
# Output: a text file with bound promoter list.
#------------------------------------------------------------ 

#!/usr/bin/perl

use strict;
use DNA::Fact_list;					     # Import Fact_list package.

#------------------------------ User input handling -----------------------------------

my($query_gene) = $ARGV[0];

if (@ARGV<2) {
  die "\n--> Error: not correct input!\n--> Usage: perl promoter_by_gene.pl gene config_file\n\n";
}

my($config_file) = &DNA::IO::validate_cf($ARGV[1]);

#----------------------------- Config. file handling ----------------------------------

my(@FIELD) = (DNA::IO::FIELD1, DNA::IO::FIELD2, DNA::IO::FIELD7, "Out_dir"); 
my(@TYPE) = (DNA::IO::TYPE1, DNA::IO::TYPE2, DNA::IO::TYPE7); 
             
my(%valid_data) = &DNA::IO::evaluate($config_file, \@FIELD, \@TYPE);  # Evaluate input validity.
   
#------------------------- Valid IO parameter assignment ------------------------------
          
print "\n--> Process started - promoter_by_gene.pl";	     # Process start.

my($out_dir) = $valid_data{Out_dir};			     # Variable assignment.  

#----------------------------- TF list object creation ---------------------------------

my($factors) = new DNA::Fact_list();	                     # Create Fact_list object.
$factors->{gene_q} = $query_gene;			     # Assign query gene to storage.
%{ $factors->{parameters} } = %valid_data;                   # Assign input data to factors.

#--------------------------------- Dataset creation -----------------------------------

$factors->gene_search();	      		             # Get list of promoters bound by given TF.

#--------------------------------- Output handling ------------------------------------

$factors->gene_search_output();	      			     # IO - output file creation.

print "\n--> Location - directory: $out_dir";
print "\n--> Process finished - promoter_by_gene.pl\n\n";
