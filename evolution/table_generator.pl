#--------------------------------------------------------------
# Author: Mirko Palla.
# Date: July 24, 2006.
# For: Task 6 at the Church Lab - Genetics Department, 
# Harvard Medical School.
# 
# Purpose: This program takes a list of S. cerevisiae genes 
# (promoters) and a list of orthologous promoters from the 
# orthogroup database [cer_orthogroups_i.txt] both specified 
# in the configuration file, and constructs a set of text files 
# containing a table of TF comparision of various yeast genomes
# for each gene as output in Perl.
#
# Input: config_file (list of genes, orthogroup database).
# 
# Output: set of tables of transcription factor comparision.
#--------------------------------------------------------------- 

#!/usr/bin/perl

use strict;
use DNA::Promoter_set;					# Import Promoter package. 

#------------------------------ User input handling -----------------------------------

if (@ARGV<1) {
  die "\n--> Error: not correct input!\n--> Usage: perl table_generator.pl config_file\n\n";
}

my($config_file) = &DNA::IO::validate_cf($ARGV[0]);

#----------------------------- Config. file handling ----------------------------------

my(@FIELD) = (DNA::IO::FIELD1, DNA::IO::FIELD2, DNA::IO::FIELD6, DNA::IO::FIELD7, "Out_dir"); 
my(@TYPE) = (DNA::IO::TYPE1, DNA::IO::TYPE2, DNA::IO::TYPE6, DNA::IO::TYPE7); 

my(%valid_data) = &DNA::IO::evaluate($config_file, \@FIELD, \@TYPE);  # Evaluate input validity.

#------------------------- Valid IO parameter assignment ------------------------------

print "\n--> Process started - orthogroup_table.pl";	              # Process start.

my($base_dir) = $valid_data{Base_dir};			              # Variable assignment.
my($out_dir) = $valid_data{Out_dir};
			     			     
my($conv_file) = $valid_data{Conv_file};
my($prom_file) = $valid_data{Prom_file};

#------------------------------- Promoter set creation --------------------------------

chdir $base_dir;
open(INFO, $prom_file) || die "--> Cannot open promoter file - $prom_file\n";          

my(%no_hits);						     # Storage hash for no hit promoters.
 
PROM: while (<INFO>) {
  my(@prom) = split(/\s/, $_);          		     # Put line member strings into an array, @prom. 
  my($prom) = substr($prom[0], 1);			     # Get promoter name.

  my($gene) = &DNA::Func::promoter_to_gene($prom, $conv_file, $base_dir); # Convert gene to promoter.

  my($prom_set) = new DNA::Promoter_set($gene, $prom);	     # Create Promoter set object.
  %{ $prom_set->{parameters} } = %valid_data;		     # Assign input data to the set.

  #------------------------------- Retrieve orthogroups -------------------------------

  my($resp) = $prom_set->get_orthogroups();

  if ($resp!=0) { 	 	                             # If promoter not in analysis or,  
    push @{ $no_hits{$resp} }, $prom_set;                    # if there are no other orthologs, skip promoter. 
    next PROM;
  }

  #---------------------------- Construct comparison table ----------------------------

  $prom_set->get_comp_table();
  $prom_set->create_table();
 
  $prom_set->list_table();				     # If this option can be commented out. 
                                                             
  #--------------------------------- Output handling ----------------------------------

  $prom_set->ortho_table_output(); 			     # IO - output file creation.
}
close PROM;

if (defined %no_hits) { &DNA::IO::excluded_proms_IO($out_dir, \%no_hits); }

print "\n--> Location - directory: $out_dir";
print "\n--> Process finished - orthogroup_table.pl\n\n";
