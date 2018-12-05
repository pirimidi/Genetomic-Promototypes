#-----------------------------------------------------------
# Author: Mirko Palla.
# Date: June 15, 2006.
# For: Task 5 at the Church Lab - Genetics Department, 
# Harvard Medical School.
# 
# Purpose: This program takes a set of TFs (genes) and a 
# configuration file, and constructs a corresponding output 
# in a text file with the calculated hyper-geometric values 
# for all gene pairs (gene1->gene2, gene2->gene1 direction 
# respectively) in Perl.
#
# Input: set of TFs (genes), config_file.
# 
# Output: a text file with hyper-geometric comparison data.
#------------------------------------------------------------ 

#!/usr/bin/perl

use strict;

use DNA::Promoter; 					# Import Promoter package.
use DNA::Fact_list; 					# Import Fact_list package.

#--------------------------------- User input handling --------------------------------

if (@ARGV<1) {
  die "\n--> Error: not correct input!\n--> Usage: perl global_hyper_g_comp.pl config_file\n\n";
}

my($config_file) = &DNA::IO::validate_cf($ARGV[0]);

#-------------------------------- Config. file handling -------------------------------

my(@FIELD) = (DNA::IO::FIELD1, DNA::IO::FIELD2, DNA::IO::FIELD7, "P_cutoff", "Out_dir"); 
my(@TYPE) = (DNA::IO::TYPE1, DNA::IO::TYPE2, DNA::IO::TYPE7, "FLOAT");  
              
my(%valid_data) = &DNA::IO::evaluate($config_file, \@FIELD, \@TYPE);  # Evaluate input validity.        
        
#---------------------------- Valid IO parameter assignment ---------------------------
                                                                               
print "\n--> Process started - global_hg_comp.pl";        # Process start.
                                                                                        
my($base_dir) = $valid_data{Base_dir};                    # Variable assignment.
my($out_dir) = $valid_data{Out_dir};
my($pssm_dir) = $valid_data{PSSMS};
          
#------------------------- TF list & hyper-geometric calculation ----------------------

my($TF_list) = &DNA::Func::TF_list($pssm_dir);		  # Get TF list of yeast.

my($combi) = new DNA::Combi();	 	  		  # Create Combi object.
%{ $combi->{parameters} } = %valid_data;  		  # Assign input data to factors.

my($hyper_hash) = $combi->global_hg_comp($TF_list);

#------------------------------- Finale - result output -------------------------------

$combi->global_hg_result($hyper_hash, 'i'); 

print "\n--> Location - directory: $out_dir";
print "\n--> Process finished - global_hg_comp.pl\n\n";
