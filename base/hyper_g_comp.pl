#-----------------------------------------------------------
# Author: Mirko Palla.
# Date: June 3, 2006.
# For: Task 4 at the Church Lab - Genetics Department, 
# Harvard Medical School.
# 
# Purpose: This program takes two genes and a configuration 
# file, and constructs a corresponding output in a textfile 
# with the calculated hyper-geometric values (gene1->gene2,
# gene2->gene1 direction respectively) in Perl.
#
# Input: gene1_name, gene2_name, config_file.
# 
# Output: a text file with hyper-geometric comparison data.
#------------------------------------------------------------ 

#!/usr/bin/perl

use strict;

use DNA::Promoter; 					# Import Promoter package.
use DNA::Fact_list; 					# Import Fact_list package.

#--------------------------------- User input handling --------------------------------

my(@genes) = ($ARGV[0], $ARGV[1]);

if (@ARGV<3) {
  die "\n--> Error: not correct input!\n--> Usage: perl hyper_g_comp.pl gene1 gene2 config_file\n\n";
}

my($config_file) = &DNA::IO::validate_cf($ARGV[2]);

#-------------------------------- Config. file handling -------------------------------

my(@FIELD) = (DNA::IO::FIELD1, DNA::IO::FIELD2, DNA::IO::FIELD7, "Out_dir"); 
my(@TYPE) = (DNA::IO::TYPE1, DNA::IO::TYPE2, DNA::IO::TYPE7);  
              
my(%valid_data) = &DNA::IO::evaluate($config_file, \@FIELD, \@TYPE);  # Evaluate input validity.        
        
#------------------------- Valid IO parameter assignment ------------------------------
                                                                               
print "\n--> Process started - hyper_g_comp.pl";          # Process start.
                                                                                        
my($base_dir) = $valid_data{Base_dir};                    # Variable assignment.
my($out_dir) = $valid_data{Out_dir};

my($conv_file) = $valid_data{Conv_file};
          
my($factors);					          # Container for Fact_list object.                                                                              
my(@sets);						  # Container for gene TF sets.

my($iter)=0;
while ($iter < @genes) {
  my($query_gene) = $genes[$iter]; 

#-------------------------------- Data set handling -----------------------------------

  $factors = new DNA::Fact_list();	   	  	  # Create Fact_list object.
  $factors->{gene_q} = $query_gene;		  	  # Assign query gene to storage.
  %{ $factors->{parameters} } = %valid_data;	  	  # Assign input data to factors.

  $factors->gene_search();	      		  	  # Get list of promoters bound by given TF.

#----------------------------- Promoter list retrieving -------------------------------

  my(@set) = keys %{ $factors->{final_list} };	  	  # Get list of promoters bound by query gene.

  $sets[$iter] = [ @set ];
  $iter++;
}

#----------------------------- Intersection calculation -------------------------------

my(@a) = @{ $sets[0] };				          # Declare target set a.
my(@b) = @{ $sets[1] };					  # Declare target set b.

my(%union);						  # Union tracker.
my(%isect);  						  # Intersection tracker.

my(@union);						  # Union container.
my(@isect);						  # Intersection container.

foreach my $e (@a, @b) { $union{$e}++ && $isect{$e}++ }
@isect = sort {$a cmp $b } (keys %isect);

if (@a>@b) { @sets = ([@b], [@a]);			  # Handle ordering s.t, for each pair (i,j), i<j. 
             @genes = ($ARGV[1], $ARGV[0]); }

#-------------------------------- Argument assignment ---------------------------------

my($combi_obj) = new DNA::Combi(\@sets, \@isect); 	  # Create Combi object.
$combi_obj->{genes} = \@genes;				  # Assign query gene pair.
%{ $combi_obj->{parameters} } = %valid_data;	  	  # Assign input data to Combi.

$combi_obj->hyper_geometric();	  	                  # Calculate hyper-geometric values.

#-------------------------------- Finale - result output ------------------------------

$combi_obj->hg_result(); 

print "\n--> Location - directory: $out_dir";
print "\n--> Process finished - hyper_g_comp.pl\n\n";
