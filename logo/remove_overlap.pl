#--------------------------------------------------------------
# Author: Mirko Palla.
# Date: April 26, 2006.
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
# Extension: Furthermore it creates a PDF file containing all
# motif mutation steps with LOGO display and p-value change.
#
# Input: gene_name, config_file.
# 
# Output: a text file with mutagenized promoter sequences
# and a PDF file with motif mutation statistics.
#--------------------------------------------------------------- 

#!/usr/bin/perl

use strict;

use DNA::Promoter;					# Import Promoter package. 
use DNA::Fact_list; 					# Import Fact_list package.
use DNA::Overlap_set; 					# Import Overlap_set package.

use PDF::API2;						# Import API2 package.
use EXT::Logo;						# Import Logo package.
use POSIX qw(ceil);

#--------------------------------------------------------------------------------------
#                            Overlap mutation set creation 				       
#--------------------------------------------------------------------------------------

#------------------------------ User input handling -----------------------------------

my($gene) = $ARGV[0];

if (@ARGV<2) {
  die "\n--> Error: not correct input!\n--> Usage: perl remove_overlap.pl gene config_file\n\n";
}

my($config_file) = &DNA::IO::validate_cf($ARGV[1]);

#----------------------------- Config. file handling ----------------------------------

my(@FIELD) = (DNA::IO::FIELD1, DNA::IO::FIELD2, DNA::IO::FIELD4, DNA::IO::FIELD7, "Logo_dir", "Out_dir"); 
my(@TYPE) = (DNA::IO::TYPE1, DNA::IO::TYPE2, DNA::IO::TYPE4, DNA::IO::TYPE7, "DIR"); 
              
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
                                            
#----------------------------- Motif overlap assignment -------------------------------

my($overlaps) = new DNA::Overlap_set($factors);         # Create Overlap_set object.
%{ $overlaps->{parameters} } = %valid_data;             # Assign input data to overlaps.

$overlaps->overlap();
$overlaps->remove_overlap();			        # Remove motif overlap pairs.

#------------------------------- Output file creation ---------------------------------

$overlaps->remove_overlap_output();

#--------------------------------------------------------------------------------------
#                                Visualization routine 				       
#--------------------------------------------------------------------------------------

#------------------------------------- Constants ------------------------------------

use constant H_HS => 700;				# Header height start. 

#----------------------------------- Page parameters ---------------------------------

my($visu) = new EXT::Logo($overlaps);		        # Create Logo object.
%{ $visu->{parameters} } = %valid_data;	                # Assign input data to overlaps.

$visu->create_PDF();
$visu->get_font();

$visu->create_header(-1);				# Create title box.
$visu->insert_config();					# Config. parameters.

my(@m_prom) = @{ $visu->{overlaps}->{removed} };	# Removed overlap (Mutant) objects.
my($pages) = ceil((scalar @m_prom)/2);			# Number of pages needed.

my($pos)=0;						# Mutant object position in list.
PG: for (my $i=0; $i<$pages; $i++) { 
  if (!defined($m_prom[$pos])) { last PG; }
  
  $visu->create_header($i);				# Create title box.

  for (my $k=0; $k<2; $k++) {				# Logo parameters.
    if (!defined($m_prom[$pos])) { last PG; }

    my($L_HS) = H_HS-(200+$k*325);			# Logo height start.

    my($pm) = $m_prom[$pos];				# Get A->B || B->A promoter mutation.
    
    my($m_A) = $pm->{overlap_pair}->{motif_A};  	# Untouched motif.
    my($m_B) = $pm->{overlap_pair}->{motif_B};  	# Mutated motif.

    my($image_1) = $m_B->{factor}.".png";		# Motif logo to mutate.
    my($image_2) = $m_A->{factor}.".png";		# Motif logo to keep.

    $visu->insert_logos($L_HS, $image_1, $image_2);     # Image parameters.
    $visu->insert_text($L_HS, $pm, $pos);		# Text parameters.

    $pos++;
  }
}

$visu->{pdf}->save;
$visu->{pdf}->end( );

#---------------------------- End of remove_overlap.pl ------------------------------

print "\n--> $visu->{file} - written."; 
print "\n--> Location - directory: $out_dir";
print "\n--> Process finished - remove_overlap.pl\n\n";
