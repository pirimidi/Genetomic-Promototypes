#-------------------------------------------------------------
# Author: Mirko Palla.
# Date: August 11, 2006.
# For: Task 7 at the Church Lab - Genetics Department,
# Harvard Medical School.
#
# Purpose: This program takes a gene specification for 'cer',
# and a configuration file, then constructs a visual represen-
# tation for the TF binding sites on the given promoter and its 
# orhologs in other yeast species. Finally, it outputs a text-
# file with the appropriate transcription factor description 
# (name, location on chromosome, P-value and its motif) in Perl. 
#
# Input: gene_name, config_file.
#
# Output: a GUI for TF binding site viewing for all orthgologs, 
# and a textfile with binded transcription factor data.
#--------------------------------------------------------------

#!/usr/local/bin/perl   
                                                                                        
use strict;

#------------------------------ Library / module import -------------------------------

use Tk; 						# For GUI display.
use Tk::ROText;						# For interactive window.

use DNA::Promoter_set;					# Import Promoter package. 
use DNA::Visu;						# Visualization class.
                              
#--------------------------------------------------------------------------------------
#                                  TF-map creation 				       
#--------------------------------------------------------------------------------------

#------------------------------ User input handling -----------------------------------

my($gene) = $ARGV[0];

if (@ARGV<2) {
  die "\n--> Error: not correct input!\n--> Usage: perl orthogroup_visualizer.pl gene config_file\n\n";
}

my($config_file) = &DNA::IO::validate_cf($ARGV[1]);

#----------------------------- Config. file handling ----------------------------------

my(@FIELD) = (DNA::IO::FIELD1, DNA::IO::FIELD2, DNA::IO::FIELD5, DNA::IO::FIELD6, DNA::IO::FIELD7, "Out_dir"); 
my(@TYPE) = (DNA::IO::TYPE1, DNA::IO::TYPE2, DNA::IO::TYPE5, DNA::IO::TYPE6, DNA::IO::TYPE7); 

my(%valid_data) = &DNA::IO::evaluate($config_file, \@FIELD, \@TYPE);  # Evaluate input validity.

#------------------------- Valid IO parameter assignment ------------------------------

print "\n--> Process started - orthogroup_visualizer.pl";             # Process start.

my($base_dir) = $valid_data{Base_dir};			              # Variable assignment.
my($out_dir) = $valid_data{Out_dir};
			     			     
my($conv_file) = $valid_data{Conv_file};

my($prom) = &DNA::Func::gene_to_promoter($gene, $conv_file, $base_dir); # Convert gene to promoter.

#------------------------------- Promoter set creation --------------------------------

my($prom_set) = new DNA::Promoter_set($gene, $prom);	     # Create Promoter set object.
%{ $prom_set->{parameters} } = %valid_data;		     # Assign input data to the set.

#------------------------------- Retrieve orthogroups ---------------------------------

$prom_set->get_orthogroups();

#---------------------------- Construct comparison table ------------------------------

$prom_set->get_comp_table();
$prom_set->get_ortholog_data();
                                                             
#--------------------------------- Output handling ------------------------------------

$prom_set->list_ortho_tfs();				     # IO - output file creation.

#--------------------------------------------------------------------------------------
#                               Visualization routine 				       
#--------------------------------------------------------------------------------------

#--------------------------- TF color / rows calculation  -----------------------------

my($f_rows) = $prom_set->TF_color_rows();			 
#$visu->TF_align($f_rows);		        	       # Alignment routine.

=pod
my($visu) = new DNA::Visu();			  	       # Create Visu object.
$visu->{gene} = $gene;				               # Assign gene to visu.
%{ $visu->{parameters} } = %valid_data;                        # Assign input data to visu.			 

#----------------------- Promoter sequence / TF row insertion --------------------------

$visu->window_setup();				               # Main window / canvas setup.

my($main) = $visu->{window};
my($canvas) = $visu->{canvas};

my($tb) = $visu->promoter_TF_row_insert(5, $canvas);     # Return ROText field.	     			 

#--------------------------- Motif coloring / TF insertion ----------------------------

#my($f_tb) = $visu->motif_color_TF_insert($tb);

#--------------------------------- Canvas update --------------------------------------

$canvas->createWindow(qw/20 20/, -window => $tb, -anchor => 'nw');
$canvas->update();

#--------------------------------- Color palette --------------------------------------

$tb->tag(qw/configure c1 -foreground/ => "#FF6600");
$tb->tag(qw/configure c2 -foreground/ => "#00CC00");
$tb->tag(qw/configure c3 -foreground/ => "#CC9933");
$tb->tag(qw/configure c4 -foreground/ => "#CC33CC");
$tb->tag(qw/configure c5 -foreground/ => "#FF0000");
$tb->tag(qw/configure c6 -foreground/ => "#0033FF");
$tb->tag(qw/configure c7 -foreground/ => "#CC0066");
$tb->tag(qw/configure c8 -foreground/ => "#A52A2A");
$tb->tag(qw/configure c9 -foreground/ => "#C71585");
$tb->tag(qw/configure c10 -foreground/ => "#5F9EA0");
$tb->tag(qw/configure c11 -foreground/ => "#663300");
$tb->tag(qw/configure c12 -foreground/ => "#FF7F50");
$tb->tag(qw/configure c13 -foreground/ => "#6495ED");
$tb->tag(qw/configure c14 -foreground/ => "#008B8B");
$tb->tag(qw/configure c15 -foreground/ => "#B8860B");
$tb->tag(qw/configure c16 -foreground/ => "#006400");
$tb->tag(qw/configure c17 -foreground/ => "#FF8C00");
$tb->tag(qw/configure c18 -foreground/ => "#E99676");
$tb->tag(qw/configure c19 -foreground/ => "#483D8B");
$tb->tag(qw/configure c20 -foreground/ => "#DC143C");

#---------------------------- End of orthogroup_visualizer.pl -------------------------

=cut

print "\n--> Location - directory: $out_dir";
print "\n--> Process finished - orthogroup_visualizer.pl\n\n";

MainLoop;
