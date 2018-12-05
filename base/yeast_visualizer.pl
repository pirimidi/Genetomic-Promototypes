#------------------------------------------------------------
# Author: Mirko Palla.
# Date: August 8, 2006.
# For: Task 3 at the Church Lab - Genetics Department,
# Harvard Medical School.
#
# Purpose: This program takes an promoter name (as listed in
# promoter sequence files), an organism name and a configura-
# tion file, then constructs a visual representation for the
# TF binding sites on the given yeast promoter. Finally, it 
# outputs a textfile with the appropriate transcription factor 
# description (name, location on chromosome, P-value and its 
# motif) in Perl. 
#
# Input: gene_name, org_name, config_file.
#
# Output: a GUI for TF binding site viewing, and a text 
# file with binded transcription factor data.
#-------------------------------------------------------------

#!/usr/local/bin/perl   
                                                                                        
use strict;

#------------------------------ Library / module import -------------------------------

use Tk; 						# For GUI display.
use Tk::ROText;						# For interactive window.

use DNA::Promoter;					# Promoter class.
use DNA::Fact_list;					# TF list object class.
use DNA::Visu;						# Visualization class.
                              
#--------------------------------------------------------------------------------------
#                                  TF-map creation 				       
#--------------------------------------------------------------------------------------

#------------------------------- User input handling ----------------------------------

my($prom) = $ARGV[0];
my($org) = $ARGV[1];

if (@ARGV!=3) {
  die "\n--> Error: not correct input!\n--> Usage: perl yeast_visualizer.pl prom_name org_name config_file\n\n";
}

my($config_file) = &DNA::IO::validate_cf($ARGV[2]);

#----------------------------- Config. file handling ----------------------------------

my(@FIELD) = (DNA::IO::FIELD5, "Base_dir", "Other_dir", "Promoters", "File_name", "Out_dir"); 
my(@TYPE) = (DNA::IO::TYPE5, "DIR", "DIR", "DIR", "STRING"); 
              
my(%valid_data) = &DNA::IO::evaluate($config_file, \@FIELD, \@TYPE);  # Evaluate input validity.
      
#------------------------- Valid IO parameter assignment ------------------------------
        
print "\n--> Process started - yeast_visualizer.pl";           # Process start.

my($base_dir) = $valid_data{Base_dir};			       # Variable assignment.
my($out_dir) = $valid_data{Out_dir};  

#------------------------------ Promoter object creation ------------------------------

my($prom_obj, $factors);				       # Declare promoter and factor list objects.

if ($org =~ /alb|cas|wal|par|mik|bay|gla|han|lac|lip|klu|kud/) { 

  $prom_obj = new DNA::Promoter($org, $prom);                  # Create promoter object.
  %{ $prom_obj->{parameters} } = %valid_data;	               # Assign input data to prom_obj.

  $prom_obj->{organism} = $org;
  $prom_obj->o_promoter_search();		               # Retrieve promoter sequence.

  if (defined $prom_obj->{no_eq}) {  
    print "\n--> Error: no such promoter in database - $prom\n\n";
    exit; 
  }
  else {
    $factors = new DNA::Fact_list($prom_obj);	               # Create Fact_list object.
    %{ $factors->{parameters} } = %valid_data;	               # Assign input data to factors.
        
    $factors->o_TF_search();     
  }
}
else { 
  print "\n--> Error: no such organism available - $org\n\n";
  exit;
}	

#--------------------------------- Output handling ------------------------------------

$factors->o_TF_search_output();				       # IO - output file creation.

#--------------------------------------------------------------------------------------
#                               Visualization routine 				       
#--------------------------------------------------------------------------------------

#--------------------------- TF color / rows calculation  -----------------------------

my($visu) = new DNA::Visu($prom_obj, $factors);		       # Create Visu object.
$visu->{gene} = $prom;					       # Assign promoter name to visu.
%{ $visu->{parameters} } = %valid_data;                        # Assign input data to visu.

my($f_rows) = $visu->TF_color_rows();			 
$visu->TF_align($f_rows);		        	       # Alignment routine.			 

#----------------------- Promoter sequence / TF row insertion --------------------------

$visu->window_setup();					       # Main window / canvas setup.

my($main) = $visu->{window};
my($canvas) = $visu->{canvas};

my($tb) = $visu->promoter_TF_row_insert($f_rows, $canvas);     # Return ROText field.	     			 

#--------------------------- Motif coloring / TF insertion ----------------------------

my($f_tb) = $visu->motif_color_TF_insert($tb);

#--------------------------------- Canvas update --------------------------------------

$canvas->createWindow(qw/20 20/, -window => $f_tb, -anchor => 'nw');
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

#---------------------------- End of promoter_visualizer ------------------------------

print "\n--> Location - directory: $out_dir";
print "\n--> Process finished - yeast_visualizer.pl\n\n";

MainLoop;
