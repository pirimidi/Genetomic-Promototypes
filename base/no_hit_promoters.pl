#-----------------------------------------------------------
# Author: Mirko Palla.
# Date: February 15, 2006.
# For: Task 3 at the Church Lab - Genetics Department, 
# Harvard Medical School.
# 
# Purpose: This program takes a list of 'cer' promoters in a 
# text file and creates an output file listing them with their 
# corresponding gene name in Perl.
#
# Input: prom_file, config_file.
# 
# Output: a text file with promoter <-> gene name ocnversion.
#------------------------------------------------------------ 

#!/usr/bin/perl

use strict;

use DNA::IO;						 # Import IO library.
use DNA::Func;						 # Import Func library.

#------------------------------ User input handling -----------------------------------

my($no_hit_file) = $ARGV[0];

if (@ARGV<2) {
  die "\n--> Error: not correct input!\n--> Usage: perl no_hit_promoters.pl no_hit_file config_file\n\n";
}

my($config_file) = &DNA::IO::validate_cf($ARGV[1]);

#----------------------------- Config. file handling ----------------------------------

my(@FIELD) = ("Base_dir", "Conv_file", "Out_dir"); 
my(@TYPE) = ("DIR", "FILE"); 

my(%valid_data) = &DNA::IO::evaluate($config_file, \@FIELD, \@TYPE);  # Evaluate input validity.

#------------------------- Valid IO parameter assignment ------------------------------

print "\n--> Process started - no_hit_promoters.pl";	 # Process start.

my($base_dir) = $valid_data{Base_dir};			 # Variable assignment.
my($out_dir) = $valid_data{Out_dir};			     
  
my($conv_file) = $valid_data{Conv_file};

#--------------------------- Promoter-to-gene conversion ------------------------------

open(INFO, $no_hit_file) || die "--> Cannot open no-promoter-hit file -> $no_hit_file\n";          

my(%prom_tf_hash);					 # Hash for prom-TF storage.

while (<INFO>) {
  chomp(my($prom) = $_);              			 # Retrieve promoter name from line.

  my($gene) = &DNA::Func::promoter_to_gene($prom, $conv_file, $base_dir); # Convert gene to promoter.
  $prom_tf_hash{$prom} = $gene;
}
close INFO;

#--------------------------------- Output handling ------------------------------------

my($out_file) = "no_hit_promoters.new";

if (!(-d $out_dir)) {		                         # Check if directory exists. 
  mkdir ($out_dir, 0777) || die "--> Error: unable to create output directory $out_dir\n\n";
  chdir $out_dir;
}
else { chdir $out_dir; }

open CONV, ">$out_file";                                 # Opens file for output.

my($i)=0;
foreach my $prom (sort {$a cmp $b} keys %prom_tf_hash) {
  print CONV ">$i\t$prom_tf_hash{$prom}\t$prom\n"; 
  $i++;
}

print "\n--> Location - directory: $out_dir";
print "\n--> Process finished - no_hit_promoters.pl\n\n";
