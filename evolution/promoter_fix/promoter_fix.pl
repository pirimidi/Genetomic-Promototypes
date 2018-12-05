#-----------------------------------------------------------
# Author: Mirko Palla.
# Date: July 17, 2006.
# For: Task 6 at the Church Lab - Genetics Department, 
# Harvard Medical School.
# 
# Purpose: This program takes either a 'cas', 'klu', or 'kud' 
# promoter file and rewrites it to standard format [convertible 
# with promoter notation in cer_orthogroups.txt] in Perl.
#
# Input: promoter_file
# 
# Output: a text file in standarized format
#------------------------------------------------------------ 

#!/usr/bin/perl

use strict;

#------------------------------ User input handling -----------------------------------

my($prom_file) = $ARGV[0];
my($out_dir) = "/root/Church_lab/Polypromoter/output/";

if (@ARGV<1) {
  die "\n--> Error: not correct input!\n--> Usage: perl promoter_fix.pl prom_file\n\n";
}

#---------------------------- File format correction ----------------------------------			              

print "\n--> Process started - promoter_fix.pl";	    # Process start.

open(INFO, $prom_file) || die "--> Cannot open file: $prom_file!\n";

my(@prom_lines);
while (<>) {  
  my(@prom) = split(/\s/, $_);

  my(@prom_name) = split(/:/, $prom[0]);
  my(@proms) = split(/_/, $prom_name[1]);

  my($tail) = substr($proms[1], 6);
  $prom[0] = '>'.$proms[0].$tail;

  push @prom_lines, [ @prom ];			            # Store updated promoter info linearly. 
}  
close INFO;

if (!(-d $out_dir)) {		                            # Check if directory exists. 
  mkdir ($out_dir, 0777) || die "--> Error: unable to create output directory $out_dir\n\n";
  chdir $out_dir;
}
else { chdir $out_dir; }

open OUTPUT, ">$prom_file";                                  # Open output file.

for (my $n=0; $n<@prom_lines; $n++) {
  my(@prom) = $prom_lines[$n];
  print OUTPUT "$prom[0][0] $prom[0][1] $prom[0][2] $prom[0][3]\t$prom[0][4]\n";
}
close OUTPUT;

print "\n--> Location - directory: $out_dir";
print "\n--> Process finished - promoter_fix.pl\n\n";
