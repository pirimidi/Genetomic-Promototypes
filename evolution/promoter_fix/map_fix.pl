#-----------------------------------------------------------
# Author: Mirko Palla.
# Date: July 20, 2006.
# For: Task 6 at the Church Lab - Genetics Department, 
# Harvard Medical School.
# 
# Purpose: This program takes either a 'cas', 'klu', or 'kud'  
# map file and rewrites it to standard format [convertible 
# with promoter notation in cer_orthogroups.txt] in Perl.
#
# Input: map_file
# 
# Output: a text file in standarized format
#------------------------------------------------------------ 

#!/usr/bin/perl

use strict;

#------------------------------ User input handling -----------------------------------

my($map_file) = $ARGV[0];
my($out_dir) = "/root/Church_lab/Polypromoter/output/";

if (@ARGV<1) {
  die "\n--> Error: not correct input!\n--> Usage: perl map_fix.pl map_file\n\n";
}

#---------------------------- File format correction ----------------------------------			              

print "\n--> Process started - map_fix.pl";	    # Process start.

open(INFO, $map_file) || die "--> Cannot open file: $map_file!\n";

my(@data);
while (<INFO>) { @data = <INFO>; } 
close INFO;

print @data;

if (!(-d $out_dir)) {		                            # Check if directory exists. 
  mkdir ($out_dir, 0777) || die "--> Error: unable to create output directory $out_dir\n\n";
  chdir $out_dir;
}
else { chdir $out_dir; }

open OUTPUT, ">$map_file";                                  # Open output file.
for (my $n=0; $n<@data; $n++) {
  $data[$n] =~ s/ORFN:|_Contig//g; 
  #$data[$n] =~ s{ORFN:|_Contig}[]gsx; 
}  
close OUTPUT;

print "\n--> Location - directory: $out_dir";
print "\n--> Process finished - map_fix.pl\n\n";
