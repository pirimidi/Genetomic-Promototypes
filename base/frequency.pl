#------------------------------------------------------------
# Author: Mirko Palla.
# Date: March 29, 2006.
# For: Task 3 at the Church Lab - Genetics Department, 
# Harvard Medical School.
# 
# Purpose: This program takes a file specification containing
# genomic data and calculates the global base frequency in 
# that data set [background] in Perl.
#
# Input: file containing genomic data.
# 
# Output: global base frequency printed on the screen.
#------------------------------------------------------------- 

#!/usr/bin/perl

use strict;

use DNA::IO;							      # Import IO package.
use DNA::Func;							      # Import Func package.

#------------------------------ User input handling -----------------------------------

my($seq_file) = $ARGV[0];

if (@ARGV<2) {
  die "\n--> Error: not correct input!\n--> Usage: perl frequency.pl seq_file config_file\n\n";
}

my($config_file) = &DNA::IO::validate_cf($ARGV[1]);

#----------------------------- Config. file handling ----------------------------------

my(@FIELD) = ("Base_dir"); my(@TYPE) = ("DIR"); 

my(%valid_data) = &DNA::IO::evaluate($config_file, \@FIELD, \@TYPE);  # Evaluate input validity.

#------------------------- Valid IO parameter assignment ------------------------------

print "\n--> Process started - frequency.pl";		              # Process start.

my($dir) = $valid_data{Base_dir};			              # Variable assignment.

#----------------------------- Frequency calculation --------------------------------

my(@frequency) = &DNA::Func::frequency($seq_file, $dir, 'g');
 
print "\n\n--> GLOBAL BASE FREQUENCY - data set: $seq_file\n";
print "\n--> BASE 'A': $frequency[0]";
print "\n--> BASE 'C': $frequency[1]";
print "\n--> BASE 'T': $frequency[2]";
print "\n--> BASE 'G': $frequency[3]\n";
print "\n--> Process finished - frequency.pl\n\n";
