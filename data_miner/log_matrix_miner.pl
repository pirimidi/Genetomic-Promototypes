#-----------------------------------------------------------
# Author: Mirko Palla.
# Date: March 15, 2006.
# For: Task 3 at the Church Lab - Genetics Department, 
# Harvard Medical School.
# 
# Purpose: This program takes raw TF matices from:
#
# http://fraenkel.mit.edu/improved_map/v1.tamo
#
# and constructs individual TF matrix textfiles in Perl.
#
# Input: TF matrix directory
# 
# Output: individual TF matrix files [ex.: FUS1.pssm] 
#------------------------------------------------------------ 

#!/usr/bin/perl

use strict;

#------------------------------ User input handling -----------------------------------

my($matrix_db) = $ARGV[0];
my($out_dir) = "/root/Church_lab/Polypromoter/pssms/log_matrices/";

if (@ARGV<1) {
  die "\n--> Error: not correct input!\n--> Usage: perl matrix_miner.pl matrix_db\n\n";
}

#---------------------------- TF matrix retrieving ---------------------------------

print "\n--> Process started - matrix_miner.pl";	    # Process start.

open(INFO, $matrix_db) || die "--> Cannot open file!\n";

my(@matrix_db);
while (<INFO>) { push @matrix_db, $_; }  

my(@matrices); my(@TFs);

for (my $k=0; $k<@matrix_db; $k++) {

  my(@pssm, @tf);
  my(@l_A, @l_C, @l_T, @l_G);

  my($line) = $matrix_db[$k]; 

  if ($line =~ "Log-odds") {

    @l_A = split(/\s+/, $matrix_db[$k+2]);
    @l_C = split(/\s+/, $matrix_db[$k+3]);
    @l_T = split(/\s+/, $matrix_db[$k+4]);
    @l_G = split(/\s+/, $matrix_db[$k+5]);

    @pssm = ([@l_A], [@l_C], [@l_T], [@l_G]);
    push @matrices, [ @pssm ];			            # Store pssm array linearly. 
  } 

  if ($line=~ "Source") { 
    @tf = split(/\s+/, $line); 
    push @TFs, $tf[1];				            # Store TF linearly. 
  } 
}  
close INFO;

if (!(-d $out_dir)) {		                            # Check if directory exists. 
  mkdir ($out_dir, 0777) || die "--> Error: unable to create output directory $out_dir\n\n";
  chdir $out_dir;
}
else { chdir $out_dir; }

my($out_file);
for (my $n=0; $n<@TFs; $n++) {
  $out_file = $TFs[$n].".pssm";		                    # Create file name - [TF.pssm]
  open OUTPUT, ">$out_file";                                # Open output file.

  my(@pssm) = @{ $matrices[$n] };			    # No. of columns in matrix. 

  for (my $i=0; $i<4; $i++) {
    for (my $j=1; $j<@{ $pssm[$i] }; $j++) { 
      print OUTPUT "@{ $pssm[$i] }[$j] ";
    }  
    print OUTPUT "\n";
  }  
}

print "\n--> Location - directory: $out_dir";
print "\n--> Process finished - matrix_miner.pl\n\n";
