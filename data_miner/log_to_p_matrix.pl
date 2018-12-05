#----------------------------------------------------------------
# Author: Mirko Palla.
# Date: April 24, 2006.
# For: Task 3 at the Church Lab - Genetics Department, 
# Harvard Medical School.
# 
# Purpose: This program takes a set of log-odds TF matrices
# and converts them into probability matrices in Perl.
# 
# Output: individual probability TF matrix files [ex.: FUS1.pssm] 
#----------------------------------------------------------------- 

#!/usr/bin/perl

use strict;
use DNA::Func;                                                        # Import Func package.

#------------------------------ User input handling -----------------------------------

my($in_dir) = "/root/Church_lab/Polypromoter/pssms/log_matrices/";
my($out_dir) = "/root/Church_lab/Polypromoter/pssms/prob_matrices/";

#---------------------------- TF matrix retrieving ---------------------------------

print "\n--> Process started - log_to_p_matrix.pl";	    # Process start.

chdir $in_dir || die "--> Error: unable to create output directory $in_dir\n\n";
my(@p_files) = <*.pssm>;

my(@bases) = ('A','C','T','G');

my(%matrices);
foreach my $f (@p_files) {
  open(INFO, $f) || die "--> Cannot open file: $f!\n";
  
  my(%p_hash);
  my(@tf) = split(/\./, $f);

  my(@pssm);
  while (<INFO>) { push(@pssm, $_); }

  my($line_A) = $pssm[0]; my($line_C) = $pssm[1];
  my($line_T) = $pssm[2]; my($line_G) = $pssm[3]; 

  my(@l_A) = split(/\s/, $line_A);  	              # Put pssm values into an array, @l_base.
  my(@l_C) = split(/\s/, $line_C);
  my(@l_T) = split(/\s/, $line_T);
  my(@l_G) = split(/\s/, $line_G);

  @pssm = ([@l_A], [@l_C], [@l_T], [@l_G]);

  for (my $i=0; $i<4; $i++) {
    my($base) = $bases[$i];

    my(@line);
    for (my $j=0; $j<@{ $pssm[$i] }; $j++) { 
   
      my($odds) = @{ $pssm[$i] }[$j];
      my($prob) = &DNA::Func::log_to_p($odds, $base);

      my($round) = sprintf("%.3f", $prob);            # Display 3 digits after decimal.
      push @line, $round; 
    } 
    $p_hash{$base} = [ @line ];
  }
  $matrices{$tf[0]} = \%p_hash;	                      # Store pssm array linearly.    
  close INFO;
}  

if (!(-d $out_dir)) {		                            # Check if directory exists. 
  mkdir ($out_dir, 0777) || die "--> Error: unable to create output directory $out_dir\n\n";
  chdir $out_dir;
}
else { chdir $out_dir; }

my($out_file);
foreach my $tf (keys %matrices) {
  $out_file = $tf.".pssm";		                    # Create file name - [TF.pssm]
  open OUTPUT, ">$out_file";                                # Open output file.

  my(%p_hash) = %{ $matrices{$tf} };			    # No. of columns in matrix. 

  for (my $i=0; $i<4; $i++) {
    my($base) = $bases[$i];

    my($columns) = scalar @{ $p_hash{$base} };

    for (my $j=0; $j<$columns; $j++) { 
      print OUTPUT "$p_hash{$base}[$j] ";
    }  
    print OUTPUT "\n";
  }  
}

print "\n--> Location - directory: $out_dir";
print "\n--> Process finished - log_to_p_matrix.pl\n\n";
