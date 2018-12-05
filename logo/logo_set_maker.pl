#----------------------------------------------------------------
# Author: Mirko Palla.
# Date: April 24, 2006.
# For: Task 3 at the Church Lab - Genetics Department,
# Harvard Medical School.
#
# Purpose: This program takes a set of probability TF matrices
# and constructs their LOGO in Perl.
#
# Output: individual LOGOs for motifs in .png file format
#-----------------------------------------------------------------

#!/usr/bin/perl

use strict;

#--------------------------------- LOGO creation ---------------------------------

print "\n--> Process started - logo_set_maker.pl";          # Process start.

my(@p_files) = <*.pssm>;
if (@p_files==0) { 
  print "\n--> Error: .pssm files must be in working directory.\n\n"; 
  exit;
}

foreach my $f (@p_files) {
  my(@tf) = split(/\./, $f);

  system("perl makefig.pl < $f > $tf[0].eps");
}

my(@e_files) = <*.eps>;
if (@e_files==0) { 
  print "\n--> Error: .eps files must be in working directory.\n\n"; 
  exit;
}

my($out_dir) = "/root/Church_lab/Polypromoter/task_4/a_sharp/logo/logos/";

if (!(-d $out_dir)) {		                            # Check if directory exists. 
  mkdir ($out_dir, 0777) || die "--> Error: unable to create output directory $out_dir\n\n";
}

foreach my $e (@e_files) { 
  my(@tf) = split(/\./, $e);

  system("perl eps2png $e");
  system("mv $tf[0].png $out_dir");
}

system("rm -r *.eps");

print "\n--> Location - directory: $out_dir";
print "\n--> Process finished - logo_set_maker.pl\n\n";
