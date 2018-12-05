#-----------------------------------------------------------
# Author: Mirko Palla.
# Date: July 21, 2006.
# For: Task 6 at the Church Lab - Genetics Department, 
# Harvard Medical School.
# 
# Purpose: This program takes the file containing 'cer' and 
# other yeast species orthologs and creates a list of yeast 
# organism names present in the database in Perl.
#
# Input: database of orthologous yeast promoters
# 
# Output: a list of yeast organism names present in analysis
#------------------------------------------------------------ 

#!/usr/bin/perl

use strict;

#------------------------------ User input handling -----------------------------------

my($o_db) = $ARGV[0];
my($out_dir) = "/root/Church_lab/Polypromoter/output/";

if (@ARGV<1) {
  die "\n--> Error: not correct input!\n--> Usage: perl org_list.pl ortholog_db\n\n";
}

#---------------------------- File format correction ----------------------------------			              

print "\n--> Process started - org_list.pl";	           # Process start.

open(DB, $o_db) || die "--> Cannot open file list: $o_db!\n";

my(@orthologs);
while (<DB>) { chomp; push @orthologs, $_; }	           # Retrieve list of orthologs to be evaluated. 
close DB;

my(%org_hash);				            	   # Hash to store conversion hashes for e.a. organism. 
foreach my $o (@orthologs) {
  if ($o =~ /^>/) {
    my(@ortho_groups) = split(/\t/, $o);                   # Put line member strings into an array, @ortho.
    my($max_iter) = scalar @ortho_groups;

    for (my $index=2; $index<$max_iter; $index++) {
      my($org_prom) = $ortho_groups[$index];
      my($org) = substr($org_prom, 0, 3);

      if ( !(exists $org_hash{$org}) ) { $org_hash{$org} = 1; }
    }
  }
}

if (!(-d $out_dir)) {		                            # Check if directory exists. 
  mkdir ($out_dir, 0777) || die "--> Error: unable to create output directory $out_dir\n\n";
  chdir $out_dir;
}
else { chdir $out_dir; }

open OUTPUT, ">org_list.txt";                               # Open output file.
print OUTPUT "--> List of yeast organisms in database:\n\n";

my($i)=0;
foreach my $key (sort keys %org_hash) {
  print OUTPUT "$i -> $key\n";
  $i++;
}
close OUTPUT;

print "\n--> Location - directory: $out_dir";
print "\n--> Process finished - org_list.pl\n\n";
