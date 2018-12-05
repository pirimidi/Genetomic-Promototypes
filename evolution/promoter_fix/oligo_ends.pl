#-----------------------------------------------------------
# Author: Mirko Palla.
# Date: July 19, 2006.
# For: BLAST task 2 at the Church Lab - Genetics Department,
# Harvard Medical School.
#
# Purpose: This program converts an oligo-database in text
# format into standard .fasta format. Then by running
# command:
#
#       formatdb -i name_of_db.nt -p F
#
# creates seven index files that Standalone BLAST needs to
# perform the searches and produce results in Perl.
#
# Input: text_file
#
# Output: .fasta file with the appropiate promoter data base
#------------------------------------------------------------

#!/usr/bin/perl

use strict;
use Cwd;

use Tie::IxHash;

#------------------------------ User input handling -----------------------------------

my($text_file) = $ARGV[0];

if (@ARGV<1) {
  die "\n--> Error: not correct input!\n--> Usage: perl oligo_ends.pl text_file\n\n";
}

#---------------------------- ORF sequence retrieving ---------------------------------

print "\n--> Process finished - oligo_ends.pl";

open(INFO, $text_file) || die "--> Cannot open file!\n";

my(@lines);
while (<>) { chomp; push @lines, $_; }                      # Retrieve list of promoter files to$close INFO;
close INFO;

tie my(%primer_hash), "Tie::IxHash";                        # Remember insertion order.

my($last) = (scalar @lines)-1;

my(@primers);
for (my $i=0; $i<@lines; $i++) {
  my($line) = $lines[$i];
  my($first) = substr($line, 0, 20);

  if($i==0) {                                               # Take care of first primer set (see$    $primer_hash{$first} = 1;
    push @primers, $first;
  }

  if(!(defined $primer_hash{$first})) {
    $primer_hash{$first} = 1;

    my($catch) = $i-1;
    my($middle) = substr($lines[$catch], 20, 100);

    $primer_hash{$primers[-1]} = $middle;
    push @primers, $first;
  }

  if($i==$last) { $primer_hash{$primers[-1]} =               # Deal with last primer set member.
                  substr($lines[$last], 20, 100); }
}

my($out_file) = "FUS1_out.txt";
open OUTPUT, ">$out_file";                                   # Open output file.

my($i)=1;
foreach my $key (keys %primer_hash) {
    print OUTPUT "$primer_hash{$key}\n";
    #print OUTPUT ">Primer set ($i) -> [$key]: $primer_hash{$key}\n";  # Check primers and middl$    $i++;
}
close OUTPUT;
print "\n--> File written: $out_file";

my($dir) = cwd;
print "\n--> Location - directory: $dir";
print "\n--> Process finished - oligo_ends.pl\n\n";
