#----------------------------------------------------------------
# Author: Mirko Palla.
# Date: May 2, 2006.
# For: Task 4 at the Church Lab - Genetics Department, 
# Harvard Medical School.
# 
# Purpose: This program takes a set of .fasta files containing 
# promoter sequence maps and converts them into standard format:
#
# >[promoter_name] (Chr. [chrom. number]: [range]) <tab> prom_seq
# in Perl.
# 
# Output: standarized individual promoter sequence files 
#----------------------------------------------------------------- 

#!/usr/bin/perl

use strict;
use DNA::Func;                                                        # Import Func package.

#------------------------------ User input handling -----------------------------------

my($in_dir) = "/root/Church_lab/Polypromoter/promoters/new_promoters/fasta/";
my($out_dir) = "/root/Church_lab/Polypromoter/promoters/new_promoters/standard/";

#---------------------------- .fasta file retrieving -----------------------------------

print "\n--> Process started - fasta_to_standard.pl";	    # Process start.

chdir $in_dir || die "--> Error: unable to change to input directory $in_dir\n\n";
my(@f_files) = <*.fasta>;

my(%standard);
foreach my $f (@f_files) {
  
  chdir $in_dir || die "--> Error: unable to change to input directory $in_dir\n\n";
  open(INFO, $f) || die "--> Cannot open file: $f!\n";
  
  my(%f_hash);
  my(@tf) = split(/\./, $f);

  my(@fasta); 
  while (<INFO>) { push(@fasta, $_); }
  chomp @fasta;

  if (!(-d $out_dir)) {		                            # Check if directory exists. 
    mkdir ($out_dir, 0777) || die "--> Error: unable to create output directory $out_dir\n\n";
    chdir $out_dir;
  }
  else { chdir $out_dir; }

  my($out_file) = $tf[0].".seq";	                    # Create file name - [.seq]
  open OUTPUT, ">$out_file";                                # Open output file.

  my($prom_seq);
  foreach my $line (@fasta) { 

    if(substr($line, 0, 1) eq '>') {
      if (defined($prom_seq)) { print OUTPUT "$prom_seq\n"; }

      my(@prom) = split(/\s/, $line);          	      # Put line member strings into an array, @prom.
      my($orf) = (split(/--/, $prom[0]))[1];          # Get ORF (promoter) name.
      my($chro) = $prom[2];			      # Get chromosome number.	 
      my(@range) = split(/-/, $prom[3]);	      # Get global promoter range.

      my($new_h) = ">$orf (Chr. $chro: [$range[0],$range[1]])\t";  
      print OUTPUT "$new_h";
      
      $prom_seq=();
    }
    else { $prom_seq = $prom_seq.$line; }
  }
  print OUTPUT "$prom_seq\n";

  close INFO;
  close OUTPUT;
}

print "\n--> Location - directory: $out_dir";
print "\n--> Process finished - fasta_to_standard.pl\n\n";
