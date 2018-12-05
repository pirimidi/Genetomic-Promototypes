#-----------------------------------------------------------
# Author: Mirko Palla.
# Date: July 17, 2006.
# For: Task 6 at the Church Lab - Genetics Department, 
# Harvard Medical School.
# 
# Purpose: This program takes the 'gla', 'han', 'lac' and 
# 'lip' promoter files and rewrites them to standard format 
# (converting gene aliases to other aliases) using the NCBI 
# [http://www.ncbi.nlm.nih.gov] web-site in Perl.
#
# Input: list of promoter_file
# 
# Output: a set of text files in standarized format
#------------------------------------------------------------ 

#!/usr/bin/perl

use strict;
use LWP::UserAgent;

#------------------------------ User input handling -----------------------------------

my($file_list) = $ARGV[0];
my($out_dir) = "/root/Church_lab/Polypromoter/output/";

if (@ARGV<1) {
  die "\n--> Error: not correct input!\n--> Usage: perl non_convertible.pl file_list\n\n";
}

#---------------------------- File format correction ----------------------------------			              

print "\n--> Process started - non_convertible.pl";	    # Process start.

open(LIST, $file_list) || die "--> Cannot open file list: $file_list!\n";

my(@files);
while (<>) { chomp; push @files, $_; }			    # Retrieve list of promoter files to be evaluated. 
close LIST;

# Create a user agent object
my($ua) = LWP::UserAgent->new;
$ua->agent("MyApp/0.1 ");
        
# Create a request
my($req) = HTTP::Request->new(POST => 'http://www.ncbi.nlm.nih.gov/entrez/query.fcgi');
$req->content_type('application/x-www-form-urlencoded');

my(%org_hash);						    # Hash to store conversion hashes for e.a. organism. 
foreach my $file (@files) {
  open(INFO, $file) || die "--> Cannot open file: $file!\n";

  print "\n--> In file: $file";

  my(@org) = split(/\_/, $file);
  my($org) = $org[0];				            # Name of organism under analysis.

  while (<INFO>) {  
    my(@prom) = split(/\s/, $_);
    my(@prom_name) = split(/\./, $prom[0]);

    my($q_gene) = substr($prom_name[0], 1);			    # Get gene (promoter) name of interest.
    my($end) = substr($q_gene, -1);

    if ($end ne 'g') {
      $req->content('term='.$q_gene.'&cmd=search&db=gene');
               
      # Pass request to the user agent and get a response back
      my($res) = $ua->request($req);
                   
      # Check the outcome of the response
      my(@content);
      if ($res->is_success) { @content = split(/\s/, $res->content); }
      else { print "Process did not work!\n"; }

      for(my $i=0; $i<@content; $i++) {
        if ($content[$i] =~ "Aliases:") {
          my(@alias) = split(/\W/, $content[++$i]);
          $org_hash{$org}{$q_gene} = $alias[3];			    # Store alias conversion in a hash.	
          #print "\n--> $q_gene\t$alias[3]"; 
        }
      }
    }
  }
  close INFO;
}  

if (!(-d $out_dir)) {		                            # Check if directory exists. 
  mkdir ($out_dir, 0777) || die "--> Error: unable to create output directory $out_dir\n\n";
  chdir $out_dir;
}
else { chdir $out_dir; }

foreach my $file (@files) {
  my(@org) = split(/\_/, $file);
  my($org) = $org[0];				            # Name of organism under analysis.

  my($out_file) = $org."_aliases.txt";
  open OUTPUT, ">$out_file";                                # Open output file.

  foreach my $key (sort keys %{ $org_hash{$org} }) {
    print OUTPUT "$key\t$org_hash{$org}{$key}\n";
  }
}
close OUTPUT;

print "\n--> Location - directory: $out_dir";
print "\n--> Process finished - non_convertible.pl\n\n";
