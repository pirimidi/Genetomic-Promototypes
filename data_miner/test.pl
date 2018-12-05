#!/usr/bin/perl

use strict;
use DNA::Fact_list;
use DNA::Promoter;
use base ("DNA::Fact", "DNA::IO", "DNA::Combi");

my($number) = 0.8874;
my($string) = 'NaN';

if ($string < $number) { print "--> HELLOKA!\n"; }

#--------------------------------------------------------------------------------------

my(@e_factors);

my($e_factor1) = new DNA::Fact("STE12");	      # Create new TF object.
$e_factor1->{chromosome} = 'I';		      	      # Assign chromosome number.
$e_factor1->{motif_start} = 0;
$e_factor1->{motif_end} = 10;
$e_factor1->{orient} = '+';
$e_factor1->{promoter} = "YHR084W|YHR083C";	      # Get promoter name list.                   

push @e_factors, $e_factor1;		              # Store Ernest's TF data.

my($e_factor2) = new DNA::Fact("FUS1");		      # Create new TF object.
$e_factor2->{chromosome} = 'II';	      	      # Assign chromosome number.
$e_factor2->{motif_start} = 100;
$e_factor2->{motif_end} = 200;
$e_factor2->{orient} = '-';
$e_factor2->{promoter} = "XHR055C|XHR054W";	      # Get promoter name list.                   
                          
push @e_factors, $e_factor2;		              # Store Ernest's TF data.

#--------------------------------------------------------------------------------------

my(@array);
foreach my $tf (@e_factors) {			      # Offset calculation.
  my($proms) = $tf->{promoter};		              # Get promoter name(s) string.
  
  my(@p) = split(/\|/, $proms);			      # Get list of promoter name(s).
  print "\n\n--> PROMOTER_LIST: @p\n";
            
  foreach my $p (@p) {
    $tf->{promoter} = $p;
    push @array, $tf;				      #	Create new member and assign it to storage.
  }
}

#--------------------------------------------------------------------------------------

my($i)=0;
foreach my $e (@array) {
  print "\n--> NO. $i";
  print "\n--> Factor: $e->{factor}";
  print "\n--> Chrom: $e->{chromosome}";
  print "\n--> Start: $e->{motif_start}";
  print "\n--> End: $e->{motif_end}";
  print "\n--> Orient: $e->{orient}";
  print "\n--> PROMOTER: $e->{promoter}\n\n";
  $i++;
} 
