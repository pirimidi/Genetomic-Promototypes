#------------------------------------------------------------------
# Author: Mirko Palla.
# Date: June 16, 2006.
# For: Task 5 at the Church Lab - Genetics Department, 
# Harvard Medical School.
# 
# Purpose: This program contains the complete code for 
# module Combi, containing mathematical functions used 
# in the hyper-geometric calculation in hyper_g_comp.pl 
# in Perl.
#
# Associated subroutines: 0. &constructor      1. &gammln	
#			  2. &ln_binom         3. &log_hg
#			  4. &hyper_geometric  5. &global_hg_comp
#------------------------------------------------------------------ 

package DNA::Combi;

use strict;
use base ("DNA::Promoter", "DNA::IO");

1;

#--------------------------------- Class Combi -------------------------------------

# 0. constructor for Combi storage object.

sub new {
  my($class, $sets, $isect) = @_;
  my($self)={};

  $self->{sets} = $sets;		 			# Pair of promoter lists. 
  $self->{isect} = $isect;                                      # List of intersecting promoters.

  $self->{genes} = ();		 			        # Query gene pair. 

  $self->{no_proms} = {};                                       # Number of promoters in yeast.
  $self->{num_set1} = ();                                       # Number of promoters bound by gene_1. 
  $self->{num_set2} = ();                                       # Number of promoters bound by gene_2. 
  $self->{num_isect} = ();                                      # Number of promoters in the intersection.

  $self->{result} = ();						# Calculated hyper-geometric value in gene1->gene2 direction. 

  $self->{parameters} = {};                                     # Input parameters.
  			                                    
  bless($self,$class);
  return($self);
}

#------------------------------------ Gammln sub. -------------------------------------

# 1. &gammln subroutine for logarithmic gamma function implementation.

sub gammln($) {
  my($xx) = @_;
  my($x, $y, $j, $tmp, $ser);

  my(@cof) = ( 76.18009172947146,
               -86.50532032941677,
               24.01409824083091,
               -1.231739572450155,
               0.1208650973866179e-2,
               -0.5395239384953e-5 );
  $x = $xx;
  $y = $x;

  $tmp = $x+5.5;
  $tmp -= ($x+0.5) * log($tmp);
  $ser = 1.000000000190015;

  for($j = 0; $j <= 5; $j++) { $ser += $cof[$j]/++$y; }
  return (-$tmp+log(2.5066282746310005 * $ser/$x));
}

#----------------------------------- Ln_binom sub. ------------------------------------

# 2. &ln_binom subroutine for natural logarithmic binomial function implementation.

sub ln_binom($$) {
  my($n, $k) = @_;
  
  if($k==$n || $k == 0) { return(1); }
  return (gammln($n+1.0)-gammln($k+1.0)-gammln($n-$k + 1));
}

#------------------------------------ Log_hg sub. -------------------------------------

# 3. &log_hg subroutine to calculate the probability to sample exactly $k ones in an m 
# subset of an n universe given that the universe have n_i 1's.

sub log_hg($$$$) {
  my($n, $n_i, $m, $k) = @_;

  if($n_i*$m/$n > $k) { return(0); }
 
  if($n_i < $k || $m < $k) {
    print STDERR "Bad HG calling params, $n $n_i $m $k\n";
    return(0);
  }

  my($tot) = -10000;

  do {
    my($add) = (ln_binom($m, $k) + ln_binom($n-$m, $n_i-$k) 
  				 - ln_binom($n, $n_i));
    if($add > $tot) {
      my($willlog) = exp($add) + exp($tot);

      if($willlog == 0) { return(-1000); }

      $tot = log($willlog);
    } 
    else { $k = $m; }
    $k++;
  } 
  while($k < $m && $k < $n_i);
    return($tot);
}

#------------------------------- Hyper_geometric sub. -------------------------------------

# 4. &hyper_geometric subroutine to calculate the hyper-geometric distribution for given gene pair.

sub hyper_geometric() {
  my($self) = shift;
  
  my(@sets) = @{ $self->{sets} };
  my(@isect) = @{ $self->{isect} };

  my($no_proms) = DNA::Promoter::PROM_COUNT;		  # Number of promoters in yeast.

  my($num_set1) = scalar @{ $sets[0] };		          # No. of gene1 targets. 
  my($num_set2) = scalar @{ $sets[1] };			  # No. of gene2 targets.
  my($num_isect) = scalar @isect;			  # Targets in intersection.

#--------------------------- Hyper-geometric distribution -----------------------------

  my($result) = exp(log_hg($no_proms, $num_set2,          # gene1->gene2 direction. 
	$num_set1, $num_isect));

  $self->{no_proms} = $no_proms;			  # Parameter assignment for IO.
  $self->{num_set1} = $num_set1;
  $self->{num_set2} = $num_set2;
  $self->{num_isect} = $num_isect;
  
  $self->{result} = $result;
}

#--------------------------------- Global_hg_comp sub. --------------------------------

# 5. &global_hg_comp subroutine to calculate the hyper-geometric distribution 
# for all gene pairs of yeast.

sub global_hg_comp($) {
  my($self, $TF_list) = @_;

  my(%hyper_hash);				  # Hash container for global hyper-geometric TF pair comparison. 
  my(%sets);					  # Container for gene TF sets.

  my(%valid_data) = %{ $self->{parameters} };	  # Get input parameters.
  my($p_thresh) = $self->{parameters}{P_cutoff};  # Get p-value cutoff.

  my($mark)=0;
  LP1: foreach my $tf1 (@$TF_list) {
    if ($tf1 eq 'YDR026C' || $tf1 eq 'YER051W') { next LP1; }
 
    LP2: foreach my $tf2 (@$TF_list) {
      if ($tf2 eq 'YDR026C' || $tf2 eq 'YER051W') { next LP2; }

      my($key1) = $tf1.'_'.$tf2;		  # Possible key combinations tf1 <-> tf2.
      my($key2) = $tf2.'_'.$tf1;

      if (($tf1 ne $tf2) && (!(defined $hyper_hash{$key1}) &&
			     !(defined $hyper_hash{$key2}))) {
        my($factors);					# Container for Fact_list object.                                                                              
        my(@genes) = ($tf1, $tf2);			# Assign current gene pair.

	my($iter)=0;
	while ($iter < @genes) {
  	  my($query_gene) = $genes[$iter]; 

  	  #----------------------------- Data set handling ----------------------------

          if (!(defined $sets{$query_gene})) {

  	    $factors = new DNA::Fact_list();		# Create Fact_list object.
 	    $factors->{gene_q} = $query_gene;		# Assign query gene to storage.
  	    %{ $factors->{parameters} } = %valid_data;	# Assign input data to factors.

  	    $factors->gene_search();	      		# Get list of promoters bound by given TF.

          #------------------------- Promoter list retrieving -------------------------

  	    my(@set) = keys %{ $factors->{final_list} };  # Get list of promoters bound by query gene.

     	    $sets{$query_gene} = [ @set ];
     	  }  
  	  $iter++;
 	}
        #--------------------------- Intersection calculation -------------------------

        my(@a) = @{ $sets{$tf1} };		        # Declare target set a.
        my(@b) = @{ $sets{$tf2} };		       	# Declare target set b.

        my(%union);		    		        # Union tracker.
        my(%isect);  					# Intersection tracker.

        my(@union);  				        # Union container.
        my(@isect);					# Intersection container.

        foreach my $e (@a, @b) { $union{$e}++ && $isect{$e}++ }
        @isect = sort {$a cmp $b } (keys %isect);

        my(@sets);
        if (@a<@b) { @sets = ([@a], [@b]); }
        else { @sets = ([@b], [@a]);                    # Handle ordering s.t, for each pair (
               @genes = ($tf2, $tf1); }

        #----------------------------- Argument assignment ----------------------------

        my($combi_obj) = new DNA::Combi(\@sets, \@isect);  # Create Combi object.
        $combi_obj->{genes} = \@genes;			   # Assign query gene pair.

        $combi_obj->hyper_geometric();  	           # Calculate hyper-geometric values.
        my($r) = $combi_obj->{result};

        if ($r < $p_thresh) {			           # If p-value enrty meets treshhold value.
          $hyper_hash{$key1} = $combi_obj;		   # Store new Combi object in hash.  
          $mark=1;
        }
      }
    }
  }
  if ($mark==0) {     
    print "\n--> Error: no p-value below THRESHOLD - $p_thresh\n";  # If no such hit: quit.
    exit;
  }
  return \%hyper_hash;					   # Return global hyper-geometric strorage hash.
}
