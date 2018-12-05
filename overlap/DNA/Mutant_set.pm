#--------------------------------------------------------------------
# Author: Mirko Palla.
# Date: February 20, 2006.
# For: Task 3 at the Church Lab - Genetics Department, 
# Harvard Medical School.
# 
# Purpose: This program contains the complete code for 
# class Mutant_set, modeling an object for a set of mutated 
# promoter sequence storage in Perl.
#
# Associated subroutines: 0. &constructor     1. &sliding_window  
#			  2. &match_search    3. &scan	    
#			  4. &scan_d  	      5. &linear_filter   
#			  6. &perm_scan       7. &get_r_kmer_library
#-------------------------------------------------------------------- 

package DNA::Mutant_set;
                                                                                        
use strict;
use vars qw($i_marker $log_sum);
use vars qw(@l_A @l_C @l_T @l_G);

use base ("DNA::Mutant", "DNA::IO");

1;

#------------------------------------ Class Mutant_set ------------------------------------

# 0. constructor for Mutant_set, an promoter mutation set.

sub new {
  my($class, $prom) = @_;  
  my($self)={};

  $self->{prom_obj} = $prom;		 # Promoter object field.
  $self->{r_kmer};			 # Rejected kmer storage.

  $self->{permuted} = [];		 # For permutate/randomize/scan result (object) storage.
  $self->{randomized} = [];
  $self->{scanned} = [];

  $self->{parameters} = {};              # Input parameter hash.
  $self->{bar_code} = ();                # Bar code field to label output files.

  bless($self,$class);
  return($self);
}

#--------------------------------- Sliding_window sub. ---------------------------------

# 1. &sliding_window subroutine to construct permuted / mutagenized sequences.
                                                                                        
sub sliding_window($) {
  my($self, $mode) = @_;

  my($kmer_l) = $self->{parameters}{Kmer_l};
  my($d) = $self->{parameters}{Window_gap};

  my($bp_diff) = $self->{parameters}{BP_diff};

  my($f_name) = $self->{parameters}{File_name};
  my($out_dir) = $self->{parameters}{Out_dir};  

  my($s_kmer); 						  # Kmer for substitution.
  if ($mode == 1) { $s_kmer = $self->match_search($kmer_l); }  # If randomize mode, generate p-value proof kmer.

  my($prom_l) = length $self->{prom_obj}{prom_seq};    	  # Lenght of gene sequence.

  my($seq);			   	  	    	  # Promoter sequence.
  my($num_seq);	   	    			    	  # Number of window iterations.
  my($last_seq);	   	    		    	  # Excess promoter ending.	
 
  if ($kmer_l > $d) { 				    	  # Number of iterations if $k > $d. 
    $num_seq = int(($prom_l-$kmer_l)/$d) + 1;
  }
  else { $num_seq = int($prom_l/$d); } 	    	    	  # Number of iterations if $k <= $d.

  my($iter)=0;				            	  # Iteration counter.
  my($offset);				            	  # Offset of kmer substitution.

  while ($iter < $num_seq) {           	            	  # Iterate for EVEN $num_seq times.
    $seq = $self->{prom_obj}{prom_seq};
  
    $offset = $iter*$d;				    	  # Current kmer offset.
    my($kmer) = substr($seq, $offset, $kmer_l);     	  # Kmer for substitution.

    if ($mode==0) {  				    	  # If permute mode:
      $s_kmer = &DNA::Func::permute($kmer, $bp_diff);     # permuted kmer.
    }

    substr($seq, $offset, $kmer_l) = $s_kmer;       	  # Insert kmer back to promoter.

    my($mutant) = new DNA::Mutant();			  # Create Mutant object.
      
    $mutant->{o_kmer} = $kmer; 				  # Store original kmer.
    $mutant->{s_kmer} = $s_kmer; 			  # Store substitition kmer.
    $mutant->{s_offset} = $offset; 			  # Store substitution offset.
    $mutant->{kmer_l} = $kmer_l;			  # Store kmer length. 
    $mutant->{mutant} = $seq;				  # Store permuted promoter sequence.

    if ($mode==0) { push @{$self->{permuted}}, $mutant; }    # Store permuted Mutant object.     
    elsif ($mode==1) { push @{$self->{randomized}}, $mutant; } # Store randomized Mutant object.     
    else { die "\n--> Error: $mode - not sliding window mode\n--> Usage of argument: 0|1!\n"; }
    $iter++;                                              
  }
  
  my($last_seq) = $prom_l - ($offset+$d); 		   # Length of last promoter segment.

  my($end_l);						   
  if ($kmer_l > $d) { 
    $end_l = $prom_l - ($offset+$kmer_l);		   # Lenght of remainder promoter segment. 	           
  }

  if (!($last_seq == 0 || $end_l == 0)) {                  # Deal with last segment.
    my(@kmer);
    
    $seq = $self->{prom_obj}{prom_seq};			   # Update to original promoter.
    my($last_offset) = $offset + $d;		           # Offset for last insertion.

    if ($last_seq > $kmer_l) { $last_seq = $kmer_l; }      # Modify parameters.
    if ($kmer_l > $d) { $kmer_l = $last_seq; }

    my($kmer) = substr($seq, $last_offset, $kmer_l);       # Present kmer for substitution.
                                                                                        
    my($mutant_l) = new DNA::Mutant();			   # Create last Mutant object.
      
    if ($mode==0) {  				    	   # If permute mode:
      $s_kmer = &DNA::Func::permute($kmer, $bp_diff);      # permuted kmer.

      substr($seq, $last_offset, $kmer_l) = $s_kmer;	   # Insert kmer into promoter.

      $mutant_l->{o_kmer} = $kmer; 		           # Store original kmer.
      $mutant_l->{s_kmer} = $s_kmer; 			   # Store substitition kmer.
      $mutant_l->{s_offset} = $last_offset;		   # Store substitution offset.
      $mutant_l->{kmer_l} = $last_seq;			   # Store kmer length. 
      $mutant_l->{mutant} = $seq;			   # Store permuted promoter sequence.

      push @{$self->{permuted}}, $mutant_l;		   # Store permuted Mutant object.     
    }
    elsif ($mode==1) {
      my($f_kmer) = substr($s_kmer, 0, $last_seq);	   # Last random kmer segment.
      substr($seq, $last_offset, $kmer_l) = $f_kmer;	   # Insert kmer into promoter.

      $mutant_l->{o_kmer} = $kmer; 		           # Store original kmer.
      $mutant_l->{s_kmer} = $f_kmer; 			   # Store substitition kmer.
      $mutant_l->{s_offset} = $last_offset;		   # Store substitution offset.
      $mutant_l->{kmer_l} = $last_seq;			   # Store kmer length. 
      $mutant_l->{mutant} = $seq;			   # Store randomized promoter sequence.

      push @{$self->{randomized}}, $mutant_l;	           # Store randomized Mutant object.     
    }
    else { die "\n--> Error: $mode - not sliding window mode\n--> Usage of argument: 0|1!\n"; }
  }
}

#---------------------------------- Match_search sub. ----------------------------------

# 2. &match_search subroutine to find p-value proof kmer.
                                                                                        
sub match_search($$) {
  my($self, $kmer_l, $kmer_p) = @_;

  my(@kmer);					      # Kmer substitution possibility.
  if (defined $kmer_p) { @kmer = split(//, $kmer_p); }
  else { @kmer = &DNA::Func::randomize($kmer_l); } 

  my($p_marker)=0;
  my($log);

  @l_A=(); 					      # Null for recursion.
  @l_C=();		
  @l_T=();
  @l_G=(); 

  chdir $self->{parameters}{PSSMS};

  use constant THRES_HOLD => 0.01;		      # Set threshold for yeast.
  my(@p_files) = <*.pssm>;

  my($matrix) = $self->{parameters}{Matrix};
  my($file_num) = scalar @p_files;		      # Number of PSSM matrices.

  MARK: foreach my $x (@p_files) {
    open(INFO, $x) || die "--> Cannot open file!\n";

    my(@pssm);
    while (<INFO>) { push(@pssm, $_); }

    my($line_A) = $pssm[0]; my($line_C) = $pssm[1];
    my($line_T) = $pssm[2]; my($line_G) = $pssm[3]; 

    @l_A = split(/\s/, $line_A);  	              # Put pssm values into an array, @l_base.
    @l_C = split(/\s/, $line_C);
    @l_T = split(/\s/, $line_T);
    @l_G = split(/\s/, $line_G);

    $log=();                                                                                       
    my($length) = (scalar @l_A)+1;            	      # Number of pssm values in a row.

    if ($length < $kmer_l) {		              # If pssm row_l less than kmer_l.
      my($iter) = ($kmer_l-$length)+1;

      my($m)=0;
      while ($m < $iter) {                             # Sliding window iteration.
        my($g)=0;
        $log = &DNA::Func::base_iterate($m, $g, $length, 
            $matrix, @kmer);

        if ($log > THRES_HOLD) {                      # Set threshold condition.
          $p_marker=1;
          last MARK;                                  # If true, go to end of MARK.
        }
        else { $log_sum = $log_sum + $log;}	      # If false, add log-values iteratively.
        $m++;
      }
    }
    elsif ($length > $kmer_l) {		              # If pssm_l greater than kmer_l.
      my($m)=0;
      my($g) = int(($length-$kmer_l)/2);              # Set starting point for sliding.

      $log = &DNA::Func::base_iterate($m, $g, $kmer_l, 
          $matrix, @kmer);

      if ($log > THRES_HOLD) {
        $p_marker=1;
        last MARK;
      }
      else { $log_sum = $log_sum + $log; }	      # If false, add log-values iteratively.
    }
    else {                		              # Else pssm_l equals kmer_l. 
      my($m)=0; my($g)=0;
      $log = &DNA::Func::base_iterate($m, $g, $kmer_l, 
          $matrix, @kmer);
    
      if ($log > THRES_HOLD) {
        $p_marker=1;
        last MARK;
      }
      else { $log_sum = $log_sum + $log; }	 # If false, add log-values iteratively.
    }
  }
  
  if ($p_marker==1) {             		 # If threshold condition holds: repeat with new kmer.

    $log_sum=0;					 # Reset log-summation.
    $i_marker=0;				 # Reset reverse complement checker.

    my($rej) = &DNA::Func::a_s(@kmer);		 # Get rejected kmer in string form.
    %{$self->{r_kmer}}->{$rej} = $log;	         # Store rejected kmer in hash.

    $self->match_search($kmer_l);		 
  }
  else {
    if ($i_marker==1) {

      my($use) = &DNA::Func::a_s(@kmer);	
      my($ave_log) = $log_sum/(2*$file_num);	 # Average log-value for substitution kmer.
	 
      push @{$self->{kmer}}, $use, $ave_log;     # Store permuted sequence.
      push @{$self->{$kmer_l}}, $use;		 # Store kmer in hash by length.

      $log_sum=0;				 # Reset log-summation.
      $i_marker=0;				 # Reset reverse complement checker.

print "\n--> Kmer to use: $use ($ave_log) [$kmer_l]";
     
      return $use; 
    }
    else {  
      my($kmer) = &DNA::Func::a_s(@kmer);	
      my($rc_kmer) = &DNA::Func::rev_complement($kmer); 

      $i_marker=1;
      $self->match_search($kmer_l, $rc_kmer);
    }
  }
}

#------------------------------------- Scan sub. --------------------------------------

# 3. &scan subroutine to keep track of kmer occurences/locations.
                                                                                        
sub scan() {
  my($self) = shift;
  my(%seq_count);					      # Kmer storage, a hash.

  my($kmer_l) = $self->{parameters}{Kmer_l};

  my($prom) = $self->{prom_obj}->{prom_seq};                  # Retrieve promoter sequence.
  my($prom_l) = length $prom;                                 # Get sequence length.
  my($iter) = ($prom_l-$kmer_l)+1;                            # Get no. of iterations.

  my($i)=0;
  while ($i<$iter) {
    my($kmer) = substr($prom, $i, $kmer_l);

    if (!exists ($seq_count{$kmer})) { 			      # If not in hash counter.
      push @{ $seq_count{$kmer}->{$kmer} }, $i;		      # Create new member. 

      my(@bases) = ('A','C','T','G');

      my($a)=0;
      while ($a<$kmer_l) {
        my($b)=0;
        my($n_kmer) = $kmer;
        while ($b<4) {
          substr($n_kmer, $a, 1) = $bases[$b];

          if (exists ($seq_count{$n_kmer})) {		                   # If in hash already.
            if (!(defined ($seq_count{$n_kmer}->{$kmer}))) {
                           @{ $seq_count{$n_kmer}->{$kmer} } = (); }

            my($mark);
            SCAN: foreach my $key (keys %{ $seq_count{$n_kmer} }) {
   	      for (my $k=0; $k < @{ $seq_count{$n_kmer}->{$key} }; $k++) {   # Eliminiate overlap.
                if (abs($seq_count{$n_kmer}->{$key}->[$k]-$i) < $kmer_l) {
                  $mark=1;
                  last SCAN;
                }
              }
            }
            if ($mark!=1) { 
              push @{ $seq_count{$n_kmer}->{$kmer} }, $i; 
              @{ $seq_count{$kmer}->{$n_kmer} } = @{ $seq_count{$n_kmer}->{$n_kmer} };  
            } 
          }
          $b++;
        }
        $a++;
      }
    } 
    else {							        # If already in hash.
      my($marker);
      SCAN: foreach my $key (keys %{ $seq_count{$kmer} }) {
   	      for (my $k=0; $k < @{ $seq_count{$kmer}->{$key} }; $k++) {   # Eliminiate overlap.
                if (abs($seq_count{$kmer}->{$key}->[$k]-$i) < $kmer_l) {
                  $marker=1;
                  last SCAN;
                }
              }
      }
      if ($marker!=1) { push @{ $seq_count{$kmer}->{$kmer} }, $i; }    	# Record location.
    }
    $i++;
  }
  return %seq_count;
}

#------------------------------------ Scan_d sub. -------------------------------------

# 4. &scan_d subroutine to keep track of kmer occurences/locations with sliding window.
    
sub scan_d() {
  my($self) = shift;
  my(%seq_count);					 # Kmer storage, a hash.

  my($kmer_l) = $self->{parameters}{Kmer_l};
  my($d) = $self->{parameters}{Window_gap};

  my($prom_l) = length $self->{prom_obj}->{prom_seq};    # Lenght of promoter sequence.

  my($prom);			   	  	    	 # Promoter sequence.
  my($num_seq);	   	    			    	 # Number of window iterations.
  my($last_seq);	   	    		    	 # Excess promoter ending.
  
  if ($kmer_l > $d) { 				    	 # Number of iterations if $k > $d. 
    $num_seq = int(($prom_l-$kmer_l)/$d) + 1;
  }
  else { $num_seq = int($prom_l/$d); } 	    	    	 # Number of iterations if $k <= $d.

  my($iter)=0;				            	 # Iteration counter.
  my($offset);				            	 # Offset of kmer substitution.

  while ($iter<$num_seq) {
    $prom = $self->{prom_obj}{prom_seq};
  
    $offset = $iter*$d;				    	 # Current kmer offset.
    my($kmer) = substr($prom, $offset, $kmer_l);

    if (!exists ($seq_count{$kmer})) { 			 # If not in hash counter.
      push @{ $seq_count{$kmer}->{$kmer} }, $offset;	 # Create new member. 

      my(@bases) = ('A','C','T','G');

      my($a)=0;
      while ($a<$kmer_l) {
        my($b)=0;
        my($n_kmer) = $kmer;
        while ($b<4) {
          substr($n_kmer, $a, 1) = $bases[$b];

          if (exists ($seq_count{$n_kmer})) {		                   # If in hash already.
            if (!(defined ($seq_count{$n_kmer}->{$kmer}))) {
                           @{ $seq_count{$n_kmer}->{$kmer} } = (); }

            my($mark);
            SCAN: foreach my $key (keys %{ $seq_count{$n_kmer} }) {
   	      for (my $k=0; $k < @{ $seq_count{$n_kmer}->{$key} }; $k++) { # Eliminiate overlap.
                if (abs($seq_count{$n_kmer}->{$key}->[$k]-$offset) < $kmer_l) {
                  $mark=1;
                  last SCAN;
                }
              }
            }
            if ($mark!=1) { 
              push @{ $seq_count{$n_kmer}->{$kmer} }, $offset; 
              @{ $seq_count{$kmer}->{$n_kmer} } = @{ $seq_count{$n_kmer}->{$n_kmer} };  
            } 
          }
          $b++;
        }
        $a++;
      }
    } 
    else {							           # If already in hash.
      my($marker);
      SCAN: foreach my $key (keys %{ $seq_count{$kmer} }) {
   	      for (my $k=0; $k < @{ $seq_count{$kmer}->{$key} }; $k++) {   # Eliminiate overlap.
                if (abs($seq_count{$kmer}->{$key}->[$k]-$offset) < $kmer_l) {
                  $marker=1;
                  last SCAN;
                }
              }
      }
      if ($marker!=1) { push @{ $seq_count{$kmer}->{$kmer} }, $offset; }   # Record location.
    }
    $iter++;
  }
  return %seq_count;
}

#---------------------------------- Linear_filter sub. --------------------------------

# 5. &linear_filter subroutine to filter out redundant kmer mutations.

sub linear_filter($) {
  my($self, %seq_count) = @_;

  my(%finalized);
  foreach my $b (sort keys %seq_count) {
    foreach my $a (sort keys %{ $seq_count{$b} }) {
      my(@sub_set) = @{ $seq_count{$b}->{$a} };
      push @{ $finalized{$b} }, @sub_set;
    }
  }

  my(@remove_kmer);						# Redundant kmer storage.

  while ( (my $f, my $v) = each %finalized ) {		        # Point-mutate kmer -> f.

    my($kmer_l) = length $f;					# Get lenght of kmer.
    my(@bases) = ('A','C','T','G');
    
    my($a)=0;
    while ($a<$kmer_l) {
      my($base) = substr($f, $a, 1);
      
      my($b)=0;
      my($n_kmer) = $f;
      LOOP: while ($b<4) {
        if ($bases[$b] ne $base){ 
          substr($n_kmer, $a, 1) = $bases[$b];

          if (exists ($finalized{$n_kmer})) {                   # If in hash already.
      
            my(@array1) = @{ $finalized{$f} };
            my(@array2) = @{ $finalized{$n_kmer} };
      
            if (@array1==@array2) {				# If array lenghts are the same.
              my(@diff); my(%count);

              foreach my $e (@array1, @array2) { $count{$e}++ }	# Calculate symmetric difference.
              foreach my $e (keys %count) {
                if ($count{$e} != 2) { push @diff, $e; }
              }
        
              if (@diff==0) { 
                delete ($finalized{$n_kmer}); 
                delete ($seq_count{$n_kmer}); 
                last LOOP;
              }
            }							# Array lenght eq. loop END.
          }							# Point mutation existance loop END.
        }
        $b++							# Base non-equality loop END.
      }
      $a++;							# Base iteration loop END.
    }								# Kmer position iteration loop END.
  }								# Kmer value iteration loop END.
  return %seq_count;   						# Return finalized kmer hash.
}                 

#----------------------------------- Perm_scan sub. -----------------------------------

# 6. &perm_scan subroutine (mutator) to permuted kmer occurences.

sub perm_scan($$) {
  my($self, $seq_count, $mode) = @_;
  my(%seq_count) = %{ $seq_count };			       # Cast hash reference into hash.
  
  my($kmer_l) = $self->{parameters}{Kmer_l};
  my($bp_diff) = $self->{parameters}{BP_diff};
  
  my(@kmers) = keys %seq_count;				       # Get kmers from hash counter.

  foreach my $o (@kmers) {
    my($mutant) = new DNA::Mutant();			       # Create Mutant object.
    $mutant->{o_kmer} = $o; 				       # Store original kmer.
    $mutant->{kmer_l} = $kmer_l;			       # Store kmer length. 

    my($scan_seq) = $self->{prom_obj}->{prom_seq};             # Retrieve promoter sequence.

    foreach my $key (keys %{ $seq_count{$o} }) {
      for (my $k=0; $k < @{ $seq_count{$o}->{$key} }; $k++) {  # Retrieve kmer from promoter.

        my($s_kmer);
        if ($mode==0) {
          $s_kmer = &DNA::Func::permute($o, $bp_diff);         # Permuted kmer.
        }

        if ($mode==1) {
          my($random) = int(rand(10));                         # Generates a random number: [0,9].
          my(%r_kmer) = DNA::IO::R_KMER;
          $s_kmer = $r_kmer{$kmer_l}[$random];                 # Retrieves a randomized kmer from lib
        }
        
        my($offset) = $seq_count{$o}->{$key}->[$k];
        substr($scan_seq, $offset, $kmer_l) = $s_kmer;         # Insert kmer back to promoter.
        
        %{ $mutant->{s_kmer} }->{$s_kmer} = $offset;           # Store substitition (kmer,offset).
        $mutant->{no_elements}++;      		     	       # Increment no. of modifications.
      }
    }
    $mutant->{mutant} = $scan_seq;			       # Store permuted promoter sequence.
    push @{ $self->{scanned} }, $mutant;	               # Store scanned Mutant object.     
  }
  @{ $self->{scanned} } = sort { $b->{no_elements} <=> 
                                $a->{no_elements} } @{ $self->{scanned} };
}

#------------------------------ Get_r_kmer_library sub. ------------------------------

# 7. &get_r_kmer_library subroutine to create randomized kmer library.

sub get_r_kmer_library() {
  my($self) = shift;

  my($pssm_dir) = $self->{parameters}{PSSMS};

  chdir $pssm_dir;
  my(@p_files) = <*.pssm>;

  LO: foreach my $f (@p_files) {
    open(INFO, $f) || die "--> Cannot open PSSMS file: $f\n";
    
    print "\n\n*** Current file: $f ***";

    my(@pssm);
    while(<INFO>) { push(@pssm, $_); }
    close INFO;

    my(@l_A) = split(/\s/, $pssm[0]); 
    my($kmer_l) = scalar @l_A;

    if (6 <= $kmer_l && $kmer_l <= 12) {     # Use kmer length restriction on algorithm: [5,12] 
      if (!(exists $self->{$kmer_l})) {	     # Generate randomized kmer library.
        for (my $i=0; $i<10; $i++) {
          my($r_kmer) = $self->match_search($kmer_l);
        }
      } 
      else { print "\n--> Kmer length library already exists: $kmer_l!"; }
    }
    else { print "\n--> $kmer_l not in kmer length range: [6,12]!"; }
  }
}
