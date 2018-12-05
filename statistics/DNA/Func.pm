#-------------------------------------------------------------------------------
# Author: Mirko Palla.
# Date: March 2, 2006.
# For: Task 1 at the Church Lab - Genetics Department, 
# Harvard Medical School.
# 
# Purpose: This program contains the complete code for 
# module Func, containing re-occurring functionality 
# subroutines in Perl.
#
# Associated subroutines: 
#			   1. &random                2. &permute            
#			   3. &randomize 	     4. &base_iterate         
#     KMER MUTATION	   5. &rev_complement        6. &a_s  
#			   7. &frequency 	     8. &validate_2bp_diff     	 	         
#   
#     PSSMS SPECIFIC       9. &log_l_ratio	    10. &log_to_p
#			  11. &pssm_mutate          12. &overlap_mutate
#	                  13. &get_pssm		    14. &overlap_motif
#
#     DATA BASE           15. &gene_to_promoter     16. &promoter_to_gene	 
#			  17. &motif_db_nir   	    	
#
#     FILTER		  18. &overlap_dup_eliminate
#
#     MUTATOR		  29. &get_global_range     20. &flip_range
#
#     GENERAL	          21. &arabic_to_roman      22. &roman_to_arabic
#			  23. &get_TF_num           24. &get_prom_name
#			  25. &get_bind_affinity    26. &get_stringency
#			  27. &get_description      28. &get_tf_cond_list
#			  29. &get_occurence        30. &prom_to_gene_library
# 		  	  31. &TF_list		    32. &create_factor 
#			  33. &get_org_file	    34. &get_conditions  
#
#     SELECTION		  35. &binding_filter	    36. &condition_filter
#			  37. &conservation_filter  38. &probability_filter
#------------------------------------------------------------------------------- 

package DNA::Func;
                                                                                        
use strict;
use POSIX qw(ceil);
use vars qw($counter $rev);

1;

#--------------------------------------------------------------------------------------#
#	  			     KMER MUTATION				       #
#--------------------------------------------------------------------------------------#

#------------------------------------ Random sub. -------------------------------------

# 1. &random subroutine for array element permutation.

sub random($) {
  my($array) = shift;
  my($i);
  for ($i = @$array; --$i;) {
    my($j) = int rand ($i+1);
    next if $i == $j;
    @$array[$i,$j] = @$array[$j,$i];
  }
}

#------------------------------------- Permute sub. -----------------------------------
                                                                                        
# 2. &permute subroutine to randomly permute a given kmer from promoter.
                                                                                        
sub permute($$) {
  my($o_kmer, $bp_thres) = @_;

  my(@kmer) = split(//, $o_kmer);		    # Convert original kmer into array form.

  my(%bases);
  while ($o_kmer =~ /(.)/g) { $bases{$1}++; }
  my($no) = scalar (keys %bases);
      
  my($p_kmer);			                    # Initialize kmer permutation.

  if ($no != 1) {
    my($base_diff);				    # Initialize base difference.
    while ($base_diff < $bp_thres) {		    # If difference < than 2.
      random(\@kmer);		                    # Randomly permute present k-mer.
      $p_kmer = a_s(@kmer);		            # Convert permuted kmer into string form.

      $base_diff = &validate_2bp_diff($o_kmer, $p_kmer); # Calculate base difference.
    }
  }
  else { $p_kmer = $o_kmer; }			    # Kmer comtains only one type of base.

  return $p_kmer;				    # Return permuted kmer string.
}

#---------------------------------- Randomize sub. ------------------------------------

# 3. &randomize subroutine for array element randomization.

sub randomize($) {
  my($k) = shift; 
  my(@r_kmer);
  my(@base) = ('G', 'A', 'C', 'T');
  for (my($z)=0; $z<$k; $z++) 
  { push(@r_kmer, $base[rand @base]); }
  return @r_kmer;
}

#--------------------------------- Base_iterate sub. ----------------------------------

# 4. &base_iterate subroutine to determine log-value.

sub base_iterate($$$$$) {
  my($m, $g, $k, $matrix, @kmer) = @_;

  my($log_val)=0;
  my($c)=0;
  my($p); 
  my($base); 

  while($c<$k) {      
    chomp($base = $kmer[$c+$m]); 

    if ($base eq 'A') { $p = $DNA::Mutant_set::l_A[$c+$g]; }	    # Get scoring matrix value.
    elsif ($base eq 'C') { $p = $DNA::Mutant_set::l_C[$c+$g]; }
    elsif ($base eq 'T') { $p = $DNA::Mutant_set::l_T[$c+$g]; }
    elsif ($base eq 'G') { $p = $DNA::Mutant_set::l_G[$c+$g]; }

    if ($matrix==0) { 
      if ($p==0) { $p = 1E-10; }
      $log_val += log_l_ratio($p, $base);             	   	    # Calculate log-likelihood value.
    }
    elsif ($matrix==1) { $log_val += $p; }
    $c++;
  }
  return $log_val;
}

#-------------------------------- Rev_complement sub. ---------------------------------

# 5. &rev_complement subroutine to reverse-complement any sequence.
                                                                                        
sub rev_complement($) {
  my($motif) = shift;
  my($length) = length $motif;
  $rev = reverse($motif);

  my($m)=0;
  while($m<$length) {
    my($b) = substr($rev, $m, 1);
    
    if ($b eq 'A') { substr($rev, $m, 1) = 'T'; }
    elsif ($b eq 'T') { substr($rev, $m, 1) = 'A'; }
    elsif ($b eq 'G') { substr($rev, $m, 1) = 'C'; }
    elsif ($b eq 'C') { substr($rev, $m, 1) = 'G'; }

    $m++;
  }
  return $rev;
}

#--------------------------------------- A_s sub. -------------------------------------

# 6. &a_s subroutine to convert an array into a string.
                                                                                        
sub a_s($) {
  my(@a_kmer) = @_;

  my($a)=0;
  my($s_kmer);
  while ($a < @a_kmer) {
    substr($s_kmer, $a, 1) = $a_kmer[$a];
    $a++;
  }
  return $s_kmer;
}

#---------------------------------- Frequency sub. ------------------------------------

# 7.a &frequency subroutine to determine base occurrence in yeast promoters/genome.

sub frequency($$$) {
  my($seq_file, $base_dir, $mode) = @_;
  chdir $base_dir;

  my($num_base)=0;

  my($a_count)=0; my($c_count)=0;
  my($t_count)=0; my($g_count)=0; 

  open(INFO, $seq_file) || die "--> Cannot open sequence file - $seq_file!\n"; 

  my($prom_seq);
  while (<INFO>) {		 
    
    if ($mode eq 'p') {
      my(@line) = split(/\s/, $_);	
      $prom_seq = $line[4];
    }
    elsif ($mode eq 'g') { $prom_seq = $_; }
    else { print "Error: $mode - not valid mode value -> usage: {p,g}\n\n"; } 	
    
    my($length) = length $prom_seq;  	    		     # Lenght of promoter sequence.

    my($i)=0;
    while ($i < $length) {
      my($base) = substr($prom_seq, $i, 1);		     # Iterate over all bases.

      if ($base eq 'A') { $a_count++; }
      elsif ($base eq 'C') { $c_count++; }
      elsif ($base eq 'T') { $t_count++; }
      elsif ($base eq 'G') { $g_count++; }

      $i++;
    }
  }
  close INFO;

  $num_base = ($a_count+$c_count+$t_count+$g_count);

  return my(@freq) = ($a_count/$num_base, 		     # Return base freq. array.
		      $c_count/$num_base,  
	              $t_count/$num_base, 
		      $g_count/$num_base);
}

#------------------------------- Frequency sub. output --------------------------------

# 7.b Constant FREQ_P, an array, to hold yeast promoter base frequency data.

  use constant FREQ_P => [0.321659652752007, 0.1820602749244,
			  0.319752394337284, 0.176527677986309];

#------------------------------- Frequency sub. output --------------------------------

# 7.c Constant FREQ_G, an array, to hold yeast genome base frequency data.

  use constant FREQ_G => [0.309013721779486, 0.191672288490661,
			  0.308012853465539, 0.191301136264315];

#------------------------------- Validate_2bp_diff sub. -------------------------------

# 8. &validate_2bp_diff to provide kmer permutation with at least 2 bp difference.

sub validate_2bp_diff($$) {
  my($p_kmer, $o_kmer) = @_;

  my($iter) = length $p_kmer;			# Number of bases in kmer.

  my($pos)=0;					# Position of base extraction.
  my($mis_match)=0;				# Number of mis-matches.

  while ($pos < $iter) {
    my($p_base) = substr($p_kmer, $pos, 1);
    my($o_base) = substr($o_kmer, $pos, 1);

    if ($p_base ne $o_base) { $mis_match++; }
    $pos++; 
  }
  return $mis_match;
}

#--------------------------------------------------------------------------------------#
#				   PSSMS SPECIFIC				       #
#--------------------------------------------------------------------------------------#

#--------------------------------- Log_l_ratio sub. -----------------------------------

# 9. &log_l_ratio subroutine to calculate log-likelihood value from p-value.

sub log_l_ratio($$) {
  my($A, $base) = @_;				      # A: probability, B: global base frequency

  my($B)=1;
  if ($base eq 'A') { $B = FREQ_G->[0]; }	      # Assign appropriate freq. value.
  elsif ($base eq 'C') { $B = FREQ_G->[1]; }
  elsif ($base eq 'T') { $B = FREQ_G->[2]; }
  elsif ($base eq 'G') { $B = FREQ_G->[3]; } 

  return log($A/$B)/log(2);			      # Return calculated log-likelihood value.
}

#----------------------------------- Log_to_p sub. ------------------------------------

# 10. &log_to_p subroutine to calculate p-value from log-likelihood entry.

sub log_to_p($$) {
  my($A, $base) = @_;				      # A: log-odds entry, B: global base frequency
  
  my($B);
  if ($base eq 'A') { $B = FREQ_G->[0]; }	      # Assign appropriate freq. value.
  elsif ($base eq 'C') { $B = FREQ_G->[1]; }
  elsif ($base eq 'T') { $B = FREQ_G->[2]; }
  elsif ($base eq 'G') { $B = FREQ_G->[3]; } 
	
  my($p) = (exp($A*log(2)))*$B;			      # Calculate p-value.
	  
  if ($p>1) { return 1; }  			      # Return corrected p-value.
  else { return $p; }
}

#----------------------------------- Pssm_mutate sub. ---------------------------------

# 11. &pssm_mutate subroutine to 'knock out' overlap (side) by selecting minimal 
# matrix entries for motif B. 
                                                                                        
sub pssm_mutate($$) {
  my($hash, $motif_B) = @_;

  my(%p_hash) = %{ $hash };		               # Get PSSM matrix.

  my($m_B_motif) = $motif_B->{motif};                  # Motif to 'knock out'.
  my($m_length) = length $m_B_motif; 	               # Lenght of motif.

  my($j)=0;
  while ($j < $m_length) {
               
    my(@rank) = sort { $p_hash{$a}[$j] <=> $p_hash{$b}[$j] } keys %p_hash;
                  
    my($sub_base) = $rank[0];
    my($sub_value) = $p_hash{$rank[0]}[$j];

    my(@same) = ($sub_base);
    my(@bases) = ('A','C','T','G');

    my($b)=0;
    while ($b<4) {
      if($bases[$b] ne $sub_base) {
        if ($p_hash{$bases[$b]}[$j] == $sub_value) { push @same, $bases[$b]; }
      }
      $b++;
    }
                                            
    if (@same>1) { $sub_base = $same[ rand @same ]; }
                                                                               
    substr($m_B_motif, $j, 1) = $sub_base;

    $j++;
  }
  
  my($orient) = $motif_B->{orient}; 	        # Orientation of motif.
  
  if ($orient eq '-') { 
    $m_B_motif = rev_complement($m_B_motif);
  } 
  return $m_B_motif;				# Return 'knocked out' motif.
}

#--------------------------------- Overlap_mutate sub. --------------------------------

# 12. &overlap_mutate subroutine to 'dominate' overlap by selecting maximum matrix 
# entries for motif A, s.t., max(motif_A - motif_B), and motif_A > motif_A*. 
                                                                                        
sub overlap_mutate($$$$$$) {
  my($o_lap, $k, $F_motif, $matrix, $keep_perc, $PSSMS) = @_;

  my($m1) = $o_lap->{motif_A};			       # Motif object in 1st position. 
  my($m2) = $o_lap->{motif_B};			       # Motif object in 2nd position.
  
  my($overlap) = $o_lap->{overlap};		       # Overlap value.

  my(%p_hash) = %{ $PSSMS };		               # Get PSSM matrix.

  #----------------------------------------------------#
  #         Motif parameters in 1st position           #
  #----------------------------------------------------#
                                                       # 
  my($m1_name) = $m1->{factor};                        # Motif name in 1st position.
  my($m1_orient) = $m1->{orient};                      # Motif orientation in ~.
  my($m1_motif) = $m1->{motif};                        # Motif sequence in ~.
  my(%m1_pssm) = %{ $p_hash{$m1_name} };               # Motif PSSM in ~.
  my($m1_length) = length $m1_motif;                   # Lenght of motif in ~.

  #----------------------------------------------------#
  #         Motif parameters in 2nd position           #
  #----------------------------------------------------#
                                                       #
  my($m2_name) = $m2->{factor};                        # Motif name in 2nd position.
  my($m2_orient) = $m2->{orient};                      # Motif orientation in ~.
  my($m2_motif) = $m2->{motif};                        # Motif sequence in ~.
  my(%m2_pssm) = %{ $p_hash{$m2_name} };               # Motif PSSM in ~.
  my($m2_length) = length $m2_motif; 	               # Lenght of motif in ~.                   				       #
  #----------------------------------------------------#

  my($m1_pos);					       # Column entry position for PSSM 1.
  my($m2_pos);					       # Column entry position for PSSM 2.

  my($iter_1);				               # Iteration number for PSSM 1 - {+1,-1}.
  my($iter_2);				               # Iteration number for PSSM 2 - {+1,-1}.

  #------------------------- PSSM column start for iteration --------------------------

  my($tail)=0;						
  if (defined ($o_lap->{tail_l})) { $tail = $o_lap->{tail_l}; }
   
  if ($m1_orient eq '+') {			       # PSSM column start for motif 1. 
    $m1_pos = $m1_length - ($overlap + $tail);
    $iter_1 = 1;
  }
  if ($m1_orient eq '-') {
    $m1_pos = ($overlap + $tail) - 1; 
    $iter_1 = -1;
  }

  if ($m2_orient eq '+') {			       # PSSM column start for motif 2. 
    $m2_pos = 0; 
    $iter_2 = 1;
  }
  if ($m2_orient eq '-') {
    $m2_pos = $m2_length - 1; 
    $iter_2 = -1;
  }

  #--------------------------------- Overlap type -------------------------------------

  my(@n_bases) = ('A','C','T','G');			# Normal base iteration.
  my(@r_bases) = ('T','G','A','C');			# Reverse complement base iteration.

  my(@bases_1);
  my(@bases_2);
  my($OL_type);
  
  SWITCH: {
    if ($m1_orient eq '+' && $m2_orient eq '+') { $OL_type = "++"; last SWITCH; }
    if ($m1_orient eq '+' && $m2_orient eq '-') { $OL_type = "+-"; last SWITCH; }
    if ($m1_orient eq '-' && $m2_orient eq '+') { $OL_type = "-+"; last SWITCH; }
    if ($m1_orient eq '-' && $m2_orient eq '-') { $OL_type = "--"; last SWITCH; }
  }

  if ($OL_type eq "++" || $OL_type eq "--") {
    @bases_1 = @n_bases; 
    @bases_2 = @n_bases;
  } 
  elsif ($OL_type eq "+-" || $OL_type eq "-+") {
    @bases_1 = @n_bases; 
    @bases_2 = @r_bases;
  } 
  else { 
    print "\n--> No such overlap type: $OL_type"; 
    exit; 
  }

  my($i)=0;					
  while ($i < $overlap) {

    my($k_base) = substr($F_motif, $i, 1);  		  # Original overlap base (keep).

    my($k_odds);
    if ($k==0) {   					  # Log-odds entry for original OL -> mutate: M1, keep: M2.
      if ($OL_type eq "+-" || $OL_type eq "--") { 
        $k_base = rev_complement($k_base);
      }
      $k_odds = $m2_pssm{$k_base}[$m2_pos]; 
    }

    if ($k==1) {     					  # Log-odds entry for original OL -> mutate: M2, keep: M1.
      if ($OL_type eq "-+" || $OL_type eq "--") { 
        $k_base = rev_complement($k_base);
      }
      $k_odds = $m1_pssm{$k_base}[$m1_pos]; 
    }  

    my($k_prob) = log_to_p($k_odds, $k_base);
    my($m_base) = substr($F_motif, $i, 1);  		  # Original overlap base (mutate).

    my($m_odds);
    if ($k==0) {   					  # Log-odds entry for original OL -> mutate: M1, keep: M2.
      if ($OL_type eq "-+" || $OL_type eq "--") { 
        $m_base = rev_complement($m_base);
      }
      $m_odds = $m1_pssm{$m_base}[$m1_pos]; 
    }

    if ($k==1) {     					  # Log-odds entry for original OL -> mutate: M2, keep: M1.
      if ($OL_type eq "+-" || $OL_type eq "--") { 
        $m_base = rev_complement($m_base);
      }
      $m_odds = $m2_pssm{$m_base}[$m2_pos]; 
    }  

    my($m_prob) = log_to_p($m_odds, $m_base);  
            
    my(%k_hash);					# Entry (keep) storage hash.
    my(%m_hash);					# Entry (mutate) storage hash.
    my(%d_hash);					# Difference storage hash.
      
    my($c)=0;
    while ($c<4) {
      my($odds_1) = $m1_pssm{$bases_1[$c]}[$m1_pos];
      my($odds_2) = $m2_pssm{$bases_2[$c]}[$m2_pos];

      my($prob_1) = log_to_p($odds_1, $bases_1[$c]);
      my($prob_2) = log_to_p($odds_2, $bases_2[$c]);

      my($p_diff);
      if ($k==0) { 
        $p_diff = $prob_2 - $prob_1;

        $k_hash{$bases_2[$c]} = $prob_2; 
        $m_hash{$bases_2[$c]} = $prob_1; 
        $d_hash{$bases_2[$c]} = $p_diff;
      } 
      
      if ($k==1) { 
        $p_diff = $prob_1 - $prob_2;

        $k_hash{$bases_1[$c]} = $prob_1; 
        $m_hash{$bases_1[$c]} = $prob_2; 
        $d_hash{$bases_1[$c]} = $p_diff;
      } 
      $c++;
    }

    my(@rank) = sort { $d_hash{$b} <=> $d_hash{$a} } keys %d_hash;
    
    END: foreach my $max_base (@rank) {
      if ($max_base ne $k_base) {        
 
        my($max_diff) = $d_hash{$max_base};
        my($max_prob) = $k_hash{$max_base};
        my($mut_prob) = $m_hash{$max_base};
  
        my($gap) = $k_prob*($keep_perc/100);
        my($inf) = $k_prob-$gap; 
        my($sup) = $k_prob+$gap;

        if (($inf < $max_prob && $max_prob < $sup) && $mut_prob < $m_prob) { 

          if ($k==0) {
            if ($OL_type eq "+-" || $OL_type eq "--") { 
              $max_base = rev_complement($max_base); 
            }
          }  
     
          if ($k==1) {
            if ($OL_type eq "-+" || $OL_type eq "--") { 
              $max_base = rev_complement($max_base); 
            }
          }

          substr($F_motif, $i, 1) = $max_base; 
        }
      } 
    } 
    $i++;
    
    $m1_pos = $m1_pos + $iter_1; 
    $m2_pos = $m2_pos + $iter_2; 
  }
  return $F_motif;						# Return mutated overlap segment.
}

#----------------------------------- Get_pssm sub. ------------------------------------

# 13. &get_pssm subroutine to retrieve PSSM matrix of given TFBS. 
                                                                                        
sub get_pssm($$) {
  my($factor, $dir) = @_;

  chdir $dir;
  if (-f $factor.".pssm") {
    open(INFO, $factor.".pssm") || die "--> Cannot open PSSMS file: $factor.pssm\n";
  }
  else { print "\n--> No PSSMS matrix for factor: $factor"; }
           
  my(@pssm);
  while(<INFO>) { push(@pssm, $_); }
  close INFO;
  
  return @pssm;
}

#---------------------------------- Overlap_motif sub. --------------------------------

# 14. &overlap_motif subroutine to return substitution motif with 
# specified ruin percentage. 

sub overlap_motif($$$$$) {
  my($o_lap, $s_motif, $k, $type, $OL_set) = @_;

  my($m_X); 
  if ($k==0) { $m_X = $o_lap->{motif_A}; } 
  if ($k==1) { $m_X = $o_lap->{motif_B}; } 

  my($m_X_name) = $m_X->{factor};
  my($m_X_orient) = $m_X->{orient};
  my($m_X_motif) = $m_X->{motif};
  my($m_X_score) = $m_X->{score};

  my($dir) = $OL_set->{parameters}{PSSMS};
  my($mat) = $OL_set->{parameters}{Matrix};
  my($ruin) = $OL_set->{parameters}{Ruin_perc};

  my($R_THRES) = $m_X_score*($ruin/100);

  if($m_X_orient eq '-') { $s_motif = rev_complement($s_motif); }
  
  my($m_length) = length $m_X_motif;

  my($new_tf) = new DNA::Fact($m_X_name);

  if ($type eq "fake" || $type eq "real->I") {  	# Fake or real->I overlap scenario.

    my($pos, $inc);				        # Position in string, increment.

    if (($k==0 && $m_X_orient eq '+') || ($k==1 && $m_X_orient eq '-')) {					# Mutate 1, keep 2.
      $pos = 0;
      $inc = 1;
    }
    
    if (($k==0 && $m_X_orient eq '-') || ($k==1 && $m_X_orient eq '+')) {
      $pos = $m_length-1;
      $inc = -1;
    }

    my($i)=0;
    while ($i < $m_length) {

      substr($m_X_motif, $pos, 1) = substr($s_motif, $pos, 1);

      $new_tf->{motif} = $m_X_motif; 
      $new_tf->get_p_value($mat, $dir);
      my($new_score) = $new_tf->{score};          

      if ($new_score <= $R_THRES) { last; }
      
      $pos = $pos + $inc;				# Increment motif position.
      $i++;						# Increment global iteration.
    }
  }
   
  if ($type eq "real->II") {

    my($overlap) = $o_lap->{overlap};
    my($head_l) = $o_lap->{head_l};
    my($tail_l) = $o_lap->{tail_l};

    if ($k==0) {					# Mutate 1, keep 2.
 
      my($mark)=0;
      
      my($iter1);					# Parameters or 1st iteration.
      my($pos, $inc);
   
      my($iter2);					# Parameters for 2nd iteration.
      my($pos1, $inc1);
      my($pos2, $inc2);

      my(@HT);
      my(@ORIENT);

      if ($head_l >= $tail_l) { 
        @HT = ($head_l, $tail_l); 
        @ORIENT = ('+', '-'); 
      }
      
      if ($head_l <  $tail_l) { 
        @HT = ($tail_l, $head_l); 
        @ORIENT = ('-', '+'); 
      }
          
      $iter1 = $HT[0] - $HT[1];

      if ($m_X_orient eq $ORIENT[0]) { ($pos, $inc) = (0, 1); }
      if ($m_X_orient eq $ORIENT[1]) { ($pos, $inc) = ($m_length-1, -1); }

      $iter2 = ceil((($overlap + $tail_l) + ($head_l - $iter1))/2);
         
      if ($m_X_orient eq $ORIENT[0]) {
        ($pos1, $inc1) = ($iter1, 1);
        ($pos2, $inc2) = ($m_length-1, -1);
      }
          
      if ($m_X_orient eq $ORIENT[1]) {
        ($pos1, $inc1) = ($pos-$iter1, -1);
        ($pos2, $inc2) = (0, 1);
      }    

      my($i)=0;
      while ($i < $iter1) {

        substr($m_X_motif, $pos, 1) = substr($s_motif, $pos, 1);

        $new_tf->{motif} = $m_X_motif; 
        $new_tf->get_p_value($mat, $dir);
        my($new_score) = $new_tf->{score};          

        if ($new_score <= $R_THRES) { $mark=1; last; }
      
        $pos = $pos + $inc;					# Increment motif position.
        $i++;							# Increment global iteration.
      }
  
      if ($mark!=1) {
          
        my($j)=0;
        while ($j < $iter2) {
    
          my($random) = int(rand(2));
            
          my(@pos);
          if ($random==0) { @pos = ($pos1, $pos2); } 
          if ($random==1) { @pos = ($pos2, $pos1); } 

          substr($m_X_motif, $pos[0], 1) = substr($s_motif, $pos[0], 1);

          $new_tf->{motif} = $m_X_motif; 
          $new_tf->get_p_value($mat, $dir);
          my($new_score) = $new_tf->{score};          

          if ($new_score <= $R_THRES) { last; }

          substr($m_X_motif, $pos[1], 1) = substr($s_motif, $pos[1], 1);

          $new_tf->{motif} = $m_X_motif; 
          $new_tf->get_p_value($mat, $dir);
          my($new_score) = $new_tf->{score};          

          if ($new_score <= $R_THRES) { last; }

          $pos1 = $pos1 + $inc1;
          $pos2 = $pos2 + $inc2;
          $j++;
        }
      }
    }

    if ($k==1) {					# Mutate 2, keep 1.
 
      my($iter1) = $m_length;				# Parameters or 1st iteration.
      my($pos, $inc);
   
      my(@ORIENT);

      if ($head_l >= $tail_l) { @ORIENT = ('-', '+'); }
      if ($head_l <  $tail_l) { @ORIENT = ('+', '-'); }
          
      if ($m_X_orient eq $ORIENT[0]) { ($pos, $inc) = (0, 1); }
      if ($m_X_orient eq $ORIENT[1]) { ($pos, $inc) = ($m_length-1, -1); }

      my($i)=0;
      while ($i < $iter1) {
            
        substr($m_X_motif, $pos, 1) = substr($s_motif, $pos, 1);

        $new_tf->{motif} = $m_X_motif; 
        $new_tf->get_p_value($mat, $dir);
        my($new_score) = $new_tf->{score};          

        if ($new_score <= $R_THRES) { last; }
      
        $pos = $pos + $inc;					# Increment motif position.
        $i++;							# Increment global iteration.
      }
    }
  }
  
  if($m_X_orient eq '-') { $new_tf->{motif} = rev_complement($m_X_motif); }
  return $new_tf;
}

#--------------------------------------------------------------------------------------#
#	                               DATA BASE				       #
#--------------------------------------------------------------------------------------#

#-------------------------------- Gene_to_promoter sub. -------------------------------

# 15. &gene_to_promoter subroutine to convert gene name to promoter name.

sub gene_to_promoter($$$) {
  my($gene, $c_file, $dir) = @_;

  if (length $gene < 4) { 
    die "\n--> Error: not correct gene input!\n\n"; 
  }
  
  chdir $dir;
  open(INFO, $c_file) || die "--> Cannot open conversion file: $c_file!\n";

  my($prom);
  my($marker);
  GENE: while (<INFO>) {
    if (/$gene/) {					# If gene matches on any line.				
      my(@line) = split(/\s/, $_);			# Parse out corresponding promoter.		
     
      if (substr($line[1],1) eq $gene) {
        $prom = $line[5];

        $marker=1;
        last GENE;
      }
    }
  } 
  close INFO;

  if($marker!=1) { 
    print "\n--> Error: no such gene - $gene\n\n"; 
    exit;
  }
  else { return $prom; }
}

#-------------------------------- Promoter_to_gene sub. -------------------------------

# 16. &promoter_to_gene subroutine to convert promoter name to gene name.

sub promoter_to_gene($$$) {
  my($prom, $c_file, $dir) = @_;

  if (length $prom < 5) { 
    die "\n--> Error: not correct promoter input - $prom\n\n"; 
  }
  
  chdir $dir;
  open(INFO, $c_file) || die "--> Cannot open conversion file - $c_file\n";

  my($gene);
  my($marker);
  PROM: while (<INFO>) {
    if (/$prom/) {					# If gene matches on any line.			
      my(@line) = split(/\s/, $_);			# Parse out corresponding gene.
      $gene = substr($line[1], 1);

      $marker=1;
      last PROM;
    }
  } 
  close INFO;

  if($marker!=1) { return $prom; }
  else { return $gene; }
}

#--------------------------------- Motif_db_nir sub. ------------------------------------

# 17. &motif_db_nir subroutine to construct TF -> motif database from cerevisiae-output.
                                                                                        
sub motif_db_nir($) {
  my($map_dir) = shift;
  
  chdir $map_dir;
  my(@files) = <*.01>;          		    # Step through a list of .01 files.

  my(%TF_hash);					    # TF -> motif database bin.

  my($marker);                                                                                         
  foreach my $w (@files) {
    my(@q) = split(/\./, $w);

    print "\n--> Processing file: $w"; 
     
    open(INFO, $w) || die "--> Cannot open map file: $w\n";   # Open the file for input.

    while (<INFO>) {		       			
      if (/^>/) {     
        my(@m) = split(/\s/, $_);

        my($tf) = new DNA::Fact($q[1]);	    	    # Create new TF object.

        $tf->{orient} = $m[5];		    	    # Assign orientation,
        $tf->{offset} = $m[4];		            # motif offset, 
        $tf->{score} = $m[14];		            # p-value (score),
        $tf->{motif} = $m[7];		  	    # and sequence.	

        chop $m[2];
        $tf->{chromosome} = $m[2];	            # Get chromosome number.
  
        my($prom_name) = substr($m[0],1);	    # Get promoter name.
        $tf->{promoter} = $prom_name;    	    # Assign it to TF object.

        push @{ $TF_hash{$q[1]} }, $tf;	            # Record new TF occurence.
      }
    }
    close INFO;
  }
  return \%TF_hash;   		 	            # Return TF -> motif database.			
}

#--------------------------------- Motif_db_nirx sub. ------------------------------------

# 17x. &motif_db_nirx subroutine to construct TF -> motif database from cerevisiae-output
# mining IME1 and FLO11 promoters binding sites also.
                                                                                        
sub motif_db_nirx($) {
  my($map_dir) = shift;
  
  chdir $map_dir;
  my(@files) = <*.01>;          		    # Step through a list of .01 files.

  my(%TF_hash);					    # TF -> motif database bin.

  my($marker);                                                                                         
  foreach my $w (@files) {
    my(@q) = split(/\./, $w);

    print "\n--> Processing file: $w"; 
     
    open(INFO, $w) || die "--> Cannot open map file: $w\n";   # Open the file for input.
   
    while (<INFO>) {		       			
      if (/^>/) {     
        my(@m) = split(/\s/, $_);

        my($tf) = new DNA::Fact($q[1]);	    	    # Create new TF object.

        $tf->{orient} = $m[5];		    	    # Assign orientation,
        $tf->{offset} = $m[4];		            # motif offset, 
        $tf->{score} = $m[14];		            # p-value (score),
        $tf->{motif} = $m[7];		  	    # and sequence.	

        chop $m[2];
        $tf->{chromosome} = $m[2];	            # Get chromosome number.
  
        my($prom_name) = substr($m[0],1);	    # Get promoter name.
        $tf->{promoter} = $prom_name;    	    # Assign it to TF object.

        push @{ $TF_hash{$q[1]} }, $tf;	            # Record new TF occurence.
        
        if ($m[2] eq 'X' || $m[2] eq 'IX') {	    # Check for IME1 & MUC1 binding sites.
          my($strand) = chop $prom_name;
        
          if ($strand eq 'C') {
            use constant IME1_S => 605569;	    # IME1 global range.
            use constant IME1_E => 604487;

            use constant MUC1_S => 393672;	    # MUC1 global range.
            use constant MUC1_E => 389569;
          
            my(@range) = split(/\D/, $m[3]);
            
            my($p_max) = $range[1];	            # Promoter start coordinate (C).
            my($fact_l) = length $m[7];		    # Get motif length.
            
            my($m_max) = $p_max-$m[4];		    # Get TF start coordinate.
            my($m_min) = $m_max-$fact_l;	    # Get TF end coordinate.
            
            if (IME1_S>=$m_max && $m_min>=IME1_E) {  # In IME1 promoter range.
              $tf->{promoter} = $prom_name;    	     # Assign it to TF object.
              push @{ $TF_hash{$q[1]} }, $tf;        # Record new TF occurence.
            }
            elsif (MUC1_S>=$m_max && $m_min>=MUC1_E) {  # In MUC1 promoter range.
              $tf->{promoter} = $prom_name;    	      
              push @{ $TF_hash{$q[1]} }, $tf;       
            }
          }					    # Loop end - C strand.
        }					    # IX or X chromosome.
      }						    # TF hit.
    }
    close INFO;
  }
  return \%TF_hash;   		 	            # Return TF -> motif database.			
}

#--------------------------------------------------------------------------------------#
#				        FILTER					       #
#--------------------------------------------------------------------------------------#

#------------------------------ Overlap_dup_eliminate sub. ----------------------------

# 18. &overlap_dup_eliminate subroutine to eliminate redundant [also same motif-motif]
# overlap objects in a list.
                                                                                        
sub overlap_dup_eliminate($@) {
  my($tf_pair, @overlaps) = @_;

  my($same)=0;                                      # Eliminate common elements in set.
  OL: foreach my $t (@overlaps) {

  #-------------------------------------------------#
  #         Overlap object list - an array          #
  #-------------------------------------------------#
                                                    #
    my($ms1) = $t->{motif_A}->{motif_start};        #
    my($ms2) = $t->{motif_B}->{motif_start};        #
                                                    #
    my($me1) = $t->{motif_A}->{motif_end};          #
    my($me2) = $t->{motif_B}->{motif_end};          #
                                                    #                       
  #-------------------------------------------------#

  #-------------------------------------------------#
  #     Current overlap object for evaluation       #
  #-------------------------------------------------#
                                                    #
    my($MS1) = $tf_pair->{motif_A}->{motif_start};  #
    my($MS2) = $tf_pair->{motif_B}->{motif_start};  #
                                                    #                                         
    my($ME1) = $tf_pair->{motif_A}->{motif_end};    #
    my($ME2) = $tf_pair->{motif_B}->{motif_end};    #
                                                    #
  #-------------------------------------------------#
                                                                                                         
    if ( ($ms1==$MS1 && $me1==$ME1 && $ms2==$MS2 && $me2==$ME2) ||  # If same overlap object pair
         ($MS1==$MS2 && $ME1==$ME2) ) { $same=1; last OL; }         # or same motif range for both
  }								    # motifs in same object - eliminate.
  return $same;
} 

#--------------------------------------------------------------------------------------#
#				       MUTATOR					       #
#--------------------------------------------------------------------------------------#
                                       
#------------------------------- Get_global_range sub. --------------------------------

# 19. &get_global_range subroutine to assign global range to factor object.   
                                                                                        
sub get_global_range($$$) {
  my($factor, $prom_start, $strand) = @_;

  my($offset) = $factor->{offset};
  my($m_length) = length $factor->{motif};

  if ($strand eq 'W' || $strand eq 'M') {   			                 
    $factor->{motif_start} = $prom_start + ($offset + 1);
    $factor->{motif_end}   = $prom_start + ($offset + $m_length);
  }
  elsif ($strand eq 'C') {   			                 
    $factor->{motif_start} = $prom_start - $offset;
    $factor->{motif_end}   = $prom_start - ($offset + $m_length) + 1;
  }
  else { die "\n--> Error: no such strand - $strand\n\n"; }
}                                                            

#----------------------------------- Flip_range sub. ----------------------------------

# 20. &flip_range subroutine to flip Crick range: [max,min] <=> [min,max].   
                                                                                        
sub flip_range($$) {
  my($min, $max) = @_;
  
  if ($min > $max) { return ($max, $min); }
  else { return ($min, $max); }
}

#--------------------------------------------------------------------------------------#
#				        GENERAL	     			               #
#--------------------------------------------------------------------------------------#

#-------------------------------- Arabic_to_roman sub. --------------------------------

# 21. &arabic_to_roman subroutine to convert an arabic numeral to roman.
                                                                                        
sub arabic_to_roman($) {
  my($a_number) = shift;

  my(@roman) = ("I", "II", "III", "IV", "V", "VI", "VII", "VIII",
                "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI");
                
  my($index) = $a_number-1;
  return $roman[$index];
}

#-------------------------------- Roman_to_arabic sub. --------------------------------

# 22. &roman_to_arabic subroutine to convert a roman numeral to arabic.
                                                                                        
sub roman_to_arabic($) {
  my($r_number) = shift;

  my(%r_a_hash) = ( I => 1, II => 2, III => 3, IV => 4,
                    V => 5, VI => 6, VII => 7, VIII => 8,
                    IX => 9, X => 10, XI => 11, XII => 12,
                    XIII => 13, XIV => 14, XV => 15, XVI => 16, );
    
  return $r_a_hash{$r_number};              
}

#-------------------------------- Get_TF_num sub. -------------------------------------

# 23. &get_TF_num subroutine to retrieve number of distinct TFs from a list.
                                                                                        
sub get_TF_num(@) {
  my(@TF_list) = @_;
  
  my(%seen);
  foreach my $tf (@TF_list) { $seen{ $tf->{factor} }++; }

  return scalar (keys %seen);
}

#------------------------------ Get_prom_name sub. -------------------------------------

# 24. &get_prom_name subroutine to retrieve proper promoter name from affinity matrix.
                                                                                        
sub get_prom_name($) {
  my($prom) = shift;        			    # Pre-promoter name.
     
  my($first) = substr($prom, 0, 1);                 # Get first character of promoter name.

  if ($first eq 'i') {
    $prom = substr($prom, 1); 
    
    if ($prom =~ /[-]/) {
      my($pre, $tail) = split(/\-/, $prom);
                              
      if ($tail =~ /\d/) { $prom = $pre; }
      else { 
        my($t) = substr($tail, 0, 1);
    
        if ($prom =~ /delta/) { $prom = $pre; }
        else { $prom = $pre."-".$t; } 
      }
    }
  }
  else {
    if ($prom =~ /[-]/) {
      my($pre, $tail) = split(/\-/, $prom);

      my($t) = substr($tail, 0, 1);
      
      if ($prom =~ /delta/) { $prom = $pre; }
      else { $prom = $pre."-".$t; } 
    }
  }
  return $prom;					    # Return justified promoter name.
}

#------------------------------ Get_bind_affinity sub. ---------------------------------

# 25. &get_bind_affinity subroutine to retrieve p-value entry for given promoter-TF pair
# from binding affinity matrix.
                                                                                        
sub get_bind_affinity($$$$) {
  my($tf, $prom, $tf_conds, $proms) = @_;           
     
  my(%tf_conds) = %{ $tf_conds };		    # Retrieve TF condition / promoter row hashes.
  my(%proms) = %{ $proms };

  my($cond, %hit);				    # Condition and p-value corresponding to minima.

  my(@p_row);
  if (defined $proms{$prom}) { 
    @p_row = @{ $proms{$prom} };                    # Retrieve promoter row containing affinity entries.

    foreach my $k (keys %{ $tf_conds{$tf} }) {
      my($pos) = $tf_conds{$tf}{$k};
      $hit{$k} = $p_row[$pos];    
    }

    my(@result) = sort { $hit{$a} <=> $hit{$b} } keys %hit; 

    my($mark)=0;
    foreach my $c (@result) {
      if ($hit{$c} ne "NaN") { 
        $cond = $c; 
        $mark=1; last; 
      }
    }
    if ($mark!=1) { $cond = 'NaN'; $hit{$cond} = '[------------]'; }
  }
  else { 
    $cond = '---'; 
    $hit{$cond} = '[------------]';	 	    # Assign - to both parameters if promoter not in affinity matrix.	
  }
  return ($cond, $hit{$cond});
}

#-------------------------------- Get_stringency sub. ---------------------------------

# 26. &get_stringency subroutine to retrieve stringency level of data base.
                                                                                        
sub get_stringency($) {
  my($d_file) = shift;           

  my(@first) = split(/\./, $d_file);
  my($str) = chop $first[0];

  if ($str =~ /\d/) { return $str; }
  else { return 0; }
}

#-------------------------------- Get_description sub. -------------------------------

# 27. &get_description subroutine to create gene <-> description library.

sub get_description($$) {
  my($dir, $c_file) = @_;

  chdir $dir;
  open(INFO, $c_file) || die "--> Cannot open conversion file - $c_file\n";

  my(%descr_hash);			# Description hash.			

  my($gene);
  while (<INFO>) {			
    my(@line) = split(/\t/, $_);	# Parse out corresponding gene.
    $gene = substr($line[1], 1);

    chomp($descr_hash{$gene} = $line[7]);
  } 
  close INFO;

  return \%descr_hash;
}

#-------------------------------- Get_tf_cond_list sub. -------------------------------

# 28. &get_tf_cond_list subroutine to create TF_condition list for promoter object.

sub get_tf_cond_list(@) {
  my(@tf_cond_array) = @_;

  my($tf_cond_list);
  foreach my $m (sort { $a cmp $b } @tf_cond_array) {
    $tf_cond_list = $tf_cond_list.$m.',';
  }
  chop $tf_cond_list;

  return $tf_cond_list;
}

#---------------------------------- Get_occurence sub. --------------------------------

# 29. &get_occurence subroutine to calculate TF occurence bound to promoter.

sub get_occurence(@) {
  my(@tf_cond_array) = @_;

  my(%tf_occurs);
  foreach my $m (@tf_cond_array) {
    my(@tf_cond) = split(/\_/, $m);
  
    $tf_occurs{$tf_cond[0]} = $m;
  }
  my($occurs) = scalar (keys %tf_occurs);

  return $occurs;
}

#------------------------------ Prom_to_gene_library sub. -----------------------------

# 30. &prom_to_gene_library subroutine to create promoter <-> gene  library.

sub prom_to_gene_library($$) {
  my($dir, $c_file) = @_;

  chdir $dir;
  open(INFO, $c_file) || die "--> Cannot open conversion file - $c_file\n";

  my(%p_g_hash);			# Promoter-to-gene hash.			

  my($gene, $prom);
  while (<INFO>) {			
    my(@line) = split(/\t/, $_);	# Parse out corresponding gene.
    $gene = substr($line[1], 1);
    $prom = $line[5];
    
    $p_g_hash{$prom} = $gene;
  } 
  close INFO;

  return \%p_g_hash;
}

#---------------------------------- TF_list sub. -------------------------------------

# 31. &TF_list subroutine to get list of all S. cerevisiae TF names.
                                                                                        
sub TF_list($) {
  my($pssm_dir) = shift;

  my(@TF_list);			         # Promoter list holder array.

  chdir $pssm_dir;
  my(@p_files) = <*.pssm>;

  foreach my $f (sort {$a cmp $b} @p_files) {
    my(@fact) = split(/\./, $f);
    push @TF_list, $fact[0];
  }
  return \@TF_list;			 # Return TF list hash reference.
}

#------------------------------- Create_factor sub. ----------------------------------

# 32. &create_factor subroutine to assign attributes to factor object.
                                                                          
sub create_factor($$) {
  my($tf_line, $prom_obj) = @_;

  my(@m) = @{ $tf_line };

  my($orient) = chop $m[0];			      # Get orientation.
  my($factor) = new DNA::Fact($m[0]);		      # Create TF object. 

  $factor->{offset} = $m[1];		       	      # Parse out field values.
  $factor->{orient} = $orient;
  $factor->{score} = $m[2]; 
  $factor->{motif} = $m[6];

  $factor->{binding} = $m[3];
  $factor->{condition} = $m[4];
  $factor->{conservation} = $m[5];

  $factor->{chromosome} = $m[7]; 
  $factor->{promoter} = $m[8];

  if (defined $prom_obj) {
    my($prom_start) = $prom_obj->{prom_start};        # Promoter start.
    my($strand) = $prom_obj->{strand};	 	      # Promoter strand.

    &DNA::Func::get_global_range($factor, $prom_start, $strand);
  }

  return $factor;
}

#--------------------------------- Get_org_file sub. ----------------------------------

# 33. &get_org_file subroutine to retrieve file of specific yeast species.
                                                                          
sub get_org_file($$) {
  my($dir, $org) = @_;

  chdir $dir;
  my(@files) = <*.*>;

  foreach my $file (@files) {
    if ($file =~ /$org/) { return $file; }
  }
}

#------------------------------- Get_conditions sub. ----------------------------------

# 34. &get_conditions subroutine to retrieve binding affinity conditions. 
                                                                                        
sub get_conditions() {
  my($self) = shift;

  my($base_dir) = $self->{parameters}{Base_dir}; 
  my($m_file) = $self->{parameters}{Affinity_file}; 

  chdir $base_dir;
  if (-f $m_file) { 
    open(INFO, $m_file) || die "--> Cannot open binding affinity file: $m_file\n"; 
  }
    
  my(@conds);       
  F: while (<INFO>) {
    chomp(my($row_1) = $_);
    @conds = split(/\,/, $row_1);        		# Put matrix column heads into an array, @row_1.
    last F;
  }

  my(%tf_conds);
  foreach my $c (sort { $a cmp $b } @conds) {
    my($tf, $cond) = split(/\_/, $c);

    if (!exists $tf_conds{$cond}) { $tf_conds{$cond} = 1; }
  }

  $self->{tf_conds} = \%tf_conds;
}

#--------------------------------------------------------------------------------------#
#				      SELECTION	     			               #
#--------------------------------------------------------------------------------------#

#-------------------------------- Binding_filter sub. ----------------------------------

# 35. &binding_filter subroutine to filter by binding affinity p-value. 
                                                                                        
sub binding_filter($$) {
  my($factors, $B_THRESH) = @_;

  my(@filtered)=();				# Array container for filtered TFs.

  for (my $i=0; $i<@$factors; $i++) {

    my($tf) = $factors->[$i]; 
    my($bind) = $tf->{binding};

    if ($bind <= $B_THRESH && !($bind =~/^\[/)) {  # If binding p-value is under a given value.
      push @filtered, $tf;			# Store TF in a list. 
    }
  }
  return @filtered;				# Return filtered TF list.
}

#------------------------------- Condition_filter sub. --------------------------------

# 36. &condition_filter subroutine to filter by binding condition. 
                                                                                        
sub condition_filter($$) {
  my($factors, $COND) = @_;

  my(@filtered)=();				# Array container for filtered TFs.

  for (my $i=0; $i<@$factors; $i++) {
    my($tf) = $factors->[$i];

    if ($tf->{condition} =~ /\b$COND\b/) {	# If binding condition of interest.
      push @filtered, $tf;			# Store TF in a list. 
    }
  }
  return @filtered;				# Return filtered TF list.
}

#------------------------------ Conservation_filter sub. ------------------------------

# 37. &conservation_filter subroutine to filter by evolutionary conservation. 
                                                                                        
sub conservation_filter($$) {
  my($factors, $CONS) = @_;

  my(@filtered)=();				# Array container for filtered TFs.

  for (my $i=0; $i<@$factors; $i++) {
    my($tf) = $factors->[$i];

    if ($tf->{conservation} >= $CONS) {		# If conserved in at least # many times in other yeast organisms.
      push @filtered, $tf;			# Store TF in a list. 
    }
  }
  return @filtered;				# Return filtered TF list.
}

#------------------------------ Probability_filter sub. -------------------------------

# 38. &probability_filter subroutine to filter by motif probability based on PSSM. 
                                                                                        
sub probability_filter($$) {
  my($factors, $PROB) = @_;

  my(@filtered)=();				# Array container for filtered TFs.

  for (my $i=0; $i<@$factors; $i++) {
    my($tf) = $factors->[$i];

    if ($tf->{score} >= $PROB) {		# If motif probability is above a given value.
      push @filtered, $tf;			# Store TF in a list. 
    }
  }
  return @filtered;				# Return filtered TF list.
}
