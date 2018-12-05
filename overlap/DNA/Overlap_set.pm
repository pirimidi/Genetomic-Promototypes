#-----------------------------------------------------------
# Author: Mirko Palla.
# Date: March 23, 2006.
# For: Task 3 at the Church Lab - Genetics Department, 
# Harvard Medical School.
# 
# Purpose: This program contains the complete code for 
# class Overlap_set, modeling a set of motif overlap objects 
# in Perl.
#
# Associated subroutines: 0. &constructor     1. &overlap
#			  2. &remove_overlap 
#------------------------------------------------------------ 

package DNA::Overlap_set;

use strict;
use base ("DNA::Overlap", "DNA::Mutant", "DNA::IO");

1;

#------------------------------- Class Overlap_set ------------------------------------
	
# 0. constructor for motif overlap set object.

sub new {
  my($class, $fact_list) = @_;
  my($self)={};

  $self->{fact_list} = $fact_list;	                # Motif list field.
  $self->{prom_obj} = $fact_list->{prom_obj};           # Promoter object field.

  $self->{bar_code} = ();                               # Bar code field to label output files.
  
  $self->{overlaps} = ();				# Overlap [TF pair] storage.
  $self->{removed} = [];				# Removed overlap (object) storage.

  $self->{pssms} = [];   				# PSSMS matrix storage.
  $self->{parameters} = {};                             # Input parameters.

  bless($self,$class);
  return($self);
}

#----------------------------------- Overlap sub. -------------------------------------

# 1. &overlap to define motif overlap as two motifs in distance 'overlap_dis' 
# from eachother.
                        
sub overlap() {
  my($self) = shift;
 
  my($dis) = $self->{parameters}{Overlap_gap};
  my($strand) = $self->{prom_obj}->{strand};			 # Get promoter strand.

  my(@tf_list) = @{ $self->{fact_list}->{final_list} };

  my(@overlaps);						 # Array container for overlap object storage.

  for (my $i=0; $i<@tf_list; $i++) {			         # Iterate through all TF pairs in promoter.
    my($TF_A) = $tf_list[$i];				         # Fixed TF in list, base of comparision.

    my($m_A_start) = $TF_A->{motif_start};		         # Get fixed TF motif start.
    my($m_A_end) = $TF_A->{motif_end};                           # Get fixed TF motif end.
  
    for (my $j=$i+1; $j<@tf_list; $j++) {
      my($TF_B) = $tf_list[$j];		         	         # TF of iterative comparision.

      my($m_B_start) = $TF_B->{motif_start};		         # Get TF motif start to compare.
      my($m_B_end) = $TF_B->{motif_end};                         # Get TF motif end to compare.

      my($range);
      if ($strand eq 'W') { $range = $m_A_end + $dis; }
      else { $range = $m_A_end - $dis; }

      if (($strand eq 'W' && $m_B_start <= $range) || 
          ($strand eq 'C' && $m_B_start >= $range)) { 
        
        my($tf_pair) = new DNA::Overlap($TF_A, $TF_B); 
          
        #-------------------------------------------------#
        #     Current overlap object for evaluation       #
        #-------------------------------------------------#
                                                          #
        my($motif_A) = $tf_pair->{motif_A};	          # 
        my($motif_B) = $tf_pair->{motif_B};		  #
                                                          #
        my($ms1) = $motif_A->{motif_start};	 	  #
        my($me1) = $motif_A->{motif_end};	 	  #
                                                          #
        my($ms2) = $motif_B->{motif_start};		  #  
        my($me2) = $motif_B->{motif_end};	 	  #
                                                          #
        #-------------------------------------------------#    
    
        my($overlap);			 
        if ($strand eq 'W' || $strand eq 'M') {		             # Watson strand.	       
          if ($me1 <= $me2) { $overlap = $me1-$ms2+1; }
          else { 
            $overlap = $me2-$ms2+1; 
            $tf_pair->{head_l} = $ms2-$ms1;			     # Motif head length.
            $tf_pair->{tail_l} = $me1-$me2;			     # Motif tail length.
          }
        }  
        else {       				                     # Crick strand.
          if ($me1 >= $me2) { $overlap = $ms2-$me1+1; }
          else { 
            $overlap = $ms2-$me2+1; 
            $tf_pair->{head_l} = $ms1-$ms2;			     # Motif head length.
            $tf_pair->{tail_l} = $me2-$me1;			     # Motif tail length.
          }
        }

        $tf_pair->{overlap} = $overlap;  			     # Assign overlap value.			       
        push @overlaps, $tf_pair; 				     # Store overlap object in an array.
      }
    }
  }
  @{ $self->{overlaps} } = @overlaps;		       		     # Store list of overlap objects in Overlap_set.
}

#--------------------------------- Remove_overlap sub. --------------------------------

# 2. &remove_overlap to remove motif overlap [motif_A vs. motif_B and vice and versa].

sub remove_overlap() {
  my($self) = shift;

  my($strand) = $self->{prom_obj}->{strand};        # Get promoter strand.
  my(%PSSMS);					    # PSSM matrix hash for all TFs in OL list.

  OL: foreach my $o_lap ( @{ $self->{overlaps} } ) {

  #-------------------------------- PSSM retrieval ----------------------------------
    
    my($motif_A) = $o_lap->{motif_A};		    # Motif in 1st position. 	    
    my($motif_B) = $o_lap->{motif_B}; 	   	    # Motif in 2nd position.

    my($m_name1) = $motif_A->{factor}; 	    
    my($m_name2) = $motif_B->{factor}; 	   

    my($mat) = $self->{parameters}{Matrix};	    # Get matrix type.
    my($dir) = $self->{parameters}{PSSMS};	    # Get PSSMS directory.
   
    my(@m_pair) = ($m_name1, $m_name2);

    my($i)=0;
    while ($i < @m_pair) {  			    # Loop twice for each TF.

      if (!exists ($PSSMS{$m_pair[$i]})) { 	    # If PSSM matrix not in hash for TF.

        my($motif) = $m_pair[$i];
        my($file) = $motif.".pssm";

        chdir $dir;				    # Change to pssms directory. 
        if (-f $file) {
          open(INFO, $file) || die "--> Cannot open PSSMS file: $file\n";

          my(@pssm);
          while(<INFO>) { push(@pssm, $_); }
          close INFO;

          my($line_A) = $pssm[0]; my($line_C) = $pssm[1];
          my($line_T) = $pssm[2]; my($line_G) = $pssm[3]; 

          my(@l_A) = split(/\s/, $line_A);          # Put pssm values into an array, @l_base.
          my(@l_C) = split(/\s/, $line_C);
          my(@l_T) = split(/\s/, $line_T);
          my(@l_G) = split(/\s/, $line_G);

          my(%pssm) = ( A => \@l_A, C => \@l_C,     # Position weight matrix, a hash.
                        G => \@l_G, T => \@l_T );

          $PSSMS{$motif} = \%pssm;                  # Store PSSM in a storage hash.
        }
        else { 
          print "\n--> No PSSMS matrix for factor: $motif\n\n"; 
          next OL;
        }
      }
      $i++
    }

    my($overlap) = $o_lap->{overlap};		    	       # Get overlap value. 	    
    my(@m_pairs) = ($motif_A, $motif_B);
    
    my($k)=0;
    REP: while ($k < @m_pairs) {  

      my($m_X) = $m_pairs[$k];		                       # Get motif to mutate.
      my($m_X_name) = $m_X->{factor};                          # Get motif name ~.
      my($m_X_offset) = $m_X->{offset};           	       # Get motif offset ~.
      my($m_X_orient) = $m_X->{orient};                	       # Get motif orientation ~.
      my($m_X_motif)  = $m_X->{motif};                 	       # Get motif sequence ~.
      my($m_X_score)  = $m_X->{score};                 	       # Get motif sequence ~.
    
      my($prom_seq) = $self->{prom_obj}->{prom_seq};	       # Get promoter sequence.
      my($m_X_length) = length $m_X_motif;		       # Motif length to delete.
      my($o_motif) = substr($prom_seq, $m_X_offset, $m_X_length);           

      my($m_D);						       # Get motif to be dominated. 
      if ($k==0) { $m_D = $m_pairs[1]; }     
      else       { $m_D = $m_pairs[0]; }     

      my($m_D_name) = $m_D->{factor};                          # Get motif name to dominate.
      my($m_D_offset) = $m_D->{offset};           	       # Get motif offset ~.
      my($m_D_orient) = $m_D->{orient};                	       # Get motif orientation ~.
      my($m_D_motif)  = $m_D->{motif};                 	       # Get motif sequence ~.
      my($m_D_score)  = $m_D->{score};                 	       # Get motif sequence ~.
    
      my($m_D_length) = length $m_D_motif;		       # Motif length to dominate.
      my($d_motif) = substr($prom_seq, $m_D_offset, $m_D_length);           

      foreach my $m (@{ $self->{removed} }) {                  # Remove repeating promoter mutations.
        if ($overlap <= 0 && 
            $m->{overlap_pair}->{sign} eq '-') { 

          my($m_tf) = $m->{overlap_pair}->{motif_B};	       # Motif to mutate (already in list).
          my($m_name) = $m_tf->{factor};		       # Motif to mutate name.
          my($m_offset) = $m_tf->{offset};	               # Motif to mutate offset.

          if ($m_X_name eq $m_name && 			       # If same occurence.
                         $m_X_offset == $m_offset) {
            push @{ $m->{untouched_tf} }, $m_D;		       # Remember motif by assigning to 'm'.
       
            $k++;  next REP;
          }
        }
      }  

      my($mutant) = new DNA::Mutant();                         # Create Mutant object.

      my($s_motif);					       # Substitution motif.
      my($sign);					       # Overlap type: {+,-}.
      my($type);					       # Overlap type: {"fake","real->I","real->II"}.

      if ($overlap <= 0) { 
        $s_motif = &DNA::Func::pssm_mutate($PSSMS{$m_X_name}, $m_X);

        $sign = '-';
        $type = "fake";
      }	
      else { 					               # Positive (real) overlap.

        my($me1) = $motif_A->{motif_end};  	 	       # Get 1st motif global end.
        my($me2) = $motif_B->{motif_end};  	      	       # Get 2nd motif global end.

        my($mat) = $self->{parameters}{Matrix};		       # Matrix type.
        my($keep_perc) = $self->{parameters}{Keep_perc};       # Maximun theshold to keep (in %).

        #----------------------------------------------------------------------------
        #
        # Sub-case I:                 ms2|ACTTTGTA|TGGGTAAAATG|me2
        #		                  ********		
        #		 ms1|ATGTTTGTGCTC|ACTTTGTA|me1			
        #
        #----------------------------------------------------------------------------

        if (($strand eq 'W' && $me1<=$me2) || 	               # Watson strand.	       
            ($strand eq 'C' && $me1>=$me2)) {	               # Crick strand.	       

          my($x_motif) = &DNA::Func::pssm_mutate($PSSMS{$m_X_name}, $m_X);
          
          my($offset_1) = $motif_A->{offset};
          my($length_1) = length $motif_A->{motif};
          my($motif_1) = substr($prom_seq, $offset_1, $length_1);           
          
          my($F_motif) = substr($motif_1, -$overlap);  		       # Original overlap motif. 
          my($I_motif) = &DNA::Func::overlap_mutate($o_lap,            # Deleted overlap "motif".
              $k, $F_motif, $mat, $keep_perc, \%PSSMS);

          my($I_offset) = $offset_1 + ($length_1 - $overlap);

          my($f_motif); 
          my($i_motif);					       # Sub. motif fragment.
          my($i_offset);                                       # Sub. motif fragment offset. 

          my($fr_length) = $m_X_length - $overlap;

          if ($k==0) { 
            $f_motif = substr($o_motif, 0, $fr_length); 
            $i_motif = substr($x_motif, 0, $fr_length); 
            $i_offset = $m_X_offset;
            $s_motif = $i_motif.$I_motif; 
          }
          elsif ($k==1) { 
            $f_motif = substr($o_motif, $overlap); 
            $i_motif = substr($x_motif, $overlap); 
            $i_offset = $m_X_offset + $overlap; 
            $s_motif = $I_motif.$i_motif; 
          }
          else { 
            print "\n--> Invalid motif iteration number: $k\n\n";
            exit;
          } 

          $sign = '+';
          $type = "real->I";
     
          @{ $mutant->{f_kmer} } = ($f_motif, $F_motif);     	   # Fragment kmer field.
          @{ $mutant->{i_kmer} } = ($i_motif, $I_motif);   	   # Intermediate kmer field.
          @{ $mutant->{i_offset} } = ($i_offset, $I_offset);       # Intermediate offset of substitution.
          @{ $mutant->{ikmer_l} } = ($fr_length, $overlap);        # Store intermediate kmer length.
        }

        #----------------------------------------------------------------------------
        #
        # Sub-case II:                  ms2|ACTTTGTA|me2
        #		                    ********		
        #		     ms1|ATGTTTGTGC|ACTTTGTA|GTG|me1			
        #
        #----------------------------------------------------------------------------

        if (($strand eq 'W' && $me1>$me2) || 	               # Watson strand.	       
            ($strand eq 'C' && $me1<$me2)) {	               # Crick strand.	       
          
          if ($k==0) {					       # Mutate 1, keep 2 case.
            my($x_motif) = &DNA::Func::pssm_mutate($PSSMS{$m_X_name}, $m_X);
          
            my($length_2) = length $m_D_motif;
            my($motif_2) = substr($prom_seq, $m_D_offset, $length_2);  # Original overlap.       

            my($head_l) = $o_lap->{head_l};
            my($tail_l) = $o_lap->{tail_l};
          
            my($o_head) = substr($o_motif, 0, $head_l);
            my($o_tail) = substr($o_motif, -$tail_l);

            my($x_head) = substr($x_motif, 0, $head_l);
            my($x_tail) = substr($x_motif, -$tail_l);
          
            my($t_offset) = ($m_X_offset + $m_X_length) - $tail_l;  
          
            my($I_motif) = &DNA::Func::overlap_mutate($o_lap,          # Deleted overlap "motif".
                $k, $motif_2, $mat, $keep_perc, \%PSSMS);

            $s_motif = $x_head.$I_motif.$x_tail;           

            @{ $mutant->{f_kmer} } = ($o_head, $motif_2, $o_tail);     # Fragment kmer field.
            @{ $mutant->{i_kmer} } = ($x_head, $I_motif, $x_tail);     # Intermediate kmer field.
            @{ $mutant->{i_offset} } = ($m_X_offset, $m_D_offset,      # Intermediate offset of substitution.
                                                       $t_offset);
            @{ $mutant->{ikmer_l} } = ($head_l, $overlap, $tail_l);    # Store intermediate kmer length.
          }
          if ($k==1) {
            $s_motif = &DNA::Func::overlap_mutate($o_lap,              # Deleted overlap "motif".
                $k, $o_motif, $mat, $keep_perc, \%PSSMS);
          }

          $sign = '+';
          $type = "real->II";
        }
      }

      #------------------------ MUTATED MOTIF PARAMETERS -----------------------------
      
      my($new_tf) = &DNA::Func::overlap_motif($o_lap, $s_motif, $k, $type, $self);

      my($final_motif) = $new_tf->{motif};
      my($final_score) = $new_tf->{score};

      substr($prom_seq, $m_X_offset, $m_X_length) = $final_motif;
      
      #----------------------- DOMINANT MOTIF PARAMETERS -----------------------------
      
      my($keep_motif) = substr($prom_seq, $m_D_offset, $m_D_length);           
      $mutant->{D_kmer} = $keep_motif;	                       # Substitution kmer field (keep).

      if($m_D_orient eq '-') { $keep_motif = &DNA::Func::rev_complement($keep_motif); }

      my($D_tf) = new DNA::Fact($m_D_name);

      $D_tf->{motif} = $keep_motif;
      $D_tf->get_p_value($mat, $dir);
      my($D_score) = $D_tf->{score};
 
      #------------------------- MUTANT OBJECT SET-UP ---------------------------------

      $mutant->{overlap_pair} = new DNA::Overlap($m_D, $m_X);  # Overlap pair: [dominate, mutate] format.

      $mutant->{overlap_pair}->{sign} = $sign;		       
      $mutant->{overlap_pair}->{type} = $type;

      $mutant->{o_kmer} = $o_motif;	                       # Original kmer from promoter (mutate).
      $mutant->{s_kmer} = $final_motif;                        # Substitution kmer field.
      $mutant->{s_score} = $final_score;	               # Substitution kmer p-value.

      $mutant->{x_kmer} = $s_motif;                            # Completely mutated kmer field.

      $mutant->{d_kmer} = $d_motif;	                       # Original kmer from promoter (keep).
      $mutant->{D_score} = $D_score;	                       # Substitution kmer p-value.
      
      $mutant->{overlap_pair}->{overlap} = $overlap;           # Overlap value.
      $mutant->{mutant} = $prom_seq;               	       # Mutated promoter sequence.

      push @{ $self->{removed} }, $mutant;		       # Store "overlap-removed" Mutant object.

      $k++;
    }
  }
}
