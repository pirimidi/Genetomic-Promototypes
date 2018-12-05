#-------------------------------------------------------------------
# Author: Mirko Palla.
# Date: March 2, 2006.
# For: Task 3 at the Church Lab - Genetics Department, 
# Harvard Medical School.
# 
# Purpose: This program contains the complete code for class 
# Fact_list, modeling storage object containing all TF objects 
# bound to a given promoter in Perl.
#
# Associated subroutines: 0. &constructor       1. &get_range_orient 
#			  2. &get_offset_motif  3. &get_p_value 
#			  4. &shift_correct   
#-------------------------------------------------------------------- 

package DNA::Fact;

use strict;
use base ("DNA::Func");

1;

#------------------------------------ Class Fact --------------------------------------

# 0. constructor for transcription factor object.

sub new {
  my($class, $f) = @_;
  my($self)={};

  $self->{factor} = $f; 					# TF field.
  $self->{score} = ();						# Log-value filed.
  $self->{offset} = ();						# Offset field.
  $self->{orient} = ();						# Orientation field.

  $self->{binding} = ();					# Binding affinity field.
  $self->{condition} = ();					# Binding condition field.
  $self->{conservation} = ();					# Evolutionary conservation field.

  $self->{motif} = (); 						# Motif field.
  $self->{motif_start} = ();					# Local (motif) start.
  $self->{motif_end} = ();					# Local (motif) end.
  $self->{promoter} = ();					# Promoter name(s).
  $self->{chromosome} = ();					# Chromosome number of TF.

  $self->{ernest} = ();						# Ernest indicator.
  $self->{nobind} = '-';					# Nobind indicator.
  $self->{p005} = '-';						# p005 indicator.

  $self->{tf_line} = ();					# TF line number.
  $self->{prom_line} = ();					# Promoter line number.

  $self->{tf_wrap} = ();                                        # TF wrap marker.
  $self->{prom_wrap_o} = ();                                    # Promoter wrap-offset.
  $self->{prom_wrap_l} = ();					# Promoter wrap-line number.

  $self->{new_offset} = ();					# TF offset on segment.
  $self->{color} = ();						# TF color.

  $self->{bound_promoters};					# No. of bound promoters.

  bless($self,$class);
  return($self);
}

#--------------------------------------------------------------------------------------#
#				       MUTATOR					       #
#--------------------------------------------------------------------------------------#

#--------------------------------- Range_orient sub. ----------------------------------

# 1. &range_orient subroutine to assign motif range and orientation to factor object.   
                                                                                        
sub range_orient($$$$) {
  my($self, $r1, $r2, $orient, $dir) = @_;
  
  if ($self->{promoter} ne "DNE") {

    my($strand) = substr($self->{promoter}, 6, 1);        # Get promoter end: W|C.
    my(@pssm) = &DNA::Func::get_pssm($self->{factor}, $dir);

    my(@column) = split(/\s/, $pssm[0]);
    my($c_length) = scalar @column;		          # PSSM column length.
    my($m_length) = abs($r1-$r2)+1;	                  # Length of inital range.

    my($diff) = abs($m_length-$c_length);		  # Difference in column/motif length. 

    if ($strand eq 'W') {                                 # Assign global motif range.
    
      if ($diff==0) {  
        $self->{motif_start} = $r1+1;
        $self->{motif_end} = $r2+1;
      }
      else {
        $self->{motif_start} = ($r1+1)+$diff;
        $self->{motif_end} = $r2+1;
      }

      $self->{orient} = $orient;
    }
    else {

      if ($diff==0) {  
        $self->{motif_start} = $r2+1;
        $self->{motif_end} = $r1+1;
      }
      else {
        $self->{motif_start} = ($r2+1)-$diff;
        $self->{motif_end} = $r1+1;
      }

      if ($orient eq '+') { $self->{orient} = '-'; }      # Assign orientation.
      else { $self->{orient} = '+'; }
    }
  }
  else {						  # Deal with 'DNE' promoter.				
    $self->{motif_start} = $r1;
    $self->{motif_end} = $r2;
    $self->{orient} = $orient;
  }
}

#------------------------------- Get_offset_motif sub. --------------------------------

# 2. &get_offset_motif subroutine to calculate offset & motif sequence of TFBS.   
                                                                                        
sub get_offset_motif($) {
  my($self, $prom_obj) = @_;
  
  my($factor) = $self->{factor};			      # Motif name.
  my($orient) = $self->{orient};			      # Motif orientation.

  my($t_start) = $self->{motif_start};                        # Motif start.
  my($t_end) = $self->{motif_end};                            # Motif end.

  my($m_length) = abs($t_start-$t_end) + 1;		      # Motif length.		      

  my($p_seq) = $prom_obj->{prom_seq};                         # Promoter sequence.
  my($p_start) = $prom_obj->{prom_start};	              # Promoter start.

  my($offset);
  if ($p_start <= $t_start) { $offset = $t_start-$p_start; }  # Motif offset.
  else { $offset = $p_start-$t_start; }

  my($motif) = substr($p_seq, $offset, $m_length); 	      # Motif sequence.

  if ($orient eq '-') { $motif = &DNA::Func::rev_complement($motif); }

  $self->{motif} = $motif;
  $self->{offset} = $offset;
}

#----------------------------------- Get_p_value sub. ---------------------------------

# 3. &get_p_value subroutine to calculate p-value for motif using PSSM matrix entries.   
                                                                                        
sub get_p_value($$) {
  my($self, $matrix, $dir) = @_;
  
  my($factor) = $self->{factor};		            # Get TF name.
  my($motif)  = $self->{motif};			            # Get motif sequence.
  my($length) = length $motif;  			    # Calculate motif length.
  
  my(@pssm) = &DNA::Func::get_pssm($factor, $dir);

  my($line_A) = $pssm[0]; my($line_C) = $pssm[1];
  my($line_T) = $pssm[2]; my($line_G) = $pssm[3];
  
  my(@l_A) = split(/\s/, $line_A);          	  # Put pssm values into an array, @l_base.
  my(@l_C) = split(/\s/, $line_C);
  my(@l_T) = split(/\s/, $line_T);
  my(@l_G) = split(/\s/, $line_G);
               
  my($p);
  my($c)=0;
  my($p_val)=1;
  
  while($c<$length) {
    chomp(my($base) = substr($motif, $c, 1));     # Check each base in motif.

    if ($base eq 'A') { $p = $l_A[$c]; }
    elsif ($base eq 'C') { $p = $l_C[$c]; }
    elsif ($base eq 'T') { $p = $l_T[$c]; }
    elsif ($base eq 'G') { $p = $l_G[$c]; }
    
    if ($matrix == 1) { $p = &DNA::Func::log_to_p($p, $base); }
                      
    if($p==0) { $p = 1E-10; }
    $p_val *= $p;                                 # Calculate p-value.
    $c++;
  }
  $self->{score} = $p_val;
}    

#--------------------------------- Shift_correct sub. ---------------------------------

# 4. &shift_correct subroutine to correct motif sequence based on PSSM.   
                                                                                        
sub shift_correct($$$) {
  my($self, $prom_obj, $dir, $mat) = @_;
 
  my(@TF_possib);					      # Strorage hash for new TF possibilties.

  my($factor) = $self->{factor};			      # Motif name.
  my($score) = $self->{score};			      	      # Motif score (p-value).
  my($offset) = $self->{offset};			      # Motif offset.
  my($orient) = $self->{orient};			      # Motif orientation.
  my($motif) = $self->{motif};				      # Motif sequence.     

  my($p_start) = $prom_obj->{p_start};                        # Promoter global start.  
  my($strand) = $prom_obj->{strand};                          # Promoter strand.  
  my($p_seq) = $prom_obj->{prom_seq};                         # Promoter sequence.  

  &DNA::Func::get_global_range($self, $p_start, $strand);     # Get new, global motif range.	

  my($t_start) = $self->{motif_start};                        # Motif start.
  my($t_end) = $self->{motif_end}; 			      # Motif end.
  my($m_length) = abs($t_start-$t_end) + 1;		      # Motif length.

  #--------- Shift right by 1 bp ---------

  my($r_tf) = new DNA::Fact($factor);         	              # Create right shifted TF object. 

  $r_tf->{offset} = $offset+1;				      # Motif offset - right shift.
  my($r_motif) = substr($p_seq, $offset+1, $m_length);        # Motif sequence - right shift.

  if ($orient eq '-') { $r_motif = &DNA::Func::rev_complement($r_motif); }
  $r_tf->{motif} = $r_motif;

  $r_tf->get_p_value($mat, $dir);			      # Get p-value of new motif.

  if ($r_tf->{score}>$score) {
    &DNA::Func::get_global_range($r_tf, $p_start, $strand);   # Get new, global motif range.	
    push @TF_possib, $r_tf;				      # Store it as possibility.
  }

  #--------- Shift left by 1 bp ---------

  my($l_tf) = new DNA::Fact($factor);         	              # Create left shifted TF object. 

  $l_tf->{offset} = $offset-1;				      # Motif offset - left shift.
  my($l_motif) = substr($p_seq, $offset-1, $m_length);        # Motif sequence - left shift.

  if ($orient eq '-') { $l_motif = &DNA::Func::rev_complement($l_motif); }
  $l_tf->{motif} = $l_motif;

  $l_tf->get_p_value($mat, $dir);			      # Get p-value of new motif.

  if ($l_tf->{score}>$score) {
    &DNA::Func::get_global_range($l_tf, $p_start, $strand);   # Get new, global motif range.	
    push @TF_possib, $l_tf;				      # Store it as possibility.
  }

  #--------- Summarize and finish shift ---------

  my($pos) = scalar @TF_possib;				      # Get new motif possibilities.

  if ($pos!=0) {
    @TF_possib = sort { $b->{score} <=> $a->{score} } @TF_possib;
    my($new_tf) = $TF_possib[0];			      # Get new motif from list.

    $self->{offset} = $new_tf->{offset};		      # Update to corrected motif parameters.
    $self->{score} = $new_tf->{score};

    $self->{motif} = $new_tf->{motif};
    $self->{motif_start} = $new_tf->{motif_start};
    $self->{motif_end} = $new_tf->{motif_end};
  }
}
