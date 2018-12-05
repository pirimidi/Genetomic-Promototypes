#-----------------------------------------------------------
# Author: Mirko Palla.
# Date: March 23, 2006.
# For: Task 3 at the Church Lab - Genetics Department, 
# Harvard Medical School.
# 
# Purpose: This program contains the complete code for 
# class Overlap, modeling a motif overlap object in Perl.
#
# Associated subroutines: 0. &constructor  
#------------------------------------------------------------ 

package DNA::Overlap;

use strict;

1;

#---------------------------------- Class Overlap -------------------------------------
	
# 0. constructor for motif overlap [TF pair] object.

sub new {
  my($class, $TF1, $TF2) = @_;
  my($self)={};

  $self->{motif_A} = $TF1;		  # Overlapping motif in position 1.
  $self->{motif_B} = $TF2;		  # Overlapping motif in position 2. 

  $self->{overlap} = ();		  # Overlap length.

  $self->{sign} = ();			  # Overlap magnitude: {+,-} -> {positive, negative}.
  $self->{type} = ();			  # Type of overlap: {"fake","real->I","real->II"}.

  $self->{head_l} = ();			  # Motif head length: only for type II.
  $self->{tail_l} = ();			  # Motif tail length: only for type II.

  bless($self,$class);
  return($self);
}
