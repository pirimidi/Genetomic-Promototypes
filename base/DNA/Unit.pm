#-----------------------------------------------------------
# Author: Mirko Palla.
# Date: March 2, 2006.
# For: Task 3 at the Church Lab - Genetics Department, 
# Harvard Medical School.
# 
# Purpose: This program contains the complete code for 
# class Unit, modeling a motif unit object in Perl.
#
# Associated subroutines: 0. &constructor  
#------------------------------------------------------------ 

package DNA::Unit;

use strict;

1;

#------------------------------------ Class Unit --------------------------------------
	
# 0. constructor for motif unit object.

sub new {
  my($class, $tf) = @_;
  my($self)={};

  push(@{ $self->{factors} }, $tf);			# Insert TF into unit.
  $self->{unit_no} = ();			        # Initialize unit number.		

  bless($self,$class);
  return($self);
}
