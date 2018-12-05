#-------------------------------------------------------------
# Author: Mirko Palla.
# Date: August 7, 2006.
# For: Task 6 at the Church Lab - Genetics Department, 
# Harvard Medical School.
# 
# Purpose: This program contains the complete code for 
# class TF_pair, modeling a motif pair object in Perl.
#
# Associated subroutines: 0. &constructor     1. &orient_type
#			  2. &calculate_dist  
#-------------------------------------------------------------- 

package DNA::TF_pair;

use strict;

1;

#------------------------------------ Class TF_pair -----------------------------------
	
# 0. constructor for motif pair object.

sub new {
  my($class, $TF_pos1, $TF_pos2) = @_;
  my($self)={};

  #      ---> 
  #      |
  #      ---------------[TF_pos1]----[TF_pos2]-------[GENE]: (Watson stand)
  #      ^
  # [Transcription start]
  #

  $self->{position_1} = $TF_pos1;			# TF in position 1, i.e., closer to transcription start.
  $self->{position_2} = $TF_pos2;			# TF in position 2, i.e., closer to gene start.

  $self->{distance} = ();				# Relative TF distance on promoter.
  $self->{orient_type} = ();				# Relative TF orientation type: {++, +-, --, -+}.

  $self->{promoter} = ();				# Promoter name.
  $self->{chromosome} = ();				# Chromosome number of TF pair.

  bless($self,$class);
  return($self);
}

#----------------------------------- Orient_type sub. ---------------------------------

# 1. &orient_type subroutine to determine and assign orientation type of TF pair.

sub orient_type() {
  my($self) = shift;

  my($strand) = $self->{strand};			# Get promoter strand: {W,C}.

  my($orient_pos1) = $self->{position_1}->{orient};	# Get orientation of TFs in order.
  my($orient_pos2) = $self->{position_2}->{orient};

  if ($strand eq 'C') {
    if ($orient_pos1 eq '+') { $orient_pos1 = '-'; }	# Take care of strand relativity.
    else { $orient_pos1 = '+'; }

    if ($orient_pos2 eq '+') { $orient_pos2 = '-'; }
    else { $orient_pos2 = '+'; }
  }

  $self->{orient_type} = $orient_pos1.$orient_pos2;	# Assign orientation type.
}

#---------------------------------- Calculate_dist sub. -------------------------------

# 2. &calculate_dist subroutine to determine relative distance of TF pair.

sub calculate_dist() {
  my($self) = shift;

  my($strand) = $self->{strand};			# Get promoter strand: {W,C}.

  my($tf_L);						# TF on the absolute left.
  my($tf_R);					        # TF on the absolute right.

  my($tf_L_length);					# Lenght of TF on the absolute left.

  if ($strand eq 'W' or $strand eq 'M') {
    $tf_L = $self->{position_1}->{offset};		
    $tf_R = $self->{position_2}->{offset};

    $tf_L_length = length $self->{position_1}->{motif};    
  }
  else {
    $tf_L = $self->{position_2}->{offset};		
    $tf_R = $self->{position_1}->{offset};

    $tf_L_length = length $self->{position_2}->{motif};    
  } 
  
  my($r_dist) = $tf_R - ($tf_L + $tf_L_length);		# Calculate relative distance of TF pair.

  $self->{distance} = $r_dist;				# Assign distance parameter to object.
}

  
