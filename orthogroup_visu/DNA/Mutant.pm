#-----------------------------------------------------------------
# Author: Mirko Palla.
# Date: February 20, 2006.
# For: Task 3 at the Church Lab - Genetics Department, 
# Harvard Medical School.
# 
# Purpose: This program contains the complete code for 
# class Mutant, modeling an object for mutated promoter sequence
# storage in Perl.
#
# Associated subroutines: 0. &constructor    
#----------------------------------------------------------------- 

package DNA::Mutant;
                                                                                        
use strict;

1;

#------------------------------------ Class Mutant ------------------------------------

# 0. constructor for Mutant, a promoter mutation.

sub new {
  my($class) = @_;  
  my($self)={};

  $self->{o_kmer} = ();                  # Original kmer from promoter (mutate).
  $self->{x_kmer} = ();                  # Completely mutated kmer field.
  $self->{s_kmer} = ();                  # Substitution kmer field.

  $self->{s_offset} = ();                # Offset of substitution.
  $self->{s_score} = ();                 # P-value of substitution.

  $self->{d_kmer} = ();                  # Original kmer from promoter (keep).
  $self->{D_kmer} = ();                  # Substitution kmer field.
  $self->{D_score} = ();                 # P-value of substitution.

  $self->{f_kmer} = ();                  # Intermediate fragment kmer field.
  $self->{i_kmer} = ();                  # Intermediate substitution kmer field.
  $self->{i_offset} = ();                # Intermediate substitution kmer offset.
  $self->{ikmer_l} = ();                 # Intermediate substitution kmer length.

  $self->{mutant} = ();   	         # Mutated promoter field.
  
  $self->{unit} = ();   	         # Unit mutant specification field.
  $self->{no_elements} = 0;   	         # Number of modifications.

  $self->{overlap_pair} = ();  	         # Overlap pair field: [motif_B, motif_A] format.
  $self->{untouched_tf} = ();  	         # List of non-mutated pair elements.

  bless($self,$class);
  return($self);
}
