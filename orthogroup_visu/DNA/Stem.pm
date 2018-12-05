#-----------------------------------------------------------
# Author: Mirko Palla.
# Date: February 13, 2006.
# For: Task 3 at the Church Lab - Genetics Department, 
# Harvard Medical School.
# 
# Purpose: This program contains the complete code for 
# class Stem, modeling any Polypromoter object available 
# in Perl.
#
# Associated subroutines: 0. &constructor
#------------------------------------------------------------ 

package DNA::Stem;

use strict;
use base ( "DNA::Fact", "DNA::Fact_list", "DNA::Mutant",   "DNA::Mutant_set",
           "DNA::Unit", "DNA::Unit_set",  "DNA::Promoter", "DNA::Unit", 
           "DNA::Visu" );
1;

#------------------------------------ Class Stem --------------------------------------

# 0. constructor for any Polypromoter object.

sub new {
  my($class) = @_;
  my($self)={};

  bless($self,$class);
  return($self);
}
