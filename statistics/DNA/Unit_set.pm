#-----------------------------------------------------------
# Author: Mirko Palla.
# Date: March 2, 2006.
# For: Task 3 at the Church Lab - Genetics Department, 
# Harvard Medical School.
# 
# Purpose: This program contains the complete code for 
# class Unit_set, modeling a set of motif unit objects 
# in Perl.
#
# Associated subroutines: 0. &constructor   1. &mode_1
#			  2. &mode_2 	    3. &mode_3	
#			  4. &remove_motif  5. &remove_units
#------------------------------------------------------------ 

package DNA::Unit_set;

use strict;
use base ("DNA::Unit", "DNA::Mutant", "DNA::IO");

1;

#------------------------------------ Class Unit_set --------------------------------------
	
# 0. constructor for motif unit set object.

sub new {
  my($class, $prom) = @_;
  my($self)={};

  $self->{prom_obj} = $prom; 		                # Promoter object field.

  $self->{units} = ();					# TF of unit storage.
  $self->{removed} = [];				# Removed unit (object) storage.

  $self->{bar_code} = ();                               # Bar code field to label output files.
  $self->{parameters} = {};                             # Input parameters.

  bless($self,$class);
  return($self);
}

#----------------------------------- Running mode 1. ----------------------------------

# 1. &mode_1 subroutine to define motif unit as every single instance of a motif.
                                                                                        
sub mode_1($) {
  my($self, $factors) = @_;

  my(@tf_list) = @{ $factors->{final_list} };

  my($u)=0;
  foreach my $tf (@tf_list) {
    my($unit) = new DNA::Unit($tf);                         # Create motif unit object.
    $unit->{unit_no} = $u;				    # Label unit number.
    push(@{ $self->{units} }, $unit);			    # Units stored in an array.
    $u++;
  }
}

#----------------------------------- Running mode 2. ----------------------------------

# 2. &mode_2 to define motif unit as all copies of a motif.
                        
sub mode_2($) {
  my($self, $factors) = @_;

  my(@tf_list) = @{ $factors->{final_list} };
  my(%units);                                                 # Unit container, a hash.

  foreach my $tf (@tf_list) {
    my($key) = $tf->{factor};				      # Get current TF factor.

    if (exists $units{$key}) {				      # If already in hash.
      push( @{ $units{$key}->{factors} }, $tf); 	      # Put in the same unit.
    }
    else { $units{$key} = new DNA::Unit($tf); }		      # Else create new unit.
  }

  my($u)=0;
  foreach my $key (sort { $units{$a}->{factors}[0]->{offset} <=> 
		          $units{$b}->{factors}[0]->{offset} } keys %units) {

    $units{$key}->{unit_no} = $u; 
    $u++; 
  }
  $self->{units} = \%units;				      # Units stored in hash.
}

#----------------------------------- Running mode 3. ----------------------------------

# 3. &mode_3 to define motif unit as all copies of a motif in distance 'd' from 
# eachother.
                        
sub mode_3($) {
  my($self, $factors) = @_;

  my($dis) = $self->{parameters}{Distance};

  my(@tf_list) = @{ $factors->{final_list} };
  my(%units);                                           # Unit container, a hash.

  foreach my $tf (@tf_list) {
    my($key) = $tf->{factor};				# Get current TF data.
    my($offset) = $tf->{offset};

    my($in_unit);					# Declare marker for "in unit".

    if (exists $units{$key}) {				# If TF already in hash.
      foreach my $unit (@{ $units{$key} }) {
        my($last_tf) = $unit->{factors}[-1];		# Get TF with highest offset.
    
        if (abs($offset-$last_tf->{offset})<=$dis) {    # If in range, put it in the
          push( @{ $unit->{factors} }, $tf);      	# unit, and iterate.
        
          $in_unit = 1;					# Set marker.
          last;						# Get out of unit iteration.
        }
      }
      if (!defined($in_unit)) { 			# If not in any units.
        push(@{ $units{$key} }, new DNA::Unit($tf));    # Create new unit.
      }
    }
    else { 						# If TF not in hash yet.
      push(@{ $units{$key} }, new DNA::Unit($tf)); 	# Create new unit.
    }      						
  }

  my(@unit_list);
  foreach my $key (keys %units) {
    foreach my $unit (@{ $units{$key} }) { push @unit_list, $unit; }
  }   

  $self->{unit_num} = scalar @unit_list;

  my(@sorted) = sort { $a->{factors}[0]->{offset} <=>  
                       $b->{factors}[0]->{offset} } @unit_list;  

  my(%ordered_units);					# Initialize finalized unit hash.

  my($u)=0;
  foreach my $unit (@sorted) { 
    $unit->{unit_no} = $u; 
    my($key) = $unit->{factors}[0]->{factor}; 

    push @{ $ordered_units{$key} }, $unit;    
    $u++; 
  } 
  $self->{units} = \%ordered_units;			# Units stored in hash.
}

#----------------------------------- Remove_motif -------------------------------------

# 4. &remove_motif to remove all copies of TF motif occuring more than once. 

sub remove_motif($) {
  my($self, $mode) = @_;

  my($units) = $self->{units};	  			  # Assign units hash to units.  
  $self->{mode} = $mode;				  # Remember 'permute' or 'randomize' selection.

  my($p_seq) = $self->{prom_obj}->{prom_seq};		  # Original promoter sequence.
  my($bp_diff) = $self->{parameters}{BP_diff};

  foreach my $key (sort keys %$units) {
    my($number) = scalar @{ $units->{$key}->{factors} };

    my($seq) = $self->{prom_obj}->{prom_seq};
    my($kmer_l) = length $units->{$key}->{factors}[0]->{motif};

    if ($number > 1) { 				    	  # If more than 1 TFs present.

      my($mutant) = new DNA::Mutant();                    # Create Mutant object.
      $mutant->{unit} = $key;				  # Assign belonging of Mutant object.
      $mutant->{no_elements} = $number;			  # Assign number of TF in unit.
            
      foreach my $tf (@{ $units->{$key}->{factors} }) {

        my($offset) = $tf->{offset};
        my($kmer) = substr($p_seq, $offset, $kmer_l);     # Kmer for substitution.
        
        my($s_kmer);
  
        if ($mode==0) { 
          $s_kmer = &DNA::Func::permute($kmer, $bp_diff); # Permuted kmer.
        }
        
        if ($mode==1) { 
          my($random) = int(rand(10));			  # Generates a random number: [0,9].
          my(%r_kmer) = DNA::IO::R_KMER;
          $s_kmer = $r_kmer{$kmer_l}[$random]; 		  # Retrieves a randomized kmer from library.
        }

        substr($seq, $offset, $kmer_l) = $s_kmer;         # Insert kmer back to promoter.

        push @{ $mutant->{o_kmer} }, $kmer;               # Store original kmer.
        push @{ $mutant->{s_kmer} }, $s_kmer;             # Store substitition kmer.
        push @{ $mutant->{s_offset} }, $offset;	          # Store substitution offset.
        push @{ $mutant->{kmer_l} }, $kmer_l;             # Store kmer length.
      }

      $mutant->{mutant} = $seq;				  # Store mutagenized sequence.
      push @{ $self->{removed} }, $mutant;		  # Store unit removal history Unit_set.
    }
  }
}

#----------------------------------- Remove_units -------------------------------------

# 5. &remove_units to remove (permutate/randomize) all selected units of promoter. 

sub remove_units($mode, $type) {
  my($self, $mode, $type) = @_;

  $self->{mode} = $type;				     # Remember 'permute' or 'randomize' selection.

  my($units) = $self->{units};   			     # Assign units hash to units.  
  my($p_seq) = $self->{prom_obj}->{prom_seq};		     # Original promoter sequence.

  my($base_dir) = $self->{parameters}{Base_dir};
  my($bp_diff) = $self->{parameters}{BP_diff};

  my($unit_file) = $self->{parameters}{Unit_file};	     # Retrieve unit file. 
  
  chdir $base_dir;  
  open(IN, $unit_file) || die "--> Cannot open unit_input file: $unit_file\n";
  
  my($line)=0;  
  while (<IN>) {
    my(@unit_list) = split(/\s/, $_);

    my($mark)=0;
    foreach my $int (@unit_list) {
      if (!($int =~ /^\d+$/)) { 
        print "\n--> Error: $int - not valid integer value in: $unit_file (line $line skipped)\n\n"; 
        $mark==1; 
      }
    }
   
    if ($mark==0) {
      my(@ordered) = sort {$a <=> $b } @unit_list;
      my($seq) = $self->{prom_obj}->{prom_seq};		     # Get promoter sequence.

      my(%num_hash);					     # Unit order hash.
      if ($mode==1) {
        foreach my $unit (@$units) {
          my($unit_no) = $unit->{unit_no};		     # Unit number of current item.
          $num_hash{$unit_no} = $unit;
        }
      }
      elsif ($mode==2) {
        foreach my $tf (keys %$units) {
          my($unit_no) = $units->{$tf}->{unit_no};	     # Current unit number.        
          $num_hash{$unit_no} = $units->{$tf};		     # Store unit in hash by number.
        }
      }                                   
      elsif ($mode==3) {
        foreach my $tf (keys %$units) {
          my($k)=0;
          while ($k < @{ $units->{$tf} }) {
            my($unit_no) = $units->{$tf}[$k]->{unit_no};
            $num_hash{$unit_no} = $units->{$tf}[$k]; 	     # Store unit in hash by number.
            $k++;
          }
        }
      }

      my(%unit_hash);
      foreach my $num (@unit_list) { $unit_hash{$num} = $num_hash{$num}; }

      my($mutant) = new DNA::Mutant();                         # Create Mutant object.
      $mutant->{unit_hash} = \%unit_hash;	 	       # Assign unit hash to mutate.
      $mutant->{unit_selection} = \@ordered;		       # Retrieve unit file. 
  
      foreach my $unit_no (sort { $a <=> $b } keys %unit_hash) {
        my($unit) = $unit_hash{$unit_no};

        my($first) = $unit->{factors}[0];		       # First unit element.
        my($kmer_l) = length $first->{motif};
        my($factor) = length $first->{factor};

        foreach my $tf (@{ $unit->{factors} }) {

          my($offset) = $tf->{offset};
          my($kmer) = substr($seq, $offset, $kmer_l);	       # Kmer for substitution.
               
          my($s_kmer);
          if ($type==0) {
            $s_kmer = &DNA::Func::permute($kmer, $bp_diff);    # Permuted kmer.
          }
          
          if ($type==1) {
            my($random) = int(rand(10));                       # Generates a random number: [0,9].
            my(%r_kmer) = DNA::IO::R_KMER;
            $s_kmer = $r_kmer{$kmer_l}[$random];	       # Retrieves a randomized kmer from lib
          }

          substr($seq, $offset, $kmer_l) = $s_kmer;            # Insert kmer back to promoter.

          push @{ $mutant->{$unit_no}->{o_kmer} }, $kmer;      # Store original kmer.
          push @{ $mutant->{$unit_no}->{s_kmer} }, $s_kmer;    # Store substitition kmer.
          push @{ $mutant->{$unit_no}->{s_offset} }, $offset;  # Store substitution offset.
          push @{ $mutant->{$unit_no}->{kmer_l} }, $kmer_l;    # Store kmer length.
        }
      }
      $mutant->{mutant} = $seq;			             # Store mutagenized sequence.
      $mutant->{input_line} = $line;		             # Store input line number.
      push @{ $self->{removed} }, $mutant;		     # Store unit removal history Unit_set.
    }
    $line++;
  }
}
