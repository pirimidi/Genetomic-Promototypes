#---------------------------------------------------------------------------
# Author: Mirko Palla.
# Date: February 22, 2006.
# For: Task 3 at the Church Lab - Genetics Department, 
# Harvard Medical School.
# 
# Purpose: This program contains the complete code for 
# class Visu, modeling an object for promoter visualization 
# in Perl.
#
# Associated subroutines: 0. &constructor         1. &TF_color_rows
#			  2. &TF_align            3. &window_setup
#			  4. &promoter_TF_insert  5. &motif_color_TF_insert
#--------------------------------------------------------------------------- 

package DNA::Visu;

use strict;

1;

#------------------------------------ Class Visu --------------------------------------

# 0. constructor for a promoter visualization object.

sub new {
  my($class, $prom, $fact) = @_;
  my($self)={};

  $self->{prom_obj} = $prom;		 		# Promoter object field.
  $self->{fact_obj} = $fact;		 		# Fact_list object field.
  $self->{prom_set} = {};		 		# Prom_set object field.

  $self->{gene} = (); 					# Gene field.
  $self->{bar_code} = ();			 	# Identification field.

  $self->{window} = ();		 			# main window field.
  $self->{canvas} = ();		 			# Canvas field.
  $self->{tb};						# ROText field.

  $self->{ROW} = ();		                        # Set # of chars in a row.

  $self->{width} = ();		 			# Canvas width field.
  $self->{height} = ();		 			# Canvas height field.
  $self->{font} = ();		 			# Font size field.

  $self->{parameters} = {};                             # Input parameters.

  bless($self,$class);
  return($self);
}

#---------------------------------- TF_color_rows sub. --------------------------------

# 1. &TF_color_rows subroutine to determine TF colors and number of empty rows for TF insertion.
                                                                                        
sub TF_color_rows() {
  my($self) = @_;

  my($c_row);                                     # Initialize number of current TF rows.
  my($f_row)=0;                                   # Initialize number of TF rows.
  my($index)=1;                                   # Initialize color list index.

  my(@ordered) = @{ $self->{fact_obj}->{final_list} };  # Assign final TF list.

  my($TF_num) = &DNA::Func::get_TF_num(@ordered); # Get number of distinct TFs.
  if ($TF_num>20) { $TF_num = 20; } 		  # Maximum color number: 20.

  foreach my $o (@ordered) {                      # Color and TF overlap calculation.
    my($x) = $o->{factor};
    my($y) = $o->{offset};

    my($l) = (length $x)+1;			  # TF name + orientation length.
    $o->{color} = 'c'.$index;                     # Assign color with given index to TF object.
     
    $c_row=0;
    foreach my $u (@ordered) {
      my($i) = $u->{factor};
      my($j) = $u->{offset};
                  
      if ($x eq $i) { $u->{color} = 'c'.$index; }            # If same TF, assign same color.

      my($diff) = $j-$y;
      if ($x ne $i && (0<=$diff && $diff<$l)) { $c_row++; }  # If overlap, but not same TF.
      
      if ($f_row<$c_row) { $f_row=$c_row; }                  # Update distance between promoter rows.
    }
    $index++;
    if ($index>$TF_num) { $index=1; }
  }
  return ($f_row+1);
}
                                
#------------------------------------ TF_align sub. -----------------------------------

# 2. &TF_align subroutine to determine promoter and TF row number.

sub TF_align($) {
  my($self, $f_rows) = @_;
  
  use constant TOP => 3;          	                  # Set gap from top.
  my($gap) = $f_rows+3;                               

  my($ROW) = $self->{parameters}{ROW};

  my(@ordered) = @{ $self->{fact_obj}->{final_list} };    # Assign final TF list.

  foreach my $p (@ordered) {                              # TF alignment calculation, assignment.
    my($f) = $p->{factor};
    my($s) = $p->{offset};

    my($l) = length $p->{motif};
    my($k) = (length $p->{factor})+1;

    my($row_count)=0;
    while (($s-($ROW*$row_count))>0) { $row_count++ };    # No. of rows to display.

    my($prom_line) = TOP+(($row_count-1)*$gap);
    my($tf_line) = (TOP+1)+(($row_count-1)*$gap);
               
    $p->{prom_line} = $prom_line;
    $p->{tf_line} = $tf_line;
                   
    if ($s>$ROW) { $p->{new_offset} = $s%$ROW; }          # New offset on specific row.
    else { $p->{new_offset} = $s; }
                        
    my($row_diff) = $ROW-$p->{new_offset};

    if ($row_diff==$ROW) {  		                  # TF row wrap handling - modulo 0.
      $p->{offset} = $ROW*$row_count;
      $p->{prom_line} = TOP+($row_count*$gap);
      push @{ $p->{tf_wrap} }, $p;	 		  # Mark TF wrap start. 
    }

    if ($row_diff<$l) {			                  # TF row wrap handling - other.
      $p->{prom_wrap_o} = $l-$row_diff;
      $p->{prom_wrap_l} = TOP+($row_count*$gap);
    }

    if ($row_diff<$k) { 
      $p->{offset} = $ROW*$row_count;
      $p->{tf_line} = (TOP+1)+($row_count*$gap);
      push @{ $p->{tf_wrap} }, $p;	 		  # Mark TF wrap start. 
    }
  }

  my($array_l) = scalar @ordered;			  # Get lenght of TF list.
  for (my $i=1; $i<$array_l; $i++) {	                  # TF overlap handling.

    my($tf_1) = $ordered[$i-1];
    my($tf_2) = $ordered[$i];

    my($o1) = $tf_1->{offset};
    my($o2) = $tf_2->{offset};

    if (defined ($tf_1->{tf_wrap}) &&
        defined ($tf_2->{tf_wrap})) {

         @{ $tf_2->{tf_wrap} } = ( @{ $tf_1->{tf_wrap} }, 
                                    @{ $tf_2->{tf_wrap} } );
    }  

    my($tf_l) = (length $tf_1->{factor})+1;      	  # Length of TF name + orientation.
    my($tf_line);

    if (($o2-$o1)<$tf_l) { 
      $tf_line = $tf_1->{tf_line};  
      $tf_2->{tf_line} = $tf_line+1;               	  # Assign updated TF line number.
    }
    elsif (defined $tf_1->{tf_wrap} && !(defined $tf_2->{tf_wrap})) {
      foreach my $m ( @{ $tf_1->{tf_wrap} } ) {
        my($o3) = $m->{offset};
        my($m_l) = (length $m->{factor})+1;      	  # Length of TF name + orientation in list.

        if (($o2-$o3)<$m_l) { 
          $tf_line = $m->{tf_line}; 
          $tf_2->{tf_line} = $tf_line+1;               	  # Assign updated TF line number.
        }
      }
    }
  }
}                    

#--------------------------------- Window_setup sub. ----------------------------------

# 3. &window_setup subroutine to create main window and canvas field for visualization.

sub window_setup() {
  my($self) = @_;

  my($gene) = $self->{prom_obj}->{gene};          # Assign gene.
  my($prom) = $self->{prom_obj}->{promoter};      # Assign promoter.
  my($bar_code) = $self->{fact_obj}->{bar_code};  # Assign bar_code.

  my($width) = $self->{parameters}{Width};        # Assign window width.
  my($height) = $self->{parameters}{Height};      # Assign window height.
  my($font) = $self->{parameters}{Font};          # Assign letter font.

  my($main) = MainWindow->new(-title => "Transcription Factor Binding Site Viewer - $gene : [$prom] - ID#: $bar_code");

  my($canvas) = $main->Scrolled('Canvas', -relief => 'raised', -width => $width,
                                          -height => $height, -cursor => 'top_left_arrow',
                                          -background => 'white', -scrollbars => 'se',
                                          -scrollregion => [-10,-10,720,5000]);

  $canvas->pack(-expand => 'yes', -fill => 'both');
  $main->fontCreate('small', -family => 'courier', -size => $font);
  
  $self->{window} = $main;					# Assign window value to Visu.
  $self->{canvas} = $canvas;					# Assign canvas value to Visu.
}                                                                          

#------------------------------ Promoter_TF_row_insert sub. ---------------------------

# 4. &promoter_TF_row_insert subroutine to insert promoter sequences and corresponding TFs.

sub promoter_TF_row_insert($$) {
  my($self, $f_rows, $canvas) = @_;

  my($prom_seq) = $self->{prom_obj}->{prom_seq};
  my($ROW) = $self->{parameters}->{ROW};

  #----------------------------- TF alignment parameters ------------------------------
                                                                                      #
  use constant TOP => 3;                            # Set gap from top.		      #
                                                                                      #
  my($seq_l) = length $prom_seq;						      #
  my($num_seq) = int($seq_l/$ROW);                  # Number of EVEN sequences.	      #
                                                                                      #
  if($seq_l%$ROW!=0) { $num_seq++; }                # If not evenly divisible by ROW. #
                                                                                      #
  #------------------------------------------------------------------------------------

  #------------------------------ ROText field set-up ---------------------------------
                                                                                      #
  my($gap) = $f_rows+3;                               				      #
  my($h)= $num_seq*$gap+(TOP-1);						      #
                                                                                      #
  my($tb) = $canvas->ROText(qw/-setgrid true/, -width => $ROW,			      #
                             -height => $h, -font => "small");			      #
                                                                                      #
  #------------------------------------------------------------------------------------

  my($empty)=" ";
  my($empty_line);

  my($a)=0;
  while ($a<$ROW) {
    substr($empty_line, $a, 1) = $empty;
    $a++;
  }
    
  $tb->insert('insert', "\n");

  my($cord1) = $self->{prom_obj}->{prom_start}+1;
  my($cord2) = $self->{prom_obj}->{prom_end}+1;
 
  my($strand) = $self->{prom_obj}->{strand};
  
  my($ROW_R);
  my($start, $end);

  if ($strand eq 'W') { 
    $ROW_R = $ROW; 
    $start = $cord1; $end = $cord2;
  } 
  elsif ($strand eq 'C') { 
    $ROW_R = $ROW*(-1); 
    $start = $cord2; $end = $cord1;
  }
    
  my($b)=0;
  while ($b<$num_seq) {
    my($row_start) = $start+($b*$ROW_R)+1;
    my($row_end) = $start+($b+1)*$ROW_R;

    if ($b==0) { $row_start = $start; }
    if (($b+1) == $num_seq) { $row_end = $end; }

    my($range) = "[".$row_start."-".$row_end."]";
          
    $tb->insert('insert', "$range\n");
    $tb->insert('insert', substr($prom_seq, $b*$ROW, $ROW));         # Promoter segment insertion.
    $tb->insert('insert', "\n");
                
    my($k)=0;                                                        # TF e_line insertion.
    my($tf_rows) = $f_rows+1;
   
    while ($k<$tf_rows) {
      $tb->insert('insert', "$empty_line");
      $tb->insert('insert', "\n");
      $k++;
    }
    $b++;
  }
  return $tb;
}

#----------------------------- Motif_color_TF_insert sub. -----------------------------

# 5. &motif_color_TF_insert subroutine to color motifs and insert corresponding TFs.

sub motif_color_TF_insert($) {
  my($self, $tb) = @_;
  
  my($prom_seq) = $self->{prom_obj}->{prom_seq};
  my($strand) = $self->{prom_obj}->{strand};

  my(@ordered) = @{ $self->{fact_obj}->{final_list} };    # Assign final TF list.

  my($seq_l) = length $prom_seq;

#print "\n\n";

  foreach my $v (@ordered) {
    my($m_length)= length $v->{motif};
    my($f_length) = (length $v->{factor})+1;
    my($offset) = $v->{offset};
   
    #print "--> $v->{factor}\t$v->{offset}\t$v->{new_offset}\t$v->{tf_line}\t$v->{color}\n";       # Print option.
     
    if (0<=$offset && $offset<=$seq_l) {
      my($sub1) = $v->{prom_line}.".".$v->{new_offset};
      my($sub2) = $v->{prom_line}.".".($v->{new_offset}+$m_length);
                  
      $tb->tagAdd($v->{color}, $sub1, $sub2);                      # Promoter segment display.
                      
      if (defined($v->{prom_wrap_o})) {
        my($sub3) = $v->{prom_wrap_l}.".".0;
        my($sub4) = $v->{prom_wrap_l}.".".$v->{prom_wrap_o};
                                      
        $tb->tagAdd($v->{color}, $sub3, $sub4);
      }
   
      if (defined($v->{tf_wrap})) { $v->{new_offset} = 0; }
 
      my($coord1) = $v->{tf_line}.".".$v->{new_offset}; 
      my($coord2) = $v->{tf_line}.".".($v->{new_offset}+$f_length);
                                  
      $tb->delete($coord1, $coord2);
      $tb->insert($coord1, $v->{factor}.$v->{orient}, $v->{color});
    }
    else { print "OUT OF RANGE - factor: $v->{factor}\toffset: $v->{offset}\n\n"; }
  }
  return $tb;
}                                                                      
