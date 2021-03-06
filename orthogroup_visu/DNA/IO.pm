#---------------------------------------------------------------------------------- 
# Author: Mirko Palla.
# Date: March 27, 2006.
# For: Task 3 at the Church Lab - Genetics Department, 
# Harvard Medical School.
# 
# Purpose: This program contains the complete code for 
# module IO, handling output for Fact, Mutant and Unit
# objects in Perl.
#
# Associated subroutines: 
#			  1.  &evaluate 	        2. &validate_cf
#	INPUT		  3.  &validate_ds
#
#			  4.  &TF_search_output         5. &results
#			  6.  &output	                7. &mode1_output  
#	OUTPUT		  8.  &mode2_output             9. &mode3_output  
#			  10. &remove_motif_output     11. &remove_overlap_output
#			  12. &db_output    	       13. &hg_output
#	                  14. &overlap_output	       15. &standard_map
#			  16. &gene_search_output      17. &get_matrix
#			  18. &final_map_output        19. &binding_summary_output
#			  20. &global_hg_result        21. &remove_units_output
#			  22. &random_kmer_output      23. R_KMER (constant value)
#			  24. &ortho_output_i	       25. &ortho_output_o
#			  26. &ortho_table_output      27. &list_table	    
#			  28. &excluded_proms_IO       29. &statistics_output
#			  30. &d_to_gene_output        31. &o_TF_search_output 
#			  32. &list_ortho_tfs          33. &response_IO
#			  34. &tf_pair_output          35. &tf_overlap_output
#			  36. &hg_ol_db_output	       37. &ol_tf_ratio_output
#			  38. &binding_ol_output       39. &kl_divergence_output
#---------------------------------------------------------------------------------- 

package DNA::IO;
                                                                                        
use strict;
use vars qw ($i $factor $offset $dir $score $binding);
use vars qw ($motif $chrom $prom $conserv $condition);

1;

#--------------------------------------------------------------------------------------#
#					 INPUT					       #
#--------------------------------------------------------------------------------------#

#------------------------------------ Constants ---------------------------------------
                                                                                        
# 0. constant arrays for IO operations.

use constant FIELD1 => ("Base_dir", "Prom_file", "Conv_file", "File_name");
use constant FIELD2 => ("TF_dir", "Finalized_db", "PSSMS", "Matrix");
use constant FIELD3 => ("Kmer_l", "Window_gap", "BP_diff");
use constant FIELD4 => ("Overlap_gap", "Keep_perc", "Ruin_perc");
use constant FIELD5 => ("ROW", "Width", "Height", "Font");
use constant FIELD6 => ("Filtered_orthogroups", "Promoters", "Other_dir", "Aliases");
use constant FIELD7 => ("Binding_THRESH", "Condition", "Conservation", "Probability", "Affinity_file");
use constant FIELD8 => ("Mode", "Distance", "Unit_file");

use constant TYPE1 => ("DIR", "FILE", "FILE", "STRING");
use constant TYPE2 => ("DIR", "FILE", "DIR", "BINARY");
use constant TYPE3 => ("KMER_L", "WIN", "INT");
use constant TYPE4 => ("INT", "INT", "INT");
use constant TYPE5 => ("INT", "INT", "INT", "INT");
use constant TYPE6 => ("FILE", "DIR", "DIR", "DIR");
use constant TYPE7 => ("FLOAT", "STRING", "INT", "FLOAT", "FILE");
use constant TYPE8 => ("INT_SET", "INT", "FILE");

#------------------------------------ Evaluate sub. -----------------------------------
                                                                                        
# 1. &evaluate to filter configuration input for validity.
                                                                                        
sub evaluate($$$) {
  my($config_file, $FIELD, $TYPE) = @_;

  #----------------------------- Config. file handling --------------------------------

  my(%config_input);                                       # Store config. input in hash.

  for (my $k=0; $k < @$FIELD; $k++) {
    my($field) = @$FIELD->[$k];
    my($type) = @$TYPE->[$k];

    open(IN, $config_file) || die "--> Cannot open config. file: $config_file\n";

    while (<IN>) {
      if (/\b$field\b/) { 
        chomp ($config_input{$field} = <IN>); 
       
        SWITCH: {

          if ($type eq "DIR") {     
            die "\n--> Error: $config_input{$field} - not valid input directory -> [$field]\n\n" 
            unless (-d $config_input{$field}); 
            last SWITCH; 
          }

          if ($type eq "FILE") {     
            die "\n--> Error: $config_input{$field} - not valid input file -> [$field]\n\n" 
            unless (-f $config_input{Base_dir}."/".$config_input{$field}); 
            last SWITCH; 
          }
              
          if ($type eq "INT") {     
            die "\n--> Error: $config_input{$field} - not valid integer value -> [$field]\n\n" 
            unless $config_input{$field} =~ /^\d+$|N\/A/;
            last SWITCH; 
          }

          if ($type eq "FLOAT") {     
            die "\n--> Error: $config_input{$field} - not valid floating point number -> [$field]\n\n" 
            unless $config_input{$field} =~ /^[+-]?\ *(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?$|N\/A/;
            last SWITCH; 
          }

          if ($type eq "INT_LIST") {
            my(@int_list) = split(/\s/, $config_input{$field});
            
            foreach my $int (@int_list) {
              die "\n--> Error: $int - not valid integer value -> [$field]\n\n" 
              unless $int =~ /^\d+$/;
            }
            last SWITCH; 
          }

          if ($type eq "INT_SET") {     
            die "\n--> Error: $config_input{$field} - not valid integer values {1,2,3} -> [$field]\n\n" 
            unless ($config_input{$field} == 1 || $config_input{$field} == 2 ||
                                                    $config_input{$field} == 3);
            last SWITCH; 
          }

          if ($type eq "BINARY") { validate_ds($config_input{$field}); }
          
          if ($type eq "KMER_L") {     
            die "\n--> Error: $config_input{$field} - not valid kmer length -> [$field]\n\n" 
            unless ($config_input{$field} =~ /^\d+$/ && $config_input{$field} != 0 && 
                                                       	    $config_input{$field} != 1);
            last SWITCH; 
          }

          if ($type eq "WIN") {     
            die "\n--> Error: $config_input{$field} - not valid sliding window size -> [$field]\n\n" 
            unless ($config_input{$field} =~ /^\d+$/ && $config_input{$field} != 0); 
            last SWITCH; 
          }
              
          if ($type eq "STRING") {     
            die "\n--> Error: $config_input{$field} - not valid string value -> [$field]\n\n" 
            unless ($config_input{$field} =~ /.*|N\/A/);
            last SWITCH; 
          }
        }
      }  
    } 
  }
  return %config_input;	     # If all configuration input were valid: return arguments.
}

#---------------------------------- Validate_ds sub. ----------------------------------
                                                                                        
# 2. &validate_ds to filter user input for {0,1} validity.
                                                                                        
sub validate_ds($) {
  my($arg) = shift;

  die "\n--> Error: $arg - not valid input value {0,1}\n\n" 
  unless ($arg eq '0' || $arg eq '1');
  return $arg,
}  

#---------------------------------- Validate_cf sub. ----------------------------------
                                                                                        
# 3. &validate_cf to filter user input for string validity.
                                                                                        
sub validate_cf($) {
  my($arg) = shift;

  die "\n--> Error: $arg - not valid input file\n\n" 
  unless (-f $arg); 
  return $arg,
}  

#--------------------------------------------------------------------------------------#
#					OUTPUT					       #
#--------------------------------------------------------------------------------------#

#-------------------------------- TF_search_output sub. -------------------------------------

# 4. &TF_search_output subroutine to output list of TFs bound to given promoter.

sub TF_search_output() {
  my($self) = shift;

  my($f_name) = $self->{parameters}{File_name};  
  my($out_dir) = $self->{parameters}{Out_dir}; 

  my(@fac) = @{ $self->{final_list} };
  my($prom) = $self->{prom_obj}->{promoter};
  my($chro) = $self->{prom_obj}->{chromosome};
  my($gene) = $self->{prom_obj}->{gene};

  my($random) = int( rand(9E5) ) + 1E5;			# Generate random number between 1E5 and 999999.
  $self->{bar_code} = $random;				# Store it in TF object.

  my($h) = "_human_".$f_name."_$random.txt";
  my($m) = "_machine_".$f_name."_$random.txt";
  my($h_name) = $gene.$h;	                	 # Create the name of the file.
  my($m_name) = $gene.$m;	                	 # Create the name of the file.

  if (!(-d $out_dir)) {		                         # Check if directory exists. 
    mkdir ($out_dir, 0777) || die "--> Error: unable to create output directory $out_dir\n\n";
    chdir $out_dir;
  }
  else { chdir $out_dir; }

  open H_RESULT, ">$h_name";                             # Opens a h_file for output.
  print H_RESULT "--> { GENE: $gene <=> PROMOTER: $prom } -> CHROMOSOME: $chro\n\n";
  print H_RESULT "--> Data set in use: $self->{parameters}{Finalized_db}\n\n";

  open M_RESULT, ">$m_name";                             # Opens a m_file for output.

  my($e)=0;
  while ($e < @fac) {
    $i = ">".$e; 
    $factor = $fac[$e]->{factor}.$fac[$e]->{orient}; 
    $offset = $fac[$e]->{offset};
    $motif = $fac[$e]->{motif}; 

    $binding = $fac[$e]->{binding}; 
    $condition = $fac[$e]->{condition}; 
    $conserv = $fac[$e]->{conservation}; 
    
    my($scr) = $fac[$e]->{score};
    #$score = sprintf("%.10f", $scr);		         # Display 10 digits after decimal. 
    $score = $scr;		  		         # Display scientific notation. 
   
    write H_RESULT;
    print M_RESULT ">$e\t$factor\t$offset\t$score\t$binding\t$condition\t$conserv\t$motif\n";
    $e++;
  }
  close H_RESULT;                                 	# Close the h_file.
  close M_RESULT;                                 	# Close the m_file.
  print "\n--> $h_name - written."; 
  print "\n--> $m_name - written.";
}

format H_RESULT =
@<<<   @<<<<<<< @<<<<  @<<<<<<<<<<<<<<<<<<<<   @<<<<<<<<<<<<<   @<<<<<<<   @<   @<<<<<<<<<<<<<<<<   
$i, $factor, $offset, $score, $binding, $condition, $conserv, $motif
.

#----------------------------------- Results sub. -------------------------------------

# 5. &results subroutine to output required results: Promoter.pm.

sub results($@) {
  my($self, $mode, @prom_list) = @_;

  my($out_dir) = $self->{parameters}{Out_dir};
  my($f_name) = $self->{parameters}{File_name};

  my($query_gene) = $self->{gene_q};			# Retrieve gene.
 
  my($random) = int( rand(9E5) ) + 1E5;			# Generate random number between 1E5 and 999999.
  $self->{bar_code} = $random;				# Store it in promorter object.

  my($h) = "_".$f_name."_".$random.".txt";
  my($h_name) = "promoter_by_".$query_gene.$h;         	# Create the name of the file.

  if (!(-d $out_dir)) {		                        # Check if directory exists. 
    mkdir ($out_dir, 0777) || die "--> Error: unable to create output directory $out_dir\n\n";
    chdir $out_dir;
  }
  else { chdir $out_dir; }

  open RESULT, ">$h_name";                              # Opens a h_file for output.
  print RESULT "--> GENE: $query_gene\n";
  print RESULT "--> Corresponding promoter list:\n\n";

  my($r)=0;
  while ($r < @prom_list) {			        # Access 'gene' / 'promoter' field in promoter object.
    my($prom) = $prom_list[$r];
    print RESULT ">$r: $prom->{chromosome}. $prom->{gene} -- ";
    print RESULT "[$prom->{promoter}]\n";
    $r++;
  }
  print RESULT "________________________________\n";

  if ($mode==0) {
    print RESULT "\n--> Promoters not in conversion file:\n\n";

    my($n)=0;
    while ($n < @{$self->{no_eq}}) {                    # Access 'no_eq' field in promoter object.
      print RESULT ">$n: $self->{no_eq}[$n]->{promoter}\n";
      $n++;
    }
  }
  
  if ($mode==1) {
    my($factor) = $self->{no_eq}[0]->{factor};
    print RESULT "\n--> TF $factor out of range of annotated ORFs:\n\n";

    my($n)=0;
    while ($n < @{$self->{no_eq}}) {                    # Access 'no_eq' field in promoter object.
      my($f) = $self->{no_eq}[$n];
      print RESULT ">$n: $f->{chromosome}. -> [$f->{motif_start},$f->{motif_end}]\n";
      $n++;
    }
  }  
  close RESULT;
  print "\n--> $h_name - written."; 
}

#------------------------------------ Output sub. -------------------------------------

# 6. &output subroutine to record required results: Mutagen_set.pm.

sub output($) {
  my($self, $mod) = @_;
  
  my($f_name) = $self->{parameters}{File_name};  
  my($out_dir) = $self->{parameters}{Out_dir}; 

  my($gene) = $self->{prom_obj}->{gene};		# Get gene name.
  my($prom) = $self->{prom_obj}->{promoter};		# Get promoter name.
  my($chro) = $self->{prom_obj}->{chromosome};          # Get chromosome no.

  my($random) = int( rand(9E5) ) + 1E5;			# Generate random number between 1E5 and 999999.
  $self->{bar_code} = $random;				# Store it in Mutagen object.
  
  if ($mod ne "permuted" && $mod ne "randomized" && $mod ne "scanned") {
    die "\n--> Error: $mod - wrong third argument\n--> Usage of parameter: permuted|randomized!\n\n"; }

  my($name) = $gene."_".$mod."_".$f_name.		# Mutant promoter sequence output file.
                             "_".$random.".txt";   

  my($d_name) = $gene."_".$mod."_".$f_name.		# Mutation documentation file.
                              "_".$random.".doc";
                    
  if (!(-d $out_dir)) {		                        # Check if directory exists. 
    mkdir ($out_dir, 0777) || die "--> Error: unable to create output directory $out_dir\n\n";
    chdir $out_dir;
  }
  else { chdir $out_dir; }
                                                                    
  open OUTPUT, ">$name";			        # Open output file.
  open DOC, ">$d_name";				        # Open documentation file.

  print OUTPUT "--> { GENE: $gene <=> PROMOTER: $prom } -> CHROMOSOME: $chro\n\n";
  print DOC "--> { GENE: $gene <=> PROMOTER: $prom } -> CHROMOSOME: $chro\n\n";

  if ($mod eq "randomized") {
    while ( (my $k, my $v) = each %{$self->{r_kmer}} ) {
      print DOC "> Random kmer $k with log-value $v rejected.\n";
    }
    print DOC "\n> Random kmer $self->{kmer}[0] with average log-value $self->{kmer}[1] for randomization.\n\n";
  }
  
  my($r)=0;
  foreach my $m_prom (@{$self->{$mod}}) {    	         # Access 'mod' field in Mutant object.
    my($kmer) = $m_prom->{o_kmer};                   	 # Assign original kmer.
    my($kmer_l) = $m_prom->{kmer_l};          		 # Assign kmer length.

    print OUTPUT ">$r: $m_prom->{mutant}\n";    

    if ($mod ne "scanned") { 
      my($s_kmer) = $m_prom->{s_kmer};         			 # Assign substitition kmer.
      my($offset) = $m_prom->{s_offset};       			 # Assign substitution offset.
      print DOC ">$r: Replace $kmer ($offset) with $mod kmer $s_kmer ($kmer_l).\n";
    }
    else {
      print DOC ">$r: Replace occurrences of $kmer at: ";
   
      foreach my $key (sort { $m_prom->{s_kmer}->{$a} <=> 
        $m_prom->{s_kmer}->{$b} } keys %{ $m_prom->{s_kmer} }) {   # Retrieve kmer from promoter.
        print DOC "($m_prom->{s_kmer}->{$key}) with $key ($kmer_l) |";
      }
      print DOC "\n";
    }
    $r++;
  }
  close OUTPUT;
  close DOC;
  
  print "\n--> $name - written."; 
  print "\n--> $d_name - written."; 
}

#---------------------------------- Mode1_output sub. ---------------------------------

# 7. &mode1_output subroutine to record required results: Unit_set.pm.

sub mode1_output() {
  my($self) = shift;

  my($f_name) = $self->{parameters}{File_name};  
  my($out_dir) = $self->{parameters}{Out_dir};  
  my($mode) = $self->{parameters}{Mode};

  my($gene) = $self->{prom_obj}->{gene};                   # Get gene name.
  my($prom) = $self->{prom_obj}->{promoter};               # Get promoter name.
  my($chro) = $self->{prom_obj}->{chromosome};             # Get chromosome no.
  
  my(@units) = @{ $self->{units} };
  my($unit_num) = scalar @units;
  
  my($random) = int( rand(9E5) ) + 1E5;			   # Generate random number between 1E5 and 999999.
  $self->{bar_code} = $random;				   # Store it in Unit object.
  
  my($name) = $gene."_units".$mode."_".$f_name.		   # Create name.
                "_".$random.".txt";
	
  if (!(-d $out_dir)) {		                           # Check if directory exists. 
    mkdir ($out_dir, 0777) || die "--> Error: unable to create output directory $out_dir\n\n";
    chdir $out_dir;
  }
  else { chdir $out_dir; }

  open OUTPUT, ">$name";					  
  print OUTPUT "--> { GENE: $gene <=> PROMOTER: $prom } -> CHROMOSOME: $chro\n\n";
  print OUTPUT "--> Data set in use: $self->{parameters}{Finalized_db}\n";
  print OUTPUT "--> Motif unit mode: $mode\n";
  print OUTPUT "--> Number of total units: $unit_num\n\n";

  foreach my $unit (@units) {				  # Iterate over unit elements.
    my($unit_n) = $unit->{unit_no};
    my($factor) = $unit->{factors}[0]->{factor};
    my($offset) = $unit->{factors}[0]->{offset};

    print OUTPUT "--> UNIT [$unit_n]\t$factor\t$offset\n";
  }
  close OUTPUT;
  print "\n--> $name - written."; 
}

#----------------------------------- Mode2_output sub. --------------------------------

# 8. &mode2_output subroutine to record required results: Unit_set.pm.

sub mode2_output() {
  my($self) = shift;

  my($f_name) = $self->{parameters}{File_name};  
  my($out_dir) = $self->{parameters}{Out_dir};  
  my($mode) = $self->{parameters}{Mode};

  my($gene) = $self->{prom_obj}->{gene};                   # Get gene name.
  my($prom) = $self->{prom_obj}->{promoter};               # Get promoter name.
  my($chro) = $self->{prom_obj}->{chromosome};             # Get chromosome no.

  my($units) = $self->{units};
  my($unit_num) = scalar (keys %{ $units });

  my($random) = int( rand(9E5) ) + 1E5;			   # Generate random number between 1E5 and 999999.
  $self->{bar_code} = $random;				   # Store it in Unit object.

  my($name) = $gene."_units".$mode."_".$f_name.
                "_".$random.".txt";

  if (!(-d $out_dir)) {		                           # Check if directory exists. 
    mkdir ($out_dir, 0777) || die "--> Error: unable to create output directory $out_dir\n\n";
    chdir $out_dir;
  }
  else { chdir $out_dir; }

  open OUTPUT, ">$name";					  
  print OUTPUT "--> { GENE: $gene <=> PROMOTER: $prom } -> CHROMOSOME: $chro\n\n";
  print OUTPUT "--> Data set in use: $self->{parameters}{Finalized_db}\n";
  print OUTPUT "--> Motif unit mode: $mode\n";
  print OUTPUT "--> Number of total units: $unit_num\n";

  foreach my $key (sort { $units->{$a}->{unit_no} <=> 
                          $units->{$b}->{unit_no} } keys %$units) {
                          
    my($val) = scalar @{ $units->{$key}->{factors} };
    my($unit_n) = $units->{$key}->{unit_no};

    print OUTPUT "\n--> UNIT [$unit_n]: $key - number of elements: $val\n";

    my($l)=0;
    while ($l < @{ $units->{$key}->{factors} }) {
      print OUTPUT "--> Element $l:\t$units->{$key}->{factors}[$l]->{offset}\n";
      $l++;
    }
    print OUTPUT "___________________________________________\n";
  }
  close OUTPUT;
  print "\n--> $name - written."; 
}

#----------------------------------- Mode3_output sub. --------------------------------

# 9. &mode3_output subroutine to record required results: Unit_set.pm.

sub mode3_output() {
  my($self) = shift;

  my($f_name) = $self->{parameters}{File_name};  
  my($out_dir) = $self->{parameters}{Out_dir};  
  my($mode) = $self->{parameters}{Mode};
  
  my($gene) = $self->{prom_obj}->{gene};                   # Get gene name.
  my($prom) = $self->{prom_obj}->{promoter};               # Get promoter name.
  my($chro) = $self->{prom_obj}->{chromosome};             # Get chromosome no.
  
  my($units) = $self->{units};
  
  my($random) = int( rand(9E5) ) + 1E5;			   # Generate random number between 1E5 and 999999.
  $self->{bar_code} = $random;				   # Store it in Unit object.
  
  my($name) = $gene."_units".$mode."_".$f_name.
                "_".$random.".txt";

  if (!(-d $out_dir)) {		                           # Check if directory exists. 
    mkdir ($out_dir, 0777) || die "--> Error: unable to create output directory $out_dir\n\n";
    chdir $out_dir;
  }
  else { chdir $out_dir; }

  open OUTPUT, ">$name";					  
  print OUTPUT "--> { GENE: $gene <=> PROMOTER: $prom } -> CHROMOSOME: $chro\n\n";
  print OUTPUT "--> Data set in use: $self->{parameters}{Finalized_db}\n";
  print OUTPUT "--> Motif unit mode: $mode\n";  
  print OUTPUT "--> Number of total units: $self->{unit_num}\n";
  print OUTPUT "--> Distance parameter: $self->{parameters}{Distance}\n";
  
  foreach my $key (sort { $units->{$a}[0]->{unit_no} <=> 
                          $units->{$b}[0]->{unit_no} } keys %$units) {
                          
    my($val) = scalar @{ $units->{$key} };
    print OUTPUT "\n--> $key - number of units: $val\n";

    my($k)=0;
    while ($k < @{ $units->{$key} }) {
      my($unit_n) = $units->{$key}[$k]->{unit_no};
      print OUTPUT "--> UNIT [$unit_n]:\n";

      my($l)=0;
      while ($l < @{ $units->{$key}[$k]->{factors} }) {
        print OUTPUT "--> Element $l: $units->{$key}[$k]->{factors}[$l]->{offset}\n";
        $l++;
      }
      $k++;
    }
    print OUTPUT "________________________________\n";
  }
  close OUTPUT;
  print "\n--> $name - written."; 
}

#------------------------------- Remove_motif_output sub. -----------------------------

# 10. &remove_motif_output subroutine to record promoter sequences where all copies
# of TF motif is removed.

sub remove_motif_output() {
  my($self) = shift;

  if ( defined @{ $self->{removed} } ) {
  
    my($f_name) = $self->{parameters}{File_name};  
    my($out_dir) = $self->{parameters}{Out_dir}; 
  
    my($gene) = $self->{prom_obj}->{gene};                 # Get gene name.
    my($prom) = $self->{prom_obj}->{promoter};             # Get promoter name.
    my($chro) = $self->{prom_obj}->{chromosome};           # Get chromosome no.

    my($random) = int( rand(9E5) ) + 1E5;		   # Generate random number between 1E5 and 999999.
    $self->{bar_code} = $random;			   # Store it in Unit object.

    my($name) = $gene."_removed_".$f_name.		   # Promoter mutation with removed units.
                          "_".$random.".txt";

    my($d_name) = $gene."_removed_".$f_name.		   # Mutation documentation file.
                          "_".$random.".doc";
              
    if (!(-d $out_dir)) {		                   # Check if directory exists. 
      mkdir ($out_dir, 0777) || die "--> Error: unable to create output directory $out_dir\n\n";
      chdir $out_dir;
    }
    else { chdir $out_dir; }

    open OUTPUT, ">$name";					  
    open DOC, ">$d_name";
			  
    print OUTPUT "--> { GENE: $gene <=> PROMOTER: $prom } -> CHROMOSOME: $chro\n\n";
    print DOC "--> { GENE: $gene <=> PROMOTER: $prom } -> CHROMOSOME: $chro\n\n";
    print DOC "--> Data set in use: $self->{parameters}{Finalized_db}]\n";

    my($mode);
    if ($self->{mode}==0) { $mode = "permuted"; }
    else { $mode = "randomized"; }

    my($r)=0; 
    foreach my $m_prom (@{ $self->{removed} }) {           # Access 'removed' field in Mutant object.
      print OUTPUT ">$r: $m_prom->{mutant}\n";

      my($unit) = $m_prom->{unit};	            	   # Assign unit name.
      my($number) = $m_prom->{no_elements};          	   # Assign number of units elements.

      print DOC "\n--> [$r] In UNIT $unit - number of elements: $number\n";

      for (my $k=0; $k < @{ $m_prom->{o_kmer} }; $k++) {     

        my($kmer) = $m_prom->{o_kmer}->[$k];            	   # Assign original kmer.
        my($kmer_l) = $m_prom->{kmer_l}->[$k];          	   # Assign kmer length.
        my($s_kmer) = $m_prom->{s_kmer}->[$k];          	   # Assign substitition kmer.
        my($offset) = $m_prom->{s_offset}->[$k];        	   # Assign substitution offset.

        print DOC ">$k: Replace $kmer ($offset) with $mode kmer $s_kmer ($kmer_l).\n";
      }
      $r++;
      print DOC "____________________________________________________________________________\n";
    }
    close OUTPUT;
    close DOC;

    print "\n--> $name - written."; 
    print "\n--> $d_name - written."; 
  }
  else { 
    print "\n--> No UNIT object present, such that, the number of TFs in a unit is greater than 1."; 
    print "\n--> Process exited - remove_motif.pl\n\n"; 
    exit;
  }
}

#----------------------------- Remove_overlap_output sub. -----------------------------

# 11. &remove_overlap_output subroutine to record promoter sequences where all motif 
# overlaps are removed in a vice and versa manner.

sub remove_overlap_output() {
  my($self) = shift;

  if ( defined @{ $self->{removed} } ) {

    my($f_name) = $self->{parameters}{File_name};  
    my($out_dir) = $self->{parameters}{Out_dir}; 

    my($gene) = $self->{prom_obj}->{gene};                 # Get gene name.
    my($prom) = $self->{prom_obj}->{promoter};             # Get promoter name.
    my($chro) = $self->{prom_obj}->{chromosome};           # Get chromosome no.

    my($random) = int( rand(9E5) ) + 1E5;		   # Generate random number between 1E5 and 999999.
    $self->{bar_code} = $random;			   # Store it in Unit object.

    my($name) = $gene."_OL-remove_".$f_name.		   # Promoter mutation with removed units.
                          "_".$random.".txt";

    my($d_name) = $gene."_OL-remove_".$f_name.		   # Mutation documentation file.
                          "_".$random.".doc";
              
    if (!(-d $out_dir)) {		                   # Check if directory exists. 
      mkdir ($out_dir, 0777) || die "--> Error: unable to create output directory $out_dir\n\n";
      chdir $out_dir;
    }
    else { chdir $out_dir; }

    open OUTPUT, ">$name";					  
    open DOC, ">$d_name";
				  
    print OUTPUT "--> { GENE: $gene <=> PROMOTER: $prom } -> CHROMOSOME: $chro\n\n";
    print DOC "--> { GENE: $gene <=> PROMOTER: $prom } -> CHROMOSOME: $chro\n\n";
    print DOC "--> Data set in use: $self->{parameters}{Finalized_db}\n";
    print DOC "--> Overlap gap parameter: $self->{parameters}{Overlap_gap}\n\n";
  
    my($r)=0; 
    foreach my $m_prom (@{ $self->{removed} }) {           # Access 'removed' field in Mutant object.
      print OUTPUT ">$r: $m_prom->{mutant}\n";

      my($m_A) = $m_prom->{overlap_pair}->{motif_A};       # Untouched motif.
      my($m_B) = $m_prom->{overlap_pair}->{motif_B}; 	   # Mutated motif.

      my($OL)  = $m_prom->{overlap_pair}->{overlap}; 	   # Overlap value.
      my($sign) = $m_prom->{overlap_pair}->{sign}; 	   # Overlap type: +/-.
      my($type) = $m_prom->{overlap_pair}->{type}; 	   # Overlap type: {"fake","real->I","real->II"}.
      
      my($m_A_motif) = $m_A->{motif}; 
      my($m_A_offset) = $m_A->{offset}; 
      my($m_A_orient) = $m_A->{orient}; 
      my($m_A_score) = $m_A->{score}; 

      my($m_B_motif) = $m_B->{motif}; 
      my($m_B_orient) = $m_B->{orient};
      my($m_B_offset) = $m_B->{offset}; 
      my($m_B_score) = $m_B->{score}; 

      my($A) = $m_A->{factor}.$m_A_orient;		   # [motif_name.orient] - [STE12+]
      my($B) = $m_B->{factor}.$m_B_orient; 

      print DOC "\n--> [$r] OVERLAP: $OL -> {$B <-> $A} => TYPE: $type\n";

      #------------------------ MUTATED MOTIF PARAMETERS -----------------------------

      my($kmer) = $m_prom->{o_kmer};            	   # Assign original kmer (mutate).
      my($s_kmer) = $m_prom->{s_kmer};          	   # Assign substitition kmer.
      my($s_score) = $m_prom->{s_score};        	   # Assign substitution p-value.

      my($x_kmer) = $m_prom->{x_kmer};          	   # Assign completely mutated kmer.
      my($kmer_l) = length $m_B_motif;          	   # Assign kmer length.

      #----------------------- DOMINANT MOTIF PARAMETERS -----------------------------

      my($d_kmer) = $m_prom->{d_kmer};            	   # Assign original kmer (keep).
      my($D_kmer) = $m_prom->{D_kmer};          	   # Assign substitition kmer.
      my($D_score) = $m_prom->{D_score};        	   # Assign substitution p-value.

      my($d_kmer_l) = length $m_A_motif;          	   # Assign kmer length.

      print DOC ">$r: MUTATED [$B] - Replace $kmer (offset: $m_B_offset, score: $m_B_score) with PSSMS knocked-out kmer $s_kmer (score: $s_score) -> total mutant: {$x_kmer} - (kmer length: $kmer_l).\n";

      if ($d_kmer eq $D_kmer) { print DOC ">   UNTOUCH [$A] - Motif $d_kmer (offset: $m_A_offset, score: $m_A_score) - (kmer length: $d_kmer_l).\n"; }
      else { print DOC ">   UNTOUCH [$A] - Replace $d_kmer (offset: $m_A_offset, score: $m_A_score) with PSSMS knocked-out kmer $D_kmer (score: $D_score) - (kmer length: $d_kmer_l).\n"; }

      if ($sign eq '+' && defined ($m_prom->{i_kmer})) {
        for (my $i=0; $i < @{ $m_prom->{i_kmer} }; $i++) {     
         
          my($f_kmer) = $m_prom->{f_kmer}[$i];         	   # Assign fragment kmer.
          my($i_kmer) = $m_prom->{i_kmer}[$i];         	   # Assign intermediate kmer.
          my($i_offset) = $m_prom->{i_offset}[$i];     	   # Assign intermediate kmer offset.
          my($ikmer_l) = $m_prom->{ikmer_l}[$i];      	   # Assign intermediate kmer length.
         
          print DOC ">$i: I. STEP [$B] - Replace $f_kmer (fregment offset: $i_offset) with PSSMS knocked-out kmer $i_kmer (fregment length: $ikmer_l).\n" unless ($ikmer_l==0);
        }
      }  

      if ( defined @{ $m_prom->{untouched_tf} } ) {
        my($count) = 1;
        foreach my $u_tf (@{ $m_prom->{untouched_tf} }) {        # Access 'untouched_tf' field in Mutant object.
         
          my($fac) = $u_tf->{factor};
          my($mot) = $u_tf->{motif};
          my($off) = $u_tf->{offset};
          my($ori) = $u_tf->{orient};
            
          print DOC ">$count: UNTOUCH [$fac$ori] - Motif $mot (offset: $off).\n";
          $count++;
        }
      }

      #----------------------- FINAL STATISTICS PARAMETERS -----------------------------

      if($m_A_orient eq '-') { $D_kmer = &DNA::Func::rev_complement($D_kmer); }
      if($m_B_orient eq '-') { $s_kmer = &DNA::Func::rev_complement($s_kmer); }

      print DOC "\n--------------------------- FINAL STATISTICS ----------------------------------\n\n";
      print DOC "--> MUTATED MOTIF:  [$B]\n";
      print DOC "--> ORIGINAL: $m_B_motif -> P-value: $m_B_score\n";
      print DOC "--> MUTATED : $s_kmer -> P-value: $s_score\n\n";
      print DOC "--> DOMINANT MOTIF: [$A]\n";
      print DOC "--> ORIGINAL: $m_A_motif -> P-value: $m_A_score\n";
      print DOC "--> MUTATED : $D_kmer -> P-value: $D_score\n";

      $r++;
      print DOC "____________________________________________________________________________\n";
    }
    close OUTPUT;
    close DOC;

    print "\n--> $name - written."; 
    print "\n--> $d_name - written."; 
  }
  else { 
    print "\n--> No OVERLAP object present in given data set: $self->{parameters}{Set_selection}";
    print "\n--> Process exited - remove_motif.pl\n\n"; 
    exit;
  }
}

#------------------------------- Motif_db_output sub. ----------------------------------

# 12. &db_output subroutine to record motif database content.

sub db_output($$$) {
  my($self, $TF_hash, $file, $m) = @_;

  my($out_dir) = $self->{parameters}{Out_dir};

  my($name);
  if (defined $file) {
    if (!defined $m) { $name = (split(/\./, $file))[0].".m_db"; }   # Promoter set by motif. 
    elsif ($m eq 'p') { $name = (split(/\./, $file))[0].".p_db"; }  # Motif set by promoter.
    else { print "--> No such output option: $m - usage 'p'.\n"; }
  }
  else { $name = "nir.m_db"; }
  
  if (!(-d $out_dir)) {		   	                            # Check if directory exists. 
    mkdir ($out_dir, 0777) || die "--> Error: unable to create output directory $out_dir\n\n";
    chdir $out_dir;
  }
  else { chdir $out_dir; }

  open OUTPUT, ">$name";
  my($header) = "--> [Transcription factor] database - ";
  
  my($tail);
  my($type);

  if ($m eq 'p') { 
    $tail = "TF list by promoter:\n"; 
    $type = "PR"; 
  }
  else { 
    $tail = "promoter list by TF:\n"; 
    $type = "TF"; 
  }
  print OUTPUT $header.$tail;

  foreach my $key (sort keys %$TF_hash) {

    if ($key ne " ") {
      my($val) = scalar @{ $TF_hash->{$key} };
      print OUTPUT "\n--> $type: $key -> number of targets: $val\n";

      if ($m eq 'p') { 
        @{ $TF_hash->{$key} } = sort { $a->{offset} <=> $b->{offset} } @{ $TF_hash->{$key} };
        @{ $TF_hash->{$key} } = sort { $a->{factor} cmp $b->{factor} } @{ $TF_hash->{$key} };
      }     
       
      my($l)=0;
      while ($l < $val) {
        my($tf) = $TF_hash->{$key}[$l];

        my($facto) = $tf->{factor};
        my($offse) = $tf->{offset};
        my($orien) = $tf->{orient};
        my($motif) = $tf->{motif};
        my($score) = $tf->{score};

        my($chrom) = $tf->{chromosome};
        my($promo) = $tf->{promoter};

        print OUTPUT ">$l\t$facto$orien\t$offse\t$score\t$motif\t$chrom\t$promo\n"; 
        $l++;
      }
      print OUTPUT "________________________________________________________________________\n";
    }
  } 
  close OUTPUT;
  print "\n--> $name - written.";
}

#---------------------------------- Hg_result sub. ------------------------------------

# 13. &hg_result subroutine to output hyper-geometric comparison results.

sub hg_result() {
  my($self) = shift;
  
  my($no_proms) = $self->{no_proms};				   # Number of promoters in yeast.
  my($num_set1) = $self->{num_set1};				   # Number of promoters bound by gene1.
  my($num_set2) = $self->{num_set2};				   # Number of promoters bound by gene2.
  my($num_isect) = $self->{num_isect}; 				   # Number of promoters in intersection.
  
  my($result) = $self->{result};				   # Hyper-geometric value in direction 1. 
  
  my(@genes) = @{ $self->{genes} };				   # Gene pair of query.
  my(@isect) = @{ $self->{isect} };			           # List of intersecting promoters.
  
  my($random) = int( rand(9E5) ) + 1E5;				   # Generate random number between 1E5 and 999999.
  $self->{bar_code} = $random;					   # Store it in Stem object.

  my($gene1) = $genes[0];					   # Query gene 1.
  my($gene2) = $genes[1];					   # Query gene 2.

  my($name) = "$gene1-$gene2"."_hg_$random.txt";

  my($out_dir) = $self->{parameters}{Out_dir};

  if (!(-d $out_dir)) {		                                   # Check if directory exists. 
    mkdir ($out_dir, 0777) || die "--> Error: unable to create output directory: $out_dir\n\n";
    chdir $out_dir;
  }
  else { chdir $out_dir; }

  open DOC, ">$name";                                              # Open output file.                                                                                        
  print DOC "--> Hyper-geometric distribution results\n";
  print DOC "--> Data set in use: $self->{parameters}{Finalized_db}\n\n";

  print DOC "1. Number of promoters in yeast: $no_proms\n";
  print DOC "2. Number of $gene1 targets: $num_set1\n";
  print DOC "3. Number of $gene2 targets: $num_set2\n";
  print DOC "4. Number of targets in intersection of $gene1 & $gene2: $num_isect\n\n";
  print DOC "--> $gene1 => $gene2: $result\n";
  close DOC;                                                       # Close output file.

  print "\n--> $name - written."; 
}

#--------------------------------- Overlap_output sub. --------------------------------

# 14. &overlap_output subroutine to record overlap occurences on promoter.

sub overlap_output() {
  my($self) = shift;
  
  my($f_name) = $self->{parameters}{File_name};  
  my($out_dir) = $self->{parameters}{Out_dir};  

  my($gene) = $self->{prom_obj}->{gene};                   # Get gene name.
  my($prom) = $self->{prom_obj}->{promoter};               # Get promoter name.
  my($chro) = $self->{prom_obj}->{chromosome};             # Get chromosome no.
  
  my(@overlaps) = @{ $self->{overlaps} };
  
  my($random) = int( rand(9E5) ) + 1E5;			   # Generate random number between 1E5 and 999999.
  $self->{bar_code} = $random;				   # Store it in Overlap object.
  
  my($name) = $gene."_overlaps_".$f_name.		   # Create name.
                "_".$random.".txt";
	
  if (!(-d $out_dir)) {		                           # Check if directory exists. 
    mkdir ($out_dir, 0777) || die "--> Error: unable to create output directory $out_dir\n\n";
    chdir $out_dir;
  }
  else { chdir $out_dir; }

  open OUTPUT, ">$name";					  
  print OUTPUT "--> { GENE: $gene <=> PROMOTER: $prom } -> CHROMOSOME: $chro\n\n";
  print OUTPUT "--> Data set in use: $self->{parameters}{Finalized_db}\n\n";

  my($k)=0;
  foreach my $overlap (@overlaps) {			  # Iterate over overlap elements.
    my($factor_A) = $overlap->{motif_A}->{factor};
    my($factor_B) = $overlap->{motif_B}->{factor};

    my($offset_A) = $overlap->{motif_A}->{offset};
    my($offset_B) = $overlap->{motif_B}->{offset};

    print OUTPUT "--> OVERLAP $k.\n--> MOTIF PAIR -> [$factor_A : $factor_B] - OFFSET -> [$offset_A : $offset_B]\n\n";
    $k++;
  }
  close OUTPUT;
  print "\n--> $name - written."; 
}

#---------------------------------- Standard_map sub. -----------------------------------

# 15. &standard_map subroutine to construct standard TF database from Ernest's maps.

sub standard_map($) {
  my($self, $prom_list) = @_;

  my(%prom_list) = %{ $prom_list };  			     # Promoter library.	
  my($base_dir) = $self->{parameters}{Base_dir};             # Base directory.	

  chdir $base_dir."ernest/new_maps/";		             # Change to raw map direcoty.
  my(@files) = <*.gff>;          		    	     # Collect all .gff files in a list.

  my($out_dir) = $base_dir."ernest/standard_maps/";

  foreach my $w (@files) {
    chdir $base_dir."ernest/new_maps/";  	  	       # Change to raw map direcoty.

    print "\n--> Processing file: $w";

    open(INFO, $w) || die "--> Cannot open *.gff file: $w\n";  # Open the file for input.

    my(@TF_hits);					     # TF hit (data) array.
    while (<INFO>) {		       			
      my(@m) = split(/\s/, $_);				     # Collect TF information.

      my($chrom) = $m[0];				     # Chromosomal location. 				     
      my($m_start) = $m[3];				     # Global motif start. 				     
      my($m_end) = $m[4];				     # Global motif end.
     
      my(@prom_set) = @{ $prom_list->{$chrom} };

      my($prom_start);
      my($prom_end);

      my($proms);					     # Promoter name list, a string.
      my($marker);						
      my($l)=0;		
      
      PROM: while ($l < @prom_set) {

        $prom_start = $prom_set[$l]->{prom_start};
        $prom_end = $prom_set[$l]->{prom_end};
      
        if ($prom_start <= $m_start && $m_end <= $prom_end) {
        
          my($prom) = $prom_set[$l]->{promoter};	     # Get promoter name TF bound to.
          $proms = $proms."$prom|";			     # Collect all suck promoter names.
          $marker=1;					     # Mark existence of such promoter(s).
        }
        $l++;
      } 
      substr($proms,-1) = ';';				     # Replace last '|' with ';'

      if ($marker!=1) { $proms = 'DNE;'; }
 
      if (/[+]/) { s/[+]/-/; }				     # Change orientation mistake.
      elsif (/[-]/) { s/[-]/+/; }

      chomp; chop;					     # Remove unnecessary characters.
      my($hit) = $_."\t$proms\n";
      push @TF_hits, $hit;				     # Store TF hit linearly in array.
    }

    if (!(-d $out_dir)) {	                             # Check if directory exists. 
      mkdir ($out_dir, 0777) || die "--> Error: unable to create output directory $out_dir\n\n";
      chdir $out_dir;
    }
    else { chdir $out_dir; }

    open RESULT, ">$w";		                             # Opens an output file.
    foreach my $t (@TF_hits) { print RESULT $t; }	     # Print standarized TF hit data.
    close RESULT;
  }  
}

#------------------------------ gene_search_output sub. -------------------------------------

# 16. &gene_search_output subroutine to output list of promoters bound by given query gene.

sub gene_search_output() {
  my($self) = @_;

  my($q_gene) = $self->{gene_q};			# Get query gene.	
  my($TF_hash) = $self->{final_list};			# Get final promoter library, a hash_ref.
	
  my($f_name) = $self->{parameters}{File_name};  
  my($out_dir) = $self->{parameters}{Out_dir}; 

  my($random) = int( rand(9E5) ) + 1E5;			# Generate random number between 1E5 and 999999.
  $self->{bar_code} = $random;				# Store it in TF object.

  my($h) = "_prom_list_".$f_name."_$random.txt";
  my($name) = $q_gene.$h;	                	# Create the name of the file.

  if (!(-d $out_dir)) {		                        # Check if directory exists. 
    mkdir ($out_dir, 0777) || die "--> Error: unable to create output directory $out_dir\n\n";
    chdir $out_dir;
  }
  else { chdir $out_dir; }

  open P_OUT, ">$name";                                 # Opens a file for output.
  print P_OUT "--> { QUERY GENE: $q_gene } -> list of promoters:\n\n";
  print P_OUT "--> Data set in use: $self->{parameters}{Finalized_db}\n\n";

  my($e)=0;
  foreach my $key (sort keys %$TF_hash) {
    my(@sorted) = sort { $a->{offset} <=> $b->{offset} } @{ $TF_hash->{$key} };       

    my($l)=0;
    while ($l < @{ $TF_hash->{$key} }) {
      my($tf) = $sorted[$l];

      $factor = $tf->{factor}.$tf->{orient};
      $offset = $tf->{offset};
      $motif = $tf->{motif};
      $score = $tf->{score};

      $binding = $tf->{binding}; 
      $condition = $tf->{condition}; 
      $conserv = $tf->{conservation}; 

      $chrom = $tf->{chromosome};
      $prom = $tf->{promoter};

      $i = ">".$e;
      write P_OUT;
     
      $e++; 
      $l++;
    }
  }
  close P_OUT;
  print "\n--> $name - written.";
}

format P_OUT =
@<<<   @<   @<<<<<<    @<<<<<<<  @<<<   @<<<<<<<<<<<<<<<<<<<<   @<<<<<<<<<<<<<   @<<<<<<<   @<   @<<<<<<<<<<<<<<<<   
$i, $chrom, $prom, $factor, $offset, $score, $binding, $condition, $conserv, $motif
.

#----------------------------------- Get_matrix sub. ------------------------------------

# 17. &get_matrix subroutine to retrieve binding affinity matrix of yeast. 
                                                                                        
sub get_matrix() {
  my($self) = shift;

  my($m_file) = $self->{parameters}{Affinity_file}; 
  my($base_dir) = $self->{parameters}{Base_dir}; 
  my($out_dir) = $self->{parameters}{Out_dir}; 

  chdir $base_dir;
  if (-f $m_file) { 
    open(INFO, $m_file) || die "--> Cannot open binding affinity file: $m_file\n"; 
  }
           
  my(@matrix);
  while(<INFO>) { chomp; push(@matrix, $_); }
  close INFO;

  my(@row_1) = split(/\,/, $matrix[0]);        		# Put matrix column heads into an array, @row_1.

  my(%tf_conds);
  for (my $j=1; $j<@row_1; $j++) {
    my($tf, $cond) = split(/\_/, $row_1[$j]);

    $tf_conds{$tf}{$cond} = $j;				# Build TF-condition-position documentation hash.
  }

  my(%proms);
  for (my $i=1; $i<@matrix; $i++) {
    my(@row_p) = split(/\s/, $matrix[$i]);     		# Split promoter row into -> promoter name & entries.
    my($prom) = $row_p[0];				# Pre-promoter name.
 
    if ($prom =~ /[+]/) {
      my(@ps) = split(/\+/, $prom);

      foreach my $p (@ps) {
        $p = &DNA::Func::get_prom_name($p); 
        @{ $proms{$p} } = split(/\,/, $row_p[1]); 
      }
    }
    else {
      my($p) = &DNA::Func::get_prom_name($prom); 
      @{ $proms{$p} } = split(/\,/, $row_p[1]);         # Create hash table of (promoter)->(binding affinity values) : (key)->(value).
    }
  }
  return (\%tf_conds, \%proms);				# Return column head and promoter row hash references. 
}

#------------------------------- Final_map_output sub. --------------------------------

# 18. &final_map_output subroutine to record finalized motif database content.

sub final_map_output($) {
  my($self, $TF_hash) = @_;

  my($name) = "finalized.db"; 			 	   	    # Finalized version.
  
  my($out_dir) = $self->{parameters}{Out_dir};

  if (!(-d $out_dir)) {		   	                            # Check if directory exists. 
    mkdir ($out_dir, 0777) || die "--> Error: unable to create output directory $out_dir\n\n";
    chdir $out_dir;
  }
  else { chdir $out_dir; }

  open OUTPUT, ">$name";
  print OUTPUT "--> [Transcription factor] database - promoter list by TF:\n";
  
  foreach my $key (sort keys %$TF_hash) {

    if ($key ne " ") {
      my($val) = scalar @{ $TF_hash->{$key} };
      print OUTPUT "\n--> PR: $key -> number of targets: $val\n";

      @{ $TF_hash->{$key} } = sort { $a->{offset} <=> $b->{offset} } @{ $TF_hash->{$key} };
      @{ $TF_hash->{$key} } = sort { $a->{factor} cmp $b->{factor} } @{ $TF_hash->{$key} };
       
      my($l)=0;
      while ($l < $val) {
        my($tf) = $TF_hash->{$key}[$l];

        my($facto) = $tf->{factor};
        my($offse) = $tf->{offset};
        my($orien) = $tf->{orient};
        my($motif) = $tf->{motif};
        my($score) = $tf->{score};

        my($chrom) = $tf->{chromosome};
        my($promo) = $tf->{promoter};

        my($condi) = $tf->{condition};
        my($bindi) = $tf->{binding};
        my($strin) = $tf->{stringency};

        print OUTPUT "$facto$orien\t$offse\t$score\t$bindi\t$condi\t$strin\t$motif\t$chrom\t$promo\n";
        $l++;
      }
      print OUTPUT "______________________________________________________________________________________________\n";
    }
  } 
  close OUTPUT;
  print "\n--> $name - written.";
}

#---------------------------- Binding_summary_output sub. -----------------------------

# 19. &binding_summary_output subroutine to record TF binding affinity statistics.

sub binding_summary_output($) {
  my($self, $TF_hash) = @_;

  my($base_dir) = $self->{parameters}{Base_dir};
  my($out_dir) = $self->{parameters}{Out_dir};
  my($c_file) = $self->{parameters}{Conv_file};

  my($description) = &DNA::Func::get_description($base_dir, $c_file);  # Create gene description map.

  if (!(-d $out_dir)) {		   	                            # Check if directory exists. 
    mkdir ($out_dir, 0777) || die "--> Error: unable to create output directory $out_dir\n\n";
    chdir $out_dir;
  }
  else { chdir $out_dir; }

  my($name) = "binding_summary.txt"; 		 	   	    # Binding summary file.

  open OUTPUT, ">$name";
  print OUTPUT "--> Binding affinity summary - TF list by promoter:\n\n";
  
  my(@proms) = sort { $TF_hash->{$a}->{gene} cmp 
                      $TF_hash->{$b}->{gene} } keys %$TF_hash;
  my($l)=0;
  foreach my $key (@proms) {
      my($prom_obj) = $TF_hash->{$key};
      
      my(@tf_cond);
      my($tf_cond);
      my($hits, $occurs);

      if (defined $prom_obj->{tf_cond_list}) {
        @tf_cond = @{ $prom_obj->{tf_cond_list} };
        $tf_cond = &DNA::Func::get_tf_cond_list(@tf_cond);
        $hits = scalar @tf_cond;
        $occurs = &DNA::Func::get_occurence(@tf_cond);
      }
      else { 
        $tf_cond = 'no_hit'; 
        $hits = 0; 
        $occurs = 0; 
      }
    
      my($prom_l) = length $prom_obj->{prom_seq};

      my($prom) = $prom_obj->{promoter};
      my($chrom) = $prom_obj->{chromosome};

      my($gene) = $prom_obj->{gene};

      my($descr);
      if ($gene==1) { $gene = 'g_DNE'; $descr = 'd_DNE'; }
      else {
        $gene = $prom_obj->{gene};
        $descr = $description->{$gene};
      }
      
      print OUTPUT "$l\t$chrom\t$gene\t[$prom]\t$prom_l\t$tf_cond\t$hits\t$occurs\t$descr\n";
      $l++;
  }
  print OUTPUT "_________________________________________________________\n";

  if (defined $self->{no_eq}) {
    print OUTPUT "\n--> Promoters in list, but not in Young's matrix:\n\n";

    my($n)=0;
    while ($n < @{$self->{no_eq}}) {                    # Access 'no_eq' field in promoter object.
      my($rej) = $self->{no_eq}[$n]; 
    
      my($prom_l) = length $rej->{prom_seq};
      my($prom) = $rej->{promoter};
      my($chrom) = $rej->{chromosome};

      my($tf_cond) = 'NDEF'; 
      my($hits) = 0; 
      my($occurs) = 0;

      my($gene) = $rej->{gene};

      my($descr);
      if ($gene==1) { $gene = 'g_DNE'; $descr = 'd_DNE'; }
      else {
        $gene = $rej->{gene};
        $descr = $description->{$gene};
      }

      print OUTPUT "$n\t$chrom\t$gene\t[$prom]\t$prom_l\t$tf_cond\t$hits\t$occurs\t$descr\n";
      $n++;
    }
  }
  close OUTPUT;
  print "\n--> $name - written.";
}

#-------------------------------- Global_hg_result sub. -------------------------------

# 20. &global_hg_result subroutine to output global hyper-geometric comparison results.

sub global_hg_result($$$) {
  my($self, $hyper_hash, $mode, $max) = @_;
  
  my($base_dir) = $self->{parameters}{Base_dir};
  my($c_file) = $self->{parameters}{Conv_file};
  my($out_dir) = $self->{parameters}{Out_dir};

  my($random) = int( rand(9E5) ) + 1E5;				   # Generate random number between 1E5 and 999999.
  $self->{bar_code} = $random;					   # Store it in Stem object.

  if (!(-d $out_dir)) {		                                   # Check if directory exists. 
    mkdir ($out_dir, 0777) || die "--> Error: unable to create output directory: $out_dir\n\n";
    chdir $out_dir;
  }
  else { chdir $out_dir; }

  my($name) = "global_hg_comp_$random.txt";

  open DOC, ">$name";                                              # Open output file.                                                                                        
  print DOC "--> Hyper-geometric comparison results\n";
  print DOC "--> Data set in use: $self->{parameters}{Finalized_db}\n";
  print DOC "--> Number of promoters in yeast: ";

  my(@ordered) = sort { $hyper_hash->{$a}->{num_set1} <=>
                        $hyper_hash->{$b}->{num_set1} } keys %$hyper_hash; 
  my($i)=0;  
  LP: foreach my $key (sort { $hyper_hash->{$a}->{result} <=>
                              $hyper_hash->{$b}->{result} } @ordered) { 
    
    my($combi) = $hyper_hash->{$key};				   # Get current Combi object.

    if ($i==0) { print DOC "$combi->{no_proms}\n\n"; }		   # Number of promoters in yeast.
    my($num_set1) = $combi->{num_set1};				   # Number of promoters bound by gene1.
    my($num_set2) = $combi->{num_set2};				   # Number of promoters bound by gene2.
    my($num_isect) = $combi->{num_isect};			   # Number of promoters in intersection.
  
    my($result) = $combi->{result};				   # Hyper-geometric value in direction 1. 
  
    my(@genes) = @{ $combi->{genes} };				   # Gene pair of query.
    my(@isect) = @{ $combi->{isect} };			           # List of intersecting promoters.
  
    my($gene1) = $genes[0];					   # Query gene 1.
    my($gene2) = $genes[1];					   # Query gene 2.

    print DOC "$i\t$gene1: $num_set1\t$gene2: $num_set2\ti: $num_isect\t$result";
    
    if ($mode eq 'i') {
      my($proms);
      foreach my $p (@isect) {
        my($g);
        if ($p ne 'DNE') {
          $g = &DNA::Func::promoter_to_gene($p, $c_file, $base_dir);  # Convert promoter to gene.
     
          if ($g==1) { $g = '-'; } 
          $proms = $proms."$p [$g], ";
        }
      }
      chop $proms; chop $proms;
      print DOC "\t$proms\n";
    }
    elsif ($mode eq '-') { 
      print DOC "\n";
      if ($i==0) { print "\n--> No promoter list in intersection requested -> 'mode': '-'"; } 
    }
    else { 
      print "--> Error: no such mode: $mode, usage: i|-\n\n";
      exit;
    }
    
    if (defined $max) { 
      if ($i < $max) { $i++; } else { last LP; }
    }
    else { $i++; }
  }
  close DOC;                                                       # Close output file.
  print "\n--> $name - written."; 
}

#------------------------------- Remove_units_output sub. -----------------------------

# 21. &remove_units_output subroutine to record promoter sequences where selected 
# units are removed.

sub remove_units_output() {
  my($self) = shift;

  if (defined @{ $self->{removed} }) {
  
    my($f_name) = $self->{parameters}{File_name};  
    my($out_dir) = $self->{parameters}{Out_dir}; 
  
    my($gene) = $self->{prom_obj}->{gene};                 # Get gene name.
    my($prom) = $self->{prom_obj}->{promoter};             # Get promoter name.
    my($chro) = $self->{prom_obj}->{chromosome};           # Get chromosome no.

    my($random) = $self->{bar_code};			   # Bar code of unit object.

    my($name) = $gene."_units_removed_".$f_name.	   # Promoter mutation with removed units.
                             "_".$random.".txt";

    my($d_name) = $gene."_units_removed_".$f_name.	   # Mutation documentation file.
                               "_".$random.".doc";
              
    if (!(-d $out_dir)) {		                   # Check if directory exists. 
      mkdir ($out_dir, 0777) || die "--> Error: unable to create output directory $out_dir\n\n";
      chdir $out_dir;
    }
    else { chdir $out_dir; }

    open OUTPUT, ">$name";					  
    open DOC, ">$d_name";
					  
    print OUTPUT "--> { GENE: $gene <=> PROMOTER: $prom } -> CHROMOSOME: $chro\n\n";
    print DOC "--> { GENE: $gene <=> PROMOTER: $prom } -> CHROMOSOME: $chro\n\n";
    print DOC "--> Data set in use: $self->{parameters}{Finalized_db}\n";

    my($mode);
    if ($self->{mode}==0) { $mode = "permuted"; }
    else { $mode = "randomized"; }

    my($r)=0; 
    foreach my $m_prom (@{ $self->{removed} }) {           # Access 'removed' field in Mutant object.
      print OUTPUT ">$r: $m_prom->{mutant}\n";
      
      my(@unit_list) = split(/\s/, $m_prom->{unit_selection});
      print DOC "\n--> Removed unit numbers: @{ $m_prom->{unit_selection} }\n";

      foreach my $key (sort { $a <=> $b } keys %$m_prom) {
        if ($key =~ /^\d+$/) {
          my($unit) = $m_prom->{$key};

          my($u_name) = $m_prom->{unit_hash}{$key}->{factors}[0]->{factor};
          my($e_numb) = scalar @{ $m_prom->{unit_hash}{$key}->{factors} };

          print DOC "\n--> UNIT [$key]: $u_name - number of elements: $e_numb\n";
        
          for (my $k=0; $k < @{ $unit->{o_kmer} }; $k++) {     

            my($kmer) = $unit->{o_kmer}->[$k];         	   # Assign original kmer.
            my($kmer_l) = $unit->{kmer_l}->[$k];       	   # Assign kmer length.
            my($s_kmer) = $unit->{s_kmer}->[$k];       	   # Assign substitition kmer.
            my($offset) = $unit->{s_offset}->[$k];     	   # Assign substitution offset.

            print DOC ">$k: Replace $kmer ($offset) with $mode kmer $s_kmer ($kmer_l).\n";
          }
        }
      }
      $r++;
      print DOC "____________________________________________________________________________\n";
    }
    close OUTPUT;
    close DOC;

    print "\n--> $name - written."; 
    print "\n--> $d_name - written."; 
  }
  else { 
    print "\n--> No UNIT object present!"; 
    print "\n--> Process exited - remove_units.pl\n\n"; 
    exit;
  }
}

#------------------------------- Random_kmer_output sub. -----------------------------

# 22. &random_kmer_output subroutine to record promoter sequences where selected 
# units are removed.

sub random_kmer_output() {
  my($self) = shift;

  my($out_dir) = $self->{parameters}{Out_dir}; 

  if (!(-d $out_dir)) {		                             # Check if directory exists. 
    mkdir ($out_dir, 0777) || die "--> Error: unable to create output directory: $out_dir\n\n";
    chdir $out_dir;
  }
  else { chdir $out_dir; }

  my($name) = "random_kmer_library.txt";
  open OUTPUT, ">$name";                                     # Open output file.

  foreach my $kmer_l (sort { $a <=> $b } keys %$self) {
    if ($kmer_l =~ /^\d+$/) {
      print OUTPUT "\n<----- Kmer length: $kmer_l ----->\n";

      my($j)=0;
      foreach my $kmer ( @{ $self->{$kmer_l} } ) {  
        print OUTPUT "$j\t$kmer\n";
        $j++;
      }
    }
  }

  close OUTPUT;  
  print "\n--> $name - written."; 
}

#------------------------------------ Constants ---------------------------------------
                                                                                        
# 23. constant hash random kmer library of various lenghts.

use constant R_KMER => ( 5  => [ "ATCGT", "ATCGA", "ATGCA", "GACGA", "AGTAC", "AATAT", "ATCGC", "ATCGG", "ATGCA", "GACGT" ], 
        		 6  => [ "ATCGTC", "ATCGAT", "ATGCAT", "GACGAT", "AGTACT", "AATATT", "ATCGTT", "ATCGAG", "ATGCAA", "GACGAG" ], 
        		 7  => [ "TAAGGCT", "ATGTCCA", "GAATACA", "AATTCAC", "CAACACA", "TGTATTC", "CATGGTC", "AATAGTT", "ATGCGAA", "GGAACTA" ],
        		 8  => [ "TATCTTTA", "AGTAGATG", "AGTAGAGT", "AGCGCTAA", "CTACGACT", "TAACGATT", "CTACGACT", "AATCTTTG", "ATACGACT", "CATCTACT" ],
        		 9  => [ "GTAAAATTG", "ACGGTTTGT", "ACCGACTGC", "CTTAAGCTG", "CAAGTCGAC", "GAATGATTG", "ATACAGTTT", "GTTTTACTC", "GTGTCATTC", "GTAGAAATG" ],
        		 10 => [ "TCCATGACCC", "GGTAGATTTA", "GATTCACTAG", "GACAGTACTA", "TAACAGTACT", "GAGACGGTTG", "AAGACGGTCA", "ACGAAGATCA", "CATCTCTTAT", "ACTACGATAT" ],
  			 11 => [ "CACCTTGCTGT", "AGTATAGACCC", "TAAAGTTTACA", "TAGAGCTATAG", "AGGGTCATTCA", "AGTAGCTCTGG", "ACGAAGATCTG", "CTAAACTTACA", "GTCGTAAGACT", "TTGAACTATGA" ],
  			 12 => [ "CTAAGATCTTGT", "CGAGATCGTAGT", "CTAGATCTATGG", "TACTATAGACCT", "TGTCTTAAGACA", "AGGTCTATAGAG", "AGATCTTTTTTT", "AGATCGTTCATG", "TATAGATCTACG", "AAGTAGATCTAT" ],
  			 13 => [ "CTAAGATCTTGTT", "CCGAGATCGTAGT", "AGCTAGATCTATG", "TACTATAGACCTT", "TGTCTTAAGACAA", "AGGTCTATAGAGG", "AGATCTTTTTTTC", "AGATCGTTCATGC", "TATAGATCTACGG", "AAGTAGATCTATT" ],
  			 14 => [ "ATCTAAGATCTTGT", "GACGAGATCGTAGT", "AACTAGATCTATGG", "GATACTATAGACCT", "ATTGTCTTAAGACA", "AAAGGTCTATAGAG", "ATAGATCTTTTTTT", "AAAGATCGTTCATG", "ATTATAGATCTACG", "GAAAGTAGATCTAT" ],
  			 15 => [ "ATCCTAAGATCTTGT", "ATGCGAGATCGTAGT", "GACCTAGATCTATGG", "AGTTACTATAGACCT", "ATGTGTCTTAAGACA", "ATCAGGTCTATAGAG", "GACAGATCTTTTTTT", "ATGAGATCGTTCATG", "ATCTATAGATCTACG", "GACAAGTAGATCTAT" ],
  			 16 => [ "ATCGCTAAGATCTTGT", "GACGCGAGATCGTAGT", "AGTACTAGATCTATGG", "AATATACTATAGACCT", "ATGCTGTCTTAAGACA", "GACGAGGTCTATAGAG", "ATGCAGATCTTTTTTT", "ATCGAGATCGTTCATG", "ATCGTATAGATCTACG", "GACGAAGTAGATCTAT" ],
  			 17 => [ "ATCGTCTAAGATCTTGT", "ATCGACGAGATCGTAGT", "ATGCACTAGATCTATGG", "GACGATACTATAGACCT", "AGTACTGTCTTAAGACA", "AATATAGGTCTATAGAG", "ATCGCAGATCTTTTTTT", "ATGCAAGATCGTTCATG", "ATGCATATAGATCTACG", "GACGTAAGTAGATCTAT" ],
);

#-------------------------------- Ortho_output_i sub. ----------------------------------

# 24. &ortho_output_i subroutine to record orthogroups of S. cerevisiae sorted by one
# promoter at a time.

sub ortho_output_i() {
  my($self) = shift;

  my($base_dir) = $self->{parameters}{Base_dir};
  my($out_dir) = $self->{parameters}{Out_dir};
  my($conv_file) = $self->{parameters}{Conv_file};

  my($name) = "cer_orthogroups_i.txt";

  if (!(-d $out_dir)) {		   	                            # Check if directory exists. 
    mkdir ($out_dir, 0777) || die "--> Error: unable to create output directory $out_dir\n\n";
    chdir $out_dir;
  }
  else { chdir $out_dir; }

  open ORTH, ">$name";
  print ORTH "--> Cerevisiae orthogroup analysis: $self->{parameters}{Ortho_file}\n\n";

  my($i)=0;
  foreach my $o_hash (@{ $self->{ortho_list} }) {
    my($cers) = scalar @{ $o_hash->{'cer'} };

    for (my $j=0; $j<$cers; $j++) {
      my($cer) = $o_hash->{'cer'}[$j];

      my(@prom) = split(/\|/, $cer);
      my($gene) = &DNA::Func::promoter_to_gene($prom[1], $conv_file, $base_dir);  # Convert promoter to gene.

      my($proms);
      foreach my $specie (sort keys %$o_hash) {
        if ($specie ne "cer") {
          my($spes) = scalar @{ $o_hash->{$specie} };
 
          for (my $k=0; $k<$spes; $k++) {
            my($other) = $o_hash->{$specie}[$k];
            $proms = $proms."$other,";
          }
          chop $proms;
          $proms = $proms."\t";
        } 
      }
      chop $proms;
      print ORTH ">$i\t$gene\t$cer\t$proms\n";
    }
    $i++;
  } 
  close ORTH;
  print "\n--> $name - written.";
}

#-------------------------------- Ortho_output_o sub. ----------------------------------

# 25. &ortho_output_o subroutine to record orthogroups of S. cerevisiae sorted by 
# appearance in original file.

sub ortho_output_o($) {
  my($self) = shift;

  my($base_dir) = $self->{parameters}{Base_dir};
  my($out_dir) = $self->{parameters}{Out_dir};
  my($conv_file) = $self->{parameters}{Conv_file};

  my($name) = "cer_orthogroups_o.txt";

  if (!(-d $out_dir)) {		   	                            # Check if directory exists. 
    mkdir ($out_dir, 0777) || die "--> Error: unable to create output directory $out_dir\n\n";
    chdir $out_dir;
  }
  else { chdir $out_dir; }

  open ORTH, ">$name";
  print ORTH "--> Cerevisiae orthogroup analysis: $self->{parameters}{Ortho_file}\n\n";

  my($i)=0;
  foreach my $o_hash (@{ $self->{ortho_list} }) {
    my($cers) = scalar @{ $o_hash->{'cer'} };

    my($genes);
    for (my $j=0; $j<$cers; $j++) {
      my($cer) = $o_hash->{'cer'}[$j];

      my(@prom) = split(/\|/, $cer);
      my($gene) = &DNA::Func::promoter_to_gene($prom[1], $conv_file, $base_dir);  # Convert promoter to gene.

      $genes = $genes."$gene,";
    }
    chop $genes;

    my($cers);
    my($proms);

    foreach my $specie (sort keys %$o_hash) {
      if ($specie ne "cer") {
        my($spes) = scalar @{ $o_hash->{$specie} };
 
        for (my $k=0; $k<$spes; $k++) {
          my($other) = $o_hash->{$specie}[$k];
          $proms = $proms."$other,";
        }
        chop $proms;
        $proms = $proms."\t";
      }
      else {
        my($spes) = scalar @{ $o_hash->{$specie} };
 
        for (my $k=0; $k<$spes; $k++) {
          my($cer) = $o_hash->{$specie}[$k];
          $cers = $cers."$cer,";
        }
        chop $cers;
        $cers = $cers."\t";
      } 
    } 
    print ORTH ">$i\t$genes\t$cers$proms\n";

    $i++;
  }
  close ORTH;
  print "\n--> $name - written.";
}


#-------------------------------- Ortho_table_output sub. -----------------------------

# 26. &ortho_table_output subroutine to output TFBS occurence comparision table of 
# gene (promoter) of interest.

sub ortho_table_output() {
  my($self) = shift;

  my($q_gene) = $self->{gene}; 	 		                    # Query gene of interest.
  my($q_prom) = $self->{promoter}; 	 		            # Equaivalent promoter name.

  my(@table) = @{ $self->{ortho_table} }; 	                    # TF comparision table.

  my($out_dir) = $self->{parameters}{Out_dir};
  my($name) = $q_gene."_orthogroups_table.txt";

  if (!(-d $out_dir)) {		   	                            # Check if directory exists. 
    mkdir ($out_dir, 0777) || die "--> Error: unable to create output directory $out_dir\n\n";
    chdir $out_dir;
  }
  else { chdir $out_dir; }

  open ORTH, ">$name";
  print ORTH "--> { Saccharomyces cerevisiae - GENE: $q_gene <=> PROMOTER: $q_prom }\n";
  print ORTH "--> Cerevisiae orthogroup analysis: $self->{parameters}{Filtered_orthogroups}\n";

  foreach my $row (@table) { print ORTH "$row\n"; }
  close ORTH;

  print "\n--> $name - written.";
}

#--------------------------------- List_table sub. ----------------------------------

# 27. &list_table subroutine to record TFBS occurences in all orthologous promoters 
# of yeast organisms in the orthogroup.

sub list_table() {
  my($self) = shift;

  my($q_gene) = $self->{gene}; 	 		                    # Query gene of interest.
  my($q_prom) = $self->{promoter}; 	 		            # Equaivalent promoter name.

  my($out_dir) = $self->{parameters}{Out_dir};
  my($name) = $q_gene."_orthogroups_table_i.txt";

  if (!(-d $out_dir)) {		   	                            # Check if directory exists. 
    mkdir ($out_dir, 0777) || die "--> Error: unable to create output directory $out_dir\n\n";
    chdir $out_dir;
  }
  else { chdir $out_dir; }

  open ORTH, ">$name";
  print ORTH "--> { Saccharomyces cerevisiae - GENE: $q_gene <=> PROMOTER: $q_prom }\n";
  print ORTH "--> Cerevisiae orthogroup listing: $self->{parameters}{Filtered_orthogroups}\n";

  my(@phylo_tree) = ('cer', 'par', 'mik', 'bay', 'gla', 
		     'cas', 'lac', 'han', 'alb', 'lip');

  foreach my $org (@phylo_tree) {
    if (exists $self->{factors_list}{$org}) {
      foreach my $factors ( @{ $self->{factors_list}{$org} } ) {
        if ( !(defined $factors->{prom_obj}->{no_eq}) ) {

          my($gene) = $factors->{prom_obj}->{gene};
          my($prom) = $factors->{prom_obj}->{promoter};

          print ORTH "\n***** ORGANISM: $org *****";
          print ORTH "\n***** PROMOTER: $prom *****";
          print ORTH "\n***** GENE: $gene *****";

          my(%o) = %{ $factors->{tf_occurence} };
          print ORTH "\n\n";
          foreach my $k (sort keys %o) { print ORTH "$k => $o{$k}\n"; }
        }
        else { print ORTH "\n--> No data for org: $org <=> prom: $factors->{prom_obj}->{no_eq}\n\n"; }
      }
    }
  }
  close ORTH;
  print "\n--> $name - written.";
}

#----------------------------- Excluded_proms_IO sub. ----------------------------------

# 28. &excluded_proms_IO subroutine to record promoters which either: a) are not in the
# promoter sequence database or b) do not have orthogroups.

sub excluded_proms_IO($$) {
  my($out_dir, $no_hits) = @_;

  if (!(-d $out_dir)) {                                     # Check if directory exists.
    mkdir ($out_dir, 0777) || die "--> Error: unable to create output directory $out_dir\n\n";
    chdir $out_dir;
  }
  else { chdir $out_dir; }

  my($out_file) = "excluded_promoters.txt";
  open NO_HIT, ">$out_file";

  foreach my $key (keys %$no_hits) {
    if ($key==1) {
      print NO_HIT "--> No data available for list of (genes <-> promoters):\n\n";
    }
    elsif ($key==2) {
      print NO_HIT "--> No orthologs available for list of (genes <-> promoters):\n\n";
    }
    else { die "\n-->Error: no such exclusion member - $key\n\n"; }

    my($i)=0;
    foreach my $e (sort { $a->{gene} <=> $b->{gene} }
                                @{ $no_hits->{$key} }) {
      my($gene) = $e->{gene};
      my($prom) = $e->{promoter};

      print NO_HIT ">$i\t$gene\t$prom\n";
    }
  }
  close NO_HIT;
}

#---------------------------- Statistics_output sub. ----------------------------------

# 29. &statistics_output subroutine to record all TF pair distances for all promoters
# in the S. cerevisiae genome.

sub statistics_output() {
  my($self) = shift;

  my($out_dir) = $self->{parameters}{Out_dir};

  my($d_file) = "relative_distances.txt";
  my($o_file) = "relative_orientation.txt";

  if (!(-d $out_dir)) {		   	                            # Check if directory exists. 
    mkdir ($out_dir, 0777) || die "--> Error: unable to create output directory $out_dir\n\n";
    chdir $out_dir;
  }
  else { chdir $out_dir; }

  open DIST, ">$d_file";
  open ORDER, ">$o_file";

  foreach my $prom_name (keys %{ $self->{tf_pair_list} }) {
    foreach my $pair (@{ $self->{tf_pair_list}{$prom_name} }) {
      print DIST "$pair->{distance}\s";

      my($orient_type) = $pair->{orient_type};

      if ($orient_type eq "++") { print ORDER "1\s"; } 
      if ($orient_type eq "+-") { print ORDER "2\s"; } 
      if ($orient_type eq "--") { print ORDER "3\s"; } 
      if ($orient_type eq "-+") { print ORDER "4\s"; } 
    }
  }
  close DIST;
  close ORDER;
}

#----------------------------- D_to_gene_output sub. ----------------------------------

# 30. &d_to_gene_output subroutine to record all TF-gene distances for all promoters
# in the S. cerevisiae genome.

sub d_to_gene_output() {
  my($self) = shift;

  my($out_dir) = $self->{parameters}{Out_dir};

  my($out_file) = "TF_gene_distances.txt";

  if (!(-d $out_dir)) {		   	                            # Check if directory exists. 
    mkdir ($out_dir, 0777) || die "--> Error: unable to create output directory $out_dir\n\n";
    chdir $out_dir;
  }
  else { chdir $out_dir; }

  open DIST, ">$out_file";

  foreach my $prom_name (keys %{ $self->{tf_gene_dist} }) {
    foreach my $dist (@{ $self->{tf_gene_dist}{$prom_name} }) {
      print DIST "$dist\s";
    }
  }
  close DIST;
}

#---------------------------- O_TF_search_output sub. ---------------------------------

# 31. &o_TF_search_output subroutine to output list of TFs bound to given yeast promoter.

sub o_TF_search_output() {
  my($self) = shift;

  my($f_name) = $self->{parameters}{File_name};  
  my($out_dir) = $self->{parameters}{Out_dir}; 

  my(@fac) = @{ $self->{final_list} };
  my($prom) = $self->{prom_obj}->{promoter};
  my($chro) = $self->{prom_obj}->{chromosome};
  my($org) = $self->{prom_obj}->{gene};

  my($random) = int( rand(9E5) ) + 1E5;			 # Generate random number between 1E5 and 999999.
  $self->{bar_code} = $random;				 # Store it in TF object.

  my($h) = "_human_".$f_name."_$random.txt";
  my($m) = "_machine_".$f_name."_$random.txt";
  my($h_name) = $prom.$h;	                	 # Create the name of the file.
  my($m_name) = $prom.$m;	                	 # Create the name of the file.

  if (!(-d $out_dir)) {		                         # Check if directory exists. 
    mkdir ($out_dir, 0777) || die "--> Error: unable to create output directory $out_dir\n\n";
    chdir $out_dir;
  }
  else { chdir $out_dir; }

  open O_RESULT, ">$h_name";                             # Opens a h_file for output.
  print O_RESULT "--> { ORGANISM: $org <-> PROMOTER: $prom } -> CHROMOSOME: $chro\n\n";

  open M_RESULT, ">$m_name";                             # Opens a m_file for output.

  my($e)=0;
  while ($e < @fac) {
    $i = ">".$e; 
    $factor = $fac[$e]->{factor}.$fac[$e]->{orient}; 
    $offset = $fac[$e]->{offset};
    $motif = $fac[$e]->{motif}; 
    
    my($scr) = $fac[$e]->{score};
    #$score = sprintf("%.10f", $scr);		         # Display 10 digits after decimal. 
    $score = $scr;		  		         # Display scientific notation. 
   
    write O_RESULT;
    print M_RESULT ">$e\t$factor\t$offset\t$score\t$motif\n";
    $e++;
  }
  close O_RESULT;                                 	# Close the h_file.
  close M_RESULT;                                 	# Close the m_file.
  print "\n--> $h_name - written."; 
  print "\n--> $m_name - written.";
}

format O_RESULT =
@<<<  @<<<<<<< @<<<<  @<<<<<<<<<<<<  @<<<<<<<<<<<<<<<<   
$i, $factor, $offset, $score, $motif
.

#--------------------------------- List_ortho_tfs sub. ----------------------------------

# 32. &list_ortho_tfs subroutine to record TFBS data in all orthologous promoters of 
# yeast organisms in the orthogroup.

sub list_ortho_tfs() {
  my($self) = shift;

  my($q_gene) = $self->{gene}; 	 		                    # Query gene of interest.
  my($q_prom) = $self->{promoter}; 	 		            # Equaivalent promoter name.

  my($out_dir) = $self->{parameters}{Out_dir};
  my($name) = $q_gene."_orthogroup_tfs.txt";

  if (!(-d $out_dir)) {		   	                            # Check if directory exists. 
    mkdir ($out_dir, 0777) || die "--> Error: unable to create output directory $out_dir\n\n";
    chdir $out_dir;
  }
  else { chdir $out_dir; }

  open ORTH, ">$name";
  print ORTH "--> { Saccharomyces cerevisiae - GENE: $q_gene <=> PROMOTER: $q_prom }\n";
  print ORTH "--> Cerevisiae orthogroup listing: $self->{parameters}{Filtered_orthogroups}\n";

  my(@phylo_tree) = ('cer', 'par', 'mik', 'bay', 'gla', 
		     'cas', 'lac', 'han', 'alb', 'lip');

  foreach my $org (@phylo_tree) {
    if (exists $self->{factors_list}{$org}) {
      foreach my $factors ( @{ $self->{factors_list}{$org} } ) {
        if ( !(defined $factors->{prom_obj}->{no_eq}) ) {

          my($gene) = $factors->{prom_obj}->{gene};
          my($prom) = $factors->{prom_obj}->{promoter};

          print ORTH "\n***** ORGANISM: $org *****";
          print ORTH "\n***** PROMOTER: $prom *****";
          print ORTH "\n***** GENE: $gene *****\n\n";

          my(@fac) = @{ $factors->{final_list} };

  	  my($e)=0;
  	  while ($e < @fac) {
   	    $i = ">".$e; 
   	    $factor = $fac[$e]->{factor}.$fac[$e]->{orient}; 
   	    $offset = $fac[$e]->{offset};
    	    $motif = $fac[$e]->{motif}; 
    
    	    my($scr) = $fac[$e]->{score};
    	    #$score = sprintf("%.10f", $scr);		         # Display 10 digits after decimal. 
  	    $score = $scr;		  		         # Display scientific notation. 

            if ($org eq 'cer') {
    	      $binding = $fac[$e]->{binding}; 
    	      $condition = $fac[$e]->{condition}; 
    	      $conserv = $fac[$e]->{conservation}; 

              print ORTH ">$e\t$factor\t$offset\t$score\t$binding\t$condition\t$conserv\t$motif\n";
            }
            else { print ORTH ">$e\t$factor\t$offset\t$score\t$motif\n"; }

            $e++;
          }
        }
        else { print ORTH "\n--> No data for org: $org <=> prom: $factors->{prom_obj}->{no_eq}\n"; }
      }
    }
  }
  close ORTH;
  print "\n--> $name - written.";
}

#---------------------------------- Response_IO sub. ----------------------------------

# 33. &response_IO subroutine to generate error message or proceed depending on ortho-
# group response for given gene name.

sub response_IO($$) {
  my($resp, $prom) = @_;

  if ($resp==1) { 
    print "\n--> Error: promoter not present in orthogroup file - $prom\n\n"; 
    exit;
  }
  elsif ($resp==2) {
    print "\n--> Error: no other orthologous promoters in any yeast organism - $prom\n\n"; 
    exit;
  }
  else { print "\n--> OK: promoter - $prom in on its way!"; } 
}

#-------------------------------- tf_pair_output sub. ---------------------------------

# 34. &tf_pair_output subroutine to output TF pairs which occure at least in 50 distinc
# promoters along with its corresponding overlap statistics.

sub tf_pair_output($) {
  my($self, $p) = @_;

  my($diff) = $self->{parameters}{Diff_prom_occurence};             # Number of promoters needed for filtering.

  my($out_dir) = $self->{parameters}{Out_dir};
  my($name) = "tf_pairs.txt";

  if (!(-d $out_dir)) {		   	                            # Check if directory exists. 
    mkdir ($out_dir, 0777) || die "--> Error: unable to create output directory $out_dir\n\n";
    chdir $out_dir;
  }
  else { chdir $out_dir; }

  open PR, ">$name";
  print PR "--> { Saccharomyces cerevisiae - all TF pairs + overlap statistics > $diff }\n";

  my(%dist_prom) = %{ $self->{dist_prom} };			    # Distinct promoter count hash [TF pair].
  my(%OL_prom) = %{ $self->{OL_prom} };			            # Distinct promoter count hash [overlap].

  my($i)=0;
  foreach my $tf_pair (@{ $self->{filtered_pairs} }) {  

    my($tf_p_occ) = scalar keys %{ $dist_prom{$tf_pair} };          # No. of different promoters containing TF pair.
    my($tf_p_max) = $self->{tf_pair}{$tf_pair};			    # No. of occurences of TF pair in genome.

    my($ol_p_o) = scalar keys %{ $OL_prom{$tf_pair} };              # No. of different promoters containing overlap. 

    my($ol_num);						    # Initialize overlap occurence number.

    if ($ol_p_o==0) { $ol_num=0; }
    else { $ol_num = $self->{tf_overlap}{$tf_pair}; }		    # No. of occurences of overlap in genome.

    my($proms);

    if ($p==1) { 
      foreach my $prom (sort {$a cmp $b} keys %{ $dist_prom{$tf_pair} }) {
        $proms = $proms.$prom.',';
      }
      chop $proms;
    }

    my($p_ol_p) = ($ol_num/$tf_p_max)*100;			    # Overlap percentage.
    my($p_num) = sprintf("%.4f", $p_ol_p);	  	            # Display 4 digits after decimal. 

    print PR "\n>$i\t$tf_pair\t$tf_p_occ\t$tf_p_max\t$ol_p_o\t$ol_num\t$p_num%\t$proms";
    $i++;
  }
  close PR;
  print "\n--> $name - written.";
}

#-------------------------------- tf_overlap_output sub. ---------------------------------

# 35. &tf_overlap_output subroutine to output overlapping TF-pairs.

sub tf_overlap_output($) {
  my($self, $p) = @_;

  my($out_dir) = $self->{parameters}{Out_dir};
  my($name) = "overlapping_pairs.txt";

  if (!(-d $out_dir)) {		   	                            # Check if directory exists. 
    mkdir ($out_dir, 0777) || die "--> Error: unable to create output directory $out_dir\n\n";
    chdir $out_dir;
  }
  else { chdir $out_dir; }

  open OL, ">$name";
  print OL "--> { Saccharomyces cerevisiae - all overlapping TF pairs }\n";

  my(%dist_prom) = %{ $self->{dist_prom} };			    # Distinct promoter count hash [TF pair].
  my(%OL_prom) = %{ $self->{OL_prom} };			            # Distinct promoter count hash.

  my($i)=0;

  foreach my $tf_pair (sort { ($self->{tf_overlap}{$b}/$self->{tf_pair}{$b}) <=>
                              ($self->{tf_overlap}{$a}/$self->{tf_pair}{$a})} keys %OL_prom) { 

  #foreach my $tf_pair (sort { keys %{ $OL_prom{$b} } <=> 
  #			      keys %{ $OL_prom{$a} } } keys %OL_prom) { 

    my($tf_p_occ) = scalar keys %{ $dist_prom{$tf_pair} };          # No. of different promoters containing TF pair.
    my($tf_p_max) = $self->{tf_pair}{$tf_pair};			    # No. of occurences of TF pair in genome.

    my($ol_p_o) = scalar keys %{ $OL_prom{$tf_pair} };              # No. of different promoters containing overlap.  
    my($ol_num) = $self->{tf_overlap}{$tf_pair};		    # No. of occurences of overlap in genome.

    my($proms);

    if ($p==1) { 
      foreach my $prom (sort {$a cmp $b} keys %{ $OL_prom{$tf_pair} }) {
        $proms = $proms.$prom.',';
      }
      chop $proms;
    }

    my($p_ol_p) = ($ol_num/$tf_p_max)*100;			    # Overlap percentage.
    my($p_num) = sprintf("%.4f", $p_ol_p);	  	            # Display 4 digits after decimal. 

    print OL "\n>$i\t$tf_pair\t$tf_p_occ\t$tf_p_max\t$ol_p_o\t$ol_num\t$p_num%\t$proms";
    $i++;
  }
  close OL;
  print "\n--> $name - written.";
}

#------------------------------- Hg_ol_db_output sub. ---------------------------------

# 36. &hg_ol_db_output subroutine to record overlap <-> hyper-geometric database content.

sub hg_ol_db_output($) {
  my($self, $hg_ol_pairs) = @_;

  my($out_dir) = $self->{parameters}{Out_dir};
  
  if (!(-d $out_dir)) {		   	                            # Check if directory exists. 
    mkdir ($out_dir, 0777) || die "--> Error: unable to create output directory $out_dir\n\n";
    chdir $out_dir;
  }
  else { chdir $out_dir; }

  my($name) = "hg_overlapping_pairs.txt";

  open OUTPUT, ">$name";
  print OUTPUT "--> Overlapping TF pairs with hyper-geometric comparision data:\n\n";

  my($l)=0;
  foreach my $key (sort { $hg_ol_pairs->{$a}->{hg_comp_value} <=> 
			  $hg_ol_pairs->{$b}->{hg_comp_value} } keys %$hg_ol_pairs) {

    my($pair) = $hg_ol_pairs->{$key};

    my($fac_1) = $pair->{position_1}->{factor};
    my($num_p1) = $pair->{position_1}->{bound_promoters};

    if ($fac_1 eq "") { next; }

    my($fac_2) = $pair->{position_2}->{factor};
    my($num_p2) = $pair->{position_2}->{bound_promoters};

    my($hg_val) = $pair->{hg_comp_value};

    print OUTPUT ">$l\t$fac_1-$fac_2\t$num_p1\t$num_p2\t$hg_val\n"; 
    $l++;
  }
  close OUTPUT;
  print "\n--> $name - written.";

  if (exists $self->{no_hg_result}) {
    my($name) = "problematic_pairs.txt";

    open PROB, ">$name";
    print PROB "--> Problematic TF pairs with no hyper-geometric comparision data:\n\n";

    my($m)=0;
    foreach my $pair (sort {$a cmp $b} @{ $self->{no_hg_result} } ) {
      print PROB ">$l\t$pair\n"; 
      $m++;
    }
    close PROB;
    print "\n--> $name - written.";
  }
}

#------------------------------ Ol_tf_ratio_output sub. -------------------------------

# 37. &ol_tf_ratio_output subroutine to calculate overlap-to-pair occurence ratio.
                                                                          
sub ol_tf_ratio_output(@) {
  my($self, @ol_tf_ratio) = @_;

  my($out_dir) = $self->{parameters}{Out_dir};
  
  if (!(-d $out_dir)) {		   	                        # Check if directory exists. 
    mkdir ($out_dir, 0777) || die "--> Error: unable to create output directory $out_dir\n\n";
    chdir $out_dir;
  }
  else { chdir $out_dir; }

  my($name) = "overlap_pair_ratio.txt";

  open OUTPUT, ">$name";
  print OUTPUT "--> Overlap-to-pair occurence comparision data:\n\n";

  for (my $i=0; $i<@ol_tf_ratio-1; $i++) {
    my($pair) = $ol_tf_ratio[$i];				# Current TF pair (overlap).

    my($pair_name) = $pair->{pair_name};
    my($d_proms) = $pair->{distinct_proms};			# No. of distinct promoters pair occurs in.
    my($o_proms) = $pair->{overlap_proms};			# No. of distinct promoters overlap occurs in.
    my($o_num) = $pair->{overlap_num};			        # No. of times overlap occurs.
    my($ratio) = $pair->{overlap_ratio};			# No. of times overlap occurs.

    print OUTPUT ">$i\t$pair_name\t$d_proms\t$o_proms\t$o_num\t[$o_proms/$d_proms]=$ratio%\n";
  }
  close OUTPUT;
  print "\n--> $name - written.";
} 

#------------------------------ Binding_ol_output sub. --------------------------------

# 38. &binding_ol_output subroutine to record overlap <-> binding strength database.

sub binding_ol_output($) {
  my($self, $binding_ol) = @_;

  my($out_dir) = $self->{parameters}{Out_dir};
  
  if (!(-d $out_dir)) {		   	                        # Check if directory exists. 
    mkdir ($out_dir, 0777) || die "--> Error: unable to create output directory $out_dir\n\n";
    chdir $out_dir;
  }
  else { chdir $out_dir; }

  my($name) = "binding_strength_overlap.txt";

  open OUTPUT, ">$name";
  print OUTPUT "--> Overlapping TF pairs with binding affinity <-> condition data:\n\n";

  my($i)=0;
  foreach my $key (sort { $binding_ol->{$b}->{ratio} <=> 
			  $binding_ol->{$a}->{ratio} } keys %$binding_ol) {

    my($pair) = $binding_ol->{$key};

    my($pair_name) = $pair->{pair_name};
    my($d_proms) = $pair->{distinct_proms};			# No. of distinct promoters pair occurs in.
    my($o_proms) = $pair->{overlap_proms};			# No. of distinct promoters overlap occurs in.
    my($o_num) = $pair->{overlap_num};			        # No. of times overlap occurs.
    my($ratio) = $pair->{overlap_ratio};			# No. of times overlap occurs.

    my($cond_1) = $pair->{position_1}->{condition};
    my($bind_1) = $pair->{position_1}->{binding};

    my($cond_2) = $pair->{position_2}->{condition};
    my($bind_2) = $pair->{position_2}->{binding};

    # ***** For same condition/binding? *****

    print OUTPUT ">$i\t$pair_name\t$d_proms\t$o_proms\t$o_num\t[$o_proms/$d_proms]=$ratio%\t$cond_1,$bind_1\t$cond_2,$bind_2\n";
    $i++;
  }
  close OUTPUT;
  print "\n--> $name - written.";
}

#----------------------------- Kl_divergence_output sub. ------------------------------

# 39. &kl_divergence_output subroutine to record all TF pair <-> KL-divergence score.

sub kl_divergence_output($) {
  my($self, $kl_hash) = @_;

  my($out_dir) = $self->{parameters}{Out_dir};
  
  if (!(-d $out_dir)) {		   	                        # Check if directory exists. 
    mkdir ($out_dir, 0777) || die "--> Error: unable to create output directory $out_dir\n\n";
    chdir $out_dir;
  }
  else { chdir $out_dir; }

  my($random) = int( rand(9E5) ) + 1E5;			        # Generate random number between 1E5 and 999999.
  $self->{bar_code} = $random;				        # Store it in promorter object.

  my($name) = "kl_divergence_".$random.".txt";

  open OUTPUT, ">$name";
  print OUTPUT "--> KL-divergence score for each TF pair in promoters:\n\n";

  my($i)=0;
  foreach my $key (sort {$a cmp $b} keys %$kl_hash) {	        # Sort by TF pair name.	
  #foreach my $key (sort { $kl_hash->{$a} <=> $kl_hash->{$b} } keys %$kl_hash) {  # Sort by increasing KL-divergence score.
    print OUTPUT ">$i\t$key\t$kl_hash->{$key}\n";
    $i++;
  }
  close OUTPUT;
  print "\n--> $name - written.";
}
