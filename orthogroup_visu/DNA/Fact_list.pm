#----------------------------------------------------------------------
# Author: Mirko Palla.
# Date: March 2, 2006.
# For: Task 3 at the Church Lab - Genetics Department, 
# Harvard Medical School.
# 
# Purpose: This program contains the complete code for 
# class Fact_list, modeling an object for TF data storage 
# in Perl.
#
# Associated subroutines: 0.  &constructor      1.  &TF_search
#			  2.  &gene_search      3.  &ordering         
#			  4.  &data_set_assign  5.  &all_tf           
#			  6.  &motif_db_ernest  7.  &promoter_db      
#			  8.  &final_map        9.  &count_TF_occurence
#			  10. &o_TF_search      11. &selection_filter
#			  12. &correct_motif    13. &union_filter
#		          14. &hg_pair_list     15. &mesh_hg_ol
#		          16. &ol_tf_ratio      17. &binding_ol
#			  18. &kl_divergence
#---------------------------------------------------------------------- 

package DNA::Fact_list;

use strict;
use vars qw($rev);
use base ("DNA::Fact", "DNA::IO", "DNA::Combi", "DNA::TF_pair");

1;

#--------------------------------- Class Fact_list -------------------------------------

# 0. constructor for transcription factor list object.

sub new {
  my($class, $prom) = @_;
  my($self)={};

  $self->{prom_obj} = $prom;		 			# Promoter object field.
  $self->{final_list} = [];                                     # Final TF list.

  $self->{tf_occurence} = {};	 			        # Hash library of TF occurence.

  $self->{parameters} = {};                                     # Input parameters.
  $self->{bar_code} = ();                                       # Bar code field to 
  			                                        # label output files.
  bless($self,$class);
  return($self);
}

#----------------------------------- TF_search sub. -----------------------------------

# 1. &TF_search subroutine to collect standarized, bound TF data given a promoter.
                                                                                        
sub TF_search($) {
  my($self, $order) = @_;
 
  my($final_db) = $self->{parameters}{Finalized_db};	         # Finalized data base name.  
   
  my($prom_obj) = $self->{prom_obj};			         # Promoter object.
  my($prom) = $prom_obj->{promoter};			         # Promoter name.
  
  chdir $self->{parameters}{Base_dir};
 
  my($marker)=0;
  my(@factors); 
                                                                                       
  open(INFO, $final_db) || die "--> Cannot open map file: $final_db\n";  # Open final map for input.

  while (<INFO>) {		
    if (/^\w/ && /$prom$/) {  
      $marker=1;        			                 # If given promoter match.

      my(@m) = split(/\s/, $_);		
      push @factors, &DNA::Func::create_factor(\@m, $prom_obj);  # Store all TF object in an array without any filtering: default.
    }	  						         # Close "hit" case.
  }							         # Close promoter match case.
  close INFO;

  if ($marker!=1) { 
    print "\n--> Error: no TF hits for promoter - $prom";        # If no promoter match: quit. 
    #exit;
  }
  else { 
    #@{ $self->{final_list} } = $self->selection_filter(@factors);  # Intersection filter. 
    @{ $self->{final_list} } = $self->union_filter(@factors); 	   # Union filter. 
    $self->ordering($order);	                                   # Order TFs as selected.
  }
}

#---------------------------------- Gene_search sub. ----------------------------------

# 2. &gene_search subroutine to collect promoter list bound by a specified gene.

sub gene_search() {
  my($self) = shift;

  my($q_gene) = $self->{gene_q};                              # Query gene.	
  my($final_db) = $self->{parameters}{Finalized_db};	      # Finalized data base name.  
  
  chdir $self->{parameters}{Base_dir};
 
  my(@factors);
  my($marker)=0;                                                                                         
  my(%TF_hash);	 
                                                                                       
  open(INFO, $final_db) || die "--> Cannot open map file: $final_db\n";   # Open final map for input.

  while (<INFO>) {		       			
    if (/^\w/ && /$q_gene/) {  
      $marker=1;        			              # If given gene match.
      my(@m) = split(/\s/, $_);

      push @factors, &DNA::Func::create_factor(\@m); 	      # Store all TF object in an array without any filtering: default.
    }							      # Close "hit" case.
  }							      # Close promoter match case.
  close INFO;

  if ($marker!=1) { 
    print "\n--> Error: no TF hits for query gene - $q_gene\n\n";  # If no promoter match: quit. 
    exit;
  }
  else { 
    #foreach my $tf ($self->selection_filter(@factors)) {      # Intersection filter.
    foreach my $tf ($self->union_filter(@factors)) {	      # Union filter. 
      push @{ $TF_hash{$tf->{promoter}} }, $tf; 
    }

    $self->{final_list} = \%TF_hash; 	                      # IO database in standard format.			    
  }
}

#----------------------------------- Ordering sub. ------------------------------------

# 3. &ordering subroutine to order TF set in array.
                                                                                        
sub ordering($) {
  my($self, $order) = @_;				
  
  my(@result) = @{ $self->{final_list} };
  my(@pre_ordered) = sort { $a->{offset} <=> $b->{offset} } @result; 
  
  my(@ordered);							     # Final order of TFs.
  if (!defined $order) { @ordered = @pre_ordered; }                  # Ordering of TF's.
  elsif ($order==1) { @ordered = sort { $a->{factor} cmp 
                                        $b->{factor} } @pre_ordered; }
  else {
    print "\n--> Usage of order parameter: 1 - alphabetical ordering\n\n";
    exit;
  }
  @{ $self->{final_list} } = @ordered;
} 

#--------------------------------- Data_set_assign sub. -------------------------------

# 4. &data_set_assign subroutine to determine data set selection for visualization.

sub data_set_assign($$$) {
  my($self, $dana, $f, $s) = @_;
  
  if ($dana==1) { $self->eliminate(); }                      # Binary input handling.
  if ($dana==0) { 
    if ($f==0 && $s==0) { 
      print "\n--> Error: empty data set was requested - [0 0 0]!\n\n"; 
      exit; 
    }
    else { @{ $self->{final_list} } = @{ $self->{ernest_list} }; }
  }  

  @{ $self->{final_list} } = sort { $a->{offset} <=> $b->{offset} }
                                          @{ $self->{final_list} };
}                                                        

#------------------------------------- All_tf sub. -----------------------------------

# 5. &all_tf subroutine to construct TF database from Ernest's maps.
                                                                          
sub all_tf($) {
  my($self, $file) = @_;

  my($pssm_dir) = $self->{parameters}{PSSMS};

  my($base_dir) = $self->{parameters}{Base_dir};
  chdir $base_dir."ernest/standard_maps/";

  my(@e_factors) = ();					      # TF storge bin from Ernest's db.

  open(INFO, $file) || die "--> Cannot open Ernest's file: $file\n";  # Open the file for input.
  while (<INFO>) {        
    my(@m) = split(/\s/, $_); 

    chop $m[10];					      # Remove ';' character at end.
    my($proms) = $m[10];				      # Get promoter name list.

    my(@p) = split(/\|/, $proms);			      # Get list of promoter name(s).

    foreach my $p (@p) {
      my($e_factor) = new DNA::Fact($m[9]);		      # Create new TF object.

      $e_factor->{chromosome} = $m[0];		      	      # Assign chromosome number.
      $e_factor->{promoter} = $p;
 
      $e_factor->range_orient($m[3], $m[4], $m[6], $pssm_dir);

      push @e_factors, $e_factor;			      # Store Ernest's TF data.
    }
  }  
  close INFO;
  return @e_factors;					      # Return non-standard TF list.
} 

#-------------------------------- Motif_db_ernest sub. -----------------------------------

# 6. &motif_db_ernest subroutine to construct TF -> motif database from Ernest's maps.

sub motif_db_ernest($@) {
  my($self, $prom_list, @e_factors) = @_;

  my(%prom_list) = %{ $prom_list };			      # Get promoter library.
  my(%TF_hash);					              # TF -> motif database bin.

  foreach my $tf (@e_factors) {				      # Offset calculation.

    my($prom) = $tf->{promoter};		              # Get promoter name(s) string.
    my($factor) = $tf->{factor};			      # Get name of current TF.

    if ($prom =~ /DNE/) {
      $tf->{promoter} = "DNE"; 
      $tf->{offset} = "[-]"; 
      $tf->{score} = "[------------------]";

      my($start) = $tf->{motif_start};
      my($end) = $tf->{motif_end};

      $tf->{motif} = "[$start:$end]"; 

      push @{ $TF_hash{$factor}}, $tf;			      #	Create new member and assign it to storage.
    }
    else {
      my($prom_obj) = $prom_list{$prom};		      # Get promoter object from library.
   
      my($dir) = $self->{parameters}{PSSMS};
      my($mat) = $self->{parameters}{Matrix};

      $tf->get_offset_motif($prom_obj);
      $tf->get_p_value($mat, $dir);

      push @{ $TF_hash{$factor}}, $tf;		              #	Create new member and assign it to storage.
    } 
  }
  return \%TF_hash;		 		 	      # Return TF -> motif database.			
}

#---------------------------------- Promoter_db sub. -----------------------------------

# 7. &promoter_db subroutine to construct promoter database from standard maps.

sub promoter_db() {
  my($self) = shift;

  my($base_dir) = $self->{parameters}{Base_dir};             # Base directory.	

  chdir $base_dir."maps/";			             # Change to map direcoty.
  my(@files) = <*.m_db>;          		    	     # Collect all .m_db files in a list.

  foreach my $w (@files) {
    chdir $base_dir."maps/";  	  	                        # Change to map direcoty.
    open(INFO, $w) || die "--> Cannot open *.m_db file: $w\n";  # Open the file for input.

    my(%TF_hash);	   				        # TF hit (data) hash.

    while (<INFO>) {
      if(/\A>/) {
        my(@m) = split(/\s/, $_);
        my($orient) = chop $m[1]; 

        # Eliminate common elements in TF hash. 

        my($same)=0;

        if (defined $TF_hash{$m[6]}) {
          LP: foreach my $t (@{ $TF_hash{$m[6]} } ) {

            my($f) = $t->{factor};
            my($o) = $t->{offset};

            if ($f eq $m[1] && $o==$m[2]) { 
              $same=1; 
              last LP;
            }
          }
        }

        if ($same!=1) {
          my($tf) = new DNA::Fact($m[1]);		     # Create TF object. 

          $tf->{orient} = $orient;	            	     # Assign orientation,
          $tf->{offset} = $m[2];		             # motif offset, 
          $tf->{score} = $m[3];		                     # p-value (score),
          $tf->{motif} = $m[4];		  	             # and sequence.	

          $tf->{chromosome} = $m[5];	            	     # Get chromosome number.
	  $tf->{promoter} = $m[6]; 		    	     # Get promoter name.

          push @{ $TF_hash{$m[6]} }, $tf;                    # Record new TF occurence.
        }
      }
    }
    close INFO;

    $self->db_output(\%TF_hash, $w, 'p');                    # IO database in standard format.			    
  }
}

#----------------------------------- Final_map sub. -----------------------------------

# 8. &final_map subroutine to construct finalized database from standard maps.

sub final_map() {
  my($self, $tf_conds, $proms) = @_;

  my(%TF_hash);	       				                # TF hit (data) hash.

  my($base_dir) = $self->{parameters}{Base_dir};                # Base directory.	
  my($d_files) = $self->{parameters}{Data_bases}; 
    
  my(@files) = split(/\s/, $d_files);

  my($i)=0; 		                                        # File iteration count.
  foreach my $w (@files) {
    chdir $base_dir."maps/";  	  	                        # Change to map direcoty.
    open(INFO, $w) || die "--> Cannot open *.p_db file: $w\n";  # Open the file for input.

    my($str) = &DNA::Func::get_stringency($w);			# Get stringency value of data set.

print "\n\nSTRINGENCY: $str.\n\n";

    while (<INFO>) {
      if(/\A>/) {
        my(@m) = split(/\s/, $_);
        my($orient) = chop $m[1]; 

        # Eliminate common elements in TF hash. 

        my($same)=0;

        if (defined $TF_hash{$m[6]}) {
          LP: foreach my $t (@{ $TF_hash{$m[6]} } ) {

            my($f) = $t->{factor};
            my($o) = $t->{offset};

            if ($m[6] eq "DNE") {

              my($r) = $t->{motif};
              my($c) = $t->{chromosome};
 
              if (($f eq $m[1]) && ($r eq $m[4]) 
                                && ($c eq $m[5])) {
 
                $t->{stringency} = $str;                     # Remember data set order.		
                $same=1; 
                last LP;
              }
            }
            elsif ($f eq $m[1] && $o==$m[2]) { 

              $t->{strincency} = $str;            		
              $same=1; 
              last LP;
            }
          }
        }

        if ($same!=1) {
          my($tf) = new DNA::Fact($m[1]);		     # Create TF object. 

          $tf->{orient} = $orient;	            	     # Assign orientation,
          $tf->{offset} = $m[2];		             # motif offset, 
          $tf->{motif} = $m[4];		  	             # and sequence.	

          $tf->{chromosome} = $m[5];	            	     # Get chromosome number.
	  $tf->{promoter} = $m[6]; 		    	     # Get promoter name.
	  
  	  $tf->{stringency} = $str;	                     # Remember data set order.		

          $tf->{score} = $m[3];             
          #if ($str==0) { $tf->{score} = $m[3]; }             # Assign calculated p-value.
          #else { $tf->{score} = "[---------]"; }    	     # In not Nir's data base: -	

          my($cond, $min_p) = &DNA::Func::get_bind_affinity($m[1], 
                                        $m[6], $tf_conds, $proms);

	  $tf->{condition} = $cond; 		    	     # Assign condition of lowest p-value.
          $tf->{binding} = $min_p; 		             # Assign lowest binding affinity value.		

          push @{ $TF_hash{$m[6]} }, $tf;                    # Record new TF occurence.
        }
      }
    }
    close INFO;
    $i++;
  }
  $self->final_map_output(\%TF_hash); 	                     # IO database in finalized format.			    
}

#------------------------------ Count_TF_occurence sub. -------------------------------

# 9. &count_TF_occurence subroutine to count occurences of each TF in the promoter.

sub count_TF_occurence() {
  my($self) = shift;

  foreach my $tf (@{ $self->{final_list} }) {
    my($tf_name) = $tf->{factor};
    $self->{tf_occurence}->{$tf_name}++;
  }
}

#---------------------------------- O_TF_search sub. ----------------------------------

# 10. &o_TF_search subroutine to collect standarized, bound TF data given a 
# promoter of chosen organism.
                                                                                        
sub o_TF_search($) {
  my($self, $order) = @_;

  my($prom_obj) = $self->{prom_obj};			      # Promoter object.
  my($prom) = $prom_obj->{promoter};			      # Promoter name.
  
  my($dir) = $self->{parameters}{Other_dir};
  my($org) = $self->{prom_obj}->{organism};  		      # Get name of yeast species.

  my($m_file) = &DNA::Func::get_org_file($dir, $org);

  open(INFO, $m_file) || die "--> Cannot open map file: $m_file\n";   # Open the file for input.
   
  my(@factors);
  while (<INFO>) {
    my(@m) = split(/\s/, $_);	

    if (/^>/ && $m[6] =~ /$prom$/) {   			      # If given promoter match.
      my($orient) = chop $m[1];			              # Get orientation.
      my($factor) = new DNA::Fact($m[1]);		      # Create TF object. 

      $factor->{offset} = $m[2];		       	      # Parse out field values.
      $factor->{orient} = $orient;
      $factor->{score} = $m[3]; 
      $factor->{motif} = $m[4];

      $factor->{chromosome} = $m[5]; 
      $factor->{promoter} = $m[6];

      my($prom_start) = $prom_obj->{prom_start};  	      # Promoter start.
      my($strand) = $prom_obj->{strand};	 	      # Promoter strand.

      &DNA::Func::get_global_range($factor, $prom_start, $strand);

      push @factors, $factor;			     	      # Store TF object in an array.
    }						      	      # Close elimination loop.
  }							      # Close "hit" case.							      # Close promoter match case.
  close INFO;

  @{ $self->{final_list} } = @factors; 
  $self->ordering($order);	                              # Order TFs as selected.
}


#------------------------------ Selection_filter sub. -----------------------------------

# 11. &selection_filter to filter TFs by binding affinity p-value threshold and/or
# binding condition and/or conservation number and/or motif probability based on PSSM. 
                                                                                        
sub selection_filter(@) {
  my($self, @factors) = @_;

  my($B_THRESH) = $self->{parameters}{Binding_THRESH};		# P-value threshold of binding affinity.
  my($COND) = $self->{parameters}{Condition};			# Condition for lowest binding p-value.
  my($CONS) = $self->{parameters}{Conservation};		# Number of yeast species TF is conserved in.
  my($PROB) = $self->{parameters}{Probability};			# Probability to be a motif based on PSSM.

  if ($B_THRESH ne "N/A") { 					# If binding affinity filter selected.
    @factors = &DNA::Func::binding_filter(\@factors, $B_THRESH);
  }  

  if ($COND ne "N/A") { 					# If binding condition filter selected.
    $self->get_conditions();					# Get all name of all experimental conditions of database.

    if (!exists $self->{tf_conds}->{$COND}) { 
      print "\n--> Error: $COND - wrong condition in field -> [Condition]\n\n";
      exit;
    }
    @factors = &DNA::Func::condition_filter(\@factors, $COND);
  }

  if ($CONS ne "N/A") { 					# If evolutionary conservation filter selected.
    @factors = &DNA::Func::conservation_filter(\@factors, $CONS);
  } 

  if ($PROB ne "N/A") { 					# If motif probability filter selected.
    @factors = &DNA::Func::probability_filter(\@factors, $PROB);
  } 

  my($num) = scalar @factors;

  if ($num==0) { 
    print "\n--> Error: no such filtered TF hits for promoter!\n\n";  # If no TF match: quit. 
    exit;
  }

  return @factors;						# Return filtered list of TFs.
}

#--------------------------------- Correct_motif sub. ---------------------------------

# 12. &correct_motif subroutine to construct finalized database with motif shift correction.

sub correct_motif($) {
  my($self, $prom_list) = @_;

  my(%TF_hash);	       				                # TF hit (data) hash.
  my(%prom_list) = %{ $prom_list };			        # Get promoter library.

  my($dir) = $self->{parameters}{PSSMS};			# Get PSSM directory.
  my($mat) = $self->{parameters}{Matrix};			# Get matrix type: 0|1.

  my($base_dir) = $self->{parameters}{Base_dir};                # Base directory.	
  my($final_db) = $self->{parameters}{Finalized_db}; 

  chdir $base_dir;  	  	                                # Change to map direcoty.
  open(INFO, $final_db) || die "--> Cannot open *.db file: $final_db\n";  # Open the file for input.

  while (<INFO>) {
    my(@m) = split(/\s/, $_);

    if(/^\w/ && defined $m[8]) {
      my($orient) = chop $m[0]; 

      my($tf) = new DNA::Fact($m[0]);		    	        # Create TF object. 

      $tf->{orient} = $orient;	            	     		# Assign orientation,
      $tf->{offset} = $m[1];		             		# motif offset, 
      $tf->{score} = $m[2];                                     # score (p-value),
      $tf->{motif} = $m[6];		  	             	# and motif sequence.	

      $tf->{binding} = $m[3]; 		                        # Assign lowest binding affinity value,
      $tf->{condition} = $m[4]; 		    	        # condition of lowest p-value,
      $tf->{stringency} = $m[5];				# and evolutionary conservation number.

      $tf->{chromosome} = $m[7];	            	     	# Get chromosome number.
      $tf->{promoter} = $m[8]; 		    	     		# Get promoter name.

      if (!($m[5]==0 || $m[8] eq "DNE")) {			# Correct only for Ernest's not Nir's TFs.
        my($prom_obj) = $prom_list{$m[8]};			# Retrieve promoter of interest.
        $tf->shift_correct($prom_obj, $dir, $mat);		# Modify motif sequence and score if needed.
      }

      push @{ $TF_hash{$m[8]} }, $tf;                           # Record new TF occurence.
    }
  }
  close INFO;

  $self->final_map_output(\%TF_hash); 	                        # IO database in finalized format.			    			    
}

#-------------------------------- Union_filter sub. -----------------------------------

# 13. &union_filter to filter TFs by binding affinity p-value threshold and/or
# binding condition and/or conservation number and/or motif probability based on PSSM. 
                                                                                        
sub union_filter(@) {
  my($self, @factors) = @_;

  my($B_THRESH) = $self->{parameters}{Binding_THRESH};		# P-value threshold of binding affinity.
  my($COND) = $self->{parameters}{Condition};			# Condition for lowest binding p-value.
  my($CONS) = $self->{parameters}{Conservation};		# Number of yeast species TF is conserved in.
  my($PROB) = $self->{parameters}{Probability};			# Probability to be a motif based on PSSM.

  if ($B_THRESH eq "N/A") { $B_THRESH=-1; }			# Assign non-sense value, so condition never true.
  if ($CONS eq "N/A") { $CONS=5; }			
  if ($PROB eq "N/A") { $PROB=5; }			

  if ($COND ne "N/A") { 					# If binding condition filter selected.
    $self->get_conditions();					# Get all name of all experimental conditions of database.

    if (!exists $self->{tf_conds}->{$COND}) { 
      print "\n--> Error: $COND - wrong condition in field -> [Condition]\n\n";
      exit;
    }
  }

  my(@union)=();						# Union array container for TFs.

  for (my $i=0; $i<@factors; $i++) {

    my($tf) = $factors[$i]; 
    my($bind) = $tf->{binding};

    if (($bind <= $B_THRESH && !($bind =~/^\[/)) ||             # If binding p-value is under a given value.
	($tf->{condition} =~ /\b$COND\b/) 	 ||		# If binding condition of interest.
        ($tf->{conservation} >= $CONS)           ||             # If conserved in at least # many times in other yeast organisms.
        ($tf->{score} >= $PROB)) {

      push @union, $tf;			                        # Store TF in a list. 
    }
  }

  my($num) = scalar @union;

  if ($num==0) { 
    print "\n--> Error: no such 'union' TF hits for promoter!\n\n";  # If no TF match: quit. 
    #exit;
  }

  return @union;						# Return union of TFs.
}

#----------------------------------- Hg_pair_list sub. --------------------------------

# 14. &hg_pair_list subroutine to construct hyper-geometric pair library.
                                                                          
sub hg_pair_list() {
  my($self) = shift;

  my(%hg_pairs)={};					        # TF pair storge bin from hyper-geometric comparision data.

  my($base_dir) = $self->{parameters}{Script_outputs};          # Data directory.	
  my($hg_file) = $self->{parameters}{Global_hg_comparison}; 

  chdir $base_dir."hyper_geometric/";			        # Change directory where data is located.
  open(INFO, $hg_file) || die "--> Cannot open global hyper-geometric comparison: $hg_file\n";  # Open the file for input.

  while (<INFO>) {  
    if (/^\d*/) {    
      my(@m) = split(/\t/, $_); 

      #------- Pair member no. 1 -------

      my(@fact_p1) = split(/\s/, $m[1]);		        # Factor name and number of bound promoters (1st pair element).
      chop $fact_p1[0];					        # Remove ':' from end of string.

      my($TF_1) = new DNA::Fact($fact_p1[0]);		        # Create new TF object.
      $TF_1->{bound_promoters} = $fact_p1[1];		        # Assign number of no. of bound promoters.

      #------- Pair member no. 2 -------

      my(@fact_p2) = split(/\s/, $m[2]);		        # Factor name and number of bound promoters (2nd pair element).
      chop $fact_p2[0];

      my($TF_2) = new DNA::Fact($fact_p2[0]);		        # Create new TF object.
      $TF_2->{bound_promoters} = $fact_p2[1];		        # Assign number of no. of bound promoters.

      my($TF_pair) = new DNA::TF_pair(sort {$a->{factor} cmp 
			       $b->{factor} } ($TF_1, $TF_2));  # Create TF pair storing mambers by alphabetical order.

      chomp($TF_pair->{hg_comp_value} = $m[4]);	                # Hyper-geometric comparision value.

      my($pair_key) = $TF_pair->{position_1}->{factor}.'-'.     # Create TF pair key for storage hash.
		      $TF_pair->{position_2}->{factor};

      $hg_pairs{$pair_key} = $TF_pair;			        # Assign TF pair with 'hg' results to library.	
    }
  }  
  close INFO;
  return \%hg_pairs;					        # Return hyper-geometric TF pair list.
} 

#------------------------------------ Mesh_hg_ol sub. ---------------------------------

# 15. &mesh_hg_ol subroutine to construct overlap pair library with hyper-geometric data.

sub mesh_hg_ol($) {
  my($self, $hg_pairs) = @_;

  my(%hg_ol_pairs)={};                                          # Meshed overlap <-> hyper-geometric comparision data library.

  my($base_dir) = $self->{parameters}{Script_outputs};          # Data directory.
  my($pair_file) = $self->{parameters}{Overlapping_pairs};

  chdir $base_dir."overlap/";                                   # Change directory where data is located.
  open(INFO, $pair_file) || die "--> Cannot open overlapping pair: $pair_file\n";  # Open the file for input.

  while (<INFO>) {
    if (/^\>/) {
      my(@m) = split(/\t/, $_);

      if (exists $hg_pairs->{$m[1]}) {
         my($TF_pair) = $hg_pairs->{$m[1]};
         $TF_pair->{distinct_proms} = $m[2];                    # Assign no. of distinct promoters pair occurs in.

         $hg_ol_pairs{$m[1]} = $TF_pair;                        # Store Tf pair with needed data.
         delete($hg_pairs->{$m[1]});                            # Free up some memory.
      }
      else { push @{ $self->{no_hg_result} }, $m[1]; }            # If no 'hg' result for pair, remember it.
    }
  }
  close INFO;
  return \%hg_ol_pairs;                                         # Return overlap <-> hyper-geometric TF pair list.
}

#------------------------------------ Ol_tf_ratio sub. --------------------------------

# 16. &ol_tf_ratio subroutine to calculate overlap-to-pair occurence ratio.
                                                                          
sub ol_tf_ratio() {
  my($self) = shift;

  my(@ratio)={};					        # TF pair storge bin for occurence ratio.

  my($base_dir) = $self->{parameters}{Script_outputs};          # Data directory.	
  my($op_file) = $self->{parameters}{Overlapping_pairs}; 

  chdir $base_dir."overlap/";			                # Change directory where data is located.
  open(INFO, $op_file) || die "--> Cannot open overlap/pair: $op_file\n";  # Open the file for input.

  while (<INFO>) {  
    if (/^\>/) {    
      my(@m) = split(/\t/, $_); 

      #------- Ratio calculation -------

      my($TF_pair) = new DNA::TF_pair();

      $TF_pair->{pair_name} = $m[1];
      $TF_pair->{distinct_proms} = $m[2];			# No. of distinct promoters pair occurs in.
      $TF_pair->{overlap_proms} = $m[4];			# No. of distinct promoters overlap occurs in.
      $TF_pair->{overlap_num} = $m[5];				# No. of times overlap occurs.

      my($ratio) = ($m[4]/$m[2])*100;
      $TF_pair->{overlap_ratio} = $ratio;			# No. of times overlap occurs.

      push @ratio, $TF_pair;					# Store pair in a list.	
    }
  }  
  close INFO;

  return (sort { $b->{overlap_ratio} <=> $a->{overlap_ratio} } @ratio);
} 

#----------------------------------- Binding_ol sub. ----------------------------------

# 17. &binding_ol subroutine to retrieve binding strength with or w/o a specific con-
# dition of overlapping pairs.

sub binding_ol() {
  my($self, $tf_conds, $proms) = @_;

  my(%binding_ol)={};	       				        # Library of pairs by pair name.

  my($base_dir) = $self->{parameters}{Script_outputs};          # Data directory.
  my($ratio_file) = $self->{parameters}{Overlap_ratio};

  chdir $base_dir."overlap/";                                   # Change directory where data is located.
  open(INFO, $ratio_file) || die "--> Cannot open overlap-ratio file: $ratio_file\n";  # Open the file for input.

  while (<INFO>) {
    if (/^\>/) {
      my(@m) = split(/\t/, $_);

      #------- Binding strenght calculation -------

      my($TF_pair) = new DNA::TF_pair();

      $TF_pair->{pair_name} = $m[1];
      $TF_pair->{distinct_proms} = $m[2];			# No. of distinct promoters pair occurs in.
      $TF_pair->{overlap_proms} = $m[3];			# No. of distinct promoters overlap occurs in.
      $TF_pair->{overlap_num} = $m[4];				# No. of times overlap occurs.
      $TF_pair->{overlap_ratio} = $m[5];			# No. of times overlap occurs.

      my(@ratio) = split(/\=/, $m[5]); chop $ratio[1];
      $TF_pair->{ratio} = $ratio[1];

      my(@pair) = split(/\-/, $m[1]);

      my($i)=1;
      foreach my $tf (@pair) {

        #my($cond, $min_p) = &DNA::Func::get_bind_affinity($tf, # ***** GET P-VAL SOMEHOW? *****
        #                                   $tf_conds, $proms);
       
        my($pos) = "position_".$i;
		             	
        #$TF_pair->{$pos}->{condition} = $cond;    	        # Assign condition of lowest p-value. 
        #$TF_pair->{$pos}->{binding} = $min_p;                  # Assign lowest binding affinity value.
     
        $i++;							# Increment position counter.
      }
      
      #my($cond, $min1, $min2) = &DNA::Func::get_bind_cond(\@pair, 	# ***** GET P-VAL/COND SOMEHOW? *****
      #                                        $tf_conds, $proms);

      $binding_ol{$m[1]} = $TF_pair;				# Store pair in a list.	
    }
  }
  close INFO;

  $self->binding_ol_output(\%binding_ol); 	                # IO database in format of interest.			    
}

#----------------------------------- Kl_divergence sub. -------------------------------

# 18. &kl_divergence subroutine to calculate KL-divergence values for each TF (gene) pair.
                                                                          
sub kl_divergence($) {
  my($self, $TF_list) = @_;

  my(%PSSMS);					    # PSSM matrix hash for all known TFs.
  my(%kl_hash);					    # Hash for KL-divergence value storage by pair.

  chdir $self->{parameters}{PSSMS}."prob_matrices"; # Change to pssms directory.

  foreach my $motif (@$TF_list) {
    my($file) = $motif.".pssm";
		     
    if (-f $file) {
      open(INFO, $file) || die "--> Cannot open PSSMS file: $file\n";

      my(@pssm);
      while(<INFO>) { push(@pssm, $_); }
      close INFO;

      my($line_A) = $pssm[0]; my($line_C) = $pssm[1];
      my($line_T) = $pssm[2]; my($line_G) = $pssm[3]; 

      my(@l_A) = split(/\s/, $line_A);              # Put pssm values into an array, @l_base.
      my(@l_C) = split(/\s/, $line_C);
      my(@l_T) = split(/\s/, $line_T);
      my(@l_G) = split(/\s/, $line_G);

      my(%pssm) = ( A => \@l_A, C => \@l_C,         # Position weight matrix, a hash.
                    G => \@l_G, T => \@l_T );

      $PSSMS{$motif} = \%pssm;                      # Store PSSM in a storage hash.
    }
  }

  my($count)=0;
 
  for (my $i=0; $i<=$#$TF_list; $i++) {	            # Create all TF pairs in promoter.
    my($fac_A) = $TF_list->[$i];	    	    # Fixed TF in list, base of comparision. 

    for (my $j=$i+1; $j<=$#$TF_list; $j++) {        # If same TF, skip to next iteration.
      my($fac_B) = $TF_list->[$j];                  # TF of iterative comparision.

      my($fac_A, $fac_B) = sort {$a cmp $b} ($fac_A, $fac_B);	 # Order TF names alphabetically.
      my($pair) = $fac_A."-".$fac_B; 	            # Create TF pair key. 

      $kl_hash{$pair} = &DNA::Func::get_kl_score($fac_A, $fac_B, \%PSSMS);

      $count++; print "\n>$count PAIR: $pair is processed!";
    }
  }
  return \%kl_hash;				    # Return storage hash reference.
}
