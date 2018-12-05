#-------------------------------------------------------------------------
# Author: Mirko Palla.
# Date: July 12, 2006.
# For: Task 6 at the Church Lab - Genetics Department, 
# Harvard Medical School.
# 
# Purpose: This program contains the complete code for 
# class Promoter_set, modeling a set of promoters in Perl.
#
# Associated subroutines: 0. &constructor           1. &get_orthogroups
#			  2. &get_comp_table        3. &create_table
#			  4. &all_tf_list_promoter  5. &relative_distance   
#			  6. &distance_to_gene      7. &get_ortholog_data   
#			  8. &tf_pair_count	    9. &get_cer_orthologs
#			  10. &correct_OL_key
#------------------------------------------------------------------------- 

package DNA::Promoter_set;

use strict;
use base ("DNA::Promoter", "DNA::Fact_list", "DNA::TF_pair", "DNA::IO");

1;

#--------------------------------- Class Promoter_set -------------------------------------

# 0. constructor for promoter set object.

sub new {
  my($class, $gene, $prom) = @_;
  my($self)={};

  $self->{gene} = $gene;		 		        # Query gene field.
  $self->{promoter} = $prom; 				        # Equivalent promoter field.

  $self->{ortho_groups} = [];					# List of orthologous promoters.
  $self->{promoter_set} = {};                                   # Promoter set (by some rule).
  $self->{factors_list} = {};					# List of TFBS set of promoters.
  $self->{ortho_table} = [];					# TF occurence table of yeast species.

  $self->{tf_pair_dist} = {};					# List of relative TF distance.
  $self->{tf_pair_orient} = {};					# List of TF pair orientation type.
  $self->{tf_gene_dist} = {};					# List of TF-gene distance.

  $self->{tf_pair} = {};					# TF pair occurence hash.
  $self->{dist_prom} = {};					# Global TF pair occurence hash.

  $self->{tf_overlap} = {};					# Overlap [TF pair] occurence hash.
  $self->{OL_prom} = {};					# Global [TF pair] overlap occurence hash.

  $self->{filtered_pairs} = [];					# Filtered TF pair list by some threshold.

  $self->{parameters} = {};                                     # Input parameters.
  $self->{bar_code} = ();                                       # Bar code field to 
  			                                        # label output files.
  bless($self,$class);
  return($self);
}

#------------------------------ Get_orthogroups sub. ----------------------------------

# 1. &get_orthogroups subroutine to retrieve orthogroups of given gene entry of 'cer'.
                                                                          
sub get_orthogroups() {
  my($self) = shift;

  my($gene) = $self->{gene}; 	 		               # Query gene of interest.
  my($prom) = $self->{promoter}; 		               # Query promoter of interest.

  my($base_dir) = $self->{parameters}{Base_dir};	       # Get base directory.
  my($comp_file) = $self->{parameters}{Filtered_orthogroups};  # Get orthogroup file name. 

  chdir $base_dir;
  open(INFO, $comp_file) || die "--> Cannot open orthogroups file -> $comp_file\n";          
 
  my($num);
  my($line_num);

  my($marker)=0;

  while (<INFO>) {
    chomp(my(@ortho) = split(/\t/, $_));       # Put line member strings into an array, @ortho.

    if ($ortho[1] eq $gene || $ortho[1] eq $prom) { 
      $num = scalar @ortho;

      @{ $self->{ortho_groups} } = @ortho; 
      $line_num = $ortho[0];

      $marker=1; 
    }   
  }
  close INFO;

  $self->get_cer_orthologs($line_num);

  if ($marker!=1) { return 1; }		       # Promoter not present in orthogroup file.
  elsif ($num<4) { return 2; }	      	       # In orthogroup file, but no other orthologous promoters in any yeast organism. 	  
  else { return 0; }
}

#------------------------------- Get_comp_table sub. ----------------------------------

# 2. &get_comp_table subroutine to construct table of comparative analysis for ortho-
# groups of given 'cer' gene entry.
                                                                          
sub get_comp_table() {
  my($self) = shift;

  my(%valid_data) = %{ $self->{parameters} };            # Get parameter hash.
  my(@ortho_groups) = @{ $self->{ortho_groups} };

  my($max_iter) = scalar @ortho_groups;

  for (my $index=2; $index<$max_iter; $index++) {
    my($org_prom) = $ortho_groups[$index]; 

    my(@o_p_s) = split(/\,/, $org_prom);
 
    foreach my $o_p (@o_p_s) {
      my(@org_prom) = split(/\|/, $o_p);
      my($org) = $org_prom[0];				 # Parse out organism name.

      if ($org =~ /gos|wal|pom|gra|cra|gri|nid/) { next; }

      my($gene, $prom);
      my($prom_obj, $factors);

      #-------- Case I: S. cerevisiae --------

      if ($org eq "cer") {
        $prom = $org_prom[1];

        my($base_dir) = $valid_data{Base_dir};	 	 # Get base directory.
        my($c_file) = $valid_data{Conv_file};            # Get conversion file.
        $gene = &DNA::Func::promoter_to_gene($prom, 
		               $c_file, $base_dir);

        $prom_obj = new DNA::Promoter($gene, $prom);     # Create promoter object.
        %{ $prom_obj->{parameters} } = %valid_data;	 # Assign input data to prom_obj.

        $prom_obj->{organism} = $org;
        $prom_obj->promoter_search();		         # Retrieve promoter sequence.

        $factors = new DNA::Fact_list($prom_obj);	 # Create Fact_list object.
        %{ $factors->{parameters} } = %valid_data;	 # Assign input data to factors.

        if ( !(defined $prom_obj->{no_eq}) ) {
          $factors->TF_search(1);			 # Retrieve all TFs bound to given promoter.
        }		                 
      }

      #------ Case II: Other yeast species ------

      else { 
        my($file);				         # Promoter file of yeast species.

        if ($org =~ /alb|cas|par|mik|bay/) { 
          $gene = "NDEF"; 

          chomp(my(@prom) = split(/\-/, $org_prom[1]));

          if (defined $prom[2]) { $prom = $prom[2]; }
          elsif (defined $prom[1]) { $prom = $prom[1]; }
          else { $prom = $prom[0]; };    		 # Parse out promoter name.

          if ($org =~ /cas/) { 
            my($end) = substr($prom, -1);
            if ($end =~ /\D/) { chop $prom; } 
          }

          if ($org =~ /par|mik|bay/) { $prom = "ORFN:".$prom; }
        }
        elsif ($org =~ /gla|han|lac|lip/) { 
          $gene = $org_prom[1];

          my($dir) = $valid_data{Aliases};     		 # Pull out converstion file data.
          $file = &DNA::Func::get_org_file($dir, $org);

          open(ALIAS, $file) || die "--> Cannot open file list: $file!\n"; 

          my($marker);
          A: while (<ALIAS>) {			
            my(@line) = split(/\t/, $_);		 # Parse out corresponding promoter.

            chomp $line[1];
            if($gene eq $line[1]) { $prom = $line[0]; 
              $marker=1; last A; 
            } 
          }  
          close ALIAS;     

          if ($marker!=1) { $prom = $gene; }    
        }	

        $prom_obj = new DNA::Promoter($gene, $prom);     # Create promoter object.
        %{ $prom_obj->{parameters} } = %valid_data;	 # Assign input data to prom_obj.

        $prom_obj->{organism} = $org;
        $prom_obj->o_promoter_search();		         # Retrieve promoter sequence.

        $factors = new DNA::Fact_list($prom_obj);	 # Create Fact_list object.
        %{ $factors->{parameters} } = %valid_data;	 # Assign input data to factors.
        
        if ( !(defined $prom_obj->{no_eq}) ) {
          $factors->o_TF_search(1);
        }     
      }	

      if ( !(defined $prom_obj->{no_eq}) ) {
        $factors->count_TF_occurence();			 # Count TF occurence for each list.
      }

      push @{ $self->{factors_list}{$org} }, $factors;
    }
  }        
}

#--------------------------------- Create_table sub. ----------------------------------

# 3. &create_table subroutine to record TFBS occurences in a specified S. cerevisiae 
# promoter and in orthologous promoters of other yeast species sorted by chronogolical 
# order according to the phylogenetic tree derived.

sub create_table() {
  my($self) = shift;

  my(@phylo_tree) = ('cer', 'par', 'mik', 'bay', 'gla', 	    # Pholygenic tree of different yeast species
		     'cas', 'lac', 'han', 'alb', 'lip');	    # in chronological order.

  my($cer_list) = $self->{factors_list}{'cer'}[0];
  my(%cer_occu) = %{ $cer_list->{tf_occurence} };
  my(@cer_TFs) = sort keys %cer_occu;				    # List of 'cer' TFs of requested gene (promoter).

  my(@table);
  my($org_header);

  my(@curr_fact);

  foreach my $org (@phylo_tree) {				    # Create header (list of yeast organisms) as present in orthology group.			
    if (exists $self->{factors_list}{$org}) {
      foreach my $factors ( @{ $self->{factors_list}{$org} } ) {
        if ( !(defined $factors->{prom_obj}->{no_eq}) ) {

          my($prom) = $factors->{prom_obj}->{promoter};
          my($org_prom) = $org.'|'.$prom;
          
          $org_header = $org_header.$org_prom."\t";

          push @curr_fact, $factors;				    # Keep track of TF list (of promoter) of orthogroup.
        }
        else { print "\n--> No data for org: $org <=> prom: $factors->{prom_obj}->{no_eq}"; }
      }
    }
  }
  chop $org_header;
  push @table, "\n[TFs]\t".$org_header;			    	    # Assign 'header' to the first row of comparison table.

  foreach my $tf (@cer_TFs) {					    # Iterate through only the existing organisms.

    my($o_line) = "$tf\t\t";
    foreach my $factors (@curr_fact) {  
  
      my($o);							    # Create list of TF occurence for each promoter. 
      if (defined $factors->{tf_occurence}->{$tf}) {
        $o = $factors->{tf_occurence}->{$tf};
      }
      else { $o = '-'; }
         
      $o_line = $o_line.$o."\t\t";
    }
    chop $o_line;
    push @table, $o_line;					    # If row completed store it in table array lineraly.	
  }
  @{ $self->{ortho_table} } = @table;			            # Assign TF occurence table to storage.
}

#-------------------------------- All_tf_list_promoter sub. ---------------------------

# 4. &all_tf_list_promoter subroutine to collect all promoter objects with bound TF list
# in S. cerevisiae into a library.
                                                                                        
sub all_tf_list_promoter() {
  my($self) = shift;
  
  my($dir) = $self->{parameters}{Base_dir};
  my(%valid_data) = %{ $self->{parameters} };

  my($c_file) = $self->{parameters}{Conv_file};
  my($prom_file) = $self->{parameters}{Prom_file};

  chdir $dir;

  my($marker);
  open(INFO, $prom_file) || die "--> Cannot open promoter file -> $prom_file\n";          
 
  my($end);
  while (<INFO>) {

      #----- Promoter object retrieval -----

      my($prom_obj) = new DNA::Promoter();       # Create Promoter object.

      my(@prom) = split(/\s/, $_);               # Put line member strings into an array, @prom.
      chomp($prom_obj->{prom_seq} = $prom[4]);   # Parse out promoter sequence and assign to Promoter object.

      $marker=1;			         # Mark existance of such Promoter.

      my(@name) = split(/\s/, $prom[0]);
      my($prom_name) = substr($name[0],1);

      $prom_obj->{promoter} = $prom_name;
      $prom_obj->{gene} = &DNA::Func::promoter_to_gene($prom_name, $c_file, $dir);

      $end = substr($prom_name,6,1);             # Get promoter end: W|C.

      if ($end eq 'W' || $end eq 'C') {
        $prom_obj->{strand} = $end;		 # Assign strand value to promoter.
      }
      else { $prom_obj->{strand} = 'M'; }	 # Mitochondrial chromosome.
      
      my(@range) = split(/\D/, $prom[3]);
      
      chop $prom[2];
      $prom_obj->{chromosome} = $prom[2];        # Assign chromosome no.     

      $prom_obj->{prom_start} = $range[1];       # Assign global range.
      $prom_obj->{prom_end} = $range[2];

      #------ Bound TF list retrieval ------

      my($factors) = new DNA::Fact_list($prom_obj);  # Create Fact_list object.
      %{ $factors->{parameters} } = %valid_data;     # Assign input data to factors.

      $factors->TF_search();		 	     # Retrieve all TFs bound to given promoter.

print "\n> Promoter [$prom_name] processed!";

      $self->{factors_list}{$prom_name} = $factors;  # Assign TF list (with promoter object) to library. 
  }
  close INFO;
}

#-------------------------------- Relative_distance sub. ------------------------------

# 5. &relative_distance subroutine to get statistics on relative distance and orien-
# tation of TF pairs in all promoters in S. cerevisiae.

sub relative_distance() {
  my($self) = shift;

  foreach my $prom_name (keys %{ $self->{factors_list} }) {

    my($fact_list) = $self->{factors_list}{$prom_name};		# Retrieve TF list of each promoter.
    my($prom_obj) = $fact_list->{prom_obj};			# Retrieve promoter object.

    my($strand) = $prom_obj->{strand};				# Get strand type,
    my($chrom) = $prom_obj->{chromosome};		        # Chromosome number of current promoter.

    my(@ordered_list) = @{ $fact_list->{final_list} };		# Get TF list ordered by increasing offset.

    for (my $i=0; $i<@ordered_list; $i++) {			# Create all TF pairs in promoter.
      my($base_TF) = $ordered_list[$i];				# Fixed TF in list, base of comparision.

      for (my $j=$i+1; $j<@ordered_list; $j++) {
						# If same TF, skip to next iteration.
        my($comp_TF) = $ordered_list[$j];			# TF of iterative comparision.

        my($TF_pair);						# Declare "empty" TF pair.
          
        if ($strand eq 'W' || $strand eq 'M') {
          $TF_pair = new DNA::TF_pair($base_TF, $comp_TF);    # Assign order depending on strand.
        }
        elsif ($strand eq 'C') {
          $TF_pair = new DNA::TF_pair($comp_TF, $base_TF);
        }
        else { die "\n--> Error: no such strand - $strand\n\n"; }

        $TF_pair->{strand} = $strand;				# Label TF pair strand.
        $TF_pair->{promoter} = $prom_name;			# Label TF pair with promoter name.
        $TF_pair->{chromosome} = $chrom;			# Label TF pair with chromosome number.

        $TF_pair->orient_type();				# Assign orientation of TF pair.
        $TF_pair->calculate_dist();				# Calculate relative distance of TFs.

        my($dis) = $TF_pair->{distance};			# Get current distance value.
        my($ori) = $TF_pair->{orient_type};			# Get current orientation value.

        DESTROY($TF_pair);					# Delete TF pair object to save up memory.
 
        push @{ $self->{tf_pair_dist}{$prom_name} }, $dis;	# Assign parameters to storage array.
        push @{ $self->{tf_pair_orient}{$prom_name} }, $ori; 
      }  
    }
    delete($self->{factors_list}{$prom_name});		        # Delete promoter object and assiciated. 
  }
}

#-------------------------------- Distance_to_gene sub. ------------------------------

# 6. &distance_to_gene subroutine to get statistics on distance to gene start of TFs
# in all promoters in S. cerevisiae.

sub distance_to_gene() {
  my($self) = shift;

  foreach my $prom_name (keys %{ $self->{factors_list} }) {

    my($fact_list) = $self->{factors_list}{$prom_name};		# Retrieve TF list of each promoter.
    my($prom_obj) = $fact_list->{prom_obj};			# Retrieve promoter object.

    my($strand) = $prom_obj->{strand};				# Get strand type,
    my($p_length) = length $prom_obj->{prom_seq};      		# promoter sequence length.

    my(@ordered_list) = @{ $fact_list->{final_list} };		# Get TF list ordered by increasing offset.

    for (my $i=0; $i<@ordered_list; $i++) {			# Create all TF pairs in promoter.
      my($base_TF) = $ordered_list[$i];				# Fixed TF in list, base of comparision.

      my($distance);						# (Shortest) distance from TF to gene start.
      my($offset) = $base_TF->{offset};			        # Relative offset of TF (from left of sequence string).

      if ($strand eq 'W' || $strand eq 'M') {
        my($tf_length) = length $base_TF->{motif}; 		# Get length of motif.
	
        $distance = $p_length - ($offset + $tf_length);		# Calculate TF-to-gene distance for 'W'.
      }
      else { $distance = $offset; }				# Calculate TF-to-gene distance for 'C'.

      push @{ $self->{tf_gene_dist}{$prom_name} }, $distance;	
    }
  }
}

#----------------------------- Get_ortholog_data sub. ---------------------------------

# 7. &get_ortholog_data subroutine to collect data for all orthogroups of given 'cer' 
# gene entry inclusive.
                                                                          
sub get_ortholog_data() {
  my($self) = shift;

  my(%valid_data) = %{ $self->{parameters} };             # Get parameter hash.
  my(@ortho_groups) = @{ $self->{ortho_groups} };

  my($max_iter) = scalar @ortho_groups;

  for (my $index=2; $index<$max_iter; $index++) {
    my($org_prom) = $ortho_groups[$index]; 

    my(@o_p_s) = split(/\,/, $org_prom);
    
    foreach my $o_p (@o_p_s) {
      my(@org_prom) = split(/\|/, $o_p);
      my($org) = $org_prom[0];				 # Parse out organism name.

      if ($org =~ /gos|wal|pom|par|mik|bay/) { next; }

      my($gene, $prom);
      my($prom_obj, $factors);

      #-------- Case I: S. cerevisiae --------

      if ($org eq "cer") {
        $prom = $org_prom[1];

        my($base_dir) = $valid_data{Base_dir};	 	 # Get base directory.
        my($c_file) = $valid_data{Conv_file};            # Get conversion file.
        $gene = &DNA::Func::promoter_to_gene($prom, 
		               $c_file, $base_dir);

        $prom_obj = new DNA::Promoter($gene, $prom);     # Create promoter object.
        %{ $prom_obj->{parameters} } = %valid_data;	 # Assign input data to prom_obj.

        $prom_obj->{organism} = $org;
        $prom_obj->promoter_search();		         # Retrieve promoter sequence.

        $factors = new DNA::Fact_list($prom_obj);	 # Create Fact_list object.
        %{ $factors->{parameters} } = %valid_data;	 # Assign input data to factors.

        if ( !(defined $prom_obj->{no_eq}) ) {
          $factors->TF_search();		         # Retrieve all TFs bound to given promoter.
        }
      }

      #------ Case II: Other yeast species ------

      else { 
        my($file);				         # Promoter file of yeast species.

        if ($org =~ /alb|cas|gos|wal|pom|par|mik|bay/) { 
          $gene = "NDEF"; 

          my(@prom) = split(/\-/, $org_prom[1]);
          if (defined $prom[1]) { $prom = $prom[1]; }
          else { $prom = $prom[0]; };    		 # Parse out promoter name.

          if ($org =~ /cas/) { 
            my($end) = substr($prom, -1);
            if ($end =~ /\D/) { chop $prom; } 
          }
        }
        elsif ($org =~ /gla|han|lac|lip/) { 
          $gene = $org_prom[1];

          my($dir) = $valid_data{Aliases};     		 # Pull out converstion file data.
          $file = &DNA::Func::get_org_file($dir, $org);

          open(ALIAS, $file) || die "--> Cannot open file list: $file!\n"; 

          my($marker);
          A: while (<ALIAS>) {			
            my(@line) = split(/\t/, $_);		 # Parse out corresponding promoter.

            chomp $line[1];
            if($gene eq $line[1]) { $prom = $line[0]; 
              $marker=1; last A; 
            } 
          }  
          close ALIAS;     

          if ($marker!=1) { $prom = $gene; }    
        }	

        $prom_obj = new DNA::Promoter($gene, $prom);     # Create promoter object.
        %{ $prom_obj->{parameters} } = %valid_data;	 # Assign input data to prom_obj.

        $prom_obj->{organism} = $org;
        $prom_obj->o_promoter_search();		         # Retrieve promoter sequence.

        $factors = new DNA::Fact_list($prom_obj);	 # Create Fact_list object.
        %{ $factors->{parameters} } = %valid_data;	 # Assign input data to factors.
        
        if ( !(defined $prom_obj->{no_eq}) ) {
          $factors->o_TF_search();
        }     
      }	

      if ( !(defined $prom_obj->{no_eq}) ) {
        $factors->count_TF_occurence();			 # Count TF occurence for each list.
      }

      push @{ $self->{factors_list}{$org} }, $factors;
    }
  }        
}

#--------------------------------- Tf_pair_count sub. -------------------------------

# 8. &tf_pair_count subroutine to get statistics on occurence number of TF pairs in 
# all promoters in S. cerevisiae.

sub tf_pair_count() {
  my($self) = shift;

  my(%valid_data) = %{ $self->{parameters} };                   # Get parameter hash.
  my($dis) = $self->{parameters}{Overlap_gap};

  foreach my $prom_name (keys %{ $self->{factors_list} }) {

    my($fact_list) = $self->{factors_list}{$prom_name};		# Retrieve TF list of each promoter.
    my(@tf_list) = @{ $fact_list->{final_list} };		# Get TF list ordered by increasing offset.

    for (my $i=0; $i<@tf_list; $i++) {			        # Iterate through all TF pairs in promoter.
      my($TF_A) = $tf_list[$i];				        # Fixed TF in list, base of comparision.
      my($fac_A) = $TF_A->{factor};				# Get name of TF.  

      my($m_A_start) = $TF_A->{motif_start};		        # Get fixed TF motif start.
      my($m_A_end) = $TF_A->{motif_end};                        # Get fixed TF motif end.

      for (my $j=$i+1; $j<@tf_list; $j++) {
        my($TF_B) = $tf_list[$j];		         	# TF of iterative comparision.
        my($fac_B) = $TF_B->{factor};				# Get name of TF.  

        #--------------------- Determine overlapping ----------------------

        my($m_B_start) = $TF_B->{motif_start};		         # Get TF motif start to compare.
        my($m_B_end) = $TF_B->{motif_end};                       # Get TF motif end to compare.

        my($strand) = substr($prom_name,6,1);                    # Get promoter end: W|C.;

        my($range);
        if ($strand eq 'C') { $range = $m_A_end - $dis; }
        else { $range = $m_A_end + $dis; }

        my($overlap);						 # Overlap signal: {0,1}.
        if ((($strand eq 'W' || $strand eq 'M') && $m_B_start <= $range) || 
             ($strand eq 'C' && $m_B_start >= $range)) { 

          $overlap=1;  						 # Overlap, assign 1.
        }
        else { $overlap=0; }					 # No overlap, assign 0.

        #--------------- Determine global/overlap occurence ---------------

        my($pair1) = $fac_A."-".$fac_B;     
        my($pair2) = $fac_B."-".$fac_A;  

        if ($pair1 ne $pair2) {					# Eliminate same TF-TF pairs.       
          if (exists $self->{tf_pair}{$pair1}) {	        # If pair in promoter already exists.
            $self->{tf_pair}{$pair1}++;

            if (!defined $self->{dist_prom}{$pair1}{$prom_name}) {
              $self->{dist_prom}{$pair1}{$prom_name}=1;         # Mark distinct promoter. 
            }
          }
          elsif (exists $self->{tf_pair}{$pair2}) {	        # Increment occurence number.          
            $self->{tf_pair}{$pair2}++;

            if (!defined $self->{dist_prom}{$pair2}{$prom_name}) {
              $self->{dist_prom}{$pair2}{$prom_name}=1;         # Mark distinct promoter. 
            }           
          }
          else { 
            $self->{tf_pair}{$pair1}=1; 			# Initialize TF pair.
            $self->{dist_prom}{$pair1}{$prom_name}=1;           # Mark distinct promoter. 
          }

          if ($overlap==1) {					# If overlapping TF pair.
            if (exists $self->{tf_overlap}{$pair1}) {	        # If pair in promoter already exists.
              $self->{tf_overlap}{$pair1}++;

              if (!defined $self->{OL_prom}{$pair1}{$prom_name}) {
                $self->{OL_prom}{$pair1}{$prom_name}=1;         # Mark distinct promoter. 
              }
            }
            elsif (exists $self->{tf_overlap}{$pair2}) {	# Increment occurence number.          
              $self->{tf_overlap}{$pair2}++;

              if (!defined $self->{OL_prom}{$pair2}{$prom_name}) {
                $self->{OL_prom}{$pair2}{$prom_name}=1;         # Mark distinct promoter. 
              }           
            }
            else { 
              $self->{tf_overlap}{$pair1}=1; 			# Initialize TF pair.
              $self->{OL_prom}{$pair1}{$prom_name}=1;           # Mark distinct promoter. 
            }
          }
        }  
      }
    }
  }
  $self->correct_OL_key;					# Correct for TF pair library key if reversed.

  my($diff) = $self->{parameters}{Diff_prom_occurence};         # Number of promoters needed for filtering.

  my(%dist_prom) = %{ $self->{dist_prom} };			# Distinct promoter count hash.
  my(%tf_pairs) = %{ $self->{tf_pair} };			# TF pair occurence count hash.

  my($i)=0;

  #foreach my $pair (sort { keys %{ $dist_prom{$b} } <=> 
  #			   keys %{ $dist_prom{$a} } } keys %dist_prom) { 

  foreach my $pair (sort { $tf_pairs{$b} <=> $tf_pairs{$a} } keys %tf_pairs) { 
 
    my($o) = scalar keys %{ $dist_prom{$pair} };                  # No. of different promoters containing TF pair.
    if ($o>=$diff) { push @{ $self->{filtered_pairs} }, $pair; }  # Filter by some occurence threshold.
  }  
}

#----------------------------- Get_cer_orthologs sub. ---------------------------------

# 9. &get_cer_orthologs subroutine to retrieve 'cer' orthologs of given gene.
                                                                          
sub get_cer_orthologs($) {
  my($self, $line_num) = @_;

  my($gene) = $self->{gene}; 	 		               # Query gene of interest.
  my($prom) = $self->{promoter}; 		               # Query promoter of interest.

  my($base_dir) = $self->{parameters}{Base_dir};	       # Get base directory.
  my($comp_file) = $self->{parameters}{Filtered_orthogroups};  # Get orthogroup file name. 

  chdir $base_dir;
  open(INFO, $comp_file) || die "--> Cannot open orthogroups file -> $comp_file\n";          
					       
  while (<INFO>) {
    chomp(my(@ortho) = split(/\t/, $_));       # Put line member strings into an array, @ortho.

    if ($ortho[0] eq $line_num && 
      !($ortho[1] eq $gene || $ortho[1] eq $prom)) { 

print "\n--> ORTHOLOG: $ortho[2]";

      push @{ $self->{ortho_groups} }, $ortho[2];
    }   
  }
  close INFO;
}

#------------------------------- Correct_OL_key sub. ----------------------------------

# 10. &correct_OL_key subroutine to correct for TF pair library key if reversed.
                                                                          
sub correct_OL_key($) {
  my($self) = shift;

  foreach my $o_p (keys %{ $self->{tf_overlap} }) {
    if (!exists $self->{tf_pair}{$o_p}) {
      my(@spl) = split(/\-/, $o_p);

      my($rev) = $spl[1].'-'.$spl[0];

      $self->{tf_overlap}{$rev} = $self->{tf_overlap}{$o_p};
      $self->{OL_prom}{$rev} = $self->{OL_prom}{$o_p};

      delete($self->{tf_overlap}{$o_p});
      delete($self->{OL_prom}{$o_p});
    }
  }
}
