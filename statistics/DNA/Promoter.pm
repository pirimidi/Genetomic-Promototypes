#-------------------------------------------------------------------
# Author: Mirko Palla.
# Date: February 3, 2006.
# For: Task 3 at the Church Lab - Genetics Department, 
# Harvard Medical School.
# 
# Purpose: This program contains the complete code for 
# class Promoter, modeling an object for promoter data 
# storage in Perl.
#
# Associated subroutines: 0. &constructor      1. &promoter_search
#			  2. &gene_search      3. &promoter_list
#		          4. &promoter_count   5. &binding_summary  
#			  6. &ortho_groups     7. &o_promoter_search
#			  8. &o_promoter_search
#-------------------------------------------------------------------- 

package DNA::Promoter;

use strict;
use base ("DNA::Func", "DNA::IO", "DNA::Fact"); 

1;

#---------------------------------- Class Promoter ------------------------------------

# 0. constructor for promoter object.

sub new {
  my($class, $gene, $prom) = @_;
  my($self)={};

  $self->{gene} = $gene;		 		# Gene equivalent field.
  $self->{promoter} = $prom; 				# Promoter field.

  $self->{organism}; 				        # Yeast organism field.

  $self->{prom_seq} = ();	 			# Promoter sequence field.
  $self->{prom_start};					# Global promoter start.
  $self->{prom_end};					# Global promoter end.

  $self->{strand};					# Strand type: W|C.
  $self->{chromosome};					# Chromosome number.

  $self->{gene_q} = ();		 		        # Query gene field.
  $self->{gene_loc} = ();		 		# Query gene location.
  $self->{gene_ori} = ();		 		# Query gene orientation.

  $self->{no_eq};		 			# No promoter match field.
  $self->{out_r};		 			# Motif out of promoter range.

  $self->{tf_cond_list} = ();	 			# List of TF-condition pairs.
  $self->{ortho_list} = [];	 			# List of orthogroups of cerevisiae.

  $self->{parameters} = {};                             # Input parameters.

  bless($self,$class);
  return($self);
}

#--------------------------------- Promoter_search sub. -------------------------------

# 1. &promoter_search subroutine to get promoter sequence.
                                                                                        
sub promoter_search() {
  my($self) = shift;

  chdir $self->{parameters}{Base_dir};

  my($prom_file) = $self->{parameters}{Prom_file};
  my($prom) = $self->{promoter};

  my($marker);
  open(INFO, $prom_file) || die "--> Cannot open promoter file - $prom_file\n";          
 
  my($end);
  PROM: while (<INFO>) {
    if (/$prom/) {                          # If given promoter matches on any line.
      my(@prom) = split(/\s/, $_);          # Put line member strings into an array, @prom.
      chomp($self->{prom_seq} = $prom[4]);  # Parse out promoter sequence and assign to Promoter object.

      $marker=1;			    # Mark existance of such Promoter.

      chop $prom[2];
      $self->{chromosome} = $prom[2];       # Assign chromosome no.     

      $end = substr($prom,6,1);	            # Get promoter end: W|C.

      if ($end eq 'W' || $end eq 'C') {
        $self->{strand} = $end;		    # Assign strand value to promoter.
      }
      else { $self->{strand} = 'M'; }	    # Mitochondrial chromosome.
      
      my(@range) = split(/\D/, $prom[3]);

      $self->{prom_start} = $range[1];      # Assign global range.
      $self->{prom_end} = $range[2];

      last PROM;
    }
  }
  close INFO;
  if(!($marker==1)) { die "\n-->Error: no such promoter - $prom\n\n"; }
}

#---------------------------------- Gene_search sub. ----------------------------------

# 2. &gene_search subroutine to collect promoter list bound by a specified gene.
                                                                                        
sub gene_search($) {
  my($self, $file) = @_;
  my($query_gene) = $self->{gene_q};			# Assign query gene (TF) name.

  my($c_file) = $self->{parameters}{Conv_file};
  my($base_dir) = $self->{parameters}{Base_dir};

  my($marker);
  my(@prom_list) = ();				        # Promoter storage array.

  chdir $base_dir."ernest/standard_maps/";

  open(INFO, $file) || die "--> Cannot open Ernest's file: $file!\n";  
  while (<INFO>) {        
    if (/$query_gene/) {				# If given gene matches on any line.
      my(@line) = split(/\s/, $_);			# Parse out promoter candidate list.

      my(@proms) = split(/\|/, $line[10]);

      foreach my $o (@proms) {				# For all promoters in list.
      	my($check) = length $o;				# Check length of candidates.
        if ($check >= 5) {				
          my($prom) = new DNA::Promoter();		# Create Promoter object.

          my($end) = substr($o, -1);
          if ($end eq ";") { chop $o; }
                                                                                        
          my($gene) = &DNA::Func::promoter_to_gene($o, $c_file, $base_dir);
 
          if ($gene==1) { 
            push @{$self->{no_eq}}, $prom;              # Store no match promoter data.
            $prom->{gene} = '-';
          }  
          else { $prom->{gene} = $gene; }		# Assign parameter values to promoter object.

          $prom->{promoter} = $o;
          $prom->{chromosome} = $line[0]; 		# Assign chromosome no.     

	  $prom->{gene_ori} = $line[6];			# Assign query gene orientation.
	  $prom->{gene_loc} = $line[3].'-'.$line[4];	# Assign query gene location.
          
          push @prom_list, $prom;    		        # Store Ernest's promoter data.
          
	  $marker=1;  
  	}
      }
    }
  }
  close INFO;
  
  if(!($marker==1)) { 					      # If no gene match: quit.
    print "\n--> Error: no such gene - $query_gene\n\n"; 
    exit;
  }
  else { return @prom_list; }				# Otherwise return list of promoters.
}

#---------------------------------- Promoter_list sub. -------------------------------

# 3. &promoter_list subroutine to get promoter sequence list in S. cerevisiae - [New_promoters].
                                                                                        
sub promoter_list($) {
  my($self, $o) = @_;
  
  my($dir) = $self->{parameters}{Base_dir};
  my($c_file) = $self->{parameters}{Conv_file};
  my($prom_file) = $self->{parameters}{Prom_file};

  my(%promoter_list);			         # Promoter list holder hash.

  chdir $dir;

  my($marker);
  open(INFO, $prom_file) || die "--> Cannot open promoter file -> $prom_file\n";          
 
  my($end);
  while (<INFO>) {
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

      if (!(defined $o)) {
        $promoter_list{$prom_name} = $prom_obj;	 # Assign promoter object to library. 
      }
      elsif ($o==0) {
        push @{ $promoter_list{$prom[2]} }, $prom_obj;
      } 
      else {
        print "\n--> Error: argument options -> -|0\n";
        exit;
      }
  }
  close INFO;
  return \%promoter_list;			 # Return promoter list hash reference.
}

#----------------------------------- Promoter_count sub. ------------------------------

# 4a. &promoter_count to determine promoter number in yeast genome.

sub promoter_count($$) {
  my($prom_file, $base_dir) = @_;
  chdir $base_dir;
    
  my($prom_count)=0;
      
  open(INFO, $prom_file) || die "--> Cannot open promoter file!\n";
       
  while (<INFO>) { $prom_count++; }
  close INFO;
            
  return $prom_count;
}
              
#------------------------------- Promoter_count sub. output ---------------------------
              
# 4b. Constant PROM_COUNT, the number of promoters in the yeast genome.
              
use constant PROM_COUNT => 5679;

#------------------------------ Binding_summary sub. ---------------------------------

# 5. &binding_summary subroutine to retrieve TF binding statistics from 
# binding affinity matrix.
                                                                                        
sub binding_summary($$) {
  my($self, $tf_conds, $proms) = @_;           

  my($p_thresh) = $self->{parameters}{P_cutoff};    # Get p-value cutoff.
  my($base_dir) = $self->{parameters}{Base_dir};
  my($c_file) = $self->{parameters}{Conv_file};

  my(%tf_conds) = %{ $tf_conds };		    # Retrieve TF condition / promoter row hashes.
  my(%proms) = %{ $proms };

  my($prom_list) = $self->promoter_list(); 	    # Get promoter set of yeast.
  my(%prom_list) = %{ $prom_list };		    # Cast it into a hash.

  my(%TF_hash);					    # Global TF binding statistics recorder.
  my(@p_row);					    # Promoter row of affinity p-values.

  foreach my $prom (sort keys %prom_list) {
    my($prom_obj) = $prom_list{$prom};              # Current promoter object.

    my(@tf_conds);				    # List of TF conditions.
    if (defined $proms{$prom}) {
      @p_row = @{ $proms{$prom} };                  # Retrieve promoter row containing affinity entries.
    
      my($mark)=0;
      foreach my $tf (keys %tf_conds) {
        foreach my $k (keys %{ $tf_conds{$tf} }) {
          my($pos) = $tf_conds{$tf}{$k};
          my($p_val) = $p_row[$pos];

          if ($p_val < $p_thresh && $p_val ne 'NaN') { # If p-value enrty meets treshhold value.
	    push @tf_conds, $tf.'_'.$k; 	       # Record TF-condition pair.
            $mark=1;
          }   
        } 	
      }
      if ($mark==1) { @{ $prom_obj->{tf_cond_list} } = @tf_conds; }
      $TF_hash{$prom} = $prom_obj;	     	    # Record new promoter occurence.
    }
    else { push @{ $self->{no_eq} }, $prom_obj; }   # Record promoters in list, but not in matrix. 
  }
  return \%TF_hash;        	                    # Return promoter database.			
}

#----------------------------------- Ortho_groups sub. -------------------------------

# 6. &ortho_groups subroutine to generate a list of orthogroups for S. cerevisiae
# and other yeast species in the phylogenic tree provided.
                                                                                        
sub ortho_groups() {
  my($self) = shift;

  my($base_dir) = $self->{parameters}{Base_dir};
  my($ortho_file) = $self->{parameters}{Ortho_file};
  
  chdir $base_dir;
  open(INFO, $ortho_file) || die "--> Cannot open orthogroups file -> $ortho_file\n";          
 
  while (<INFO>) {
    if (/^\d+/ && /cer/) {                       # If orthologous promoters on line.
      my(@ortho) = split(/\t/, $_);              # Put line member strings into an array, @ortho.
      
      my(%ortho_list);			         # Orthologous promoter list holder hash.
      foreach my $o (@ortho) {
        if ($o =~ /\|/) { 
          my(@prom) = split(/\|/, $o);
          push @{ $ortho_list{$prom[0]} }, $o;
        }
      }
      push @{ $self->{ortho_list} }, { %ortho_list }; 
    }
  }
  close INFO;
}

#------------------------------- O_promoter_search sub. -------------------------------

# 7. &o_promoter_search subroutine to get promoter sequence of selected yeast organism.
                                                                                        
sub o_promoter_search() {
  my($self) = shift;

  my($dir) = $self->{parameters}{Promoters};
  my($org) = $self->{organism};		    # Get name of yeast species.
	
  my($p_file) = &DNA::Func::get_org_file($dir, $org);

  chdir $dir;
  open(INFO, $p_file) || die "--> Cannot open other promoter file - $p_file\n";          

  my($prom) = $self->{promoter};	    # Get promoter name of yeast organism.

  my($marker);
  PROM: while (<INFO>) {
    my(@prom) = split(/\s/, $_);            # Put line member strings into an array, @prom.
    if ($prom[0] =~ /$prom/) {              # If given promoter matches on any line.
      chomp($self->{prom_seq} = $prom[4]);  # Parse out promoter sequence and assign to Promoter object.

      $marker=1;			    # Mark existance of such Promoter.

      chop $prom[2];
      $self->{chromosome} = $prom[2];       # Assign chromosome no.     
      
      my(@range) = split(/\D/, $prom[3]);
      $self->o_global_range(@range);        # Assign global range.

      last PROM;
    }
  }
  close INFO;
  if(!($marker==1)) { $self->{no_eq} = $prom; }  # Requested promoter does not exists in sequence database.
}

#-------------------------------- O_global range sub. ---------------------------------
              
# 8. &o_global_range subroutine to get promoter range (of not 'cer') according to its
# strand type.
                                                                                        
sub o_global_range(@) {
  my($self, @range) = @_;

  if ($range[1] < $range[2]) { 
    $self->{strand} = 'W';			# Watson strand (+).
    $self->{prom_start} = $range[1];
    $self->{prom_end} = $range[2];
  }  
  else {
    $self->{strand} = 'C';			# Crick strand (-).
    $self->{prom_start} = $range[2];
    $self->{prom_end} = $range[1];
  }  
}
