#---------------------------------------------------------------------------
# Author: Mirko Palla.
# Date: April 25, 2006.
# For: Task 3 at the Church Lab - Genetics Department, 
# Harvard Medical School.
# 
# Purpose: This program contains the complete code for 
# class Logo, modeling an object for motif overlap statistics in Perl.
#
# Associated subroutines:   0. &constructor	  1. &create_PDF
#			    2. &get_font	  3. &create_header
#			    4. &insert_logos	  5. &insert_text
#		            6. &insert_config
#
# {Note: Make sure [Logo_dir] field is specified in configuration file.}
#--------------------------------------------------------------------------- 

package EXT::Logo;
use base ("EXT::Func");

use strict;

1;

#------------------------------------ Class Logo --------------------------------------

# 0. constructor for a overlap visualization object (PDF).

sub new {
  my($class, $o_set) = @_;
  my($self)={};

  $self->{overlaps} = $o_set; 				# Overlap_set field.

  $self->{pdf} = ();		 			# PDF object field.
  $self->{info} = ();		 			# Information field.
  $self->{file} = ();		 			# File name field.

  $self->{page} = ();		 			# PDF page field.
  $self->{text} = ();		 			# PDF text field.
  $self->{font} = ();		 			# Font type field.

  $self->{bar_code} = ();		 		# Bar code field.
  $self->{parameters} = ();		 		# Parameter list field.

  bless($self,$class);
  return($self);
}

#---------------------------------- PDF constants ----------------------------------

use constant H_HS => 700;                               # Header height start.

use constant L_LS1 => 60;                               # 1st logo left start.
use constant L_LS2 => 335;                              # 2nd logo left start.
use constant L_H => 120;                                # Logo height.
use constant L_L => 200;                                # Logo lenght.

use constant S_S => 5;                                  # Separator start.
use constant S_E => 590;                                # Separator end.

#---------------------------------- Create_PDF sub. --------------------------------

# 1. &create_PDF subroutine to construct PDF object with its appropiate attributes.
                                                                                        
sub create_PDF() {
  my($self) = @_;

  my($f_name) = $self->{parameters}{File_name};  
  my($out_dir) = $self->{parameters}{Out_dir}; 
       
  my($o_laps) = $self->{overlaps};
  my($random) = $o_laps->{bar_code};			   # Get bar code.
  my($gene) = $o_laps->{prom_obj}->{gene};                 # Get gene name.
              
  my($title) = $gene."_OL-remove_".$f_name.		   # Promoter mutation with removed units.
                               "_".$random.".pdf";
  $self->{file} = $title;				   # Assign file name to Logo object.

  chdir $out_dir;
  my($pdf) = PDF::API2->new(-file => $title);
  $pdf->mediabox(595,842);

  my($theTime) = &EXT::Func::get_time();

  my(%info) = $pdf->info( 'Author' => "Mirko Palla",
                          'CreationDate' => $theTime );

  $self->{pdf} = $pdf;
  $self->{info} = \%info;
}

#-------------------------------- Get_font sub. ----------------------------------

# 2. &get_font subroutine to assign different font values available to PDF.

sub get_font() {
  my($self) = @_;
  my($pdf) = $self->{pdf};

  my(%font) = (
    Arial => {
        Bold => $pdf->corefont('Arial-Bold', -encoding => 'latin1'),
        Roman => $pdf->corefont('Arial', -encoding => 'latin1'),
    },

    Courier => {
        Bold => $pdf->corefont('Courier-Bold', -encoding => 'latin1'),
        Roman => $pdf->corefont('Courier', -encoding => 'latin1'),
    },

    Times => {
        Bold => $pdf->corefont('Times-Bold', -encoding => 'latin1'),
        Roman => $pdf->corefont('Times', -encoding => 'latin1'),
        Italic => $pdf->corefont('Times-Italic', -encoding => 'latin1'),
    },
  );

  $self->{font} = \%font;
}

#-------------------------------- Create_header sub. ----------------------------------

# 3. &create_header subroutine to construct header of page.
                                                                                        
sub create_header() {
  my($self, $i) = @_;

  my($pdf) = $self->{pdf};				  # PDF object.
  my(%font) = %{ $self->{font} };			  # Font palette.
  my(%info) = %{ $self->{info} };			  # Doc. information.

  my($page) = $pdf->page;
  $self->{page} = $page;

  $pdf->preferences( -onecolumn => 1,
                     -afterfullscreenoutlines => 1,
                     -firstpage => [ $page, -fit => 1],
                   );

  my($title_box) = $page->gfx;
  $title_box->fillcolor('lightblue');
  $title_box->rect(S_S, H_HS, 585, 132);
  $title_box->fill;

  &EXT::Func::draw_line($page, 'red', S_S, S_E, H_HS);
  &EXT::Func::draw_line($page, 'red', S_S, S_E, 832);

  #---------------------------------- Text parameters --------------------------------

  my($txt) = $page->text;
  $self->{text} = $txt;

  $txt->fillcolor('black');

  #--------------------------------- Header parameters -------------------------------

  $txt->textstart;
  $txt->font( $font{'Arial'}{'Roman'}, 14 );

  $txt->translate(297,20);
  $txt->text($i+1);

  $txt->translate(30,800);
  $txt->text( "Created by: $info{'Author'}" );

  $txt->translate(570,800);
  $txt->text_right( "On: $info{'CreationDate'}" );

  $txt->font( $font{'Arial'}{'Bold'}, 20 );
  $txt->translate(300,746);
  $txt->text_center("Remove_overlap.pl statistics");

  $txt->textend;
}

#-------------------------------- Insert_logos sub. ----------------------------------

# 4. &insert_logos subroutine to insert logos into page.
                                                                                        
sub insert_logos($$$) {
  my($self, $L_HS, $image_1, $image_2) = @_;

  my($pdf) = $self->{pdf};				  # PDF object.
  my($page) = $self->{page};

  #-------------------------------- Image parameters ----------------------------------

  chdir $self->{parameters}{Logo_dir};

  my($logo_1) = $pdf->image_png($image_1);
  my($logo_2) = $pdf->image_png($image_2);

  my($image) = $page->gfx;

  $image->image($logo_1, L_LS1, $L_HS, L_L, L_H);
  $image->image($logo_2, L_LS2, $L_HS, L_L, L_H);
}

#---------------------------------- Insert_text sub. ----------------------------------

# 5. &insert_text subroutine to insert motif related text into page.
                                                                                        
sub insert_text($$$) {
  my($self, $L_HS, $m_prom, $pos) = @_;

  #------------------------- Parameters -----------------------------

  my($m_A) = $m_prom->{overlap_pair}->{motif_A};   # Untouched motif.
  my($m_B) = $m_prom->{overlap_pair}->{motif_B};   # Mutated motif.
            
  my($OL)  = $m_prom->{overlap_pair}->{overlap};   # Overlap value.
  my($sign) = $m_prom->{overlap_pair}->{sign}; 	   # Overlap type: +/-.
  my($type) = $m_prom->{overlap_pair}->{type}; 	   # Overlap type: {"fake","real->I","real->II"}.
                                   
  #------------------ Motif A & B -----------------#
                                                   #
         my($m_A_factor) = $m_A->{factor};    	   #
         my($m_A_motif) = $m_A->{motif};    	   #
         my($m_A_offset) = $m_A->{offset};    	   #
         my($m_A_orient) = $m_A->{orient};   	   #
         my($m_A_score) = $m_A->{score};     	   #
                                                   #
         my($m_B_factor) = $m_B->{factor};         #
         my($m_B_motif) = $m_B->{motif};           #
         my($m_B_offset) = $m_B->{offset};         #
         my($m_B_orient) = $m_B->{orient};   	   # 
         my($m_B_score) = $m_B->{score};     	   #
                                                   #                                  
  #------------------------------------------------#
  
  my($A) = $m_A_factor.$m_A_orient;        	   # [motif_name.orient] - [STE12+]
  my($B) = $m_B_factor.$m_B_orient; 

  my($text_h) = "[$pos] OVERLAP: $OL -> {$B ($m_B_offset) <-> $A ($m_A_offset)} => TYPE: $type";

  #----------------------- Mutated motif ----------------------------
      
  my($s_kmer) = $m_prom->{s_kmer};          	   # Assign substitition kmer.
  my($s_score) = $m_prom->{s_score};        	   # Assign substitution p-value.

  if($m_B_orient eq '-') { $s_kmer = &DNA::Func::rev_complement($s_kmer); }
                      
  my($un_l) = "Mutated motif: [$B]";
  my($se_l) = "O: $m_B_motif";
  my($th_l) = "P: $m_B_score";
  my($fo_l) = "M: $s_kmer";
  my($fi_l) = "P: $s_score";

  #---------------------- Dominant motif -----------------------------
      
  my($D_kmer) = $m_prom->{D_kmer};          	   # Assign substitition kmer.
  my($D_score) = $m_prom->{D_score};        	   # Assign substitution p-value.
                      
  if($m_A_orient eq '-') { $D_kmer = &DNA::Func::rev_complement($D_kmer); }

  my($un_s) = "Untouch motif: [$A]";
  my($se_s) = "O: $m_A_motif";
  my($th_s) = "P: $m_A_score";
  my($fo_s) = "M: $D_kmer";
  my($fi_s) = "P: $D_score";

  my($repet);
  if ( defined @{ $m_prom->{untouched_tf} } ) {
   
    $repet = "U: ";
    my($count) = 1;
    
    foreach my $u_tf (@{ $m_prom->{untouched_tf} }) {  # Access 'untouched_tf' field in Mutant object.
                        
      my($fac) = $u_tf->{factor};
      my($off) = $u_tf->{offset};
      my($ori) = $u_tf->{orient};
                                                                      
      $repet = $repet."[$fac$ori] ($off) ";
      $count++;
    }
  }

  #-------------------------- LOGO top ------------------------------
    
  my($txt) = $self->{text};
  my(%font) = %{ $self->{font} };
            
  $txt->textstart;
  $txt->font( $font{'Arial'}{'Bold'}, 16 );
           
  $txt->translate(300, $L_HS+170);
  $txt->text_center($text_h);
                            
  $txt->translate(163, $L_HS+130);
  $txt->text_center($m_B_factor);
                                   
  $txt->translate(443, $L_HS+130);
  $txt->text_center($m_A_factor);

  #---------------------- LOGO bottom (left) ------------------------
    
  my($t_start) = $L_HS-30;

  $txt->font( $font{'Courier'}{'Bold'}, 16 );
  $txt->translate(L_LS1, $t_start);
  $txt->text($un_l);
                   
  $txt->font( $font{'Courier'}{'Roman'}, 16 );
                       
  $txt->translate(L_LS1, $t_start-18);
  $txt->text($se_l);
                                
  $txt->translate(L_LS1, $t_start-32);
  $txt->text($th_l);
                                       
  $txt->translate(L_LS1, $t_start-48);
  $txt->text($fo_l);
                                               
  $txt->translate(L_LS1, $t_start-62);
  $txt->text($fi_l);

  #---------------------- LOGO bottom (right) ------------------------

  $txt->font( $font{'Courier'}{'Bold'}, 16 );
  $txt->translate(L_LS2, $t_start);
  $txt->text($un_s);

  $txt->font( $font{'Courier'}{'Roman'}, 16 );

  $txt->translate(L_LS2, $t_start-18);
  $txt->text($se_s);

  $txt->translate(L_LS2, $t_start-32);
  $txt->text($th_s);

  $txt->translate(L_LS2, $t_start-48);
  $txt->text($fo_s);

  $txt->translate(L_LS2, $t_start-62);
  $txt->text($fi_s);

  $txt->translate(L_LS2, $t_start-78);
  $txt->text($repet);
  $txt->textend;

  #--------------------------------- Red separator ------------------------------------

  my($page) = $self->{page};
  my($l_height) = $t_start-95;

  &EXT::Func::draw_line($page, 'red', S_S, S_E, $l_height);
}                                           

#-------------------------------- Insert_config sub. ----------------------------------

# 6. &insert_config subroutine to insert configuration parameters into zeroth page.
                                                                                        
sub insert_config($) {
  my($self) = @_;

  use constant L_S => 60; 
  use constant H_S => 500; 
  #use constant L_S => 300; 

  #------------------------- Parameters -----------------------------

  my($o_lap) = $self->{overlaps};

  my($gene) = $o_lap->{prom_obj}->{gene};                 # Get gene name.
  my($prom) = $o_lap->{prom_obj}->{promoter};             # Get promoter name.
  my($chro) = $o_lap->{prom_obj}->{chromosome};           # Get chromosome no.

  my($header) = "{ GENE: $gene <=> PROMOTER: $prom } -> CHROMOSOME: $chro\n";

  my($txt) = $self->{text};
  my(%font) = %{ $self->{font} };
            
  $txt->textstart;

  $txt->font( $font{'Arial'}{'Bold'}, 16 );
  $txt->translate(300,600);
  $txt->text_center($header);

  $txt->font( $font{'Arial'}{'Bold'}, 14 );
  $txt->translate(300,550);
  $txt->text_center("Configuration file parameters");

  $txt->font( $font{'Arial'}{'Roman'}, 14 );

  my($iter)=0;
  my(%param) = %{ $self->{parameters} };
  foreach my $k (keys %param) {
  
    my($no) = $iter+1;
    my($text) = "$no. $k => $param{$k}"; 
  
    $txt->translate(L_S, H_S-($iter*20));
    $txt->text($text);
    #$txt->text_center($text);
    
    $iter++;
  }
}                                           
