#-------------------------------------------------------------------------------
# Author: Mirko Palla.
# Date: April 25, 2006.
# For: Task 3 at the Church Lab - Genetics Department, 
# Harvard Medical School.
# 
# Purpose: This program contains the complete code for 
# module Func, containing re-occurring functionality 
# subroutines in Perl.
#
# Associated subroutines:   1. &get_time	2. &draw_line
#------------------------------------------------------------------------------- 

package EXT::Func;
                                                                                        
use strict;

1;

#------------------------------------ Get_time sub. -------------------------------------

# 1. &get_time subroutine to get current time.

sub get_time() {

  my(@months) = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
  my(@weekDays) = qw(Sun Mon Tue Wed Thu Fri Sat Sun);

  my($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, 
              $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();

  my($year) = 1900 + $yearOffset;

  if($minute<10) { $minute = '0'.$minute; }
  if($second<10) { $second = '0'.$second; }

  my($theTime) = "$hour:$minute:$second, $weekDays[$dayOfWeek], $months[$month] $dayOfMonth, $year";
  return $theTime;
}

#------------------------------------ draw_line sub. -------------------------------------

# 2. &draw_line subroutine to draw line as specified.

sub draw_line($$$$$) {
  my($page, $color, $S_S, $S_E, $H_HS) = @_;

  my($red_line) = $page->gfx;
  $red_line->strokecolor($color);

  $red_line->move($S_S, $H_HS);                           # Start coordinate (origin in x,y cords) of pen.
  $red_line->line($S_E, $H_HS);                           # End poit of sliding (in x, y cords).

  $red_line->stroke;  
}
