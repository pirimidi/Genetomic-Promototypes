#-----------------------------------------------------------
# Author: Mirko Palla.
# Date: February 28, 2006.
# For: Task 3 at the Church Lab - Genetics Department, 
# Harvard Medical School.
# 
# Purpose: This program removes .~ files in a specified 
# directory in Perl.
#------------------------------------------------------------ 

#!/usr/bin/perl

use strict;

my(@director) = ('base', 'data_miner', 'logo', 'evolution', 'statistics');
my($base_dir) = "/root/Church_lab/Polypromoter/task_7/a_sharp/";

my($i)=0;
foreach my $dir (@director) {
  chdir $base_dir.$director[$i];
  system("rm *~");
  system("rm DEADJOE");

  chdir $base_dir.$director[$i]."/DNA/";
  system("rm *~");
  system("rm DEADJOE");

  if ($dir eq 'logo') {
    chdir $base_dir.$director[$i]."/EXT/";
    system("rm *~");
    system("rm DEADJOE");
  }
  
  $i++;
}

chdir $base_dir;
system("clear");
system("pwd");
