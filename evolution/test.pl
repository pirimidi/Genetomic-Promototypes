#!/usr/bin/perl

use strict;

# ----- LINE[2] ------
my($num) = 0.0000048948;
if($num =~ /^\d/) { print "\n--> THIS IS A MULT_DIG: $num\n"; }

my($non) = '[-----------]';
if($non =~ /^\d+$/) { print "\n--> THIS IS A MULT_DIG: $non\n"; }

# ----- LINE[5] ------
my(@two) = split(//, '32');
if($two[0] =~ /^\d/) { print "\n--> THIS IS IN c3: $two[0]\n"; }
if($two[1] =~ /^\d/) { print "\n--> THIS IS IN c2: $two[1]\n\n"; }

my(@o_p_s1) = split(/\,/, "cer|YBL021C");
my(@o_p_s2) = split(/\,/, "alb|YBL021C, alb|YGH098W");

print "--> I: @o_p_s1\n";
print "--> II: @o_p_s2\n\n";

my($multi_tab) = "STE12\t\t4\t\t1\t\t2\t\t1\n";
print "\nMulti tab line: $multi_tab";

chomp(my(@tab_line) = split(/\t\t/, $multi_tab));
print "Multi line as array: @tab_line\n";

print "Multi line[0]: $tab_line[0]\n";
print "Multi line[1]: $tab_line[1]\n";
print "Multi line[2]: $tab_line[2]\n";
print "Multi line[3]: $tab_line[3]\n";
print "Multi line[4]: $tab_line[4]\n";

my(%hash)=[];
if(defined %hash) { print "\n--> HASH IS DEFINED!\n\n"; }
