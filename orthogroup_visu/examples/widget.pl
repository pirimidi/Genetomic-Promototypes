#!/usr/local/bin/perl -w
                                                                                        
use strict;
use Tk;

my $mw = MainWindow->new;
my $scrollbar = $mw->Scrollbar(-orient => 'horizontal');
my $entry = $mw->Entry(-width => 30, 
                       -xscrollcommand => ['set' , $scrollbar]);
$scrollbar->configure(-command => ['xview', $entry]);
$scrollbar->pack(-side => 'bottom', -fill => 'x');
$entry->pack(-side => 'bottom', -fill => 'x');

MainLoop;
