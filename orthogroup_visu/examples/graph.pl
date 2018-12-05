#!/usr/local/bin/perl -w
use strict;
use Tk;
my($o, $s) = (250, 20);
my($pi, $x, $y) = (3.1415926, 0);
my $mw = MainWindow->new;
my $c = $mw->Canvas(-width => 500, -height => 500);
$c->pack;
$c->createLine(50, 250, 450, 250);
$c->createText(10, 250, -fill => 'blue', -text => 'X');
$c->createLine(250, 50, 250, 450);
$c->createText(250, 10, -fill => 'blue', -text => 'Y');

for ($x = -(3*$pi); $x <= +(3*$pi); $x += 0.1) {
    $y = sin($x);
    $c->createText( $x*$s+$o, $y*$s+$o, -fill => 'red', -text => '.');
    $y = cos($x);
    $c->createText( $x*$s+$o, $y*$s+$o, -fill => 'green', -text => '.');
}

MainLoop;
