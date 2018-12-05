#! /usr/bin/perl5 -w

use strict;
use Tk;

my $main = new MainWindow;
my $box = $main->Listbox(-relief => 'ridge',
                         -width => 75, 
                         -height => 20,
                         -setgrid => 'yes');

my @items = qw(One Two Three Four Five Six Seven Eight Nine Ten Eleven Twelve);

foreach (@items) {
  $box->insert('end', $_);
}
        
my $scroll = $main->Scrollbar(-command => ['yview', $box],
                              -orient => 'vertical');

  $box->configure(-yscrollcommand => ['set', $scroll]);
  $box->pack(-side => 'left', -fill => 'both', -expand => 'yes');
  $scroll->pack(-side => 'right', -fill => 'y');

MainLoop;
