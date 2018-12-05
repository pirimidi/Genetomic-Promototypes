#! /usr/bin/perl5 -w

use Tk;
my($mw) = new MainWindow;

my($txt) = $mw -> Scrolled('Text',-scrollbars=>"os") -> pack;

$txt -> insert('0.0',
"Arthur: \"It's at times like this I wish I'd listened to my mother.\"
Ford : \"Why, what did she say?\"
Arthur: \"I don't know, I never listened.\"
		Douglas Adams");

$txt -> insert('10.8', "HELLO MIRKO!!!"); 

MainLoop;
