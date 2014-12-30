#!/usr/bin/perl

BEGIN {
	use Term::ANSIColor qw(:constants);
	$Term::ANSIColor::AUTORESET = 1;
	print YELLOW "\n$0"; print WHITE ": Running PreCheck.pl to check library consistencies\n";
	system("perl PreCheck.pl -v") == 0 or die "\n";
	# Required libraries
	my @libraries = qw(dripc_parser.pm dripc_database.pm);
	
	# Get current directory and get library directory
	my ($dir) = `pwd` =~ /^(.+)\n$/;
	$dir .= "/lib";

	# Push library directory to @INC
	push (@INC, $dir);
}

############# USAGE ############# 

# Basics
use strict; use warnings;

# Libraries
use Exporter qw(import);
use Cache::FileCache;
use FAlite;
use Statistics::Basic qw(:all);
use Term::ANSIColor qw(:constants);
$Term::ANSIColor::AUTORESET = 1;
use Thread;
use Thread::Queue;
use Getopt::Std;
use vars qw($opt_h);
getopts("h");

# Custom Libraries
use mitochy;
use R_toolbox;
use dripc_misc;
use dripc_misc qw(print0);
use dripc_database;
use dripc_parser;
############# MAIN ############# 

print_help() if ($opt_h);
my ($input) = @ARGV;
pre_check();

# Input files...

############# SUBS ############# 
die "END OF MAIN\n";

sub pre_check {
	if (@ARGV == 0) {
		errorMsg();
		print WHITE "\t- "; print0(); print WHITE ": Argument missing\n\n";
		print RED "SOLUTION:\n";
		print MAGENTA "Specify input file as below:\n";
		print_help();
	}
}
sub print_help {
	print GREEN "\n########### USAGE ############\n";
	print WHITE "\nUsage: $0 <BED FILE>\n";
	print GREEN "\n##############################\n";
	die "\n";
}
