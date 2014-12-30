#!/usr/bin/perl

BEGIN {
	use Term::ANSIColor qw(:constants);
	$Term::ANSIColor::AUTORESET = 1;
	use Cwd;
	my ($dir) = getcwd;
	if ($dir !~ /dripc_bin$/) {
		print RED "\n#####################################################\n"; print RED "FATAL ERROR on "; print YELLOW "$0"; print RED "\nREASON:\n";
		print WHITE "\t- Please run "; print YELLOW "$0 "; print WHITE "from dripc_bin folder!\n";
		die "\n";
	}
	print YELLOW "\n$0"; print WHITE ": Running PreCheck.pl to check library consistencies\n";
	system("perl PreCheck.pl -v") == 0 or die "\n";
	
	# Get current directory and get library directory
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
use Cwd;

# Custom Libraries
use mitochy;
use mitochy qw(getFilename);
use R_toolbox;
use class::file;
use functions;
use functions qw(print0);
use database;
use parser;
use Carp;
############# MAIN ############# 

print_help() if ($opt_h);
my ($input) = @ARGV;
pre_check();

# Input files...
my %att;
$att{name}	= getFilename($input);
$att{folder}	= getcwd($input);
$att{md5}	= getmd5($input);
$att{ext}	= getext($input);
my $file = file->new (%att);

my $name   	= $file->name;
my $folder 	= $file->folder;
my $ext		= $file->ext;
my $md5		= $file->md5;
print "NAME $name\nFOLDER $folder\nEXT $ext\nMD5 $md5\n";



############# SUBS ############# 
print GREEN "\nEND OF MAIN\n\n";

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
