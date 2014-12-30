package parser;
# Parse stuff

############# USAGE ############# 

# Basics
use strict; use warnings;

# Libraries
use Exporter qw(import);
use vars qw(@EXPORT);
use Cache::FileCache;
use FAlite;
use Statistics::Basic qw(:all);
use Term::ANSIColor qw(:constants);
$Term::ANSIColor::AUTORESET = 1;
use Thread;
use Thread::Queue;
use Getopt::Std;

# Custom Libraries
use mitochy;
use R_toolbox;
use functions;

my $caller = getsubname();

# Put all functions into @EXPORT
my @functions = `grep -P "^sub " lib/$caller.pm | grep -v getsubname | grep -v print0`;
foreach my $function (@functions) {
	my ($function_name) = $function =~ /^sub (\w+) {/;
	if (not defined($function_name)) {
		errorMsg();
		print WHITE "\t- "; print0(); print WHITE ": Not defined function name of function "; print RED "$function"; print WHITE "\n";
		die "\n";
	}
	push (@EXPORT, $function_name);
}

############# SUBS ############# 

sub getsubname {
	return(caller());
}

sub print0 {
	print YELLOW $caller;
}

sub parse {
	

}

1;
