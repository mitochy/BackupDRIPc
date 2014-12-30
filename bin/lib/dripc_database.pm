package dripc_database;

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

my $caller = getsubname();

############# SUBS ############# 

sub getsubname {
        return(caller());
}

sub print0 {
	print YELLOW "$caller";
}
1;
