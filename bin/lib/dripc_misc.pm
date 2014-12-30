package dripc_misc;

############# USAGE ############# 

# Basics
use strict; use warnings;

# Libraries
use Exporter qw(import);
use vars qw(@EXPORT @EXPORT_OK);
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

@EXPORT    = qw(errorMsg);
@EXPORT_OK = qw(print0);
############# SUBS ############# 

sub getsubname {
        return(caller());
}

sub errorMsg {
        print RED "\n#####################################################\n"; print RED "FATAL ERROR on "; print0(); print RED "\nREASON:\n";
}

sub print0 {
        print YELLOW "$0";
}

1;
