package dripc_file;
# Class for storing data about a file

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

my $caller = caller();

############# SUBS ############# 

sub getsubname {
        return(caller());
}

sub new {
	my $self = {};
	bless ($self, "dripc_file");
	return($self);
}


1;
