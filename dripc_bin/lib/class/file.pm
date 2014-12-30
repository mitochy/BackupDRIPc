package file;
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

my $caller = getsubname();

############# SUBS ############# 

sub getsubname {
        return(caller());
}

sub new {
	my $class = shift;
	my $self = {@_};
	bless ($self, $class);
	return($self);
}

sub name {
	my $self = shift;
	my $var = $self->{name};
	return($var) if defined($var);
	return("NA");
}

sub ext {
	my $self = shift;
	my $var = $self->{ext};
	return($var) if defined($var);
	return("NA");
}

sub folder {
	my $self = shift;
	my $var = $self->{folder};
	return($var) if defined($var);
	return("NA");
}

sub md5 {
	my $self = shift;
	my $var = $self->{md5};
	return($var) if defined($var);
	return("NA");
}

1;
