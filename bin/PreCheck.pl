#!/usr/bin/perl

#-------------Usage------------
# Basics
use strict; use warnings;

# Libraries
use Cache::FileCache;
use FAlite;
use Statistics::Basic qw(:all);
use Term::ANSIColor qw(:constants);
$Term::ANSIColor::AUTORESET = 1;
use Thread;
use Thread::Queue;
use Getopt::Std;
use vars qw($opt_v);
getopts("v");

# Required libraries
my @libraries = qw(dripc_parser.pm dripc_database.pm dripc_misc.pm);
	
# Get current directory and get library directory
my ($dir) = `pwd` =~ /^(.+)\n$/;
$dir .= "/lib";

# Push library directory to @INC
push (@INC, $dir);

# Check if library directory $dir is in perl library @INC
if ($opt_v) {
	print "\n"; print0();
	print WHITE ": Checking if library "; print BLUE "$dir "; print WHITE "is in perl library "; print BLUE "\@INC\n";
	if (not grep(/^$dir$/, @INC)) {
		print RED "\tERROR:   "; print BLUE "$dir "; print WHITE "is not in "; print BLUE "\@INC!\n";
	}
	else {
	print GREEN "\tSUCCESS: "; print BLUE "$dir "; print WHITE "is in "; print BLUE "\@INC!\n";
	}
}
die if not grep(/^$dir$/, @INC);

# Check if all @libraries is in $dir
if ($opt_v) {print "\n"; print0(); print WHITE ": Checking existance of libraries in"; print BLUE " $dir\n";}
my @libLibraryCheck;
for (my $i = 0; $i < @libraries; $i++) {
	my $library = $libraries[$i];
	if (not -e "$dir/$library") {
		push(@libLibraryCheck, "$library");
		if ($opt_v) {print RED "\tERROR:   "; print WHITE "$i. Library "; print BLUE "$library "; print WHITE "does not exist!\n";}
	}
	else {
		if ($opt_v) {print GREEN "\tSUCCESS: "; print WHITE "$i. Library "; print BLUE "$library "; print WHITE "exist!\n";}
	}
}
print WHITE "\n" if $opt_v;

# Otherwise, ask user to put all libraries into their perl library themselves
if (@libLibraryCheck != 0) {
	print RED "\n-----------------------------------------------------\n"; print RED "FATAL ERROR on "; print0(); print RED "\nREASON: "; print WHITE "Libraries missing or failed to put some libraries to "; print BLUE "$dir:\n"; 
	for (my $i = 0; $i < @libLibraryCheck; $i++) {
		print BLUE "\t$i. $libLibraryCheck[$i]\n";
	}
	print WHITE "\nSOLUTION: Please manually copy the missing libraries from "; print BLUE "$dir "; print WHITE "to one of these Perl library folders:\n";
	foreach my $INCLIB (@INC) {
		chomp($INCLIB);
		next if $INCLIB !~ /^\/usr/;
		next if $INCLIB =~ /^\/usr/ and $INCLIB !~ /local/;
		print BLUE "\t$INCLIB\n";
	}
	foreach my $INCLIB (@INC) {
		chomp($INCLIB);
		next if $INCLIB =~ /^\/usr/;
		print BLUE "\t$INCLIB\n";
	}
	print WHITE "\n";
	die "\n";
}

sub print0 {
	print YELLOW "$0";
}
