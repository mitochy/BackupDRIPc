package functions;
# Misc functions (directly exported to main script)

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

sub getmd5 {
	my ($input) = @_;
	my ($md5) = `md5sum $input` =~ /^(\w+) /;
	return($md5);
}

sub getext {
	my ($input) = @_;
	my ($ext) = $input =~ /\.(\w+)$/;
	return($ext) if defined($ext);
	return("NA");
}
1;
