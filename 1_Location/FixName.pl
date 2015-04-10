#!/usr/bin/perl

use strict; use warnings;

my ($input) = @ARGV;
exit if $input =~ /All/;
die "usage: $0 <dripc_promoter.bed or genomic_promoter.bed>\n" unless @ARGV == 1;
my $check = 0;
fix_names($input);
print "$input was fixed\n" if $check == 1;
print "$input was okay\n" if $check == 0;

sub fix_names {
        # This fix multiple same names on the name column (e.g. ENST00000379694.4 dripc_both has 51x re$
        my ($input) = @_;
        open (my $in, "<", $input) or die "fix_names: Failed to open $input: $!\n";
        open (my $out, ">", "$input\_FixedNames.temp") or die "fix_names: Failed to open $input\_FixedNames.temp: $!\n";
        while (my $line = <$in>) {
                chomp($line);
                if ($line =~ /^#/ or $line =~ /^track/) {
                        print $out "$line\n";
                }
                my ($chr, $start, $end, $names, @others) = split("\t", $line);
                die "fix_names: Died due to not defined names at line $line\n" unless defined($names);
                my $others = join("\t", @others); #just in case others is not BED6

                # Fix multiple same names
                my @names = split(";", $names);
                my @fixed;
		my $check2 = 0;
                for (my $i = 0; $i < @names; $i++) {
			print "$line\n" if $names[$i] !~ /^ENS/ and $input !~ /drip/i and $input !~ /intergenic/;
			$check = 1 if (grep(/^$names[$i]$/, @fixed));
			#print "$line\n" if (grep(/^$names[$i]$/, @fixed)) and $check2 == 0;
			$check2 = 1 if (grep(/^$names[$i]$/, @fixed));
                        push(@fixed, $names[$i]) if not grep(/^$names[$i]$/, @fixed);
			
                }
                my $fixed = join(";", @fixed);

                # Print out
                print $out "$chr\t$start\t$end\t$fixed\t$others\n";
        }
        close $in;
        close $out;
}
