#!/usr/bin/perl

use strict; use warnings;

=command
#promoter
`bedtools intersect -a E14_promoter.bed -b 3T3_promoter.bed > promoter_E143T3.BED`;
`bedtools intersect -a NT2_promoter.bed -b 3T3_promoter.bed > promoter_NT23T3.BED`;
#`bedtools intersect -a Fib_promoter.bed -b 3T3_promoter.bed > promoter_Fib3T3.BED`;
`bedtools intersect -a E14_promoter.bed -b NT2_promoter.bed > promoter_E14NT2.BED`;
#`bedtools intersect -a Fib_promoter.bed -b NT2_promoter.bed > promoter_FibNT2.BED`;
`bedtools intersect -a 3T3_promoter.bed -b NT2_promoter.bed > promoter_3T3NT2.BED`;
#`bedtools intersect -a E14_promoter.bed -b Fib_promoter.bed > promoter_E14Fib.BED`;
#`bedtools intersect -a NT2_promoter.bed -b Fib_promoter.bed > promoter_NT2Fib.BED`;
#`bedtools intersect -a 3T3_promoter.bed -b Fib_promoter.bed > promoter_3T3Fib.BED`;
`bedtools intersect -a NT2_promoter.bed -b E14_promoter.bed > promoter_NT2E14.BED`;
#`bedtools intersect -a Fib_promoter.bed -b E14_promoter.bed > promoter_FibE14.BED`;
`bedtools intersect -a 3T3_promoter.bed -b E14_promoter.bed > promoter_3T3E14.BED`;

#promoter_ext
`bedtools intersect -a E14_promoter_ext.bed -b 3T3_promoter_ext.bed > promoter_ext_E143T3.BED`;
`bedtools intersect -a NT2_promoter_ext.bed -b 3T3_promoter_ext.bed > promoter_ext_NT23T3.BED`;
#`bedtools intersect -a Fib_promoter_ext.bed -b 3T3_promoter_ext.bed > promoter_ext_Fib3T3.BED`;
`bedtools intersect -a E14_promoter_ext.bed -b NT2_promoter_ext.bed > promoter_ext_E14NT2.BED`;
#`bedtools intersect -a Fib_promoter_ext.bed -b NT2_promoter_ext.bed > promoter_ext_FibNT2.BED`;
`bedtools intersect -a 3T3_promoter_ext.bed -b NT2_promoter_ext.bed > promoter_ext_3T3NT2.BED`;
#`bedtools intersect -a E14_promoter_ext.bed -b Fib_promoter_ext.bed > promoter_ext_E14Fib.BED`;
#`bedtools intersect -a NT2_promoter_ext.bed -b Fib_promoter_ext.bed > promoter_ext_NT2Fib.BED`;
#`bedtools intersect -a 3T3_promoter_ext.bed -b Fib_promoter_ext.bed > promoter_ext_3T3Fib.BED`;
`bedtools intersect -a NT2_promoter_ext.bed -b E14_promoter_ext.bed > promoter_ext_NT2E14.BED`;
#`bedtools intersect -a Fib_promoter_ext.bed -b E14_promoter_ext.bed > promoter_ext_FibE14.BED`;
`bedtools intersect -a 3T3_promoter_ext.bed -b E14_promoter_ext.bed > promoter_ext_3T3E14.BED`;

#terminal
`bedtools intersect -a E14_terminal.bed -b 3T3_terminal.bed > terminal_E143T3.BED`;
`bedtools intersect -a NT2_terminal.bed -b 3T3_terminal.bed > terminal_NT23T3.BED`;
#`bedtools intersect -a Fib_terminal.bed -b 3T3_terminal.bed > terminal_Fib3T3.BED`;
`bedtools intersect -a E14_terminal.bed -b NT2_terminal.bed > terminal_E14NT2.BED`;
#`bedtools intersect -a Fib_terminal.bed -b NT2_terminal.bed > terminal_FibNT2.BED`;
`bedtools intersect -a 3T3_terminal.bed -b NT2_terminal.bed > terminal_3T3NT2.BED`;
#`bedtools intersect -a E14_terminal.bed -b Fib_terminal.bed > terminal_E14Fib.BED`;
#`bedtools intersect -a NT2_terminal.bed -b Fib_terminal.bed > terminal_NT2Fib.BED`;
#`bedtools intersect -a 3T3_terminal.bed -b Fib_terminal.bed > terminal_3T3Fib.BED`;
`bedtools intersect -a NT2_terminal.bed -b E14_terminal.bed > terminal_NT2E14.BED`;
#`bedtools intersect -a Fib_terminal.bed -b E14_terminal.bed > terminal_FibE14.BED`;
`bedtools intersect -a 3T3_terminal.bed -b E14_terminal.bed > terminal_3T3E14.BED`;

#terminal_ext
`bedtools intersect -a E14_terminal_ext.bed -b 3T3_terminal_ext.bed > terminal_ext_E143T3.BED`;
`bedtools intersect -a NT2_terminal_ext.bed -b 3T3_terminal_ext.bed > terminal_ext_NT23T3.BED`;
#`bedtools intersect -a Fib_terminal_ext.bed -b 3T3_terminal_ext.bed > terminal_ext_Fib3T3.BED`;
`bedtools intersect -a E14_terminal_ext.bed -b NT2_terminal_ext.bed > terminal_ext_E14NT2.BED`;
#`bedtools intersect -a Fib_terminal_ext.bed -b NT2_terminal_ext.bed > terminal_ext_FibNT2.BED`;
`bedtools intersect -a 3T3_terminal_ext.bed -b NT2_terminal_ext.bed > terminal_ext_3T3NT2.BED`;
#`bedtools intersect -a E14_terminal_ext.bed -b Fib_terminal_ext.bed > terminal_ext_E14Fib.BED`;
#`bedtools intersect -a NT2_terminal_ext.bed -b Fib_terminal_ext.bed > terminal_ext_NT2Fib.BED`;
#`bedtools intersect -a 3T3_terminal_ext.bed -b Fib_terminal_ext.bed > terminal_ext_3T3Fib.BED`;
`bedtools intersect -a NT2_terminal_ext.bed -b E14_terminal_ext.bed > terminal_ext_NT2E14.BED`;
#`bedtools intersect -a Fib_terminal_ext.bed -b E14_terminal_ext.bed > terminal_ext_FibE14.BED`;
`bedtools intersect -a 3T3_terminal_ext.bed -b E14_terminal_ext.bed > terminal_ext_3T3E14.BED`;

#genebody
`bedtools intersect -a E14_genebody.bed -b 3T3_genebody.bed > genebody_E143T3.BED`;
`bedtools intersect -a NT2_genebody.bed -b 3T3_genebody.bed > genebody_NT23T3.BED`;
#`bedtools intersect -a Fib_genebody.bed -b 3T3_genebody.bed > genebody_Fib3T3.BED`;
`bedtools intersect -a E14_genebody.bed -b NT2_genebody.bed > genebody_E14NT2.BED`;
#`bedtools intersect -a Fib_genebody.bed -b NT2_genebody.bed > genebody_FibNT2.BED`;
`bedtools intersect -a 3T3_genebody.bed -b NT2_genebody.bed > genebody_3T3NT2.BED`;
#`bedtools intersect -a E14_genebody.bed -b Fib_genebody.bed > genebody_E14Fib.BED`;
#`bedtools intersect -a NT2_genebody.bed -b Fib_genebody.bed > genebody_NT2Fib.BED`;
#`bedtools intersect -a 3T3_genebody.bed -b Fib_genebody.bed > genebody_3T3Fib.BED`;
`bedtools intersect -a NT2_genebody.bed -b E14_genebody.bed > genebody_NT2E14.BED`;
#`bedtools intersect -a Fib_genebody.bed -b E14_genebody.bed > genebody_FibE14.BED`;
`bedtools intersect -a 3T3_genebody.bed -b E14_genebody.bed > genebody_3T3E14.BED`;

#intergenic
`bedtools intersect -a E14_intergenic.bed -b 3T3_intergenic.bed > intergenic_E143T3.BED`;
`bedtools intersect -a NT2_intergenic.bed -b 3T3_intergenic.bed > intergenic_NT23T3.BED`;
#`bedtools intersect -a Fib_intergenic.bed -b 3T3_intergenic.bed > intergenic_Fib3T3.BED`;
`bedtools intersect -a E14_intergenic.bed -b NT2_intergenic.bed > intergenic_E14NT2.BED`;
#`bedtools intersect -a Fib_intergenic.bed -b NT2_intergenic.bed > intergenic_FibNT2.BED`;
`bedtools intersect -a 3T3_intergenic.bed -b NT2_intergenic.bed > intergenic_3T3NT2.BED`;
#`bedtools intersect -a E14_intergenic.bed -b Fib_intergenic.bed > intergenic_E14Fib.BED`;
#`bedtools intersect -a NT2_intergenic.bed -b Fib_intergenic.bed > intergenic_NT2Fib.BED`;
#`bedtools intersect -a 3T3_intergenic.bed -b Fib_intergenic.bed > intergenic_3T3Fib.BED`;
`bedtools intersect -a NT2_intergenic.bed -b E14_intergenic.bed > intergenic_NT2E14.BED`;
#`bedtools intersect -a Fib_intergenic.bed -b E14_intergenic.bed > intergenic_FibE14.BED`;
`bedtools intersect -a 3T3_intergenic.bed -b E14_intergenic.bed > intergenic_3T3E14.BED`;

# three way intersect
`bedtools intersect -a promoter_E143T3.BED -b NT2_promoter.bed > promoter_E143T3NT2.BED`;
`bedtools intersect -a promoter_ext_E143T3.BED -b NT2_promoter_ext.bed > promoter_ext_E143T3NT2.BED`;
`bedtools intersect -a terminal_E143T3.BED -b NT2_terminal.bed > terminal_E143T3NT2.BED`;
`bedtools intersect -a terminal_ext_E143T3.BED -b NT2_terminal_ext.bed > terminal_ext_E143T3NT2.BED`;
`bedtools intersect -a genebody_E143T3.BED -b NT2_genebody.bed > genebody_E143T3NT2.BED`;
`bedtools intersect -a intergenic_E143T3.BED -b NT2_intergenic.bed > intergenic_E143T3NT2.BED`;
=cut

my @sample = qw(E14 3T3 NT2);
my @feature = qw(promoter promoter_ext terminal terminal_ext genebody intergenic);

# Get total from each original .bed
my %total;
print "Original\n";
foreach my $feature (@feature) {
	foreach my $sample (@sample) {
		my $count = `bedtools_bed_stats.pl $sample\_$feature.bed`;
		($count) = $count =~ /^$sample\_$feature.bed\s+(\d+) bp/;
		$count /= 1000000;
		$total{$sample}{$feature} = $count;
		print "\t$feature\t$sample\t$count\n";
	}
}

print "Two-way\n";
foreach my $feature (@feature) {			
	my @used;
	foreach my $sample (@sample) {
		foreach my $sample2 (@sample) {
			next if $sample eq $sample2;
			next if grep(/^$sample\_$sample2$/, @used);
			next if grep(/^$sample2\_$sample$/, @used);
			my $count = `bedtools_bed_stats.pl $feature\_$sample$sample2.BED`;
			($count) = $count =~ /^$feature\_$sample$sample2.BED\s+(\d+) bp/;
			$count /= 1000000;
			my $total = $total{$sample}{$feature};
			#my $perc  = int(1000*$count / $total)/10;
			print "\t$feature\t$sample\t$sample2\t$count\n";
			push(@used, "$sample\_$sample2");
			push(@used, "$sample2\_$sample");
		}
	}
}

print "Three-way\n";
foreach my $feature (@feature) {
	my $count = `bedtools_bed_stats.pl $feature\_E143T3NT2.BED`;
	($count) = $count =~ /^$feature\_E143T3NT2.BED\s+(\d+) bp/;
			$count /= 1000000;
	print "\t$feature\tTHREE\t$count\n";
}

__END__

GOOGLE APPS VENN:
https://chart.googleapis.com/chart?cht=v&chs=300x100&chd=t:8.232851,14.258268,17.940902,5.328091,4.844798,7.591592,3.578145&chdl=E14|3T3|NT2|Promoter
https://chart.googleapis.com/chart?cht=v&chs=300x100&chd=t:&chdl=E14|3T3|NT2|Promoter_ext
https://chart.googleapis.com/chart?cht=v&chs=300x100&chd=t:&chdl=E14|3T3|NT2|Terminal
https://chart.googleapis.com/chart?cht=v&chs=300x100&chd=t:&chdl=E14|3T3|NT2|Terminal_ext
https://chart.googleapis.com/chart?cht=v&chs=300x100&chd=t:&chdl=E14|3T3|NT2|Genebody
https://chart.googleapis.com/chart?cht=v&chs=300x100&chd=t:&chdl=E14|3T3|NT2|Intergenic
