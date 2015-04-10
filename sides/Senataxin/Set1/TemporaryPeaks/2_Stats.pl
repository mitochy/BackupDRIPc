#!/usr/bin/perl
# This script do anova on Senataxin KD DRIP vs scramble DRIP to find genes with significantly lower/higher DRIP
# This produces NOSIG/SIGUP/SIGDOWN.BED
use strict; use warnings; use mitochy;
use Math::CDF qw(:all);

print "This script do anova on Senataxin KD DRIP

Input are *_cut.txt
";

my (@input) = <./*_cut.txt>;
my %data;
foreach my $input1 (sort @input) {
	my ($folder1, $fileName1) = mitochy::getFilename($input1, "folder");
	my ($sample) = $fileName1 =~ /(\w)_OldPeaks_cut/;
	open (my $in1, "<", $input1) or die "Cannot read from $input1: $!\n";
	while (my $line = <$in1>) {
		chomp($line);
		next if $line =~ /#/;
		my ($chr, $start, $end, $value, $names, $zero, $strand, $chr2, $start2, $end2) = split("\t", $line);
		if ($chr2 eq "\.") {
			$chr2 = "NA";
			$start2 = "NA";
			$end2 = "NA";
			$names = "Intergenic";
		}
		my $namez = "$chr\t$start\t$end\t0\t$strand\t$chr2\t$start2\t$end2";
		$data{$namez}{coor} = $namez;
		$data{$namez}{name} = $names;
		$data{$namez}{strand} = $strand;
		$data{$namez}{$sample} = $value;
		$data{$namez}{chr} = $chr;
		$data{$namez}{start} = $start;
		$data{$namez}{end} = $end;
		$data{$namez}{chr2} = $chr2;
		$data{$namez}{start2} = $start2;
		$data{$namez}{end2} = $end2;
	}
	close $in1;
}

print "Output:

1. SIGUP.BED
2. SIGDOWN.BED
3. NOSIG.BED

UCSC-formatted:
4. SIGUP_UCSC.BED
5. SIGDOWN_UCSC.BED
6. NOSIG_UCSC.BED
7. ALLPEAK_UCSC.BED
";
open (my $out1, ">", "SIGUP.BED") or die "Cannot write to SIGUP.BED: $!\n";
open (my $out2, ">", "SIGDOWN.BED") or die "Cannot write to SIGDOWN.BED: $!\n";
open (my $out3, ">", "NOSIG.BED") or die "Cannot write to NOSIG.BED: $!\n";
open (my $out4, ">", "SIGUP_UCSC.BED") or die "Cannot write to SIGUP_UCSC.BED: $!\n";
open (my $out5, ">", "SIGDOWN_UCSC.BED") or die "Cannot write to SIGDOWN_UCSC.BED: $!\n";
open (my $out6, ">", "NOSIG_UCSC.BED") or die "Cannot write to NOSIG_UCSC.BED: $!\n";
open (my $out7, ">", "ALLPEAK_UCSC.BED") or die "Cannot write to ALLPEAK_UCSC.BED: $!\n";
print $out7 "track name=Combined_DRIP_Peak itemRgb=On description=\"Value 99 = Inf, -99 = -Inf. Color code: Green = >2x up in KD (Knockdown, p ANOVA < 0.05), Lightgreen = >2x up in KD (Not Significant), Red = >2x down in KD (p ANOVA < 0.05), Light Red = >2x down in KD (NS), Light Grey = not >2x down or up in KD\"\n";
foreach my $namez (sort {$data{$a}{chr} cmp $data{$b}{chr} || $data{$a}{start} <=> $data{$b}{start}} keys %data) {
	my $chr = $data{$namez}{chr};
	my $start = $data{$namez}{start};
	my $end = $data{$namez}{end};
	my $chr2 = $data{$namez}{chr2}; # GENE CHR
	my $start2 = $data{$namez}{start2}; # GENE START
	my $end2 = $data{$namez}{end2}; # GENE END
	my $name = $data{$namez}{name};
	my $strand = $data{$namez}{strand};
	my $C = $data{$namez}{C};
	my $D = $data{$namez}{D};
	my $E = $data{$namez}{E};
	my $F = $data{$namez}{F};
	my $pval = anova($C, $D, $E, $F);
	my $meanCD = ($C + $D) / 2;
	my $meanEF = ($E + $F) / 2;
	my $fold = $meanCD == 0 ? "Infpos" : $meanEF == 0 ? "Infneg" : log($meanEF / $meanCD) / log(2);
	my $rank = $fold =~ /Inf/ ? $meanCD : $meanCD * $fold;
	my $fold2 = $fold =~ /Infneg/ ? -99 : $fold =~ /Infpos/ ? 99 : 0;
	$strand = "+" if $strand eq ".";
	if ($pval <= 0.05 and ($fold =~ /Infpos/ or ($fold !~ /Inf/ and $fold >= 1)) and ($meanCD > 1.5 or $meanEF > 1.5)) {
		print $out1 "$chr2\t$start2\t$end2\t$name\t$fold\t$strand\t$C\t$D\t$E\t$F\t$meanCD\t$meanEF\t$pval\t$chr\t$start\t$end\t$rank\n";
		print $out4 "$chr\t$start\t$end\t$name\t$fold2\t$strand\t$start\t$end\t0,155,0\n";
		print $out7 "$chr\t$start\t$end\t$name\t$fold2\t$strand\t$start\t$end\t0,155,0\n";
	}
	elsif ($pval <= 0.05 and ($fold =~ /Infneg/ or ($fold !~ /Inf/ and $fold <=- 1)) and ($meanCD > 1.5 or $meanEF > 1.5)) {
		print $out2 "$chr2\t$start2\t$end2\t$name\t$fold\t$strand\t$C\t$D\t$E\t$F\t$meanCD\t$meanEF\t$pval\t$chr\t$start\t$end\t$rank\n";
		print $out5 "$chr\t$start\t$end\t$name\t$fold2\t$strand\t$start\t$end\t255,0,0\n";
		print $out7 "$chr\t$start\t$end\t$name\t$fold2\t$strand\t$start\t$end\t255,0,0\n";
	}
	else {
		print $out3 "$chr2\t$start2\t$end2\t$name\t$fold\t$strand\t$C\t$D\t$E\t$F\t$meanCD\t$meanEF\t$pval\t$chr\t$start\t$end\t$rank\n";

		if ($fold =~ /Infneg/) {
			print $out6 "$chr\t$start\t$end\t$name\t-99\t$strand\t$start\t$end\t253,174,97\n";
			print $out7 "$chr\t$start\t$end\t$name\t-99\t$strand\t$start\t$end\t253,174,97\n";
		}
		elsif ($fold =~ /Infpos/) {
			print $out6 "$chr\t$start\t$end\t$name\t99\t$strand\t$start\t$end\t116,196,118\n";
			print $out7 "$chr\t$start\t$end\t$name\t99\t$strand\t$start\t$end\t116,196,118\n";
		}
		elsif ($fold <= -1) {
			print $out6 "$chr\t$start\t$end\t$name\t$fold\t$strand\t$start\t$end\t253,174,97\n";
			print $out7 "$chr\t$start\t$end\t$name\t$fold\t$strand\t$start\t$end\t253,174,97\n";
		}
		elsif ($fold >= 1) {
			print $out6 "$chr\t$start\t$end\t$name\t$fold\t$strand\t$start\t$end\t116,196,118\n";
			print $out7 "$chr\t$start\t$end\t$name\t$fold\t$strand\t$start\t$end\t116,196,118\n";
		}
		else {
			print $out6 "$chr\t$start\t$end\t$name\t$fold\t$strand\t$start\t$end\t200,200,200\n";
			print $out7 "$chr\t$start\t$end\t$name\t$fold\t$strand\t$start\t$end\t200,200,200\n";
		}
	}

}
close $out1;


sub anova {
	my ($C, $D, $E, $F) = @_;
	my $meanCD = ($C + $D) / 2;
	my $SSwithinCD = ($meanCD - $C)**2 + ($meanCD - $D)**2;
	my $meanEF = ($E + $F) / 2;
	my $SSwithinEF = ($meanEF - $E)**2 + ($meanEF - $F)**2;
	my $dfwithin = (2-1)+(2-1); # Sum of sample number in each group - 1, group 1 = 2 (C and D) etc
	my $MSwithin = ($SSwithinCD + $SSwithinEF) / $dfwithin;
	
	my $grandmean = ($meanCD + $meanEF) / 2;
	my $SSbetween = 2*($meanCD - $grandmean)**2 + 2*($meanEF - $grandmean)**2; # 2*(mean is from number of sample per group
	my $dfbetween = (2-1); # two groups: C/D and E/F - 1
	my $MSbetween = $SSbetween / $dfbetween;
	
	my $Fval = $MSwithin == 0 ? "Inf" : $MSbetween / $MSwithin;
	return($Fval) if $Fval =~ /Inf/;
	my $pval = 1 - pf($Fval, $dfbetween, $dfwithin);
	return($pval);
	die "$C $D $E $F $Fval $pval\n";
}


__END__
P value = [ 1 / Β(ndf/2,ddf/2) ] × [ (ndf × x) / (ndf × x + ddf) ]ndf/2 × [1 - (ndf × x) / (ndf × x + ddf) ]ddf/2 × x-1
