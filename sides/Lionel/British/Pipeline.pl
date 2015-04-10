#!/usr/bin/perl
# This script find 5' and 3' that's high in DRIP
use strict; use warnings; use mitochy; use Term::ANSIColor qw(:constants);
$Term::ANSIColor::AUTORESET = 1;

my ($input1) = @ARGV;
die "usage: $0 <input1>\n" unless @ARGV;


my $PeakFile = "3T3_DRIP_Peaks_mm9.bed"; 
my $WIGFile  = "/data/mitochi/Work/Project/DRIPc/wig/3T3_mm9DRIP.wig";
my $highDRIPFile = "3T3_high_DRIP_Peaks_mm9.bed"; # This is outputted later
my $highDRIPValue = 10;

print GREEN "
###################################

Input gene annotation: $input1
Peak File: $PeakFile
WIG File: $WIGFile
DRIP Value Threshold: 20


###################################
";
###################################
# Names
my ($folder1, $fileName1) = mitochy::getFilename($input1, "folder");
my ($peakName) = mitochy::getFilename($PeakFile);
my ($wigName)  = mitochy::getFilename($WIGFile);

# 1. Get >30k gene long
system("bedLength.pl $input1 30000 > $fileName1\_30k.temp") == 0 or die "Failed at bedLength.pl: $!\n";

# 2. Map wig to bed and get high DRIP peaks
system("map_wig_to_bed_BIG.pl -w $WIGFile -m $PeakFile") if not -e "$wigName\_$peakName.txt";
open (my $in1, "<", "$wigName\_$peakName.txt") or die;
open (my $out1, ">", $highDRIPFile) or die;
while (my $line = <$in1>) {
	chomp($line);
	my ($chr, $start, $end, $val) = split("\t", $line);
	print $out1 "$chr\t$start\t$end\t3T3_DRIP_mm9\t$val\t\+\n" if $val > $highDRIPValue;
}
close $out1;
close $in1;

# 3. Get clean promoter and terminal that intersect with mm9DRIP
# 3a. Get clean promoter
system("bedtools_bed_change.pl -a -x -10000 -y -1 -i $fileName1\_30k.temp -o $fileName1\_prom.bed");
system("bedtools_bed_change.pl -b -x 1 -y 10000 -i $fileName1\_30k.temp -o $fileName1\_term.bed");
system("bedtools intersect -v -a $fileName1\_prom.bed -b $input1 > $fileName1\_PROM && mv $fileName1\_PROM $fileName1\_prom.bed");
system("bedtools intersect -v -a $fileName1\_term.bed -b $input1 > $fileName1\_TERM && mv $fileName1\_TERM $fileName1\_term.bed");

# 3b. Intersect wih DRIP
system("bedtools_bed_change.pl -b -x -1999 -y 2001 -i $fileName1\_prom.bed -o $fileName1\_prom\.temp");
system("bedtools_bed_change.pl -a -x -2001 -y 1999 -i $fileName1\_term.bed -o $fileName1\_term\.temp");
system("bedtools intersect -wb -a $fileName1\_prom\.temp -b $highDRIPFile > $fileName1\_prom.highDRIP");
system("bedtools intersect -wb -a $fileName1\_term\.temp -b $highDRIPFile > $fileName1\_term.highDRIP");
process_wb("$fileName1\_prom.highDRIP");
process_wb("$fileName1\_term.highDRIP");
# 4. Get promoter/terminal which gene body doesn't have large DRIP
my %prom = %{process_bed("$fileName1\_prom.highDRIP")};
my %term = %{process_bed("$fileName1\_term.highDRIP")};
my %gene = %{process_bed("$input1")};

my %final;
open (my $out3, ">", "$fileName1\_genebody.bed") or die;
foreach my $name (sort {$gene{$a}{chr} cmp $gene{$b}{chr} and $gene{$a}{start} <=> $gene{$b}{start}} keys %gene) {
	if (defined($prom{$name}) and defined($term{$name})) {
		my ($promVal, $termVal) = ($prom{$name}{value}, $term{$name}{value});
		my $promstart 	= $prom{$name}{start};
		my $promend 	= $prom{$name}{end};
		my $termstart 	= $term{$name}{start};
		my $termend 	= $term{$name}{end};

		# Next if prom/term intersect with each other (gene too short)
		next if (intersect($promstart, $promend, $termstart, $termend) == 1);

		my $chr 	= $gene{$name}{chr};
		my $start 	= $gene{$name}{start};
		my $end 	= $gene{$name}{end};
		my $strand 	= $gene{$name}{strand};
		my $combValue	= sqrt((1+$promVal) * (1+$termVal));
		$final{$name} = "$chr\t$start\t$end\t$combValue\t$name\t$strand";
		die if not defined($start);

		# print genebody
		print $out3 "$chr\t$promend\t$termstart\t$name\t$combValue\t$strand\n" if $strand eq "+";
		print $out3 "$chr\t$termend\t$promstart\t$name\t$combValue\t$strand\n" if $strand eq "-";
	}
}
close $out3;

# 5. Intersect genebody with high DRIP and find big region
my @genebody = `bedtools intersect -wao -a $fileName1\_genebody.bed -b $highDRIPFile`;
my %genebody;
foreach my $line (@genebody) {
	chomp($line);
	my ($chr1, $start1, $end1, $name1, $value1, $strand1, $chr2, $start2, $end2, $name2, $value2,$strand2, $length) = split("\t", $line);
	$genebody{$name1} = $length if $length != 0;
}

open (my $out2, ">", "$fileName1\_out.out") or die "Cannot write to $fileName1\_out.out: $!\n";
foreach my $line (keys %final) {
	my ($chr, $start, $end, $value, $name, $strand) = split("\t", $final{$line});
	die "Died at name $name\n" if not defined($start);
	my $length = $end - $start;
	print $out2 "$chr\t$start\t$end\t$name\t$value\t$strand\n" and next if not defined($genebody{$name});
	my $percentile = (1 - ($genebody{$name} / $length))**8 * $value;
	#my $percentile = int($genebody{$name} / $length * 100)/100;
	print $out2 "$chr\t$start\t$end\t$name\t$percentile\t$strand\n";
	#print "$value becomes $percentile at $name $chr\t$start\t$end\n";

}
close $out2;

# 5. Get unique columns
system("unique_column.pl $fileName1\_out.out 5 | sort -k5,5rn > $fileName1.FINAL");

# 6. Post process
mkdir "$fileName1\_TEMPFiles";
system("mv $fileName1\_*.* $fileName1\_TEMPFiles");
print "Output: $fileName1.FINAL\n";

sub intersect {
	my ($start1, $end1, $start2, $end2) = @_;
	return 1 if ($start1 >= $start2 and $start1 <= $end2);
	return 1 if ($start2 >= $start1 and $start2 <= $end1);
	return 0;
}
sub process_bed {
	my ($input) = @_;
	my %bed;
	open (my $in, "<", $input) or die "Cannot read from $input; $!\n";
	while (my $line = <$in>) {
		chomp($line);
		my ($chr, $start, $end, $name, $value, $strand) = split("\t", $line);
		$bed{$name}{value} = $value;
		$bed{$name}{strand} = $strand;
		$bed{$name}{chr} = $chr;
		$bed{$name}{start} = $start;
		$bed{$name}{end} = $end;
	}
	close $in;
	return(\%bed);
}

sub process_wb {
	my ($input) = @_;
	open (my $in, "<", $input) or die "Cannot read from $input; $!\n";
	open (my $out, ">", "$input.temp") or die "Cannot read from $input.temp; $!\n";
	while (my $line = <$in>) {
		chomp($line);
		my ($chr, $start, $end, $name, $value, $strand, $chr2, $start2, $end2, $name2, $value2, $strand2) = split("\t", $line);
		print $out "$chr\t$start\t$end\t$name\t$value2\t$strand\n";
	}
	close $in;
	close $out;
	system("mv $input.temp $input");
}
__END__










open (my $in1, "<", $input1) or die "Cannot read from $input1: $!\n";
while (my $line = <$in1>) {
	chomp($line);
	next if $line =~ /#/;
	my @arr = split("\t", $line);
}
close $in1;

