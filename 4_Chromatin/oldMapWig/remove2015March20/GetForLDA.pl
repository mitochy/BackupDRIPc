#!/usr/bin/perl
# Version: 15 Oct 2014
# - Parse Original
# - Parse Shuffled
# - Get the first 150 shuffled (otherwise discard)
use strict; use warnings; use mitochy; use R_toolbox; use Statistics::Test::WilcoxonRankSum; use Math::CDF qw (:all); use Statistics::Multtest qw(:all); use Statistics::Normality qw(:all);
my ($origInput) = @ARGV;
die "usage: $0 <Histone_origInput.txt>\n" unless @ARGV == 1;

my ($folder, $fileName) = mitochy::getFilename($origInput, "folder");
my $shuffleInput = "MappedShuffled/$fileName.txt"; die "Shuffled $shuffleInput does not exist\n" if not -e $shuffleInput;
my ($histone, $feature) = $fileName =~ /^(\w+)_drip\w?_(\w+)$/; 
($histone, $feature) = $fileName =~ /^(\w+)_enh\w*_(\w+)$/ if not defined($histone);
die "Histone and/or Feature undefined\n" unless defined($histone) and defined($feature);
my ($rnaFile) = "/data/mitochi/Work/Project/DRIPc/data/NT2.rpkm";
my %rna = %{parse_rna($rnaFile)};
my $usd = 0;
##########################
# 1. Parse Shuffled
print "1. Processing $shuffleInput\n";
my %data;
my %names;
open (my $in1, "<", $shuffleInput) or die "Cannot read from $shuffleInput: $!\n";
my $linecount = 0;
while (my $line = <$in1>) {
	chomp($line);
	next if $line =~ /#/;
	$linecount ++;
	my ($chr, $start, $end, $value, $name, $val, $strand, $info) = split("\t", $line);
	my ($gene) = $info =~ /TWIN=(\w+\.\d+),chr/;
	($gene) = $info =~ /TWIN=(\w+\.\d+)/ if not defined($gene);
	my ($orig) = $info =~ /ORIG=(\w+\.\d+),chr/;
	($orig) = $info =~ /ORIG=(\w+\.\d+)/ if not defined($orig);
	$names{$name}{orig} = $orig;
	$names{$name}{gene} = $gene;
	die "Undefined RNA seq for gene $gene shuffled\n" unless defined($rna{$gene});
	my $rna = $rna{$orig};
	next if $rna < 10;
	#next if $rna <95 or $rna > 105;
	my @arr = split("\t", $line);
	#push(@{$data{$name}{random}}, $value);
	$value = int($value * 100)/100;
	push(@{$data{$orig}{random}}, "$line\t$linecount");
}
close $in1;
##########################

##########################
# 2. Parse Original
print "2. Processing $origInput\n";
open (my $in2, "<", $origInput) or die "Cannot read from $origInput: $!\n";
my ($totalDRIPc, $undefined) = (0,0);
while (my $line = <$in2>) {
	chomp($line);
	next if $line =~ /#/;
	my ($chr, $start, $end, $value, $name) = split("\t", $line);
	# If the name isn't the same as orig
	$name = $names{$name}{orig} if $name !~ /ENS/;
	print "Undefined orig for name $name\n" and next if (not defined($name));
	$value = int($value * 100)/100;
	my @arr = split("\t", $line);	
	$totalDRIPc ++;
	#next if $rna{$name} >10;# and $rna{$name} > 10;
	# Next if shuffled of this data does not exist
	# Because DRIPc too low or number of twin is less than 50 or less than 100 shuffles
	if (not defined($data{$name})) {
		$undefined ++;
		next;
	}
	else {
		$data{$name}{orig} = $line;
	}
}
close $in2;
printf "DRIPc peak not used due to no shuffles: $undefined / $totalDRIPc (%.2f %%)\n", 100 * $undefined / $totalDRIPc;
##########################

##########################
# 3. 
print "3. Print Out Shuffles\n";
my %namez;
if (-e "RANDOM_SHUFFLES.data") {
	open (my $ins, "<", "RANDOM_SHUFFLES.data");
	while (my $line = <$ins>) {
		chomp($line);
		my ($origName, $shufName, $lineCount) = split("\t", $line);
		$namez{$lineCount}{orig} = $origName;
		$namez{$lineCount}{shuf} = $shufName;
	}	
}
else {
	open (my $outz, ">", "RANDOM_SHUFFLES.data") or die;
	my %check;
	foreach my $name (keys %data) {
		my $valueOrig = $data{$name}{orig};
		next if not defined($valueOrig);
		next if @{$data{$name}{random}} < 5;#125;
		$usd++;
	
		my @random = shuffle(@{$data{$name}{random}});
		for (my $i = 0; $i < @random; $i++) {
			my ($chr, $start, $end, $value, $name, $val, $strand, $info, $lineCount) = split("\t", $random[$i]);
			my ($shuf) = $info =~ /TWIN=(\w+\.\d+),chr/;
			($shuf) = $info =~ /TWIN=(\w+\.\d+)/ if not defined($shuf);
			my ($orig) = $info =~ /ORIG=(\w+\.\d+),chr/;
			($orig) = $info =~ /ORIG=(\w+\.\d+)/ if not defined($orig);
			if (keys %{$check{$orig}} == 5) {
				last;
			}
			elsif (not defined($check{$orig}{$shuf})) {
				$check{$orig}{$shuf} = $lineCount;
				$namez{$lineCount}{orig} = $orig;
				$namez{$lineCount}{shuf} = $shuf;
				print $outz "$orig\t$shuf\t$lineCount\n";
			}
			else {
				next;
			}
		}
		
	}
}


open (my $outShuf, ">", "$fileName\_shuf5.txt") or die;
open (my $outOrig, ">", "$fileName\_orig5.txt") or die;
foreach my $name (keys %data) {
	my $check = 0;
	for (my $i = 0; $i < @{$data{$name}{random}}; $i++) {
		my ($chr, $start, $end, $value, $name, $val, $strand, $info, $lineCount) = split("\t", $data{$name}{random}[$i]);
		if (defined($namez{$lineCount}{orig}) and $namez{$lineCount}{orig} =~ /\w+/) {
			my ($shuf) = $info =~ /TWIN=(\w+\.\d+),chr/;
			($shuf) = $info =~ /TWIN=(\w+\.\d+)/ if not defined($shuf);
			my ($orig) = $info =~ /ORIG=(\w+\.\d+),chr/;
			($orig) = $info =~ /ORIG=(\w+\.\d+)/ if not defined($orig);
			die "NAME OF LINECOUNT $lineCount ORIG $namez{$lineCount}{orig} isn't the same as name in shuffle file $orig\n" if $namez{$lineCount}{orig} ne $orig;
			die "NAME OF LINECOUNT $lineCount shuf $namez{$lineCount}{shuf} isn't the same as name in shuffle file $shuf\n" if $namez{$lineCount}{shuf} ne $shuf;
			print $outShuf "$chr\t$start\t$end\t$value\t$name\t$val\t$strand\t$info\n";
			$check ++;
		}
	}
	print $outOrig "$data{$name}{orig}\n" if $check > 0;
}

die "There is no peak with more than 125 shuffles\n" if $usd == 0;

sub mean {
	my (@value) = @_;
	my $mean;
	foreach my $value (@value) {
		print "MEAN: UNDEFINED VALUE\n" and next if not defined($value);
		if ($value eq "NA" or $value =~ /inf/i) {
			$value = 0;
			print "NON DIGIT VALUE FOUND AT $value\n";
		}
		$mean += $value / @value;
	}
	return($mean);
}

sub median {
	my (@value) = @_;
	@value = sort {$a <=> $b} (@value);
	my $median = $value[int(@value/2)];
	return($median);
}

sub shuffle {
	my (@value) = @_;
	#print "Before: @value\n";
	for (my $i = 0; $i < 10000; $i++) {
		my $rand1 = int(rand(@value));
		my $rand2 = int(rand(@value));
		my $val1 = $value[$rand1];
		my $val2 = $value[$rand2];
		$value[$rand1] = $val2;
		$value[$rand2] = $val1;
	}
	#print "After: @value\n";
	return(@value);
}

sub parse_rna {
        my ($input) = @_;
        my %data;
        open (my $in, "<", $input) or die "Cannot read from $input: $!\n";
        while (my $line = <$in>) {
                chomp($line);
                next if $line =~ /#/;
                my ($gene, $val) = split("\t", $line);
                $data{$gene} = $val;
        }
        close $in;
        return(\%data);
}
__END__
