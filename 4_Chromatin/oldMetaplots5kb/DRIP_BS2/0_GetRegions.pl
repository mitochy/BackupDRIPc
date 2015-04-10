#!/usr/bin/perl
# This script take SHUFFLED files, reduce the shuffled into 15 (randomly), and change coordinate to +/- 5kb
# For promoter/antisense/terminal, can be centered on TSS or the peak

use strict; use warnings; use mitochy;

# Get shuffled regions and turns them into +/- 5kb
my @files = @ARGV;
die "usage: $0 <input.shuffled>\n" unless @ARGV;

foreach my $file (@files) {
	print "Processing file $file\n";
	# Promoter/term/antisense need to get fixed (+/-2kb of TSS/TTS) and +/- 5kb
	# Genebody just take 5 shuffles and print out
	# Intergenic just print out 
	#next if $file =~ /both/;
	#next if $file !~ /promoter\./ and $file !~ /antisense\./ and $file !~ /terminal\./;
	getCoorFixed($file);
}

sub getCoorFixed {
	my ($input) = @_;
	my ($folder, $fileName) = mitochy::getFilename($input, "folder");
	open (my $out1, ">", "$fileName\_orig.tmp") or die "Cannot write to $fileName\_orig.tmp: $!\n";
	open (my $out2, ">", "$fileName\_shuf.tmp") or die "Cannot write to $fileName\_shuf.tmp: $!\n";

	##############################################################

	# GET GENE COORDINATES #
	# DEFAULT (DONT USE THIS THOUGH)
	my $gencode = "/data/mitochi/Work/Project/DRIPc/bed/hg19_gencode.bed";

	# PROMOTER/TERMINAl/ANTISENSE CENTERED on TSS/TTS
	$gencode = "/data/mitochi/Work/Project/DRIPc/bed/hg19_gencode_promoter.bed"  if $input =~ /promoter/;
	$gencode = "/data/mitochi/Work/Project/DRIPc/bed/hg19_gencode_terminal.bed"  if $input =~ /terminal/;
	$gencode = "/data/mitochi/Work/Project/DRIPc/bed/hg19_gencode_antisense.bed" if $input =~ /antisense/;

	# PROMOTER/TERMINAl/ANTISENSE CENTERED ON DRIPBS PEAK
	#$gencode = "/data/mitochi/Work/Project/DRIPc/1_Location/4_CreateRegion_DRIP/Region/NT2_DRIPBS_Peak_region_promoter.bed" if $input =~ /promoter\./;
	#$gencode = "/data/mitochi/Work/Project/DRIPc/1_Location/4_CreateRegion_DRIP/Region/NT2_DRIPBS_Peak_region_terminal.bed" if $input =~ /terminal\./;
	#$gencode = "/data/mitochi/Work/Project/DRIPc/1_Location/4_CreateRegion_DRIP/Region/NT2_DRIPBS_Peak_region_antisense.bed" if $input =~ /antisense\./;
	
	# GENEBODY
	$gencode = "/data/mitochi/Work/Project/DRIPc/1_Location/4_CreateRegion_DRIP/Region/NT2_DRIPBS_Peak_region_genebody.bed" if $input =~ /genebody/;

	# INTERGENIC
	$gencode = "/data/mitochi/Work/Project/DRIPc/1_Location/4_CreateRegion_DRIP/Region/NT2_DRIPBS_Peak_region_intergenic.bed" if $input =~ /intergenic/;
	
	my %data;
	open (my $in1, "<", $gencode) or die "Cannot read from $gencode: $!\n";
	while (my $line = <$in1>) {
		chomp($line);
		next if $line =~ /#/;
		next if $line =~ /track/;
		my ($chr, $start, $end, $name, $val, $strand) = split("\t", $line);
		$data{$name} = $line;
	}
	close $in1;

	##############################################################
	
	# Process input file so we can get random 15 shuffles	
	my %lines;
	open (my $in2, "<", $input) or die "Cannot read from $input: $!\n";
	while (my $line = <$in2>) {
		chomp($line);
		next if $line =~ /#/;
		next if $line =~ /track/;
		my ($chr, $start, $end, $name, $val, $strand, $info) = split("\t", $line);
		my ($orig, $shuf) = $info =~ /ORIG=(\w+\.?\d*),chr.+TWIN=(\w+\.?\d*),chr/;
		#next if defined($bad{$orig}) and $bad{$orig} == 1;
		#next if defined($bad{$shuf}) and $bad{$shuf} == 1;
		next if $chr eq "chrM";
		push(@{$lines{$name}}, $line) if not defined($lines{$name});
	}
	close $in2;

	my %usedOrig;
	foreach my $name (keys %lines) {
		my @val = @{$lines{$name}};
		# Get only 5 shuffles
		#my $shuffle_times = 5;
		#my $iter = 0;
		# my @line;
		#while (@line < $shuffle_times) {
		#	my $rand = int(rand(@val));
		#	last if $iter > 10;
		#	$iter ++ and next if $val[$rand] =~ /chrM/;
		#	push(@line, $rand) if not grep(/^$rand$/, @line);
		#}
		#print "Name $name nexted\n" and next if $iter == 10;
		for (my $i = 0; $i < @val; $i++) {
			my $line = $val[$i];
			my ($chr, $start, $end, $name, $val, $strand, $info) = split("\t", $line);
			die "Died at line $line\n" if not defined($info);
			my ($orig, $shuf) = $info =~ /ORIG=(\w+\.?\d*),chr.+TWIN=(\w+\.?\d*),chr/;
			die "Died at getting original/shuffled gene name at file $input: Failed at getCoorFixed\n" if not defined($orig) or not defined($shuf);
	
			#if ($input =~ /promoter\./ or $input =~ /terminal\./ or $input =~ /antisense\./) {
			#	my $origCoor = $data{$orig}; my ($chrOrig, $startOrig, $endOrig, $nameOrig, $valOrig, $strandOrig) = split("\t", $origCoor);
			#	my $shufCoor = $data{$shuf}; my ($chrShuf, $startShuf, $endShuf, $nameShuf, $valShuf, $strandShuf) = split("\t", $shufCoor);
			#	die "Died at getting original/shuffled gene coordinates at file $input: Failed at getCoorFixed\n" if not defined($origCoor) or not defined($shufCoor);
			#	print $out1 "$chrOrig\t$startOrig\t$endOrig\t$orig\t$valOrig\t$strandOrig\t$info\n" if not defined($usedOrig{$orig});
			#	print $out2 "$chrShuf\t$startShuf\t$endShuf\t$shuf\t$valShuf\t$strand\t$info\n";
			#	$usedOrig{$orig} = 1;
			#}
			#else {
				my $origCoor = $data{$name}; 
				$origCoor = $data{$orig} if not defined($origCoor);
				my ($chrOrig, $startOrig, $endOrig, $nameOrig, $valOrig, $strandOrig) = split("\t", $origCoor);
				die "Died at getting original/shuffled gene coordinates at file $input\n$line\n: Failed at getCoorFixed\n" if not defined($origCoor);
				print $out1 "$chrOrig\t$startOrig\t$endOrig\t$orig\t$valOrig\t$strandOrig\t$info\n" if not defined($usedOrig{$name});
				print $out2 "$chr\t$start\t$end\t$shuf\t$val\t$strand\t$info\n";
				$usedOrig{$name} = 1;
			#}
		}
	}
	
	runbash("bedtools_bed_change.pl -c -x -5000 -y 5000 -i $fileName\_shuf.tmp -o $fileName\_shuf.input");
	runbash("bedtools_bed_change.pl -c -x -5000 -y 5000 -i $fileName\_orig.tmp -o $fileName\_orig.input");
}

sub runbash {
        my ($cmd) = @_;
        print "\t$cmd\n";
        system($cmd) == 0 or die "Failed to run $cmd: $!\n";
}














__END__
my ($input) = @ARGV;
die "usage: $0 <input>\n" unless @ARGV;

my ($folder, $fileName) = mitochy::getFilename($input, "folder");

open (my $out, ">", "$fileName.out") or die "Cannot write to $fileName.out: $!\n";

open (my $in, "<", $input) or die "Cannot read from $input: $!\n";
while (my $line = <$in>) {
	chomp($line);
	next if $line =~ /#/;
	my @arr = split("\t", $line);
}

close $in;
close $out;
