#!/usr/bin/perl
# This script find random genes of each dripc peak
# 1. Parse RNA-seq (/data/mitochi/Work/Project/DRIPc/data/NT2.rpkm)
# 2. Parse and Cache Twin Gene Expression Set (NT2_RNAzero.rpkm and NT2_RNAtwin.rpkm)
#	Index by expression value
#	NT2_RNAtwin.$index
# 3. Parse Genomic Files
# 4. Parse DRIPc Files
# 5. For each DRIPc peak, get RNA seq 
# 6. Original: $dripcname.bed, Output: $dripcname.shuffled
# Run map wig to original and output
use strict; use warnings; use mitochy; use Cache::FileCache;
use Term::ANSIColor qw(:constants);
$Term::ANSIColor::AUTORESET = 1;

my ($input) = @ARGV;
die "usage: $0 <dripc_*.bed>\n" unless @ARGV == 1;
my ($folder, $dripcName) = mitochy::getFilename($input, "folder");

####################################################
# 1. Parse RNA-seq
print BLUE "A. Parsing RNA-seq /data/mitochi/Work/Project/DRIPc/data/NT2.rpkm\n";
my %rna;
open (my $inRNA, "<", "/data/mitochi/Work/Project/DRIPc/data/NT2.rpkm") or die "Cannot read from /data/mitochi/Work/Project/DRIPc/data/NT2.rpkm: $!\n";
while (my $line = <$inRNA>) {
	chomp($line);
	next if $line =~ /#/;

	# Format: GENE\tVALUE\n and GENE is only 1
	my ($gene, $value) = split("\t", $line);
	$rna{$gene} = $value;
}
close $inRNA;
####################################################

####################################################
# 2. Parse RNA-seq Twin into indexed cache hash (It's big!)
print BLUE "B. Parsing and Caching Twin Gene Exprsesion Set NT2_RNAtwin and NT2_RNAzero.rpkm\n";
my $cache = new Cache::FileCache; $cache->set_cache_root("/data/mitochi/Work/Cache/");
my $data  = $cache->get("NT2_RNAtwin");
parse_twin() if (not defined($data)); # If not exist, cache the twin gene
####################################################

####################################################
# 3. Put gene name to column 7 of each DRIPc peak except intergenic by intersecting with its genomic file
my ($feature)   = $dripcName =~ /dripc_(\w+)$/; die "Undefined Feature from dripcname $dripcName (is it dripc_\\w+?)\n" unless defined($feature);
my $genomicFile = "/data/mitochi/Work/Project/DRIPc/1_Location/genomic_$feature.bed"; die "Died undefined genomic $genomicFile for $dripcName\n" unless -e $genomicFile;
runbash("bedtools intersect -s -wb -a $input -b $genomicFile | cut -f1-6,13 > $dripcName.name");

# Antisense_other intersect with antisense of all genic region (genomic_GENIC, which is genic minus antisense (promoter* terminal* genebody both*))
if ($input =~ /antisense_other/) {
	$genomicFile = "/data/mitochi/Work/Project/DRIPc/1_Location/genomic_GENIC.bed";
	runbash("bedtools intersect -S -wb -a $input -b $genomicFile | cut -f1-6,13 > $dripcName.name");
}
$input = "$dripcName.name";
####################################################

####################################################
# 4. Parse Genomic Location
print BLUE "C. Parsing Genomic File $genomicFile Coordinates\n";
my %genomic;
open (my $in2, "<", $genomicFile) or die "Cannot read from $genomicFile: $!\n";
while (my $line = <$in2>) {
	chomp($line);
	my ($chr, $start, $end, $names, $zero, $strand) = split("\t", $line);
	die "DIED AT $line\n" if not defined($strand);
	my @names = split(";", $names);
	foreach my $name (@names) {
		if ($input =~ /antisense_other/) {
			if ($strand eq "+") {
				$strand = "-";
			}
			else {
				$strand = "+";
			}
		}
		$genomic{$name} = "$chr\t$start\t$end\t$name\t0\t$strand\n";
	}
}
close $in2;
####################################################

# 5. Parse DRIPc data and sort by RNA-seq
print BLUE "D. Parsing DRIPc File $input Coordinates\n";
open (my $in3, "<", $input) or die "Cannot read from $input: $!\n";
my %gene;
while (my $line = <$in3>) {
	chomp($line);
	next if $line =~ /#/;
	my ($chr, $start, $end, $name, $value, $strand, $genes) = split("\t", $line);
	my @names = split(";", $genes);

	# Since a genomic region might contain overlappnig genes, use gene with highest RNA to be used to find its twin
	my $rna = 0;
	my $gene;
	$gene = "NA" if $input =~ /intergenic/;
	foreach my $genename (@names) {
		if ($input !~ /intergenic/ and $rna <= $rna{$genename}) {
			next if not defined($genomic{$genename});
			$rna = $rna{$genename};
			$gene = $genename;
		}
	}
	next if not defined($gene);
	$gene{$name}{rna}     = $rna;
	$gene{$name}{genomic} = $line;
	$gene{$name}{gene}    = $gene;
}
####################################################

# 6. Shuffle
# This is separate from above because we need to sort by rna-seq for cache purpose otherwise it'll be helllla slow
print BLUE "E. Shuffling\n";
open (my $out, ">", "$dripcName\_Shuffled.bed") or die "Cannot write to $dripcName\_Shuffled.bed: $!\n";
if ($input =~ /intergenic/) {

	# Get intergenic genomic coordinates as twin
	my @twin;
	foreach my $gene (keys %genomic) {
		push(@twin, $gene);
	}
	
	# Shuffle
	my $count_Done = 0;
	foreach my $dripName (keys %gene) {
		my $total = (keys %gene);
		printf "Intergenic Done = $count_Done / $total (%.2f %%)\n", int($count_Done / $total * 10000) / 100 if $count_Done % 200 == 0;
		my $iter = 0;
		my ($chr, $start, $end, $name, $value, $strand, $genes) = split("\t", $gene{$dripName}{genomic});
		next if $value < 5;
		for (my $i = 0; $i < 5000; $i++) {
			my $randGene = $twin[int(rand(@twin))];
			my $genomic     = $genomic{$randGene};
			if (not defined($randGene) or not defined($genomic)) {
				next;
			}
			chomp($genomic);
			my ($chr2, $start2, $end2, $name2, $zero2, $strand2) = split("\t", $genomic);
			die "COOR $genomic CHR $chr2 START $start2 END $end2 NAME $name2 ZERO $zero2 STRAND $strand2\n" if not defined($strand2);
			die "Died at gene $randGene\n" if not defined($chr2);
			my $length  = $end - $start;
			my $length2 = $end2 - $start2;
	
			my $spacer = $length2 - $length;
			my $randStart = $start2 + int(rand($spacer));
			my $randEnd   = $randStart + $length;
			if ($spacer < 1) {
				if ($iter == 10) {
					$iter = 0;
					$randStart = $start2;
					$randEnd   = $end2;
				}
				$iter ++;
				$i --;		
				next;
			}
	
			# Put random strand (intergenic can be positive or negative strand)
			$strand2 = "-" if rand() < 0.5;
			print $out "$chr2\t$randStart\t$randEnd\t$name\t$value\t$strand2\tORIG=$name,TWIN=$name2\n";
		}
		
		$count_Done ++;
	}
}
else {
	my %data;
	my $last_index = "INIT";
	my $count_notused_lowdripc = 0;
	my $count_notused_lowtwin = 0;
	my $twin_notused = 0;
	my %twin_notused;
	foreach my $dripName (sort {$gene{$a}{rna} <=> $gene{$b}{rna}} keys %gene) {
		my ($chr, $start, $end, $name, $value, $strand, $genes) = split("\t", $gene{$dripName}{genomic});
		
		$count_notused_lowdripc ++ and next if $value < 0;
		my $gene = $gene{$dripName}{gene};
		my $length = $end - $start;
	
		# Define Index and get the cache
		my $index = 0;
		my $rna = $gene{$dripName}{rna};
		die "Undefined rna of $dripName of gene $gene\n" if not defined($rna);
		if ($rna >= 10) {
			$index = get_index($rna);
		}
		
		if ($index ne $last_index) {
			if ($index == 0) {
				print "Doing index $index\n";
				my $data = $cache->get("NT2_RNAzero");
				@{$data{zero}} = @{$data};
			}
			else {
				print "Doing index $index\n";
				my $data = $cache->get("NT2_RNAtwin.$index");
				die "INDEX $index not found\n" if not defined($data);
				%data = %{$data};
			}
		}
		else {
			#print "DRIPNAME $dripName GENE $gene index $index\n";
		}
		$last_index = $index;
	
		# Get its gene exp twin (or zero)	
		my @twin = $rna >= 10 ? @{$data{$gene}} : @{$data{zero}};
		$count_notused_lowtwin ++ and next if @twin < 20;
		my $iter = 0;
		my $twin_iter = 0;
		for (my $i = 0; $i < 5000; $i++) {
			my $randGene = $twin[int(rand(@twin))];
			my $genomic     = $genomic{$randGene};
			if (not defined($randGene) or not defined($genomic)) {
				if ($twin_iter == 10) {
					$twin_iter = 0;
					next;
				}
				$twin_iter ++;
				$i --;
				next;
			}
			chomp($genomic);
			my ($chr2, $start2, $end2, $name2, $zero2, $strand2) = split("\t", $genomic);
			die "Died at gene $randGene\n" if not defined($chr2);
			die "Undefined strand2 at genomic $genomic CHR $chr2 START $start2 END $end2 NAME $name2 ZERO $zero2 STRAND $strand2\n" if not defined($strand2);
			my $length2 = $end2 - $start2;
	
			my $spacer = $length2 - $length;
			my $randStart = $start2 + int(rand($spacer));
			my $randEnd   = $randStart + $length;
			if ($spacer < 1) {
				if ($iter == 10) {
					$iter = 0;
					$randStart = $start2;
					$randEnd   = $end2;
				}
				$iter ++;
				$i --;		
				next;
			}
			$twin_notused{$dripName} ++;
			print $out "$chr2\t$randStart\t$randEnd\t$name\t$value\t$strand2\tORIG=$gene,TWIN=$randGene\n";
		}
		
	}
	close $in3;

	my $count_notused_lowtwin2 = 0;
	foreach my $dripName (keys %twin_notused) {
		$count_notused_lowtwin2 ++ if $twin_notused{$dripName} < 20;
	}

	print "Not Used DRIPc < 5: $count_notused_lowdripc, Number of Twin < 50: $count_notused_lowtwin, Twin < 50: $count_notused_lowtwin2\n";
}
close $out;
####################################################

# 7. Below is to filter shuffled peaks that has too high DRIPc and genes that has than 100 shuffles
open (my $in10, "<", "$dripcName\_Shuffled.bed") or die;
my %data;
while (my $line = <$in10>) {
	chomp($line);
	my ($chr, $start, $end, $name) = split("\t", $line);
	$data{$name}{count}++;
	push(@{$data{$name}{line}}, $line);
}
close $in10;
open (my $out5, ">", "$dripcName.shuffled") or die "Cannot write to $dripcName.shuffled: $!\n";
foreach my $name (keys %data) {
	if ($data{$name}{count} >= 1000) {
		for (my $i = 0; $i < 1000; $i++) {
			print $out5 "$data{$name}{line}[$i]\n";
		}
	}
}
close $out5;
runbash("cat $dripcName.shuffled | sort -k1,1 -k2,2n > $dripcName.tmp && mv $dripcName.tmp $dripcName.shuffled");


###############
# SUBROUTINES #
###############

sub parse_twin {
	# Data Hash
	my %data;
	my @data;

	# 1. Processing NT2 RNA lses than 10
	print "\t1. Processing NT2_RNAzero.rpkm\n";
	my $inputZero = "NT2_RNAzero.rpkm";
	open (my $in0, "<", $inputZero) or die "Cannot read from $inputZero: $!\n";
	while (my $line = <$in0>) {
		chomp($line);
		next if $line =~ /#/;
		push(@data, $line);
	}
	close $in0;
	$cache->set("NT2_RNAzero", \@data);
	@data = ();
	
	print "\t2. Processing NT2Twin.rpkm\n";
	my $inputtwin = "NT2_RNAtwin.rpkm";
	open (my $in1, "<", $inputtwin) or die "Cannot read from $inputtwin: $!\n";
	my $last_index = "INIT";
	while (my $line = <$in1>) {
		chomp($line);
		next if $line =~ /#/;
		my ($gene, @twin) = split("\t", $line);

		# If number of gene twin is less than 20, don't use it
		my $rna = $rna{$gene};

		# Get index based on gene rpkm
		my $index = get_index($rna);
		print "\tProcessing index $index\n" if $last_index ne $index;
		if ($index ne $last_index and $last_index ne "INIT") {
			$cache->set("NT2_RNAtwin.$last_index", \%data);
			%data = ();
		}
		@{$data{$gene}} = @twin;
		$last_index = $index;
	}
	$cache->set("NT2_RNAtwin.$last_index", \%data);
	%data = ();
	close $in1;
	
	$cache->set("NT2_RNAtwin", 1);
}

sub runbash {
        my ($cmd) = @_;
        print "\t$cmd\n";
        system($cmd) == 0 or die "Failed to run $cmd: $!\n";
}

sub get_index {
	my ($rna) = @_;
	my $index = 0;
	$index = 10 + int($rna / 10)    if $rna >= 10   and $rna < 100;
	$index = 20 + int($rna / 100)   if $rna >= 100  and $rna < 1000;
	$index = 30 + int($rna / 1000)  if $rna >= 1000 and $rna < 10000;
	$index = 40 + int($rna / 10000) if $rna >= 10000;
	return($index);
}
