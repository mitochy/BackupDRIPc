#!/usr/bin/perl
# This script find random genes of each dripc peak
# 1. Parse RNA-seq (/data/mitochi/Work/Project/DRIPc/data/$RNAFile)
# 2. Parse and Cache Twin Gene Expression Set (E14_RNAzero.rpkm and E14_RNAtwin.rpkm)
#	Index by expression value
#	E14_RNAtwin.$index
# 3. Parse Genomic Files
# 4. Parse DRIPc Files
# 5. For each DRIPc peak, get RNA seq 
# 6. Original: $dripcname.bed, Output: $dripcname.shuffled
# Run map wig to original and output
use strict; use warnings; use mitochy; use Cache::FileCache; use Getopt::Std;
use Term::ANSIColor qw(:constants);
$Term::ANSIColor::AUTORESET = 1;

my ($input) = @ARGV;
die "usage: $0 <../1_Location/drip3T3mm9_*.bed>\n" unless @ARGV == 1;
my ($folder, $dripcName) = mitochy::getFilename($input, "folder");

my $shuffleNumber = 500;
####################################################
# 1. Parse RNA-seq
my $RNAFile = $input =~ /E14/ ? "/data/mitochi/Work/Project/DRIPc/data/E14.rpkm" : $input =~ /3T3/ ? "/data/mitochi/Work/Project/DRIPc/data/3T3.rpkm" : die "Input must be 3T3 or E14!\n";
print BLUE "A. Parsing RNA-seq /data/mitochi/Work/Project/DRIPc/data/$RNAFile\n";
my %rna;
open (my $inRNA, "<", "$RNAFile") or die "Cannot read from $RNAFile: $!\n";
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
print BLUE "B. Parsing and Caching Twin Gene Exprsesion Set E14_RNAtwin and E14_RNAzero.rpkm\n";
my $cache = new Cache::FileCache; $cache->set_cache_root("/data/mitochi/Work/Cache/");
my $data  = $cache->get("E14_RNAtwin") if $input =~ /E14/;
$data  = $cache->get("3T3_RNAtwin") if $input =~ /3T3/;
parse_twin() if (not defined($data)); # If not exist, cache the twin gene
####################################################

####################################################
# 3. Put gene name to column 7 of each DRIPc peak except intergenic by intersecting with its genomic file
my ($feature)   = $dripcName =~ /drip[0-9a-zA-Z]*_(\w+)$/; die "Undefined Feature from dripcname $dripcName (is it dripc_\\w+?)\n" unless defined($feature);
my $genomicFile = "/data/mitochi/Work/Project/DRIPc/Mouse/1_Location/genomic_$feature.bed"; die "Died undefined genomic $genomicFile for $dripcName\n" unless -e $genomicFile;
runbash("bedtools intersect -wb -a $input -b $genomicFile | cut -f1-6,13 > $dripcName.name");
#runbash("./getStrand.pl $dripcName.name ../../bed/mm9_gencode.bed $RNAFile");
# DEPRECATED: Antisense_other intersect with antisense of all genic region (which is genic minus antisense and promoter (terminal* genebody both*))
#if ($input =~ /antisense_other/) {
#	$genomicFile = "/data/mitochi/Work/Project/DRIPc/1_Location/genomic_antisense_other.bed";
#	runbash("bedtools intersect -s -wb -a $input -b $genomicFile | cut -f1-6,13 > $dripcName.name");
#}
$input = "$dripcName.name";
####################################################

####################################################
# 4. Parse Genomic Location, including mm9_gencode_promoter and terminal
# For promoter, antisense, and terminal, we match shuffles to the same location as DRIPc peak relative to TSS/TTS
print BLUE "C. Parsing Genomic File $genomicFile Coordinates\n";
my %genomic;
my %prom;
my %term;
my %anti;
##my %temp;
if ($input !~ /promoter\./ and $input !~ /antisense\./ and $input !~ /terminal\./) {
	open (my $in2, "<", $genomicFile) or die "Cannot read from $genomicFile: $!\n";
	while (my $line = <$in2>) {
		chomp($line);
		my ($chr, $start, $end, $names, $zero, $strand) = split("\t", $line);
		die "Parsing genomcFile $genomicFile error: Undefined strand at line $line\n" if not defined($strand);
		my @names = split(";", $names);
		foreach my $name (@names) {
			##$temp{exist}{0}{$name} = 1 if not defined($temp{exist}{0}{$name});
			##$temp{exist}{1}{$name} = 1 if defined($genomic{$name});
			##$temp{used}{$name} ++;
	
			push(@{$genomic{$name}}, "$chr\t$start\t$end\t$name\t0\t$strand");
		}
	}
	close $in2;
}

if ($input =~ /promoter\./) {
	open (my $inprom, "<", "/data/mitochi/Work/Project/DRIPc/bed/mm9_gencode_promoter.bed") or die "Cannot read from /data/mitochi/Work/Project/DRIPc/bed/mm9_gencode_promoter.bed: $!\n";
	while (my $line = <$inprom>) {
		chomp($line);
		my ($chr, $start, $end, $names, $zero, $strand) = split("\t", $line);
		die "DIED AT $line\n" if not defined($strand);
		print "\tWARNING: More than 1 promoter location for gene $names\n" if defined($prom{$names});
		$prom{$names} = "$chr\t$start\t$end\t$names\t0\t$strand\n";
	}
	close $inprom;
}

elsif ($input =~ /antisense\./) {
	open (my $inanti, "<", "/data/mitochi/Work/Project/DRIPc/bed/mm9_gencode_antisense.bed") or die "Cannot read from /data/mitochi/Work/Project/DRIPc/bed/mm9_gencode_antisense.bed: $!\n";
	while (my $line = <$inanti>) {
		chomp($line);
		my ($chr, $start, $end, $names, $zero, $strand) = split("\t", $line);
		die "DIED AT $line\n" if not defined($strand);
		print "\tWARNING: More than 1 antisense location for gene $names\n" if defined($anti{$names});
		$anti{$names} = "$chr\t$start\t$end\t$names\t0\t$strand\n";
	}
	close $inanti;
}

elsif ($input =~ /terminal\./) {
	open (my $interm, "<", "/data/mitochi/Work/Project/DRIPc/bed/mm9_gencode_terminal.bed") or die "Cannot read from /data/mitochi/Work/Project/DRIPc/bed/mm9_gencode_terminal.bed: $!\n";
	while (my $line = <$interm>) {
		chomp($line);
		my ($chr, $start, $end, $names, $zero, $strand) = split("\t", $line);
		die "DIED AT $line\n" if not defined($strand);
		print "\tWARNING: More than 1 terminal location for gene $names\n" if defined($term{$names});
		$term{$names} = "$chr\t$start\t$end\t$names\t0\t$strand\n";
	}
	close $interm;
}
####################################################
# DRIP ONLY: GET STRAND
my %stranddata;
open (my $in1000, "<", "../../bed/mm9_gencode.bed") or die "Cannot read from ../../bed/mm9_gencode.bed: $!\n";
while (my $line = <$in1000>) {
        chomp($line);
        next if $line =~ /#/;
        my ($chr, $start, $end, $name, $zero, $strand) = split("\t", $line);
        $stranddata{$name} = $strand;
}
close $in1000;

my %rnaseqs;
open (my $in00, "<", $RNAFile) or die "Cannot read from $RNAFile: $!\n";
while (my $line = <$in00>) {
        chomp($line);
        next if $line =~ /#/;
        my ($name, $val) = split("\t", $line);
        $rnaseqs{$name} = $val;
}
close $in00;

# 5. Parse DRIPc data and sort by RNA-seq
print BLUE "D. Parsing DRIPc File $input Coordinates\n";
open (my $in3, "<", $input) or die "Cannot read from $input: $!\n";
my %gene;
while (my $line = <$in3>) {
	chomp($line);
	next if $line =~ /#/;

	# DRIPc Peaks
	my ($chr, $start, $end, $name, $value, $strand, $genes) = split("\t", $line);

	# Since a genomic region might contain overlapping genes, use gene with highest RNA to be used to find its twin
	my @names = split(";", $genes);
	my $rna = 0;
	my $gene;
	$gene = $genes if $input =~ /intergenic/;

	# 5a. Non-promoter/terminal/antisense/intergenic
	if ($input !~ /promoter\./ and $input !~ /terminal\./ and $input !~ /antisense\./) {
		my ($chrGen, $startGen, $endGen) = (0,0,0);
		foreach my $genename (@names) {
			my $currrna = defined($rna{$genename}) ? $rna{$genename} : 0;
			if ($input !~ /intergenic/ and $rna <= $currrna) {
				foreach my $genomic_gene (@{$genomic{$genename}}) {
					next if not defined($genomic_gene);
					my ($chr2, $start2, $end2) = split("\t", $genomic_gene);
					next if ($start < $start2 and $end < $start2);
					next if ($start > $end2   and $end > $end);
					next unless ($start >= $start2 and $start <= $end2 and $end >= $start2 and $end <= $end2);
					$rna = $currrna;
					$gene = $genename;
				}
			}
		}
	}
	
	# 5b If it's promoter/terminal then use gene which start/end is closest to dripc peak
	else {
		my ($diffstart, $diffend) = (99999999999,99999999999);
		my ($chrGen, $startGen, $endGen) = (0,0,0);
		$rna = -1;
		foreach my $genename (@names) {
			my ($chr2, $start2, $end2) = $input =~ /promoter/ ? split("\t", $prom{$genename}) : $input =~ /antisense/ ? split("\t", $anti{$genename}) :split("\t", $term{$genename});
			my $currrna = defined($rna{$genename}) ? $rna{$genename} : 0;
			die "Died at $genename\n" if not defined($start2);
			my $currdiffstart = ($start >= $start2 and $start <= $end2  ) ? 0 : abs($start2 - $start);
			my $currdiffend   = ($end   >= $start2 and $end   <= $start2) ? 0 : abs($end - $end2);
			#print "\t$genename\t$chr2\t$start2\t$end2\tRNA $currrna\n" if $name eq "PROMOTER310";
			next if $start > $end2;
			next if $end   < $start2;
			if ($start >= $start2 and $start <= $end2 and $end >= $start2 and $end <= $end2) {
				if ($diffstart > 0 or $diffend > 0) {
					($chrGen, $startGen, $endGen) = ($chr2, $start2, $end2);
					$gene = $genename;
					$rna = $currrna;
					#print "\t\t$genename $diffstart > 0 or $diffend > 0: $genename rna $rna\n" if $name eq "PROMOTER310";
					$diffstart = $currdiffstart;
					$diffend = $currdiffend;
				}
				elsif ($rna < $currrna) {
					#print "\t\t$genename $rna < $currrna: $genename YES\n" if $name eq "PROMOTER310";
					($chrGen, $startGen, $endGen) = ($chr2, $start2, $end2);
					$gene = $genename;
					$rna = $currrna;
					$diffstart = $currdiffstart;
					$diffend = $currdiffend;
				}

			}
			elsif ($currdiffstart + $currdiffend < $diffstart + $diffend) {
				#print "\t\t$genename ELSE ($currdiffstart + $currdiffend < $diffstart + $diffend $genename YES\n" if $name eq "PROMOTER310";
				($chrGen, $startGen, $endGen) = ($chr2, $start2, $end2);
				$gene = $genename;
				$rna = $currrna;
				$diffstart = $currdiffstart;
				$diffend = $currdiffend;
			}
		}
	}
	
	print "$name ($chr $start $end) nexted: $genes does not overlap...\n" and next if (not defined($gene));

	my $beststrand = $stranddata{$gene} if $input !~ /intergenic/;
	die "Died at gene $gene line $line\n" if not defined($beststrand);
	$line = "$chr\t$start\t$end\t$name\t$value\t$beststrand\t$genes" if $input !~ /intergenic/;
	$gene{$name}{gene}    = $gene; # Gene name to be used to find RNA twin
	$gene{$name}{rna}     = $rna; # RNA-seq of the gene
	$gene{$name}{genomic} = $line; # DRIP peak coordinate
}
####################################################

# 6. Shuffle
# This is separate from above because we need to sort by rna-seq for cache purpose otherwise it'll be helllla slow
print BLUE "E. Shuffling $input\n";
my $minValue = $input =~ /E14/ ? 2.5 : $input =~ /3T3/ ? 5 : die "Cannot process other than E14 or 3T3\n";
open (my $out, ">", "$dripcName\_Shuffled.bed") or die "Cannot write to $dripcName\_Shuffled.bed: $!\n";
if ($input =~ /intergenic/) {
	# Get intergenic genomic coordinates as twin
	my @twin;
	foreach my $gene (keys %genomic) {
		push(@twin, $gene);
	}
	@twin = shuffle(@twin);
	# Shuffle
	my $count_Done = 0;
	foreach my $dripName (keys %gene) {
		my $total = (keys %gene);
		printf "Intergenic Done = $count_Done / $total (%.2f %%)\n", int($count_Done / $total * 10000) / 100 if $count_Done % 200 == 0;
		my $iter = 0;
		my ($chr, $start, $end, $name, $value, $strand, $genes) = split("\t", $gene{$dripName}{genomic});

		# Next if dripc value is too low (less than 5)
		next if $value < $minValue;
		for (my $i = 0; $i < 1000; $i++) {
			last if not defined($twin[$i]);
			my $randGene = $twin[$i];
			my @genomic     = @{$genomic{$randGene}};
			my $genomic;
			#print "$i $dripName $randGene\n";
			if (not defined($randGene) or not defined($genomic[0])) {
				$i --;
				#print "\t$i\t$dripName\t$chr $start $end\ttwin $randGene\t$genomic\tUNDEF i -1 = $i\n";
				next;
			}
			my $randomRegion = int(rand(@genomic));
			chomp($genomic[$randomRegion]);
			my ($chr2, $start2, $end2, $name2, $zero2, $strand2) = split("\t", $genomic[$randomRegion]);

			# If the region is intergenic, find huge region since small region might be close to genes
			if ($end2 - $start2 < 10000) {
				$i--;
				next;
			}

			die "COOR $genomic CHR $chr2 START $start2 END $end2 NAME $name2 ZERO $zero2 STRAND $strand2\n" if not defined($strand2);
			if ($end2 - $start2 < $end - $start) {
				$i --;
				my $dist = $end - $start;
				my $dist2 = $end2 - $start2;
				#print "\t$i\t$dist > $dist2: $dripName\t$chr $start $end\ttwin $randGene\t$genomic[$randomRegion]\tTOOSMALL i -1 = $i\n";
				next;
			}
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
			my $rand = rand();
			$strand2 = $rand < 0.5 ? "+" : "-";
			print $out "$chr2\t$randStart\t$randEnd\t$name\t$value\t$strand2\tDRIP=$chr,$start,$end,$strand;ORIG=$dripName,$chr,$start,$end;TWIN=$randGene,$chr2,$start2,$end2\n";
		}
		
		$count_Done ++;
	}
}
elsif ($input =~ /promoter\.name/ or $input =~ /terminal\.name/ or $input =~ /antisense\.name/) {
	my %data;
	my $last_index = "INIT";
	my $count_notused_lowdripc = 0;
	my $count_notused_lowtwin = 0;
	my $twin_notused = 0;
	my %twin_notused;
	my ($countFixed, $totalFixed) = (0,0);
	open (my $outLog, ">", "$input\_LOG.txt") or die;
	foreach my $dripName (sort {$gene{$a}{rna} <=> $gene{$b}{rna}} keys %gene) {
		my ($chr, $start, $end, $genename, $value, $strand, $genes) = split("\t", $gene{$dripName}{genomic});
		$count_notused_lowdripc ++ and next if $value < $minValue;
		

		my $rna = $gene{$dripName}{rna};
		my $gene = $gene{$dripName}{gene};
		my $length = $end - $start;

		# Get genomic coordinates (real promoter)
		my ($chrGen, $startGen, $endGen) = $input =~ /promoter/ ? split("\t", $prom{$gene}) : $input =~ /antisense/ ? split("\t", $anti{$gene}) :split("\t", $term{$gene});
		my $startpos = $start - $startGen;# < 0 ? 0 : $start - $startGen;
		my $endpos   = $endGen - $end;# < 0 ? 0 : $endGen - $end;
		# Define Index and get the cache
		my $index = 0;
		die "Undefined rna of $dripName of gene $gene\n" if not defined($rna);
		if ($rna >= 10) {
			$index = get_index($rna);
		}
		
		if ($index ne $last_index) {
			if ($index == 0) {
				next;
				print "Doing index $index\n";
				my $data = $cache->get("E14_RNAzero") if $input =~ /E14/;
				$data = $cache->get("3T3_RNAzero") if $input =~ /3T3/;
				@{$data{zero}} = @{$data};
			}
			else {
				print "Doing index $index\n";
				my $data = $cache->get("E14_RNAtwin.$index") and print "E14 data\n" if $input =~ /E14/;
				$data = $cache->get("3T3_RNAtwin.$index") if $input =~ /3T3/;
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
		if (@twin < 500) {
			$count_notused_lowtwin ++;
			next;
		}
		my $iter = 0;
		my $twin_iter = 0;
		@twin = shuffle(@twin);
		for (my $i = 0; $i < 1000; $i++) {
			my $randGene = $twin[$i]; last if not defined($randGene);
			my $genomic     = $input =~ /promoter/ ? $prom{$randGene} : $input =~ /terminal/ ? $term{$randGene} : $anti{$randGene};
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
			my $randStart;
			my $randEnd;
			if ($input !~ /antisense/) {
				$randStart = $strand eq $strand2 ? $start2 + $startpos : $start2 + $endpos;
				$randEnd   = $strand eq $strand2 ? $end2   - $endpos   : $end2 - $startpos;
			}
			else {
				$randStart = $strand ne $strand2 ? $start2 + $endpos     : $start2 + $startpos;
				$randEnd   = $strand ne $strand2 ? $end2   - $startpos   : $end2 - $endpos;
			}
			my $tempstart = $strand eq "+" ? $start2 - $randStart : $endGen - $randEnd;
			my $tempend = $strand eq "+" ? $end2 - $randEnd : $startGen - $randStart;
			my $tempstart2 = $strand eq "+" ? $start2 - $randStart : $endGen - $randEnd;
			my $tempend2 = $strand eq "+" ? $end2 - $randEnd : $startGen - $randStart;
			my $diff = $end - $start;
			my $diff2 = $randEnd - $randStart;
			die "Died at $chr $start $end $dripName gene $chrGen $startGen $endGen: randstart > randedn ($randStart > $randEnd, start2 = $start2 startpos = $startpos end2 = $end2 endpos = $endpos\n" if $randEnd < $randStart;
			die "Died at: DRIP:$chr $start $end $dripName DRIPGENE: $gene $chrGen $startGen $endGen TWIN: $chr2 $randStart $randEnd TWINGENG $chr2 $start2 $end2, diff: $diff vs $diff2, STARTPOS/ENDPOS = $tempstart/$tempend vs $tempstart2/$tempend2\n" if $randEnd - $randStart != $end - $start;
			$twin_notused{$dripName} ++;
			my $newstrand = $strand2;
			#if ($input =~ /antisense/) {
			#	$newstrand = $strand2 eq "+" ? "-" : "+";
			#}
			print $out "$chr2\t$randStart\t$randEnd\t$randGene\t$value\t$newstrand\tDRIP=$chr,$start,$end,$strand;ORIG=$gene,$chrGen,$startGen,$endGen;TWIN=$randGene,$chr2,$start2,$end2\n";
			#print "$chr2\t$randStart\t$randEnd\t$genename\t$value\t$newstrand\tDRIP=$chr,$start,$end;ORIG=$gene,$chrGen,$startGen,$endGen;TWIN=$randGene,$chr2,$start2,$end2\n";
		}
		
	}
	close $in3;

	my $count_notused_lowtwin2 = 0;
	foreach my $dripName (keys %twin_notused) {
		$count_notused_lowtwin2 ++ if $twin_notused{$dripName} < 500;
	}
	print $outLog "Not Used DRIPc < $minValue: $count_notused_lowdripc, Number of Twin < 500: $count_notused_lowtwin, Shuffles less than <500: $count_notused_lowtwin2\n";
	#runbash("bedtools intersect -v -s -a $dripcName\_Shuffled.bed -b /data/mitochi/Work/Project/DRIPc/bed/dripc.bed | sort -T /data/mitochi/Work/SortTMP/ -k1,1 -k2,2n > $dripcName\_Shuffled_noDRIPc.bed");
}
else {
	my %data;
	my $last_index = "INIT";
	my $count_notused_lowdripc = 0;
	my $count_notused_lowtwin = 0;
	my $twin_notused = 0;
	my %twin_notused;
	my $count_Done = 0;
	my $total = (keys %gene);
	open (my $outLog, ">", "$input\_LOG.txt") or die;
	foreach my $dripName (sort {$gene{$a}{rna} <=> $gene{$b}{rna}} keys %gene) {
		printf "$input Done = $count_Done / $total (%.2f %%)\n", int($count_Done / $total * 10000) / 100 if $count_Done % 200 == 0;
		my ($chr, $start, $end, $name, $value, $strand, $genes) = split("\t", $gene{$dripName}{genomic});
		
		$count_notused_lowdripc ++ and next if $value < $minValue;
		my $gene = $gene{$dripName}{gene};
		my $length = $end - $start;

		my ($chrGen, $startGen, $endGen) = $genomic{$gene};
	
		# Define Index and get the cache
		my $index = 0;
		my $rna = $gene{$dripName}{rna};
		die "Undefined rna of $dripName of gene $gene\n" if not defined($rna);
		if ($rna >= 10) {
			$index = get_index($rna);
		}
		
		if ($index ne $last_index) {
			if ($index == 0) {
				next;
				#print $outLog "Doing index $index\n";
				my $data = $cache->get("E14_RNAzero") if $input =~ /E14/;
				$data = $cache->get("3T3_RNAzero") if $input =~ /3T3/;
				@{$data{zero}} = @{$data};
			}
			else {
				#print $outLog "Doing index $index\n";
				my $data = $cache->get("E14_RNAtwin.$index") if $input =~ /E14/;
				$data = $cache->get("3T3_RNAtwin.$index") if $input =~ /3T3/;
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
		if (@twin < 500) {
			$count_notused_lowtwin ++;
			next;
		}
		my $iter = 0;
		my $iter2 = 0;
		my $twin_iter = 0;
		@twin = shuffle(@twin);
		for (my $i = 0; $i < 1000; $i++) {
			my $randGene = $twin[$i]; last if not defined($twin[$i]);
			#print $outLog "$i $dripName $randGene\n";
			my $randGeneGenomic = $genomic{$randGene};
			if ($iter2 == 100) {
				$iter2 = 0;
				#print $outLog "\t$i$dripName nexted\n";
				next;
			}
			$iter2 ++;

			if (not defined($randGene) or not defined($randGeneGenomic)) {
				$i --;
				#print $outLog "\t$i\t$dripName\t$chr $start $end\ttwin $randGene\tUNDEF i -1 = $i\n";
				next;
			}
			my @genomic     = @{$genomic{$randGene}};
			my $randomRegion = int(rand(@genomic));
			chomp($genomic[$randomRegion]);
			my ($chr2, $start2, $end2, $name2, $zero2, $strand2) = split("\t", $genomic[$randomRegion]);
			# Except promoter and terminal, we can't use region that's too short
			# Too short region will cause sucky shuffling therefore ignore genomic
			# regions less than 1kb since dripc peak median are ~1kb
			if ($end2 - $start2 < 100) {
				my $dist2 = $end2 - $start2;
				#print $outLog "\t$i\t$dist2: $dripName\t$chr $start $end\ttwin $randGene\t$genomic[$randomRegion]\tLESSTHAN100 i -1 = $i\n";
				$i --;
				next;
			}

			die "COOR $genomic[$randomRegion] CHR $chr2 START $start2 END $end2 NAME $name2 ZERO $zero2 STRAND $strand2\n" if not defined($strand2);
			if ($end2 - $start2 < $end - $start) {
				$i --;
				my $dist = $end - $start;
				my $dist2 = $end2 - $start2;
				#print $outLog "\t$i\t$dist > $dist2: $dripName\t$chr $start $end\ttwin $randGene\t$genomic[$randomRegion]\tTOOSMALL i -1 = $i\n";
				next;
			}
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
				#print $outLog "\t$i\t$spacer < 1: $dripName\t$chr $start $end\ttwin $randGene\t$genomic[$randomRegion]\tTOOSMALL i -1 = $i\n";
				next;
			}

			$twin_notused{$dripName} ++;	
			my $newstrand = $strand2;
			print $out "$chr2\t$randStart\t$randEnd\t$randGene\t$value\t$newstrand\tDRIP=$chr,$start,$end,$strand;ORIG=$gene,$chr,$start,$end;TWIN=$randGene,$chr2,$start2,$end2\n";
			$iter2 = 0;
		}
		$count_Done++;
		
	}
	close $in3;

	my $count_notused_lowtwin2 = 0;
	foreach my $dripName (keys %twin_notused) {
		$count_notused_lowtwin2 ++ if $twin_notused{$dripName} < 500;
	}

	print $outLog "Not Used DRIPc < $minValue: $count_notused_lowdripc, Number of Twin < 500: $count_notused_lowtwin, Shuffles less than < 500: $count_notused_lowtwin2\n";
	#runbash("bedtools intersect -v -s -a $dripcName\_Shuffled.bed -b /data/mitochi/Work/Project/DRIPc/bed/dripc.bed | sort -T /data/mitochi/Work/SortTMP/ -k1,1 -k2,2n > $dripcName\_Shuffled_noDRIPc.bed");
}
close $out;

# Make sure +/- 2kb of the region doesn't have DRIP and DRIPc peaks
# Change to +/- 2kb of region
runbash("bedtools_bed_change.pl -x -2000 -y 2000 -i $dripcName\_Shuffled.bed -o $dripcName\_TEMP.bed");
# Filter out those with DRIPc
#runbash("bedtools intersect -v -a $dripcName\_TEMP.bed -b /data/mitochi/Work/Project/DRIPc/bed/dripc.bed > $dripcName\_TEMP2.bed");
# Filter out those with DRIP
runbash("bedtools intersect -v -a $dripcName\_TEMP.bed -b /data/mitochi/Work/Project/DRIPc/bed/E14_DRIP_Peaks_mm9.bed > $dripcName\_TEMP2.bed") if $input =~ /E14/;
runbash("bedtools intersect -v -a $dripcName\_TEMP.bed -b /data/mitochi/Work/Project/DRIPc/bed/3T3_DRIP_Peaks_mm9.bed > $dripcName\_TEMP2.bed") if $input =~ /3T3/;
# Change width back
runbash("bedtools_bed_change.pl -x 2000 -y -2000 -i $dripcName\_TEMP2.bed -o $dripcName\_TEMP.bed");
# Sort and remove temp files
runbash("sort -T /data/mitochi/Work/SortTMP/ -k1,1 -k2,2n $dripcName\_TEMP.bed > $dripcName.txt && rm $dripcName\_TEMP*.bed");
####################################################




###############
# SUBROUTINES #
###############

sub parse_twin {
	# Data Hash
	my %data;
	my @data;

	# 1. Processing NT2 RNA lses than 10
	print "\t1. Processing E14/3T3_RNAzero.rpkm\n";
	my $inputZero = "E14_RNAzero.rpkm" if $input =~ /E14/;
	$inputZero = "3T3_RNAzero.rpkm" if $input =~ /3T3/;
	open (my $in0, "<", $inputZero) or die "Cannot read from $inputZero: $!\n";
	while (my $line = <$in0>) {
		chomp($line);
		next if $line =~ /#/;
		push(@data, $line);
	}
	close $in0;
	$cache->set("E14_RNAzero", \@data) if $input =~ /E14/;
	$cache->set("3T3_RNAzero", \@data) if $input =~ /3T3/;
	@data = ();
	
	print "\t2. Processing E14/3T3Twin.rpkm\n";
	my $inputtwin = "E14_RNAtwin.rpkm" if $input =~ /E14/;
	$inputtwin = "3T3_RNAtwin.rpkm" if $input =~ /3T3/;
	open (my $in1, "<", $inputtwin) or die "Cannot read from $inputtwin: $!\n";
	my $last_index = "INIT";
	while (my $line = <$in1>) {
		chomp($line);
		next if $line =~ /#/;
		my ($gene, @twin) = split("\t", $line);
		print "ENSMUST00000165223.1 exist\n" if $gene eq "ENSMUST00000165223.1";
		# If number of gene twin is less than 20, don't use it
		my $rna = $rna{$gene};

		# Get index based on gene rpkm
		my $index = get_index($rna);
		print "\tProcessing index $index\n" if $last_index ne $index;
		if ($index ne $last_index and $last_index ne "INIT") {
			$cache->set("E14_RNAtwin.$last_index", \%data) if $input =~ /E14/;
			$cache->set("3T3_RNAtwin.$last_index", \%data) if $input =~ /3T3/;
			%data = ();
		}
		@{$data{$gene}} = @twin;
		print "ENSMUST00000165223.1 data gene $gene: exist = @twin\n" if $gene eq "ENSMUST00000165223.1";
		$last_index = $index;
	}
	$cache->set("E14_RNAtwin.$last_index", \%data) if $input =~ /E14/;
	$cache->set("3T3_RNAtwin.$last_index", \%data) if $input =~ /3T3/;
	%data = ();
	close $in1;
	
	$cache->set("E14_RNAtwin", 1) if $input =~ /E14/;
	$cache->set("3T3_RNAtwin", 1) if $input =~ /3T3/;
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

