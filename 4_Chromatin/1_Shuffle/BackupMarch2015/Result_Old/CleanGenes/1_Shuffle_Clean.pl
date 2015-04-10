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
die "usage: $0 <../../1_Location/dripc_*.bed>\n" unless @ARGV == 1;
my ($folder, $dripcName) = mitochy::getFilename($input, "folder");

my $shuffleNumber = 1000;
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
my $genomicFile = "/data/mitochi/Work/Project/DRIPc/1_Location/4_CreateRegion_Clean/genomic_$feature.bed"; die "Died undefined genomic $genomicFile for $dripcName\n" unless -e $genomicFile;
runbash("bedtools intersect -s -wb -a $input -b $genomicFile | cut -f1-6,13 > $dripcName.name");

# DEPRECATED: Antisense_other intersect with antisense of all genic region (which is genic minus antisense and promoter (terminal* genebody both*))
#if ($input =~ /antisense_other/) {
#	$genomicFile = "/data/mitochi/Work/Project/DRIPc/1_Location/genomic_antisense_other.bed";
#	runbash("bedtools intersect -s -wb -a $input -b $genomicFile | cut -f1-6,13 > $dripcName.name");
#}
$input = "$dripcName.name";
####################################################

####################################################
# 4. Parse Genomic Location, including hg19_gencode_promoter and terminal
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
	open (my $inprom, "<", "/data/mitochi/Work/Project/DRIPc/bed/hg19_gencode_promoter_clean.bed") or die "Cannot read from /data/mitochi/Work/Project/DRIPc/bed/hg19_gencode_promoter_clean.bed: $!\n";
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
	open (my $inanti, "<", "/data/mitochi/Work/Project/DRIPc/bed/hg19_gencode_antisense.bed") or die "Cannot read from /data/mitochi/Work/Project/DRIPc/bed/hg19_gencode_antisense.bed: $!\n";
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
	open (my $interm, "<", "/data/mitochi/Work/Project/DRIPc/bed/hg19_gencode_terminal_clean.bed") or die "Cannot read from /data/mitochi/Work/Project/DRIPc/bed/hg19_gencode_terminal_clean.bed: $!\n";
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
			if ($input !~ /intergenic/ and $rna <= $rna{$genename}) {
				foreach my $genomic_gene (@{$genomic{$genename}}) {
					next if not defined($genomic_gene);
					my ($chr2, $start2, $end2) = split("\t", $genomic_gene);
					next if ($start < $start2 and $end < $start2);
					next if ($start > $end2   and $end > $end);
					next unless ($start >= $start2 and $start <= $end2 and $end >= $start2 and $end <= $end2);
					$rna = $rna{$genename};
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
			my $currrna = $rna{$genename};
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

	$gene{$name}{gene}    = $gene; # Gene name to be used to find RNA twin
	$gene{$name}{rna}     = $rna; # RNA-seq of the gene
	$gene{$name}{genomic} = $line; # DRIP peak coordinate
}
####################################################

# 6. Shuffle
# This is separate from above because we need to sort by rna-seq for cache purpose otherwise it'll be helllla slow
print BLUE "E. Shuffling $input\n";
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

		# Next if dripc value is too low (less than 5)
		next if $value < 5;
		for (my $i = 0; $i < $shuffleNumber; $i++) {
			my $randGene = $twin[int(rand(@twin))];
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
			print $out "$chr2\t$randStart\t$randEnd\t$name\t$value\t$strand2\tDRIP=$chr,$start,$end;ORIG=$dripName,$chr,$start,$end;TWIN=$randGene,$chr2,$start2,$end2\n";
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
	foreach my $dripName (sort {$gene{$a}{rna} <=> $gene{$b}{rna}} keys %gene) {
		my ($chr, $start, $end, $genename, $value, $strand, $genes) = split("\t", $gene{$dripName}{genomic});
		$count_notused_lowdripc ++ and next if $value < 5;
		

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
		$count_notused_lowtwin ++ and next if @twin < 50;
		my $iter = 0;
		my $twin_iter = 0;
		for (my $i = 0; $i < $shuffleNumber; $i++) {
			my $randGene = $twin[int(rand(@twin))];
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
			print $out "$chr2\t$randStart\t$randEnd\t$genename\t$value\t$newstrand\tDRIP=$chr,$start,$end;ORIG=$gene,$chrGen,$startGen,$endGen;TWIN=$randGene,$chr2,$start2,$end2\n";
			#print "$chr2\t$randStart\t$randEnd\t$genename\t$value\t$newstrand\tDRIP=$chr,$start,$end;ORIG=$gene,$chrGen,$startGen,$endGen;TWIN=$randGene,$chr2,$start2,$end2\n";
		}
		
	}
	close $in3;

	my $count_notused_lowtwin2 = 0;
	foreach my $dripName (keys %twin_notused) {
		$count_notused_lowtwin2 ++ if $twin_notused{$dripName} < 50;
	}
	print "Not Used DRIPc < 5: $count_notused_lowdripc, Number of Twin < 50: $count_notused_lowtwin, Twin < 50: $count_notused_lowtwin2\n";
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
	#open (my $outLog, ">", "$input\_LOG.txt") or die;
	foreach my $dripName (sort {$gene{$a}{rna} <=> $gene{$b}{rna}} keys %gene) {
		printf "$input Done = $count_Done / $total (%.2f %%)\n", int($count_Done / $total * 10000) / 100 if $count_Done % 200 == 0;
		my ($chr, $start, $end, $name, $value, $strand, $genes) = split("\t", $gene{$dripName}{genomic});
		
		$count_notused_lowdripc ++ and next if $value < 5;
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
				#print $outLog "Doing index $index\n";
				my $data = $cache->get("NT2_RNAzero");
				@{$data{zero}} = @{$data};
			}
			else {
				#print $outLog "Doing index $index\n";
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
		$count_notused_lowtwin ++ and next if @twin < 50;
		my $iter = 0;
		my $iter2 = 0;
		my $twin_iter = 0;
		for (my $i = 0; $i < $shuffleNumber; $i++) {
			my $randGene = $twin[int(rand(@twin))];
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
			print $out "$chr2\t$randStart\t$randEnd\t$name\t$value\t$newstrand\tDRIP=$chr,$start,$end;ORIG=$gene,$chr,$start,$end;TWIN=$randGene,$chr2,$start2,$end2\n";
			$iter2 = 0;
		}
		$count_Done++;
		
	}
	close $in3;

	my $count_notused_lowtwin2 = 0;
	foreach my $dripName (keys %twin_notused) {
		$count_notused_lowtwin2 ++ if $twin_notused{$dripName} < 50;
	}

	#print $outLog "Not Used DRIPc < 5: $count_notused_lowdripc, Number of Twin < 50: $count_notused_lowtwin, Twin < 50: $count_notused_lowtwin2\n";
	#runbash("bedtools intersect -v -s -a $dripcName\_Shuffled.bed -b /data/mitochi/Work/Project/DRIPc/bed/dripc.bed | sort -T /data/mitochi/Work/SortTMP/ -k1,1 -k2,2n > $dripcName\_Shuffled_noDRIPc.bed");
}
close $out;

# Make sure +/- 2kb of the region doesn't have DRIP and DRIPc peaks
# Change to +/- 2kb of region
runbash("bedtools_bed_change.pl -x -2000 -y 2000 -i $dripcName\_Shuffled.bed -o $dripcName\_TEMP.bed");
# Filter out those with DRIPc
runbash("bedtools intersect -v -a $dripcName\_TEMP.bed -b /data/mitochi/Work/Project/DRIPc/bed/dripc.bed > $dripcName\_TEMP2.bed");
# Filter out those with DRIP (Merged 1 and 2)
runbash("bedtools intersect -v -a $dripcName\_TEMP2.bed -b /data/mitochi/Work/Project/DRIPc/bed/DRIP_Intersect_Peak.bed > $dripcName\_TEMP.bed");
# Change width back
runbash("bedtools_bed_change.pl -x 2000 -y -2000 -i $dripcName\_TEMP.bed -o $dripcName\_TEMP2.bed");
# Sort and remove temp files
runbash("sort -T /data/mitochi/Work/SortTMP/ -k1,1 -k2,2n $dripcName\_TEMP2.bed > $dripcName.txt && rm $dripcName\_TEMP*.bed");

=comment
# 7. Below is to filter shuffled peaks that has too high DRIPc (DRIP is not done coz it's intersect)
# a. Get positive and negative strands and run DRIPc
my $outputName = "$dripcName\_Shuffled.bed";
runbash("grep + $outputName | sort -T /data/mitochi/Work/SortTMP/ -k1,1 -k2,2n > $outputName.POS");
runbash("grep - $outputName | sort -T /data/mitochi/Work/SortTMP/ -k1,1 -k2,2n > $outputName.NEG");
# b. Run map_wig_to_bed_BIG.pl to DRIPc wig file to filter 
runbash("map_wig_to_bed_BIG.pl -w /data/mitochi/Work/Project/DRIPc/wig/NT2_DRIPc_pos.wig -m -r /data/mitochi/Work/Cache/ $outputName.POS");
runbash("map_wig_to_bed_BIG.pl -w /data/mitochi/Work/Project/DRIPc/wig/NT2_DRIPc_neg.wig -m -r /data/mitochi/Work/Cache/ $outputName.NEG");
# c. Cat the result togetehr
runbash("cat NT2_DRIPc_pos_$dripcName\_Shuffled.txt NT2_DRIPc_neg_$dripcName\_Shuffled.txt > $dripcName.txt");
runbash("rm NT2_DRIPc_pos_$dripcName\_Shuffled.txt && rm NT2_DRIPc_neg_$dripcName\_Shuffled.txt");
runbash("rm $outputName.POS $outputName.NEG $dripcName.name");
=cut
print "\n\n$input DONE!\n$input output = $dripcName.txt\n\n";
####################################################




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
