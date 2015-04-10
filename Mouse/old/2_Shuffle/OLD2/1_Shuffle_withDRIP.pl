#!/usr/bin/perl
# This script find random genes of each dripc peak
# 1. Parse RNA-seq ($RNAInput)
# 2. Parse and Cache Twin Gene Expression Set ($zeroInput and $twinInput)
#	Index by expression value
#	$twinInput.$index
# 3. Parse Genomic Files
# 4. Parse DRIPc Files
# 5. For each DRIPc peak, get RNA seq 
# 6. Original: $regionName.bed, Output: $regionName.shuffled
# Run map wig to original and output
use strict; use warnings; use mitochy; use Cache::FileCache; use Getopt::Std;
use Term::ANSIColor qw(:constants);
$Term::ANSIColor::AUTORESET = 1;
use vars qw($opt_i $opt_r $opt_z $opt_t $opt_g $opt_d $opt_w);
getopt("i:r:z:t:g:d:w:");

my ($regionInput, $RNAInput, $zeroInput, $twinInput, $genomeInput, $dripInputWig, $dripInput) = ($opt_i, $opt_r, $opt_z, $opt_t, $opt_g, $opt_w, $opt_d);
die "usage: $0 -r <RNA.rpkm> -z <RNA_RNAzero.id> -t <RNA_RNAtwin.id> -g <Region_Genomic.bed> -w <DRIP Wig> -d <DRIP Peak> -i <Region.bed>\n" unless defined($opt_i) and defined($opt_r) and defined($opt_z) and defined($opt_t) and defined($opt_g) and defined($opt_d);
my ($folder1, $zeroInputName) = mitochy::getFilename($zeroInput, "folder");
my ($folder2, $twinInputName) = mitochy::getFilename($twinInput, "folder");
my ($folder, $regionName)      = mitochy::getFilename($regionInput, "folder");

my $shuffleNumber = 1000;
####################################################
# 1. Parse RNA-seq
print BLUE "A. Parsing RNA-seq $RNAInput\n";
my %rna;
open (my $inRNA, "<", "$RNAInput") or die "Cannot read from $RNAInput: $!\n";
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
print BLUE "B. Parsing and Caching Twin Gene Exprsesion Set $twinInputName and $zeroInputName\n";
my $cache = new Cache::FileCache; $cache->set_cache_root("/data/mitochi/Work/Cache/");
my $data  = $cache->get("$twinInputName");
parse_twin() if (not defined($data)); # If not exist, cache the twin gene
####################################################

####################################################
# 3. Put gene name to column 7 of each DRIPc peak except intergenic by intersecting with its genomic file
runbash("bedtools intersect -wb -a $regionInput -b $genomeInput | cut -f1-6,13 > $regionName.name");
$regionInput = "$regionName.name";
####################################################

####################################################
# 4. Parse Genomic Location, including mm10_gencode_promoter and terminal
print BLUE "C. Parsing Genomic File $genomeInput Coordinates\n";
my %genomic;
my %prom;
my %term;
my %anti;
##my %temp;
open (my $in2, "<", $genomeInput) or die "Cannot read from $genomeInput: $!\n";
while (my $line = <$in2>) {
	chomp($line);
	my ($chr, $start, $end, $names, $zero, $strand) = split("\t", $line);
	die "Parsing genomcFile $genomeInput error: Undefined strand at line $line\n" if not defined($strand);
	my @names = split(";", $names);
	foreach my $name (@names) {
		push(@{$genomic{$name}}, "$chr\t$start\t$end\t$name\t0\t$strand");
	}
}
close $in2;

open (my $inprom, "<", $genomeInput) or die "Cannot read from $genomeInput: $!\n";
while (my $line = <$inprom>) {
	chomp($line);
	my ($chr, $start, $end, $names, $zero, $strand) = split("\t", $line);
	die "DIED AT $line\n" if not defined($strand);
	print "\tWARNING: More than 1 promoter location for gene $names\n" if defined($prom{$names});
	$prom{$names} = "$chr\t$start\t$end\t$names\t0\t$strand\n";
}
close $inprom;
open (my $interm, "<", $genomeInput) or die "Cannot read from $genomeInput: $!\n";
while (my $line = <$interm>) {
	chomp($line);
	my ($chr, $start, $end, $names, $zero, $strand) = split("\t", $line);
	die "DIED AT $line\n" if not defined($strand);
	print "\tWARNING: More than 1 terminal location for gene $names\n" if defined($term{$names});
	$term{$names} = "$chr\t$start\t$end\t$names\t0\t$strand\n";
}
close $interm;
####################################################

# 5. Parse DRIPc data and sort by RNA-seq
print BLUE "D. Parsing DRIPc File $regionInput Coordinates\n";
open (my $in3, "<", $regionInput) or die "Cannot read from $regionInput: $!\n";
my %gene;
while (my $line = <$in3>) {
	chomp($line);
	next if $line =~ /#/;

	# DRIPc Peaks
	my ($chr, $start, $end, $name, $value, $strand, $genes) = split("\t", $line);

	print "Processing $name ($chr $start $end)\n" if $name eq "PROMOTER310";
	# Since a genomic region might contain overlappnig genes, use gene with highest RNA to be used to find its twin
	my @names = split(";", $genes);
	my $rna = 0;
	my $gene;
	$gene = $genes if $regionInput =~ /intergenic/;
	if ($regionInput !~ /promoter\./ and $regionInput !~ /terminal\./ and $regionInput !~ /antisense\./) {
		my ($chrGen, $startGen, $endGen) = (0,0,0);
		if (@names == 1) {
			$rna = $rna{$names[0]};
			$gene = $names[0];
		}
		else {
			foreach my $genename (@names) {
				if ($regionInput !~ /intergenic/ and $rna <= $rna{$genename}) {
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
	}
	
	# If it's promoter/terminal then use gene which start/end is closest to dripc peak
	else {
		my ($diffstart, $diffend) = (99999999999,99999999999);
		my ($chrGen, $startGen, $endGen) = (0,0,0);
		$rna = -1;
		foreach my $genename (@names) {
			my ($chr2, $start2, $end2) = $regionInput =~ /promoter/ ? split("\t", $prom{$genename}) : $regionInput =~ /antisense/ ? split("\t", $anti{$genename}) : split("\t", $term{$genename});
			my $currrna = $rna{$genename};
			die "Died at $genename\n" if not defined($start2);
			my $currdiffstart = ($start >= $start2 and $start <= $end2  ) ? 0 : abs($start2 - $start);
			my $currdiffend   = ($end   >= $start2 and $end   <= $start2) ? 0 : abs($end - $end2);
			print "\t$genename\t$chr2\t$start2\t$end2\tRNA $currrna\n" if $name eq "PROMOTER310";
			next if $start > $end2;
			next if $end   < $start2;
			if ($start >= $start2 and $start <= $end2 and $end >= $start2 and $end <= $end2) {
				if ($diffstart > 0 or $diffend > 0) {
					($chrGen, $startGen, $endGen) = ($chr2, $start2, $end2);
					$gene = $genename;
					$rna = $currrna;
					print "\t\t$genename $diffstart > 0 or $diffend > 0: $genename rna $rna\n" if $name eq "PROMOTER310";
					$diffstart = $currdiffstart;
					$diffend = $currdiffend;
				}
				elsif ($rna < $currrna) {
					print "\t\t$genename $rna < $currrna: $genename YES\n" if $name eq "PROMOTER310";
					($chrGen, $startGen, $endGen) = ($chr2, $start2, $end2);
					$gene = $genename;
					$rna = $currrna;
					$diffstart = $currdiffstart;
					$diffend = $currdiffend;
				}

			}
			elsif ($currdiffstart + $currdiffend < $diffstart + $diffend) {
				print "\t\t$genename ELSE ($currdiffstart + $currdiffend < $diffstart + $diffend $genename YES\n" if $name eq "PROMOTER310";
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
print BLUE "E. Shuffling $regionInput\n";
open (my $out, ">", "$regionName\_Shuffled.bed") or die "Cannot write to $regionName\_Shuffled.bed: $!\n";
if ($regionInput =~ /intergenic/) {

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
	runbash("bedtools intersect -v -a $regionName\_Shuffled.bed -b $dripInput | sort -T /data/mitochi/Work/SortTMP/ -k1,1 -k2,2n > $regionName\_Shuffled_noDRIP.bed");
}
elsif ($regionInput =~ /promoter\./ or $regionInput =~ /terminal\./ or $regionInput =~ /antisense\./) {
	my %data;
	my $last_index = "INIT";
	my $count_notused_lowdripc = 0;
	my $count_notused_lowtwin = 0;
	my $twin_notused = 0;
	my %twin_notused;
	my ($countFixed, $totalFixed) = (0,0);
	foreach my $dripName (sort {$gene{$a}{rna} <=> $gene{$b}{rna}} keys %gene) {
		my ($chr, $start, $end, $genename, $value, $strand, $genes) = split("\t", $gene{$dripName}{genomic});
		$count_notused_lowdripc ++ and next if $value < 5 and $regionInput =~ /drip/i;
		

		my $rna = $gene{$dripName}{rna};
		my $gene = $gene{$dripName}{gene};
		my $length = $end - $start;

		# Get genomic coordinates (real promoter)
		my ($chrGen, $startGen, $endGen) = $regionInput =~ /promoter/ ? split("\t", $prom{$gene}) : $regionInput =~ /antisense/ ? split("\t", $anti{$gene}) :split("\t", $term{$gene});
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
				my $data = $cache->get("$zeroInputName");
				@{$data{zero}} = @{$data};
			}
			else {
				print "Doing index $index\n";
				my $data = $cache->get("$twinInputName.$index");
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
		for (my $i = 0; $i < 1000; $i++) {
			my $randGene = $twin[int(rand(@twin))];
			my $genomic     = $regionInput =~ /promoter/ ? $prom{$randGene} : $term{$randGene};
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
			my $randStart = $strand eq $strand2 ? $start2 + $startpos : $start2 + $endpos;
			my $randEnd   = $strand eq $strand2 ? $end2   - $endpos   : $end2 - $startpos;
			my $tempstart = $strand eq "+" ? $start2 - $randStart : $endGen - $randEnd;
			my $tempend = $strand eq "+" ? $end2 - $randEnd : $startGen - $randStart;
			my $tempstart2 = $strand eq "+" ? $start2 - $randStart : $endGen - $randEnd;
			my $tempend2 = $strand eq "+" ? $end2 - $randEnd : $startGen - $randStart;
			my $diff = $end - $start;
			my $diff2 = $randEnd - $randStart;
			die "Died at $chr $start $end $dripName gene $chrGen $startGen $endGen: randstart > randedn ($randStart > $randEnd, start2 = $start2 startpos = $startpos end2 = $end2 endpos = $endpos\n" if $randEnd < $randStart;
			die "Died at: DRIP:$chr $start $end $dripName DRIPGENE: $gene $chrGen $startGen $endGen TWIN: $chr2 $randStart $randEnd TWINGENG $chr2 $start2 $end2, diff: $diff vs $diff2, STARTPOS/ENDPOS = $tempstart/$tempend vs $tempstart2/$tempend2\n" if $randEnd - $randStart != $end - $start;
			$twin_notused{$dripName} ++;
			print $out "$chr2\t$randStart\t$randEnd\t$genename\t$value\t$strand2\tDRIP=$chr,$start,$end;ORIG=$gene,$chrGen,$startGen,$endGen;TWIN=$randGene,$chr2,$start2,$end2\n";
		}
		
	}
	close $in3;

	my $count_notused_lowtwin2 = 0;
	foreach my $dripName (keys %twin_notused) {
		$count_notused_lowtwin2 ++ if $twin_notused{$dripName} < 50;
	}
	print "Not Used DRIPc < 5: $count_notused_lowdripc, Number of Twin < 50: $count_notused_lowtwin, Twin < 50: $count_notused_lowtwin2\n";
	runbash("bedtools intersect -v -a $regionName\_Shuffled.bed -b $dripInput | sort -T /data/mitochi/Work/SortTMP/ -k1,1 -k2,2n > $regionName\_Shuffled_noDRIP.bed");
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
	open (my $outLog, ">", "$regionInput\_LOG.txt") or die;
	foreach my $dripName (sort {$gene{$a}{rna} <=> $gene{$b}{rna}} keys %gene) {
		printf "$regionInput Done = $count_Done / $total (%.2f %%)\n", int($count_Done / $total * 10000) / 100 if $count_Done % 200 == 0;
		my ($chr, $start, $end, $name, $value, $strand, $genes) = split("\t", $gene{$dripName}{genomic});
		
		$count_notused_lowdripc ++ and next if $value < 5 and $regionInput =~ /drip/i;
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
				print $outLog "Doing index $index\n";
				my $data = $cache->get("$zeroInputName");
				@{$data{zero}} = @{$data};
			}
			else {
				print $outLog "Doing index $index\n";
				my $data = $cache->get("$twinInputName.$index");
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
		if (@twin < 50) {
			$count_notused_lowtwin ++;
		#	print "$gene (RNA $rna)\n";
			next;
		}
		my $iter = 0;
		my $iter2 = 0;
		my $twin_iter = 0;
		for (my $i = 0; $i < $shuffleNumber; $i++) {
			my $randGene = $twin[int(rand(@twin))];
			print $outLog "$i $dripName $randGene\n";
			my $randGeneGenomic = $genomic{$randGene};
			if ($iter2 == 100) {
				$iter2 = 0;
				print $outLog "\t$i$dripName nexted\n";
				next;
			}
			$iter2 ++;

			if (not defined($randGene) or not defined($randGeneGenomic)) {
				$i --;
				print $outLog "\t$i\t$dripName\t$chr $start $end\ttwin $randGene\tUNDEF i -1 = $i\n";
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
				print $outLog "\t$i\t$dist2: $dripName\t$chr $start $end\ttwin $randGene\t$genomic[$randomRegion]\tLESSTHAN100 i -1 = $i\n";
				$i --;
				next;
			}

			die "COOR $genomic[$randomRegion] CHR $chr2 START $start2 END $end2 NAME $name2 ZERO $zero2 STRAND $strand2\n" if not defined($strand2);
			if ($end2 - $start2 < $end - $start) {
				$i --;
				my $dist = $end - $start;
				my $dist2 = $end2 - $start2;
				print $outLog "\t$i\t$dist > $dist2: $dripName\t$chr $start $end\ttwin $randGene\t$genomic[$randomRegion]\tTOOSMALL i -1 = $i\n";
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
				print $outLog "\t$i\t$spacer < 1: $dripName\t$chr $start $end\ttwin $randGene\t$genomic[$randomRegion]\tTOOSMALL i -1 = $i\n";
				next;
			}

			$twin_notused{$dripName} ++;	
			print $out "$chr2\t$randStart\t$randEnd\t$name\t$value\t$strand2\tDRIP=$chr,$start,$end;ORIG=$gene,$chr,$start,$end;TWIN=$randGene,$chr2,$start2,$end2\n";
			$iter2 = 0;
		}
		$count_Done++;
		
	}
	close $in3;

	my $count_notused_lowtwin2 = 0;
	foreach my $dripName (keys %twin_notused) {
		$count_notused_lowtwin2 ++ if $twin_notused{$dripName} < 100;
	}

	print $outLog "Not Used DRIPc < 5: $count_notused_lowdripc, Number of Twin < 50: $count_notused_lowtwin, Shuffled Less than 100: $count_notused_lowtwin2\n";
	runbash("bedtools intersect -v -a $regionName\_Shuffled.bed -b $dripInput | sort -T /data/mitochi/Work/SortTMP/ -k1,1 -k2,2n > $regionName\_Shuffled_noDRIP.bed");
}
close $out;

runbash("map_wig_to_bed_BIG.pl -w $dripInputWig -m -r /data/mitochi/Work/Cache/ $regionName\_Shuffled_noDRIP.bed");
# b. Run map_wig_to_bed

###############
# SUBROUTINES #
###############

sub parse_twin {
	# Data Hash
	my %data;
	my @data;

	# 1. Processing NT2 RNA lses than 10
	print "\t1. Processing $zeroInput\n";
	my $regionInputZero = "$zeroInput";
	open (my $in0, "<", $regionInputZero) or die "Cannot read from $regionInputZero: $!\n";
	while (my $line = <$in0>) {
		chomp($line);
		next if $line =~ /#/;
		push(@data, $line);
	}
	close $in0;
	$cache->set("$zeroInput", \@data);
	@data = ();
	
	print "\t2. Processing NT2Twin.rpkm\n";
	my $regionInputtwin = "$twinInput";
	open (my $in1, "<", $regionInputtwin) or die "Cannot read from $regionInputtwin: $!\n";
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
			$cache->set("$twinInputName.$last_index", \%data);
			%data = ();
		}
		@{$data{$gene}} = @twin;
		$last_index = $index;
	}
	$cache->set("$twinInputName.$last_index", \%data);
	%data = ();
	close $in1;
	
	$cache->set("$twinInputName", 1);
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
