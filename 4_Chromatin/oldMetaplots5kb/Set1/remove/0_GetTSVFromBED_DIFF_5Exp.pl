#!/usr/bin/perl
# dripc_promoter.shuffled is in /data/mitochi/Work/Project/DRIPc/4_Chromatin/1_Shuffle/Result_125/dripc_promoter.shuffled
use strict; use warnings; use mitochy; use Getopt::Std; use R_toolbox;
use vars qw($opt_r $opt_i $opt_a $opt_b);
getopts("r:i:a:b:");

die "
usage: $0 -r <RNAseq> -i <BED dripc_promoter.shuf> -a <orig.tsv e.g. MCF_groseq_promoter_orig.tsv> -b <shuf.tsv>

E.g. :

./0_GetTSVFromBED.pl -r ../../../data/NT2.rpkm -a MCF_groseq_promoter_orig.tsv -b MCF_groseq_promoter_shuf.tsv
./0_GetTSVFromBED.pl -r ../../../data/NT2.rpkm -a MCF_groseq_terminal_orig.tsv -b MCF_groseq_terminal_shuf.tsv
./0_GetTSVFromBED.pl -r ../../../data/NT2.rpkm -a MCF_groseq_genebody_orig.tsv -b MCF_groseq_genebody_shuf.tsv
./0_GetTSVFromBED.pl -r ../../../data/NT2.rpkm -a MCF_groseq_antisense_orig.tsv -b MCF_groseq_antisense_shuf.tsv


" unless defined($opt_r) and defined($opt_a) and defined($opt_b);
my ($rnaFile, $origFile, $shufFile) = ($opt_r, $opt_a, $opt_b);
#my ($bedFile, $rnaFile, $origFile, $shufFile) = ($opt_i, $opt_r, $opt_a, $opt_b);
my ($folder1, $fileNameOrig) = mitochy::getFilename($origFile, "folder");
my ($folder2, $fileNameShuf) = mitochy::getFilename($shufFile, "folder");

# 0. Parse NT2.rpkm
print "1. Parsing RNA file $rnaFile\n";
my %rna = %{parse_rna($rnaFile)};

#my %data = %{parse_peak($bedFile)};
# 2. Divide orig tsv file into into expression quartiles
open (my $outsuper1,  ">", "$fileNameOrig\_super.temp")  or die "Cannot write to $fileNameOrig\_super.temp: $!\n";
open (my $outhigh1,   ">", "$fileNameOrig\_high.temp")  or die "Cannot write to $fileNameOrig\_high.temp: $!\n";
open (my $outmed1,    ">", "$fileNameOrig\_med.temp")  or die "Cannot write to $fileNameOrig\_med.temp: $!\n";
open (my $outlow1,    ">", "$fileNameOrig\_low.temp")  or die "Cannot write to $fileNameOrig\_low.temp: $!\n";
open (my $outzero1,   ">", "$fileNameOrig\_zero.temp")  or die "Cannot write to $fileNameOrig\_zero.temp: $!\n";

print "3. Dividing orig tsv $origFile\n";
open (my $in, "<", $origFile) or die "Cannot read from $origFile: $!\n";
while (my $line = <$in>) {
	chomp($line);
	next if $line =~ /#/;
	my ($name) = split("\t", $line);
	#if ($name !~ /^ENS/) {
	#	$name = $data{name}{$name};
	#}
	my $rna = defined($rna{$name}) ? $rna{$name} : 0;
	#next if not defined($data{orig}{$name}{count});
	#for (my $i = 0; $i < $data{orig}{$name}{count}; $i++) {
		$line =~ s/\tNA/\t0/g;
		print $outsuper1  "$line\n" if $rna > 200;
		print $outhigh1   "$line\n" if $rna <= 200 and $rna > 100;
		print $outmed1    "$line\n" if $rna >= 50 and $rna < 100;
		print $outlow1    "$line\n" if $rna >= 10 and $rna < 50;
		print $outzero1   "$line\n" if $rna < 10;
		
	#}
	#$data{orig}{$name}{count} = 0;
}
close $in;

print "4. Dividing Shuffled tsv $shufFile\n";
# 2. Divide shuffled tsv file into into expression quartiles
open (my $outsuper2,  ">", "$fileNameShuf\_super.temp")  or die "Cannot write to $fileNameShuf\_super.temp: $!\n";
open (my $outhigh2,   ">", "$fileNameShuf\_high.temp")  or die "Cannot write to $fileNameShuf\_high.temp: $!\n";
open (my $outmed2,    ">", "$fileNameShuf\_med.temp")  or die "Cannot write to $fileNameShuf\_med.temp: $!\n";
open (my $outlow2,    ">", "$fileNameShuf\_low.temp")  or die "Cannot write to $fileNameShuf\_low.temp: $!\n";
open (my $outzero2,   ">", "$fileNameShuf\_zero.temp")  or die "Cannot write to $fileNameShuf\_zero.temp: $!\n";

open (my $in2, "<", $shufFile) or die "Cannot read from $shufFile: $!\n";
while (my $line = <$in2>) {
	chomp($line);
	next if $line =~ /#/;
	my ($name, @val) = split("\t", $line);
	if ($shufFile =~ /antisense/) {
		@val = reverse(@val);
	}
	my $val = join("\t", @val);
	$line = "$name\t$val";
	#if ($name !~ /^ENS/) {
	#	$name = $data{name}{$name};
	#}
	#else {
		#next if not defined($data{shuf}{$name});
		#foreach my $origName (keys %{$data{shuf}{$name}}) {
		#	for (my $i = 0; $i < $data{shuf}{$name}{$origName}; $i++) {
		#		$line =~ s/\tNA/\t0/g;
				#my $rna = defined($rna{$origName}) ? $rna{$origName} : 0;
				my $rna = defined($rna{$name}) ? $rna{$name} : 0;
				print $outsuper2  "$line\n" if $rna > 200;
				print $outhigh2   "$line\n" if $rna <= 200 and $rna > 100;
				print $outmed2    "$line\n" if $rna >= 50 and $rna < 100;
				print $outlow2    "$line\n" if $rna >= 10 and $rna < 50;
				print $outzero2   "$line\n" if $rna < 10;
		#	}
		#	$data{shuf}{$name}{$origName} = 0;
		#}
	#}
}
close $in2;

# Run Graph.pl
runbash("./Graph_DIFF_5Exp.pl $fileNameOrig\_high.temp");

# run Rscript
my $RscriptName = $fileNameShuf;
$RscriptName =~ s/_shuf//;
runbash("run_Rscript.pl $RscriptName.R > $RscriptName.Rlog 2>&1");

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

sub parse_peak {
	my ($input) = @_;
	my %data;
	my %used;
	# Parse the bedfile
	open (my $in1, "<", $input) or die "Cannot read from $input: $!\n";
	while (my $line = <$in1>) {
		chomp($line);
		next if $line =~ /#/;
		my ($chr, $start, $end, $name, $val, $strand, $info) = split("\t", $line);
		my ($orig, $shuf) = $info =~ /ORIG=(\w+\.\d+),chr.+TWIN=(\w+\.\d+),chr/;
		die "Undefined twin gene at line: $line\n" if not defined($shuf);
		push(@{$data{orig}{$orig}{temp}}, $shuf);
		next if defined($used{$name}) and $used{$name} == 1;
		$data{orig}{$orig}{count} ++;
		$data{name}{$name} = $orig;
		$used{$name} = 1;
	}
	
	close $in1;

	# randomize array and take first 10
	foreach my $name (keys %{$data{orig}}) {
		my @value = shuffle(@{$data{orig}{$name}{temp}});
		@{$data{orig}{$name}{temp}} = ();
		for (my $i = 0; $i < 10; $i++) {
			$data{shuf}{$value[$i]}{$name} ++;
		}
		@value = ();
	}
	return(\%data);
}

sub shuffle {
        my (@value) = @_;
	my $shuffleTimes = @value < 1000 ? 1000 : @value;
        for (my $i = 0; $i < $shuffleTimes; $i++) {
                my $rand1 = int(rand(@value));
                my $rand2 = int(rand(@value));
                my $val1 = $value[$rand1];
                my $val2 = $value[$rand2];
                $value[$rand1] = $val2;
                $value[$rand2] = $val1;
        }
        return(@value);
}

sub runbash {
        my ($cmd) = @_;
        print "\t$cmd\n";
        system($cmd) == 0 or die "Failed to run $cmd: $!\n";
}

__END__

		#print $outzero1  "$line\n" if $rna == 0;
		#print $outlow1   "$line\n" if $rna > 10 and $rna <= 50;
		#print $outmed1   "$line\n" if $rna > 50 and $rna <= 100;
		#print $outhigh1  "$line\n" if $rna > 100 and $rna <= 200;
		#print $outsuper1 "$line\n" if $rna > 200;


./0_GetTSVFromBED.pl -r ../../../data/NT2.rpkm -i ../../1_Shuffle/Result_125/dripc_promoter.shuffled -a MCF_groseq_promoter_orig.tsv -b MCF_groseq_promoter_shuf.tsv
./0_GetTSVFromBED.pl -r ../../../data/NT2.rpkm -i ../../1_Shuffle/Result_125/dripc_terminal.shuffled -a MCF_groseq_terminal_orig.tsv -b MCF_groseq_terminal_shuf.tsv
./0_GetTSVFromBED.pl -r ../../../data/NT2.rpkm -i ../../1_Shuffle/Result_125/dripc_genebody.shuffled -a MCF_groseq_genebody_orig.tsv -b MCF_groseq_genebody_shuf.tsv
./0_GetTSVFromBED.pl -r ../../../data/NT2.rpkm -i ../../1_Shuffle/Result_125/dripc_antisense.shuffled -a MCF_groseq_antisense_orig.tsv -b MCF_groseq_antisense_shuf.tsv

