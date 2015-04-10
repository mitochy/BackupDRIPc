#!/usr/bin/perl

use strict; use warnings; use mitochy; use Thread; use Thread::Queue;

my @histone = qw(
CTCF
DNASE
DNASE1
Faire
H3K27ac
H3K27me3
H3K36me3
H3K4me1
H3K4me3
H3K79me2
H3K9me3
H4K20me1
MNase
NT2_DRIPc
K562_PROseq
);

my @feature = qw(
promoter
promoter_ext
terminal
terminal_ext
antisense
antisense_ext
intergenic
antisense_other
both_ext
both
genebody
);

my $Q = new Thread::Queue;
for (my $i = 0; $i < @histone; $i++) {
	my $histone = $histone[$i];
	for (my $j = 0; $j < @feature; $j++) {
		my $feature = $feature[$j];
		my $orig = "$histone\_dripc\_$feature\_orig.tsv";
		my $shuf = "$histone\_dripc\_$feature\_shuf.tsv";
		print "PDF EXIST: $histone $feature\n" and next if -e "$histone\_dripc\_$feature.pdf";
		print "NOT MADE: $histone $feature not found\n" and next if not -e $orig and not -e $shuf;
		if (not -e $orig or not -e $shuf) {
			print "MISSING: $orig not found\n" if not -e $orig;
			print "MISSING: $shuf not found\n" if not -e $shuf;
			next;
		}
		my $command = "./0_GetTSVFromBED_5Exp.pl -r ../../../data/NT2.rpkm -a $orig -b $shuf";
		$Q->enqueue($command);
	}
}
$Q->end();
my @threads;

for (my $i = 0; $i < 5; $i++) {
	$threads[$i] = threads->create(\&worker, $i, $Q);
}
for (my $i = 0; $i < 5; $i++) {
	$threads[$i]->join();
}

sub worker {
        my ($thread, $queue) = @_;
        my $tid = threads->tid;

        while ($queue->pending) {
                my $command = $queue->dequeue;
                print "processing $command in thread $tid\n";
                `$command`;
        }
    return;
}
