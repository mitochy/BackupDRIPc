#!/usr/bin/sh

grep -P "\texon\t" hg19_gencode.gtf | cut -f1,4,5,7,9 | perl -pi -e 's/gene_id \"(.+)\"\; transcript_id.+$/$1/' | perl -pi -e 's/^(.+)\t(.+)\t(.+)\t(.+)\t(.+)$/$1\t$2\t$3\t$5\t0\t$4/' > hg19_gencode_exon.bed
bedtools merge -s -nms -i hg19_gencode_exon.bed > hg19_gencode_exon_merged.bed
./1_MergeGene.pl hg19_gencode_exon_merged.bed
bed2gtf hg19_gencode_exon_merged.out
