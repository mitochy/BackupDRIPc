#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($input) = @ARGV;
die "usage: $0 <input>\n" unless @ARGV;

my ($folder, $filename) = mitochy::getFilename($input, "folder");
system("echo \"df = numeric(0)\" > $input.TEMP");
system("cat $input >> $input.TEMP");
system("mv $input.TEMP $input");
system("echo \"library(ggplot2);pdf(\'$filename.pdf\');ggplot(df,aes(histone,val)) + geom_boxplot(aes(fill=histone),outlier.shape=NA) + coord_cartesian(ylim=c(0,2.5));dev.off()\" >> $input");
