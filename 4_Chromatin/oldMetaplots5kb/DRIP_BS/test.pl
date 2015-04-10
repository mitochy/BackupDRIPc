#!/usr/bin/perl

use strict; use warnings; use mitochy; use Getopt::Std; use Cache::FileCache;
use vars qw($opt_a $opt_b $opt_i $opt_f);
getopt("a:b:fi");

die "$opt_a B $opt_b\n";
