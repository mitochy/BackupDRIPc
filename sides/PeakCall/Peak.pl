#!/usr/bin/perl

use strict; use warnings; use mitochy;

my ($input1) = @ARGV;
die "usage: $0 <input1>\n" unless @ARGV;

my ($folder1, $fileName1) = mitochy::getFilename($input1, "folder");
my $strand = $input1 =~ /pos/i ? "+" : "-";
my $name   = $input1 =~ /pos/i ? "pos" : "neg";
my @curr;
my $curr_val = "INIT";
my $curr_chr = 0;
my $curr_dir = 0;
my $curr_start = 0;
my $curr_end = 0;
my $freeze_end;
my $freeze_start;
my $check_dir = 0;
my @val;
my $count = 0;
my $prev_mean;
my $prev_line;
open (my $in1, "<", $input1) or die "Cannot read from $input1: $!\n";
while (my $line = <$in1>) {
	chomp($line);
	next if $line =~ /#/;
	my ($chr, $start, $end, $val) = split("\t", $line);
	next if $val < 4;
	if ($curr_val eq "INIT") {
		$curr_chr   = $chr;
		$curr_start = $start;
		$curr_end   = $end;
		$curr_val   = $val;
		$curr_dir   = 0;
		push(@val, $val);
	}
	elsif ($chr ne $curr_chr) {
		my $mean = mean(@val); $count++;
		my $color = getcolor($mean);
		my $line = "$curr_chr\t$curr_start\t$curr_end\t1$name\_$count\t$mean\t$strand\t$curr_start\t$curr_end\t$color";
		print "$line\n";
		$check_dir = 0;
		@val = ();
		$curr_chr   = $chr;
		$curr_start = $start;
		$curr_end   = $end;
		$curr_val   = $val;
		$curr_dir   = 0;
		push(@val, $val);
	}
	elsif ($start != $curr_end + 1) {
		my $mean = mean(@val); $count++;
		my $color = getcolor($mean);
		my $line = "$curr_chr\t$curr_start\t$curr_end\t2$name\_$count\t$mean\t$strand\t$curr_start\t$curr_end\t$color";
		if (defined($prev_mean) and $prev_mean < 5) {
			my ($chr2, $start2, $end2, $name2, $val2, $strand2) = split("\t", $prev_line);
			if ($mean < 5 and $end2 > 0 and $curr_start - $end2 < 300) {
				$count -= 1;
				$prev_mean = mean($prev_mean, $mean);
				$color = getcolor($prev_mean);
				$prev_line = "$chr2\t$start2\t$curr_end\tprev1$name2\t$prev_mean\t$strand\t$start2\t$curr_end\t$color";
			}
			else {
				print "$prev_line\n" if $end2 > 0;
				$check_dir = 0;
				if ($mean < 5) {
					$prev_mean = $mean;
					$prev_line = $line;
				}
				else {
					print "$line\n";
					$check_dir = 0;
					$prev_mean = 0;
					$prev_line = "chr\t-1000\t-1000";
				}
			}
		}
		else {
			if ($mean < 5) {
				$prev_mean = $mean;
				$prev_line = $line;
			}
			else {
				$prev_mean = 0;
				$prev_line = "chr\t-1000\t-1000";
				$check_dir = 0;
				print "$line\n";
			}
		}
		@val = ();
		$curr_chr   = $chr;
		$curr_start = $start;
		$curr_end   = $end;
		$curr_val   = $val;
		$curr_dir   = 0;
		push(@val, $val);
	}
	elsif ($curr_dir == 0) {
		$curr_end   = $end;
		$curr_dir   = $val - $curr_val > 0 ? 1 : -1;
		$curr_val   = $val;		
		$curr_chr   = $chr;
		push(@val, $val);
	}
	elsif ($curr_dir != 0) {
		if ($check_dir == 2) {
			$check_dir = 0;
			my $mean = mean(@val); $count++;
			my $color = getcolor($mean);
			my $line = "$curr_chr\t$freeze_start\t$freeze_end\t3$name\_$count\t$mean\t$strand\t$freeze_start\t$freeze_end\t$color";
			if (defined($prev_mean) and $prev_mean < 5) {
				my ($chr2, $start2, $end2, $name2, $val2, $strand2) = split("\t", $prev_line);
				if ($mean < 5 and $end2 > 0 and $curr_start - $end2 < 300) {
					$count -= 1;
					$prev_mean = mean($prev_mean, $mean);
					$color = getcolor($prev_mean);
					$prev_line = "$chr2\t$start2\t$freeze_end\tprev2$name2\t$prev_mean\t$strand\t$start2\t$freeze_end\t$color";
				}
				else {
					print "$prev_line\n" if $end2 > 0;
					if ($mean < 5) {
						$prev_mean = $mean;
						$prev_line = $line;
					}
					else {
						$prev_mean = 0;
						$prev_line = "chr\t-1000\t-1000";
						print "$line\n";
					}
				}
			}
			else {
				if ($mean < 5) {
					$prev_mean = $mean;
					$prev_line = $line;
				}
				else {
					$prev_mean = 0;
					$prev_line = "chr\t-1000\t-1000";
					print "$line\n";
				}
			}
			$curr_start = $start;
			$curr_end   = $end;
			@val = ();
			push(@val, $val);
		}
		my $dir     = $val - $curr_val > 0 ? 1 : -1;
		if ($dir != $curr_dir) {
			if (($dir == -1 and $val/$curr_val < 0.5) or ($dir == 1 and $curr_val/$val < 0.5)) {
				my $mean = mean(@val); $count++;
				my $color = getcolor($mean);
				my $line = "$curr_chr\t$curr_start\t$curr_end\t4$name\_$count\t$mean\t$strand\t$curr_start\t$curr_end\t$color";
				if (defined($prev_mean) and $prev_mean < 5) {
					my ($chr2, $start2, $end2, $name2, $val2, $strand2) = split("\t", $prev_line);
					if ($mean < 5 and $end2 > 0 and $curr_start - $end2 < 300) {
						$count -= 1;
						$prev_mean = mean($prev_mean, $mean);
						$color = getcolor($prev_mean);
						$prev_line = "$chr2\t$start2\t$curr_end\tprev3$name2\t$prev_mean\t$strand\t$start2\t$curr_end\t$color";
					}
					else {
						print "$prev_line\n" if $end2 > 0;
						if ($mean < 5) {
							$prev_mean = $mean;
							$prev_line = $line;
						}
						else {
							$prev_mean = 0;
							$prev_line = "chr\t-1000\t-1000";
							print "$line\n";
						}
					}
				}
				else {
					if ($mean < 5) {
						$prev_mean = $mean;
						$prev_line = $line;
					}
					else {
						$prev_mean = 0;
						$prev_line = "chr\t-1000\t-1000";
						print "$line\n";
					}
				}
				@val = ();
				$curr_chr   = $chr;
				$curr_start = $start;
				$curr_end   = $end;
				$curr_val   = $val;
				$curr_dir   = 0;
				$check_dir = 0;
				push(@val, $val);
			}
			else {
				$freeze_start = $curr_start if $check_dir == 0;
				$freeze_end = $curr_end if $check_dir == 0;
				$check_dir ++;
				push(@val, $val);
			}
		}
		elsif ($dir == $curr_dir) {
			undef $freeze_end;
			$check_dir = 0;
			push(@val, $val);
		}
		$curr_end = $end;
	}
}
close $in1;

sub mean {
	my @val = @_;
	my $mean = 0;
	foreach my $val (@val) {
		$mean += $val / @val;
	}
	return($mean);

}

sub getcolor {
	my ($val) = @_;
	return("150,0,13")    if $val > 15 and $strand eq "+";
	return("8,69,148")    if $val > 15 and $strand eq "-";
	return("203,24,29")   if $val <= 15 and $val > 10 and $strand eq "+";
	return("33,113,181")  if $val <= 15 and $val > 10 and $strand eq "-";
	return("239,59,44")   if $val <= 10 and $val > 7 and $strand eq "+";
	return("66,146,198")  if $val <= 10 and $val > 7 and $strand eq "-";
	return("251,106,74")  if $val <= 7 and $strand eq "+";
	return("107,174,214") if $val <= 7 and $strand eq "-";

}
__END__
open (my $out1, ">", "$fileName1.out") or die "Cannot write to $fileName1.out: $!\n";
close $out1;
