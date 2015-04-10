# Get only unique reads without headers
run_script_in_paralel2.pl "grep -v -P '^\@\w\w\t' SX014C.sam | grep -v -P '\t4\t\*\t0\t' > FILENAME.out" . sam 4
# Count reads in each sample
wc -l *.sam.out
   14687670 SX014C.sam.out
   20841866 SX014D.sam.out
   17001011 SX014E.sam.out
   10664011 SX014F.sam.out
# Take lowest read and get that amount from each sample (10000000)
run_script_in_paralel2.pl "random_line.pl FILENAME 10000000 > FILENAME.final" . sam.out 4
