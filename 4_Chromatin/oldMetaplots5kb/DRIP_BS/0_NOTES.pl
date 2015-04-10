# Note: Only processing promoter DRIP on ORIGINAL ONLY (This doesn't use SHUFFLED coz we're comparing DRIPBS and whole BS)
# 1. Go to DRIPc/Shuffle to shuffle 125x
./1_Shuffle.pl -r ../data/NT2.rpkm -z ../data/NT2_RNAzero.rpkm -t ../data/NT2_RNAtwin.rpkm -g ../bed/hg19_gencode_promoter.bed -i ../1_Location/4_CreateRegion_DRIP/Region/NT2_DRIPBS_Peak_region_promoter.bed
./2_PostProcess.pl NT2_DRIPBS_Peak_region_promoter_Shuffled.bed
Output: NT2_DRIPBS_Peak_region_promoter_Shuffled.shuffled

# 2. Use the shuffled and get .orig and .shuf regions and change to +/- 5kb of TSS or peak
#	a. Take input.shuffled (not used but it can reduce the number of shuffled peaks into whatever instead of 125)
#	b. For promoter/terminal/antisense can be centered on TSS or peak itself
#	c. Print out into original (orig.input) and shuffled (shuf.input)
#	d. Change coordinate into +/- 5kb of peak or TSS (Use TSS for this)
./0_GetRegions.pl ../../../Shuffle/result/NT2_DRIPBS_Peak_region_promoter_Shuffled.shuffled
rm *.tmp

# 3. Run map_wig_to_bed_BIG.pl -w BSMETH_INTERSECT.wig -p input_orig.input
map_wig_to_bed_BIG.pl  -w ../../../wig/BSMETH_intersect.wig -p NT2_DRIPBS_Peak_region_promoter_Shuffled_orig.input 
map_wig_to_bed_BIG.pl  -w ../../../wig/DRIPBSMETH_intersect.wig -p NT2_DRIPBS_Peak_region_promoter_Shuffled_orig.input 

# 4. To get diffmeth: run 2_DRIPBS_GetHighUnmeth.pl

# 5. To graph:
# Run calculate_heatmap_meth.pl -w 200 -s 100 input_orig.txt
# 0b_GetTSVFromBED_meth_5Exp.pl: Divide into 5 clusters for meth only. The diff is meth doesn't turn NA into 0
# create figure for promoter/terminal etc categorized by expression levels: 
# - Super: RNA >= 200
# - High : RNA >= 100 & $NA < 200
# - Med  : RNA >= 50 & RNA < 100
# - Low  : RNA >= 10 & RNA < 50
# - Zero : RNA < 10
# Output pdf
