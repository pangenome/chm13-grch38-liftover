#Notes to set up this on octopus
#guix install openjdk
#wget -c https://github.com/broadinstitute/picard/releases/download/2.26.11/picard.jar && mv picard.jar picard.2.26.11.jar && chmod +x picard.2.26.11.jar


# USAGE: bash liftover_evaluation.sh grch38_onto_chm13_autosome_X_Y_wfmash-0.7.0+a36ab5f_p95_s5k_l25k.trim.chain out



CHAIN_FILE_GRCH38_ONTO_CHM13=$1
DIR_OUTPUT=$2

bash lift_from_chm13vs2_to_grch38.sh /lizardfs/guarracino/liftover/data/CHM13.combined.v4.unique.bed $CHAIN_FILE_GRCH38_ONTO_CHM13 $DIR_OUTPUT 2> /dev/null
cd $DIR_OUTPUT

# intersect gencode lifted chm13v2 locations with CAT/liftoff predicted locations:
bedtools intersect -f 1.0 -F 1.0 -b CHM13.combined.v4.unique.grch38.picard.sort.bed -a /lizardfs/guarracino/liftover/data/gencodeV35.bed.gz -u > gencodeV35.primary.unique.ucsc.inCHM13.combined.v4.unique.grch38.sort.bed
bedtools intersect -f 1.0 -F 1.0 -a CHM13.combined.v4.unique.grch38.picard.sort.bed -b /lizardfs/guarracino/liftover/data/gencodeV35.primary.unique.bed.gz -u > CHM13.combined.v4.unique.grch38.ucsc.ingencodeV35.bed

# count number of CHM13.combined.v4 (CAT/liftoff) intervals that exactly match gencodeV35 annotations when lifted to GRCh38:
#
# Original count of CHM13.combined.v4 intervals (in chm13v1.1)
n=$(zgrep -v '#' /lizardfs/guarracino/liftover/data/CHM13.combined.v4.unique.bed.gz -c)
echo -e "Total intervals in CHM13.combined.v4: $n"


# Count of intervals that were unmappable
echo "Unmappable intervals according to UCSC liftOver:"
grep '#' CHM13.combined.v4.unique.grch38.ucsc.unmapped.bed | sort | uniq -c

# Count of intervals that were successfully lifted
n=$(wc -l CHM13.combined.v4.unique.grch38.ucsc.sort.bed | cut -f 1 -d ' ')
echo -e "Intervals successfully created in chm13v2.0: $n"

# Count of intervals that lifted to intervals that match the (GRCh38) gencodeV35 intervals
n=$(wc -l CHM13.combined.v4.unique.grch38.ucsc.ingencodeV35.bed | cut -f 1 -d ' ')
echo -e "Lifted intervals which match the (chm13) CAT/liftoff intervals: $n"
