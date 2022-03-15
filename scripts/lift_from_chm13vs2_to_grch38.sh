#!/bin/bash

# USAGE: sh.liftover_bedfile <BED file path> <CHAIN file path grch38_onto_chm13_autosome_X_Y_wfmash-0.7.0+a36ab5f_p95_s5k_l25k.trim.chain>
#bash lift_from_chm13vs2_to_grch38.sh /lizardfs/guarracino/liftover/data/CHM13.combined.v4.unique.bed.gz grch38_onto_chm13_autosome_X_Y_wfmash-0.7.0+a36ab5f_p95_s5k_l25k.trim.chain out

# The script converts the BED file to a Picard interval list file, using
# the fourth BED file as an identifier, or creating one if it doesn't exist. Any fields
# beyond the fourth field in the input bed file will be ignored in the picard run

PATH_LIFTOVER=/home/guarracino/tools/liftOver
PATH_PICARD_JAR=/home/guarracino/tools/picard.2.26.11.jar

export OLDREFDICT=/lizardfs/guarracino/liftover/data/chm13v2.0.dict  #/data/Phillippy/t2t-share/team-liftover/liftover_genes/chm13v2.0.dict
export NEWREFDICT=/lizardfs/guarracino/liftover/data/hg38.chrom.dict #/data/Phillippy/t2t-share/team-variants/grch38_t2t_liftover/ref_files/hg38.chrom.dict

export BED=$1
export CHAINFILE=$2
export DIR_OUTPUT=$3


export BED4=`echo $BED | sed 's:.*/::' | sed 's/.bed/.bed4/'`
export UCSCOUTPUT=`echo $BED | sed 's:.*/::' | sed 's/.bed/.grch38.ucsc.bed/'`
export UCSCUNMAPPED=`echo $BED | sed 's:.*/::' | sed 's/.bed/.grch38.ucsc.unmapped.bed/'`
export UCSCSORTED=`echo $BED | sed 's:.*/::' | sed 's/.bed/.grch38.ucsc.sort.bed/'`

export PICARDINTERVAL=`echo $BED | sed 's:.*/::' | sed 's/.bed/.picard.interval_list/'`
export PICARDOUTPUT=`echo $BED | sed 's:.*/::' | sed 's/.bed/.grch38.picard.interval_list/'`
export PICARDOUTPUTBED=`echo $PICARDOUTPUT | sed 's/.interval_list/.bed/'`
export PICARDSORTED=`echo $PICARDOUTPUT | sed 's/.interval_list/.sort.bed/'`
export PICARDUNMAPPED=`echo $BED | sed 's:.*/::' | sed 's/.bed/.grch38.picard.unmapped.interval_list/'`

mkdir -p $DIR_OUTPUT

# create an "interval_list" file for Picard:
cat $OLDREFDICT > $DIR_OUTPUT/$PICARDINTERVAL
zcat $BED | awk '(NF<4 || $4=="") {OFS="\t"; print $1, $2+1, $3, "+", "INTERVAL_"NR} $4!="" {OFS="\t"; print $1, $2+1, $3, "+", $4}' | grep -v '#chrom' >> $DIR_OUTPUT/$PICARDINTERVAL

## create a bed4 file (with identifiers for regions) for liftOver:
zcat $BED | awk '(NF<4 || $4=="") {OFS="\t"; print $1, $2, $3, "INTERVAL_"NR} $4!="" {OFS="\t"; print $1, $2, $3, $4}' > $DIR_OUTPUT/$BED4

## run UCSC liftOver:
$PATH_LIFTOVER $DIR_OUTPUT/$BED4 $CHAINFILE $DIR_OUTPUT/$UCSCOUTPUT $DIR_OUTPUT/$UCSCUNMAPPED -multiple
sort -k1,1 -k2,2n -k3,3n $DIR_OUTPUT/$UCSCOUTPUT > $DIR_OUTPUT/$UCSCSORTED

# run picard LiftOverIntervalList:
java -Xmx6g -jar $PATH_PICARD_JAR LiftOverIntervalList -I $DIR_OUTPUT/$PICARDINTERVAL -O $DIR_OUTPUT/$PICARDOUTPUT --CHAIN $CHAINFILE --REJECT $DIR_OUTPUT/$PICARDUNMAPPED -SD $NEWREFDICT


grep -v '^@' $DIR_OUTPUT/$PICARDOUTPUT | awk '{OFS="\t"; print $1, $2-1, $3, $5}' > $DIR_OUTPUT/$PICARDOUTPUTBED
sort -k1,1 -k2,2n -k3,3n $DIR_OUTPUT/$PICARDOUTPUTBED > $DIR_OUTPUT/$PICARDSORTED
