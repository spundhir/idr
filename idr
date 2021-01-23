#!/bin/bash
#PBS -l nodes=1:ppn=4

## initialize variables with default values
MACS_P_VALUE="1e-3"
IDR_THRESHOLD="0.01"
PROGDIR="/localhome/bric/xfd783/software/idrCode/"
GENOME="mm"

#### usage ####
usage() {
	echo Program: "idr (run IDR analysis for ChIP-seq data)"
	echo Author: BRIC, University of Copenhagen, Denmark
	echo Version: 1.0
	echo Contact: pundhir@binf.ku.dk
    echo "Usage: idr -i <files> -c <file> -o <dir> [OPTIONS]"
	echo "Options:"
    echo " -i <file>   [mapped tag (sample) files in BAM format (separated by comma)]"
    echo "             [format: <identical file name>_Rep[1|2].bam]"
    echo " -c <file>   [mapped tag (control) files in BAM format (separated by comma)]"
    echo "             [format: <identical file name>_Rep[1|2].bam]"
    echo "             [both control and real samples should be in same directory]"
    echo " -o <dir>    [output directory (****should be ABSOLUTE path****)]"
    echo "[OPTIONS]"
    echo " -p <dir>    [path to dependent R scripts (default: ~/software/idrCode)]"
    echo " -t <float>  [IDR threshold (default: 0.01)]"
    echo " -g <string> [effective genome size, required by macs2 (default: mm)]"
    echo "             [availble: hs|mm|ce|dm]"
	echo " -h          [help]"
	echo
	exit 0
}

## function to compute the factorial of a number
factorial(){
    fact=1
    n=$1
    while [ $n -ge 1 ]; do
        fact=`expr $fact \* $n`
        n=`expr $n - 1`
    done
    echo "$fact"
}

#### parse options ####
while getopts i:c:o:p:t:g:h ARG; do
	case "$ARG" in
        i) CHIP_INPUT=$OPTARG;;
        c) CONTROL_INPUT=$OPTARG;;
        o) OUTDIR=$OPTARG;;
        p) PROGDIR=$OPTARG;;
        t) IDR_THRESHOLD=$OPTARG;;
        g) GENOME=$OPTARG;;
		h) HELP=1;;
	esac
done

## usage, if necessary file and directories and not given/exist
if [ ! "$CHIP_INPUT" -o ! "$CONTROL_INPUT" -o ! "$OUTDIR" -o "$HELP" ]; then
    usage
fi

## create appropriate directory structure
echo
echo -n "Create appropriate directory structure.. "
if [ ! -d "$OUTDIR/homer" ]; then
    mkdir -p $OUTDIR/homer
    mkdir -p $OUTDIR/tagAlign
    mkdir -p $OUTDIR/macs
    mkdir -p $OUTDIR/quality
    mkdir -p $OUTDIR/logs
    mkdir -p $OUTDIR/consistency
    mkdir -p $OUTDIR/consistency/reps
    mkdir -p $OUTDIR/consistency/selfPseudoReps
    mkdir -p $OUTDIR/consistency/pooledPseudoReps
fi
echo "done"

###############################################
## determine number of input ChIP samples
IFS=","
TMP=($CHIP_INPUT)
CHIP_COUNT=${#TMP[@]}
IFS=" "

## determine common id of the input ChIP samples
CHIP_ID=`echo $TMP[0] | sed 's/^.*\///g' | sed 's/Rep.*//g'`

## determine directory of the input ChIP/control sample files
INDIR=`echo $TMP[0] | perl -ane 'if($_=~/\//) { $_=~s/\/[^\/]+$//g; print $_; } else { print "."; }'`

echo "Total ChIP samples: $CHIP_COUNT ($CHIP_ID)"

## determine number of input control samples
IFS=","
TMP=($CONTROL_INPUT)
CONTROL_COUNT=${#TMP[@]}
IFS=" "

## determine common id of the input control samples
CONTROL_ID=`echo $TMP[0] | sed 's/^.*\///g' | sed 's/Rep.*//g'`

echo "Total control samples: $CONTROL_COUNT ($CONTROL_ID)"

<<"COMMENT1"
COMMENT1
## convert BAM into tagAlign format for ChIP samples
echo -n "Convert ChIP and control samples from BAM to tagAlign format.. "
for (( i=1; i<=$CHIP_COUNT; i++ )); do
    samtools view -b -F 1548 -q 30 $INDIR/$CHIP_ID"Rep"$i.bam | bamToBed -i stdin | awk 'BEGIN{S="\t";OFS="\t"}{$4="N";print $0}' | gzip -c > $OUTDIR/tagAlign/$CHIP_ID"Rep"$i.tagAlign.gz
done

## convert BAM into tagAlign format for control samples
for (( i=1; i<=$CONTROL_COUNT; i++ )); do
    samtools view -b -F 1548 -q 30 $INDIR/$CONTROL_ID"Rep"$i.bam | bamToBed -i stdin | awk 'BEGIN{S="\t";OFS="\t"}{$4="N";print $0}' | gzip -c > $OUTDIR/tagAlign/$CONTROL_ID"Rep"$i.tagAlign.gz
done
echo "done"

## pool all control sample files into one file
echo -n "Pool all control samples into one.. "
COMMAND="" 
for (( i=1; i<=$CONTROL_COUNT; i++ )); do
    COMMAND="$COMMAND $OUTDIR/tagAlign/$CONTROL_ID"Rep"$i.tagAlign.gz";
done

gunzip -c $COMMAND | gzip -c > $OUTDIR/tagAlign/$CONTROL_ID"Rep0.tagAlign.gz"
echo "done"

##########################################################################
############ CALL PEAKS ON INDIVIDUAL REPLICATES
##########################################################################

## determine fragment length for all ChIP samples using macs
echo -n "Determine fragment length of each ChIP sample.. "
FRAGMENT_LENGTH=(0)
for (( i=1; i<=$CHIP_COUNT; i++ )); do
    macs2 predictd -i $OUTDIR/tagAlign/$CHIP_ID"Rep"$i.tagAlign.gz -g $GENOME --outdir $OUTDIR/logs 2> $OUTDIR/logs/$CHIP_ID"Rep"$i.predictd
    FRAGMENT_LENGTH[$i]=`grep "predicted fragment length" $OUTDIR/logs/$CHIP_ID"Rep"$i.predictd | perl -ane 'print $F[scalar(@F)-2];'`
done
echo "done"

## peak calling on all the ChIP samples using macs
echo -n "Call peaks on each ChIP sample.. "
for (( i=1; i<=$CHIP_COUNT; i++ )); do
    #SHIFT_SIZE=$(echo "scale=0; ${FRAGMENT_LENGTH[$i]}/2" | bc)
    SHIFT_SIZE=`echo $(( FRAGMENT_LENGTH[$i] ))`
    macs2 callpeak -t $OUTDIR/tagAlign/$CHIP_ID"Rep"$i.tagAlign.gz -c $OUTDIR/tagAlign/$CONTROL_ID"Rep0".tagAlign.gz -f BED -n $OUTDIR/macs/$CHIP_ID"Rep"$i"_Vs_"$CONTROL_ID"Rep0" -g $GENOME -p $MACS_P_VALUE --nomodel --extsize $SHIFT_SIZE 2>$OUTDIR/logs/$CHIP_ID"Rep"$i"_Vs_"$CONTROL_ID"Rep0.log"
done
echo "done"

## retrieve top 100,000 peaks sorted by p-value
echo -n "Retrieve top 100,000 peaks sorted by p-value.. "
for (( i=1; i<=$CHIP_COUNT; i++ )); do
    sort -k 8nr,8nr $OUTDIR/macs/$CHIP_ID"Rep"$i"_Vs_"$CONTROL_ID"Rep0_peaks.narrowPeak" | head -n 100000 | gzip -c > $OUTDIR/macs/$CHIP_ID"Rep"$i"_Vs_"$CONTROL_ID"Rep0.regionPeak.gz"
done
echo "done"

##########################################################################
############ CALL PEAKS ON POOLED REPLICATES
##########################################################################

## pool all ChIP sample files into one file
echo -n "Pool all ChIP samples into one file.. "
COMMAND=""
for (( i=1; i<=$CHIP_COUNT; i++ )); do
    COMMAND="$COMMAND $OUTDIR/tagAlign/$CHIP_ID""Rep"$i".tagAlign.gz"
done

gunzip -c $COMMAND | gzip -c > $OUTDIR/tagAlign/$CHIP_ID"Rep0.tagAlign.gz"
echo "done"

## determine fragment length for pooled ChIP sample using macs
echo -n "Determine fragment length for pooled ChIP sample.. "
macs2 predictd -i $OUTDIR/tagAlign/$CHIP_ID"Rep0.tagAlign.gz" -g $GENOME --outdir $OUTDIR/logs 2> $OUTDIR/logs/$CHIP_ID"Rep0.predictd"
FRAGMENT_LENGTH=`grep "predicted fragment length" $OUTDIR/logs/$CHIP_ID"Rep0.predictd" | perl -ane 'print $F[scalar(@F)-2];'`
echo "done"

## peak calling on pooled ChIP sample using macs
echo -n "Call peaks on pooled ChIP sample.. "
SHIFT_SIZE=`echo $(( FRAGMENT_LENGTH ))`
macs2 callpeak -t $OUTDIR/tagAlign/$CHIP_ID"Rep0.tagAlign.gz" -c $OUTDIR/tagAlign/$CONTROL_ID"Rep0.tagAlign.gz" -f BED -n $OUTDIR/macs/$CHIP_ID"Rep0_Vs_"$CONTROL_ID"Rep0" -g $GENOME -p $MACS_P_VALUE --nomodel --extsize $SHIFT_SIZE 2>$OUTDIR/logs/$CHIP_ID"Rep0_Vs_"$CONTROL_ID"Rep0.log"
echo "done"

## retrieve top 100,000 peaks sorted by p-value
echo -n "Retrieve top 100,000 peaks sorted by p-value.. "
sort -k 8nr,8nr $OUTDIR/macs/$CHIP_ID"Rep0_Vs_"$CONTROL_ID"Rep0_peaks.narrowPeak" | head -n 100000 | gzip -c > $OUTDIR/macs/$CHIP_ID"Rep0_Vs_"$CONTROL_ID"Rep0.regionPeak.gz"
echo "done"

## determine fragment length for pooled ChIP sample using macs (control)
echo -n "Determine fragment length for pooled ChIP sample (control).. "
macs2 predictd -i $OUTDIR/tagAlign/$CONTROL_ID"Rep0.tagAlign.gz" -g $GENOME --outdir $OUTDIR/logs 2> $OUTDIR/logs/$CONTROL_ID"Rep0.predictd"
FRAGMENT_LENGTH=`grep "predicted fragment length" $OUTDIR/logs/$CONTROL_ID"Rep0.predictd" | perl -ane 'print $F[scalar(@F)-2];'`
echo "done"

## peak calling on pooled ChIP sample using macs (control)
echo -n "Call peaks on pooled ChIP sample (control).. "
SHIFT_SIZE=`echo $(( FRAGMENT_LENGTH ))`
macs2 callpeak -t $OUTDIR/tagAlign/$CONTROL_ID"Rep0.tagAlign.gz" -c $OUTDIR/tagAlign/$CHIP_ID"Rep0.tagAlign.gz" -f BED -n $OUTDIR/macs/$CONTROL_ID"Rep0_Vs_"$CHIP_ID"Rep0" -g $GENOME -p $MACS_P_VALUE --nomodel --extsize $SHIFT_SIZE 2>$OUTDIR/logs/$CONTROL_ID"Rep0_Vs_"$CHIP_ID"Rep0.log"
echo "done"

##########################################################################
############ CALL PEAKS ON PSEUDOREPLICATES OF INDIVIDUAL CHIP REPLICATES
##########################################################################

## randomly split input ChIP samples into two parts
echo -n "Randomly split ChIP samples into two parts.. "
for (( i=1; i<=$CHIP_COUNT; i++ )); do
    filename=$OUTDIR/tagAlign/$CHIP_ID"Rep"$i".tagAlign.gz"
    nlines=$( gunzip -c $filename | wc -l );
    nlines=$(( (nlines + 1) / 2 ))
    gunzip -c ${filename} | shuf | split -l ${nlines} - $OUTDIR/tagAlign/$CHIP_ID"Rep"$i
    gzip $OUTDIR/tagAlign/$CHIP_ID"Rep"$i"aa"
    gzip $OUTDIR/tagAlign/$CHIP_ID"Rep"$i"ab"
    mv $OUTDIR/tagAlign/$CHIP_ID"Rep"$i"aa.gz" $OUTDIR/tagAlign/$CHIP_ID"Rep"$i".pr1.tagAlign.gz"
    mv $OUTDIR/tagAlign/$CHIP_ID"Rep"$i"ab.gz" $OUTDIR/tagAlign/$CHIP_ID"Rep"$i".pr2.tagAlign.gz"
done
echo "done"

## determine fragment length for all randomly split ChIP samples using macs
echo -n "Determine fragment length for all randomly split ChIP samples.. "
FRAGMENT_LENGTH=(0)
j=0
for (( i=1; i<=$CHIP_COUNT; i++ )); do
    j=$(( j + 1 ))
    macs2 predictd -i $OUTDIR/tagAlign/$CHIP_ID"Rep"$i.pr1.tagAlign.gz -g $GENOME --outdir $OUTDIR/logs 2> $OUTDIR/logs/$CHIP_ID"Rep"$i.pr1.predictd
    FRAGMENT_LENGTH[$j]=`grep "predicted fragment length" $OUTDIR/logs/$CHIP_ID"Rep"$i.pr1.predictd | perl -ane 'print $F[scalar(@F)-2];'`
    j=$(( j + 1 ))
    macs2 predictd -i $OUTDIR/tagAlign/$CHIP_ID"Rep"$i.pr2.tagAlign.gz -g $GENOME --outdir $OUTDIR/logs 2> $OUTDIR/logs/$CHIP_ID"Rep"$i.pr2.predictd
    FRAGMENT_LENGTH[$j]=`grep "predicted fragment length" $OUTDIR/logs/$CHIP_ID"Rep"$i.pr2.predictd | perl -ane 'print $F[scalar(@F)-2];'`
done
echo "done"

## peak calling on all the randomly split ChIP samples using macs
echo -n "Call peaks on each of the randomly split ChIP samples.. "
j=0
for (( i=1; i<=$CHIP_COUNT; i++ )); do
    j=$(( j + 1 ))
    SHIFT_SIZE=`echo $(( FRAGMENT_LENGTH[$j] ))`
    macs2 callpeak -t $OUTDIR/tagAlign/$CHIP_ID"Rep"$i.pr1.tagAlign.gz -c $OUTDIR/tagAlign/$CONTROL_ID"Rep0".tagAlign.gz -f BED -n $OUTDIR/macs/$CHIP_ID"Rep"$i".pr1_Vs_"$CONTROL_ID"Rep0" -g $GENOME -p $MACS_P_VALUE --nomodel --extsize $SHIFT_SIZE 2>$OUTDIR/logs/$CHIP_ID"Rep"$i".pr1_Vs_"$CONTROL_ID"Rep0.log"
    j=$(( j + 1 ))
    SHIFT_SIZE=`echo $(( FRAGMENT_LENGTH[$j] ))`
    macs2 callpeak -t $OUTDIR/tagAlign/$CHIP_ID"Rep"$i.pr2.tagAlign.gz -c $OUTDIR/tagAlign/$CONTROL_ID"Rep0".tagAlign.gz -f BED -n $OUTDIR/macs/$CHIP_ID"Rep"$i".pr2_Vs_"$CONTROL_ID"Rep0" -g $GENOME -p $MACS_P_VALUE --nomodel --extsize $SHIFT_SIZE 2>$OUTDIR/logs/$CHIP_ID"Rep"$i".pr2_Vs_"$CONTROL_ID"Rep0.log" 
done
echo "done"

## retrieve top 100,000 peaks sorted by p-value
echo -n "Retrieve top 100,000 peaks sorted by p-value.. "
for (( i=1; i<=$CHIP_COUNT; i++ )); do
    sort -k 8nr,8nr $OUTDIR/macs/$CHIP_ID"Rep"$i".pr1_Vs_"$CONTROL_ID"Rep0_peaks.narrowPeak" | head -n 100000 | gzip -c > $OUTDIR/macs/$CHIP_ID"Rep"$i".pr1_Vs_"$CONTROL_ID"Rep0.regionPeak.gz"
    sort -k 8nr,8nr $OUTDIR/macs/$CHIP_ID"Rep"$i".pr2_Vs_"$CONTROL_ID"Rep0_peaks.narrowPeak" | head -n 100000 | gzip -c > $OUTDIR/macs/$CHIP_ID"Rep"$i".pr2_Vs_"$CONTROL_ID"Rep0.regionPeak.gz"
done
echo "done"

##########################################################################
############ CALL PEAKS ON PSEUDOREPLICATES OF POOLED CHIP SAMPLES
##########################################################################

## randomly split pooled input ChIP samples into two parts
echo -n "Randomly split pooled ChIP samples into two parts.. "
filename=$OUTDIR/tagAlign/$CHIP_ID"Rep0.tagAlign.gz"
nlines=$( gunzip -c $filename | wc -l );
nlines=$(( (nlines + 1) / 2 ))
gunzip -c ${filename} | shuf | split -l ${nlines} - $OUTDIR/tagAlign/$CHIP_ID"Rep0"
gzip $OUTDIR/tagAlign/$CHIP_ID"Rep0aa"
gzip $OUTDIR/tagAlign/$CHIP_ID"Rep0ab"
mv $OUTDIR/tagAlign/$CHIP_ID"Rep0aa.gz" $OUTDIR/tagAlign/$CHIP_ID"Rep0.pr1.tagAlign.gz"
mv $OUTDIR/tagAlign/$CHIP_ID"Rep0ab.gz" $OUTDIR/tagAlign/$CHIP_ID"Rep0.pr2.tagAlign.gz"
echo "done"

## determine fragment length for all randomly split pooled ChIP samples using macs
echo -n "Determine fragment length for all randomly splot pooled ChIP samples.. "
FRAGMENT_LENGTH=(0)
j=0
j=$(( j + 1 ))
macs2 predictd -i $OUTDIR/tagAlign/$CHIP_ID"Rep0.pr1.tagAlign.gz" -g $GENOME --outdir $OUTDIR/logs 2> $OUTDIR/logs/$CHIP_ID"Rep0.pr1.predictd"
FRAGMENT_LENGTH[$j]=`grep "predicted fragment length" $OUTDIR/logs/$CHIP_ID"Rep0.pr1.predictd" | perl -ane 'print $F[scalar(@F)-2];'`
j=$(( j + 1 ))
macs2 predictd -i $OUTDIR/tagAlign/$CHIP_ID"Rep0.pr2.tagAlign.gz" -g $GENOME --outdir $OUTDIR/logs 2> $OUTDIR/logs/$CHIP_ID"Rep0.pr2.predictd"
FRAGMENT_LENGTH[$j]=`grep "predicted fragment length" $OUTDIR/logs/$CHIP_ID"Rep0.pr2.predictd" | perl -ane 'print $F[scalar(@F)-2];'`
echo "done"

## peak calling on the randomly split pooled ChIP samples using macs
echo -n "Call peaks on randomly split pooled ChIP samples.. "
j=0
j=$(( j + 1 ))
SHIFT_SIZE=`echo $(( FRAGMENT_LENGTH[$j] ))`
macs2 callpeak -t $OUTDIR/tagAlign/$CHIP_ID"Rep0.pr1.tagAlign.gz" -c $OUTDIR/tagAlign/$CONTROL_ID"Rep0".tagAlign.gz -f BED -n $OUTDIR/macs/$CHIP_ID"Rep0.pr1_Vs_"$CONTROL_ID"Rep0" -g $GENOME -p $MACS_P_VALUE --nomodel --extsize $SHIFT_SIZE 2>$OUTDIR/logs/$CHIP_ID"Rep0.pr1_Vs_"$CONTROL_ID"Rep0.log" 
j=$(( j + 1 ))
SHIFT_SIZE=`echo $(( FRAGMENT_LENGTH[$j] ))`
macs2 callpeak -t $OUTDIR/tagAlign/$CHIP_ID"Rep0.pr2.tagAlign.gz" -c $OUTDIR/tagAlign/$CONTROL_ID"Rep0".tagAlign.gz -f BED -n $OUTDIR/macs/$CHIP_ID"Rep0.pr2_Vs_"$CONTROL_ID"Rep0" -g $GENOME -p $MACS_P_VALUE --nomodel --extsize $SHIFT_SIZE 2>$OUTDIR/logs/$CHIP_ID"Rep0.pr2_Vs_"$CONTROL_ID"Rep0.log"
echo "done"

## retrieve top 100,000 peaks sorted by p-value
echo -n "Retrieve top 100,000 peaks sorted by p-value.. "
sort -k 8nr,8nr $OUTDIR/macs/$CHIP_ID"Rep0.pr1_Vs_"$CONTROL_ID"Rep0_peaks.narrowPeak" | head -n 100000 | gzip -c > $OUTDIR/macs/$CHIP_ID"Rep0.pr1_Vs_"$CONTROL_ID"Rep0.regionPeak.gz"
sort -k 8nr,8nr $OUTDIR/macs/$CHIP_ID"Rep0.pr2_Vs_"$CONTROL_ID"Rep0_peaks.narrowPeak" | head -n 100000 | gzip -c > $OUTDIR/macs/$CHIP_ID"Rep0.pr2_Vs_"$CONTROL_ID"Rep0.regionPeak.gz"
echo "done"


##########################################################################
############ IDR ANALYSIS ON ORIGINAL CHIP SAMPLES
##########################################################################

## move to idrCode directory
cd $PROGDIR

## IDR analysis on the original ChIP samples
echo -n "IDR analysis on the original ChIP samples.. "
COMMAND=""
for (( i=1; i<=$CHIP_COUNT; i++ )); do
    for (( j=i+1; j<=$CHIP_COUNT; j++ )); do
        if [ -f $OUTDIR/macs/$CHIP_ID"Rep"$i"_Vs_"$CONTROL_ID"Rep0.regionPeak.gz" ]; then
            gunzip $OUTDIR/macs/$CHIP_ID"Rep"$i"_Vs_"$CONTROL_ID"Rep0.regionPeak.gz" -f
            gunzip $OUTDIR/macs/$CHIP_ID"Rep"$j"_Vs_"$CONTROL_ID"Rep0.regionPeak.gz" -f
        fi
        Rscript batch-consistency-analysis.r $OUTDIR/macs/$CHIP_ID"Rep"$i"_Vs_"$CONTROL_ID"Rep0.regionPeak" $OUTDIR/macs/$CHIP_ID"Rep"$j"_Vs_"$CONTROL_ID"Rep0.regionPeak" -1 $OUTDIR/consistency/reps/$CHIP_ID"Rep"$i"_Vs_"$CHIP_ID"Rep"$j 0 F p.value
        COMMAND="$COMMAND "$OUTDIR/consistency/reps/$CHIP_ID"Rep"$i"_Vs_"$CHIP_ID"Rep"$j
    done
done

## plot IDR plots
numerator=$( factorial $CHIP_COUNT );
denominator=$(( 2 * $( factorial `expr $CHIP_COUNT - 2` ) ));
plot_argument=`expr $numerator / $denominator`;
#echo -e "$numerator\t$denominator\t$plot_argument";
Rscript batch-consistency-plot.r $plot_argument $OUTDIR/consistency/reps/chipSampleAllReps $COMMAND 2>/dev/null
echo "done"

## IDR analysis on the pseudo-replicates of the individual ChIP samples"
echo -n "IDR analysis on pseudo-replicates of the individual ChIP samples.. "
COMMAND=""
for (( i=1; i<=$CHIP_COUNT; i++ )); do
    if [ -f $OUTDIR/macs/$CHIP_ID"Rep"$i".pr1_Vs_"$CONTROL_ID"Rep0.regionPeak.gz" ]; then
        gunzip $OUTDIR/macs/$CHIP_ID"Rep"$i".pr1_Vs_"$CONTROL_ID"Rep0.regionPeak.gz" -f
        gunzip $OUTDIR/macs/$CHIP_ID"Rep"$i".pr2_Vs_"$CONTROL_ID"Rep0.regionPeak.gz" -f
    fi
    Rscript batch-consistency-analysis.r $OUTDIR/macs/$CHIP_ID"Rep"$i".pr1_Vs_"$CONTROL_ID"Rep0.regionPeak" $OUTDIR/macs/$CHIP_ID"Rep"$i".pr2_Vs_"$CONTROL_ID"Rep0.regionPeak" -1 $OUTDIR/consistency/selfPseudoReps/$CHIP_ID"Rep"$i".pr1_Vs_"$CHIP_ID"Rep"$i".pr2" 0 F p.value
    COMMAND="$COMMAND "$OUTDIR/consistency/selfPseudoReps/$CHIP_ID"Rep"$i".pr1_Vs_"$CHIP_ID"Rep"$i".pr2"
done

## plot IDR plots
Rscript batch-consistency-plot.r $plot_argument $OUTDIR/consistency/selfPseudoReps/chipSampleAllSelfReps $COMMAND 2>/dev/null
echo "done"

## IDR analysis on the pseudo-replicates of the pooled ChIP samples"
echo -n "IDR analysis on pseudo-replicates of the pooled ChIP samples.. "
if [ -f $OUTDIR/macs/$CHIP_ID"Rep0.pr1_Vs_"$CONTROL_ID"Rep0.regionPeak.gz" ]; then
    gunzip $OUTDIR/macs/$CHIP_ID"Rep0.pr1_Vs_"$CONTROL_ID"Rep0.regionPeak.gz" -f
    gunzip $OUTDIR/macs/$CHIP_ID"Rep0.pr2_Vs_"$CONTROL_ID"Rep0.regionPeak.gz" -f
fi

Rscript batch-consistency-analysis.r $OUTDIR/macs/$CHIP_ID"Rep0.pr1_Vs_"$CONTROL_ID"Rep0.regionPeak" $OUTDIR/macs/$CHIP_ID"Rep0.pr2_Vs_"$CONTROL_ID"Rep0.regionPeak" -1 $OUTDIR/consistency/pooledPseudoReps/$CHIP_ID"Rep0.pr1_Vs_"$CHIP_ID"Rep0.pr2" 0 F p.value

## plot IDR plots
Rscript batch-consistency-plot.r 1 $OUTDIR/consistency/pooledPseudoReps/chipSamplePooledReps $OUTDIR/consistency/pooledPseudoReps/$CHIP_ID"Rep0.pr1_Vs_"$CHIP_ID"Rep0.pr2" 2>/dev/null
echo "done"


##########################################################################
############ GETTING THRESHOLDS TO TRUNCATE PEAK LISTS
##########################################################################

## compute threshold on ChIP samples
echo -n "Compute threshold on ChIP samples.. "
peakReps=()
for (( i=1; i<=$CHIP_COUNT; i++ )); do
    for (( j=i+1; j<=$CHIP_COUNT; j++ )); do
        #peakReps+=( $(awk '$11 <= 0.01 {print $0}' $OUTDIR/consistency/reps/$CHIP_ID"Rep"$i"_Vs_"$CHIP_ID"Rep"$j"-overlapped-peaks.txt" | wc -l) )
        peakReps+=( $(awk '$11 <= '$IDR_THRESHOLD' {print $0}' $OUTDIR/consistency/reps/$CHIP_ID"Rep"$i"_Vs_"$CHIP_ID"Rep"$j"-overlapped-peaks.txt" | wc -l) )
    done
done
echo "done"

## compute threshold on pseudo-replicates of ChIP samples
echo -n "Compute threshold on pseudo-replicates of ChIP samples.. "
peakSelfPseudoReps=()
for (( i=1; i<=$CHIP_COUNT; i++ )); do
    peakSelfPseudoReps+=( $(awk '$11 <= 0.02 {print $0}' $OUTDIR/consistency/selfPseudoReps/$CHIP_ID"Rep"$i".pr1_Vs_"$CHIP_ID"Rep"$i".pr2-overlapped-peaks.txt" | wc -l) )
done
echo "done"

## compute threshold on pooled pseudo-replicates of ChIP samples
echo -n "Compute threshold on pooled pseudo-replicates of ChIP samples.. "
peakPooledPseudoReps=( $(awk '$11 <= 0.0025 {print $0}' $OUTDIR/consistency/pooledPseudoReps/$CHIP_ID"Rep0.pr1_Vs_"$CHIP_ID"Rep0.pr2-overlapped-peaks.txt" | wc -l) )
echo "done"

## compute conservative and optimal number of peaks
max_numPeaks_Rep=${peakReps[0]}
for i in "${peakReps[@]}"; do
    if [ "$i" -gt "$max_numPeaks_Rep" ]; then
        max_numPeaks_Rep=$i
    fi
done

echo "max_numPeaks_Rep: $max_numPeaks_Rep";

echo -n "self-consistency thresholds: "
for i in "${peakSelfPseudoReps[@]}"; do
    echo -n "$i ";
done
echo

echo "numPeaks_Rep0: $peakPooledPseudoReps"

## determine final set of peaks (conservative)
zcat $OUTDIR/macs/$CHIP_ID"Rep0_Vs_"$CONTROL_ID"Rep0.regionPeak.gz" | sort -k8nr,8nr | head -n $max_numPeaks_Rep | gzip -c > $OUTDIR/macs/"conservative."$CHIP_ID"Rep0_Vs_"$CONTROL_ID"Rep0.regionPeak.gz"
## determine final set of summits (conservative)
cat $OUTDIR/macs/$CHIP_ID"Rep0_Vs_"$CONTROL_ID"Rep0_summits.bed" | sort -k5nr,5nr | head -n $max_numPeaks_Rep | gzip -c > $OUTDIR/macs/"conservative."$CHIP_ID"Rep0_Vs_"$CONTROL_ID"Rep0.summits.bed.gz"

## determine final set of peaks (optimal)
if [ "$max_numPeaks_Rep" -gt "$peakPooledPseudoReps" ]; then
    zcat $OUTDIR/macs/$CHIP_ID"Rep0_Vs_"$CONTROL_ID"Rep0.regionPeak.gz" | sort -k8nr,8nr | head -n $max_numPeaks_Rep | gzip -c > $OUTDIR/macs/"optimal."$CHIP_ID"Rep0_Vs_"$CONTROL_ID"Rep0.regionPeak.gz"
    cat $OUTDIR/macs/$CHIP_ID"Rep0_Vs_"$CONTROL_ID"Rep0_summits.bed" | sort -k5nr,5nr | head -n $max_numPeaks_Rep | gzip -c > $OUTDIR/macs/"optimal."$CHIP_ID"Rep0_Vs_"$CONTROL_ID"Rep0.summits.bed.gz"
else
    zcat $OUTDIR/macs/$CHIP_ID"Rep0_Vs_"$CONTROL_ID"Rep0.regionPeak.gz" | sort -k8nr,8nr | head -n $peakPooledPseudoReps | gzip -c > $OUTDIR/macs/"optimal."$CHIP_ID"Rep0_Vs_"$CONTROL_ID"Rep0.regionPeak.gz"
    cat $OUTDIR/macs/$CHIP_ID"Rep0_Vs_"$CONTROL_ID"Rep0_summits.bed" | sort -k5nr,5nr | head -n $peakPooledPseudoReps | gzip -c > $OUTDIR/macs/"optimal."$CHIP_ID"Rep0_Vs_"$CONTROL_ID"Rep0.summits.bed.gz"
fi

## compute enrichment and quality measure for input ChIP-seq data
echo -n "Compute enrichment and quality measure for input ChIP-seq data... "
if [ -e "$PROGDIR/phantompeakqualtools/run_spp.R" ]; then
    for (( i=1; i<=$CHIP_COUNT; i++ )); do
        ## details about the output file are here: https://github.com/kundajelab/phantompeakqualtools#typical-usage
        Rscript $PROGDIR/phantompeakqualtools/run_spp.R -c=$OUTDIR/tagAlign/$CHIP_ID"Rep"$i.tagAlign.gz -fdr=$IDR_THRESHOLD -savp -odir=$OUTDIR/quality -out=$OUTDIR/quality/quality.txt &> $OUTDIR/logs/quality.log
    done
fi
echo "done"

#for (( i=1; i<=$CHIP_COUNT; i++ )); do
#    makeTagDirectory $OUTDIR/homer/$CHIP_ID"Rep"$i $INDIR/$CHIP_ID"Rep"$i.bam
#done

exit
