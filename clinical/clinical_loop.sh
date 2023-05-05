#!/bin/bash
while getopts i:c:t:b:j:o: flag
do
    case "${flag}" in
        i) input=${OPTARG};; #input bamlist
        c) cov=${OPTARG};; #minimum coverage
        t) trim=${OPTARG};; #end trim
        b) refbg=${OPTARG};; #reference name (hg19/hg38)
        j) thred=${OPTARG};; #thread
        o) odir=${OPTARG};; #output directory
    esac
done
rm -rf tmp
mkdir tmp
basename=$( cat $input | cut -d, -f1 )
ochk="T"
if [ -z $input ]; then
	echo "Input bamlist was not given, exiting"
	exit 1
fi

if [ -z $cov ]; then
	echo "No minimum coverage was set, default is 1"
	cov=1
elif [ $cov -lt 0 ]; then
	cov=1
fi

if [ -z $trim ]; then
	echo "No endtrim was set, default is 2"
	trim=2
elif [ $bqthr -lt 0 ]; then
	trim=0
fi

if [ -z $refbg ]; then
	echo "No reference was given, assuming hg19"
	refbg="hg19"
fi

if [ -z $thred ]; then
	echo "No thread was given, default is 12"
	thred=12
elif [ $bqthr -lt 0 ]; then
	thred=1
fi

if [ -z $odir ]; then
	echo "No output directory was given, outputting to source directory"
	ochk="F"
fi
samtools idxstats $( cat $input | head -1 ) | cut -f 1 | uniq | head -25 > ./tmp/refnames.txt
for j in $basename; do
	file=$( echo ${j##*/} )
	bbn=$( echo "$file" | cut -d'.' -f1 )
	if [ $ochk == "F" ]; then odir=$( dirname $j ); fi
	echo "$odir"
	Rscript --vanilla clinical1.r full_panel.csv F ./tmp/ ./tmp/refnames.txt $refbg 
	angsd -gl 2 -nThreads $thred -doGlf 2 -minQ 30 -minMapQ 25 -i $j -rf ./tmp/snps.bed -doMajorMinor 1 -doCounts 1 -dumpCounts 3 -trim $trim -out "$odir"/"$bbn"_clin
	Rscript --vanilla clinical2.r "$odir"/"$bbn"_clin.counts.gz "$odir"/"$bbn"_clin.pos.gz $cov $cov full_panel.csv $refbg ./tmp/refnames.txt "$odir"/"$bbn"_clin.beagle.gz "$odir/" "$bbn"
done
