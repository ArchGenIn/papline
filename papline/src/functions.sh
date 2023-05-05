#!/bin/bash

source ./tmp/paramwork

function defaults {
	#check if parameters aren set, if not, set defaults
	#number of threads
	if [ -z $thred ]; then
		thred=24
	elif [ $thred -lt 1 ]; then
		thred=24
	fi
	#BQ
	if [ -z $bqthr ]; then
		bqthr=30
	elif [ $bqthr -lt 0 ]; then
		bqthr=0
	fi
	#seed
	fi [ -z $seeds ]; then
		seeds=2
	elif [ $seeds -lt 1 ]; then
		seeds=0
	fi
	#overlapping base number
	if [ -z $overl ]; then
		overl=5
	elif [ $overl -lt 1 ]; then
		overl=0
	fi
	#minimum read length
	if [ -z $minle ]; then
		minle=15
	elif [ $minle -lt 1 ]
		minle=0
	fi
	#bases trimmed off from end of the reads
	if [ -z $ctrim ]; then
		ctrim=2
	elif [ $ctrim -lt 1 ]; then
		ctrim=0
	fi
	#MQ
	if [ -z $mapqt ]; then
		mapqt=30
	elif [ $mapqt -lt 0 ]; then
		mapqt=0
	fi
	#mtDNA minimum coverage
	if [ -t $mtcov ]; then
		mtcov=2
	elif [ $mtcov -lt 1 ]; then
		mtcov=0
	fi
	#disease mincov
	if [ -z $disis ]; then
		disis=3
	elif [ $disis -lt 1 ]; then
		disis=0
	fi
	#pigmentation mincov
	if [ -z $pigme ]; then
		pigme=1
	elif [ $pigme -lt 1 ]; then
		pigme=1
	fi
	#Variant discovery background reference
	if [ -z $refbg ]; then
		refbg="hg19"
	fi
	#whether run fastqc T/F
	if [ -z $fstqc ]; then
		fstqc=TRUE
	fi
	#illumina adapter trim T/F
	if [ -z $atrim ]; then
		atrim=FALSE
	fi
	#inner barcode trim T/F
	if [ -z $btrim ]; then
		btrim=TRUE
	fi
	#discard barcodeless reads T/F
	if [ -z $bdisc ]; then
		bdisc=TRUE
	fi
	#set default aligner type
	if [ -z $align ]; then
		align=Bowtie2_EVF
	fi
	#keep raw bam T/F
	if [ -z $krbam ]; then
		krbam=FALSE
	fi
	#remove duplicates T/F
	if [ -z $rmdup ]; then
		rmdup=TRUE
	fi
	#run mapdamage T/F
	if [ -z $mapda ]; then
		mapda=TRUE
	fi
	#run contammix T/F
	if [ -z $cntmx ]; then
		cntmx="Inclusive"
	fi
	#mtDNA fasta creating method
	if [ -z $mtfas ]; then
		mtfas="Intelligent"
	fi
	#call only transversions T/F
	if [ -z $cotrv ]; then
		cotrv=FALSE
	fi
	#call pseudohaploide genome T/F
	if [ -z $psudo ]; then
		psudo="None"
	fi
	#create eigenstrat format output T/F
	if [ -z $eigen ]; then
		eigen=FALSE
	fi
	#create plink format output T/F
	if [ -z $plink ]; then
		plink=FALSE
	fi
}

#function for creating directories
function dirmake {
	#if papline structure exists T/F
	if [ $workf == "PAPLINE" ] || [ $workf == "papline" ] || [ $workf == "PAPline" ] || [ $workf == "Papline" ]; then
		ppp=TRUE
		erendir=FALSE
		erindir=FALSE
	else
		ppp=FALSE
		erendir=FALSE
		erindir=FALSE
	fi
	#if papline structure does not exist create endir
	if [ $ppp == FALSE ]; then
		#if endir is given, create endir within
		if [ ! -z $endir ]; then
			rm -rf "$endir"/"$runam"
			mkdir "$endir"/"$runam"
			if [ ! -d $endir ]; then
				echo "Given end directory can not be created"
				exit 1
			fi
			basename=$( cat $listf | cut -d, -f1 )
			for i in basename; do
				echo "Create directories"
				mkdir -p "$endir"/"$runam"/"$i"
				mkdir "$endir"/"$runam"/"$i"/RUNLOGS
			done 2> "$endir"/error_dirmake.log
			erendir=TRUE
		#if endir is not given, create endir within indir
		else
			rm -rf "$indir"/"$runam"
			mkdir "$indir"/"$runam"
			basename=$( cat $listf | cut -d, -f1 )
			for i in $basename; do
				echo "Create directories"
				mkdir -p "$indir"/"$runam"/"$i"
				mkdir "$indir"/"$runam"/"$i"/RUNLOGS
			done 2> "$indir"/error_dirmake.log
			erindir=TRUE
		fi
	#if papline structure does exist
	elif [ $ppp == TRUE ]; then
		#but papli_log file is missing, give error
		if [ ! -s "$indir"/papli_log ]
			echo "papli_log file does not exist or is renamed or is empty despite workf=PAPLINE, check your input or PARAMETERFILE"
			exit 1
		fi
	fi
}

#function for running fastqc analysis FIN
function fastqc {
	if ! command -v fastqc &> /dev/null; then
		echo "FASTQC is not installed despite: fastqc=TRUE, try sudo apt install fastqc -y, also check dependencies"
		exit 1
	fi
	basename=$( cat $listf | cut -d, -f1 )
	echo "Quality check of FASTQ files"
	#if endir was not given, set indir as endir
	if [ $erindir == TRUE ]; then
		endir=$( echo $indir )
	fi
	for i in $basename; do
		shopt -s nullglob
		echo "FastQC run for $i"
		fastqc "$indir"/"$i"*.{fastq,fastq.gz,fq,fq.gz} 
		shopt -u nullglob
		rm -rf "$endir"/"$runam"/"$i"/FASTQC
		mkdir "$endir"/"$runam"/"$i"/FASTQC
		mv "$indir"/"$i"*.html "$endir"/"$runam"/"$i"/FASTQC/
		mv "$indir"/"$i"*.zip "$endir"/"$runam"/"$i"/FASTQC/
	done 2> x
	mv x "$endir"/"$runam"/error_fastqc.log
}

#function for trimming low quality bases FIN
function basequalitytrim {
	basename=$( cat $listf | cut -d, -f1 )
	echo "BQ trim of FASTQ files"
	if ! command -v cutadapt &> /dev/null; then
		echo "Cutadapt is not installed despite bqthr is set, try sudo apt install cutadapt -y, also check dependencies"
		exit 1
	fi
	if [ $erindir == TRUE ]; then
		endir=$( echo $indir )
	fi
	for i in $basename; do	
		rm -rf "$endir"/"$runam"/"$i"/BQTRIM
		mkdir "$endir"/"$runam"/"$i"/BQTRIM
		rm -rf "$endir"/"$runam"/"$i"/BQTRIM1
		mkdir "$endir"/"$runam"/"$i"/BQTRIM1
		if [ ! -d "$endir"/"$runam"/"$i"/RUNLOGS ]; then
			mkdir  "$endir"/"$runam"/"$i"/RUNLOGS
		fi
		shopt -s nullglob
		for j in "$indir"/"$i"*.{fastq,fastq.gz,fq,fq.gz}; do
			echo "$j"
			file=$( echo ${j##*/} )
			bbn=$( echo "$file" | cut -d'.' -f1 )
			cutadapt -j $thred -q $bqthr -o "$endir"/"$runam"/"$i"/BQTRIM1/"$bbn"_"$bqthr"bq.fastq.gz $j > "$endir"/"$runam"/"$i"/RUNLOGS/"$bbn".qualtrim.full.txt
			echo $( awk '{for(i=1;i<=NF;i++) if ($i=="basepairs") print $(i+2)}' "$endir"/"$runam"/"$i"/RUNLOGS/"$bbn".qualtrim.full.txt ) | tr --delete "," > "$endir"/"$runam"/"$i"/RUNLOGS/"$bbn".totalbp.txt
			echo $( awk '{for(i=1;i<=NF;i++) if ($i=="Quality-trimmed:") print $(i+1)}' "$endir"/"$runam"/"$i"/RUNLOGS/"$bbn".qualtrim.full.txt ) | tr --delete "," > "$endir"/"$runam"/"$i"/RUNLOGS/"$bbn".qualtrim.txt
			expr $( zcat "$endir"/"$runam"/"$i"/BQTRIM1/"$bbn".bq"$bqthr".fastq.gz | wc -l ) / 4 &> "$endir"/"$runam"/"$i"/RUNLOGS/"$bbn".readno1.txt
		done
		shopt -u nullglob
		array1=($(cat "$endir"/"$runam"/"$i"/RUNLOGS/*R1*readno1.txt))
		sum1=$(IFS=+; echo "$((${array1[*]}))")
		echo $sum1 > "$endir"/"$runam"/"$i"/RUNLOGS/"$i".Starting_read_number.txt
		array2=($(cat "$endir"/"$runam"/"$i"/RUNLOGS/*totalbp.txt))
		sum2=$(IFS=+; echo "$((${array2[*]}))")
		echo $sum2 > "$endir"/"$runam"/"$i"/RUNLOGS/"$i".Starting_basepairs.txt
		array3=($(cat "$endir"/"$runam"/"$i"/RUNLOGS/*qualtrim.txt))
		sum3=$(IFS=+; echo "$((${array3[*]}))")
		echo $sum3 > "$endir"/"$runam"/"$i"/RUNLOGS/"$i".Trimmed_off_basepairs.txt
		rm -rf "$endir"/"$runam"/"$i"/RUNLOGS/*.totalbp.txt
		rm -rf "$endir"/"$runam"/"$i"/RUNLOGS/*.qualtrim.txt
		rm -rf "$endir"/"$runam"/"$i"/RUNLOGS/*.readno1.txt
		rm -rf "$endir"/"$runam"/"$i"/RUNLOGS/*.qualtrim.full.txt
		shopt -s nullglob
		countr1=$( ls -1 "$endir"/"$runam"/"$i"/BQTRIM1/"$i"*{r1,R1}.{fastq,fq,fastq.gz,fq.gz} 2>/dev/null | wc -l )
		countr2=$( ls -1 "$endir"/"$runam"/"$i"/BQTRIM1/"$i"*{r2,R2}.{fastq,fq,fastq.gz,fq.gz} 2>/dev/null | wc -l )
		if [ $countr1 -eq 0 ] && [ $countr2 -eq 0 ]; then #if SE
			zcat "$endir"/"$runam"/"$i"/BQTRIM1/"$i"*fastq.gz| gzip -f > "$endir"/"$runam"/"$i"/BQTRIM/"$i"_"$bqthr"bp.fastq.gz
			rm -rf "$endir"/"$runam"/"$i"/BQTRIM1
		elif [ $countr1 -eq $countr2 ]; then #if PE
			zcat rm -rf "$endir"/"$runam"/"$i"/BQTRIM1/"$i"*{r1,R1}*fastq.gz | gzip -f > rm -rf "$endir"/"$runam"/"$i"/BQTRIM/"$i"_R1_"$bqthr"bp.fastq.gz
			zcat rm -rf "$endir"/"$runam"/"$i"/BQTRIM1/"$i"*{r2,R2}*fastq.gz | gzip -f > rm -rf "$endir"/"$runam"/"$i"/BQTRIM/"$i"_R2_"$bqthr"bp.fastq.gz
			rm -rf rm -rf "$endir"/"$runam"/"$i"/BQTRIM1
		else
			echo "Uneven number of PE reads, check input data. If SE, delete R1 and R2 tags in filenames."
			exit 1
		fi
		shopt -u nullglob
	done 2> "$endir"/"$runam"/error_fastqc.log
}

#illumina adapter trimming function FIN
function adaptertrim {
	basename=$( cat $listf | cut -d, -f1 )
	echo "ADAPTER trim of FASTQ files"
	if [ $erindir == TRUE ]; then
		endir=$( echo $indir )
	fi
	if ! command -v cutadapt &> /dev/null; then
		echo "Cutadapt is not installed despite bqthr is set, try sudo apt install cutadapt -y, also check dependencies"
		exit 1
	fi
	rm -rf "$endir"/"$runam"/"$i"/ADAPTERTRIM
	mkdir "$endir"/"$runam"/"$i"/ADAPTERTRIM
	rm -rf "$endir"/"$runam"/"$i"/ADAPTERTRIM1
	mkdir "$endir"/"$runam"/"$i"/ADAPTERTRIM1
	#setting difference ratio
	r1e=$( echo "0.$( echo $((1000*$seeds/${#p5ada})) )" )
	r2e=$( echo "0.$( echo $((1000*$seeds/${#p7ada})) )" )
	for i in $basename; do
		if [ -d "$endir"/"$runam"/"$i"/BQTRIM ] && [ `ls -1 "$endir"/"$runam"/"$i"/BQTRIM/*{fastq.gz,fq.gz,fastq,fq} 2>/dev/null | wc -l` != 0 ]; then #papli structure
			cr1=$( ls -l "$endir"/"$runam"/"$i"/BQTRIM/*{R1,r1}*{fastq.gz,fq.gz,fastq,fq} 2>dev/null | wc -l )
			cr2=$( ls -l "$endir"/"$runam"/"$i"/BQTRIM/*{R2,r2}*{fastq.gz,fq.gz,fastq,fq} 2>dev/null | wc -l )
			if [ $cr1 == $cr2 ] && [ $cr1 != 0 ]; then #if PE
				shopt -s nullglob
				for j in "$endir"/"$runam"/"$i"/BQTRIM/*{R1,r1}*{fastq.gz,fq.gz,fastq,fq}; do
					echo "$j"
					file=$( echo ${j##*/} )
					bbn=$( echo "$file" | cut -d'.' -f1 )
					cutadapt -j $thred -a $p5ada -e $r1e -o "$endir"/"$runam"/"$i"/ADAPTERTRIM1/"$bbn".atrim.fastq.gz $j
				done
				for j in "$endir"/"$runam"/"$i"/BQTRIM/*{R2,r2}*{fastq.gz,fq.gz,fastq,fq}; do
					echo "$j"
					file=$( echo ${j##*/} )
					bbn=$( echo "$file" | cut -d'.' -f1 )
					cutadapt -j $thred -a $p7ada -e $r2e -o "$endir"/"$runam"/"$i"/ADAPTERTRIM1/"$bbn".atrim.fastq.gz $j
				done
				shopt -u nullglob
			elif [ $cr1 == 0 ] && [ $cr2 == 0 ]; then #if SE
				shopt -s nullglob
				for j in "$endir"/"$runam"/"$i"/BQTRIM/*{fastq.gz,fq.gz,fastq,fq}; do
					echo "$j"
					file=$( echo ${j##*/} )
					bbn=$( echo "$file" | cut -d'.' -f1 )
					cutadapt -j $thred -g $p5ada -a $p7ada -e $r2e -o "$endir"/"$runam"/"$i"/ADAPTERTRIM1/"$bbn".atrim.fastq.gz $j
				shopt -u nullglob
				done
			else #if uneven
				echo "Uneven number of PE reads, check input data. If SE, delete R1 and R2 tags in filenames."
				exit 1
			fi
		elif [ `ls -1 "$indir"/*{fastq.gz,fq.gz,fastq,fq} 2>/dev/null | wc -l` != 0 ]; then #fq outside
			cr1=$( ls -l "$indir"/"$i"*{R1,r1}*{fastq.gz,fq.gz,fastq,fq} 2>dev/null | wc -l )
			cr2=$( ls -l "$indir"/"$i"*{R2,r2}*{fastq.gz,fq.gz,fastq,fq} 2>dev/null | wc -l )
			if [ $cr1 == $cr2 ] && [ $cr1 != 0 ]; then #if PE
				shopt -s nullglob
				for j in "$indir"/"$i"*{R1,r1}*{fastq.gz,fq.gz,fastq,fq}; do
					echo "$j"
					file=$( echo ${j##*/} )
					bbn=$( echo "$file" | cut -d'.' -f1 )
					cutadapt -j $thred -a $p5ada -e $r1e -o "$endir"/"$runam"/"$i"/ADAPTERTRIM1/"$bbn".atrim.fastq.gz $j
				done
				for j in "$indir"/"$i"*{R2,r2}*{fastq.gz,fq.gz,fastq,fq}; do
					echo "$j"
					file=$( echo ${j##*/} )
					bbn=$( echo "$file" | cut -d'.' -f1 )
					cutadapt -j $thred -a $p7ada -e $r2e -o "$endir"/"$runam"/"$i"/ADAPTERTRIM1/"$bbn".atrim.fastq.gz $j
				done
				shopt -u nullglob
			elif [ $cr1 == 0 ] && [ $cr2 == 0 ]; then #if SE
				shopt -s nullglob
				for j in "$indir"/"$i"*{fastq.gz,fq.gz,fastq,fq}; do
					echo "$j"
					file=$( echo ${j##*/} )
					bbn=$( echo "$file" | cut -d'.' -f1 )
					cutadapt -j $thred -g $p5ada -a $p7ada -e $r2e -o "$endir"/"$runam"/"$i"/ADAPTERTRIM1/"$bbn".atrim.fastq.gz $j
				done
				shopt -u nullglob
			else #if uneven
				echo "Uneven number of PE reads, check input data. If SE, delete R1 and R2 tags in filenames."
				exit 1
			fi
		else
			echo "FASTQ files do not exist in given directories, exiting 1"
			exit 1
		fi
		shopt -s nullglob
		countr1=$( ls -1 "$endir"/"$runam"/"$i"/ADAPTERTRIM1/"$i"*{r1,R1}.{fastq,fq,fastq.gz,fq.gz} 2>/dev/null | wc -l )
		countr2=$( ls -1 "$endir"/"$runam"/"$i"/ADAPTERTRIM1/"$i"*{r2,R2}.{fastq,fq,fastq.gz,fq.gz} 2>/dev/null | wc -l )
		if [ $countr1 -eq 0 ] && [ $countr2 -eq 0 ]; then #if SE
			zcat "$endir"/"$runam"/"$i"/ADAPTERTRIM1/"$i"*fastq.gz| gzip -f > "$endir"/"$runam"/"$i"/ADAPTERTRIM/"$i".atrim.fastq.gz
			rm -rf "$endir"/"$runam"/"$i"/ADAPTERTRIM1
		elif [ $countr1 -eq $countr2 ]; then #if PE
			zcat rm -rf "$endir"/"$runam"/"$i"/ADAPTERTRIM1/"$i"*{r1,R1}*fastq.gz | gzip -f > rm -rf "$endir"/"$runam"/"$i"/ADAPTERTRIM/"$i"_R1.atrim.fastq.gz
			zcat rm -rf "$endir"/"$runam"/"$i"/ADAPTERTRIM1/"$i"*{r2,R2}*fastq.gz | gzip -f > rm -rf "$endir"/"$runam"/"$i"/ADAPTERTRIM/"$i"_R2.atrim.fastq.gz
			rm -rf rm -rf "$endir"/"$runam"/"$i"/ADAPTERTRIM1
		else
			echo "Uneven number of PE reads, check input data. If SE, delete R1 and R2 tags in filenames."
			exit 1
		fi
		shopt -u nullglob
	done 2> "$endir"/"$runam"/error_adapter.log
}

#trim inner adapter sequences FIN
function barcodetrim {
	list=$( cat $listf )
	echo "BARCODE trim of FASTQ files"
	if [ $erindir == TRUE ]; then
		endir=$( echo $indir )
	fi
	if ! command -v cutadapt &> /dev/null; then
		echo "Cutadapt is not installed despite bqthr is set, try sudo apt install cutadapt -y, also check dependencies"
		exit 1
	fi
	for k in $list; do
		i=$( echo $k | cut -d , -f1 )
		p5=$( echo $k | cut -d , -f2 )
		p7=$( echo $k | cut -d , -f3 )
		barcode_error=$( echo "0.$( echo $((1000*$seeds/${#p5})) )" )
		echo "Barcode trim for $i"
		rm -rf "$endir"/"$runam"/"$i"/BARCODETRIM
		mkdir "$endir"/"$runam"/"$i"/BARCODETRIM
		if [ ! -d "$endir"/"$runam"/"$i"/RUNLOGS ]; then
			mkdir  "$endir"/"$runam"/"$i"/RUNLOGS
		fi
		shopt -s nullglob
		if [ -d "$endir"/"$runam"/"$i"/ADAPTERTRIM ] && [ `ls -1 "$endir"/"$runam"/"$i"/ADAPTERTRIM/*{fastq.gz,fq.gz} 2>/dev/null | wc -l` != 0 ]; then
			cr1=$( ls -l "$endir"/"$runam"/"$i"/ADAPTERTRIM/*{r1,R1}*{fastq.gz,fq.gz} 2>dev/null | wc -l )
			cr1=$( ls -l "$endir"/"$runam"/"$i"/ADAPTERTRIM/*{r2,R2}*{fastq.gz,fq.gz} 2>dev/null | wc -l )
			if [ $cr1 == $cr2 ] && [ $cr1 != 0 ]; then #if PE
				zcat "$endir"/"$runam"/"$i"/ADAPTERTRIM/*{r1,R1}*{fastq.gz,fq.gz} | gzip -f > "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R1_input.fastq.gz
				zcat "$endir"/"$runam"/"$i"/ADAPTERTRIM/*{r2,R2}*{fastq.gz,fq.gz} | gzip -f > "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R2_input.fastq.gz
				ispe=T
			elif [ $cr1 == 0 ] && [ $cr2 == 0 ]; then #if SE
				zcat "$endir"/"$runam"/"$i"/ADAPTERTRIM/*{fastq.gz,fq.gz} | gzip -f > "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_input.fastq.gz
				ispe=F
			else #if uneven
				echo "Uneven number of PE reads, check input data. If SE, delete R1 and R2 tags in filenames."
				exit 1
			fi
		elif [ -d "$endir"/"$runam"/"$i"/BQTRIM ] && [ `ls -1 "$endir"/"$runam"/"$i"/BQTRIM/*{fastq.gz,fq.gz} 2>/dev/null | wc -l` != 0 ]; then
			cr1=$( ls -l "$endir"/"$runam"/"$i"/BQTRIM/*{r1,R1}*{fastq.gz,fq.gz} 2>dev/null | wc -l )
			cr1=$( ls -l "$endir"/"$runam"/"$i"/BQTRIM/*{r2,R2}*{fastq.gz,fq.gz} 2>dev/null | wc -l )
			if [ $cr1 == $cr2 ] && [ $cr1 != 0 ]; then #if PE
				zcat "$endir"/"$runam"/"$i"/BQTRIM/*{r1,R1}*{fastq.gz,fq.gz} | gzip -f > "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R1_input.fastq.gz
				zcat "$endir"/"$runam"/"$i"/BQTRIM/*{r2,R2}*{fastq.gz,fq.gz} | gzip -f > "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R2_input.fastq.gz
				ispe=T
			elif [ $cr1 == 0 ] && [ $cr2 == 0 ]; then #if SE
				zcat "$endir"/"$runam"/"$i"/BQTRIM/*{fastq.gz,fq.gz} | gzip -f > "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_input.fastq.gz
				ispe=F
			else #if uneven
				echo "Uneven number of PE reads, check input data. If SE, delete R1 and R2 tags in filenames."
				exit 1
			fi
		elif [ `ls -1 "$indir"/*{fastq.gz,fq.gz,fastq,fq} 2>/dev/null | wc -l` != 0 ]; then
			countr1=$( ls -1 "$indir"/"$i"*{r1,R1}.{fastq,fq,fastq.gz,fq.gz} 2>/dev/null | wc -l )
			countr2=$( ls -1 "$indir"/"$i"*{r2,R2}.{fastq,fq,fastq.gz,fq.gz} 2>/dev/null | wc -l )
			if [ $countr1 -eq 0 ] && [ $countr2 -eq 0 ]; then #if SE
				countgz=$( ls -1 "$indir"/"$i"*{fastq.gz,fq.gz} 2>/dev/null | wc -l )
				countuz=$( ls -1 "$indir"/"$i"*{fastq,fq} 2>/dev/null | wc -l )
				ispe=F
				if [ $countgz -gt 0 ] && [ $countuz -eq 0 ]; then
					zcat "$indir"/"$i"*{fastq.gz,fq.gz} | gzip -f > "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_input.fastq.gz
				elif [ $countuz -gt 0 ] && [ $countgz -eq 0 ]; then
					cat "$indir"/"$i"*{fastq,fq} | gzip -f > "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_input.fastq.gz
				else
					echo "both gzipped and unzipped files are exist for $i, these will processed, but may cause distortion in output if file duplicates are present"
					zcat "$indir"/"$i"*{fastq.gz,fq.gz} | gzip -f > "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_input1.fastq.gz
					cat "$indir"/"$i"*{fastq,fq} | gzip -f > "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_input.fastq2.gz
					zcat "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_input.fastq1.gz "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_input.fastq2.gz | gzip -f > "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_input.fastq.gz
					rm -rf "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_input.fastq1.gz
					rm -rf "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_input.fastq2.gz
				fi
			elif [ $countr1 -eq $countr2 ]; then #if PE
				countgz=$( ls -1 "$indir"/"$i"*{r1,R1}.{fastq.gz,fq.gz} 2>/dev/null | wc -l )
				countuz=$( ls -1 "$indir"/"$i"*{r1,R1}.{fastq,fq} 2>/dev/null | wc -l )
				ispe=T
				if [ $countgz -gt 0 ] && [ $countuz -eq 0 ]; then
					zcat "$indir"/"$i"*{r1,R1}.{fastq.gz,fq.gz} | gzip -f > "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R1_input.fastq.gz
					zcat "$indir"/"$i"*{r2,R2}.{fastq.gz,fq.gz} | gzip -f > "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R2_input.fastq.gz
				elif [ $countuz -gt 0 ] && [ $countgz -eq 0 ]; then
					cat "$indir"/"$i"*{r1,R1}.{fastq,fq} | gzip -f > "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R1_input.fastq.gz
					cat "$indir"/"$i"*{r2,R2}.{fastq,fq} | gzip -f > "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R2_input.fastq.gz
				else
					echo "Both gzipped and unzipped files are exist for a $i, please gzip or unzip all input"
					exit 1
				fi
			else
				echo "Uneven number of PE reads, check input data. If SE, delete R1 and R2 tags in filenames."
				exit 1
			fi
		else
			echo "FASTQ files do not exist in given directories, or are not gzipped within PAPline file structure, exiting 1"
			exit 1
		fi
		shopt -u nullglob
		if [ $bdisc == TRUE ] || [ $bdisc == true ] || [ $bdisc == T ] || [ $bdisc == t ]; then
			if [ $ispe == "T" ]; then #if PE
				if [ -z $p5 ]; then #if missing p5
					revc7=$( Rscript --vanilla ./src/reverse_complement.r $p7 )
					cutadapt \
					-j $thred \
					-a ${revc7}X \
					-e $barcode_error \
					-o "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R1_inputX.fastq.gz \
					"$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R1_input.fastq.gz > x
					rm -rf x
					rm -rf "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R1_input.fastq.gz
					cutadapt \
					-j $thred \
					--discard-untrimmed \
					-G ^$p7 \
					-m $minle \
					-e $barcode_error \
					-o "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R1_bartrim.fastq.gz \
					-p "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R2_bartrim.fastq.gz \
					"$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R1_inputX.fastq.gz \
					"$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R2_input.fastq.gz > "$endir"/"$runam"/"$i"/RUNLOGS/"$i"_barcode_trim.txt
					rm -rf "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R1_inputX.fastq.gz
					rm -rf "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R2_input.fastq.gz
				elif [ -z $p7 ]; then #if missing p7
					revc5=$( Rscript --vanilla ./src/reverse_complement.r $p5 )
					cutadapt \
					-j $thred \
					-a ${revc5}X \
					-e $barcode_error \
					-o "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R2_inputX.fastq.gz \
					"$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R2_input.fastq.gz > x
					rm -rf x
					rm -rf "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R2_input.fastq.gz
					cutadapt \
					-j $thred \
					--discard-untrimmed \
					-g ^$p5 \
					-m $minle \
					-e $barcode_error \
					-o "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R1_bartrim.fastq.gz \
					-p "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R2_bartrim.fastq.gz \
					"$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R1_input.fastq.gz \
					"$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R2_inputX.fastq.gz > "$endir"/"$runam"/"$i"/RUNLOGS/"$i"_barcode_trim.txt
					rm -rf "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R1_input.fastq.gz
					rm -rf "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R2_inputX.fastq.gz
				elif [ -z $p5 ] && [ -z $p7 ]; then #if both barcodes are missing
					mv "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R1_input.fastq.gz "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R1_bartrim.fastq.gz
					mv "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R2_input.fastq.gz "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R2_bartrim.fastq.gz
					expr $( zcat "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R1_bartrim.fastq.gz | wc -l ) / 4 &> "$endir"/"$runam"/"$i"/RUNLOGS/"$i"_barcode_trim.txt
				else #if both barcodes are given
					revc5=$( Rscript --vanilla ./src/reverse_complement.r $p5 )
					revc7=$( Rscript --vanilla ./src/reverse_complement.r $p7 )
					cutadapt \
					-j $thred \
					-a ${revc7}X \
					-A ${revc5}X \
					-e $barcode_error \
					-o "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R1_bartrimX.fastq.gz \
					-p "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R2_bartrimX.fastq.gz \
					"$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R1_input.fastq.gz \
					"$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R2_input.fastq.gz > x
					rm -rf x
					rm -rf "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R1_input.fastq.gz
					rm -rf "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R2_input.fastq.gz
					cutadapt \
					-j $thred \
					--discard-untrimmed \
					-g ^$p5 \
					-G ^$p7 \
					-m $minle \
					-e $barcode_error \
					-o "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R1_bartrim.fastq.gz \
					-p "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R2_bartrim.fastq.gz \
					"$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R1_inputX.fastq.gz \
					"$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R2_inputX.fastq.gz > "$endir"/"$runam"/"$i"/RUNLOGS/"$i"_barcode_trim.txt
					rm -rf "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R1_inputX.fastq.gz
					rm -rf "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R2_inputX.fastq.gz
				fi
			elif [ $ispe == "F" ]; then
				if [ -z $p5 ]; then #if missing p5
					echo "For SE data, P7 barcode sequence will be converted to reverse complement"
					revc7=$( Rscript --vanilla ./src/reverse_complement.r $p7 )
					cutadapt \
					-j $thred \
					--discard-untrimmed \
					-a ${revc7}X \
					-m $minle \
					-e $barcode_error \
					-o "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_bartrim.fastq.gz \
					"$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_input.fastq.gz > "$endir"/"$runam"/"$i"/RUNLOGS/"$i"_barcode_trim.txt
					rm -rf "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_input.fastq.gz
				elif [ -z $p7 ]; then #if missing p7
					cutadapt \
					-j $thred \
					--discard-untrimmed \
					-g ^$p5 \
					-m $minle \
					-e $barcode_error \
					-o "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_bartrim.fastq.gz \
					"$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_input.fastq.gz > "$endir"/"$runam"/"$i"/RUNLOGS/"$i"_barcode_trim.txt
					rm -rf "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_input.fastq.gz
				elif [ -z $p5 ] && [ -z $p7 ]; then #if both barcodes are missing
					mv "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_input.fastq.gz "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_bartrim.fastq.gz
					expr $( zcat "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_bartrim.fastq.gz | wc -l ) / 4 &> "$endir"/"$runam"/"$i"/RUNLOGS/"$i"_barcode_trim.txt
				else #if both barcodes are given
					echo "For SE data, P7 barcode sequence will be converted to reverse complement"
					revc7=$( Rscript --vanilla ./src/reverse_complement.r $p7 )
					cutadapt \
					-j $thred \
					--discard-untrimmed \
					-g ^$p5 \
					-a ${revc7}X \
					-m $minle \
					-e $barcode_error \
					-o "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_bartrim.fastq.gz \
					"$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_input.fastq.gz > "$endir"/"$runam"/"$i"/RUNLOGS/"$i"_barcode_trim.txt
					rm -rf "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_input.fastq.gz
				fi
			fi
		else
			if [ $ispe == "T" ]; then #if PE
				if [ -z $p5 ]; then #if missing p5
					revc7=$( Rscript --vanilla ./src/reverse_complement.r $p7 )
					cutadapt \
					-j $thred \
					-a ${revc7}X \
					-e $barcode_error \
					-o "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R1_inputX.fastq.gz \
					"$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R1_input.fastq.gz > x
					rm -rf x
					rm -rf "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R1_input.fastq.gz
					cutadapt \
					-j $thred \
					-G ^$p7 \
					-m $minle \
					-e $barcode_error \
					-o "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R1_bartrim.fastq.gz \
					-p "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R2_bartrim.fastq.gz \
					"$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R1_inputX.fastq.gz \
					"$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R2_input.fastq.gz > "$endir"/"$runam"/"$i"/RUNLOGS/"$i"_barcode_trim.txt
					rm -rf "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R1_inputX.fastq.gz
					rm -rf "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R2_input.fastq.gz
				elif [ -z $p7 ]; then #if missing p7
					revc5=$( Rscript --vanilla ./src/reverse_complement.r $p5 )
					cutadapt \
					-j $thred \
					-a ${revc5}X \
					-e $barcode_error \
					-o "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R2_inputX.fastq.gz \
					"$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R2_input.fastq.gz > x
					rm -rf x
					rm -rf "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R2_input.fastq.gz
					cutadapt \
					-j $thred \
					-g ^$p5 \
					-m $minle \
					-e $barcode_error \
					-o "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R1_bartrim.fastq.gz \
					-p "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R2_bartrim.fastq.gz \
					"$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R1_input.fastq.gz \
					"$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R2_inputX.fastq.gz > "$endir"/"$runam"/"$i"/RUNLOGS/"$i"_barcode_trim.txt
					rm -rf "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R1_input.fastq.gz
					rm -rf "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R2_inputX.fastq.gz
				elif [ -z $p5 ] && [ -z $p7 ]; then #if both barcodes are missing
					mv "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R1_input.fastq.gz "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R1_bartrim.fastq.gz
					mv "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R2_input.fastq.gz "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R2_bartrim.fastq.gz
					expr $( zcat "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R1_bartrim.fastq.gz | wc -l ) / 4 &> "$endir"/"$runam"/"$i"/RUNLOGS/"$i"_barcode_trim.txt
				else #if both barcodes are given
					revc5=$( Rscript --vanilla ./src/reverse_complement.r $p5 )
					revc7=$( Rscript --vanilla ./src/reverse_complement.r $p7 )
					cutadapt \
					-j $thred \
					-a ${revc7}X \
					-A ${revc5}X \
					-e $barcode_error \
					-o "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R1_bartrimX.fastq.gz \
					-p "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R2_bartrimX.fastq.gz \
					"$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R1_input.fastq.gz \
					"$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R2_input.fastq.gz > x
					rm -rf x
					rm -rf "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R1_input.fastq.gz
					rm -rf "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R2_input.fastq.gz
					cutadapt \
					-j $thred \
					-g ^$p5 \
					-G ^$p7 \
					-m $minle \
					-e $barcode_error \
					-o "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R1_bartrim.fastq.gz \
					-p "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R2_bartrim.fastq.gz \
					"$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R1_inputX.fastq.gz \
					"$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R2_inputX.fastq.gz > "$endir"/"$runam"/"$i"/RUNLOGS/"$i"_barcode_trim.txt
					rm -rf "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R1_inputX.fastq.gz
					rm -rf "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_R2_inputX.fastq.gz
				fi
			elif [ $ispe == "F" ]; then
				if [ -z $p5 ]; then #if missing p5
					echo "For SE data, P7 barcode sequence will be converted to reverse complement"
					revc7=$( Rscript --vanilla ./src/reverse_complement.r $p7 )
					cutadapt \
					-j $thred \
					-a ${revc7}X \
					-m $minle \
					-e $barcode_error \
					-o "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_bartrim.fastq.gz \
					"$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_input.fastq.gz > "$endir"/"$runam"/"$i"/RUNLOGS/"$i"_barcode_trim.txt
					rm -rf "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_input.fastq.gz
				elif [ -z $p7 ]; then #if missing p7
					cutadapt \
					-j $thred \
					-g ^$p5 \
					-m $minle \
					-e $barcode_error \
					-o "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_bartrim.fastq.gz \
					"$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_input.fastq.gz > "$endir"/"$runam"/"$i"/RUNLOGS/"$i"_barcode_trim.txt
					rm -rf "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_input.fastq.gz
				elif [ -z $p5 ] && [ -z $p7 ]; then #if both barcodes are missing
					mv "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_input.fastq.gz "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_bartrim.fastq.gz
					expr $( zcat "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_bartrim.fastq.gz | wc -l ) / 4 &> "$endir"/"$runam"/"$i"/RUNLOGS/"$i"_barcode_trim.txt
				else #if both barcodes are given
					echo "For SE data, P7 barcode sequence will be converted to reverse complement"
					revc7=$( Rscript --vanilla ./src/reverse_complement.r $p7 )
					cutadapt \
					-j $thred \
					-g ^$p5 \
					-a ${revc7}X \
					-m $minle \
					-e $barcode_error \
					-o "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_bartrim.fastq.gz \
					"$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_input.fastq.gz > "$endir"/"$runam"/"$i"/RUNLOGS/"$i"_barcode_trim.txt
					rm -rf "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"_input.fastq.gz
				fi
			fi
		fi
		echo $( awk '{for(i=1;i<=NF;i++) if ($i=="filters):") print $(i+1)}' "$endir"/"$runam"/"$i"/RUNLOGS/"$i"_barcode_trim.txt ) | tr --delete "," > "$endir"/"$runam"/"$i"/RUNLOGS/"$i".Reads_passed_barcode_trim.txt
		rm -rf "$endir"/"$runam"/"$i"/RUNLOGS/"$i"_barcode_trim.txt
	done 2> "$endir"/"$runam"/error_barcode.log
}

#merge PE reads to SE FIN
function mergeread {
	basename=$( cat $listf | cut -d, -f1 )
	echo "MERGE PE reads"
	if [ $erindir == TRUE ]; then
		endir=$( echo $indir )
	fi
	if ! command -v seqprep &> /dev/null; then
		echo "SeqPrep is not installed despite overl is set, try sudo apt install seqprep -y, also check dependencies"
		exit 1
	fi
	for i in $basename; do
		if [ ! -d "$endir"/"$runam"/"$i"/RUNLOGS ]; then
			mkdir  "$endir"/"$runam"/"$i"/RUNLOGS
		fi
		echo "Merging reads for $i"
		rm -rf "$endir"/"$runam"/"$i"/MERGE
		mkdir "$endir"/"$runam"/"$i"/MERGE
		shopt -s nullglob
		#if input from barcodetrim
		if [ -d "$endir"/"$runam"/"$i"/BARCODETRIM ] && [ `ls -1 "$endir"/"$runam"/"$i"/BARCODETRIM/*{fastq.gz,fq.gz,fastq,fq} 2>/dev/null | wc -l` != 0 ]; then
			mkdir "$endir"/"$runam"/"$i"/tmp
			temppath=$( echo "$endir"/"$runam"/"$i"/tmp )
			countr1=$( ls -1 "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"*{r1,R1}.{fastq,fq,fastq.gz,fq.gz} 2>/dev/null | wc -l )
			countr2=$( ls -1 "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"*{r2,R2}.{fastq,fq,fastq.gz,fq.gz} 2>/dev/null | wc -l )
			if [ $countr1 -eq $countr2 ]; then #if PE
				countgz=$( ls -1 "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"*{r1,R1}.{fastq.gz,fq.gz} 2>/dev/null | wc -l )
				countuz=$( ls -1 "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"*{r1,R1}.{fastq,fq} 2>/dev/null | wc -l )
				if [ $countgz -gt 0 ] && [ $countuz -eq 0 ]; then
					zcat "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"*{r1,R1}.{fastq.gz,fq.gz} | gzip -f > "$temppath"/"$i"_R1_minput.fastq.gz
					zcat "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"*{r2,R2}.{fastq.gz,fq.gz} | gzip -f > "$temppath"/"$i"_R2_minput.fastq.gz
				elif [ $countuz -gt 0 ] && [ $countgz -eq 0 ]; then
					cat "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"*{r1,R1}.{fastq,fq} | gzip -f > "$temppath"/"$i"_R1_minput.fastq.gz
					cat "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"*{r2,R2}.{fastq,fq} | gzip -f > "$temppath"/"$i"_R2_minput.fastq.gz
				else
					echo "Both gzipped and unzipped files are exist for $i, please gzip or unzip all input"
					exit 1
				fi
			elif [ $countr1 -eq 0 ] && [ $countr2 -eq 0 ]; then #if SE
				echo "No PE reads are present for $i, copying and concatenating all fastq files to MERGE directory"
				countgz=$( ls -1 "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"*{fastq.gz,fq.gz} 2>/dev/null | wc -l )
				countuz=$( ls -1 "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"*{fastq,fq} 2>/dev/null | wc -l )
				if [ $countgz -gt 0 ] && [ $countuz -eq 0 ]; then
					zcat "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"*{fastq.gz,fq.gz} | gzip -f > "$endir"/"$runam"/"$i"/MERGE/"$i".merged.fastq.gz
				elif [ $countuz -gt 0 ] && [ $countgz -eq 0 ]; then
					cat "$endir"/"$runam"/"$i"/BARCODETRIM/"$i"*{fastq,fq} | gzip -f > "$endir"/"$runam"/"$i"/MERGE/"$i".merged.fastq.gz
				else
					echo "Both gzipped and unzipped files are exist for $i, please gzip or unzip all input"
					exit 1
				fi
			else
				echo "Uneven number of PE reads for $i, check input data."
				exit 1
			fi
		#if input from adaptertrim
		elif [ -d "$endir"/"$runam"/"$i"/ADAPTERTRIM ] && [ `ls -1 "$endir"/"$runam"/"$i"/ADAPTERTRIM/*{fastq.gz,fq.gz,fastq,fq} 2>/dev/null | wc -l` != 0 ]; then
			mkdir "$endir"/"$runam"/"$i"/tmp
			temppath=$( echo "$endir"/"$runam"/"$i"/tmp )
			countr1=$( ls -1 "$endir"/"$runam"/"$i"/ADAPTERTRIM/"$i"*{r1,R1}.{fastq,fq,fastq.gz,fq.gz} 2>/dev/null | wc -l )
			countr2=$( ls -1 "$endir"/"$runam"/"$i"/ADAPTERTRIM/"$i"*{r2,R2}.{fastq,fq,fastq.gz,fq.gz} 2>/dev/null | wc -l )
			if [ $countr1 -eq $countr2 ]; then #if PE
				countgz=$( ls -1 "$endir"/"$runam"/"$i"/ADAPTERTRIM/"$i"*{r1,R1}.{fastq.gz,fq.gz} 2>/dev/null | wc -l )
				countuz=$( ls -1 "$endir"/"$runam"/"$i"/ADAPTERTRIM/"$i"*{r1,R1}.{fastq,fq} 2>/dev/null | wc -l )
				if [ $countgz -gt 0 ] && [ $countuz -eq 0 ]; then
					zcat "$endir"/"$runam"/"$i"/ADAPTERTRIM/"$i"*{r1,R1}.{fastq.gz,fq.gz} | gzip -f > "$temppath"/"$i"_R1_minput.fastq.gz
					zcat "$endir"/"$runam"/"$i"/ADAPTERTRIM/"$i"*{r2,R2}.{fastq.gz,fq.gz} | gzip -f > "$temppath"/"$i"_R2_minput.fastq.gz
				elif [ $countuz -gt 0 ] && [ $countgz -eq 0 ]; then
					cat "$endir"/"$runam"/"$i"/ADAPTERTRIM/"$i"*{r1,R1}.{fastq,fq} | gzip -f > "$temppath"/"$i"_R1_minput.fastq.gz
					cat "$endir"/"$runam"/"$i"/ADAPTERTRIM/"$i"*{r2,R2}.{fastq,fq} | gzip -f > "$temppath"/"$i"_R2_minput.fastq.gz
				else
					echo "Both gzipped and unzipped files are exist for $i, please gzip or unzip all input"
					exit 1
				fi
			elif [ $countr1 -eq 0 ] && [ $countr2 -eq 0 ]; then #if SE
				echo "No PE reads are present for $i, copying and concatenating all fastq files to MERGE directory"
				countgz=$( ls -1 "$endir"/"$runam"/"$i"/ADAPTERTRIM/"$i"*{fastq.gz,fq.gz} 2>/dev/null | wc -l )
				countuz=$( ls -1 "$endir"/"$runam"/"$i"/ADAPTERTRIM/"$i"*{fastq,fq} 2>/dev/null | wc -l )
				if [ $countgz -gt 0 ] && [ $countuz -eq 0 ]; then
					zcat "$endir"/"$runam"/"$i"/ADAPTERTRIM/"$i"*{fastq.gz,fq.gz} | gzip -f > "$endir"/"$runam"/"$i"/MERGE/"$i".merged.fastq.gz
				elif [ $countuz -gt 0 ] && [ $countgz -eq 0 ]; then
					cat "$endir"/"$runam"/"$i"/ADAPTERTRIM/"$i"*{fastq,fq} | gzip -f > "$endir"/"$runam"/"$i"/MERGE/"$i".merged.fastq.gz
				else
					echo "Both gzipped and unzipped files are exist for $i, please gzip or unzip all input"
					exit 1
				fi
			else
				echo "Uneven number of PE reads for $i, check input data."
				exit 1
			fi
		#if input from bqtrim
		elif [ -d "$endir"/"$runam"/"$i"/BQTRIM ] && [ `ls -1 "$endir"/"$runam"/"$i"/BQTRIM/*{fastq.gz,fq.gz,fastq,fq} 2>/dev/null | wc -l` != 0 ]; then
			mkdir "$endir"/"$runam"/"$i"/tmp
			temppath=$( echo "$endir"/"$runam"/"$i"/tmp )
			countr1=$( ls -1 "$endir"/"$runam"/"$i"/BQTRIM/"$i"*{r1,R1}.{fastq,fq,fastq.gz,fq.gz} 2>/dev/null | wc -l )
			countr2=$( ls -1 "$endir"/"$runam"/"$i"/BQTRIM/"$i"*{r2,R2}.{fastq,fq,fastq.gz,fq.gz} 2>/dev/null | wc -l )
			if [ $countr1 -eq $countr2 ]; then #if PE
				countgz=$( ls -1 "$endir"/"$runam"/"$i"/BQTRIM/"$i"*{r1,R1}.{fastq.gz,fq.gz} 2>/dev/null | wc -l )
				countuz=$( ls -1 "$endir"/"$runam"/"$i"/BQTRIM/"$i"*{r1,R1}.{fastq,fq} 2>/dev/null | wc -l )
				if [ $countgz -gt 0 ] && [ $countuz -eq 0 ]; then
					zcat "$endir"/"$runam"/"$i"/BQTRIM/"$i"*{r1,R1}.{fastq.gz,fq.gz} | gzip -f > "$temppath"/"$i"_R1_minput.fastq.gz
					zcat "$endir"/"$runam"/"$i"/BQTRIM/"$i"*{r2,R2}.{fastq.gz,fq.gz} | gzip -f > "$temppath"/"$i"_R2_minput.fastq.gz
				elif [ $countuz -gt 0 ] && [ $countgz -eq 0 ]; then
					cat "$endir"/"$runam"/"$i"/BQTRIM/"$i"*{r1,R1}.{fastq,fq} | gzip -f > "$temppath"/"$i"_R1_minput.fastq.gz
					cat "$endir"/"$runam"/"$i"/BQTRIM/"$i"*{r2,R2}.{fastq,fq} | gzip -f > "$temppath"/"$i"_R2_minput.fastq.gz
				else
					echo "Both gzipped and unzipped files are exist for $i, please gzip or unzip all input"
					exit 1
				fi
			elif [ $countr1 -eq 0 ] && [ $countr2 -eq 0 ]; then #if SE
				echo "No PE reads are present for $i, copying and concatenating all fastq files to MERGE directory"
				countgz=$( ls -1 "$endir"/"$runam"/"$i"/BQTRIM/"$i"*{fastq.gz,fq.gz} 2>/dev/null | wc -l )
				countuz=$( ls -1 "$endir"/"$runam"/"$i"/BQTRIM/"$i"*{fastq,fq} 2>/dev/null | wc -l )
				if [ $countgz -gt 0 ] && [ $countuz -eq 0 ]; then
					zcat "$endir"/"$runam"/"$i"/BQTRIM/"$i"*{fastq.gz,fq.gz} | gzip -f > "$endir"/"$runam"/"$i"/MERGE/"$i".merged.fastq.gz
				elif [ $countuz -gt 0 ] && [ $countgz -eq 0 ]; then
					cat "$endir"/"$runam"/"$i"/BQTRIM/"$i"*{fastq,fq} | gzip -f > "$endir"/"$runam"/"$i"/MERGE/"$i".merged.fastq.gz
				else
					echo "Both gzipped and unzipped files are exist for $i, please gzip or unzip all input"
					exit 1
				fi
			else
				echo "Uneven number of PE reads for $i, check input data."
				exit 1
			fi
		#if input directly from input directory
		elif [ `ls -1 "$indir"/*{fastq.gz,fq.gz,fastq,fq} 2>/dev/null | wc -l` != 0 ]; then
			mkdir "$endir"/"$runam"/"$i"/tmp
			temppath=$( echo "$endir"/"$runam"/"$i"/tmp )
			countr1=$( ls -1 "$indir"/"$i"*{r1,R1}.{fastq,fq,fastq.gz,fq.gz} 2>/dev/null | wc -l )
			countr2=$( ls -1 "$indir"/"$i"*{r2,R2}.{fastq,fq,fastq.gz,fq.gz} 2>/dev/null | wc -l )
			if [ $countr1 -eq $countr2 ]; then #if PE
				countgz=$( ls -1 "$indir"/"$i"*{r1,R1}.{fastq.gz,fq.gz} 2>/dev/null | wc -l )
				countuz=$( ls -1 "$indir"/"$i"*{r1,R1}.{fastq,fq} 2>/dev/null | wc -l )
				if [ $countgz -gt 0 ] && [ $countuz -eq 0 ]; then
					zcat "$indir"/"$i"*{r1,R1}.{fastq.gz,fq.gz} | gzip -f > "$temppath"/"$i"_R1_minput.fastq.gz
					zcat "$indir"/"$i"*{r2,R2}.{fastq.gz,fq.gz} | gzip -f > "$temppath"/"$i"_R2_minput.fastq.gz
				elif [ $countuz -gt 0 ] && [ $countgz -eq 0 ]; then
					cat "$indir"/"$i"*{r1,R1}.{fastq,fq} | gzip -f > "$temppath"/"$i"_R1_minput.fastq.gz
					cat "$indir"/"$i"*{r2,R2}.{fastq,fq} | gzip -f > "$temppath"/"$i"_R2_minput.fastq.gz
				else
					echo "Both gzipped and unzipped files are exist for $i, please gzip or unzip all input"
					exit 1
				fi
			elif [ $countr1 -eq 0 ] && [ $countr2 -eq 0 ]; then #if SE
				echo "No PE reads are present for $i, copying and concatenating all fastq files to MERGE directory"
				countgz=$( ls -1 "$indir"/"$i"*{fastq.gz,fq.gz} 2>/dev/null | wc -l )
				countuz=$( ls -1 "$indir"/"$i"*{fastq,fq} 2>/dev/null | wc -l )
				if [ $countgz -gt 0 ] && [ $countuz -eq 0 ]; then
					zcat "$indir"/"$i"*{fastq.gz,fq.gz} | gzip -f > "$endir"/"$runam"/"$i"/MERGE/"$i".merged.fastq.gz
				elif [ $countuz -gt 0 ] && [ $countgz -eq 0 ]; then
					cat "$indir"/"$i"*{fastq,fq} | gzip -f > "$endir"/"$runam"/"$i"/MERGE/"$i".merged.fastq.gz
				else
					echo "Both gzipped and unzipped files are exist for $i, please gzip or unzip all input"
					exit 1
				fi
			else
				echo "Uneven number of PE reads for $i, check input data."
				exit 1
			fi
		else
			echo "FASTQ files do not exist in given directories, exiting 1"
			exit 1
		fi
		seqprep \
		-f "$temppath"/*R1*.{fastq.gz,fq.gz,fastq,fq} \
		-r "$temppath"/*R2*.{fastq.gz,fq.gz,fastq,fq} \
		-1 "$endir"/"$runam"/"$i"/MERGE/"$i"_unmerged_R_1.fastq.gz \
		-2 "$endir"/"$runam"/"$i"/MERGE/"$i"_unmerged_R_2.fastq.gz \
		-L $minle \
		-s "$endir"/"$runam"/"$i"/MERGE/"$i".merged.fastq.gz \
		-o $overl 
		expr $( zcat "$endir"/"$runam"/"$i"/MERGE/"$i".merged.fastq.gz | wc -l ) / 4 &> "$endir"/"$runam"/"$i"/RUNLOGS/"$i".Reads_passed_merging.txt
		shopt -u nullglob
		rm -rf "$endir"/"$runam"/"$i"/tmp
	done | pv 2> "$endir"/"$runam"/error_merge.log
}
#minimum read length trim FIN
function minlentrim {
	basename=$( cat $listf | cut -d, -f1 )
	echo "LENGTH trim of FASTQ files"
	if [ $erindir == TRUE ]; then
		endir=$( echo $indir )
	fi
	if ! command -v cutadapt &> /dev/null; then
		echo "Cutadapt is not installed despite bqthr is set, try sudo apt install cutadapt -y, also check dependencies"
		exit 1
	fi
	shopt -s nullglob
	for i in $basename; do	
		rm -rf "$endir"/"$runam"/"$i"/LENGTHTRIM
		mkdir "$endir"/"$runam"/"$i"/LENGTHTRIM
		if [ ! -d "$endir"/"$runam"/"$i"/RUNLOGS ]; then
			mkdir  "$endir"/"$runam"/"$i"/RUNLOGS
		fi
		if [ -d "$endir"/"$runam"/"$i"/MERGE ] && [ `ls -1 "$endir"/"$runam"/"$i"/MERGE/*{fastq.gz,fq.gz,fastq,fq} 2>/dev/null | wc -l` != 0 ]; then
			temppath=$( echo "$endir"/"$runam"/"$i"/MERGE )
		elif [ -d "$endir"/"$runam"/"$i"/BARCODETRIM ] && [ `ls -1 "$endir"/"$runam"/"$i"/BARCODETRIM/*{fastq.gz,fq.gz,fastq,fq} 2>/dev/null | wc -l` != 0 ]; then
			temppath=$( echo "$endir"/"$runam"/"$i"/BARCODETRIM )
		elif [ -d "$endir"/"$runam"/"$i"/ADAPTERTRIM ] && [ `ls -1 "$endir"/"$runam"/"$i"/ADAPTERTRIM/*{fastq.gz,fq.gz,fastq,fq} 2>/dev/null | wc -l` != 0 ]; then
			temppath=$( echo "$endir"/"$runam"/"$i"/ADAPTERTRIM )
		elif [ -d "$endir"/"$runam"/"$i"/BQTRIM ] && [ `ls -1 "$endir"/"$runam"/"$i"/BQTRIM/*{fastq.gz,fq.gz,fastq,fq} 2>/dev/null | wc -l` != 0 ]; then
			temppath=$( echo "$endir"/"$runam"/"$i"/BQTRIM )
		elif [ `ls -1 "$indir"/*{fastq.gz,fq.gz,fastq,fq} 2>/dev/null | wc -l` != 0 ]; then
			temppath="$indir"
		else
			echo "FASTQ files do not exist in given directories, exiting 1"
			exit 1
		fi
		for j in "$temppath"/"$i"*.{fastq,fastq.gz,fq,fq.gz}; do
			echo "$j"
			file=$( echo ${j##*/} )
			bbn=$( echo "$file" | cut -d'.' -f1 )
			cutadapt -j $thred -q $bqthr -o "$endir"/"$runam"/"$i"/LENGTHTRIM/"$bbn"_"$minle"minlength.fastq.gz $j > "$endir"/"$runam"/"$i"/RUNLOGS/"$bbn".minlen.full.txt
			echo $( awk '{for(i=1;i<=NF;i++) if ($i=="filters):") print $(i+2)}' "$endir"/"$runam"/"$i"/RUNLOGS/"$bbn".minlen.full.txt ) | tr --delete "," > "$endir"/"$runam"/"$i"/RUNLOGS/"$bbn".passedminlen.txt
			expr $( zcat "$endir"/"$runam"/"$i"/LENGTHTRIM/"$bbn".bq"$bqthr".fastq.gz | wc -l ) / 4 &> "$endir"/"$runam"/"$i"/RUNLOGS/"$bbn".readno1.txt
		done
		array1=($(cat "$endir"/"$runam"/"$i"/RUNLOGS/*R1*readno1.txt))
		sum1=$(IFS=+; echo "$((${array1[*]}))")
		echo $sum1 > "$endir"/"$runam"/"$i"/RUNLOGS/"$i".Starting_read_number.txt
		array2=($(cat "$endir"/"$runam"/"$i"/RUNLOGS/*passedminlen.txt))
		sum2=$(IFS=+; echo "$((${array2[*]}))")
		echo $sum2 > "$endir"/"$runam"/"$i"/RUNLOGS/"$i".Minlength_passed_reads.txt
		rm -rf "$endir"/"$runam"/"$i"/RUNLOGS/*.passedminlen.txt
		rm -rf "$endir"/"$runam"/"$i"/RUNLOGS/*.readno1.txt
		rm -rf "$endir"/"$runam"/"$i"/RUNLOGS/*.minlen.full.txt
	done 2> "$endir"/"$runam"/error_minlength.log
	shopt -u nullglob
}

#mapping to ref FIN
function mapping {
	basename=$( cat $listf | cut -d, -f1 )
	echo "MAPPING of reads to reference"
	if [ $erindir == TRUE ]; then
		endir=$( echo $indir )
	fi
	shopt -s nullglob
	befer=$( echo "${refer%.*}" )
	for i in $basename; do
		rm -rf "$endir"/"$runam"/"$i"/ALNFILES
		mkdir "$endir"/"$runam"/"$i"/ALNFILES
		if [ ! -d "$endir"/"$runam"/"$i"/RUNLOGS ]; then
			mkdir  "$endir"/"$runam"/"$i"/RUNLOGS
		fi
		echo "Mapping $i to reference"
		#set input dir
		if [ -d "$endir"/"$runam"/"$i"/LENGTHTRIM ] && [ `ls -1 "$endir"/"$runam"/"$i"/LENGTHTRIM/*{fastq.gz,fq.gz,fastq,fq} 2>/dev/null | wc -l` != 0 ]; then
			temppath1=$( echo "$endir"/"$runam"/"$i"/MERGE )
		elif [ -d "$endir"/"$runam"/"$i"/MERGE ] && [ `ls -1 "$endir"/"$runam"/"$i"/MERGE/*{fastq.gz,fq.gz,fastq,fq} 2>/dev/null | wc -l` != 0 ]; then
			temppath1=$( echo "$endir"/"$runam"/"$i"/MERGE )
		elif [ -d "$endir"/"$runam"/"$i"/BARCODETRIM ] && [ `ls -1 "$endir"/"$runam"/"$i"/BARCODETRIM/*{fastq.gz,fq.gz,fastq,fq} 2>/dev/null | wc -l` != 0 ]; then
			temppath1=$( echo "$endir"/"$runam"/"$i"/BARCODETRIM )
		elif [ -d "$endir"/"$runam"/"$i"/ADAPTERTRIM ] && [ `ls -1 "$endir"/"$runam"/"$i"/ADAPTERTRIM/*{fastq.gz,fq.gz,fastq,fq} 2>/dev/null | wc -l` != 0 ]; then
			temppath1=$( echo "$endir"/"$runam"/"$i"/ADAPTERTRIM )
		elif [ -d "$endir"/"$runam"/"$i"/BQTRIM ] && [ `ls -1 "$endir"/"$runam"/"$i"/BQTRIM/*{fastq.gz,fq.gz,fastq,fq} 2>/dev/null | wc -l` != 0 ]; then
			temppath1=$( echo "$endir"/"$runam"/"$i"/BQTRIM )
		elif [ `ls -1 "$indir"/*{fastq.gz,fq.gz,fastq,fq} 2>/dev/null | wc -l` != 0 ]; then
			temppath1="$indir"
		else
			echo "FASTQ files do not exist in given directories, exiting 1"
			exit 1
		fi
		#prepare input data
		mkdir "$endir"/"$runam"/"$i"/tmp
		temppath=$( echo "$endir"/"$runam"/"$i"/tmp )
		countr1=$( ls -1 "$indir"/"$i"*{r1,R1}.{fastq,fq,fastq.gz,fq.gz} 2>/dev/null | wc -l )
		countr2=$( ls -1 "$indir"/"$i"*{r2,R2}.{fastq,fq,fastq.gz,fq.gz} 2>/dev/null | wc -l )
		if [ $countr1 -eq $countr2 ]; then #if PE
			countgz=$( ls -1 "$temppath1"/"$i"*{r1,R1}.{fastq.gz,fq.gz} 2>/dev/null | wc -l )
			countuz=$( ls -1 "$temppath1"/"$i"*{r1,R1}.{fastq,fq} 2>/dev/null | wc -l )
			ispe=T
			if [ $countgz -gt 0 ] && [ $countuz -eq 0 ]; then
				zcat "$temppath1"/"$i"*{r1,R1}.{fastq.gz,fq.gz} | gzip -f > "$temppath"/"$i"_R1_input.fastq.gz
				zcat "$temppath1"/"$i"*{r2,R2}.{fastq.gz,fq.gz} | gzip -f > "$temppath"/"$i"_R2_input.fastq.gz
			elif [ $countuz -gt 0 ] && [ $countgz -eq 0 ]; then
				cat "$temppath1"/"$i"*{r1,R1}.{fastq,fq} | gzip -f > "$temppath"/"$i"_R1_input.fastq.gz
				cat "$temppath1"/"$i"*{r2,R2}.{fastq,fq} | gzip -f > "$temppath"/"$i"_R2_input.fastq.gz
			else
				echo "Both gzipped and unzipped files are exist for $i, please gzip or unzip all input"
				exit 1
			fi
		elif [ $countr1 -eq 0 ] && [ $countr2 -eq 0 ]; then #if SE
			countgz=$( ls -1 "$temppath1"/"$i"*{fastq.gz,fq.gz} 2>/dev/null | wc -l )
			countuz=$( ls -1 "$temppath1"/"$i"*{fastq,fq} 2>/dev/null | wc -l )
			ispe=F
			if [ $countgz -gt 0 ] && [ $countuz -eq 0 ]; then
				zcat "$temppath1"/"$i"*{fastq.gz,fq.gz} | gzip -f > "$temppath"/"$i"_input.fastq.gz
			elif [ $countuz -gt 0 ] && [ $countgz -eq 0 ]; then
				cat "$temppath1"/"$i"*{fastq,fq} | gzip -f > "$temppath"/"$i"_input.fastq.gz
			else
				echo "Both gzipped and unzipped files are exist for $i, please gzip or unzip all input"
				exit 1
			fi
		else
			echo "Uneven number of PE reads for $i, check input data."
			exit 1
		fi
		#run mapping
		if [ $ispe == "F" ]; then #if SE
			if [ $align -eq BWA_ALN ]; then
				if ! command -v bwa &> /dev/null; then
					echo "BWA is not installed despite bqthr is set, try sudo apt install bwa -y, also check dependencies"
					exit 1
				fi
				bwa aln -t $thred -k $seeds $refer "$temppath"/"$i"*{fastq.gz,fq.gz,fastq,fq} > "$endir"/"$runam"/"$i"/ALNFILES/"$i".bwa_aln_se.sai
				bwa samse $refer "$endir"/"$runam"/"$i"/ALNFILES/"$i".bwa_aln_se.sai "$temppath"/"$i"*{fastq.gz,fq.gz,fastq,fq} > "$endir"/"$runam"/"$i"/ALNFILES/"$i".bwa_aln_se.sam
				samtools view -b "$endir"/"$runam"/"$i"/ALNFILES/"$i".bwa_aln_se.sam > "$endir"/"$runam"/"$i"/ALNFILES/"$i".bwa_aln_se.bam
				samtools sort -o "$endir"/"$runam"/"$i"/ALNFILES/"$i".bwa_aln_se.sorted.bam "$endir"/"$runam"/"$i"/ALNFILES/"$i".bwa_aln_se.bam
				if [ $krbam == FALSE ] || [ $krbam == false ] || [ $krbam == F ] || [ $krbam == f ]; then
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bwa_aln_se.bam
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bwa_aln_se.sam
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bwa_aln_se.sai
				else
					mkdir "$endir"/"$runam"/"$i"/ALNFILES/RAW
					mv "$endir"/"$runam"/"$i"/ALNFILES/"$i".bwa_aln_se.bam "$endir"/"$runam"/"$i"/ALNFILES/RAW/"$i".bwa_aln_se.raw.bam
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bwa_aln_se.sam
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bwa_aln_se.sai
				fi
				if [ $rmdup == FALSE ] || [ $rmdup == false ] || [ $rmdup == F ] || [ $rmdup == f ]; then
					samtools index "$endir"/"$runam"/"$i"/ALNFILES/"$i".bwa_aln_se.sorted.bam
				fi
				samtools view -F 0x904 "$endir"/"$runam"/"$i"/ALNFILES/"$i".bwa_aln_se.sorted.bam | cut -f 1 | sort | uniq | wc -l > "$endir"/"$runam"/"$i"/RUNLOGS/"$i".Reads_mapped.txt
			elif [ $align -eq BWA_MEM ]; then
				if ! command -v bwa &> /dev/null; then
					echo "BWA is not installed despite bqthr is set, try sudo apt install bwa -y, also check dependencies"
					exit 1
				fi
				bwa mem -t $thred $refer "$temppath"/"$i"*{fastq.gz,fq.gz,fastq,fq} > "$endir"/"$runam"/"$i"/ALNFILES/"$i".bwa_mem_se.sam
				samtools view -b "$endir"/"$runam"/"$i"/ALNFILES/"$i".bwa_mem_se.sam > "$endir"/"$runam"/"$i"/ALNFILES/"$i".bwa_mem_se.bam
				samtools sort -o "$endir"/"$runam"/"$i"/ALNFILES/"$i".bwa_mem_se.sorted.bam "$endir"/"$runam"/"$i"/ALNFILES/"$i".bwa_mem_se.bam
				if [ $krbam == FALSE ] || [ $krbam == false ] || [ $krbam == F ] || [ $krbam == f ]; then
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bwa_mem_se.bam
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bwa_mem_se.sam
				else
					mkdir "$endir"/"$runam"/"$i"/ALNFILES/RAW
					mv "$endir"/"$runam"/"$i"/ALNFILES/"$i".bwa_mem_se.bam "$endir"/"$runam"/"$i"/ALNFILES/RAW/"$i".bwa_mem_se.raw.bam
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bwa_mem_se.sam
				fi
				if [ $rmdup == FALSE ] || [ $rmdup == false ] || [ $rmdup == F ] || [ $rmdup == f ]; then
					samtools index "$endir"/"$runam"/"$i"/ALNFILES/"$i".bwa_mem_se.sorted.bam
				fi
				samtools view -F 0x904 "$endir"/"$runam"/"$i"/ALNFILES/"$i".bwa_mem_se.sorted.bam | cut -f 1 | sort | uniq | wc -l > "$endir"/"$runam"/"$i"/RUNLOGS/"$i".Reads_mapped.txt
			elif [ $align -eq Bowtie2_EVF ]; then
				if ! command -v bowtie2 &> /dev/null; then
					echo "Bowtie2 is not installed despite bqthr is set, try sudo apt install bowtie2 -y, also check dependencies"
					exit 1
				fi
				bowtie2 --end-to-end --very-fast --un-gz "$endir"/"$runam"/"$i"/ALNFILES/"$i".unmapped_reads.fastq.gz -t --met-file "$endir"/"$runam"/"$i"/ALNFILES/"$i".metrics.log -p $thred -x $befer -U "$temppath"/"$i"*{fastq.gz,fq.gz,fastq,fq} -S "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_evf_se.sam > "$endir"/"$i"/ALNFILES/"$i".bt2_evf_se.log
				samtools view -b "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_evf_se.sam > "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_evf_se.bam
				samtools sort -o "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_evf_se.sorted.bam "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_evf_se.bam
				if [ $krbam == FALSE ] || [ $krbam == false ] || [ $krbam == F ] || [ $krbam == f ]; then
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_evf_se.bam
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_evf_se.sam
				else
					mkdir "$endir"/"$runam"/"$i"/ALNFILES/RAW
					mv "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_evf_se.bam "$endir"/"$runam"/"$i"/ALNFILES/RAW/"$i".bt2_evf_se.raw.bam
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_evf_se.sam
				fi
				if [ $rmdup == FALSE ] || [ $rmdup == false ] || [ $rmdup == F ] || [ $rmdup == f ]; then
					samtools index "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_evf_se.sorted.bam
				fi
				samtools view -F 0x904 "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_evf_se.sorted.bam | cut -f 1 | sort | uniq | wc -l > "$endir"/"$runam"/"$i"/RUNLOGS/"$i".Reads_mapped.txt
			elif [ $align -eq Bowtie2_EF ]; then
				if ! command -v bowtie2 &> /dev/null; then
					echo "Bowtie2 is not installed despite bqthr is set, try sudo apt install bowtie2 -y, also check dependencies"
					exit 1
				fi
				bowtie2 --end-to-end --fast --un-gz "$endir"/"$runam"/"$i"/ALNFILES/"$i".unmapped_reads.fastq.gz -t --met-file "$endir"/"$runam"/"$i"/ALNFILES/"$i".metrics.log -p $thred -x $befer -U "$temppath"/"$i"*{fastq.gz,fq.gz,fastq,fq} -S "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_ef_se.sam > "$endir"/"$i"/ALNFILES/"$i".bt2_ef_se.log
				samtools view -b "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_ef_se.sam > "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_ef_se.bam
				samtools sort -o "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_ef_se.sorted.bam "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_ef_se.bam
				if [ $krbam == FALSE ] || [ $krbam == false ] || [ $krbam == F ] || [ $krbam == f ]; then
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_ef_se.bam
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_ef_se.sam
				else
					mkdir "$endir"/"$runam"/"$i"/ALNFILES/RAW
					mv "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_ef_se.bam "$endir"/"$runam"/"$i"/ALNFILES/RAW/"$i".bt2_ef_se.raw.bam
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_ef_se.sam
				fi
				if [ $rmdup == FALSE ] || [ $rmdup == false ] || [ $rmdup == F ] || [ $rmdup == f ]; then
					samtools index "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_ef_se.sorted.bam
				fi
				samtools view -F 0x904 "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_ef_se.sorted.bam | cut -f 1 | sort | uniq | wc -l > "$endir"/"$runam"/"$i"/RUNLOGS/"$i".Reads_mapped.txt
			elif [ $align -eq Bowtie2_EVS ]; then
				if ! command -v bowtie2 &> /dev/null; then
					echo "Bowtie2 is not installed despite bqthr is set, try sudo apt install bowtie2 -y, also check dependencies"
					exit 1
				fi
				bowtie2 --end-to-end --very-sensitive --un-gz "$endir"/"$runam"/"$i"/ALNFILES/"$i".unmapped_reads.fastq.gz -t --met-file "$endir"/"$runam"/"$i"/ALNFILES/"$i".metrics.log -p $thred -x $befer -U "$temppath"/"$i"*{fastq.gz,fq.gz,fastq,fq} -S "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_evs_se.sam > "$endir"/"$i"/ALNFILES/"$i".bt2_evs_se.log
				samtools view -b "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_evs_se.sam > "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_evs_se.bam
				samtools sort -o "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_evs_se.sorted.bam "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_evs_se.bam
				if [ $krbam == FALSE ] || [ $krbam == false ] || [ $krbam == F ] || [ $krbam == f ]; then
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_evs_se.bam
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_evs_se.sam
				else
					mkdir "$endir"/"$runam"/"$i"/ALNFILES/RAW
					mv "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_evs_se.bam "$endir"/"$runam"/"$i"/ALNFILES/RAW/"$i".bt2_evs_se.raw.bam
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_evs_se.sam
				fi
				if [ $rmdup == FALSE ] || [ $rmdup == false ] || [ $rmdup == F ] || [ $rmdup == f ]; then
					samtools index "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_evs_se.sorted.bam
				fi
				samtools view -F 0x904 "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_evs_se.sorted.bam | cut -f 1 | sort | uniq | wc -l > "$endir"/"$runam"/"$i"/RUNLOGS/"$i".Reads_mapped.txt
			elif [ $align -eq Bowtie2_ES ]; then
				if ! command -v bowtie2 &> /dev/null; then
					echo "Bowtie2 is not installed despite bqthr is set, try sudo apt install bowtie2 -y, also check dependencies"
					exit 1
				fi
				bowtie2 --end-to-end --sensitive --un-gz "$endir"/"$runam"/"$i"/ALNFILES/"$i".unmapped_reads.fastq.gz -t --met-file "$endir"/"$runam"/"$i"/ALNFILES/"$i".metrics.log -p $thred -x $befer -U "$temppath"/"$i"*{fastq.gz,fq.gz,fastq,fq} -S "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_es_se.sam > "$endir"/"$i"/ALNFILES/"$i".bt2_es_se.log
				samtools view -b "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_es_se.sam > "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_es_se.bam
				samtools sort -o "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_es_se.sorted.bam "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_es_se.bam
				if [ $krbam == FALSE ] || [ $krbam == false ] || [ $krbam == F ] || [ $krbam == f ]; then
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_es_se.bam
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_es_se.sam
				else
					mkdir "$endir"/"$runam"/"$i"/ALNFILES/RAW
					mv "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_es_se.bam "$endir"/"$runam"/"$i"/ALNFILES/RAW/"$i".bt2_es_se.raw.bam
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_es_se.sam
				fi
				if [ $rmdup == FALSE ] || [ $rmdup == false ] || [ $rmdup == F ] || [ $rmdup == f ]; then
					samtools index "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_es_se.sorted.bam
				fi
				samtools view -F 0x904 "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_es_se.sorted.bam | cut -f 1 | sort | uniq | wc -l > "$endir"/"$runam"/"$i"/RUNLOGS/"$i".Reads_mapped.txt
			elif [ $align -eq Bowtie2_LVF ]; then
				if ! command -v bowtie2 &> /dev/null; then
					echo "Bowtie2 is not installed despite bqthr is set, try sudo apt install bowtie2 -y, also check dependencies"
					exit 1
				fi
				bowtie2 --local --very-fast-local --un-gz "$endir"/"$runam"/"$i"/ALNFILES/"$i".unmapped_reads.fastq.gz -t --met-file "$endir"/"$runam"/"$i"/ALNFILES/"$i".metrics.log -p $thred -x $befer -U "$temppath"/"$i"*{fastq.gz,fq.gz,fastq,fq} -S "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lvf_se.sam > "$endir"/"$i"/ALNFILES/"$i".bt2_lvf_se.log
				samtools view -b "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lvf_se.sam > "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lvf_se.bam
				samtools sort -o "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lvf_se.sorted.bam "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lvf_se.bam
				if [ $krbam == FALSE ] || [ $krbam == false ] || [ $krbam == F ] || [ $krbam == f ]; then
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lvf_se.bam
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lvf_se.sam
				else
					mkdir "$endir"/"$runam"/"$i"/ALNFILES/RAW
					mv "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lvf_se.bam "$endir"/"$runam"/"$i"/ALNFILES/RAW/"$i".bt2_lvf_se.raw.bam
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lvf_se.sam
				fi
				if [ $rmdup == FALSE ] || [ $rmdup == false ] || [ $rmdup == F ] || [ $rmdup == f ]; then
					samtools index "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lvf_se.sorted.bam
				fi
				samtools view -F 0x904 "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lvf_se.sorted.bam | cut -f 1 | sort | uniq | wc -l > "$endir"/"$runam"/"$i"/RUNLOGS/"$i".Reads_mapped.txt
			elif [ $align -eq Bowtie2_LF ]; then
				if ! command -v bowtie2 &> /dev/null; then
					echo "Bowtie2 is not installed despite bqthr is set, try sudo apt install bowtie2 -y, also check dependencies"
					exit 1
				fi
				bowtie2 --local --fast-local --un-gz "$endir"/"$runam"/"$i"/ALNFILES/"$i".unmapped_reads.fastq.gz -t --met-file "$endir"/"$runam"/"$i"/ALNFILES/"$i".metrics.log -p $thred -x $befer -U "$temppath"/"$i"*{fastq.gz,fq.gz,fastq,fq} -S "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lf_se.sam > "$endir"/"$i"/ALNFILES/"$i".bt2_lf_se.log
				samtools view -b "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lf_se.sam > "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lf_se.bam
				samtools sort -o "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lf_se.sorted.bam "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lf_se.bam
				if [ $krbam == FALSE ] || [ $krbam == false ] || [ $krbam == F ] || [ $krbam == f ]; then
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lf_se.bam
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lf_se.sam
				else
					mkdir "$endir"/"$runam"/"$i"/ALNFILES/RAW
					mv "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lf_se.bam "$endir"/"$runam"/"$i"/ALNFILES/RAW/"$i".bt2_lf_se.raw.bam
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lf_se.sam
				fi
				if [ $rmdup == FALSE ] || [ $rmdup == false ] || [ $rmdup == F ] || [ $rmdup == f ]; then
					samtools index "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lf_se.sorted.bam
				fi
				samtools view -F 0x904 "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lf_se.sorted.bam | cut -f 1 | sort | uniq | wc -l > "$endir"/"$runam"/"$i"/RUNLOGS/"$i".Reads_mapped.txt
			elif [ $align -eq Bowtie2_LVS ]; then
				if ! command -v bowtie2 &> /dev/null; then
					echo "Bowtie2 is not installed despite bqthr is set, try sudo apt install bowtie2 -y, also check dependencies"
					exit 1
				fi
				bowtie2 --local --very-sensitive-local --un-gz "$endir"/"$runam"/"$i"/ALNFILES/"$i".unmapped_reads.fastq.gz -t --met-file "$endir"/"$runam"/"$i"/ALNFILES/"$i".metrics.log -p $thred -x $befer -U "$temppath"/"$i"*{fastq.gz,fq.gz,fastq,fq} -S "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lvs_se.sam > "$endir"/"$i"/ALNFILES/"$i".bt2_lvs_se.log
				samtools view -b "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lvs_se.sam > "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lvs_se.bam
				samtools sort -o "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lvs_se.sorted.bam "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lvs_se.bam
				if [ $krbam == FALSE ] || [ $krbam == false ] || [ $krbam == F ] || [ $krbam == f ]; then
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lvs_se.bam
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lvs_se.sam
				else
					mkdir "$endir"/"$runam"/"$i"/ALNFILES/RAW
					mv "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lvs_se.bam "$endir"/"$runam"/"$i"/ALNFILES/RAW/"$i".bt2_lvs_se.raw.bam
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lvs_se.sam
				fi
				if [ $rmdup == FALSE ] || [ $rmdup == false ] || [ $rmdup == F ] || [ $rmdup == f ]; then
					samtools index "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lvs_se.sorted.bam
				fi
				samtools view -F 0x904 "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lvs_se.sorted.bam | cut -f 1 | sort | uniq | wc -l > "$endir"/"$runam"/"$i"/RUNLOGS/"$i".Reads_mapped.txt
			elif [ $align -eq Bowtie2_LS ]; then
				if ! command -v bowtie2 &> /dev/null; then
					echo "Bowtie2 is not installed despite bqthr is set, try sudo apt install bowtie2 -y, also check dependencies"
					exit 1
				fi
				bowtie2 --local --sensitive-local --un-gz "$endir"/"$runam"/"$i"/ALNFILES/"$i".unmapped_reads.fastq.gz -t --met-file "$endir"/"$runam"/"$i"/ALNFILES/"$i".metrics.log -p $thred -x $befer -U "$temppath"/"$i"*{fastq.gz,fq.gz,fastq,fq} -S "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_ls_se.sam > "$endir"/"$i"/ALNFILES/"$i".bt2_ls_se.log
				samtools view -b "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_ls_se.sam > "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_ls_se.bam
				samtools sort -o "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_ls_se.sorted.bam "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_ls_se.bam
				if [ $krbam == FALSE ] || [ $krbam == false ] || [ $krbam == F ] || [ $krbam == f ]; then
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_ls_se.bam
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_ls_se.sam
				else
					mkdir "$endir"/"$runam"/"$i"/ALNFILES/RAW
					mv "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_ls_se.bam "$endir"/"$runam"/"$i"/ALNFILES/RAW/"$i".bt2_ls_se.raw.bam
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_ls_se.sam
				fi
				if [ $rmdup == FALSE ] || [ $rmdup == false ] || [ $rmdup == F ] || [ $rmdup == f ]; then
					samtools index "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_ls_se.sorted.bam
				fi
				samtools view -F 0x904 "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_ls_se.sorted.bam | cut -f 1 | sort | uniq | wc -l > "$endir"/"$runam"/"$i"/RUNLOGS/"$i".Reads_mapped.txt
			else
				echo "align option was incorrectly set, check PARAMETERFILE"
				exit 1
			fi
		elif [ $ispe == "T" ]; then #if PE
			if [ $align -eq BWA_ALN ]; then
				if ! command -v bwa &> /dev/null; then
					echo "BWA is not installed despite bqthr is set, try sudo apt install bwa -y, also check dependencies"
					exit 1
				fi
				bwa aln -t $thred -k $seeds $refer "$temppath"/"$i"*{r1,R1}*{fastq.gz,fq.gz,fastq,fq} > "$endir"/"$runam"/"$i"/ALNFILES/"$i".R1.bwa_aln_pe.sai
				bwa aln -t $thred -k $seeds $refer "$temppath"/"$i"*{r2,R2}*{fastq.gz,fq.gz,fastq,fq} > "$endir"/"$runam"/"$i"/ALNFILES/"$i".R2.bwa_aln_pe.sai
				bwa sampe $refer "$endir"/"$runam"/"$i"/ALNFILES/"$i".R1.bwa_aln_pe.sai "$endir"/"$runam"/"$i"/ALNFILES/"$i".R2.bwa_aln_pe.sai "$temppath"/"$i"*{r1,R1}*{fastq.gz,fq.gz,fastq,fq} "$temppath"/"$i"*{r2,R2}*{fastq.gz,fq.gz,fastq,fq} > "$endir"/"$runam"/"$i"/ALNFILES/"$i".bwa_aln_pe.sam
				samtools view -b "$endir"/"$runam"/"$i"/ALNFILES/"$i".bwa_aln_pe.sam > "$endir"/"$runam"/"$i"/ALNFILES/"$i".bwa_aln_pe.bam
				samtools sort -o "$endir"/"$runam"/"$i"/ALNFILES/"$i".bwa_aln_pe.sorted.bam "$endir"/"$runam"/"$i"/ALNFILES/"$i".bwa_aln_pe.bam
				if [ $krbam == FALSE ] || [ $krbam == false ] || [ $krbam == F ] || [ $krbam == f ]; then
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bwa_aln_pe.bam
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bwa_aln_pe.sam
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".R1.bwa_aln_pe.sai
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".R2.bwa_aln_pe.sai
				else
					mkdir "$endir"/"$runam"/"$i"/ALNFILES/RAW
					mv "$endir"/"$runam"/"$i"/ALNFILES/"$i".bwa_aln_pe.bam "$endir"/"$runam"/"$i"/ALNFILES/RAW/"$i".bwa_aln_pe.raw.bam
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bwa_aln_pe.sam
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".R1.bwa_aln_pe.sai
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".R2.bwa_aln_pe.sai
				fi
				if [ $rmdup == FALSE ] || [ $rmdup == false ] || [ $rmdup == F ] || [ $rmdup == f ]; then
					samtools index "$endir"/"$runam"/"$i"/ALNFILES/"$i".bwa_aln_pe.sorted.bam
				fi
				samtools view -F 0x904 "$endir"/"$runam"/"$i"/ALNFILES/"$i".bwa_aln_pe.sorted.bam | cut -f 1 | sort | uniq | wc -l > "$endir"/"$runam"/"$i"/RUNLOGS/"$i".Reads_mapped.txt
			elif [ $align -eq BWA_MEM ]; then
				if ! command -v bwa &> /dev/null; then
					echo "BWA is not installed despite bqthr is set, try sudo apt install bwa -y, also check dependencies"
					exit 1
				fi
				bwa mem -t $thred $refer "$temppath"/"$i"*{r1,R1}*{fastq.gz,fq.gz,fastq,fq} "$temppath"/"$i"*{r2,R2}*{fastq.gz,fq.gz,fastq,fq} > "$endir"/"$runam"/"$i"/ALNFILES/"$i".bwa_mem_pe.sam
				samtools view -b "$endir"/"$runam"/"$i"/ALNFILES/"$i".bwa_mem_pe.sam > "$endir"/"$runam"/"$i"/ALNFILES/"$i".bwa_mem_pe.bam
				samtools sort -o "$endir"/"$runam"/"$i"/ALNFILES/"$i".bwa_mem_pe.sorted.bam "$endir"/"$runam"/"$i"/ALNFILES/"$i".bwa_mem_pe.bam
				if [ $krbam == FALSE ] || [ $krbam == false ] || [ $krbam == F ] || [ $krbam == f ]; then
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bwa_mem_pe.bam
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bwa_mem_pe.sam
				else
					mkdir "$endir"/"$runam"/"$i"/ALNFILES/RAW
					mv "$endir"/"$runam"/"$i"/ALNFILES/"$i".bwa_mem_pe.bam "$endir"/"$runam"/"$i"/ALNFILES/RAW/"$i".bwa_mem_pe.raw.bam
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bwa_mem_pe.sam
				fi
				if [ $rmdup == FALSE ] || [ $rmdup == false ] || [ $rmdup == F ] || [ $rmdup == f ]; then
					samtools index "$endir"/"$runam"/"$i"/ALNFILES/"$i".bwa_mem_pe.sorted.bam
				fi
				samtools view -F 0x904 "$endir"/"$runam"/"$i"/ALNFILES/"$i".bwa_mem_pe.sorted.bam | cut -f 1 | sort | uniq | wc -l > "$endir"/"$runam"/"$i"/RUNLOGS/"$i".Reads_mapped.txt
			elif [ $align -eq Bowtie2_EVF ]; then
				if ! command -v bowtie2 &> /dev/null; then
					echo "Bowtie2 is not installed despite bqthr is set, try sudo apt install bowtie2 -y, also check dependencies"
					exit 1
				fi
				bowtie2 --end-to-end --very-fast --un-gz "$endir"/"$runam"/"$i"/ALNFILES/"$i".unmapped_reads.fastq.gz -t --met-file "$endir"/"$runam"/"$i"/ALNFILES/"$i".metrics.log -p $thred -x $befer -1 "$temppath"/"$i"*{r1,R1}*{fastq.gz,fq.gz,fastq,fq} -2 "$temppath"/"$i"*{r2,R2}*{fastq.gz,fq.gz,fastq,fq} -S "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_evf_pe.sam > "$endir"/"$i"/ALNFILES/"$i".bt2_evf_pe.log
				samtools view -b "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_evf_pe.sam > "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_evf_pe.bam
				samtools sort -o "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_evf_pe.sorted.bam "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_evf_pe.bam
				if [ $krbam == FALSE ] || [ $krbam == false ] || [ $krbam == F ] || [ $krbam == f ]; then
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_evf_pe.bam
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_evf_pe.sam
				else
					mkdir "$endir"/"$runam"/"$i"/ALNFILES/RAW
					mv "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_evf_pe.bam "$endir"/"$runam"/"$i"/ALNFILES/RAW/"$i".bt2_evf_pe.raw.bam
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_evf_pe.sam
				fi
				if [ $rmdup == FALSE ] || [ $rmdup == false ] || [ $rmdup == F ] || [ $rmdup == f ]; then
					samtools index "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_evf_pe.sorted.bam
				fi
				samtools view -F 0x904 "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_evf_pe.sorted.bam | cut -f 1 | sort | uniq | wc -l > "$endir"/"$runam"/"$i"/RUNLOGS/"$i".Reads_mapped.txt
			elif [ $align -eq Bowtie2_EF ]; then
				if ! command -v bowtie2 &> /dev/null; then
					echo "Bowtie2 is not installed despite bqthr is set, try sudo apt install bowtie2 -y, also check dependencies"
					exit 1
				fi
				bowtie2 --end-to-end --fast --un-gz "$endir"/"$runam"/"$i"/ALNFILES/"$i".unmapped_reads.fastq.gz -t --met-file "$endir"/"$runam"/"$i"/ALNFILES/"$i".metrics.log -p $thred -x $befer -1 "$temppath"/"$i"*{r1,R1}*{fastq.gz,fq.gz,fastq,fq} -2 "$temppath"/"$i"*{r2,R2}*{fastq.gz,fq.gz,fastq,fq} -S "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_ef_pe.sam > "$endir"/"$i"/ALNFILES/"$i".bt2_ef_pe.log
				samtools view -b "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_ef_pe.sam > "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_ef_pe.bam
				samtools sort -o "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_ef_pe.sorted.bam "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_ef_pe.bam
				if [ $krbam == FALSE ] || [ $krbam == false ] || [ $krbam == F ] || [ $krbam == f ]; then
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_ef_pe.bam
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_ef_pe.sam
				else
					mkdir "$endir"/"$runam"/"$i"/ALNFILES/RAW
					mv "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_ef_pe.bam "$endir"/"$runam"/"$i"/ALNFILES/RAW/"$i".bt2_ef_pe.raw.bam
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_ef_pe.sam
				fi
				if [ $rmdup == FALSE ] || [ $rmdup == false ] || [ $rmdup == F ] || [ $rmdup == f ]; then
					samtools index "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_ef_pe.sorted.bam
				fi
				samtools view -F 0x904 "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_ef_pe.sorted.bam | cut -f 1 | sort | uniq | wc -l > "$endir"/"$runam"/"$i"/RUNLOGS/"$i".Reads_mapped.txt
			elif [ $align -eq Bowtie2_EVS ]; then
				if ! command -v bowtie2 &> /dev/null; then
					echo "Bowtie2 is not installed despite bqthr is set, try sudo apt install bowtie2 -y, also check dependencies"
					exit 1
				fi
				bowtie2 --end-to-end --very-sensitive --un-gz "$endir"/"$runam"/"$i"/ALNFILES/"$i".unmapped_reads.fastq.gz -t --met-file "$endir"/"$runam"/"$i"/ALNFILES/"$i".metrics.log -p $thred -x $befer -1 "$temppath"/"$i"*{r1,R1}*{fastq.gz,fq.gz,fastq,fq} -2 "$temppath"/"$i"*{r2,R2}*{fastq.gz,fq.gz,fastq,fq} -S "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_evs_pe.sam > "$endir"/"$i"/ALNFILES/"$i".bt2_evs_pe.log
				samtools view -b "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_evs_pe.sam > "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_evs_pe.bam
				samtools sort -o "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_evs_pe.sorted.bam "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_evs_pe.bam
				if [ $krbam == FALSE ] || [ $krbam == false ] || [ $krbam == F ] || [ $krbam == f ]; then
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_evs_pe.bam
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_evs_pe.sam
				else
					mkdir "$endir"/"$runam"/"$i"/ALNFILES/RAW
					mv "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_evs_pe.bam "$endir"/"$runam"/"$i"/ALNFILES/RAW/"$i".bt2_evs_pe.raw.bam
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_evs_pe.sam
				fi
				if [ $rmdup == FALSE ] || [ $rmdup == false ] || [ $rmdup == F ] || [ $rmdup == f ]; then
					samtools index "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_evs_pe.sorted.bam
				fi
				samtools view -F 0x904 "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_evs_pe.sorted.bam | cut -f 1 | sort | uniq | wc -l > "$endir"/"$runam"/"$i"/RUNLOGS/"$i".Reads_mapped.txt
			elif [ $align -eq Bowtie2_ES ]; then
				if ! command -v bowtie2 &> /dev/null; then
					echo "Bowtie2 is not installed despite bqthr is set, try sudo apt install bowtie2 -y, also check dependencies"
					exit 1
				fi
				bowtie2 --end-to-end --sensitive --un-gz "$endir"/"$runam"/"$i"/ALNFILES/"$i".unmapped_reads.fastq.gz -t --met-file "$endir"/"$runam"/"$i"/ALNFILES/"$i".metrics.log -p $thred -x $befer -1 "$temppath"/"$i"*{r1,R1}*{fastq.gz,fq.gz,fastq,fq} -2 "$temppath"/"$i"*{r2,R2}*{fastq.gz,fq.gz,fastq,fq} -S "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_es_pe.sam > "$endir"/"$i"/ALNFILES/"$i".bt2_es_pe.log
				samtools view -b "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_es_pe.sam > "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_es_pe.bam
				samtools sort -o "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_es_pe.sorted.bam "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_es_pe.bam
				if [ $krbam == FALSE ] || [ $krbam == false ] || [ $krbam == F ] || [ $krbam == f ]; then
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_es_pe.bam
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_es_pe.sam
				else
					mkdir "$endir"/"$runam"/"$i"/ALNFILES/RAW
					mv "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_es_pe.bam "$endir"/"$runam"/"$i"/ALNFILES/RAW/"$i".bt2_es_pe.raw.bam
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_es_pe.sam
				fi
				if [ $rmdup == FALSE ] || [ $rmdup == false ] || [ $rmdup == F ] || [ $rmdup == f ]; then
					samtools index "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_es_pe.sorted.bam
				fi
				samtools view -F 0x904 "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_es_pe.sorted.bam | cut -f 1 | sort | uniq | wc -l > "$endir"/"$runam"/"$i"/RUNLOGS/"$i".Reads_mapped.txt
			elif [ $align -eq Bowtie2_LVF ]; then
				if ! command -v bowtie2 &> /dev/null; then
					echo "Bowtie2 is not installed despite bqthr is set, try sudo apt install bowtie2 -y, also check dependencies"
					exit 1
				fi
				bowtie2 --local --very-fast-local --un-gz "$endir"/"$runam"/"$i"/ALNFILES/"$i".unmapped_reads.fastq.gz -t --met-file "$endir"/"$runam"/"$i"/ALNFILES/"$i".metrics.log -p $thred -x $befer -1 "$temppath"/"$i"*{r1,R1}*{fastq.gz,fq.gz,fastq,fq} -2 "$temppath"/"$i"*{r2,R2}*{fastq.gz,fq.gz,fastq,fq} -S "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lvf_pe.sam > "$endir"/"$i"/ALNFILES/"$i".bt2_lvf_pe.log
				samtools view -b "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lvf_pe.sam > "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lvf_pe.bam
				samtools sort -o "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lvf_pe.sorted.bam "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lvf_pe.bam
				if [ $krbam == FALSE ] || [ $krbam == false ] || [ $krbam == F ] || [ $krbam == f ]; then
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lvf_pe.bam
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lvf_pe.sam
				else
					mkdir "$endir"/"$runam"/"$i"/ALNFILES/RAW
					mv "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lvf_pe.bam "$endir"/"$runam"/"$i"/ALNFILES/RAW/"$i".bt2_lvf_pe.raw.bam
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lvf_pe.sam
				fi
				if [ $rmdup == FALSE ] || [ $rmdup == false ] || [ $rmdup == F ] || [ $rmdup == f ]; then
					samtools index "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lvf_pe.sorted.bam
				fi
				samtools view -F 0x904 "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lvf_pe.sorted.bam | cut -f 1 | sort | uniq | wc -l > "$endir"/"$runam"/"$i"/RUNLOGS/"$i".Reads_mapped.txt
			elif [ $align -eq Bowtie2_LF ]; then
				if ! command -v bowtie2 &> /dev/null; then
					echo "Bowtie2 is not installed despite bqthr is set, try sudo apt install bowtie2 -y, also check dependencies"
					exit 1
				fi
				bowtie2 --local --fast-local --un-gz "$endir"/"$runam"/"$i"/ALNFILES/"$i".unmapped_reads.fastq.gz -t --met-file "$endir"/"$runam"/"$i"/ALNFILES/"$i".metrics.log -p $thred -x $befer -1 "$temppath"/"$i"*{r1,R1}*{fastq.gz,fq.gz,fastq,fq} -2 "$temppath"/"$i"*{r2,R2}*{fastq.gz,fq.gz,fastq,fq} -S "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lf_pe.sam > "$endir"/"$i"/ALNFILES/"$i".bt2_lf_pe.log
				samtools view -b "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lf_pe.sam > "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lf_pe.bam
				samtools sort -o "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lf_pe.sorted.bam "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lf_pe.bam
				if [ $krbam == FALSE ] || [ $krbam == false ] || [ $krbam == F ] || [ $krbam == f ]; then
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lf_pe.bam
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lf_pe.sam
				else
					mkdir "$endir"/"$runam"/"$i"/ALNFILES/RAW
					mv "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lf_pe.bam "$endir"/"$runam"/"$i"/ALNFILES/RAW/"$i".bt2_lf_pe.raw.bam
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lf_pe.sam
				fi
				if [ $rmdup == FALSE ] || [ $rmdup == false ] || [ $rmdup == F ] || [ $rmdup == f ]; then
					samtools index "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lf_pe.sorted.bam
				fi
				samtools view -F 0x904 "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lf_pe.sorted.bam | cut -f 1 | sort | uniq | wc -l > "$endir"/"$runam"/"$i"/RUNLOGS/"$i".Reads_mapped.txt
			elif [ $align -eq Bowtie2_LVS ]; then
				if ! command -v bowtie2 &> /dev/null; then
					echo "Bowtie2 is not installed despite bqthr is set, try sudo apt install bowtie2 -y, also check dependencies"
					exit 1
				fi
				bowtie2 --local --very-sensitive-local --un-gz "$endir"/"$runam"/"$i"/ALNFILES/"$i".unmapped_reads.fastq.gz -t --met-file "$endir"/"$runam"/"$i"/ALNFILES/"$i".metrics.log -p $thred -x $befer -1 "$temppath"/"$i"*{r1,R1}*{fastq.gz,fq.gz,fastq,fq} -2 "$temppath"/"$i"*{r2,R2}*{fastq.gz,fq.gz,fastq,fq} -S "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lvs_pe.sam > "$endir"/"$i"/ALNFILES/"$i".bt2_lvs_pe.log
				samtools view -b "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lvs_pe.sam > "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lvs_pe.bam
				samtools sort -o "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lvs_pe.sorted.bam "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lvs_pe.bam
				if [ $krbam == FALSE ] || [ $krbam == false ] || [ $krbam == F ] || [ $krbam == f ]; then
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lvs_pe.bam
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lvs_pe.sam
				else
					mkdir "$endir"/"$runam"/"$i"/ALNFILES/RAW
					mv "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lvs_pe.bam "$endir"/"$runam"/"$i"/ALNFILES/RAW/"$i".bt2_lvs_pe.raw.bam
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lvs_pe.sam
				fi
				if [ $rmdup == FALSE ] || [ $rmdup == false ] || [ $rmdup == F ] || [ $rmdup == f ]; then
					samtools index "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lvs_pe.sorted.bam
				fi
				samtools view -F 0x904 "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_lvs_pe.sorted.bam | cut -f 1 | sort | uniq | wc -l > "$endir"/"$runam"/"$i"/RUNLOGS/"$i".Reads_mapped.txt
			elif [ $align -eq Bowtie2_LS ]; then
				if ! command -v bowtie2 &> /dev/null; then
					echo "Bowtie2 is not installed despite bqthr is set, try sudo apt install bowtie2 -y, also check dependencies"
					exit 1
				fi
				bowtie2 --local --sensitive-local --un-gz "$endir"/"$runam"/"$i"/ALNFILES/"$i".unmapped_reads.fastq.gz -t --met-file "$endir"/"$runam"/"$i"/ALNFILES/"$i".metrics.log -p $thred -x $befer -1 "$temppath"/"$i"*{r1,R1}*{fastq.gz,fq.gz,fastq,fq} -2 "$temppath"/"$i"*{r2,R2}*{fastq.gz,fq.gz,fastq,fq} -S "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_ls_pe.sam > "$endir"/"$i"/ALNFILES/"$i".bt2_ls_pe.log
				samtools view -b "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_ls_pe.sam > "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_ls_pe.bam
				samtools sort -o "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_ls_pe.sorted.bam "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_ls_pe.bam
				if [ $krbam == FALSE ] || [ $krbam == false ] || [ $krbam == F ] || [ $krbam == f ]; then
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_ls_pe.bam
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_ls_pe.sam
				else
					mkdir "$endir"/"$runam"/"$i"/ALNFILES/RAW
					mv "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_ls_pe.bam "$endir"/"$runam"/"$i"/ALNFILES/RAW/"$i".bt2_ls_pe.raw.bam
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_ls_pe.sam
				fi
				if [ $rmdup == FALSE ] || [ $rmdup == false ] || [ $rmdup == F ] || [ $rmdup == f ]; then
					samtools index "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_ls_pe.sorted.bam
				fi
				samtools view -F 0x904 "$endir"/"$runam"/"$i"/ALNFILES/"$i".bt2_ls_pe.sorted.bam | cut -f 1 | sort | uniq | wc -l > "$endir"/"$runam"/"$i"/RUNLOGS/"$i".Reads_mapped.txt
			else
				echo "align option was incorrectly set, check PARAMETERFILE"
				exit 1
			fi
		else
			echo "Sumfins wrooong here"
			exit 1
		fi
		rm -rf "$endir"/"$runam"/"$i"/tmp
		if [ $presq == TRUE ] || [ $presq == true ] || [ $presq == T ] || [ $presq == t ]; then
			preseq c_curve -o "$endir"/"$runam"/"$i"/RUNLOGS/"$i".extant_complexity -B "$endir"/"$runam"/"$i"/ALNFILES/"$i"*bam #for extant complexity
			preseq lc_extrap -o "$endir"/"$runam"/"$i"/RUNLOGS/"$i".future_yield -B "$endir"/"$runam"/"$i"/ALNFILES/"$i"*bam #for future yield
		fi
	done 2> "$endir"/"$runam"/error_mapping.log
	shopt -u nullglob
}

#remove pcr duplicates FIN
function duplicate {
	basename=$( cat $listf | cut -d, -f1 )
	echo "PCR DUPLICATE removal"
	if [ $erindir == TRUE ]; then
		endir=$( echo $indir )
	fi
	shopt -s nullglob
	for i in $basename; do
		rm -rf "$endir"/"$runam"/"$i"/ALNFILES/DEDUP
		mkdir "$endir"/"$runam"/"$i"/ALNFILES/DEDUP
		if [ ! -d "$endir"/"$runam"/"$i"/RUNLOGS ]; then
			mkdir  "$endir"/"$runam"/"$i"/RUNLOGS
		fi
		echo "Mapping $i to reference"
		#set input dir
		if [ -d "$endir"/"$runam"/"$i"/ALNFILES ] && [ `ls -1 "$endir"/"$runam"/"$i"/ALNFILES/*bam 2>/dev/null | wc -l` != 0 ]; then
			temppath1=$( echo "$endir"/"$runam"/"$i"/ALNFILES )
			sortede="T"
		elif [ `ls -1 "$indir"/*bam 2>/dev/null | wc -l` != 0 ]; then
			temppath1="$indir"
		else
			echo "BAM file does not exist in given directory, exiting 1"
			exit 1
		fi
		mkdir "$endir"/"$runam"/"$i"/tmp
		temppath=$( echo "$endir"/"$runam"/"$i"/tmp )
		countbam=$( ls -1 "$temppath1"/"$i"*bam 2>/dev/null | wc -l )
		if [ $countbam -gt 1 ]; then
			samtools merge "$temppath"/"$i".merged.bam "$temppath1"/"$i"*bam
			samtools sort -o "$temppath"/"$i".merged.sorted.bam "$temppath"/"$i".merged.bam
			rm -rf "$temppath"/"$i".merged.bam
			tmpset="T"
		elif [ $countbam -eq 1 ]; then
			if [ $sortede == "T" ]; then
				temppath=$( echo "$temppath1" )
			else
				mv "$temppath1"/"$i"*bam "$temppath"/"$i".bam
				samtools sort -o "$temppath"/"$i".sorted.bam "$temppath"/"$i".bam
				rm -rf "$temppath"/"$i".bam
				tmpset="T"
			fi
		fi
		#run stuff
		samtools view -F 0x904 "$temppath"/"$i"*bam | cut -f 1 | sort | uniq | wc -l > "$endir"/"$runam"/"$i"/RUNLOGS/"$i".Reads_mapped.txt
		if [ $presq == TRUE ] || [ $presq == true ] || [ $presq == T ] || [ $presq == t ]; then
			if [ ! -f "$endir"/"$runam"/"$i"/RUNLOGS/*extant_complexity ] && [ ! -f "$endir"/"$runam"/"$i"/RUNLOGS/*future_yield ]; then
				preseq c_curve -o "$endir"/"$runam"/"$i"/RUNLOGS/"$i".extant_complexity -B "$temppath"/*bam #for extant complexity
				preseq lc_extrap -o "$endir"/"$runam"/"$i"/RUNLOGS/"$i".future_yield -B "$temppath"/*bam #for future yield
			fi
		fi
		samtools rmdup -S "$temppath"/"$i"*bam "$endir"/"$runam"/"$i"/ALNFILES/DEDUP/"$i".dedup.bam
		samtools index "$endir"/"$runam"/"$i"/ALNFILES/DEDUP/"$i".dedup.bam
		samtools view -F 0x904 "$endir"/"$runam"/"$i"/ALNFILES/DEDUP/"$i".dedup.bam | cut -f 1 | sort | uniq | wc -l > "$endir"/"$runam"/"$i"/RUNLOGS/"$i".Reads_mapped_dedup.txt
		if [ $krbam == FALSE ] || [ $krbam == false ] || [ $krbam == F ] || [ $krbam == f ]; then
			rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i"*bam
			rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i"*bai
		else
			if [ ! -d "$endir"/"$runam"/"$i"/ALNFILES/RAW ]; then
				mkdir "$endir"/"$runam"/"$i"/ALNFILES/RAW
			fi
			if [ $tmpset == "T" ]; then
				mv "$endir"/"$runam"/"$i"/tmp/"$i"*sorted.bam "$endir"/"$runam"/"$i"/ALNFILES/RAW/"$i".sorted.bam
			else
				cp "$endir"/"$runam"/"$i"/ALNFILES/"$i"*bam "$endir"/"$runam"/"$i"/ALNFILES/RAW/
				cp "$endir"/"$runam"/"$i"/ALNFILES/"$i"*bai "$endir"/"$runam"/"$i"/ALNFILES/RAW/
				rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i"*bam
				rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i"*bai
			fi
		fi
		rm -rf "$endir"/"$runam"/"$i"/tmp
	done 2> "$endir"/"$runam"/error_dedup.log
	shopt -u nullglob
}

#make inference on genetic sex and aneuploidies, estimating contamination
function sexdeter_contam {
	basename=$( cat $listf | cut -d, -f1 )
	echo "Inferring GENETIC SEX, ANEUPLOIDIES and levels of CONTAMINATION (if set) assuming reference is human"
	if [ $erindir == TRUE ]; then
		endir=$( echo $indir )
	fi
	shopt -s nullglob
	for i in $basename; do
		echo "Running $i"
		if [ ! -d "$endir"/"$runam"/"$i"/RUNLOGS ]; then
			mkdir  "$endir"/"$runam"/"$i"/RUNLOGS
		fi
		if [ -d "$endir"/"$runam"/"$i"/ALNFILES/DEDUP ] && [ `ls -1 "$endir"/"$runam"/"$i"/ALNFILES/DEDUP/*bam 2>/dev/null | wc -l` != 0 ]; then
			temppath1=$( echo "$endir"/"$runam"/"$i"/ALNFILES/DEDUP )
		elif [ -d "$endir"/"$runam"/"$i"/ALNFILES ] && [ `ls -1 "$endir"/"$runam"/"$i"/ALNFILES/*bam 2>/dev/null | wc -l` != 0 ]; then
			temppath1=$( echo "$endir"/"$runam"/"$i"/ALNFILES )
			echo "BAM file(s) may not went under PCR duplication removal, assuming they did, but if not, results will probably be highly distorted ...just sayin"
		elif [ `ls -1 "$indir"/*bam 2>/dev/null | wc -l` != 0 ]; then
			temppath1="$indir"
		else
			echo "BAM file does not exist in given directory, exiting 1"
			exit 1
		fi
		mkdir "$endir"/"$runam"/"$i"/tmp
		temppath=$( echo "$endir"/"$runam"/"$i"/tmp )
		countbam=$( ls -1 "$temppath1"/"$i"*bam 2>/dev/null | wc -l )
		if [ $countbam -gt 1 ]; then
			samtools merge "$temppath"/"$i".merged.bam "$temppath1"/"$i"*bam
			samtools sort -o "$temppath"/"$i".merged.sorted.bam "$temppath"/"$i".merged.bam
			rm -rf "$temppath"/"$i".merged.bam
			samtools index "$temppath"/"$i".merged.sorted.bam
			tmpset="T"
		elif [ $countbam -eq 1 ]; then
			temppath=$( echo "$temppath1" )
		fi
		countbai=$( ls -1 "$temppath"/"$i"*bai 2>/dev/null | wc -l )
		if [ $countbai -lt 1 ]; then
			echo "BAM index does not exists for $i, creating one (may take some time)"
			samtools index "$temppath"/"$i"*bam
		fi
		if [ $workf == FQ ] && [ $overl != 0 ] && [ ! -z $overl ]; then
			if [ $mapqt -gt 0 ]; then
				samtools view -q $mapqt -b "$temppath"/"$i"*bam > "$temppath"/"$i"_mq_bam
				genomeCoverageBed -ibam "$temppath"/"$i"_mq_bam > "$endir"/"$runam"/"$i"/RUNLOGS/"$i".coverage_hist.txt
				samtools stats "$temppath"/"$i"_mq_bam | grep ^RL | cut -f 2- > "$endir"/"$runam"/"$i"/RUNLOGS/"$i".readlen_distr.txt
			else
				genomeCoverageBed -ibam "$temppath"/"$i"*bam > "$endir"/"$runam"/"$i"/RUNLOGS/"$i".coverage_hist.txt
				samtools stats "$temppath"/"$i"*bam | grep ^RL | cut -f 2- > "$endir"/"$runam"/"$i"/RUNLOGS/"$i".readlen_distr.txt
			fi
		elif [ $workf == FQ ] && [ $overl == 0 ]; then
			if [ $mapqt -gt 0 ]; then
				samtools view -q $mapqt -b "$temppath"/"$i"*bam > "$temppath"/"$i"_mq_bam
				genomeCoverageBed -pc -ibam "$temppath"/"$i"_mq_bam > "$endir"/"$runam"/"$i"/RUNLOGS/"$i".coverage_hist.txt
				echo "For PE alignments, read length distribution is calculated from the first fragments"
				samtools stats "$temppath"/"$i"_mq_bam | grep ^FRL | cut -f 2- > "$endir"/"$runam"/"$i"/RUNLOGS/"$i".readlen_distr.txt
			else
				genomeCoverageBed -pc -ibam "$temppath"/"$i"*bam > "$endir"/"$runam"/"$i"/RUNLOGS/"$i".coverage_hist.txt
				samtools stats "$temppath"/"$i"*bam | grep ^RL | cut -f 2- > "$endir"/"$runam"/"$i"/RUNLOGS/"$i".readlen_distr.txt
			fi
		else
			echo "SE alignment is considered for coverage and read length calculations"
			if [ $mapqt -gt 0 ]; then
				samtools view -q $mapqt -b "$temppath"/"$i"*bam > "$temppath"/"$i"_mq_bam
				genomeCoverageBed -ibam "$temppath"/"$i"_mq_bam > "$endir"/"$runam"/"$i"/RUNLOGS/"$i".coverage_hist.txt
				samtools stats "$temppath"/"$i"_mq_bam | grep ^RL | cut -f 2- > "$endir"/"$runam"/"$i"/RUNLOGS/"$i".readlen_distr.txt
			else
				genomeCoverageBed -ibam "$temppath"/"$i"*bam > "$endir"/"$runam"/"$i"/RUNLOGS/"$i".coverage_hist.txt
				samtools stats "$temppath"/"$i"*bam | grep ^RL | cut -f 2- > "$endir"/"$runam"/"$i"/RUNLOGS/"$i".readlen_distr.txt
			fi
		fi
		xarray=("chrX" "chrx" "x" "X" "23" "chr23")
		yarray=("chrY" "chry" "y" "Y" "24" "chr24")
		for j in ${xarray[@]}; do
			check=$( samtools view -H "$temppath"/"$i"*bam | grep '^@SQ' | cut -f2 | grep $j | wc -l )
			if [ $check == 1 ]; then
				xref=$( echo "$j" )
				break
			elif [ $check -gt 1 ]; then
				echo "The reference genome used for BAM making has more than one of the following sequence names: chrX, chrx, x, X, 23, chr23. Please make sure that your reference does not contain duplicates."
				exit 1
		done
		if [ $check == 0 ]; then
			skipping=T
		fi
		check=""
		for j in ${yarray[@]}; do
			check=$( samtools view -H "$temppath"/"$i"*bam | grep '^@SQ' | cut -f2 | grep $j | wc -l )
			if [ $check == 1 ]; then
				yref=$( echo "$j" )
				break
			elif [ $check -gt 1 ]; then
				echo "The reference genome used for BAM making has more than one of the following sequence names: chrY, chry, y, Y, 24, chr24. Please make sure that your reference does not contain duplicates."
				exit 1
		done
		if [ $check == 0 ]; then
			skipping=T
		fi
		check=""
		mitoarray=("MT" "M" "MITO" "chrM" "chrMT" "chrMITO" "mitochondrial" "NC_012920.1" "rCRS" "rcrs" "CRS" "crs" "RSRS" "rsrs" "mt" "m" "chrmt" "chrm")
		for j in ${mitoarray[@]}; do
			check=$( samtools view -H "$temppath"/"$i"*bam | grep '^@SQ' | cut -f2 | grep $j | wc -l )
			if [ $check == 1 ]; then
				mitorefname=$( echo "$j" )
				break
			elif [ $check -gt 1 ]; then
				echo "The reference genome used for BAM making has more than one of the following sequence names: m, M, chrm, chrM, mt, MT, chrmt, chrMT, mito, MITO, chrmito, chrMITO, rcrs, rCRS, rsrs, RSRS, mitochondrial, NC_012920.1. Please make sure that your reference does not contain duplicates."
				exit 1
			fi
		done
		if [ $check == 0 ]; then
			echo "The reference does not contain mtDNA sequence name, aneuploidy will be performed, but check your reference for further steps"
			mitorefname="NA"
		fi
		if [ $skipping == T ]; then
			echo "The reference genome does not contain any of the following chromosome names: chrX, chrx, x, X, 23, chr23, chrY, chry, y, Y, 24, chr24, assuming not human reference, skipping ploidy test"
		else
			xxx=$( samtools view -F 4 "$temppath"/"$i"*bam $xref | wc -l )
			yyy=$( samtools view -F 4 "$temppath"/"$i"*bam $yref | wc -l )
			echo $xxx  > "$endir"/"$runam"/"$i"/RUNLOGS/"$i".chrx_readNo.txt
			echo $yyy  > "$endir"/"$runam"/"$i"/RUNLOGS/"$i".chry_readNo.txt
			Rscript --vanilla ./src/gensex.r "$endir"/"$runam"/"$i"/RUNLOGS/"$i".coverage_hist.txt "$endir"/"$runam"/"$i"/RUNLOGS/ $xxx $yyy $i $xref $yref ./references/basedata.csv $mitorefname
		fi
		sex=$( cat "$endir"/"$runam"/"$i"/RUNLOGS/*sex.txt )
		if [ $sex == "MALE" ]; then
			if ! command -v Yleaf.py &> /dev/null; then
				echo "Yleaf is not installed or not in your path, skipping Y chromosome haplogroup detection, skipping"
			elif []; then
			else
				echo "For chrY haplogroup inference PAPline considers hg19 as reference genome"
				if [ $cotrv == TRUE ] || [ $cotrv == T ] || [ $cotrv == True ] || [ $cotrv == true ] || [ $cotrv == t ]; then
					posfiley=$( echo ./references/WGS_hg19_TRV.txt )
				else
					posfiley=$( echo ~/Position_files/WGS_hg19.txt )
				fi
				if [ $ctrim != 0 ]; then
					samtools view -b "$temppath"/"$i"*bam $yref > "$endir"/"$runam"/"$i"/RUNLOGS/"$i".chrY.bam
					bam trimBam "$endir"/"$runam"/"$i"/RUNLOGS/"$i".chrY.bam "$endir"/"$runam"/"$i"/RUNLOGS/"$i".chrY_trim.bam $ctrim
					Yleaf.py -bam "$endir"/"$runam"/"$i"/RUNLOGS/"$i".chrY_trim.bam -pos $posfiley -out "$endir"/"$runam"/"$i"/RUNLOGS/"$i"_chrY_hg -r 1 -q $mapqt -b 90 -t $thred
					rm -rf "$endir"/"$runam"/"$i"/RUNLOGS/"$i".chrY.bam
					rm -rf "$endir"/"$runam"/"$i"/RUNLOGS/"$i".chrY_trim.bam
				else
					Yleaf.py -bam "$temppath"/"$i"*bam -pos $posfiley -out "$endir"/"$runam"/"$i"/RUNLOGS/"$i"_chrY_hg -r 1 -q $mapqt -b 90 -t $thred
				fi
			fi
		fi
		if [ $cntmx != "None" ]; then
			if [ $cntmx != "Male_chrX" ] && [ $cntmx != "mtDNA" ] && [ $cntmx != "Joint" ] && [ $cntmx != "Exclusive" ] && [ $cntmx != "Inclusive" ]; then
				echo "Contamination estimation method is unrecogniseable given by user, please take care of letter cases and typos, now calculating with Inclusive method"
				cntmx="Inclusive"
			fi
			function contamix {
				if [ ! -f "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".mito.bam ]; then
					echo "Creating mtDNA BAM file"
					samtools view -b "$temppath"/"$i"*bam $mitorefname > "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".mito.bam
					samtools view -H "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".mito.bam | grep HD > hd.txt
					samtools view -H "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".mito.bam | grep $mitorefname > chrm.txt
					samtools view -H "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".mito.bam | grep PG > pg.txt
					cat hd.txt chrm.txt pg.txt > "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".mito.sam
					rm hd.txt
					rm chrm.txt
					rm pg.txt
					samtools view -F 4 "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".mito.bam >> "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".mito.sam
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".mito.bam
					samtools view -hbS "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".mito.sam > "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".mito.bam
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".mito.sam
					samtools index "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".mito.bam
				fi
				echo "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".mito.bam > path.txt
				if [ $ctrim -gt 0 ]; then
					bam trimBam "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".mito.bam "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".trimmed_"$ctrim".mito.bam $ctrim
					samtools index "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".trimmed_"$ctrim".mito.bam
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".mito.bam
					echo "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".trimmed_"$ctrim".mito.bam > path.txt
				else
					samtools index "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".mito.bam
					echo "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".mito.bam > path.txt
				fi
				if [ $mtfas != Intelligent ] && [ $mtfas != Majority_rule ] && [ $mtfas != Random ]; then
					mtfas2=Intelligent
					echo "FASTA making was not set or not set correctly, for contamination estimation it is set to Intelligent"
				else
					mtfas2=$( echo $mtfas )
				fi
				if [ $mtfas2 == Intelligent ]; then
					angsd -doCounts 1 -minQ $baseq -dumpCounts 3 -bam path.txt -out "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i"
					rm path.txt
					gunzip -f "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".counts.gz
					gunzip -f "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".pos.gz
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".counts.gz
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".pos.gz
					Rscript --vanilla ./src/fastamaker_intelligent.r "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/ "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".pos "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".counts "$i" "$mtcov"
				elif [ $mtfas2 == Majority_rule ]; then
					angsd -doCounts 1 -minQ $baseq -dumpCounts 3 -bam path.txt -out "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i"
					rm path.txt
					gunzip -f "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".counts.gz
					gunzip -f "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".pos.gz
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".counts.gz
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".pos.gz
					Rscript --vanilla ./src/fastamaker_majority.r "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/ "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".pos "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".counts "$i" "$mtcov"
				elif [ $mtfas2 == Random ]; then
					angsd -doCounts 1 -minQ $baseq -dumpCounts 3 -bam path.txt -out "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i"
					rm path.txt
					gunzip -f "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".counts.gz
					gunzip -f "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".pos.gz
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".counts.gz
					rm -rf "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".pos.gz
					Rscript --vanilla ./src/fastamaker_random.r "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/ "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".pos "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".counts "$i" "$mtcov"
				fi
				echo "Align $i consensus to 311 reference"
				cat "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".mito.fasta ./references/311.fasta > ./references/311.tmp.fasta
				mafft --auto --thread $thred ./references/311.tmp.fasta > ./references/311.tmp.aln.fasta
				rm -rf ./references/311.tmp.fasta
				samtools view -h -o "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".nsc.sam "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".mito.bam
				awk 'BEGIN {OFS="\t"} {split($6,C,/[0-9]*/); split($6,L,/[SMDIN]/); if (C[2]=="S") {$10=substr($10,L[1]+1); $11=substr($11,L[1]+1)}; if (C[length(C)]=="S") {L1=length($10)-L[length(L)-1]; $10=substr($10,1,L1); $11=substr($11,1,L1); }; gsub(/[0-9]*S/,"",$6); print}' "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".nsc.sam > "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".sc.sam
				rm -rf "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".nsc.sam
				samtools fastq "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".sc.sam > "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".fq
				bwa index "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".mito.fasta
				echo "Remapping $i to consensus"
				bwa mem -t $thred "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".mito.fasta "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".fq > "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".remap.sam
				rm -rf "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".fq
				samtools view -q $mapqt -bT "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".mito.fasta "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".remap.sam > "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".remap.bam
				rm -rf "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".remap.sam
				samtools sort -o "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".sort.bam "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".remap.bam
				rm -rf "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".remap.bam
				samtools index "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".sort.bam
				echo "Estimating contamination"
				if [ $ctrim -lt 1 ]; then
					ctrim2=2
				else
					ctrim2=$( echo $ctrim )
				fi
				estimate.R --samFn "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".sort.bam --malnFn ./references/311.tmp.aln.fasta --trimBases $ctrim2 --figure "$endir"/"$runam"/"$i"/RUNLOGS/"$i".out.pdf > "$endir"/"$runam"/"$i"/RUNLOGS/"$i".contamix.out
			}
			function chromoxm {
				if [ $ctrim -lt 1 ]; then
					ctrim2=2
				else
					ctrim2=$( echo $ctrim )
				fi
				angsd -i "$temppath"/"$i"*bam -r "$xref":5000000-154900000 -doCounts 1 -iCounts 1 -minMapQ $mapqt -minQ $baseq -trim $ctrim2 -out "$endir"/"$runam"/"$i"/RUNLOGS/"$i".chrXcontam
				Rscript ~/R/contamination.R mapFile="~/RES/chrX.unique.gz" hapFile="~/RES/HapMapChrX.gz" countFile="$endir"/"$runam"/"$i"/RUNLOGS/"$i".chrXcontam.icnts.gz mc.cores="$thred" > "$endir"/"$runam"/"$i"/RUNLOGS/"$i".chrXcontam.est
			}
			xcov=$( cat "$endir"/"$runam"/"$i"/RUNLOGS/*chrx_cov.txt )
			mcov=$( cat "$endir"/"$runam"/"$i"/RUNLOGS/*mtdna_cov.txt )
			sex=$( cat "$endir"/"$runam"/"$i"/RUNLOGS/*sex.txt )
			if [ $cntmx == "Male_chrX" ]; then
				if [ $xcov -gt 0.5 ] && [ $sex == "MALE" ]; then
					if ! command -v ~/R/contamination.R &> /dev/null; then
						echo "~/R/contamination.R is not installed to path, skipping"
					else
						echo "Running Male_chrX contamination estimation for $i"
						if $(chromoxm); then
							echo "Contamination estimation is finished"
						else
							echo "Something went wrong, exciting"
						fi
					fi
				else
					echo "Male_chrX contamination estimation for $i is skipped"
					echo "NA" > "$endir"/"$runam"/"$i"/RUNLOGS/"$i".chrXcontam.est
				fi
				echo "NA" > "$endir"/"$runam"/"$i"/RUNLOGS/"$i".contamix.out
			elif [ $cntmx == "mtDNA" ]; then
				if ! command -v estimate.R &> /dev/null; then
					echo "ContamMix is not installed or is not in your path, skipping"
				else
					echo "Running mtDNA contamination estimation for $i"
					if $(contamix); then
						echo "Contamination estimation is finished"
						echo "NA" > "$endir"/"$runam"/"$i"/RUNLOGS/"$i".chrXcontam.est
					else
						echo "Something went wrong, exciting"
					fi
				fi
			elif [ $cntmx == "Joint" ]; then
				if ! command -v estimate.R &> /dev/null || ! command -v ~/R/contamination.R &> /dev/null; then
					echo "ContamMix and/or ANGSD contamination estimator is not installed or is not in your path, skipping"
				else
					if [ $xcov -gt 0.5 ]; then
						echo "Running Male_chrX contamination estimation for $i"
						if $(chromoxm); then
							echo "Contamination estimation is finished"
						else
							echo "Something went wrong, exciting"
						fi
					else
						echo "Male_chrX contamination estimation for $i is skipped for low coverage"
						echo "NA" > "$endir"/"$runam"/"$i"/RUNLOGS/"$i".chrXcontam.est
					fi
					echo "Running mtDNA contamination estimation for $i"
					if $(contamix); then
						echo "Contamination estimation is finished"
					else
						echo "Something went wrong, exciting"
					fi
				fi
			elif [ $cntmx == "Exclusive" ]; then
				if ! command -v estimate.R &> /dev/null || ! command -v ~/R/contamination.R &> /dev/null; then
					echo "ContamMix and/or ANGSD contamination estimator is not installed or is not in your path, skipping"
				else
					if [ $xcov -gt 0.5 ] && [ $sex == "MALE" ]; then
						echo "Running Male_chrX contamination estimation for $i"
						if $(chromoxm); then
							echo "Contamination estimation is finished"
						else
							echo "Something went wrong, exciting"
						fi
					elif [ $xcov -lt 0.5 ] || [ $sex != "MALE" ]; then
						echo "NA" > "$endir"/"$runam"/"$i"/RUNLOGS/"$i".chrXcontam.est
						echo "Running mtDNA contamination estimation for $i"
						if $(contamix); then
							echo "Contamination estimation is finished"
						else
							echo "Something went wrong, exciting"
						fi
					fi
				fi
			elif [ $cntmx == "Inclusive" ]; then
				if ! command -v estimate.R &> /dev/null || ! command -v ~/R/contamination.R &> /dev/null; then
					echo "ContamMix and/or ANGSD contamination estimator is not installed or is not in your path, skipping"
				else
					if [ $xcov -gt 0.5 ] && [ $sex == "MALE" ]; then
						echo "Running Male_chrX contamination estimation for $i"
						if $(chromoxm); then
							echo "Contamination estimation is finished"
						else
							echo "Something went wrong, exciting"
						fi
					else
						echo "NA" > "$endir"/"$runam"/"$i"/RUNLOGS/"$i".chrXcontam.est
					fi
					echo "Running mtDNA contamination estimation for $i"
					if $(contamix); then
						echo "Contamination estimation is finished"
					else
						echo "Something went wrong, exciting"
					fi
				fi
			fi
		fi
		rm -rf "$endir"/"$runam"/"$i"/tmp
	done 2> "$endir"/"$runam"/error_sexdet_contam.log
	shopt -u nullglob
}

#inferring damage patterns FIN
function mapdamage {
	basename=$( cat $listf | cut -d, -f1 )
	echo "Measuring DNA DAMAGE"
	if [ $erindir == TRUE ]; then
		endir=$( echo $indir )
	fi
	shopt -s nullglob
	for i in $basename; do
		if [ ! -d "$endir"/"$runam"/"$i"/RUNLOGS ]; then
			mkdir  "$endir"/"$runam"/"$i"/RUNLOGS
		fi
		if [ -d "$endir"/"$runam"/"$i"/ALNFILES/DEDUP ] && [ `ls -1 "$endir"/"$runam"/"$i"/ALNFILES/DEDUP/*bam 2>/dev/null | wc -l` != 0 ]; then
			temppath1=$( echo "$endir"/"$runam"/"$i"/ALNFILES/DEDUP )
		elif [ -d "$endir"/"$runam"/"$i"/ALNFILES ] && [ `ls -1 "$endir"/"$runam"/"$i"/ALNFILES/*bam 2>/dev/null | wc -l` != 0 ]; then
			temppath1=$( echo "$endir"/"$runam"/"$i"/ALNFILES )
		elif [ `ls -1 "$indir"/*bam 2>/dev/null | wc -l` != 0 ]; then
			temppath1="$indir"
		else
			echo "BAM file does not exist in given directory, exiting 1"
			exit 1
		fi
		mkdir "$endir"/"$runam"/"$i"/tmp
		temppath=$( echo "$endir"/"$runam"/"$i"/tmp )
		countbam=$( ls -1 "$temppath1"/"$i"*bam 2>/dev/null | wc -l )
		if [ $countbam -gt 1 ]; then
			samtools merge "$temppath"/"$i".merged.bam "$temppath1"/"$i"*bam
			samtools sort -o "$temppath"/"$i".merged.sorted.bam "$temppath"/"$i".merged.bam
			rm -rf "$temppath"/"$i".merged.bam
			tmpset="T"
		elif [ $countbam -eq 1 ]; then
			temppath=$( echo "$temppath1" )
		fi
		mkdir "$endir"/"$runam"/"$i"/RUNLOGS/MAPDAM
		echo "DNA damage pattern inference for $i"
		mapDamage -n 100000 -v -t "$i" -i "$temppath"/"$i"*bam -r $refer
		cp ./results_"$i"*/* "$endir"/"$runam"/"$i"/RUNLOGS/MAPDAM
		rm -rf ./results_"$i"*/
		echo $( awk '{for(i=1;i<=NF;i++) if ($i=="1") print $(i+1)}' "$endir"/"$i"/5pCtoT_freq.txt ) > "$endir"/"$runam"/"$i"/RUNLOGS/"$i".mapdam.txt #mapdamage report
		if [ $krbam == FALSE ] || [ $krbam == false ] || [ $krbam == F ] || [ $krbam == f ]; then
			if [ ! -d "$endir"/"$runam"/"$i"/ALNFILES/RAW ]; then
				mkdir "$endir"/"$runam"/"$i"/ALNFILES/RAW
			fi
			if [ $tmpset == "T" ]; then
				mv "$endir"/"$runam"/"$i"/tmp/"$i"*sorted.bam "$endir"/"$runam"/"$i"/ALNFILES/RAW/"$i".sorted.bam
			else
				cp "$endir"/"$runam"/"$i"/ALNFILES/"$i"*bam "$endir"/"$runam"/"$i"/ALNFILES/RAW/
				cp "$endir"/"$runam"/"$i"/ALNFILES/"$i"*bai "$endir"/"$runam"/"$i"/ALNFILES/RAW/
				rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i"*bam
				rm -rf "$endir"/"$runam"/"$i"/ALNFILES/"$i"*bai
			fi
		fi
		rm -rf "$endir"/"$runam"/"$i"/tmp
	done 2> "$endir"/"$runam"/error_mapdam.log
	shopt -u nullglob
}

#make mtDNA fasta and run haplogrep FIN
function fastamaker {
	basename=$( cat $listf | cut -d, -f1 )
	echo "Making mtDNA FASTA files"
	if [ $erindir == TRUE ]; then
		endir=$( echo $indir )
	fi
	shopt -s nullglob
	for i in $basename; do
		if [ ! -d "$endir"/"$runam"/"$i"/RUNLOGS ]; then
			mkdir  "$endir"/"$runam"/"$i"/RUNLOGS
		fi
		if [ -d "$endir"/"$runam"/"$i"/ALNFILES/DEDUP ] && [ `ls -1 "$endir"/"$runam"/"$i"/ALNFILES/DEDUP/*bam 2>/dev/null | wc -l` != 0 ]; then
			temppath1=$( echo "$endir"/"$runam"/"$i"/ALNFILES/DEDUP )
		elif [ -d "$endir"/"$runam"/"$i"/ALNFILES ] && [ `ls -1 "$endir"/"$runam"/"$i"/ALNFILES/*bam 2>/dev/null | wc -l` != 0 ]; then
			temppath1=$( echo "$endir"/"$runam"/"$i"/ALNFILES )
		elif [ `ls -1 "$indir"/*bam 2>/dev/null | wc -l` != 0 ]; then
			temppath1="$indir"
		else
			echo "BAM file does not exist in given directory, exiting 1"
			exit 1
		fi
		mkdir "$endir"/"$runam"/"$i"/tmp
		temppath=$( echo "$endir"/"$runam"/"$i"/tmp )
		countbam=$( ls -1 "$temppath1"/"$i"*bam 2>/dev/null | wc -l )
		if [ $countbam -gt 1 ]; then
			samtools merge "$temppath"/"$i".merged.bam "$temppath1"/"$i"*bam
			samtools sort -o "$temppath"/"$i".merged.sorted.bam "$temppath"/"$i".merged.bam
			rm -rf "$temppath"/"$i".merged.bam
			tmpset="T"
		elif [ $countbam -eq 1 ]; then
			temppath=$( echo "$temppath1" )
		fi
		mitoarray=("MT" "M" "MITO" "chrM" "chrMT" "chrMITO" "mitochondrial" "NC_012920.1" "rCRS" "rcrs" "CRS" "crs" "RSRS" "rsrs" "mt" "m" "chrmt" "chrm")
		for j in ${mitoarray[@]}; do
			check=$( samtools view -H "$temppath"/"$i"*bam | grep '^@SQ' | cut -f2 | grep $j | wc -l )
			if [ $check == 1 ]; then
				mitorefname=$( echo "$j" )
				break
			elif [ $check -gt 1 ]; then
				echo "The reference genome used for BAM making has more than one of the following sequence names: m, M, chrm, chrM, mt, MT, chrmt, chrMT, mito, MITO, chrmito, chrMITO, rcrs, rCRS, rsrs, RSRS, mitochondrial, NC_012920.1. Please make sure that your reference does not contain duplicates."
				exit 1
			fi
		done
		if [ $check == 0 ]; then
			echo "The reference genome does not contain any of the following chromosome names: m, M, chrm, chrM, mt, MT, chrmt, chrMT, mito, MITO, chrmito, chrMITO, rcrs, rCRS, rsrs, RSRS, mitochondrial, NC_012920.1. Check your reference genome!"
			exit 1
		fi
		rm -rf "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL
		mkdir "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL
		echo "Making MT FASTA for $i"
		countbai=$( ls -1 "$temppath"/"$i"*bai 2>/dev/null | wc -l )
		if [ $countbai -lt 1 ]; then
			echo "BAM index does not exists for $i, creating one (may take some time)"
			samtools index "$temppath"/"$i"*bam
		fi 
		echo "Creating mtDNA BAM file"
		samtools view -b "$temppath"/"$i"*bam $mitorefname > "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".mito.bam
		samtools view -H "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".mito.bam | grep HD > hd.txt
		samtools view -H "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".mito.bam | grep $mitorefname > chrm.txt
		samtools view -H "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".mito.bam | grep PG > pg.txt
		cat hd.txt chrm.txt pg.txt > "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".mito.sam
		rm hd.txt
		rm chrm.txt
		rm pg.txt
		samtools view -F 4 "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".mito.bam >> "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".mito.sam
		rm -rf "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".mito.bam
		samtools view -hbS "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".mito.sam > "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".mito.bam
		rm -rf "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".mito.sam
		samtools index "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".mito.bam
		echo "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".mito.bam > path.txt
		if [ $ctrim -gt 0 ]; then
			bam trimBam "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".mito.bam "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".trimmed_"$ctrim".mito.bam $ctrim
			samtools index "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".trimmed_"$ctrim".mito.bam
			rm -rf "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".mito.bam
			echo "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".trimmed_"$ctrim".mito.bam > path.txt
		else
			samtools index "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".mito.bam
			echo "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".mito.bam > path.txt
		fi
		if [ $mtfas != Intelligent ] && [ $mtfas != Majority_rule ] && [ $mtfas != Random ]; then
			mtfas=Intelligent
			echo "FASTA making was not set or not set correctly, for contamination estimation it is set to Intelligent"
		fi
		if [ $mtfas == Intelligent ]; then
			angsd -doCounts 1 -minQ $baseq -dumpCounts 3 -bam path.txt -out "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i"
			rm path.txt
			gunzip -f "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".counts.gz
			gunzip -f "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".pos.gz
			rm -rf "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".counts.gz
			rm -rf "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".pos.gz
			Rscript --vanilla ./src/fastamaker_intelligent.r "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/ "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".pos "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".counts "$i" "$mtcov"
		elif [ $mtfas == Majority_rule ]; then
			angsd -doCounts 1 -minQ $baseq -dumpCounts 3 -bam path.txt -out "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i"
			rm path.txt
			gunzip -f "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".counts.gz
			gunzip -f "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".pos.gz
			rm -rf "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".counts.gz
			rm -rf "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".pos.gz
			Rscript --vanilla ./src/fastamaker_majority.r "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/ "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".pos "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".counts "$i" "$mtcov"
		elif [ $mtfas == Random ]; then
			angsd -doCounts 1 -minQ $baseq -dumpCounts 3 -bam path.txt -out "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i"
			rm path.txt
			gunzip -f "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".counts.gz
			gunzip -f "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".pos.gz
			rm -rf "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".counts.gz
			rm -rf "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".pos.gz
			Rscript --vanilla ./src/fastamaker_random.r "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/ "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".pos "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".counts "$i" "$mtcov"
		fi
		if ! command -v haplogrep &> /dev/null; then
			echo "Haplogrep is not installed despite $mtfas method is set, check dependencies"
			exit 1
		fi
		rm path.txt
		echo "Classifying mtDNA"
		haplogrep classify --in "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".mito.fasta --format fasta --output "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".haplogrep_class_report
		echo $( echo $( sed '2q;d' "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".haplogrep_class_report | awk '{print $2}' ) | tr -d '"' ) > "$endir"/"$runam"/"$i"/ALNFILES/RUNLOGS/"$i".hg.txt
		echo $( echo $( sed '2q;d' "$endir"/"$runam"/"$i"/ALNFILES/MITOCHONDRIAL/"$i".haplogrep_class_report | awk '{print $4}' ) | tr -d '"' ) > "$endir"/"$runam"/"$i"/ALNFILES/RUNLOGS/"$i".hgq.txt
	done 2> "$endir"/"$runam"/error_mitoref.log
	shopt -u nullglob
}


function pigment_disease {
	basename=$( cat $listf | cut -d, -f1 )
	echo "Analysing variants, assuming HG19 as reference genome"
	if [ $erindir == TRUE ]; then
		endir=$( echo $indir )
	fi
	shopt -s nullglob
	for i in $basename; do
		if [ ! -d "$endir"/"$runam"/"$i"/RUNLOGS ]; then
			mkdir  "$endir"/"$runam"/"$i"/RUNLOGS
		fi
		if [ -d "$endir"/"$runam"/"$i"/ALNFILES/DEDUP ] && [ `ls -1 "$endir"/"$runam"/"$i"/ALNFILES/DEDUP/*bam 2>/dev/null | wc -l` != 0 ]; then
			temppath1=$( echo "$endir"/"$runam"/"$i"/ALNFILES/DEDUP )
		elif [ -d "$endir"/"$runam"/"$i"/ALNFILES ] && [ `ls -1 "$endir"/"$runam"/"$i"/ALNFILES/*bam 2>/dev/null | wc -l` != 0 ]; then
			temppath1=$( echo "$endir"/"$runam"/"$i"/ALNFILES )
		elif [ `ls -1 "$indir"/*bam 2>/dev/null | wc -l` != 0 ]; then
			temppath1="$indir"
		else
			echo "BAM file does not exist in given directory, exiting 1"
			exit 1
		fi
		mkdir "$endir"/"$runam"/"$i"/tmp
		temppath=$( echo "$endir"/"$runam"/"$i"/tmp )
		countbam=$( ls -1 "$temppath1"/"$i"*bam 2>/dev/null | wc -l )
		if [ $countbam -gt 1 ]; then
			samtools merge "$temppath"/"$i".merged.bam "$temppath1"/"$i"*bam
			samtools sort -o "$temppath"/"$i".merged.sorted.bam "$temppath"/"$i".merged.bam
			rm -rf "$temppath"/"$i".merged.bam
			tmpset="T"
		elif [ $countbam -eq 1 ]; then
			temppath=$( echo "$temppath1" )
		fi
		reflist=$( samtools idxstats "$temppath"/"$i"*bam | cut -f 1 | uniq | head -25 )
		echo $reflist > "$endir"/"$runam"/"$i"/RUNLOGS/refnames.txt
		checked=""
		for i in {1..22}; do 
			if [[ $reflist != *$i* ]]; then
				echo "Your reference does not contain the chromosome set of a human genome, variant calling will be skipped"
				checked=F
			fi
		done
		if [[ $reflist != *25* ]]; then
			echo "Your reference does not contain the chromosome set of a human genome, variant calling will be skipped"
			checked=F
		fi
		if [ $checked != F ]; then
			Rscript --vanilla ./src/clinical1.r ./references/full_panel.csv $cotrv ./tmp/ "$endir"/"$runam"/"$i"/RUNLOGS/refnames.txt $refbg 
			angsd -gl 2 -nThreads $thred -doGlf 2 -minQ $bqthr -minMapQ $mapqt -i "$temppath"/"$i"*bam -rf ./tmp/snps.bed -doMajorMinor 1 -doCounts 1 -dumpCounts 3 -trim $ctrim -out "$endir"/"$runam"/"$i"/RUNLOGS/"$i"_clin
			Rscript --vanilla ./src/clinical2.r "$endir"/"$runam"/"$i"/RUNLOGS/"$i"_clin.counts.gz "$endir"/"$runam"/"$i"/RUNLOGS/"$i"_clin.pos.gz $disis $pigme ./references/full_panel.csv $refbg "$endir"/"$runam"/"$i"/RUNLOGS/refnames.txt "$endir"/"$runam"/"$i"/RUNLOGS/"$i"_clin.beagle.gz "$endir"/"$runam"/"$i"/RUNLOGS/ "$i"
		fi
		rm -rf "$endir"/"$runam"/"$i"/tmp
	done 2> "$endir"/"$runam"/error_mapdam.log
	shopt -u nullglob
}

function genocall {
	basename=$( cat $listf | cut -d, -f1 )
	echo "Analysing variants, assuming HG19 as reference genome"
	if [ $erindir == TRUE ]; then
		endir=$( echo $indir )
	fi
	if [ -f "$endir"/"$runam"/bamlist.txt ]; then
		rm -rf "$endir"/"$runam"/bamlist.txt
	fi
	if [ -z $eigen ]; then
		echo "SNP file was not provided for genotype call, skipping"
	elif [ $psudo != "Pseudohaploid_Random" ] && [ $psudo != "Pseudohaploid_Majority" ] && [ $psudo != "Diploid_Random" ] && [ $psudo != "Diploid_GL" ]; then
		echo "psudo option was incorrectly given, please check for typos and letter cases, now skipping"
	else
		shopt -s nullglob
		for i in $basename; do
			if [ ! -d "$endir"/"$runam"/"$i"/RUNLOGS ]; then
				mkdir  "$endir"/"$runam"/"$i"/RUNLOGS
			fi
			if [ -d "$endir"/"$runam"/"$i"/ALNFILES/DEDUP ] && [ `ls -1 "$endir"/"$runam"/"$i"/ALNFILES/DEDUP/*bam 2>/dev/null | wc -l` != 0 ]; then
				temppath1=$( echo "$endir"/"$runam"/"$i"/ALNFILES/DEDUP )
			elif [ -d "$endir"/"$runam"/"$i"/ALNFILES ] && [ `ls -1 "$endir"/"$runam"/"$i"/ALNFILES/*bam 2>/dev/null | wc -l` != 0 ]; then
				temppath1=$( echo "$endir"/"$runam"/"$i"/ALNFILES )
			elif [ `ls -1 "$indir"/*bam 2>/dev/null | wc -l` != 0 ]; then
				temppath1="$indir"
			else
				echo "BAM file does not exist in given directory, exiting 1"
				exit 1
			fi
			mkdir "$endir"/"$runam"/"$i"/tmp
			temppath=$( echo "$endir"/"$runam"/"$i"/tmp )
			countbam=$( ls -1 "$temppath1"/"$i"*bam 2>/dev/null | wc -l )
			if [ $countbam -gt 1 ]; then
				samtools merge "$temppath"/"$i".merged.bam "$temppath1"/"$i"*bam
				samtools sort -o "$temppath"/"$i".merged.sorted.bam "$temppath"/"$i".merged.bam
				rm -rf "$temppath"/"$i".merged.bam
				tmpset="T"
			elif [ $countbam -eq 1 ]; then
				temppath=$( echo "$temppath1" )
			fi
			if [ $ctrim -gt 0 ]; then
				if [ ! -d "$temppath"/ENDTRIM ] || [ `ls -1 "$temppath"/ENDTRIM/*bam 2>/dev/null | wc -l` == 0 ]; then
					echo "Trimming read ends for $i"
					mkdir $temppath/ENDTRIM
					bam trimBam "$temppath1"/"$i"*bam "$temppath"/ENDTRIM/"$i".trim"$ctrim".bam $ctrim
					samtools index "$temppath"/ENDTRIM/"$i".trim"$ctrim".bam
					temppath=$( echo "$temppath"/ENDTRIM )
				else if [ `ls -1 "$temppath"/ENDTRIM/*bam 2>/dev/null | wc -l` -gt 1 ]; then
					echo "ENDTRIM library contains more than one BAM file for $i, quitting"
					exit 1
				else if [ `ls -1 "$temppath"/ENDTRIM/*bam 2>/dev/null | wc -l` == 1 ] && [ `ls -1 "$temppath"/ENDTRIM/*bai 2>/dev/null | wc -l` ==0 ]; then
					"BAM indexing $i"
					samtools index "$temppath"/ENDTRIM/"$i"*bam
					temppath=$( echo "$temppath"/ENDTRIM )
				fi
			fi
			if [ `ls -1 "$temppath"/"$i"*bai 2>/dev/null | wc -l` == 0 ]; then
				echo "BAM indexing $i"
				samtools index "$temppath"/"$i"*bam
			fi
			if [ -f "$endir"/"$runam"/bamlist.txt ]; then
				echo "$temppath"/"$i"*bam >> "$endir"/"$runam"/bamlist.txt
			else
				echo "$temppath"/"$i"*bam > "$endir"/"$runam"/bamlist.txt
			fi
			if [ -f "$endir"/"$runam"/samplelist.txt ]; then
				echo "$i" >> "$endir"/"$runam"/samplelist.txt
			else
				echo "$i" > "$endir"/"$runam"/samplelist.txt
			fi
			rm -rf "$endir"/"$runam"/"$i"/tmp
		done 2> "$endir"/"$runam"/error_mapdam.log
		shopt -u nullglob
		rm -rf "$endir"/"$runam"/GENOFILES
		mkdir "$endir"/"$runam"/GENOFILES
		echo "Preparing BED file"
		samtools idxstats "$temppath"/"$i"*bam | cut -f 1 | uniq | head -25 > "$endir"/"$runam"/"$i"/RUNLOGS/refnames.txt
		Rscript --vanilla ./src/bedmaker.r $eigen "$endir"/"$runam"/GENOFILES/ $psudo $cotrv "$endir"/"$runam"/"$i"/RUNLOGS/refnames.txt
		if [ $psudo == "Pseudohaploid_Random" ]; then
			echo "Running genotype call: pseudohaploid genome by random allele selection"
			samtools mpileup -R -B -q $mapqt -Q $bqthr -l "$endir"/"$runam"/GENOFILES/genocall.bed -f $refer -b "$endir"/"$runam"/bamlist.txt | sed 's/chr//' | pileupCaller --randomHaploid --sampleNameFile "$endir"/"$runam"/samplelist.txt --samplePopName "$runam" -f "$eigen" -e "$endir"/"$runam"/GENOFILES/"$runam"_geno
		else if [ $psudo == "Pseudohaploid_Majority" ]; then
			echo "Running genotype call: pseudohaploid genome by majority rule allele selection"
			samtools mpileup -R -B -q $mapqt -Q $bqthr -l "$endir"/"$runam"/GENOFILES/genocall.bed -f $refer -b "$endir"/"$runam"/bamlist.txt | sed 's/chr//' | pileupCaller --majorityCall --sampleNameFile "$endir"/"$runam"/samplelist.txt --samplePopName "$runam" -f "$eigen" -e "$endir"/"$runam"/GENOFILES/"$runam"_geno
		else if [ $psudo == "Diploid_Random" ]; then
			echo "Running genotype call: diploid genome by random reads selection"
			samtools mpileup -R -B -q $mapqt -Q $bqthr -l "$endir"/"$runam"/GENOFILES/genocall.bed -f $refer -b "$endir"/"$runam"/bamlist.txt | sed 's/chr//' | pileupCaller --randomDiploid --sampleNameFile "$endir"/"$runam"/samplelist.txt --samplePopName "$runam" -f "$eigen" -e "$endir"/"$runam"/GENOFILES/"$runam"_geno
		else if [ $psudo == "Diploid_GL" ]; then
			echo "Running genotype call: diploid genome by genotype likelihoods"
			angsd -gl 2 -out "$endir"/"$runam"/GENOFILES/genocall_GL -nThreads $thred -doGlf 2 -minQ $bqthr -minMapQ $mapqt -b "$endir"/"$runam"/bamlist.txt -rf "$endir"/"$runam"/GENOFILES/genocall.bed -doMajorMinor 1 -doCounts 1 -trim $ctrim
			Rscript --vanilla ./src/genocall_gl.r "$endir"/"$runam"/GENOFILES/genocall_GL.beagle.gz "$endir"/"$runam"/GENOFILES/genocall.bed "$endir"/"$runam"/samplelist.txt "$eigen" "$endir"/"$runam"/GENOFILES/genocall_altref.bed "$eigen" "$runam" "$endir"/"$runam"/GENOFILES/
		fi
		if [ $plink == TRUE ] || [ $plink == T ] || [ $plink == true ] || [ $plink == t ] || [ $plink == True ]; then
			echo "genotypename: $endir/$runam/GENOFILES/$runam_geno.geno" > ./tmp/paramfile
			echo "snpname: $endir/$runam/GENOFILES/$runam_geno.snp" >> ./tmp/paramfile
			echo "indivname: $endir/$runam/GENOFILES/$runam_geno.ind" >> ./tmp/paramfile
			echo "outputformat: PACKEDPED" >> ./tmp/paramfile
			echo "genotypeoutname: $endir/$runam/GENOFILES/$runam_geno.bed" >> ./tmp/paramfile
			echo "snpoutname: $endir/$runam/GENOFILES/$runam_geno.bim" >> ./tmp/paramfile
			echo "indivoutname: $endir/$runam/GENOFILES/$runam_geno.fam" >> ./tmp/paramfile
			convertf -p ./tmp/paramfile
		fi
	fi
}

function reportmake {
	
}
