#!/usr/bin/env bash

# Version 2.6 for use with large genome pooled resequencing projects

export LC_ALL=en_US.UTF-8
export SHELL=bash

#mtDNA Chromosome
mtDNA="NC_007175.2"

# This script uses the individual raw.bcf files for filtering instead of a single VCF file.  This enables parallelization and much faster processing of larger files

# Usage is `PVCF_filter PREFIX PROCESSORS DEPTH_CUTFOFF`

if [[ -z "$2" ]]; then
echo "PVCF_filter PREFIX PROCESSORS DEPTH_CUTFOFF"
exit 1
fi

if [[ -z "$3" ]]; then
echo "PVCF_filter PREFIX PROCESSORS DEPTH_CUTFOFF"
exit 1
fi

NumProc=$2
PREFIX=$1
SProc=$(($NumProc / 5))
DEPTH=$3

NumInd=$(ls raw.*.vcf.gz 2> /dev/null | wc -l)
NumInd=$(($NumInd - 0))


if [ "$NumInd" -gt 0 ];then
	GZIP="TRUE"
fi

if [ "$GZIP" == "TRUE" ]; then
	ls raw.*.vcf.gz | sed 's/raw.//g' | sed 's/.vcf.gz//g' > list
	#cat list | parallel --no-notice -j $NumProc "vcftools --gzvcf raw.{}.vcf.gz --minQ 30 --recode --recode-INFO-all --max-missing 0.25 --out $PREFIX.TRS.{} 2> $PREFIX.filter.errors"
	#cat list | parallel --no-notice -j $NumProc "bcftools view -i 'F_MISSING<0.75 && QUAL > 29' raw.{}.vcf.gz -O v | vcftools --vcf - --minDP 5 --recode --recode-INFO-all --stdout 2> $PREFIX.filter.errors | vcftools --vcf - --max-missing 0.75 --maf 0.0001 --recode --recode-INFO-all --stdout 2>> $PREFIX.filter.errors| bgzip -c > $PREFIX.TRSdp5g75.{}.recode.vcf.gz 2>> $PREFIX.filter.errors"
	cat list | bcftools view -i 'F_MISSING<0.75 && QUAL > 29' raw.larvcado.vcf.gz -O v | bcftools +setGT -- -t q -n . -i 'FORMAT/DP<5' 2> $PREFIX.filter.errors | bcftools view -i 'F_MISSING<0.25 && MAF > 0.0001' -O z -o $PREFIX.TRSdp.$DEPTH.g75.larvcado.recode.vcf.gz 2>> $PREFIX.filter.errors

else
	NumInd=$(ls Collated.raw.*.bcf 2> /dev/null | wc -l)
	NumInd=$(($NumInd - 0))
	NumK=$(($NumInd/1000))
	ls Collated.raw.*.bcf | sed 's/raw.//g' | sed 's/.bcf//g' > list
	#cat list | parallel --no-notice --no-notice -j $NumProc "vcftools --vcf raw.{}.vcf --minQ 30 --recode --recode-INFO-all --max-missing 0.25 --out $PREFIX.TRS.{} 2> /dev/null"
	cat list | parallel --no-notice -j $NumProc "bcftools view -i 'F_MISSING<0.75' Collated.raw.{}.bcf  | bcftools +setGT -- -t q -n . -i 'FORMAT/DP<"$DEPTH"' 2> $PREFIX.filter.errors | bcftools view -i 'F_MISSING<0.25 && INFO/AO/INFO/RO > 0.015' -O z -o $PREFIX.TRSdp.$DEPTH.g75.{}.recode.vcf.gz 2>> $PREFIX.filter.errors"
	
fi

if [ ! -d "raw" ]; then
	mkdir raw
fi

mv raw.larvcado.vcf.gz Collated.raw.*.bcf ./raw 2> /dev/null



cat list | zgrep -v $mtDNA $PREFIX.TRSdp.$DEPTH.g75.larvcado.recode.vcf.gz | bgzip -c > $PREFIX.TRSdp.$DEPTH.g75.nDNA.larvcado.vcf.gz


if [ ! -d "filtered" ]; then
	mkdir filtered
fi


if [ ! -d "TRSdp20g75" ]; then
        mkdir TRSdp.$DEPTH.g75
fi

mv $PREFIX.TRSdp.$DEPTH.g75.larvcado.recode.vcf.gz*  TRSdp.$DEPTH.g75
   
#nDNA processing
echo "This script will automatically filter a FreeBayes generated VCF file using criteria related to site depth,"
echo "quality versus depth, allelic balance at heterzygous individuals, and paired read representation."
echo -e "The script assumes that loci and individuals with low call rates (or depth) have already been removed. \n"
echo -e "Contact Jon Puritz (jpuritz@gmail.com) for questions and see script comments for more details on particular filters \n"

NEW_PREFIX="$PREFIX.TRSdp.$DEPTH.g75.nDNA"

#Creates a file with the original site depth and qual for each locus

AWK1='!/#/ {print $1 "\t" $2 "\t" $6}'
AWK2='!/#/ {print $1}'
AWK3='!/NP/ && !/#/'
AWK4='!/#/ {print $1 "\t" $2}'

#seq 0 $NumK | parallel --no-notice -k -j $NumProc "zcat $NEW_PREFIX.{}.vcf.gz | cut -f8  | grep -P -oe 'DP=[0-9]*' | sed -s 's/DP=//g' " > TEMP.{}.DEPTH
#cat TEMP.*.DEPTH > $NEW_PREFIX.DEPTH
#rm TEMP*.DEPTH
zcat $NEW_PREFIX.larvcado.vcf.gz | grep -v "^#" | cut -f8 | grep -oe 'DP=[0-9]*' | sed -s 's/DP=//g' > $NEW_PREFIX.DEPTH

#seq 0 $NumK  | parallel --no-notice -k -j $NumProc "zcat $NEW_PREFIX.{}.vcf.gz | mawk '$AWK1' " > TEMP.{}.loci.qual
#cat TEMP.*.loci.qual > $NEW_PREFIX.loci.qual
#rm TEMP.*.loci.qual
zcat $NEW_PREFIX.larvcado.vcf.gz | mawk '$AWK1' > $NEW_PREFIX.loci.qual


OLD=$( cat list | zcat $NEW_PREFIX.larvcado.vcf.gz | mawk '$AWK2' | wc -l | mawk '{sum = sum + $1} END {print sum}' )

if [ "$PE" != "yes" ]; then
	FILTERED2=$(cat list | zcat $NEW_PREFIX.larvcado.vcf.gz | vcffilter -f 'PAIRED > 0.05 & PAIREDR > 0.05 & PAIREDR / PAIRED < 5 & PAIREDR / PAIRED > 0.05 | PAIRED < 0.05 & PAIREDR < 0.05' -s | mawk '$AWK2' | wc -l | mawk '{sum = sum + $1} END {print sum}' )
	NUMFIL2=$(($OLD - $FILTERED2))
	echo -e "Number of sites filtered based on properly paired status\n" $NUMFIL2 "of" $OLD "\n"
	echo -e "Number of sites filtered based on properly paired status\n" $NUMFIL2 "of" $OLD "\n" >> $NEW_PREFIX.filterstats
else	
	FILTERED2=$(cat list | zcat $NEW_PREFIX.larvcado.vcf.gz | vcffilter -f 'PAIRED < 0.005 & PAIREDR > 0.005 | PAIRED > 0.005 & PAIREDR < 0.005' -t NP -F PASS -A | mawk '$AWK3' | wc -l | mawk '{sum = sum + $1} END {print sum}' )
	NUMFIL2=$(($OLD - $FILTERED2))
	echo -e "Number of sites filtered based on properly paired status\n" $NUMFIL2 "of" $OLD "\n"
	echo -e "Number of sites filtered based on properly paired status\n" $NUMFIL2 "of" $OLD "\n" >> $NEW_PREFIX.filterstats

fi

wait
wait
wait

#Calculates the average depth and three times the square root of mean depth
DEPTH=$(mawk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' $NEW_PREFIX.DEPTH)
MD=$(python -c "import math; print( round(math.sqrt($DEPTH)))")
MD=$(python -c "import math; print( round( $MD * 3))")
DEPTH=$(python -c "import math; print( round($DEPTH + $MD))")

#Filters loci above the mean depth + 1 standard deviation that have quality scores that are less than 1*DEPTH
paste $NEW_PREFIX.loci.qual $NEW_PREFIX.DEPTH | mawk -v x=$DEPTH '$4 > x'| mawk -v x=$DEPTH '$3 < 1 * x' > $NEW_PREFIX.lowQDloci

SITES=$(cat list | vcftools --gzvcf $NEW_PREFIX.larvcado.vcf.gz --exclude-positions $NEW_PREFIX.lowQDloci 2>&1 | grep Sites | cut -f4 -d " " | mawk '{sum = sum + $1} END {print sum}' )
LQDL=$(( $OLD - $SITES ))

echo -e "Number of sites filtered based on high depth and lower than 2*DEPTH quality score\n" $LQDL "of" $OLD "\n"
echo -e "Number of sites filtered based on high depth and lower than 2*DEPTH quality score\n" $LQDL "of" $OLD "\n" >> $NEW_PREFIX.filterstats

#Recalculates site depth for sites that have not been previously filtered
if [ "$PE" != "yes" ]; then

	cat list | zcat $NEW_PREFIX.larvcado.vcf.gz | vcffilter -f 'PAIRED > 0.05 & PAIREDR > 0.05 & PAIREDR / PAIRED < 5 & PAIREDR / PAIRED > 0.05 | PAIRED < 0.05 & PAIREDR < 0.05' -s | vcftools --vcf - --remove-filtered NP --exclude-positions $NEW_PREFIX.lowQDloci --site-depth --out $NEW_PREFIX.larvcado 2>> $NEW_PREFIX.filter.errors
else
	cat list | zcat $NEW_PREFIX.larvcado.vcf.gz | vcffilter -f 'PAIRED < 0.005 & PAIREDR > 0.005 | PAIRED > 0.005 & PAIREDR < 0.005' -t NP -F PASS -A  | vcftools --vcf - --remove-filtered NP --exclude-positions $NEW_PREFIX.lowQDloci --site-depth --out $NEW_PREFIX.larvcado 2>> $NEW_PREFIX.filter.errors
fi

cut -f3 $NEW_PREFIX.larvcado.ldepth | mawk '!/SUM_DEPTH/' > $NEW_PREFIX.site.depth


cat list | gzip $NEW_PREFIX.larvcado.ldepth

DP=$(mawk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' $NEW_PREFIX.site.depth)
SD=$(mawk '{delta = $1 - avg; avg += delta / NR; mean2 += delta * ($1 - avg); } END { print sqrt(mean2 / NR); }' $NEW_PREFIX.site.depth)

#Calculates actual number of individuals in VCF file
#This is important because loci will now be filtered by mean depth calculated with individuals present in VCF
FFILE=$(head -1 list)
IND=$(zcat $NEW_PREFIX.$FFILE.vcf.gz |head -1000 | mawk '/#/' | tail -1 | wc -w)
IND=$(($IND - 9))

mawk '!/D/' $NEW_PREFIX.site.depth | mawk -v x=$IND '{print $1/x}' > meandepthpersite

#Calculates a mean depth cutoff to use for filtering
DP=$(perl -e "print ($DP+ 1.645*$SD) / $IND")
PP=$(mawk '!/SUM/' $NEW_PREFIX.site.depth | sort -rn | perl -e '$d=.05;@l=<>;print $l[int($d*$#l)]' )
PP=$(perl -e "print int($PP / $IND)")
GP=$(perl -e "print int($PP * 1.25)")
export GP

gnuplot << \EOF >> $NEW_PREFIX.filterstats
set terminal dumb size 120, 30
set autoscale
high=system("echo $GP")
set xrange [10:high]
unset label
set title "Histogram of mean depth per site"
set ylabel "Number of Occurrences"
set xlabel "Mean Depth"
#set yr [0:100000]
binwidth=1
bin(x,width)=width*floor(x/width) + binwidth/2.0
set xtics floor(high/20)
plot 'meandepthpersite' using (bin($1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF

gnuplot << \EOF
set terminal dumb size 120, 30
set autoscale
high=system("echo $GP")
set xrange [10:high]
unset label
set title "Histogram of mean depth per site"
set ylabel "Number of Occurrences"
set xlabel "Mean Depth"
#set yr [0:100000]
binwidth=1
bin(x,width)=width*floor(x/width) + binwidth/2.0
set xtics floor(high/20)
plot 'meandepthpersite' using (bin($1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF


if [[ -z "$4" ]]; then
echo "The 95% cutoff would be" $PP
echo "Would you like to use a different maximum mean depth cutoff than "$PP", yes or no"

read NEWCUTOFF
else
NEWCUTOFF=$4
fi

if [ "$NEWCUTOFF" != "yes" ]; then
echo -e "Maximum mean depth cutoff is" $PP
echo -e "Maximum mean depth cutoff is" $PP >> $NEW_PREFIX.filterstats

else
	if [[ -z "$5" ]]; then
		echo "Please enter new cutoff"
		read PP
	else
		PP=$5
	fi
echo -e "Maximum mean depth cutoff is" $PP >> $NEW_PREFIX.filterstats
fi

#Combines all filters to create filtered VCF files
if [ "$PE" != "yes" ]; then

	cat list | zcat $NEW_PREFIX.larvcado.vcf.gz | vcffilter -f 'PAIRED > 0.05 & PAIREDR > 0.05 & PAIREDR / PAIRED < 5 & PAIREDR / PAIRED > 0.05 | PAIRED < 0.05 & PAIREDR < 0.05' -s | vcftools --vcf - --remove-filtered NP --max-meanDP $PP --recode --exclude-positions $NEW_PREFIX.lowQDloci --recode-INFO-all --stdout 2>> $NEW_PREFIX.filter.errors | vcfstreamsort -a | bcftools view -O b -o $NEW_PREFIX.larvcado.FIL.bcf 2>> $NEW_PREFIX.filter.errors
else
	cat list | zcat $NEW_PREFIX.larvcado.vcf.gz | vcffilter -f 'PAIRED < 0.005 & PAIREDR > 0.005 | PAIRED > 0.005 & PAIREDR < 0.005' -t NP -F PASS -A  | vcftools --vcf - --remove-filtered NP --max-meanDP $PP --recode --exclude-positions $NEW_PREFIX.lowQDloci --recode-INFO-all --stdout 2>> $NEW_PREFIX.filter.errors | vcfstreamsort -a | bcftools view -O b -o $NEW_PREFIX.larvcado.FIL.bcf 2>> $NEW_PREFIX.filter.errors

fi


FILTERED3=$(cat list | bcftools view $NEW_PREFIX.larvcado.FIL.bcf -O v | mawk '$AWK2' | wc -l | mawk '{sum = sum + $1} END {print sum}' ) && OLD2=$(cat $NEW_PREFIX.site.depth | wc -l)

NUMFIL3=$(($FILTERED2 - $FILTERED3))

echo -e "Number of sites filtered based on maximum mean depth\n" $NUMFIL3 "\n"
echo -e "Number of sites filtered based on maximum mean depth\n" $NUMFIL3 "\n" >> $NEW_PREFIX.filterstats

NUMFIL4=$(($OLD - $FILTERED3))

echo -e "Total number of sites filtered\n" $NUMFIL4 "of" $OLD "\n"
echo -e "Total number of sites filtered\n" $NUMFIL4 "of" $OLD "\n" >> $NEW_PREFIX.filterstats

echo -e "Remaining sites\n" $FILTERED3 "\n"
echo -e "Remaining sites\n" $FILTERED3 "\n" >> $NEW_PREFIX.filterstats


if [ ! -d "nDNA" ]; then
	mkdir nDNA
	mkdir metrics
	mkdir SNPs
	mkdir nDNA.INDels
	mkdir nDNA.FIL
fi

mv $NEW_PREFIX.larvcado.vcf.gz* ./nDNA

echo -e "Variants will now be composed into SNPs and INDels\n"
echo -e "Variants will now be composed into SNPs and INDels\n" >> $NEW_PREFIX.filterstats


ln -s ../reference.fasta* .	

cat list | bcftools norm -m - -f reference.fasta $NEW_PREFIX.larvcado.FIL.bcf | vcfallelicprimitives -k -g | vcfstreamsort | bcftools norm -m + -f reference.fasta  | bcftools +fill-tags | bcftools view --threads 2 -V indels,other -O z -o SNP.$NEW_PREFIX.larvcado.FIL.bcf 2>> $NEW_PREFIX.filter.errors
cat list | bcftools norm -m - -f reference.fasta $NEW_PREFIX.larvcado.FIL.bcf | vcfallelicprimitives -k -g | vcfstreamsort | bcftools norm -m + -f reference.fasta  | bcftools +fill-tags | bcftools view --threads 2 -v indels,other -O z -o INDELS.$NEW_PREFIX.larvcado.bcf 2>> $NEW_PREFIX.filter.errors

#bcftools norm -m - -f reference.fasta CASE.TRSdp20g1.nDNA.vcf.gz | vcfallelicprimitives -k -g | vcffilter -f "TYPE = snp" | vcfstreamsort | bcftools norm -m + -f reference.fasta | vcftools --vcf - --max-alleles 2 --recod 

echo "Numk is " $NumK

ls $NEW_PREFIX.larvcado.FIL.bcf > bcf.larvcado.list
bcftools concat -n -f bcf.larvcado.list -O b -o $NEW_PREFIX.Collated.larvcado.FIL.bcf

ls $NEW_PREFIX.Collated.*.FIL.bcf > collated.bcf.list
bcftools concat -n -f collated.bcf.list -O b | bcftools view -O z --threads $NumProc -o $NEW_PREFIX.FIL.vcf.gz
rm $NEW_PREFIX.Collated.*.FIL.bcf


mv $NEW_PREFIX.*.FIL.bcf ./nDNA.FIL/ 2>> $NEW_PREFIX.filter.errors


FILTERED4=$(cat list | bcftools view SNP.$NEW_PREFIX.larvcado.FIL.bcf | mawk '$AWK2' | wc -l | mawk '{sum = sum + $1} END {print sum}' )
FILTERED5=$(cat list | bcftools view INDELS.$NEW_PREFIX.larvcado.bcf | mawk '$AWK2' | wc -l | mawk '{sum = sum + $1} END {print sum}' )

echo -e "Number of SNPs before filtering for call rate and minor allele frequency\n" $FILTERED4 "\n"
echo -e "Number of SNPs before filtering for call rate and minor allele frequency\n" $FILTERED4 "\n" >> $NEW_PREFIX.filterstats

echo -e "Number of INDels\n" $FILTERED5 "\n"
echo -e "Number of INDels\n" $FILTERED5 "\n" >> $NEW_PREFIX.filterstats

mv INDELS.$NEW_PREFIX.*.bcf ./nDNA.INDels 2>> $NEW_PREFIX.filter.errors

# Filter by genotype and minor allele frequency
cat list | bcftools view SNP.$NEW_PREFIX.larvcado.FIL.bcf 2>/dev/null | bcftools view -i 'F_MISSING = 0' -O b -o SNP.$NEW_PREFIX.g1.larvcado.FIL.bcf 2>> $NEW_PREFIX.filter.errors

FILTERED10=$(cat list | bcftools view SNP.$NEW_PREFIX.g1.larvcado.FIL.bcf | mawk '$AWK2' | wc -l | mawk '{sum = sum + $1} END {print sum}' )
echo -e "Number of SNPs with 100% call rate\n" $FILTERED10 "\n"
echo -e "Number of SNPs with 100% call rate\n" $FILTERED10 "\n" >> $NEW_PREFIX.filterstats

# Collate
ls SNP.$NEW_PREFIX.larvcado.FIL.bcf > bcf.larvcado.list
bcftools concat -n -f bcf.larvcado.list -O b -o SNP.$NEW_PREFIX.Collated.larvcado.FIL.bcf

ls SNP.$NEW_PREFIX.Collated.larvcado.FIL.bcf > collated.bcf.list
bcftools concat -n -f collated.bcf.list -O b | bcftools view -O z --threads $NumProc -o SNP.$NEW_PREFIX.FIL.vcf.gz

rm SNP.$NEW_PREFIX.Collated.larvcado.FIL.bcf

mv SNP.$NEW_PREFIX.larvcado.FIL.bcf SNP.$NEW_PREFIX.FIL.vcf.gz ./SNPs


#ls SNP.$NEW_PREFIX.g1.larvcado.FIL.bcf > bcf.larvcado.list
#bcftools concat -n -f bcf.larvcado.list -O b -o SNP.$NEW_PREFIX.g1.Collated.larvcado.FIL.bcf
#rm SNP.$NEW_PREFIX.g1.larvcado.FIL.bcf
#ls SNP.$NEW_PREFIX.g1.Collated.*.FIL.bcf > collated.bcf.list
#bcftools concat -n -f collated.bcf.list -O b 2>> $NEW_PREFIX.filter.errors | bcftools view -O z -o SNP.$NEW_PREFIX.g1.FIL.vcf.gz
#rm SNP.$NEW_PREFIX.g1.Collated.*.FIL.bcf
	
#mv SNP.$NEW_PREFIX.*.*.FIL.vcf.gz ./filtered

ls $NEW_PREFIX.DEPTH $NEW_PREFIX.lo* $NEW_PREFIX.ldepth* $NEW_PREFIX.site.depth meandepthpersite | gzip {}
mv $NEW_PREFIX.DEPTH.gz $NEW_PREFIX.lo* $NEW_PREFIX.site.depth.gz meandepthpersite.gz  ./metrics
cat $NEW_PREFIX.larvcado.ldepth.gz > $NEW_PREFIX.ldepth.gz 
mv $NEW_PREFIX.ldepth.gz ./metrics
rm $NEW_PREFIX.larvcado.ldepth.gz

mv $NEW_PREFIX.filterstats ./filtered

echo -e "Filter stats and filtered VCF files stored in $NEW_PREFIX.filterstats\n"
echo -e "Both are stored in the filtered directory\n"

