#!/bin/bash
# 20220525

echo "run_plink2.sh -h"
func() {
	echo "Usage:"
	echo "run_plink2.sh -i <trio vcf> -o <output dir> -p <prefix>"
	echo "Description"
	echo "trio vcf, vcf with multi samples"
	echo "output dir, output directory"
	echo "prefix, output prefix"
	exit -1
}
while getopts "i:o:p:" OPT; do
	case $OPT in
		i) Input_Vcf="$OPTARG";;
		o) Output_Dir="$OPTARG";;
		p) Sample_ID="$OPTARG";;
		h) func;;
		?) func;;
	esac
done

#
# Input_Vcf=$1
# Output_Dir=$2
# Sample_ID=$3
#

###########################
plink2="/yk/apps/biosoft/plink_linux/plink2"
###########################

if [ ! -d $Output_Dir ]; then
    mkdir $Output_Dir
fi

# Kinship #
$plink2 \
	--vcf $Input_Vcf \
	-out $Output_Dir/$Sample_ID \
	--make-bed --make-king-table

rm -f $Output_Dir/${Sample_ID}.kin.temp
touch $Output_Dir/${Sample_ID}.kin.temp
echo "Kinship" > $Output_Dir/${Sample_ID}.kin.temp
cut -f 6 $Output_Dir/${Sample_ID}.kin0 | sed '1d' | while read line;
do
	if [ `echo "$line >= 0.354"|bc` -eq 1 ]; then
		echo "Duplicates/Monozygotic-Twins" >> $Output_Dir/${Sample_ID}.kin.temp
	elif [ `echo "$line >= 0.177"|bc` -eq 1 ]; then
		echo "Parent-Child/Full-Siblings" >> $Output_Dir/${Sample_ID}.kin.temp
	elif [ `echo "$line >= 0.0884"|bc` -eq 1 ]; then
		echo "2nd-Degree" >> $Output_Dir/${Sample_ID}.kin.temp
	elif [ `echo "$line >= 0.0442"|bc` -eq 1 ]; then
		echo "3rd-Degree" >> $Output_Dir/${Sample_ID}.kin.temp
	else
		echo "Unrelated" >> $Output_Dir/${Sample_ID}.kin.temp
	fi
done
paste $Output_Dir/${Sample_ID}.kin0 $Output_Dir/${Sample_ID}.kin.temp > $Output_Dir/${Sample_ID}.kinship
rm -f $Output_Dir/${Sample_ID}.kin.temp
