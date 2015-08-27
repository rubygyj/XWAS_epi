#!/bin/sh

# runs full qc for XWAS and autosomal association analysis
# assumes bed format [might change this later]
# adding default male female seperate QC

paramfile=$1

## Program will now parse relevant arguments from the param file ###
#------------------------------------------------------------------#
echo "1_Parsing arguments"

fname=$( awk '$1=="filename" {print $2}' $paramfile ) # name of dataset file (without the extension)
xwasloc=$( awk '$1=="xwasloc" {print $2}' $paramfile ) # location of extended plink executable
eigstratloc=$( awk '$1=="eigstratloc" {print $2}' $paramfile ) # smartpca and convertf file locations. 
exclind=$( awk '$1=="exclind" {print $2}' $paramfile ) # 0 if no file exists, 1 if file exists
xchrpca=$( awk '$1=="excludexchrPCA" {print $2}' $paramfile ) # YES to exclude xchr data when calculating PCA
bld=$( awk '$1=="build" {print $2}' $paramfile ) # build of the dataset (18-19)
alpha=$( awk '$1=="alpha" {print $2}' $paramfile ) # alpha threshold to set
binaryplink=$( awk '$1=="plinkformat" {print $2}' $paramfile ) # ped = ped/map format, bed = binary plink format
minmaf=$( awk '$1=="maf" {print $2}' $paramfile ) # remove variants with MAF < minmaf
mindthresh=$( awk '$1=="missindiv" {print $2}' $paramfile ) # filter for missingness per individual
genothresh=$( awk '$1=="missgeno" {print $2}' $paramfile ) # filter for missingness per genotype
numpc=$( awk '$1=="numpc" {print $2}' $paramfile ) # number of PC's 
pithresh=$( awk '$1=="related" {print $2}' $paramfile ) # threshold at which to remove related individuals
quant=$( awk '$1=="quant" {print $2}' $paramfile ) #is this a quantitative trait (either 0 or 1)

if [ $eigstratloc != "." ]
then
	cp $eigstratloc/smartpca $eigstratloc/convertf ./
fi
if [ $exclind -eq 0 ] 
then
	echo -n "" > ${fname}_exclind.remove
	cp ${fname}_exclind.remove ${fname}_male_exclind.remove
	cp ${fname}_exclind.remove ${fname}_female_exclind.remove
fi

if [ $binaryplink = "ped" ] 
then
	${xwasloc}/xwas --file ${fname} --make-bed --out ${fname} --silent
fi

#------------------------------------------------------------------#
#------------------------------------------------------------------#
totsnps=$(wc -l ${fname}.bim | awk '{print $1}')
bonf=$(echo "scale=20;$alpha/$totsnps" | bc)

# Setting PAR locations
if [ $bld -eq 19 ]
then

	echo  '23 60001 2699520 par1\n' > pars.txt
	echo  '23 154931044 155260560 par2\n' >> pars.txt

fi

if [ $bld -eq 18 ]
then
	echo  '23 1 2709520 par1\n' > pars.txt
	echo  '23 154584238 154913754 par2\n' >> pars.txt

fi

#------------------- Quality Controlling ---------------------------#
#-------------------------------------------------------------------#

# remove pseudo-autosomal regions
echo "Removing pseudo-autosomal regions"
${xwasloc}/xwas --bfile ${fname} --make-bed --out ${fname}_psar --exclude pars.txt --range --silent

if [ $quant -eq 1 ]
then 
	echo "Quality Control for Quantitative Traits"
	echo "Quality Control for Male"
	# extract male file
	${xwasloc}/xwas --bfile ${fname}_psar --filter-males --make-bed --out ${fname}_male --silent
	
	# HWEf
	echo "HWE"
	${xwasloc}/xwas --bfile ${fname}_male --hardy --out ${fname}_male_hwe --filter-controls --silent 
	cat ${fname}_male_hwe.hwe | grep ALL | awk -v bf=$bonf '$9<bf {print $2}' > ${fname}_male_snp.exclude

	# MAF, missingness per snp/ per individual
	echo "MAF, missingness per SNP, missingness per individual"
	${xwasloc}/xwas --bfile ${fname}_male --make-bed --out ${fname}_male_qc1 --mind ${mindthresh} --maf ${minmaf} --geno ${genothresh} --exclude ${fname}_male_snp.exclude --silent

	# check for wrong-sex individuals (do automatically if > 1000 snps on the X chromosome) Otherwise do this manually
	echo "Checking for wrong-sex individuals"
	xtot=`awk '$1==23 {print $0}' ./${fname}_male_qc1.bim | wc -l`
	if [ $xtot -gt 999 ] 
	then
		${xwasloc}/xwas --bfile ${fname}_male_qc1 --check-sex --out ${fname}_male_sexcheck --silent
		cat ${fname}_male_sexcheck.sexcheck | grep PROBLEM | awk '{print $1,$2,"sex"}' > ${fname}_male_exclind.remove
	fi
	
	# ld prune dataset
	echo "Checking relatedness"
	echo  '6 20000000 80000000 HLA\n' > hla.excl # generous buffer region for the HLA region
	${xwasloc}/xwas --bfile ${fname}_male_qc1 --indep-pairwise 50 5 0.3 --out ${fname}_male_qc1 --remove ${fname}_male_exclind.remove  --exclude hla.excl --range  --silent 
	allsnps=$(wc -l ${fname}_male_qc1.bim | awk '{print $1}') 
	thinamt=$(echo "scale=3;150000/$allsnps" | bc)

	# check for related individuals 
	${xwasloc}/xwas --bfile ${fname}_male_qc1 --thin $thinamt --genome --min ${pithresh} --extract ${fname}_male_qc1.prune.in --out ${fname}_male_rel --remove ${fname}_male_exclind.remove --silent

	# pihat > 0.12 remove
	awk 'NR>1 {print $0}' ${fname}_male_rel.genome | awk  '$10>$pithresh {print $1,$2,"related"}' >> ${fname}_male_exclind.remove
	${xwasloc}/xwas --bfile ${fname}_male_qc1 --make-bed --out ${fname}_male_qc2 --remove ${fname}_male_exclind.remove --silent

	echo "Quality Control for Female"
	# extract male file
	${xwasloc}/xwas --bfile ${fname}_psar --filter-females --make-bed --out ${fname}_female --silent
	
	# HWEf
	echo "HWE"
	${xwasloc}/xwas --bfile ${fname}_female --hardy --out ${fname}_female_hwe --filter-controls --silent 
	cat ${fname}_female_hwe.hwe | grep ALL | awk -v bf=$bonf '$9<bf {print $2}' > ${fname}_female_snp.exclude

	# MAF, missingness per snp/ per individual
	echo "MAF, missingness per SNP, missingness per individual"
	${xwasloc}/xwas --bfile ${fname}_female --make-bed --out ${fname}_female_qc1 --mind ${mindthresh} --maf ${minmaf} --geno ${genothresh} --exclude ${fname}_female_snp.exclude --silent

	# check for wrong-sex individuals (do automatically if > 1000 snps on the X chromosome) Otherwise do this manually
	echo "Checking for wrong-sex individuals"
	xtot=`awk '$1==23 {print $0}' ./${fname}_female_qc1.bim | wc -l`
	if [ $xtot -gt 999 ] 
	then
		${xwasloc}/xwas --bfile ${fname}_female_qc1 --check-sex --out ${fname}_female_sexcheck --silent
		cat ${fname}_female_sexcheck.sexcheck | grep PROBLEM | awk '{print $1,$2,"sex"}' > ${fname}_female_exclind.remove
	fi
	
	# ld prune dataset
	echo "Checking relatedness"
	echo '6 20000000 80000000 HLA\n' > hla.excl # generous buffer region for the HLA region
	${xwasloc}/xwas --bfile ${fname}_female_qc1 --indep-pairwise 50 5 0.3 --out ${fname}_female_qc1 --remove ${fname}_female_exclind.remove  --exclude hla.excl --range  --silent 
	allsnps=$(wc -l ${fname}_female_qc1.bim | awk '{print $1}') 
	thinamt=$(echo "scale=3;150000/$allsnps" | bc)

	# check for related individuals 
	${xwasloc}/xwas --bfile ${fname}_female_qc1 --thin $thinamt --genome --min ${pithresh} --extract ${fname}_female_qc1.prune.in --out ${fname}_female_rel --remove ${fname}_female_exclind.remove --silent

	# pihat > 0.12 remove
	awk 'NR>1 {print $0}' ${fname}_female_rel.genome | awk  '$10>$pithresh {print $1,$2,"related"}' >> ${fname}_female_exclind.remove
	${xwasloc}/xwas --bfile ${fname}_female_qc1 --make-bed --out ${fname}_female_qc2 --remove ${fname}_female_exclind.remove --silent

else 
	echo "Quality Control for Case/Control Traits"
	echo "Quality Control for Male"

	# extract male file
	${xwasloc}/xwas --bfile ${fname}_psar --filter-males --make-bed --out ${fname}_male --silent
	
	echo "HWE"
	${xwasloc}/xwas --bfile ${fname}_male --hardy --out ${fname}_male_hwe --filter-controls --silent 
	cat ${fname}_male_hwe.hwe | grep ALL | awk -v bf=$bonf '$9<bf {print $2}' > ${fname}_male_snp.exclude
	
	# Missingness correlated with phenotype
	echo "Correlation between missingness and phenotype"
	${xwasloc}/xwas --bfile ${fname}_male --test-missing --out ${fname}_male_mcc --silent  
	awk -v bf=$bonf '$5<bf {print $2}' ${fname}_male_mcc.missing >> ${fname}_male_snp.exclude

	# MAF, missingness per snp/ per individual
	echo "MAF, missingness per SNP, missingness per individual"
	${xwasloc}/xwas --bfile ${fname}_male --make-bed --out ${fname}_male_qc1 --mind ${mindthresh} --maf ${minmaf} --geno ${genothresh} --exclude ${fname}_male_snp.exclude --silent

	# check for wrong sex individuals (do automatically if > 1000 snps on the X chromosome) Otherwise do this manually
	echo "Checking for wrong-sex individuals"
	xtot=`awk '$1==23 {print $0}' ./${fname}_male_qc1.bim | wc -l`
	if [ $xtot -gt 999 ] 
	then
#		echo "Sex check (# X-snps > 999)" 
		${xwasloc}/xwas --bfile ${fname}_male_qc1 --check-sex --out ${fname}_male_sexcheck --silent
		cat ${fname}_male_sexcheck.sexcheck | grep PROBLEM | awk '{print $1,$2,"sex"}' > ${fname}_male_exclind.remove
	fi

	# ld prune dataset
	echo "Checking relatedness"
	echo '6 20000000 80000000 HLA\n' > hla.excl # generous buffer region for the HLA region
	${xwasloc}/xwas --bfile ${fname}_male_qc1 --indep-pairwise 50 5 0.3 --out ${fname}_male_qc1 --remove ${fname}_male_exclind.remove  --exclude hla.excl --range  --silent 
	allsnps=$(wc -l ${fname}_male_qc1.bim | awk '{print $1}') 
	thinamt=$(echo "scale=3;150000/$allsnps" | bc)

	# check for related individuals 
	#	echo "check for related individuals"
	${xwasloc}/xwas --bfile ${fname}_male_qc1 --thin $thinamt --genome --min ${pithresh} --extract ${fname}_male_qc1.prune.in --out ${fname}_male_rel --remove ${fname}_male_exclind.remove --silent
	
	# pihat > 0.12 remove
	awk 'NR>1 {print $0}' ${fname}_male_rel.genome | awk  '$10>$pithresh {print $1,$2,"related"}' >> ${fname}_male_exclind.remove
	${xwasloc}/xwas --bfile ${fname}_male_qc1 --make-bed --out ${fname}_male_qc2 --remove ${fname}_male_exclind.remove --silent
	
	#	echo "ends for man qualitative"

	echo "Quality Control for Female"
	# extract male file
	#	echo "extract female file"
	${xwasloc}/xwas --bfile ${fname}_psar --filter-females --make-bed --out ${fname}_female --silent
	
	echo "HWE"
	${xwasloc}/xwas --bfile ${fname}_female --hardy --out ${fname}_female_hwe --filter-controls --silent 
	cat ${fname}_female_hwe.hwe | grep ALL | awk -v bf=$bonf '$9<bf {print $2}' > ${fname}_female_snp.exclude
	
	# Missingness correlated with phenotype
	echo "Correlation between missingness and phenotype"
	${xwasloc}/xwas --bfile ${fname}_female --test-missing --out ${fname}_female_mcc --silent  
	awk -v bf=$bonf '$5<bf {print $2}' ${fname}_female_mcc.missing >> ${fname}_female_snp.exclude

	# MAF, missingness per snp/ per individual
	echo "MAF, missingness per SNP, missingness per individual"
	${xwasloc}/xwas --bfile ${fname}_female --make-bed --out ${fname}_female_qc1 --mind ${mindthresh} --maf ${minmaf} --geno ${genothresh} --exclude ${fname}_female_snp.exclude --silent

	# check for wrong sex individuals (do automatically if > 1000 snps on the X chromosome) Otherwise do this manually
	echo "Checking wrong-sex individuals"
	xtot=`awk '$1==23 {print $0}' ./${fname}_female_qc1.bim | wc -l`
	if [ $xtot -gt 999 ] 
	then
#		echo "Sex check (# X-snps > 999)" 
		${xwasloc}/xwas --bfile ${fname}_female_qc1 --check-sex --out ${fname}_female_sexcheck --silent
		cat ${fname}_female_sexcheck.sexcheck | grep PROBLEM | awk '{print $1,$2,"sex"}' > ${fname}_female_exclind.remove
	fi

	# ld prune dataset
	echo "Checking relatedness"
	echo  '6 20000000 80000000 HLA\n' > hla.excl # generous buffer region for the HLA region
	${xwasloc}/xwas --bfile ${fname}_female_qc1 --indep-pairwise 50 5 0.3 --out ${fname}_female_qc1 --remove ${fname}_female_exclind.remove  --exclude hla.excl --range  --silent 
	allsnps=$(wc -l ${fname}_female_qc1.bim | awk '{print $1}') 
	thinamt=$(echo "scale=3;150000/$allsnps" | bc)

	# check for related individuals 
#	echo "check for related individuals"
	${xwasloc}/xwas --bfile ${fname}_female_qc1 --thin $thinamt --genome --min ${pithresh} --extract ${fname}_female_qc1.prune.in --out ${fname}_female_rel --remove ${fname}_female_exclind.remove --silent
	
	# pihat > 0.12 remove
	awk 'NR>1 {print $0}' ${fname}_female_rel.genome | awk  '$10>$pithresh {print $1,$2,"related"}' >> ${fname}_female_exclind.remove
	${xwasloc}/xwas --bfile ${fname}_female_qc1 --make-bed --out ${fname}_female_qc2 --remove ${fname}_female_exclind.remove --silent
	
#	echo "ends for female qualitative"
fi

# Merge male and female qc1
echo "Merge male and female qc1"
${xwasloc}/xwas --bfile ${fname}_female_qc1 --bmerge ${fname}_male_qc1.bed ${fname}_male_qc1.bim ${fname}_male_qc1.fam --merge-mode 6 --out ${fname}_sex_diff_snps --silent

${xwasloc}/xwas --bfile ${fname}_female_qc1 --exclude ${fname}_sex_diff_snps.diff --make-bed --out ${fname}_female_qc1_tmp --silent
${xwasloc}/xwas --bfile ${fname}_male_qc1 --exclude ${fname}_sex_diff_snps.diff --make-bed --out ${fname}_male_qc1_tmp --silent

${xwasloc}/xwas --bfile ${fname}_female_qc1_tmp --bmerge ${fname}_male_qc1_tmp.bed ${fname}_male_qc1_tmp.bim ${fname}_male_qc1_tmp.fam --make-bed --out ${fname}_qc1 --silent

# Merge male and female qc2
echo "Merge male and female qc2"
${xwasloc}/xwas --bfile ${fname}_female_qc2 --bmerge ${fname}_male_qc2.bed ${fname}_male_qc2.bim ${fname}_male_qc2.fam --merge-mode 6 --out ${fname}_sex_diff_snps --silent

${xwasloc}/xwas --bfile ${fname}_female_qc2 --exclude ${fname}_sex_diff_snps.diff --make-bed --out ${fname}_female_qc2_tmp --silent
${xwasloc}/xwas --bfile ${fname}_male_qc2 --exclude ${fname}_sex_diff_snps.diff --make-bed --out ${fname}_male_qc2_tmp --silent

${xwasloc}/xwas --bfile ${fname}_female_qc2_tmp --bmerge ${fname}_male_qc2_tmp.bed ${fname}_male_qc2_tmp.bim ${fname}_male_qc2_tmp.fam --make-bed --out ${fname}_qc2 --silent

# Merge male female exclind.remove and prune.in
cat ${fname}_male_exclind.remove ${fname}_female_exclind.remove > ${fname}_exclind.remove
cat ${fname}_male_qc1.prune.in ${fname}_female_qc1.prune.in > ${fname}_qc1.prune.in

# PC
echo "Population stratification"
${xwasloc}/xwas --bfile ${fname}_qc2 --extract ${fname}_qc1.prune.in --recode12 --out ${fname}_qc2 --silent
# convert to eigenstrat format
${eigstratloc}/convertf.perl ${fname}_qc2 
# now get the PC's using eigenstrat
${eigstratloc}/smartpca_loadings_mod.perl -i ${fname}_qc2.eigenstratgeno -a ${fname}_qc2.snp -b ${fname}_qc2.ind -k ${numpc} -o ${fname}_qc2.pca -p ${fname}_qc2.plot -e ${fname}_qc2.eval -l ${fname}_qc2.log -d ${fname}_qc2.load -m 5 -t 10 -s 6.0 -g ${xchrpca} > ${fname}_qc2.eig.log
# make the covariate file
awk 'BEGIN {OFS="\t"} FNR==NR {f1[$2]=$1;next} ($1 in f1)  {$NF=""} {print f1[$1],$0}' ${fname}_qc1.fam ${fname}_qc2.pca.evec | grep -v eigvals: > ${fname}.covar
# copying outlier individuals to the individuals to remove file
cat ${fname}_qc2.log | grep REMOVED | awk '{print $3}' > temp.a
awk 'FNR==NR {f1[$0]; next} ($2 in f1) {print $1,$2,"pca"}' ./temp.a ${fname}_qc2.fam >> ${fname}_exclind.remove


# Significant diff in MAF between males and females
echo "Significant difference in MAF between males and females"
${xwasloc}/xwas --bfile ${fname}_qc1 --remove ${fname}_exclind.remove --out ${fname}_final --make-bed --silent
totx=$(awk '$1=="X" || $1==23 {print $0}' ${fname}_final.bim | wc -l)
bonfx=$(echo "scale=20;$alpha/$totx" | bc)
${xwasloc}/xwas --bfile ${fname}_final --xwas --make-bed --out ${fname}_final_maf_x --chr X --freqdiff-x ${bonfx} --silent

echo "Significant difference in Missingness between males and females"
totx=$(awk '$1=="X" || $1==23 {print $0}' ${fname}_final_maf_x.bim | wc -l)
bonfx=$(echo "scale=20;$alpha/$totx" | bc)
${xwasloc}/xwas --bfile ${fname}_final_maf_x --xwas --make-bed --out ${fname}_final_x --chr X --missdiff-x ${bonfx} --silent

# ask about removing more if possible, because this still leaves a lot more in there
rm *eigenstratgeno *ind *snp *temp* *_qc1.bed *_qc1.fam *_qc1.bim *_qc2.bed *_qc2.bim *_qc2.fam pars.txt *male_qc1_tmp* *male_qc2_tmp* *psar.* *_final_maf_x.*
if [ $eigstratloc != "." ]
then
	rm smartpca convertf
fi
