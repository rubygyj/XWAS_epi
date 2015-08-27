#!/bin/bash

# This will generate a series of datasets for the XWAS functions to be tested on
# These datasets will all be on the X Chromosome and will contain both 
# quantitative and case/control phenotypes. 

# Each of the datasets will contain approximately 500 individuals and 1000 SNPs on
# on the X Chromosome

# Bash function to recode the sex of the individuals.
recode_sex(){
	awk '{ if (rand() < 0.5) print $1"\t"$2"\t"$3"\t"$4"\t"1"\t"$6; else print $1"\t"$2"\t"$3"\t"$4"\t"2"\t"$6; }' $1
}

# Recode the chromosome to be 23
recode_as_x(){
	awk '{print "23\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' $1
}

# Just generates random covariates in bash
generate_covars(){
	awk '{print $1"\t"$2"\t"rand()}' $1 > $2
}


echo "Beginning Generation of Datasets!"
# In later versions we will ask the user to point this to a compiled Plink 1.07 (or 1.9) Binary, 
# since we will likely not include the simulation functions in the XWAS package in further releases
XWAS_ORIG=../bin/xwas

# Simulating quantitative datasets
echo "Generating dummy_quant1!"
touch dummy_quant1.sim
echo "1000 qtl 0.05 0.95  0.001 0" >> dummy_quant1.sim
${XWAS_ORIG} --simulate-qt dummy_quant1.sim --simulate-n 500 --make-bed --out dummy_quant1 --silent
generate_covars dummy_quant1.fam dummy_quant1.covar
recode_as_x dummy_quant1.bim > test.bim
mv test.bim dummy_quant1.bim
recode_sex dummy_quant1.fam > test.fam 
mv test.fam dummy_quant1.fam

echo "Generating dummy_quant2!"
touch dummy_quant2.sim
echo "500 qtl 0.10 0.90 0.002 1" >> dummy_quant2.sim
${XWAS_ORIG} --simulate-qt dummy_quant2.sim --simulate-n 500 --chr 23 --make-bed --out dummy_quant2 --silent
generate_covars dummy_quant2.fam dummy_quant2.covar
recode_as_x dummy_quant2.bim > test.bim
mv test.bim dummy_quant2.bim
recode_sex dummy_quant2.fam > test.fam 
mv test.fam dummy_quant2.fam

# Simulating case/control phenotypes
echo "Generating dummy_case1!"
touch dummy_case1.sim
echo "1000 null    0.05 1.00 1.00 1.00" >> dummy_case1.sim 
echo "100  disease 0.05 1.00 2.00 mult" >> dummy_case1.sim
${XWAS_ORIG} --simulate dummy_case1.sim --simulate-ncases 250 --simulate-ncontrols 1000  --chr 23 --make-bed --out dummy_case1 --silent
generate_covars dummy_case1.fam dummy_case1.covar
recode_as_x dummy_case1.bim > test.bim
mv test.bim dummy_case1.bim
recode_sex dummy_case1.fam > test.fam
mv test.fam dummy_case1.fam

# much lower frequencies 
echo "Generating dummy_case2!"
touch dummy_case2.sim
echo "100 null 0.05 0.95 1.00 1.00" >> dummy_case2.sim
echo "10 disease 0.001 0.01 2.00 mult" >> dummy_case2.sim
${XWAS_ORIG} --simulate dummy_case2.sim --simulate-ncases 100 --simulate-ncontrols 1000 --chr 23 --make-bed --out dummy_case2 --silent
generate_covars dummy_case2.fam dummy_case2.covar
recode_as_x dummy_case2.bim > test.bim
mv test.bim dummy_case2.bim
recode_sex dummy_case2.fam > test.fam
mv test.fam dummy_case2.fam

echo "Generating dummy_case3!"
touch dummy_case3.sim
echo "100 control 0.05 0.95 1.00 1.00" >> dummy_case3.sim
echo "10 case 0.05 0.50 2.00 mult" >> dummy_case3.sim
${XWAS_ORIG} --simulate dummy_case3.sim --simulate-ncases 100 --simulate-ncontrols 800 --chr 23 --make-bed --out dummy_case3 --silent
generate_covars dummy_case3.fam dummy_case3.covar
recode_as_x dummy_case3.bim > test.bim 
mv test.bim dummy_case3.bim
recode_sex dummy_case3.fam > test.fam
mv test.fam dummy_case3.fam

echo "Finished Generating Datasets!"
rm -f *.sim *.simfreq *.hh
