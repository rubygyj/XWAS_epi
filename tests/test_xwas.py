#!/usr/bin/python

import subprocess as sp
import sys

# PLINK XWAS LOCATION
plink_location='../bin/xwas'

# Defining datasets that we have after running setup script
quant_data = ['dummy_quant1', 'dummy_quant2']
case_data = ['dummy_case1', 'dummy_case2', 'dummy_case3']

# Main testing loop for Case/Control Data
print 'Testing Case/Control Datasets\n'
for cd in case_data:
	cmdstr = ""
	cmdstr = plink_location + ' --bfile ' + cd + ' --xwas --freq-x --out ' + cd + ' --silent'
	print cmdstr
	retval = sp.call(cmdstr, shell=True)
	if retval	!= 0:
		print 'Unexpected error in \'freq-x\' option.'
		sys.exit(1)
	cmdstr = plink_location + ' --bfile ' + cd + ' --xwas --freqdiff-x 0.05 --out ' + cd + ' --silent'
	print cmdstr
	retval = sp.call(cmdstr, shell=True)
	if retval	!= 0:
		print 'Unexpected error in \'freqdiff-x\' option.'
		sys.exit(1)
	cmdstr = plink_location + ' --bfile ' + cd + ' --xwas --strat-sex --fishers --out ' + cd + ' --silent'
	print cmdstr
	retval = sp.call(cmdstr, shell=True)
	if retval != 0:
		print 'Unexpected error in \'fishers\' option.'
		sys.exit(1)
	cmdstr = plink_location + ' --bfile ' + cd + ' --xwas --strat-sex --stouffers --out ' + cd + ' --silent'
	print cmdstr
	retval = sp.call(cmdstr, shell=True)
	if retval != 0:
		print 'Unexpected error in \'stouffers\' option.'
		sys.exit(1)
	cmdstr = plink_location + ' --bfile ' + cd + ' --xwas --strat-sex --fishers --stouffers --out ' + cd + ' --silent'
	print cmdstr
	retval = sp.call(cmdstr, shell=True)
	if retval != 0:
		print 'Unexpected error in \'fishers\' and \'stouffers\' option.'
		sys.exit(1)
	# TODO : interesting that the below line needs 'stouffers' 
	cmdstr = plink_location + ' --bfile ' + cd + ' --xwas --sex-weight --strat-sex --stouffers --out ' + cd + ' --silent'
	print cmdstr
	retval = sp.call(cmdstr, shell=True)
	if retval != 0:
		print 'Unexpected error in \'sex-weight\' option.'
		sys.exit(1)
	cmdstr = plink_location + ' --bfile ' + cd + ' --xwas --xbeta --strat-sex --fishers --stouffers --out ' + cd + ' --silent'
	print cmdstr
	retval = sp.call(cmdstr, shell=True)
	if retval != 0:
		print 'Unexpected error in \'xbeta\' or \'strat-sex\' or \'fishers\' or \'stouffers\' option.'
		sys.exit(1) 
	retval = sp.call(plink_location + ' --bfile ' + cd + ' --xwas --var-het --out ' + cd + ' --silent 2> /dev/null', shell=True)
	if retval == 0:
		# It should fail because --var-het should be only applicable to quantitative data
		print 'Unexpected error in \'var-het\' option.'
		sys.exit(1)

# Main Testing Loop for Quantitative Data 
print 'Testing Quantitative Datasets!\n'
for qd in quant_data:
	cmdstr = ''
	cmdstr = plink_location + ' --bfile ' + qd + ' --xwas --freq-x --out ' + qd + ' --silent'
	print cmdstr
	retval = sp.call(cmdstr, shell=True)
	if retval	!= 0:
		print 'Unexpected error in \'freq-x\' option.'
		sys.exit(1)
	cmdstr = plink_location + ' --bfile ' + qd + ' --xwas --strat-sex --fishers --out ' + qd + ' --silent'
	print cmdstr
	retval = sp.call(cmdstr, shell=True)
	if retval != 0:
		print 'Unexpected error in \'fishers\' option.'
		sys.exit(1)
	cmdstr = plink_location + ' --bfile ' + qd + ' --xwas --strat-sex --stouffers --out ' + qd + ' --silent'
	print cmdstr
	retval = sp.call(cmdstr, shell=True)
	if retval != 0:
		print 'Unexpected error in \'stouffers\' option.'
		sys.exit(1)
	cmdstr = plink_location + ' --bfile ' + qd + ' --xwas --strat-sex --fishers --stouffers --out ' + qd + ' --silent'
	print cmdstr
	retval = sp.call(cmdstr, shell=True)
	if retval != 0:
		print 'Unexpected error in \'fishers\' and \'stouffers\' option.'
		sys.exit(1)
	cmdstr = plink_location + ' --bfile ' + qd + ' --xwas --sex-weight --strat-sex --stouffers --out ' + qd + ' --silent'
	print cmdstr
	retval = sp.call(cmdstr, shell=True)
	if retval != 0:
		print 'Unexpected error in \'strat-sex\' option.'
		sys.exit(1)
	cmdstr = plink_location + ' --bfile ' + qd + ' --xwas --sex-diff --out ' + qd + ' --silent'
	print cmdstr
	retval = sp.call(cmdstr, shell=True)
	if retval != 0:
		print 'Unexpected error in \'sex-diff\' option.'
		sys.exit(1)
	cmdstr = plink_location + ' --bfile ' + qd + ' --xwas --var-het --out ' + qd + ' --silent'
	print cmdstr
	retval = sp.call(cmdstr, shell=True)
	if retval != 0:
		print 'Unexpected error in \'var-het\' option.'
		sys.exit(1)
	cmdstr = plink_location + ' --bfile ' + qd + ' --xwas --var-het --covar ' + (qd+'.covar') +' --out ' + qd + ' --silent'
	print cmdstr
	retval = sp.call(cmdstr, shell=True)
	if retval != 0:
		print 'Unexpected error in \'var-het\' option regarding covariates.'
		sys.exit(1)
	cmdstr = plink_location + ' --bfile ' + qd + ' --xwas --xbeta --strat-sex --fishers --stouffers --out ' + qd + ' --silent'
	print cmdstr
	retval = sp.call(cmdstr, shell=True)
	if retval != 0:
		print 'Unexpected error in \'xbeta\' or \'strat-sex\' or \'fishers\' or \'stouffers\' option.'
		sys.exit(1)
	retval = sp.call(plink_location + ' --bfile ' + qd + ' --xwas --freqdiff-x 0.05 --out ' + qd + ' --silent 2> /dev/null', shell=True) 
	if retval == 0:
		# This should not work for quantitative datasets
		print 'Unexpected error in \'freqdiff-x\' option.'
		sys.exit(1)


print 'XWAS testing completed! '
