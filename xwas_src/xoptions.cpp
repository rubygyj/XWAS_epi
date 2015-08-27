

///////////////////////////////////////////////////////////////////////////
//       Adapted from plink source codes by Purcell et al. 2007          //
//       By Feng Gao and Diana Chang                                     //
///////////////////////////////////////////////////////////////////////////


#include "options.h"
#include "helper.h"
#include "xoptions.h"
#include "plink.h"

bool xpar::xwas = false; //use xwas?
bool xpar::strat_sex = false; //stratify sex?
bool xpar::sex_diff = false; //sex difference test?
bool xpar::fishers = false; // use fishers method?
bool xpar::stouffers = false; // use stouffers method?
bool xpar::sex_weight = false; // weight females as double?
bool xpar::af_sex = false; //output sex-stratified MAF report?
bool xpar::af_sex_test = false; //allele frequency between sexes?
bool xpar::am_sex_test = false; // missingness between sexes?
bool xpar::xreturn_beta = false; //return beta instead of odds ratio for logistic regression
bool xpar::var_het = false; //test for variance of heterogeneity in females?
double xpar::af_sex_test_limit = 0.01; //p-value for af_sex_test
double xpar::am_sex_test_limit = 0.01; //p-value for am_sex_test
bool xpar::xepi = false;// calc SNP x SNP interaction
double xpar::xepi_alpha1 = 1; // p_value for recording interaction pairs
double xpar::xepi_alpha2 = 1; // p_value for recording significant interaction of a SNP in the summary file
bool xpar::Set_by_Set = true;

void xSetOptions(CArgs & a){ //fully parse the command for xwas -- this part is in parse.cpp for previous pLINK code. setOptions() will go first, and at that time all previous pLINK options will have been set.
	if (a.find("--xwas")){
		xpar::xwas = true;
	}
	if (a.find("--freq-x")){
		if (!xpar::xwas) error("Cannot specify XWAS-related command --freq-x without --xwas");
		if (par::af_read) error("Cannot specify --freq-x and --read-freq together");
		xpar::af_sex = true;
	}
	if (a.find("--freqdiff-x")){
		if (!xpar::xwas) error("Cannot specify XWAS-related command --freqdiff without --xwas"); //can --freqdiff-x and --filter_females, --filter_males appear at the same time?
		if (xpar::af_sex) error("Cannot specify --freq-x and --freqdiff-x together");
		if (par::af_read) error("Cannot specify --freqdiff-x and --read-freq together");
		xpar::af_sex_test = true;      		
		xpar::af_sex_test_limit = a.value_double("--freqdiff-x");
      	if (xpar::af_sex_test_limit < 0 || xpar::af_sex_test_limit > 1) error("--freqdiff-x must have a value between 0 and 1");
	}
	if (a.find("--missdiff-x")){
		if (!xpar::xwas) error("Cannot specify XWAS-related command --missdiff without --xwas");
		xpar::am_sex_test = true;
		xpar::am_sex_test_limit = a.value_double("--missdiff-x");
			if (xpar::am_sex_test_limit < 0 || xpar::am_sex_test_limit > 1) error("--missdiff-x must have a value between 0 and 1");
	}
	if (a.find("--strat-sex")){
		if (!xpar::xwas) error("Cannot specify XWAS-related command --strat-sex without --xwas");
		if (par::filter_females) error("Cannot specify --strat-sex and --filter_females together");
		if (par::filter_males) error("Cannot specify --strat-sex and --filter_males together");
		if (a.find("--fishers")){
			xpar::fishers = true;
		}
		if (a.find("--stouffers")){
			xpar::stouffers = true;
		}
		if (!xpar::fishers && !xpar::stouffers) error("Cannot specify sex-stratified test without either --fishers or --stouffers ");
		xpar::strat_sex = true;
	}
	if (a.find("--sex-diff")){ //sex difference test flag
		if (!xpar::xwas) error("Cannot specify XWAS-related command --sex-diff without --xwas");
		if (par::filter_females) error("Cannot specify --sex-diff and --filter_females together");
		if (par::filter_males) error("Cannot specify --sex-diff and --filter_males together");
		xpar::sex_diff = true;
	}
	if (a.find("--xbeta")){
		if (!xpar::xwas) error("Cannot specify XWAS-related command --xbeta without --xwas");
		if (!xpar::strat_sex && !xpar::sex_diff) error("Cannot specify --xbeta without --strat-sex or --sex-diff");
		xpar::xreturn_beta = true;
	}
	if (a.find("--sex-weight")){
		if (!xpar::xwas) error("Cannot specify XWAS-related command --sex-weight without --xwas");
		if (!xpar::strat_sex) error("Cannot specify --sex-weight without --sex-strat");
		if (!xpar::stouffers) error("Cannot specify --sex-weight without --stouffers");
		xpar::sex_weight = true;
	}
	if (a.find("--xepi"))
	{
		if (!xpar::xwas) error("Cannot specify XWAS-related command --xepi without --xwas");
		xpar::xepi = true; 
		cout << "came here";
		if (a.find("--set")) par::set_test = true;
	}
	if (a.find("--set-by-all")) xpar::Set_by_Set = par::drop_sets = false;
	if (a.find("--xepi1")) xpar::xepi_alpha1 = a.value_double("--xepi1");
	if (a.find("--xepi2")) xpar::xepi_alpha2 = a.value_double("--xepi2");
 	if (a.find("--var-het")){
		if (!xpar::xwas) error("Cannot specify XWAS-related command --var-het without --xwas");
		xpar::var_het = true;
	}
	if (xpar::xwas && !xpar::strat_sex && !xpar::af_sex_test && !xpar::am_sex_test && !xpar::af_sex && !xpar::sex_diff && !xpar::var_het && !xpar::xepi){
		error("No command associated with --xwas is specified");
	}
}

// TODO: include miss-diff, and var_het in there as well
void xParseArgs(CArgs& a) {
	if (a.find("--xhelp")){
		cout << "\n"
		<< "Please consult the XWAS manual for more details on functionality.\n"
		<< "\n"
		<< "We have provided a listing of the XWAS options below:\n"
		<< "\n";
		cout 
		<< "--xwas          Specify XWAS functionality to use any of the functions below\n"
		<< "--freq-x        Output allele frequencies for males and females separately\n"
		<< "--freqdiff-x    Differential MAF filtering between males and females\n"
		<< "--missdiff-x    Differential Missingness filtering between males and femalse\n"
		<< "--strat-sex     Sex stratified test of association\n"
		<< "--fishers       Use Fishers method to combine p-values\n"
		<< "--stouffers     Use Stouffers method to combine p-values\n"
		<< "--sex-weight    Weigh females as twice that of males\n"
		<< "--sex-diff      Test the difference in the effect size between males and females\n"
		<< "--var-het       Apply the variance heterogeneity test for females\n"
		<< "--xbeta         Output beta coefficient instead of odds-ratio\n"
		<< "--xepi          Calculate SNP x SNP interaction test\n\n"
		<< "Options from the original PLINK v1.07 that are commonly used in PLINK_XWAS:\n\n"
		<< "--logistic      Perform logistic regression for case-control phenotypes\n"
		<< "--linear        Perform linear regression for quantitative phenotypes\n"
		<< "--xchr-model 2  Code the genotypes of males as 02\n\n";
		cout << "Please see the XWAS manual and manuscript for further details regarding these functions.\n\n";
		shutdown();
	}
}
	
	
