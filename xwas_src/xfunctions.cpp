

///////////////////////////////////////////////////////////////////////////
//       Adapted from plink source codes by Purcell et al. 2007          //
//       By Feng Gao and Diana Chang                                     //
///////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>
#include <iterator>
#include <algorithm>
#include <cmath>
#include "options.h"
#include "helper.h"
#include "stats.h"
#include "crandom.h"
#include "linear.h"
#include "logistic.h"
#include "fisher.h"
#include "plink.h"
#include "whap.h"
#include "phase.h"
#include "xoptions.h"
#include "dcdflib.h"
#include "stats.h"

extern ofstream LOG;

#define MISSING(i,l) ( P.SNP[l]->one[i] && ( ! P.SNP[l]->two[i] ) )

void xfilterSNPs(Plink& P)
{ 
	//////////////////////////////////////////////////////
	// This functions applies the following filters and
	// functions:
	// Counts of number of founders, nonfounders
	// Per-individual genotyping rate
	// Read in, or calculate, allele frequencies (and save these)
	// Optionally write allele frequencies, then close
	// Exclude SNPs with too many missing genotypes
	// Identify/correct heterozygote haploid
	// Identify SNPs with no founder genotypes
	// Calculate/report/filter on HWE tests
	// Calculate/report genotyping rate per SNP/per individual
	// Filter on MAF
	// Remove filtered-out SNPs


	bool original_SNP_major = par::SNP_major;

	if ( ! par::SNP_major ) P.Ind2SNP();

	// Which SNPs to delete
	vector<bool> del(P.locus.size(),false);

	// Which individuals to delete
 	vector<bool> indel(P.sample.size(),false);

 	// Male and female controls counts
 	int male_cont = 0, female_cont = 0;

	//////////////////////////////////////////
	// Display number of founders/nonfounders
	P.cnt_f=0;
	vector<Individual*>::iterator person = P.sample.begin();
	while ( person != P.sample.end() )
	{
		if ( (*person)->founder ) P.cnt_f++;
		// Signals they are controls
		if ((*person)->phenotype == 1){
			(*person)->sex ? male_cont++: female_cont++;
		}
		person++;
	}
	P.printLOG(int2str(P.cnt_f)+" founders and "+int2str(P.n-P.cnt_f)+" non-founders found\n");
	if (P.cnt_f<P.n) par::has_nonfounders = true;
 
	////////////////////////////////////////////
	// If we have obligatory missing genotypes:
	// ensure they really are missng
	if ( par::oblig_missing ){
		set<int2>::iterator p = P.oblig_missing.begin();
		while ( p != P.oblig_missing.end() ){
			int l = p->p1;
			int k = p->p2;
			for (int i = 0; i < P.sample.size(); i++){
				Individual * person = P.sample[i];
				if ( person->sol == k ){
					P.SNP[l]->one[i] = true;
					P.SNP[l]->two[i] = false;
				}
			}
			++p;
		}
	}

	/////////////////////////////////////////////////
	// Remove individuals with too many missing calls

	double total_genotyping = 0;

	if ( par::MAX_IND_MISSING < 1 ){
		int n_removed = 0;
		int n_orig = P.n;
		// Consider each individual
		if ( ! par::oblig_missing ){
			for (int i = 0;i < P.sample.size(); i ++){
				bool female = !(P.sample[i]->sex);
				// Sum missingness over all SNPs
				int m=0;       // Missing SNPs
				int nsnps=0;   // All non-obligatory missing SNPs
				for (int l = 0; l < P.locus.size(); l ++){
					// Skip female Y chromosomes
					if ( female && par::chr_Y[P.locus[l]->chr] ) continue;
		  			++nsnps;
					if ( MISSING(i,l) ) m++;
				}
				// Too much missingness?
				if ( (double)m/(double)nsnps > par::MAX_IND_MISSING ){
					indel[i] = true;
					n_removed++;
				}
			} // next individual
		}
		else{ // ... allow oblig missing values
			for (int i=0;i<P.sample.size();i++)
			{
				bool female = ! (P.sample[i]->sex);
				// Sum missingness over all SNPs
				int m=0;       // Missing SNPs
				int nsnps=0;   // All non-obligatory missing SNPs
				for (int l=0; l<P.locus.size();l++){
		  			// Skip female Y chromosomes
					if (female && par::chr_Y[P.locus[l]->chr]) continue;
					if ( ! (P.obligMissing(i,l)) ){
						if (MISSING(i,l)) ++m;
						++nsnps;
					}
				}
				// Too much missingness?
				if ( (double)m/(double)nsnps > par::MAX_IND_MISSING ){
					indel[i] = true;
					n_removed++;
				}
			} // next individual
		} // end if oblig-missing section
		
		////////////////////////////////////////
		// Save list of any removed individuals
		
		if (n_removed>0){
			string f = par::output_file_name + ".irem";
			P.printLOG("Writing list of removed individuals to [ " + f + " ]\n");
			ofstream REM;
			REM.open(f.c_str(), ifstream::out);
			for (int i=0;i<P.sample.size();i++)
				if (indel[i]) REM << P.sample[i]->fid << "\t" << P.sample[i]->iid << "\n";
			REM.close();
			
			// And now remove these individuals, so that
			// SNP-based statistics are calculated with
			// these samples already excluded
			n_removed = P.deleteIndividuals(indel);
		}

		P.printLOG(int2str(n_removed)+" of "+int2str(n_orig));
		P.printLOG(" individuals removed for low genotyping ( MIND > ");
		P.printLOG(dbl2str(par::MAX_IND_MISSING)+" )\n");
	} // end of remove people conditional


	/////////////////////////////////
	// Calculate or read from file?

	if (par::af_read){
		checkFileExists(par::af_file);
		P.printLOG( "Reading allele frequencies from [ " + par::af_file + " ] \n");

		// Make hash of original SNP names
		map<string,int> mlocus;
		map<string,int>::iterator ilocus;

		vector<Locus*>::iterator loc = P.locus.begin();
		int l=0;
		while ( loc != P.locus.end() ){
			mlocus.insert(make_pair( (*loc)->name,l));
			loc++;
			l++;
		}

		// Read allele frequencies
		ifstream FRQ;
		FRQ.open(par::af_file.c_str());
		FRQ.clear();

		string dum1, dum2, dum3, dum4, dum5, dum6;
		string snpname;
		double freq;
		int nm;
		loc = P.locus.begin();
		while ( loc != P.locus.end() ){
			(*loc)->freq = -1;
			(*loc)->nm = 0;
			loc++;
		}
		// Skip header line
		FRQ >> dum1 >> dum2 >> dum3 >> dum4 >> dum5 >> dum6;
		while(!FRQ.eof()){
			vector<string> tokens = tokenizeLine(FRQ);
			if (tokens.size() == 0) continue;
			else if (tokens.size() != 6){
				string sline="";
				for (int i=0; i<tokens.size(); i++)
					sline += tokens[i] + " ";
					error("Problem with allele frequency line: 6 fields required:\n" + sline + "\n");
			}
			ilocus = mlocus.find(tokens[1]);
			if (ilocus != mlocus.end()){
				string rareAllele = tokens[2];
				string commonAllele = tokens[3];
				Locus * loc = P.locus[ilocus->second]; //not very sure
				if( ! from_string<double>( loc->freq, tokens[4],std::dec)){
					loc->freq = 0;
					loc->nm = 0;
				}
				else if( ! from_string<int>(loc->nm,tokens[5],std::dec)){
					loc->freq = 0;
					loc->nm = 0;
				}
				// But was that pointing to the correct allele?
				if ( rareAllele == loc->allele2 && rareAllele != par::missing_genotype && loc->allele2 != par::missing_genotype )
					loc->freq = 1 - loc->freq;
				else if ( commonAllele == loc->allele1 && commonAllele != par::missing_genotype && loc->allele1 != par::missing_genotype )
					loc->freq = 1 - loc->freq;
			}
		}
		FRQ.clear();
		FRQ.close();
	}

	/////////////////////////////////
	// Calculate allele frequencies

	vector<string> hetlist(0);
	vector<bool>::iterator d = del.begin();
	vector<Locus*>::iterator loc = P.locus.begin();
	vector<CSNP*>::iterator s = P.SNP.begin();
	int l = 0; // Main locus counter
	int exc_maf = 0;
	int exc_miss = 0;

	vector<Locus*> no_founders_found_list;
	vector< double > freq_f(P.locus.size(), 0), freq_m(P.locus.size(), 0); //added by Feng Gao, Diana Chang
	vector< int > miss_f(P.locus.size(), 0), miss_m(P.locus.size(), 0); //counter for missingness
	vector< int > nm_f(P.locus.size(), 0), nm_m(P.locus.size(), 0); //added by Feng Gao, Diana Chang, for allele frequencies of chrX for males and females separately.
	vector< string > allele_1(P.locus.size()), allele_2(P.locus.size()); //added by Feng Gao & Diana Chang. Record the two alleles when the freq calculation is done to avoid the swap of alleles.

	vector< string >::iterator a1_it = allele_1.begin(); //added by Feng Gao, Diana Chang
	vector< string >::iterator a2_it = allele_2.begin(); //added by Feng Gao, Diana Chang
	vector< double >::iterator freq_f_it = freq_f.begin(); //added by Feng Gao, Diana Chang
	vector< double >::iterator freq_m_it = freq_m.begin(); //added by Feng Gao, Diana Chang
	vector< int >::iterator miss_f_it = miss_f.begin(); //added by Arjun Biddanda (missingness in females per locus)
	vector< int >::iterator miss_m_it = miss_m.begin(); // added by Arjun Biddanda (missingness in males per locus)

	vector< int >::iterator nm_f_it = nm_f.begin(); //added by Feng Gao, Diana Chang
	vector< int >::iterator nm_m_it = nm_m.begin(); //added by Feng Gao, Diana Chang
	int exc_sexdiff = 0; //added by Feng Gao, Diana Chang
	int exc_missdiff = 0; //added by Arjun Biddanda
	while ( loc != P.locus.end() ){//go over each loci
		if (!par::af_read){
			(*loc)->freq = 0;
			// count 1 per allele, for frequency
			(*loc)->nm = 0;
			*freq_f_it = 0; //added by Feng Gao, Diana Chang
			*freq_m_it = 0; //added by Feng Gao, Diana Chang
			*miss_f_it = 0; //added by Arjun Biddanda
			*miss_m_it = 0; // added by Arjun Biddanda
			*nm_f_it = 0; //added by Feng Gao, Diana Chang
			*nm_m_it = 0; //added by Feng Gao, Diana Chang
		}
		*a1_it = (*loc)->allele1; //added by Feng Gao & Diana Chang
		*a2_it = (*loc)->allele2; //added by Feng Gao & Diana Chang
		// count 1 per genotype, for missingness
		int geno_nm = 0;
		// int male_controls = 0;
		// int female_controls = 0;

		// count 1 per non-obligatory missing genotype
		// (or set to N individuals)
		int geno_real = 0;
		bool X = false;
		bool haploid = false;

		// Determine type of SNP
		if (par::chr_sex[(*loc)->chr]) X=true; //Is this loc a X chr?
		else if (par::chr_haploid[(*loc)->chr]) haploid=true;

		///////////////////////////////
		// Iterate over each individual

		vector<bool>::iterator i1 = (*s)->one.begin();
		vector<bool>::iterator i2 = (*s)->two.begin();
		vector<Individual*>::iterator person = P.sample.begin();
		int i = 0;
		while ( person != P.sample.end() ){	
			bool s1 = *i1;
			bool s2 = *i2;
			// Check female Y genotypes
			if ( par::chr_Y[(*loc)->chr] && ! (*person)->sex ){
				// Set to missing, unless in a RECODE mode
				if ( ! par::preserve_all_genotypes ){
					s1 = *i1 = true;
					s2 = *i2 = false;
				}
				// But in any case, do not include this marker in
				// any genotype counts: skip to next person
				++person;
				++i;
				++i1;
				++i2;
				continue;
			}

			// For haploid heterozygosity check, also consider all individuals
			if ( haploid || ( X && (*person)->sex ) ){
				if ( (!s1) && s2 ){
					hetlist.push_back( (*person)->fid + "\t" + (*person)->iid + "\t" + (*loc)->name );
					// Set to missing, unless in a RECODE mode
					if ( ! par::preserve_all_genotypes ){
						s1 = *i1 = true;
						s2 = *i2 = false;
					}
				}	
			}
		
			// For missing genotypes (Arjun Biddanda)
			if ( !(s1 && (!s2)) ){
			 	geno_nm++;
			 	// Only calculate missingness in controls
			 	if ((*person)->phenotype == 1){
			 		if ((*person)->sex){
			 			*miss_m_it += 1;
			 		}
			 		else{
			 			*miss_f_it += 1;
			 		}
			 	}
			}
			// But is this a real genotype in any case?
			if ( par::oblig_missing ){
				if ( ! (P.obligMissing(i,l)) ) ++geno_real;
			}
			else ++geno_real;
			 



			// Do not recount alleles if we have read in allele frequencies
			if (!par::af_read){
				// For allele frequencies
				// only consider founders
				if ( par::summ_nonfounders || (*person)->founder ){
					// haploid or male(1) on X
					if ( haploid || ( X && (*person)->sex ) ){//chrX for males should be treated differently
						//////////////////
						// Haploid counts
						// "1" allele count

						if ( (!s1) && (!s2) ){   //  FF = hom(11), false = allele1 according to input.cpp
							(*loc)->freq++;
							(*loc)->nm++;
							if (!(*person)->aff) { //added by Feng Gao & Diana Chang 
								*freq_m_it += 1;
								*nm_m_it += 1;
							}
						}
						else if ( s1 && s2 ){   //  TT = hom(22)
							(*loc)->nm++;
							if (!(*person)->aff) { //added by Feng Gao & Diana Chang
								*nm_m_it += 1;
							}
						}
					}
					else{
						//////////////////
						// Autosomal count 
						// "1" allele count

						if (!s1){
							if (!s2){ //   00 = hom(11)
								(*loc)->freq += 2;
								(*loc)->nm += 2;
								if ((*person)->sex) { // added by Feng Gao & Diana Chang, unaffected male
									if (!(*person)->aff) {
										*freq_m_it += 2;
										*nm_m_it += 2;
									}
								} 
								else { // added by Feng Gao & Diana Chang, female
									if (!(*person)->aff) { 
										*freq_f_it += 2;
										*nm_f_it += 2;
									}
								}
							}
							else{ //   01 = het(12)
								(*loc)->freq += 1;
								(*loc)->nm += 2;
								if ((*person)->sex){ // added by Feng Gao & Diana Chang
									if(!(*person)->aff) { 
										*freq_m_it += 1;
										*nm_m_it += 2;
									}
								} 
								else { // added by Feng Gao & Diana Chang
									if(!(*person)->aff) { // unaffected
										*freq_f_it += 1;
										*nm_f_it += 2;
									}
								}
							}
						}
						else if ( s2 ){ // 11 = hom(22)
							(*loc)->nm+=2;
							if ((*person)->sex && !(*person)->aff){ // added by Feng Gao & Diana Chang
								*nm_m_it += 2;
							}
							else if (!(*person)->sex && !(*person)->aff){ // added by Feng Gao & Diana Chang
								*nm_f_it += 2;
							}
						}
					}
				}
			}

			// Print the number of missing controls..
			// Next individual
			++person;
			++i;
			++i1;
			++i2;
		}


		////////////////////////////////
		// Calculate allele frequencies

		if (!par::af_read){
			if ((*loc)->nm>0) (*loc)->freq /= (double)(*loc)->nm;
			else{
				(*loc)->freq = 1;
				// If we aren't getting rid of it anyway
				if ( (double)geno_nm/(double)geno_real >= (1-par::MAX_GENO_MISSING)) no_founders_found_list.push_back(*loc);
			}
		}
		
		//////////////////////////////////////////
		// added by Feng Gao & Diana Chang. Calculating sex stratified X AF
		if (xpar::af_sex) {
			// Females
 			if (*nm_f_it > 0) *freq_f_it /= (double)*nm_f_it;
			else *freq_f_it = 1;
			// Males
			if (*nm_m_it > 0) *freq_m_it /= (double)*nm_m_it;
			else *freq_m_it = 1;
		}

		//////////////////////////////////////////
		// Record total proportion of missingness

		double snp_genotyping = P.n>0 ? (double)geno_nm/(double)geno_real : 0;
		if (snp_genotyping < 1.0){
			// cout<"SOME MISSING!\n";
		}
		total_genotyping += snp_genotyping;

		/////////////////////////////////////////////////
		// Exclude if SNP has too many missing genotypes
		if ( snp_genotyping < (1-par::MAX_GENO_MISSING) ){
			*d = true;
			exc_miss++;
		}

		////////////////////////////////////////////////
		// Make allele1 always the least common allele

		if ( par::make_minor_allele && (!par::af_count) && (*loc)->freq > 0.5 ){

			// then we need to swap alleles
			(*loc)->freq = 1 - (*loc)->freq;			
			string tmp = (*loc)->allele2;
			(*loc)->allele2 = (*loc)->allele1;
			(*loc)->allele1 = tmp;

			vector<bool>::iterator i1 = (*s)->one.begin();
			vector<bool>::iterator i2 = (*s)->two.begin();

			while ( i1 != (*s)->one.end() ){
				if ( (*i1) == (*i2) ){
					*i1 = ! (*i1);
					*i2 = ! (*i2);
				}
				i1++;
				i2++;
			}
		}

		// Next SNP
		++d;
		++loc;
		++l;
		++s;
		++a1_it;
		++a2_it;
		++freq_f_it;
		++freq_m_it;
		++miss_f_it;
		++miss_m_it;
		++nm_f_it;
		++nm_m_it;
	}

	/////////////////////////////////////////////////
	// Save list of any heterozygous haploid alleles

	if (hetlist.size()>0){
		P.printLOG(int2str( hetlist.size()) + " heterozygous haploid genotypes; set to missing\n");
		string f = par::output_file_name + ".hh";
		P.printLOG("Writing list of heterozygous haploid genotypes to [ " + f + " ]\n");
		ofstream REM;
		REM.open(f.c_str(), ifstream::out);
		for (int i=0; i<hetlist.size(); i++) REM << hetlist[i] << "\n";
		REM.close();
	}
	hetlist.clear();

	/////////////////////////////////////////////////
	// Save list of SNPs with no founders observed

	if (no_founders_found_list.size()>0){
		P.printLOG(int2str( no_founders_found_list.size()) + " SNPs with no founder genotypes observed\n");
		P.printLOG("Warning, MAF set to 0 for these SNPs (see --nonfounders)\n");
		string f = par::output_file_name + ".nof";
		P.printLOG( "Writing list of these SNPs to [ " + f + " ]\n");
		ofstream NOF;
		NOF.open(f.c_str(), ifstream::out);
		for (int i=0; i<no_founders_found_list.size(); i++) NOF << no_founders_found_list[i]->name << "\n";
		NOF.close();
	}
	no_founders_found_list.clear();

	//////////////////////////
	// Write allele freq file

	if (par::af_write){
		if (par::include_cluster_from_file) P.calcStratifiedAlleleFreqs();
		else{
			ofstream FRQ;
			string f = par::output_file_name + ".frq";
			if (par::af_count) f += ".count";
			if (par::summ_nonfounders) P.printLOG("Writing allele frequencies (all individuals) to [ " + f + " ] \n");
			else P.printLOG("Writing allele frequencies (founders-only) to [ " + f + " ] \n");
			if (par::af_count) P.printLOG("Display counts rather than frequencies\n");
			FRQ.open(f.c_str(), ifstream::out);
			FRQ.precision(4);
			FRQ << setw(4) << "CHR" << " "
			<< setw(par::pp_maxsnp) << "SNP" << " "
			<< setw(4) << "A1" << " "
			<< setw(4) << "A2" << " ";
			if (par::af_count)
				FRQ << setw(6) << "C1" << " "
				<< setw(6) << "C2" << " "
				<< setw(6) << "G0" << "\n";
			else
				FRQ << setw(12) << "MAF" << " "
				<< setw(8) << "NCHROBS"
				<< "\n";

			vector<Locus*>::iterator loc = P.locus.begin();
			while (loc != P.locus.end() ){
				string a1 = (*loc)->allele1;
				string a2 = (*loc)->allele2;
				if (a1=="") a1="0";
				if (a2=="") a2="0";
				FRQ << setw(4)  << (*loc)->chr  << " "
				<< setw(par::pp_maxsnp) << (*loc)->name  << " "
				<< setw(4)  << a1  << " "
				<< setw(4)  << a2  << " ";

				if (par::af_count){
					FRQ << setw(6) << int( (*loc)->freq ) << " "
					<< setw(6) << int( (*loc)->bp   ) << " "
					<< setw(6) << int( (*loc)->pos  ) << "\n";
				}
				else{
					if ( (*loc)->nm > 0 )	FRQ << setw(12) << (*loc)->freq << " ";
					else	FRQ << setw(12) << "NA" << " ";
					FRQ << setw(8) << (*loc)->nm << "\n";
				}
				loc++;
			}
			FRQ.close();
		}

		// Close after we've done alle freqs,
		shutdown();
	}

	//////////////////////////
	// added by Feng Gao & Diana Chang. Write allele freq file for X alleles

	if (xpar::af_sex && !par::af_read){
		ofstream FRQ;
		string f = par::output_file_name + ".xfrq";
		if (par::summ_nonfounders)	P.printLOG("Writing X sex stratified allele frequencies (all individuals) to [ " + f + " ] \n");
		else	P.printLOG("Writing X sex stratified allele frequencies (founders-only) to [ " + f + " ] \n");
		FRQ.open(f.c_str(), ifstream::out);
		FRQ.precision(4);
		FRQ << setw(4) << "CHR" << " "
		<< setw(par::pp_maxsnp) << "SNP" << " "
		<< setw(4) << "A1_F" << " "
		<< setw(4) << "A2_F" << " "
		<< setw(12) << "MAF_F" << " "
		<< setw(8) << "NCHROBS_F" << " "
		<< setw(4) << "A1_M" << " "
		<< setw(4) << "A2_M" << " "
		<< setw(12) << "MAF_M" << " "
		<< setw(8) << "NCHROBS_M" << "\n";
		vector<Locus*>::iterator loc = P.locus.begin();
		vector< string >::iterator a1_it = allele_1.begin();
		vector< string >::iterator a2_it = allele_2.begin();
		vector< double >::iterator freq_f_it = freq_f.begin(); 
		vector< double >::iterator freq_m_it = freq_m.begin(); 
		vector< int >::iterator nm_f_it = nm_f.begin(); 
		vector< int >::iterator nm_m_it = nm_m.begin(); 
		while (loc != P.locus.end() ){
			if (par::chr_sex[(*loc)->chr]){ //only if this locus is on chrX
				string a1 = *a1_it;
				string a2 = *a2_it;
				if (a1=="") a1="0";
				if (a2=="") a2="0";
				string a1_f = (*freq_f_it <= (1. - *freq_f_it)) ? a1 : a2;
				string a2_f = (*freq_f_it <= (1. - *freq_f_it)) ? a2 : a1;
				string a1_m = (*freq_m_it <= (1. - *freq_m_it)) ? a1 : a2;
				string a2_m = (*freq_m_it <= (1. - *freq_m_it)) ? a2 : a1;
				FRQ << setw(4)  << (*loc)->chr  << " "
				<< setw(par::pp_maxsnp) << (*loc)->name  << " "
				<< setw(4)  << a1_f  << " "
				<< setw(4)  << a2_f << " ";
				if ( *nm_f_it > 0 )	FRQ << setw(12) << min(*freq_f_it, 1. - *freq_f_it) << " ";
				else	FRQ << setw(12) << "NA" << " ";
				FRQ << setw(8) << *nm_f_it << " "
				<< setw(4)  << a1_m  << " "
				<< setw(4)  << a2_m << " ";
				if ( *nm_m_it > 0 )	FRQ << setw(12) << min(*freq_m_it, 1. - *freq_m_it) << " ";
				else	FRQ << setw(12) << "NA" << " ";
				FRQ << setw(8) << *nm_m_it << "\n";
			}
			loc ++;
			a1_it ++;
			a2_it ++;
			nm_f_it ++;
			nm_m_it ++;
			freq_f_it ++;
			freq_m_it ++;
		}
		FRQ.close();
		// Close after we've done alle freqs,
		shutdown();
	}

	//////////////////////////
	// added by Feng Gao & Diana Chang.
	// Check for significant differences between male / female allele counts
	// using fisher's exact test, then write results to a file.
	// This is testing for sig gender differences in controls only for qual traits.
	if (xpar::af_sex_test && !par::af_read){
		if (par::qt) 
			error("Frequency difference test between sexes is for qualitative phenotype only. Quantitative phenotype is detected and command --freqdiff-x is not valid." );
		ofstream FRQ;
		string f = par::output_file_name + ".xfreqtest";
		if (par::summ_nonfounders)
			P.printLOG("Writing X allele frequencies differences (all individuals) to [ " + f + " ] \n");
		else
			P.printLOG("Writing X allele frequencies differences (founders-only) to [ " + f + " ] \n");
		FRQ.open(f.c_str(), ifstream::out);
		FRQ.precision(4);
		FRQ << setw(4) << "CHR" << " "
		<< setw(par::pp_maxsnp) << "SNP" << " "
		<< setw(4) << "A1" << " "
		<< setw(4) << "A2" << " "
		<< setw(8) << "M1" << " "
		<< setw(8) << "M2" << " "
		<< setw(8) << "F1" << " "
		<< setw(8) << "F2" << " "
		<< setw(12) << "Pvalue"
		<< "\n";

		vector<bool>::iterator d = del.begin();
		vector<Locus*>::iterator loc = P.locus.begin();
		vector< string >::iterator a1_it = allele_1.begin();
		vector< string >::iterator a2_it = allele_2.begin();
		vector< double >::iterator freq_f_it = freq_f.begin(); 
		vector< double >::iterator freq_m_it = freq_m.begin(); 
		vector< int >::iterator nm_f_it = nm_f.begin(); 
		vector< int >::iterator nm_m_it = nm_m.begin(); 
		while (loc != P.locus.end() ){
			if (par::chr_sex[(*loc)->chr]){ //only if this locus is on chrX
				string a1 = *a1_it;
				string a2 = *a2_it;
				if (a1=="") a1="0";
				if (a2=="") a2="0";
				FRQ << setw(4)  << (*loc)->chr  << " "
				<< setw(par::pp_maxsnp) << (*loc)->name  << " "
				<< setw(4)  << a1  << " "
				<< setw(4)  << a2  << " "
				<< setw(8)  << *freq_m_it  << " "
				<< setw(8)  << *nm_m_it-*freq_m_it  << " "
				<< setw(8)  << *freq_f_it  << " "
				<< setw(8)  << *nm_f_it-*freq_f_it  << " ";
				table_t t;
				double pvalue;
				sizeTable(t,2,2);
				t[0][0] = *freq_m_it;
				t[0][1] = *nm_m_it - *freq_m_it;
				t[1][0] = *freq_f_it;
				t[1][1] = *nm_f_it - *freq_f_it;
				pvalue = fisher(t);
				if ( pvalue > -1 )	FRQ << setw(12) << pvalue << "\n";
				else	FRQ << setw(12) << "NA" << "\n";
				if (pvalue <= xpar::af_sex_test_limit && pvalue > -1){
					exc_sexdiff++;	
					*d = true;
				}
			}
			loc ++;
			d ++;
			a1_it ++;
			a2_it ++;
			nm_f_it ++;
			nm_m_it ++;
			freq_f_it ++;
			freq_m_it ++;
		}
		FRQ.close();
	}

	///////////////////////////
	// DIFFERENTIAL MISSINGNESS
	// Check for significant differences in genotype missingness rates between 
	// Male and female controls for qualitative traits using Fishers exact test
	if (xpar::am_sex_test){
		// TODO : make this for quantitative traits
		if (par::qt) error("Missingness difference test between sexes is for case/control phenotype only. Quantitative phenotype is detected and command --missdiff-x is not valid." );
		ofstream MIS;
		string m = par::output_file_name + ".xmisstest";
		if (par::summ_nonfounders)
			P.printLOG("Writing X missingness differences (all individuals) to [ " + m + " ] \n");
		else
			P.printLOG("Writing X missingness differences (founders-only) to [ " + m + " ] \n");
		MIS.open(m.c_str(), ifstream::out);
		MIS.precision(4);
		MIS << setw(4) << "CHR" << " "
		<< setw(par::pp_maxsnp) << "SNP" << " "
		<< setw(4) << "A1" << " "
		<< setw(4) << "A2" << " "
		<< setw(8) << "MISS_M" << " "
		<< setw(8) << "TOTAL_M" << " "
		<< setw(8) << "MISS_F" << " "
		<< setw(8) << "TOTAL_F" << " "
		<< setw(12) << "Pvalue"
		<< "\n";

		vector<bool>::iterator d = del.begin();
		vector<Locus*>::iterator loc = P.locus.begin();
		vector< string >::iterator a1_it = allele_1.begin();
		vector< string >::iterator a2_it = allele_2.begin();
		vector< int >::iterator miss_f_it = miss_f.begin(); 
		vector< int >::iterator miss_m_it = miss_m.begin(); 
		vector< int >::iterator nm_f_it = nm_f.begin(); 
		vector< int >::iterator nm_m_it = nm_m.begin(); 
		while (loc != P.locus.end() ){
			if (par::chr_sex[(*loc)->chr]){ //only if this locus is on chrX
				string a1 = *a1_it;
				string a2 = *a2_it;
				if (a1=="") a1="0";
				if (a2=="") a2="0";

				MIS << setw(4)             << (*loc)->chr  << " "
				<< setw(par::pp_maxsnp)    << (*loc)->name  << " "
				<< setw(4)  << a1          << " "
				<< setw(4)  << a2          << " "
				<< setw(8)  << *miss_m_it  << " "
				<< setw(8)  << male_cont   << " "
				<< setw(8)  << *miss_f_it  << " "
				<< setw(8)  << female_cont << " ";
				table_t t;
				double pvalue;
				sizeTable(t,2,2);
				t[0][0] = *miss_m_it;
				t[0][1] = male_cont;
				t[1][0] = *miss_f_it;
				t[1][1] = female_cont;
				pvalue = fisher(t);
				if ( pvalue > -1 ) MIS << setw(12) << pvalue << "\n";
				else	MIS << setw(12) << "NA" << "\n";
				if (pvalue <= xpar::am_sex_test_limit && pvalue > -1){
					exc_sexdiff++;	
					*d = true;
				}
			}
			loc ++;
			d ++;
			a1_it ++;
			a2_it ++;
			nm_f_it ++;
			nm_m_it ++;
			miss_f_it ++;
			miss_m_it ++;
		}

		MIS.close();

		//Close after we have done the differential missingness
		shutdown();
	}

	/////////////////////////
	// Write HWE statistics

	if (par::HWD_test || par::HWD_report){
		ofstream HWD;
		if (par::HWD_report){
			if (par::summ_nonfounders)
				P.printLOG("Writing Hardy-Weinberg tests (all individuals) to [ " + par::output_file_name + ".hwe ] \n");
			else
				P.printLOG("Writing Hardy-Weinberg tests (founders-only) to [ " + par::output_file_name + ".hwe ] \n");
			string f = par::output_file_name + ".hwe";
			HWD.open(f.c_str(), ifstream::out);
			HWD.precision(4);
			HWD << setw(4) << "CHR" << " "
			<< setw(par::pp_maxsnp) << "SNP" << " "
			<< setw(8) << "TEST" << " "
			<< setw(4) << "A1" << " "
			<< setw(4) << "A2" << " "
			<< setw(20) << "GENO" << " "
			<< setw(8) << "O(HET)" << " "
			<< setw(8) << "E(HET)" << " "
			<< setw(12) << "P" << " "
			<< "\n";
		}

		int cnt=0, cnt_a=0, cnt_u=0;
		////////////////////////
		// Consider each locus

		vector<bool>::iterator d = del.begin();
		vector<Locus*>::iterator loc = P.locus.begin();

		for ( int l = 0 ; l < P.locus.size() ; l++ ){
			// Compute p-values for HWE test in cases, controls & all
			// Only consider founders
			int a11, a12, a22;
			int u11, u12, u22;
			int b11, b12, b22;
			a11=a12=a22=0;
			u11=u12=u22=0;
			b11=b12=b22=0;

			bool X = false, haploid = false;
			if (par::chr_sex[(*loc)->chr])	X=true;
			else if (par::chr_haploid[(*loc)->chr])	haploid=true;
			///////////////////////////////////////////////
			// Iterate over each individual, founders only

			for ( int i = 0 ; i < P.sample.size() ; i++ ){
				Individual * person = P.sample[i];
				///////////////////////////////////////////////
				// Only consider founders, & diploid genotypes
				if ( par::summ_nonfounders || person->founder )
					if ( ! ( haploid || ( X && person->sex ) ) ){
						bool s1 = P.SNP[l]->one[i];
						bool s2 = P.SNP[l]->two[i];
						// Consider everybody, irrespective of phenotype
						// (QT, C/C or missing)
						if (!s1){
							if (!s2) b11++;   //   00 = hom(11)
							else b12++;       //   01 = het(12)
						}
						else if ( s2 )	b22++; // 11 = hom(22)
						if (par::bt){  // for binary trait, separately for cases/controls
							if (person->phenotype == 1){
								if (!s1){
									if (!s2) u11++;   //   00 = hom(11)
									else u12++;       //   01 = het(12)
								}
								else if ( s2 ) u22++; //   11 = hom(22)
							}
							else if (person->phenotype == 2){
								if (!s1){
									if (!s2) a11++;   //   00 = hom(11)
									else a12++;         //   01 = het(12)
								}
								else if ( s2 ) a22++; //   11 = hom(22)
							}
						}
					}
			// Next individual
			}

			// Allele frequencies
			double afreq = 0, ufreq = 0, freq = 0;
			bool include_cases = true;
			bool include_controls = true;
			if (par::qt)	freq = ( b11 + (double)b12/2.0 ) / (double)( b11+b12+b22 );
			else{
				afreq = ( a11 + (double)a12/2.0 ) / (double)( a11+a12+a22 );
				ufreq = ( u11 + (double)u12/2.0 ) / (double)( u11+u12+u22 );
				freq =  ( b11 + (double)b12/2.0 ) / (double)( b11+b12+b22 );
				if ( a11+a12+a22 == 0 ) include_cases = false;
				if ( u11+u12+u22 == 0 ) include_controls = false;
			}
			if (par::qt){
				double p;
				if (par::HWD_standard){
				double tot = b11 + b12 + b22;
				double exp_11 = freq * freq * tot;
				double exp_12 = 2 * freq * (1-freq) * tot;
				double exp_22 = (1-freq) * (1-freq) * tot;
				double chisq = ( (b11-exp_11)*(b11-exp_11) ) / exp_11
							+ ( (b12-exp_12)*(b12-exp_12) ) / exp_12
							+ ( (b22-exp_22)*(b22-exp_22) ) / exp_22 ;
				p = chiprobP(chisq,1);
				}
				else	p = SNPHWE( b12, b11, b22 );
				if (par::HWD_report){
					HWD << setw(4) << (*loc)->chr << " "
						<< setw(par::pp_maxsnp) << (*loc)->name << " "
						<< setw(8) << "ALL(QT)" << " "
						<< setw(4) << (*loc)->allele1 << " "
						<< setw(4) << (*loc)->allele2 << " "
						<< setw(20) << (int2str(b11)+"/"+int2str(b12)+"/"+int2str(b22)) << " "
						<< setw(8) << (double)b12/(double)(b11+b12+b22) << " "
						<< setw(8) << 2 * freq * (1-freq)  << " ";
					if ( realnum(p) )	HWD << setw(12) << p << "\n";
					else	HWD << setw(12) << "NA" << "\n";
				}
				if ( p <= par::HWD_limit && p > -1 ){
					cnt++;
					*d = true;
				}
			}
			else{
				// For case/control data
				double p, p_a, p_u;
				if (par::HWD_standard){
					double exp_a11 = afreq * afreq * (a11+a12+a22);
					double exp_a12 = 2 * afreq * (1-afreq) * (a11+a12+a22);
					double exp_a22 = (1-afreq) * (1-afreq) * (a11+a12+a22);
					double exp_u11 = ufreq * ufreq * (u11+u12+u22);
					double exp_u12 = 2 * ufreq * (1-ufreq) * (u11+u12+u22);
					double exp_u22 = (1-ufreq) * (1-ufreq) * (u11+u12+u22);
					double exp_11 = freq * freq * (b11+b12+b22);
					double exp_12 = 2 * freq * (1-freq) * (b11+b12+b22);
					double exp_22 = (1-freq) * (1-freq) * (b11+b12+b22);
					double chisq_a = ( (a11-exp_a11)*(a11-exp_a11) ) / exp_a11
									+ ( (a12-exp_a12)*(a12-exp_a12) ) / exp_a12
									+ ( (a22-exp_a22)*(a22-exp_a22) ) / exp_a22 ;
					double chisq_u = ( (u11-exp_u11)*(u11-exp_u11) ) / exp_u11
									+ ( (u12-exp_u12)*(u12-exp_u12) ) / exp_u12
									+ ( (u22-exp_u22)*(u22-exp_u22) ) / exp_u22 ;
					double chisq = ( (b11-exp_11)*(b11-exp_11) ) / exp_11
									+ ( (b12-exp_12)*(b12-exp_12) ) / exp_12
									+ ( (b22-exp_22)*(b22-exp_22) ) / exp_22 ;
					p = chiprobP(chisq,1);
					p_a = chiprobP(chisq_a,1);
					p_u = chiprobP(chisq_u,1);
				}
				else{
					p = SNPHWE( b12, b11, b22 );
					p_a = SNPHWE( a12, a11, a22 );
					p_u = SNPHWE( u12, u11, u22 );
				}
				if (par::HWD_report){
					HWD << setw(4) << (*loc)->chr << " "
						<< setw(par::pp_maxsnp) << (*loc)->name  << " "
						<< setw(8) << "ALL" << " "
						<< setw(4) << (*loc)->allele1 << " "
						<< setw(4) << (*loc)->allele2 << " "
						<< setw(20)
						<< int2str(b11)+"/"+int2str(b12)+"/"+int2str(b22) << " "
						<< setw(8) << (double)b12/(double)(b11+b12+b22) << " "
						<< setw(8) << 2 * freq * (1-freq)  << " ";
					if ( p > -1 )	HWD << setw(12) << p  << "\n";
					else	HWD << setw(12) << "NA"  << "\n";
					HWD << setw(4) << (*loc)->chr << " "
						<< setw(par::pp_maxsnp) << (*loc)->name  << " "
						<< setw(8) << "AFF" << " "
						<< setw(4) << (*loc)->allele1 << " "
						<< setw(4) << (*loc)->allele2 << " "
						<< setw(20)
						<< int2str(a11)+"/"+int2str(a12)+"/"+int2str(a22) << " "
						<< setw(8) << (double)a12/(double)(a11+a12+a22) << " "
						<< setw(8) << 2 * afreq * (1-afreq)  << " ";
					if (include_cases && p_a > -1 )	HWD << setw(12) << p_a  << "\n";
					else	HWD << setw(12) << "NA" << "\n";
					HWD << setw(4) << (*loc)->chr << " "
						<< setw(par::pp_maxsnp) << (*loc)->name  << " "
						<< setw(8) << "UNAFF" << " "
						<< setw(4) << (*loc)->allele1 << " "
						<< setw(4) << (*loc)->allele2 << " "
						<< setw(20)
						<< int2str(u11)+"/"+int2str(u12)+"/"+int2str(u22) << " "
						<< setw(8) << (double)u12/(double)(u11+u12+u22) << " "
						<< setw(8) << 2 * ufreq * (1-ufreq)  << " ";
					if (include_controls && p_u > -1 )	HWD << setw(12) << p_u  << "\n";
					else	HWD << setw(12) << "NA" << "\n";
				}
				// Increase counts: in cases
				if ( include_cases && p_a < par::HWD_limit && p_a > -1 ) cnt_a++;
				// Controls (and, if possible, exclude on this value)
				if ( include_controls && p_u < par::HWD_limit && p_u > -1 ){
					cnt_u++;
					if ( ! par::HWD_filter_on_all ){
						*d = true;
						cnt++;
					}
				}
				// In total sample, and if needed, exclude here
				if ( p < par::HWD_limit && p>-1 ){
					if ( par::HWD_filter_on_all || ! include_controls ){
						*d = true;
						cnt++;
					}
				}
	    	}
	    	// next locus
			++loc;
			++d;
		}
		// Finish the report...
		if (par::HWD_report)	HWD.close();
		// ...or finish pruning
		P.printLOG( int2str(cnt) + " markers to be excluded based on HWE test ( p <= " + dbl2str(par::HWD_limit) + " )\n");
		if (par::bt){
			P.printLOG("\t" + int2str(cnt_a) + " markers failed HWE test in cases\n");
			P.printLOG("\t" + int2str(cnt_u) + " markers failed HWE test in controls\n");
		}
	}
	
	///////////////////////////////////////////////////
	// Summary statistics for genotyping/missing rates
	if (par::report_missing){
		///////////////////////////////////////////
		// Report by genotyping rate by individual
		// possibly allowing for obligatory missingness
		P.printLOG( "Writing individual missingness information to [ " + par::output_file_name + ".imiss ] \n");
		ofstream MIS;
		string f = par::output_file_name + ".imiss";
		MIS.open(f.c_str(), ifstream::out);
		MIS.precision(4);
		MIS << setw(par::pp_maxfid) << "FID" << " "
			<< setw(par::pp_maxiid) << "IID" << " "
			<< setw(10) << "MISS_PHENO" << " "
			<< setw(8) << "N_MISS" << " ";
		MIS << setw(8) << "N_GENO" << " ";
		MIS << setw(8) << "F_MISS" << "\n";
		for (int i=0; i<P.n; i++){
			MIS << setw(par::pp_maxfid) << P.sample[i]->fid << " "
				<< setw(par::pp_maxiid) << P.sample[i]->iid << " ";
			if (P.sample[i]->missing) MIS << setw(10) << "Y" << " ";
			else MIS << setw(10) << "N" << " " ;
			// Sum missingness over all SNPs
			int m=0;       // Missing SNPs
			int nsnps=0;   // All non-obligatory missing SNPs
			bool female = ! (P.sample[i]->sex);
			if ( ! par::oblig_missing ){
				for (int l=0; l<P.locus.size();l++){
					// Skip female Y chromosomes
					if ( female && par::chr_Y[P.locus[l]->chr] )	continue;
					if ( MISSING(i,l) ) ++m;
					++nsnps;
				}
			}
			else{ // ... allow oblig missing values
				for (int l=0; l<P.locus.size();l++){
					// Skip female Y chromosomes
					if ( female && par::chr_Y[P.locus[l]->chr] )	continue;
					if ( ! (P.obligMissing(i,l)) ){
						if ( MISSING(i,l) )	++m;
						++nsnps;
					}
				}
			}
			MIS << setw(8) << m << " ";
			MIS << setw(8) << nsnps << " ";
			MIS << setw(8) << (double)m/(double)nsnps << "\n";
		}
		MIS.close();

		///////////////////////////////////////////
		// Report by genotyping rate by locus
		// possibly allowing for sample strata
		// possibly allowing for obligatory missingness
		P.printLOG("Writing locus missingness information to [ " + par::output_file_name +".lmiss ] \n");
		f = par::output_file_name + ".lmiss";
		MIS.open(f.c_str(), ifstream::out);
		MIS.clear();
		MIS.precision(4);

		MIS << setw(4) << "CHR" << " " << setw(par::pp_maxsnp) << "SNP" << " ";
		if (par::include_cluster_from_file)	MIS << setw(10) << "CLST" << " ";
		MIS << setw(8) << "N_MISS" << " ";
		MIS << setw(8) << "N_GENO" << " ";
		if (par::include_cluster_from_file)	MIS << setw(8) << "N_CLST" << " ";
		MIS << setw(8) << "F_MISS" << "\n";
		for (int l=0; l<P.locus.size(); l++){
			Locus * loc = P.locus[l];
			bool chrY = par::chr_Y[P.locus[l]->chr];
			// nk==1 for basic missingness (i.e. not stratified by cluster)
			for (int k=0; k<P.nk; k++){
				MIS << setw(4) << loc->chr << " " << setw(par::pp_maxsnp) << loc->name << " ";
				if (par::include_cluster_from_file)	MIS << setw(10) << P.kname[k] << " ";
				int m=0;     // Number of missing genotypes
				int c=0;     // Number of people in cluster
				int nsnps=0; // Number of actual genotypes in cluster
				for ( int i=0; i<P.sample.size(); i++){
					// Skip female Y chromosome calls
					if ( chrY && ! P.sample[i]->sex )	continue;
					if (par::include_cluster_from_file){
						if ( P.sample[i]->sol == k ){
							if ( ( ! par::oblig_missing ) || ( ! (P.obligMissing(i,l)) ) ){
								if ( MISSING(i,l) ) ++m;
								++nsnps;
							}
							++c;
						}
					}
					else{ // ... ignore cluster strata
						if ( ( ! par::oblig_missing ) || ( ! (P.obligMissing(i,l)) ) ){
							if ( MISSING(i,l) ) ++m;
							++nsnps;
						}
					}
				// Next individual
				}
				MIS << setw(8) << m << " ";
				if (par::include_cluster_from_file)	MIS << setw(8) << c << " ";
				MIS << setw(8) << nsnps << " ";
				MIS << setw(8) << (double)m / (double)nsnps << "\n";
			}
		// Next SNP
		}
		MIS.close();
	}

	/////////////////////////////////
	// Remove rare SNPs

	loc = P.locus.begin();
	d = del.begin();
	while ( loc != P.locus.end() ){
		// Note epsilon correction for MAF, due to floating point
		// issues: only apply to the lower MAF range
		if ( (*loc)->freq < 0 || (*loc)->freq + par::epsilon < par::min_af || (*loc)->freq > par::max_af ){
			*d = true;
			exc_maf++;
		}
		d++;
		loc++;
	}

	/////////////////////////////////////////
	// Remove SNPs based on thresholds

	if ( P.locus.size() > 0 ) P.printLOG("Total genotyping rate in remaining individuals is " + dbl2str(total_genotyping/(double)P.locus.size())+"\n");
	P.printLOG(int2str(exc_miss) + " SNPs failed missingness test ( GENO > " + dbl2str(par::MAX_GENO_MISSING)+" )\n");
	P.printLOG(int2str(exc_maf)+" SNPs failed frequency test ( MAF < "+dbl2str(par::min_af));
	if (par::max_af < 0.5 ) P.printLOG(" or MAF > " + dbl2str(par::max_af));
  	P.printLOG(" )\n");
  	if (xpar::af_sex_test) //added by Feng Gao & Diana Chang. Output this sentence only if the test is specified.
		P.printLOG(int2str(exc_sexdiff)+" SNPs failed sex frequency difference test ( p-value < " + dbl2str(xpar::af_sex_test_limit)+" )\n");
  	
	int tmp = P.deleteSNPs(del);
	//////////////////////////////////////////
	// Need to make back to individual major?

	if ( ! original_SNP_major ) P.SNP2Ind();
	return;
}

void xSexRelatedTests(Plink& P) {
	matrix_t resultpv_male, resultpv_female;
	xglmAssocSexSep(P, resultpv_male, resultpv_female);
	if (xpar::strat_sex)	xStratSexAssocTest(P, resultpv_male, resultpv_female);
	if (xpar::sex_diff)	xSexDiffTest(P, resultpv_male, resultpv_female);
	
}
void xStratSexAssocTest(Plink& P, matrix_t& resultpv_male, matrix_t& resultpv_female){
	string f = par::output_file_name;
	if (par::bt){
		f += ".xstrat.logistic";
		P.printLOG("Writing sex-stratified logistic model association results to [ " + f + " ] \n");
	} 
	else{
		f += ".xstrat.linear";
		P.printLOG("Writing sex-stratified linear model association results to [ " + f + " ] \n");
	}
	ofstream ASC;
	ASC.open(f.c_str(),ios::out);
	ASC << setw(4) << "CHR" << " "
	<< setw(par::pp_maxsnp) << "SNP" << " "
	<< setw(10) << "BP" << " "
	<< setw(4) << "A1" << " "
	<< setw(10) << "TEST" << " ";

	if (par::bt && !xpar::xreturn_beta)	ASC << setw(10) << "OR_M" << " ";
	else ASC << setw(10) << "BETA_M" << " ";
	ASC << setw(12) << "P_M" << " ";
	if (par::bt && !xpar::xreturn_beta)	ASC << setw(10) << "OR_F" << " ";
	else ASC << setw(10) << "BETA_F" << " ";
	ASC << setw(12) << "P_F" << " ";
	
	// If-cases for determining headers
	if (xpar::fishers && !xpar::stouffers){
		ASC << setw(12) << "Fisher_Chi_Squared" << " ";
		ASC << setw(12) << "P_comb_Fisher" << " " << "\n";
	}
	if (!xpar::fishers && xpar::stouffers){
		ASC << setw(12) << "Stouffers_Z" << " ";
		ASC << setw(12) << "P_comb_Stouffer" << " " << "\n";
	}
	if (xpar::fishers && xpar::stouffers){
		ASC << setw(12) << "Fisher_Chi_Squared" << " " 
		<< setw(12) << "P_comb_Fisher" << " ";
		ASC << setw(12) << "Stouffers_Z" << " "
		<< setw(12) << "P_comb_Stouffer" << " " << "\n";
	}
	ASC.precision(4);

	// combine the p-values
	vector<Locus*>::iterator loc = P.locus.begin();
	matrix_t::iterator pv_f_it = resultpv_female.begin();
	matrix_t::iterator pv_m_it = resultpv_male.begin();

	//only if there is something for output. If the user specifies --xchr-model 0, recessive, etc., it won't work for chrX.
	if ((resultpv_male.size() > 0) && (resultpv_female.size() > 0)){ 
		while ( loc != P.locus.end() ){
			if (par::chr_sex[(*loc)->chr]){
				ASC << setw(4)  << (*loc)->chr << " "
				<< setw(par::pp_maxsnp) << (*loc)->name << " "
				<< setw(10) << (*loc)->bp << " "
				<< setw(4)  << (*loc)->allele1 << " "
				<< setw(10) << "SexStrat" << " ";
				if (par::bt && !xpar::xreturn_beta){
					if (*(pv_m_it->begin() + 1) > 0)	
						ASC << setw(10) << exp(*(pv_m_it->begin())) << " " 
						<< setw(12) << *(pv_m_it->begin() + 1) << " ";
					else
						ASC << setw(10) << "NA" << " " 
						<< setw(12) << "NA" << " ";
					if (*(pv_f_it->begin() + 1) > 0)
						ASC << setw(10) << exp(*(pv_f_it->begin())) << " " 
						<< setw(12) << *(pv_f_it->begin() + 1) << " ";
					else
						ASC << setw(10) << "NA" << " " 
						<< setw(12) << "NA" << " ";
				}
				else{
					if (*(pv_m_it->begin() + 1) > 0)	
						ASC << setw(10) << *(pv_m_it->begin()) << " " 
						<< setw(12) << *(pv_m_it->begin() + 1) << " ";
					else
						ASC << setw(10) << "NA" << " " 
						<< setw(12) << "NA" << " ";
					if (*(pv_f_it->begin() + 1) > 0)
						ASC << setw(10) << *(pv_f_it->begin()) << " " 
						<< setw(12) << *(pv_f_it->begin() + 1) << " ";
					else
						ASC << setw(10) << "NA" << " " 
						<< setw(12) << "NA" << " ";
				}

				// *(pv_m_it->begin() + 1) is the p-value we are looking for
				double chisq = ((*(pv_f_it->begin() + 1) > 0) && (*(pv_m_it->begin() + 1) > 0)) ? -2*(log(*(pv_f_it->begin() + 1))+log(*(pv_m_it->begin() + 1))) : -1.;

				// We consider weighting by individual, not chromosome here
				int nmales = 0, nfemales = 0;
				for (int i=0; i<P.n; i++) {
					if ( ! P.sample[i]->missing ){
						if ( P.sample[i]->sex )	nmales++;
						else nfemales++;
					}
				}

				// Weighting the females more heavily (if the flag is there)
				if (xpar::sex_weight){
					nfemales = 2 * nfemales;
				}

				// Sample size weighting method (maybe implement inverse variance method sometime?)
				double wm = sqrt(nmales);
				double wf = sqrt(nfemales);	
				double pvm = (*(pv_m_it->begin() + 1));
				double pvf = (*(pv_f_it->begin() + 1));

				// Direction of effect from the odds ratio
				double efm = *(pv_m_it->begin()) > 0 ? 1 : -1;
				double eff = *(pv_f_it->begin()) > 0 ? 1 : -1;

				//Individual test statistics
				double zm = fabs(ltqnorm((pvm / 2.0))) * efm;
				double zf = fabs(ltqnorm((pvf / 2.0))) * eff;

				double stouffer_z = ((wm * zm) + (wf * zf)) / sqrt(nmales + nfemales);
				double pv_stouffer = 2 * (1 - normdist(fabs(stouffer_z)));

				//Weighting by Inverse variance (TODO)





				if (xpar::fishers && !xpar::stouffers){
					if (chisq > 0)	ASC << setw(12) << chisq << " " << setw(12) << chiprobP(chisq,4);
					else	ASC << setw(12) << "NA" << " " << setw(12) << "NA";
				}
				// Here we have set the parameters on this 
				if (!xpar::fishers && xpar::stouffers){
					if ((pvf > 0) && (pvf > 0)) ASC << setw(12) << stouffer_z << " " << setw(12) << pv_stouffer;
			  	else ASC << setw(12) << "NA" << " " << setw(12) << "NA";
				}
				if (xpar::fishers && xpar::stouffers){
					// Need to establish other field for stouffers z-score
					if (chisq > 0) ASC << setw(12) << chisq << " " << setw(12) << chiprobP(chisq, 4);
					else ASC << setw(12) << "NA" << " " << setw(12) << "NA";
					if ((pvf > 0) && (pvf > 0)) ASC << setw(12) << stouffer_z << " " << setw(12) << pv_stouffer;
					else ASC << setw(12) << "NA" << " " << setw(12) << "NA";
				}
				ASC << "\n";
				pv_f_it ++;
				pv_m_it ++;
			}
			loc ++;
		}
	}
}

// ----------- Borrowed functions from PLINK ------

// static inline uint32_t realnum(double dd) {
//   return (dd == dd) && (dd != INFINITY) && (dd != -INFINITY);
// }

double calc_tprob(double tt, double df) {  
  int32_t st = 0;
  int32_t ww = 1;
  double bnd = 1;
  double pp;
  double qq;
  if (!realnum(tt)) {
    return -9;
  }
  tt = fabs(tt);
  cdft(&ww, &pp, &qq, &tt, &df, &st, &bnd);
  if (st != 0) {
    return -9;
  }
  return 2 * qq;
}

double inverse_tprob(double dbl_qq, double df) {
  double qq = dbl_qq * 0.5;
  double pp = 1 - qq;
  int32_t st = 0;
  int32_t ww = 2;
  double bnd = 1;
  double tt;
  cdft(&ww, &pp, &qq, &tt, &df, &st, &bnd);
  if (st != 0) {
    return -9;
  }
  return tt;
}

// ------ End borrowed functions

//Below are helper functions for sex difference test

inline double square(double x) {
	return x * x;
}

int makeRanks(vector< beta >& b) { //b is already sorted by value
	vector< beta >::iterator it = b.begin(), tieSt = it, tmpIt = it;
	double tmpValue = it -> value;
	int nTie = 1, tmpRank = 1, sumRanks = 1;
	it ++;
	while(1) {
		if ((it == b.end()) || (it -> value != tmpValue)) {
			double avgRank = 1. * sumRanks / nTie;
			for(tmpIt = tieSt; tmpIt != it; tmpIt ++) {
				tmpIt -> rank = avgRank;
			}
			if (it == b.end()) {
				break;
			}
			nTie = 1;
			tmpRank ++;
			tmpValue = it -> value;
			sumRanks = tmpRank;
			tieSt = it;
		}
		else { //tie case
			nTie ++;
			tmpRank ++;
			sumRanks += tmpRank;
		}
		it ++;
	}
	return 0; 
}

double corrRank(const vector< beta >& x, const vector< beta >& y) {
	double sumX = 0., sumY = 0., sumX2 = 0., sumY2 = 0., sumXY = 0., n = x.size();
	vector< beta >::const_iterator itX = x.begin(), itY = y.begin();
	for(; itX != x.end(); itX ++, itY ++) {
		sumX += itX -> rank;
		sumY += itY -> rank;
		sumX2 += square(itX -> rank);
		sumY2 += square(itY -> rank);
		sumXY += itX -> rank * itY -> rank;
	}
	return (n * sumXY - sumX * sumY) / sqrt(n * sumX2 - square(sumX)) 
			/ sqrt(n * sumY2 - square(sumY));
}

double spearman(vector< beta >& betaM, vector< beta >& betaF) {
	sort(betaM.begin(), betaM.end(), byValue());
	sort(betaF.begin(), betaF.end(), byValue());
	makeRanks(betaM); makeRanks(betaF);
	sort(betaM.begin(), betaM.end(), byIndex());
	sort(betaF.begin(), betaF.end(), byIndex());
	return corrRank(betaM, betaF);
}

double sexDiffPVal(double bM, double bF, double nM, double nF, double pM, double pF, double corr) { //using # of individuals
	// boost::math::normal ndist(0., 1.);
	// Possibly use ltqnorm?
	// double test_ltq = -ltqnorm(0.5 * pM);
	// double original = quantile(complement(ndist, 0.5 * pM));
	// cout << "Testing ltqnorm:  "<< test_ltq <<" Original:  " << original << "\n";
	// double SEM = bM / quantile(complement(ndist, 0.5 * pM));
	// double SEF = bF / quantile(complement(ndist, 0.5 * pF));
	double SEM = bM / (-ltqnorm(0.5 * pM));
	double SEF = bM / (-ltqnorm(0.5 * pF));
	double SEM2 = square(SEM), SEF2 = square(SEF);
	double SEMF = sqrt(SEM2 + SEF2 - 2. * corr * SEM * SEF);
	double stat = (bM - bF) / SEMF;
	double df = square(SEM2 / nM + SEF2 / nF) / (square(SEM2 / nM) / (nM - 1.) + 
				square(SEF2 / nF) / (nF - 1.));
	// boost::math::students_t tdist(df);
	// double pval = calc_tprob(stat, df);
	// double original_pval = 2 * cdf(complement(tdist, fabs(stat)));
	// cout << "Testing calc_tprob:  "<< test_pval << " Original:  " << original_pval<<"\n";
	// return pval;
	return calc_tprob(stat, df);
}

//sex difference test
void xSexDiffTest(Plink& P, matrix_t& resultpv_male, matrix_t& resultpv_female){
	string f = par::output_file_name;
	if (par::bt){
		f += ".xdiff.logistic";
		P.printLOG("Writing logistic model sex difference test results to [ " + f + " ] \n");
	} 
	else{
		f += ".xdiff.linear";
		P.printLOG("Writing linear model sex difference test results to [ " + f + " ] \n");
	}
	ofstream ASC;
	ASC.open(f.c_str(),ios::out);
	ASC << setw(4) << "CHR" << " "
	<< setw(par::pp_maxsnp) << "SNP" << " "
	<< setw(10) << "BP" << " "
	<< setw(4) << "A1" << " "
	<< setw(10) << "TEST" << " ";

	// We will actually need to start the parsing from here

	if (par::bt && !xpar::xreturn_beta)	ASC << setw(10) << "OR_M" << " ";
	else ASC << setw(10) << "BETA_M" << " ";
	ASC << setw(12) << "P_M" << " ";
	if (par::bt && !xpar::xreturn_beta)	ASC << setw(10) << "OR_F" << " ";
	else ASC << setw(10) << "BETA_F" << " ";
	ASC << setw(12) << "P_F" << " ";
	ASC << setw(12) << "P_DIFF" << "\n";
	
	ASC.precision(4);

	// combine the p-values
	vector<Locus*>::iterator loc = P.locus.begin();
	matrix_t::iterator pv_f_it = resultpv_female.begin();
	matrix_t::iterator pv_m_it = resultpv_male.begin();

	//only if there is something for output. If the user specifies --xchr-model 0, recessive, etc., it won't work for chrX.
	if ((resultpv_male.size() > 0) && (resultpv_female.size() > 0)){ 
		int n = 0;
		vector< beta > betaF, betaM;
		while ( loc != P.locus.end() ){
			if (par::chr_sex[(*loc)->chr]){
				if ((*(pv_m_it->begin() + 1) > 0) && (*(pv_f_it->begin() + 1) > 0)) {
					beta tBM, tBF;
					tBM.index = n; tBF.index = n;
					tBM.value = *(pv_m_it->begin()); tBF.value = *(pv_f_it->begin());
					betaM.push_back(tBM);
					betaF.push_back(tBF);
					n ++;
				} 
				pv_f_it ++;
				pv_m_it ++;
			}
			loc ++;
		}
		double r = spearman(betaM, betaF);
		loc = P.locus.begin();
		pv_f_it = resultpv_female.begin();
		pv_m_it = resultpv_male.begin();
		while ( loc != P.locus.end() ){
			if (par::chr_sex[(*loc)->chr]){
				ASC << setw(4)  << (*loc)->chr << " "
				<< setw(par::pp_maxsnp) << (*loc)->name << " "
				<< setw(10) << (*loc)->bp << " "
				<< setw(4)  << (*loc)->allele1 << " "
				<< setw(10) << "SexDiff" << " ";
				if (par::bt && !xpar::xreturn_beta){
					if (*(pv_m_it->begin() + 1) > 0)	
						ASC << setw(10) << exp(*(pv_m_it->begin())) << " " 
						<< setw(12) << *(pv_m_it->begin() + 1) << " ";
					else
						ASC << setw(10) << "NA" << " " 
						<< setw(12) << "NA" << " ";
					if (*(pv_f_it->begin() + 1) > 0)
						ASC << setw(10) << exp(*(pv_f_it->begin())) << " " 
						<< setw(12) << *(pv_f_it->begin() + 1) << " ";
					else
						ASC << setw(10) << "NA" << " " 
						<< setw(12) << "NA" << " ";
				}
				else{
					if (*(pv_m_it->begin() + 1) > 0)	
						ASC << setw(10) << *(pv_m_it->begin()) << " " 
						<< setw(12) << *(pv_m_it->begin() + 1) << " ";
					else
						ASC << setw(10) << "NA" << " " 
						<< setw(12) << "NA" << " ";
					if (*(pv_f_it->begin() + 1) > 0)
						ASC << setw(10) << *(pv_f_it->begin()) << " " 
						<< setw(12) << *(pv_f_it->begin() + 1) << " ";
					else
						ASC << setw(10) << "NA" << " " 
						<< setw(12) << "NA" << " ";
				}
				
				// We consider weighting by individual, not chromosome here
				int nmales = 0, nfemales = 0;
				for (int i=0; i<P.n; i++) {
					if ( ! P.sample[i]->missing ){
						if ( P.sample[i]->sex )	nmales++;
						else nfemales++;
					}
				}

				//should we initiate some check in here? to make sure they are below 1?
				if ((*(pv_m_it->begin() + 1) > 0) && (*(pv_f_it->begin() + 1) > 0) && (*(pv_m_it->begin() + 1) < 1) && (*(pv_f_it->begin() + 1) < 1)) {
					double bM = *(pv_m_it->begin()), bF = *(pv_f_it->begin());
					double pM = *(pv_m_it->begin() + 1), pF = *(pv_f_it->begin() + 1);
					ASC << setw(12) << sexDiffPVal(bM, bF, nmales, nfemales, pM, pF, r) << " ";
				}
				else {
					ASC << setw(12) << "NA" << " ";
				}
				
				ASC << "\n";
				pv_f_it ++;
				pv_m_it ++;
			}
			loc ++;
		}
	}
}

// Note that if you only have one heterozygous female you 
// cannot do this test
double calcVar(vector_t& z, double mu){
	double n = z.size();
	double var = 0;
	for (vector<double>::iterator it = z.begin(); it != z.end(); ++it){
		var += square((*it - mu));
	}
	var /= (n - 1.);
	return var;
}

double calcMean(vector_t& z){
	double n = z.size();
	double m = 0;
	for (vector<double>::iterator it = z.begin(); it != z.end(); ++it){
		m += (*it);
	}
	m /= n;
	return m;
}

double calcT(double z1, double z02, double s1, double s02, double n1, double n02){
	double t = (abs(z1 - z02)) / sqrt(s1/n1 + s02/n02);
	return t;
}

double calcDF(double s1, double s02, double n1, double n02){
	double num = (s1 / n1) + (s02 / n02);
	num = square(num);
	double denom = square((s1 / n1)) / (n1 - 1.0) + square((s02 / n02)) / (n02 - 1.0);
	double df = num / denom;
	return df;
}

void xVarianceHeterogeneity(Plink& P){
	if (!par::qt) error("--var-het for quantitative traits only");
	string f = par::output_file_name;
	f += ".xvarhet"	;
	P.printLOG("Writing variance heterogeneity test results to [ " + f + " ] \n");
	ofstream ASC;
	ASC.open(f.c_str(),ios::out);
	ASC << setw(4) << "CHR" << " "
	<< setw(par::pp_maxsnp) << "SNP" << " "
	<< setw(10) << "BP" << " "
	<< setw(4) << "A1" << " "
	<< setw(10) << "T-STAT" << " "
	<< setw(10) << "DF" << " "
	<< setw(10) << "PVAL\n";

	vector<Locus *>::iterator loc1 = P.locus.begin();
	boolmatrix_t ones;
	boolmatrix_t twos;
	for (int l=0; l < P.locus.size(); l++){
		if (par::chr_sex[(*loc1)->chr]){
			ones.push_back((P.SNP[l]->one));
			twos.push_back((P.SNP[l]->two));
		}
		++loc1;
	}

	// Get all of the regression residuals for this PLINK fileset
	// TODO : The residuals function changes the lengths of P.SNP[l]->one and P.SNP[l]->two (TODO : Need to resolve this...)
	matrix_t residuals = xglmAssocResidualsVarHet(P); 

	int xsnpindex = 0;
	vector<Locus *>::iterator loc = P.locus.begin();
	// vector<CSNP *>::iterator s = P.SNP.begin();
	for (int l = 0; l < P.locus.size(); l++){
		if (par::chr_sex[(*loc)->chr]){
			// Need to write output here

			ASC << setw(4)  << (*loc)->chr << " "
			<< setw(par::pp_maxsnp) << (*loc)->name << " "
			<< setw(10) << (*loc)->bp << " "
			<< setw(4)  << (*loc)->allele1 << " ";

			int i = 0;
			int j = 0;
			vector_t z1, z02;

			while (j < P.n) {
				Individual * person = P.sample[j];
				if (!(person->missing) && !(person->sex)) { // not male && not missing
					bool s1 = ones[l][j];
					bool s2 = twos[l][j];

					if (residuals[xsnpindex][i] > 0){
						if (!s1){ 
							if (!s2){ //00 == hom(11)
								z02.push_back (residuals[xsnpindex][i]);
							}
							else{ //10 || 01 == het
								z1.push_back (residuals[xsnpindex][i]);
							}
						}	
						else if (s2){ //11 = hom(22)
							z02.push_back (residuals[xsnpindex][i]);
						}
					}
					i++;
				}
				j++;
			}

			//Calculate the actual statistic
			double z1_mean = calcMean(z1);
			double z02_mean = calcMean(z02);
			double s1 = calcVar(z1, z1_mean);
			double s02 = calcVar(z02, z02_mean);
			double n1 = z1.size();
			double n02 = z02.size();
			
			if (n1 > 1 && n02 > 2){
				double t = calcT(z1_mean, z02_mean, s1, s02, n1, n02);
				double df = calcDF(s1, s02, n1, n02);
				double pval = calc_tprob(t, df);
				ASC << setw(10) << t << " "
				<< setw(10) << df << " " 
				<< setw(10) << pval << " ";
				ASC << "\n";
			}else{
				ASC << setw(10) << "N/A" << " "
				<< setw(10) << "N/A" << " " 
				<< setw(10) << "N/A" << " ";
				ASC << "\n";
			}
			xsnpindex++;
		}
		// Next SNP
		++loc;
	}
}

void xglmAssocSexSep(Plink& P, matrix_t& resMale, matrix_t& resFemale){
	if (par::bt)	P.printLOG("Fitting logistic model to males only \n");
	else	P.printLOG("Fitting linear model to males only \n");
//	Plink P_male = P;
//	P_male.filterOnMale();
	vector< bool > missingStatus(P.sample.size());
	vector<Individual*>::iterator person = P.sample.begin();
	vector< bool >::iterator missing_it = missingStatus.begin();
	while(person != P.sample.end()){ //retain a copy for the real missing status
		*missing_it = (*person)->missing;
		person ++;
		missing_it ++;
	}
	person = P.sample.begin();
	while(person != P.sample.end()){ //Set all females to missing
		if ((*person)->sexcode != "1")	(*person)->missing = true;
		person ++;
	}
	resMale = xglmAssoc(P);
//	Plink* PP_male = &P_male;
//	PP_male -> cleanup();
	if (par::bt)	P.printLOG("Fitting logistic model to females only \n");
	else	P.printLOG("Fitting linear model to females only \n");
//	Plink P_female = P;
//	P_female.filterOnFemale();
	person = P.sample.begin();
	missing_it = missingStatus.begin();
	while(person != P.sample.end()){ //Set all females to missing
		if ((*person)->sexcode != "2")	(*person)->missing = true;
		else (*person)->missing = *missing_it;
		person ++;
		missing_it ++;
	}
	resFemale = xglmAssoc(P);
//	Plink* PP_female = &P_female;
//	PP_female -> cleanup();
	person = P.sample.begin();
	missing_it = missingStatus.begin();
	while(person != P.sample.end()){ //Set all females to missing
		(*person)->missing = *missing_it;
		person ++;
		missing_it ++;
	}//recover P
}

matrix_t xglmAssoc(Plink& P){ //basically the same as Plink::glmAssoc(), except that this function returns p-values
	matrix_t pval; //stores both pval and beta
	Plink* Q = &P;
	P.SNP2Ind();
	vector<Locus*>::iterator loc = P.locus.begin();
	/////////////////////////////
	// Determine sex distribution
	int nmales = 0, nfemales = 0;
	for (int i=0; i<P.n; i++)
		if ( ! P.sample[i]->missing ){
			if ( P.sample[i]->sex )	nmales++;
			else	nfemales++;
		}
	bool variationInSex = nmales > 0 && nfemales > 0;
	bool orig_without_main_snp = par::assoc_glm_without_main_snp;
	bool orig_standard_beta = par::standard_beta;
	par::assoc_glm_without_main_snp = false;
	par::standard_beta = false;
	if (par::xchr_model != 0){ //only if chrX is considered
		int l = 0;
		while ( loc != P.locus.end() ){
     		if (par::chr_sex[(*loc)->chr] || par::chr_haploid[(*loc)->chr]){ //only if chrX SNP or haplotype to save running time
				bool X=true;
				bool automaticSex=false;
				Model * lm;
				if (par::bt){
					LogisticModel * m = new LogisticModel(Q);
					lm = m;
				}
				else{
					LinearModel * m = new LinearModel(Q);
					lm = m;
				}
				lm->setMissing();
				if ( par::glm_dominant )	lm->setDominant();
				else if ( par::glm_recessive || par::twoDFmodel_hethom )
					lm->setRecessive();
				string mainEffect = "";
				bool genotypic = false;
				if ( ! par::assoc_glm_without_main_snp ) {//need to assume that this test is for all loci?
					genotypic = par::chr_haploid[(*loc)->chr] ? false : par::twoDFmodel ;
					if ( par::glm_recessive )	mainEffect = "REC";
					else if ( par::glm_dominant )	mainEffect = "DOM";
					else if ( par::twoDFmodel_hethom )	mainEffect = "HOM";
					else	mainEffect = "ADD";	  
					lm->addAdditiveSNP(l); 
					lm->label.push_back(mainEffect);
					if ( genotypic ){
						lm->addDominanceSNP(l);	      
						if ( par::twoDFmodel_hethom )	lm->label.push_back("HET");
						else	lm->label.push_back("DOMDEV");
					}
				}  
				//////////////////////////////////////////////////////////
				// Haplotype test
				if ( par::chap_test ){
					for (int h=1; h < P.whap->current->group.size(); h++){
						lm->addHaplotypeDosage( P.whap->current->group[h] );	    
						lm->label.push_back( "WHAP"+int2str(h+1) );
					}
				}     
				if ( par::proxy_glm ){
					set<int> t1 = P.haplo->makeSetFromMap(P.haplo->testSet);
					lm->addHaplotypeDosage( t1 );
					lm->label.push_back( "PROXY" );	    
				}
				if ( par::test_hap_GLM ){
					set<int>::iterator i = P.haplo->sets.begin();
					while ( i != P.haplo->sets.end() ){
						set<int> t;
						t.insert(*i);
						lm->addHaplotypeDosage( t );
						lm->label.push_back( P.haplo->haplotypeName( *i ) );
						++i;
					}
				}
				//////////////////////////////////////////////////////////
				// Conditioning SNPs?      
				if (par::conditioning_snps){
					if ( par::chap_test ){
						for (int c=0; c<P.conditioner.size(); c++){
							if ( P.whap->current->masked_conditioning_snps[c] ){
								lm->addAdditiveSNP(P.conditioner[c]); 
								lm->label.push_back(P.locus[P.conditioner[c]]->name);
							}
						}
					}
					else{
						for (int c=0; c<P.conditioner.size(); c++){
							lm->addAdditiveSNP(P.conditioner[c]); 
							lm->label.push_back(P.locus[P.conditioner[c]]->name);
						}
					}
				}
				//////////////////////////////////////////////////////////      
				// Sex-covariate (necessary for X chromosome models, unless
				// explicitly told otherwise). For male/female-only tests,
				// clearly there is no variation in sex
				if ( ( par::glm_sex_effect || ( X && !par::glm_no_auto_sex_effect ) ) && variationInSex ){
					automaticSex = true;
					lm->addSexEffect();
					lm->label.push_back("SEX");	  
				}
				//////////////////////////////////////////////////////////
				// Covariates?
				if (par::clist){
					for (int c=0; c<par::clist_number; c++){
						lm->addCovariate(c);
						lm->label.push_back(P.clistname[c]);
					}
				}
				//////////////////////////////////////////////////////////
				// Interactions
				// addInteraction() takes parameter numbers
				// i.e. not covariate codes      
				// 0 intercept
				// 1 {A}
				//   {D}
				//   {conditioning SNPs}
				//   {sex efffect}
				//   {covariates}
				// Allow for interactions between conditioning SNPs, sex, covariates, etc
      	
				////////////////////////////////////////
				// Basic SNP x covariate interaction?       
				// Currently -- do not allow interactions if no main effect 
				// SNP -- i.e. we need a recoding of things here.
				if ( par::simple_interaction && ! par::assoc_glm_without_main_snp ){
	  				// A, D and haplotypes by conditioning SNPs, sex, covariates	  
	  				int cindex = 2;
	  				if ( genotypic )	cindex = 3;
	  				for (int c=0; c<P.conditioner.size(); c++){
						lm->addInteraction(1,cindex);
						lm->label.push_back(mainEffect+"xCSNP"+int2str(c+1));
						if ( genotypic ){
							lm->addInteraction(2,cindex);
							if ( par::twoDFmodel_hethom )	lm->label.push_back("HETxCSNP"+int2str(c+1));	  
							else	lm->label.push_back("DOMDEVxCSNP"+int2str(c+1));	  
						}	      
						cindex++;
					}
					if ( automaticSex ){
						lm->addInteraction(1,cindex);
						lm->label.push_back(mainEffect+"xSEX");	  
						if ( genotypic ){
							lm->addInteraction(2,cindex);
							if ( par::twoDFmodel_hethom )	lm->label.push_back("HETxSEX");
							else	lm->label.push_back("DOMDEVxSEX");
						}	      
						cindex++;
					}
					for (int c=0; c<par::clist_number; c++){
						lm->addInteraction(1,cindex);
						lm->label.push_back(mainEffect+"x"+P.clistname[c]);	  	      
						if ( genotypic ){
							lm->addInteraction(2,cindex);		  
							if ( par::twoDFmodel_hethom )		  
								lm->label.push_back("HETx"+P.clistname[c]); 
							else
								lm->label.push_back("DOMDEVx"+P.clistname[c]); 
						}	      
						cindex++;	      
					}
				}

				//////////////////////////////
				// Fancy X chromosome models      
				if ( X && automaticSex && par::xchr_model > 2 ){ //obviously these fancy models are not used with only females or males
					// Interaction between allelic term and sex (i.e. allow scale of male effect to vary)
					int sindex = 2;
					if ( genotypic )	sindex++;
					sindex += P.conditioner.size(); 
					lm->addInteraction(2,sindex);
					lm->label.push_back("XxSEX");	  	  
					// xchr model 3 : test ADD + XxSEX
					// xchr model 4 : test ADD + DOM + XxSEX
				}
				//////////////////////////////
				// Build design matrix
				lm->buildDesignMatrix();
				//////////////////////////////
				// Clusters specified?
				if ( par::include_cluster ){
					lm->setCluster();
				}
				//////////////////////////////////////////////////
				// Fit linear or logistic model (Newton-Raphson)      
				lm->fitLM();
				lm->validParameters();
				vector_t pval_temp(2, -1.);
				//cout << setprecision(25) << lm->getVar()[0] << " " << setprecision(25) << lm->getVar()[1] << " " << lm->isValid() << " " << realnum(lm->getVar()[1]) << endl; //for test purpose only
				if (lm->isValid() && !(lm->getVar()[1] < 1e-20) && realnum(lm->getVar()[1])){ //the model or result is OK
					pval_temp[0] = lm->getCoefs()[1];
					pval_temp[1] = lm->getPVals()[0];
				}
				pval.push_back(pval_temp);
				/************I don't think I need this any longer***********
				cout << "var: " << setprecision(22) << lm->getVar()[1] << endl;
				cout << (!(lm->isValid())) << " " << (lm->getVar()[1] < 1e-20) << " " << !realnum(lm->getVar()[1]) << endl;
				cout << "pval: " << lm->getPVals()[0] << endl;
				vector< bool > missing = lm->getMissing();
				int nAllele = 0;
				int allele = -9;
				for (int i=0; i < P.n; i++){
					if (missing[i])	continue;
					bool i1 = P.sample[i]->one[l];
					bool i2 = P.sample[i]->two[l];
					if (i1 != allele){	
						allele = i1;
						nAllele ++;
					}
					if (i2 != allele){
						allele = i2;
						nAllele ++;
					}
					if (nAllele >= 2)	break;
				}
				if (nAllele < 2)	pval.push_back(-1.); //if there is only one allele
				else	pval.push_back(lm->getPVals()[0]);
				***********************************************************/
				delete lm; //clear linear model
			}
			loc ++;
			l ++; 
		}
	}
	par::assoc_glm_without_main_snp = orig_without_main_snp;
	par::standard_beta = orig_standard_beta;
	return pval;
}

matrix_t xglmAssocResidualsVarHet(Plink& P){ //basically the same as Plink::glmAssoc(), except that this function returns residuals
	matrix_t residuals; //stores the residuals 
	Plink* Q = &P;
	P.SNP2Ind();
	vector<Locus*>::iterator loc = P.locus.begin();
	
	// vector<Locus *>::iterator localt = P.locus.begin();
	// for (int l=0; l < P.locus.size(); l++){
	// 	if (par::chr_sex[(*localt)->chr]){
	// 		cout << "This is the length of s1 before : " << (P.SNP[l]->one).size() << "\n";
	// 		cout << "This is the length of s2 before : " << (P.SNP[l]->two).size() << "\n";
	// 		cout << "\n";
	// 	}
	// 	++localt;
	// }

	// Make sure to reset missing status for real afterwards..
	vector< bool > missingStatus(P.sample.size());
	vector<Individual*>::iterator person = P.sample.begin();
	vector< bool >::iterator missing_it = missingStatus.begin();
	while(person != P.sample.end()){ //retain a copy for the real missing status
		*missing_it = (*person)->missing;
		person ++;
		missing_it ++;
	}



	// Set all males as missing
	for (int i=0; i<P.n; i++){
		if ( ! P.sample[i]->missing ){
			if ( P.sample[i]->sex )	P.sample[i]->missing = true;
		}
	}

	bool orig_without_main_snp = par::assoc_glm_without_main_snp;
	bool orig_standard_beta = par::standard_beta;
	par::assoc_glm_without_main_snp = false;
	par::standard_beta = false;
	if (par::xchr_model != 0){ //only if chrX is considered
		int l = 0;
		while ( loc != P.locus.end() ){
     		if (par::chr_sex[(*loc)->chr] || par::chr_haploid[(*loc)->chr]){ //only if chrX SNP or haplotype to save running time
				bool X=true;
				bool automaticSex=false;
				Model * lm;
				if (par::bt){
					LogisticModel * m = new LogisticModel(Q);
					lm = m;
				}
				else{
					LinearModel * m = new LinearModel(Q);
					lm = m;
				}
				lm->setMissing();

				if ( par::glm_dominant )	lm->setDominant();
				else if ( par::glm_recessive || par::twoDFmodel_hethom )
					lm->setRecessive();
				string mainEffect = "";
				bool genotypic = false;
				if ( ! par::assoc_glm_without_main_snp ) {//need to assume that this test is for all loci?
					genotypic = par::chr_haploid[(*loc)->chr] ? false : par::twoDFmodel ;
					if ( par::glm_recessive )	mainEffect = "REC";
					else if ( par::glm_dominant )	mainEffect = "DOM";
					else if ( par::twoDFmodel_hethom )	mainEffect = "HOM";
					else	mainEffect = "ADD";	  
					lm->addAdditiveSNP(l); 
					lm->label.push_back(mainEffect);
					if ( genotypic ){
						lm->addDominanceSNP(l);	      
						if ( par::twoDFmodel_hethom )	lm->label.push_back("HET");
						else	lm->label.push_back("DOMDEV");
					}
				}  
				//////////////////////////////////////////////////////////
				// Haplotype test
				if ( par::chap_test ){
					for (int h=1; h < P.whap->current->group.size(); h++){
						lm->addHaplotypeDosage( P.whap->current->group[h] );	    
						lm->label.push_back( "WHAP"+int2str(h+1) );
					}
				}     
				if ( par::proxy_glm ){
					set<int> t1 = P.haplo->makeSetFromMap(P.haplo->testSet);
					lm->addHaplotypeDosage( t1 );
					lm->label.push_back( "PROXY" );	    
				}
				if ( par::test_hap_GLM ){
					set<int>::iterator i = P.haplo->sets.begin();
					while ( i != P.haplo->sets.end() ){
						set<int> t;
						t.insert(*i);
						lm->addHaplotypeDosage( t );
						lm->label.push_back( P.haplo->haplotypeName( *i ) );
						++i;
					}
				}
				//////////////////////////////////////////////////////////
				// Conditioning SNPs?      
				if (par::conditioning_snps){
					if ( par::chap_test ){
						for (int c=0; c<P.conditioner.size(); c++){
							if ( P.whap->current->masked_conditioning_snps[c] ){
								lm->addAdditiveSNP(P.conditioner[c]); 
								lm->label.push_back(P.locus[P.conditioner[c]]->name);
							}
						}
					}
					else{
						for (int c=0; c<P.conditioner.size(); c++){
							lm->addAdditiveSNP(P.conditioner[c]); 
							lm->label.push_back(P.locus[P.conditioner[c]]->name);
						}
					}
				}

				//////////////////////////////////////////////////////////
				// Covariates?
				if (par::clist){
					for (int c=0; c<par::clist_number; c++){
						lm->addCovariate(c);
						lm->label.push_back(P.clistname[c]);
					}
				}
				//////////////////////////////////////////////////////////
				// Interactions
				// addInteraction() takes parameter numbers
				// i.e. not covariate codes      
				// 0 intercept
				// 1 {A}
				//   {D}
				//   {conditioning SNPs}
				//   {sex efffect}
				//   {covariates}
				// Allow for interactions between conditioning SNPs, sex, covariates, etc
      	
				////////////////////////////////////////
				// Basic SNP x covariate interaction?       
				// Currently -- do not allow interactions if no main effect 
				// SNP -- i.e. we need a recoding of things here.
				if ( par::simple_interaction && ! par::assoc_glm_without_main_snp ){
	  				// A, D and haplotypes by conditioning SNPs, sex, covariates	  
	  				int cindex = 2;
	  				if ( genotypic )	cindex = 3;
	  				for (int c=0; c<P.conditioner.size(); c++){
						lm->addInteraction(1,cindex);
						lm->label.push_back(mainEffect+"xCSNP"+int2str(c+1));
						if ( genotypic ){
							lm->addInteraction(2,cindex);
							if ( par::twoDFmodel_hethom )	lm->label.push_back("HETxCSNP"+int2str(c+1));	  
							else	lm->label.push_back("DOMDEVxCSNP"+int2str(c+1));	  
						}	      
						cindex++;
					}
					if ( automaticSex ){
						lm->addInteraction(1,cindex);
						lm->label.push_back(mainEffect+"xSEX");	  
						if ( genotypic ){
							lm->addInteraction(2,cindex);
							if ( par::twoDFmodel_hethom )	lm->label.push_back("HETxSEX");
							else	lm->label.push_back("DOMDEVxSEX");
						}	      
						cindex++;
					}
					for (int c=0; c<par::clist_number; c++){
						lm->addInteraction(1,cindex);
						lm->label.push_back(mainEffect+"x"+P.clistname[c]);	  	      
						if ( genotypic ){
							lm->addInteraction(2,cindex);		  
							if ( par::twoDFmodel_hethom )		  
								lm->label.push_back("HETx"+P.clistname[c]); 
							else
								lm->label.push_back("DOMDEVx"+P.clistname[c]); 
						}	      
						cindex++;	      
					}
				}

				//////////////////////////////
				// Build design matrix
				lm->buildDesignMatrix();
				//////////////////////////////
				// Clusters specified?
				if ( par::include_cluster ){
					lm->setCluster();
				}
				//////////////////////////////////////////////////
				// Fit linear or logistic model (Newton-Raphson)      
				lm->fitLM();
				lm->validParameters();

				vector< vector<double> > Xdes = lm->X;
				int nind = lm->Ysize();
				int np = lm->getNP();

				vector_t resid_temp(nind, -1);
				if (lm->isValid() && !(lm->getVar()[1] < 1e-20) && realnum(lm->getVar()[1])){ //the model or result is OK
					// i = index in nind , j = index in P.n
					int i = 0;
					int j = 0;
					while (i < nind && j < P.n){	
						if (!(P.sample[j]->missing)){
							double e_i = P.sample[j]->pperson->phenotype;
							for (int p=0; p<np; p++){
								e_i -= (lm->getCoefs())[p] * Xdes[i][p];
							}
							e_i = abs(e_i);
							resid_temp[i] = e_i;
							i++;
						}
						j++;
					}
				}

				residuals.push_back(resid_temp);

				delete lm; //clear linear model
			}
			loc ++;
			l ++; 
		}
	}

	person = P.sample.begin();
	missing_it = missingStatus.begin();
	while(person != P.sample.end()){ //Set all females to missing
		(*person)->missing = *missing_it;
		person ++;
		missing_it ++;
	}//recover P

	// vector<Locus *>::iterator loc1 = P.locus.begin();
	// for (int l=0; l < P.locus.size(); l++){
	// 	if (par::chr_sex[(*loc1)->chr]){
	// 		cout << "This is the length of s1 : " << (P.SNP[l]->one).size() << "\n";
	// 		cout << "This is the length of s2 : " << (P.SNP[l]->two).size() << "\n";
	// 		cout << "\n";
	// 	}
	// 	++loc1;
	// }

	par::assoc_glm_without_main_snp = orig_without_main_snp;
	par::standard_beta = orig_standard_beta;
	return residuals;
}

void xAssocFunctions(Plink& P){ //only one function now. May add others in the future.
	if (xpar::var_het) xVarianceHeterogeneity(P);
	if (xpar::strat_sex || xpar::sex_diff)	xSexRelatedTests(P);
}


void calcXepistasis(Plink& P){
	Plink* Q = &P;
	if ( par:: SNP_major) P.SNP2Ind();

	////////////////////////////////////////////
	// Set up results files
	ofstream XEPI;
	string f = par::output_file_name;
	if (par::qt)
		f += ".xepi.qt";
	else {
		f += ".xepi.cc";
	}
	XEPI.open (f.c_str(),ios::out);
	P.printLOG("Writing xepistasis pairwise results to [ " + f + " ] \n");

	XEPI.precision(4);

	XEPI << setw(4) << "CHR1" << " "
	<< setw(par::pp_maxsnp) << "SNP1" << " "
	<< setw(4) << "CHR2" << " "
	<< setw(par::pp_maxsnp) << "SNP2" << " "; 

	if (par::bt)
		XEPI << setw(12) << "OR_INT" << " ";
	else {
		XEPI << setw(12) << "BETA_INT" << " ";
	}

	XEPI << setw(12) << "STAT" << " "
	<< setw(12) << "P" << " "
	<< "\n";

	/////////////////////////////////////////////////////////////////
	// xepi1 and xepi2 thresholds were given in terms of 1(two-sides)
	// calculate appropriate absolute Z score

	P.printLOG("threshold for displaying xepistatic result (--xepi1) : p <= " +dbl2str(xpar::xepi_alpha1)+"\n");
	P.printLOG("threshold for counting epistatic result (--xepi2) : p <= " +dbl2str(xpar::xepi_alpha2)+ "\n");

	xpar::xepi_alpha1 = fabs (ltqnorm(xpar::xepi_alpha1 / 2));
	xpar::xepi_alpha2 = fabs (ltqnorm(xpar::xepi_alpha2 / 2));

	// Regression based test: case/control or quantitative trait
	
	// Take a list of SNPs, or all SNPs (vector<bool> epi1) ??????????????????
	// Test these against either themselves, or all SNPs (vector<bool> epi2)

	// A    x    B
	// ALL  x    ALL   skip e1 > e2
	// SET1 x    ALL  
	// SET1 x    SET1  skip e1 > e2
	// SET1 x    SET2

	bool skip_symm = false;

	// Only output epistatic tests that have p < xpar::xepi_alpha1;
	// Do not even attempt to save any epistatic results

	// Also present summary results for all SNPs
	// (i.e. average / proportion of significant epistatic tests)
	// at a certain alpha level, xpar::xepi_alpha2

	vector<bool> sA(P.nl_all,false);
	vector<bool> sB(P.nl_all,false);

	// Are we using a test set? If so, construct now
	if (par::set_test)
	{
		if (P.snpset.size() > 2)
			error ("Can only specify one or two SETs when testing for xepistasis");
		if (P.snpset.size() == 0)
			error ("There is no valid sets specifed");

		for (int e=0;e<P.snpset[0].size();e++)
			sA [P.snpset[0][e]] = true;

		// Has a second set been specified?

		if (P.snpset.size()==2){
			P.printLOG("SET1 x SET2 xepistasis mode\n");
			for (int e=0;e<P.snpset[1].size();e++)
				sB[P.snpset[1][e]] = true;
		} else if (xpar::Set_by_Set)
		{
			P.printLOG("SET1 x SET1 xepistasis mode\n");
			skip_symm = true;
			for (int e=0;e<P.snpset[0].size();e++)
				sB[P.snpset[0][e]] = true;
		} else // ALL SNPs in the second set
		{
			P.printLOG("SET1 x ALL xepistasis mode\n");
			for (int e=0;e<P.nl_all;e++)
				sB[e] = true;
		}

	}  else 
	{
		P.printLOG ("ALL x ALL xepistasis mode\n");
		skip_symm = true;
		for (int e=0;e<P.nl_all;e++)
		{
			sA[e] = true;
			sB[e] = true; 
		}
	}

	// Use fast aff coding 

	if (par::bt)
		affCoding(P);

	// Count how many items in the SET1
	int epc = 0;
	for (vector<bool>::iterator e1 = sA.begin();
		e1 != sA.end(); e1++)
		if (*e1) epc++;

	int epcc = 0;

	// Keep track of how many epistatic tests actually performed 
	long int nepi = 0;

	vector<int> summary_sig(P.nl_all,0);
	vector<int> summary_good(P.nl_all,0);
	vector<double> best_score(P.nl_all,0);
	vector<int> best_partner(P.nl_all);

	//////////////////////////////////////////
	// Begin interating over pairs : SET x SET

	for (int e1=0;e1<P.nl_all;e1++)
	{
		if (sA[e1])
		{
			if (!par::silent) 
			{
				cout << "Performing tests of xepistasis: group "
				<< ++epcc << " of " << epc << "              \n";
				cout.flush();
			}

			for (int e2=0;e2<P.nl_all;e2++)
			{
				//////////////////////////////////////////
				// Skip this test under certain conditions

				// The SNP is not in the set
				if (!sB[e2]) { cout << "skipping...\n"; continue; }

				// We've already performed this test
				if (e1 >= e2 && skip_symm) continue;

				// Same SNP
				if (e1 == e2) continue;

				//////////////////////////////////
				// Perform test of epistasis here

				Model * lm;

				if (par::bt)
				{
					LogisticModel * m = new LogisticModel(Q);
					lm = m;
				}
				else
				{
					LinearModel * m = new LinearModel(Q);
					lm = m;
				}

				// Set missing data
				lm -> setMissing();

				// Main effect of SNP 1
				lm -> addAdditiveSNP(e1);
				lm -> label.push_back("ADD1");

				// Main effect of SNP 2
				lm -> addAdditiveSNP(e2);
				lm -> label.push_back("ADD2");

				// Epistasis
				lm -> addInteraction(1,2);
				lm -> label.push_back("XEPI");

				// Add Covariates?
				if (par::clist) 
				{
					for (int c=0;c<par::clist_number;c++)
					{
						lm -> addCovariate(c);
						lm -> label.push_back(P.clistname[c]);
					}
				}

				// Build design matrix

				lm -> buildDesignMatrix();

				// Fit linear model
				lm -> fitLM();

				// Did model fit ok?
				lm -> validParameters();

				// Obtain estimates and statistic

				lm -> testParameter = 3; // interaction
				vector_t b = lm -> getCoefs();
				double chisq = lm -> getStatistic();
				double pvalue = chiprobP(chisq,1);
				double z = sqrt(chisq); //??????????????????????

				// Is this result worth displaying?
				if (lm -> isValid())
				{
					// One more valid test performed
					nepi++;

					// Count as a good result

					summary_good[e1]++;
					if (sA[e2]) summary_good[e2]++;

					// Do we want to record this as part of the summary for the first set?
					if ( z >= xpar::xepi_alpha2)
					{
						// first variable will always be in A set
						summary_sig[e1]++;

						// but the second may also be in A set
						if (sA[e2]) summary_sig[e2]++;
					}

					// Is this result the best score yet for marker in set A?

					if ( z >= best_score[e1])
					{
						best_score[e1] = z;
						best_partner[e1] = e2;
					}

					// The second marker might also be in set A

					if (sA[e2])
					{
						if ( z >= best_score[e2])
						{
							best_score[e2] = z;
							best_partner[e2] = e1;
						}
					}

					// Is this result worth displaying?

					if ( z >= xpar::xepi_alpha1)
					{
						XEPI << setw(4) << P.locus[e1]->chr << " "
						<< setw(par::pp_maxsnp) << P.locus[e1]->name << " "
						<< setw(4) << P.locus[e2]->chr << " "
						<< setw(par::pp_maxsnp) << P.locus[e2]->name << " ";

						if (lm -> isValid())
						{
							if (par::bt) 
								XEPI << setw(12) << exp(b[3]) << " "
								<< setw(12) << chisq << " "
								<< setw(12) << pvalue << " "
								<< "\n";
							else
								XEPI << setw(12) << b[3] << " "
								<< setw(12) << chisq << " "
								<< setw(12) << pvalue << " "
								<< "\n" ;
						}
						else 
							XEPI << setw(12) << "NA" << " "
							<< setw(12) << "NA" << " "
							<< setw(12) << "NA" << " "
							<< "\n";

						XEPI.flush();
					}

					// Clean up
					delete lm;

				}


			} // Next pair of SNPs
		}
	}

	if (!par::silent)
		cout << "\n";

	XEPI.close();

	/////////////////////
	// Summary of results

	if (true)
	{
		f += ".summary";
		XEPI.open(f.c_str(), ios::out);
		XEPI.clear();

		P.printLOG("Performed a total of "+int2str(nepi)+" valid SNP x SNP tests\n");

		P.printLOG("Writing epsitasis summary results to [ "+ f + " ]\n");

		XEPI.precision(4);
		XEPI << setw(4) << "CHR" << " "
		<< setw(par::pp_maxsnp) << "SNP" << " "
		<< setw(12) << "N_SIG" << " "
		<< setw(12) << "N_TOT" << " "
		<< setw(12) << "PROP" << " "
		<< setw(12) << "BEST_CHISQ" << " "
		<< setw(12) << "BEST_CHR" << " "
		<< setw(par::pp_maxsnp) << "BEST_SNP" << " "
		<< "\n";

		int c=0;
		for (int e1=0;e1 <P.nl_all;e1++)
		{
			if (sA[e1])
			{
				XEPI << setw(4) << P.locus[e1]->chr << " "
				<< setw(par::pp_maxsnp) << P.locus[e1]->name << " "
				<< setw(12) << summary_sig[e1] << " "
				<< setw(12) << summary_good[e1] << " "
				<< setw(12) << (double)summary_sig[e1] / (double) summary_good[e1] << " "
				<< setw(12) << best_score[e1] * best_score[e1] << " "
				<< setw(4) << P.locus[best_partner[e1]]->chr << " "
				<< setw(par::pp_maxsnp) << P.locus[best_partner[e1]]->name << " "
				<< "\n";
			}
		}

		XEPI.close();

	}






}
	
