////////////////////
////////////   COMPUTE SNV GENOTYPE
////////////////////
//
//
// 
// In this script, the likelihood for each genotype is computed and the most probable genotype
// (given the observed data) is reported as the final genotype.
//  The probability for four possible genotypes is computed:
// Heterozygous reference (R/NR), Homozygous reference (R/R),
// Heterozygous non-reference (NR1/NR2), Homozygous non-reference (NR/NR).
// For The child, this script does not print the homozygous reference-sites.
//
//
//              COMMAND LINE.
// Compute_Genotype.out out.perbase FilePerbaseType IndividualType undetermined.out > genotype.out
//
//
// PARAMETERS.
//	out.perbase 		File with single nucleotide polymorphism data
//	FilePerbaseType 	P if total number of alignments per allele is in column 9 
//	IndividualType	For individual SNV discovery set to CHILD
//	undetermined.out	Output file with undetermined genotypes
//	genotype.out 		File with genotype information per SNV 
//
// OUTPUT.
// A file with assigned genotypes per basepair per chromosome with 11 columns:
//	1-9 will contain the information from the file out.perbase
//	A numeric classificarion for the most probable genotype:
//		1 - Heterozygous reference (Ref/NoRef)
//		2 - Homozygous reference (Ref/Ref)
//		3 - Heterozygous non-reference (NoRef1/NoRef2)
//		5 - Homozygous non-reference (NoRef1/NoRef1) 
//	The different alleles for the assigned genotype
//
//
////////////////////
////////////////////
////////////////////


#include <stdlib.h> 
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm> 
#include "util.cc"

using std::string;
using std::endl;
using std::cout;

struct EvCov{
        std::string nt;
        int cov;
        int plus_total, plus_partial, minus_total, minus_partial;
};

struct by_cov{
        bool operator() (EvCov const &a, EvCov const &b){
                return a.cov > b.cov;
        }
};


// OBTAIN THE TWO HIGHEST VALUES
// g = 0; HIGHEST PROBABILITY HOMOZYGOUS NON-REFERENCE
// g = 1; HIGHEST PROBABILITY HETEROZYGOUS
// g = 2; HIGHEST PROBABILITY HOMOZYGOUS REFERENCE
void two_max( long double a, long double b, long double c , long double  &max, long double &sec_max, unsigned int &g)
{
        if (a > b){
                if (a > c){
                        max = a;
                        g = 0;
                        sec_max = ( b < c ) ? c : b;
                }else{
                        max = c;
                        g = 2;
                        sec_max = a;
                }
        }else{
                if (b > c){
                        max = b;
                        g = 1;
                        sec_max = ( a < c) ? c : a;
                }else{
                        max = c;
                        g = 2;
                        sec_max = a;
                }
        }
}





void gen_max( long double g1, long double g2, long double g4, long double g5, long double  &max, long double &sec_max)
{
	std::vector<long double> g_vector;
	long double gs[] = {g1, g2, g4, g5};

	g_vector.assign(gs, gs+4);
	std::sort(g_vector.begin(), g_vector.end());

	max = g_vector[g_vector.size() -1];
	sec_max = g_vector[g_vector.size() -2];

}

// COMPUTE GENOTYPE LIKELIHOOD
void compute_likelihood(string ref, std::vector<string> nt_split, std::vector<string> ratio_split, long double &g0, long double &g1, long double &g2, double pl, unsigned int &g)
{
	float error = 0.01;
	int ref_reads = 0;
	int noref_reads = 0;
	int total_reads = 0;
	int max_noref = 0;
	int hete_reads = 0;
	int cov;

	long double g0_ref, g0_noref, g2_ref, g2_noref, max, sec_max;
	
	// ALL READS ARE DIVIDED INTO TWO CATEGORIES REFERENCE READS AND NON-REFERENCE READS
	// GET NUMBER OF READS FOR REFERENCE AND NON-REFERENCE ALLELES
	// IT CAN BE MORE THAN ONE NON-REFERENCE ALLELE
	for (unsigned int i = 0; i < nt_split.size(); i++){
		cov = atoi(ratio_split[i].c_str());
		if (nt_split[i].compare(ref) == 0 ){
			ref_reads = cov;
		}else{
			if (cov > max_noref){
				max_noref = cov;
			}
			noref_reads += atoi(ratio_split[i].c_str());  
		}
	} 
	
	total_reads = ref_reads + noref_reads;
	
	// PROBABILITY OF GENOTYPE BEING HOMOZYGOUS NON-REFERENCE
	g0_ref = pow(error, ref_reads);
	g0_noref = pow(1-error, noref_reads);
	g0 = g0_ref*g0_noref;
	
	// PROBABILITY OF GENOTYPE BEING HETEROZYGOUS
	hete_reads = ref_reads + max_noref;
	g1 = pow(pl, hete_reads) * pow(error, (total_reads - hete_reads)) ;
	
	// PROBABILITY OF GENOTYPE BEING HOMOZYGOUS REFERENCE
	g2_ref = pow(1-error, ref_reads);
	g2_noref = pow(error, noref_reads);
	g2 = g2_ref * g2_noref;

	two_max(g0, g1, g2, max, sec_max, g);
}

//STRUCTRES NEED TO DO THE ALLELE ASSIGNATION


// TYPE_FILE = 'T' ALL ALN INFO IS IN COLUMN 6
// TYPE_FILE = 'P' ALL ALN INFO IS IN COLUMN 8
int main(int argc, char *argv[])
{

  int c;
  string ind, undetermined;
  char *var_file, *type_file;

 while ((c = getopt (argc, argv, "hp:t:i:u:")) != -1)
  {

 switch(c)
 {
  case 'p':
   var_file = optarg;
   break;
  case 't':
   type_file = optarg;
   break;
  case 'i':
   ind = std::string(optarg);
   break;
  case 'u':
   undetermined = std::string(optarg);
   break;
  case 'h':
     fprintf (stderr,  "USAGE:\n\tCompute_Genotype.out -p out.perbase -t FilePerbaseType -i IndividualType -u undetermined.out > genotype.out");
     fprintf (stderr, "\nPARAMETERS.\n\tp [out.perbase]\t\t\tFile with single nucleotide polymorphism data");
     fprintf (stderr, "\n\tt [FilePerbaseType]\t\tP if total number of alignments per allele is in column 9");
     fprintf (stderr, "\n\ti [IndividualType]\t\tFor individual SNV discovery set to CHILD\n\t\t\t\t\tSet to PARENT for SNV discovery in the parents");
     fprintf (stderr, "\n\tu [undetermined.out]\t\tOutput file with undetermined genotypes");
     fprintf (stderr, "\n\tgenotype.out\t\t\tFile with genotype information per SNV\n\n");
    return (0);
  }
 }

   
	std::ofstream un_file;
	string line, ref, nt, new_ref, alleles;
	std::vector<string> var_info;
	long double g0, g1, g2, g3, g4, g5, g6, g7, g8, lh, lh_noref;
	unsigned int g, g_temp, max_i, coverage, al_size;
	long double  max_lh, sec_max_lh;
    std::ifstream var;
	std::vector<string> ratio_split, ratio_hete, nt_split, nt_hete;
	EvCov coverage_record;
    std::vector<EvCov> coverage_struct;

	int ploidy = 2;
	string ratio = "";
	unsigned int max=0;
	unsigned int homo_ref=0;



       un_file.open(undetermined);
	

    var.open(var_file);
	double pl = double(1)/double(ploidy);

	for (line; std::getline(var, line);){
		g0 = 0; g1 = 0; g2 = 0; g3 = 0; g4 = 0; g5 =0; max = 0;
 		var_info.clear();
		split(line, '\t', var_info);
		ref = var_info[4];
		nt = var_info[5];
		if (type_file[0] == 'T')
			ratio = var_info[6];
		else
			ratio = var_info[8];
		
		
		ratio_split.clear();
		nt_split.clear();
		split(ratio, '/', ratio_split);
		split(nt, '/', nt_split);

		compute_likelihood(ref, nt_split, ratio_split, g0, g1, g2, pl, g);
		// MOST PROBABLE GENOTYPE HOMOZYGOUS NON-REFERENCE
		
		if (g == 0){
			// CONGRUENCY TEST
			// THE NEW REFERENCE IS THE MOST FREQUENT ALLELE
			// g = 0 (g3) = HomoNoRef (the most common is not common enough ... non informative site)
			// g = 1 (g4) = Heterozygous  (NR1/NR2)
			// g = 2 (g5) = HomoRef (NR/NR))

			for (unsigned int i = 0; i < nt_split.size(); i++){
				coverage = atoi(ratio_split[i].c_str());
					if (max < coverage ){
					max = coverage;
					max_i = i;
				}
			}
			
			new_ref = nt_split[max_i];
			compute_likelihood(new_ref, nt_split, ratio_split, g3, g4, g5, pl, g);
			// IF MOST MOST PROBABLE GENOTYPE IS AGAIN HOMOZYGOUS NON-REFERENCE
			// ASSIGN GENOTYPE TO NA =1000
			g = g + 3;				
			if (g == 3){
				g = 1000;
			}
		}	

		// IF THE GENOTYPE IS HETEROZYGOUS, THE SECOND ALLELE SHOULD BE WELL DEFINED
		// ALLELES WITH LESS THAN TWO READS ARE DISCARDED
		// ELIMINATE REFERENCE 
		// ASSIGN AS REFERENCE THE MOST COMMON ALLELE, IF IT IS NOT COMMON ENOUGH TO BE HOMOZYGOUS REFERENCE ASSIGN GENOTYPE TO NA=1000
		if ( g == 1 || g == 4){
			max = 0; max_i = 0;
			nt_hete.clear(); ratio_hete.clear();
			for (unsigned int i = 0; i < nt_split.size(); i++){
				if (nt_split.at(i).compare(ref) == 0){
					continue;
				}
				if (atoi(ratio_split[i].c_str()) <= 2){
					continue;
				}
				nt_hete.push_back(nt_split[i]);
				ratio_hete.push_back(ratio_split[i]);
                if (max < atoi(ratio_split[i].c_str()) ){
					max = atoi(ratio_split[i].c_str());
                    max_i = i;
                }
            }
			new_ref = nt_split[max_i];
			compute_likelihood(new_ref, nt_hete, ratio_hete, g6, g7, g8, pl, g_temp);
			if (g_temp != 2){
				g = 1000;
			}
		}
	
		//DOING ALLELE ASIGNATION
		alleles = "";
        coverage_struct.clear();
		
        for (unsigned int i = 0; i < nt_split.size(); i++){
            coverage_record.nt = nt_split[i];
            coverage_record.cov = atoi(ratio_split[i].c_str());
            coverage_struct.push_back(coverage_record);
    	}
		
		//SORT COVERAGE STRUCTURE BY COVERAGE, THE MOST ABUNDANT ALLELES ARE THE ONES INVOLVED IN THE GENOTYPE
        std::sort(coverage_struct.begin(), coverage_struct.end(), by_cov());

		// NUMBER OF ALLELES PER GENOTYPE
        
          switch(g){
            case 2:
            	al_size = 1;
                break;
			case 5:
				al_size = 1;
                break;
            case 1:
                al_size = 2;
                break;
			case 4:
				al_size = 2;
                break;
            default :
                al_size = 0;
                alleles = "NA";
         }

		// SOME EVENTS WITH ONLY ONE KIND OF READS WILL DE CLASSIFIED AS HETEROZYGOTES BECAUSE OF THE LOW COVERAGE
		if (al_size > nt_split.size()){
			alleles = "NA-GEN";
			g = 1000;
		}else{
		// CREATE ALLELES STRUCTURE WITH THE ALLELES INVOLVED IN THE GENOTYPE
			for (unsigned int i = 0; i < al_size; i++){
                alleles += coverage_struct.at(i).nt;
                if (i != (al_size -1)){
                    alleles += "/";
                }
            }
		}

	
		// IF THE GENOTYPE CAN NOT BE DETERMINED ( g = 1000)
		// FOR THE CHILD DOES NOT PRINT THE HOMOZYGOUS REFERENCE POSITIONS
		if (g >= 1000){
			un_file << line << "\t" << "NA"  << "\t" << alleles << endl;
		}else if (ind.compare("CHILD") == 0 and g == 2){
			homo_ref = homo_ref +1;
			continue;	
		}else{
			gen_max(g1, g2, g4, g5, max_lh, sec_max_lh);
			lh = max_lh/sec_max_lh;
			cout << line << "\t" << g  << "\t" << alleles << endl;
		}
	}
//	std::cout << "HOMO-REF "<< homo_ref << std::endl;
}
