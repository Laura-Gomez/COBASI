
////////////////
////////////////// GET SIGNATURE READS
//////////////////
////
//// This script will obtain the sequences and some useful information from the reads that contain either the PreCS or both SignatureCSs.
//// This script will process one FASTA or FASTQ file (with the read sequences) per run.
////
////        Retrieve_SignatureReads.out  SignatureCS ReadType ReadFileX kmer_size > outX
////
////  PARAMETERS.
////	SignatureCS 	A multi-fasta file with the sequences of all Signature CSs
////	ReadType 	FASTA or FASTQ formats are supported [FASTA|FASTQ]
////	ReadFileX	A fasta or fastq file with the reads of the sequencing project
////	kmer_size	CSs size = Jellyfish database kmer size
////	outX 		Output file
////
////OUTPUT.
//// IT WILL DEPEND ON THE KIND OF DISCOVERY BEING MADE
////
//////////////////
////////////////// 
//////////////////
//
//

#include <unordered_map>
#include <fstream>
#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <unistd.h>


#include "read.cc"

using std::string;
using std::ifstream;
using std::cout;
using std::endl;


// REVERSE COMPLEMENT FUNCTION

string reverse_complement(string &seq)
{
	string rc;
	for ( string::iterator it=seq.end(); it != seq.begin()-1; --it)
	{ 
	char letter = *it;
		switch (letter)
		{
		// try to find the complement of every char of the source string
		case 'A':
			rc.append("T"); break;
		case 'T':
			rc.append("A"); break;
              	case 'G':
			rc.append("C");	break;
		case 'C':
			rc.append("G");	break;
              // special characters
              	case 'N':
			rc.append("N");	break;
		default:
			rc.append("N");	break;
		}
	}
	return(rc);                               
}

// SPLIT BY DEMILITER FUNCTION

void split(string &str,char delim, std::vector<string> &splitstr) {
    string buf = "";
    int i = 0;
    while (i < str.length()) {
        if (str[i] != delim)
            buf += str[i];
        else if (buf.length() > 0) {
            splitstr.push_back(buf);
            buf = "";
        }
        i++;
    }
    if (!buf.empty())
        splitstr.push_back(buf);
}

struct RetrReadRow{
	string chr, kmer, id_read, seq_read, qual_read;
	unsigned int pos, st_ev, end_ev, st_ev_parent, end_ev_parent;
	unsigned short pos_read;
	char strand;
};

// SORT RetrReadRow DATA STRUCTURE
// SORT FIRST BY CHR
// ONLY IF POS, ST_EV AND END_EV ARE THE SAME THE RECORD WILL BE SORTED TOGETHER

bool sort_RRR(const RetrReadRow &first,  const RetrReadRow &second)
{
	unsigned int compare_first, compare_second;
	if (first.chr.compare(second.chr) == 0){
		compare_first = first.pos + first.st_ev + first.end_ev;
		compare_second = second.pos + second.st_ev + second.end_ev;
		if (compare_first < compare_second)
		{
			return true;
		}else{
			return false;
		}
	}else{
		if (first.chr.compare(second.chr) < 0 ){
			return true;;
		}else{
			return false;
		}
	}
}

// CONVERSION OF DATA FORMAT
// CONVERSION TO FIT RetrReadRow DATA STRUCTURE
// THE STRING REF HAST MULTIPLE PIECES OF INFORMATION
void add_read_kmer_row(char strand, int pos, string seq, string id, string qual, string kmer, string ref, std::vector<RetrReadRow> &match, string ind){
	RetrReadRow read_kmer_row;
 	read_kmer_row.pos_read = pos;
        read_kmer_row.seq_read = seq;
        read_kmer_row.id_read = id;
        read_kmer_row.qual_read = qual; //SAVING QUALITY VALUES                    
        read_kmer_row.strand = strand;
        read_kmer_row.kmer = kmer;
       	
	string ref_info = ref;
        
	std::vector<string> splitstr;
        split(ref_info, ':', splitstr);
        
	read_kmer_row.chr = splitstr[0];
        read_kmer_row.pos = atoi(splitstr[1].c_str());

	if (ind.compare("CHILD") == 0){
        	read_kmer_row.st_ev = atoi(splitstr[2].c_str());
        	read_kmer_row.end_ev = atoi(splitstr[3].c_str());
        }else{
        	read_kmer_row.st_ev_parent = atoi(splitstr[2].c_str());
        	read_kmer_row.end_ev_parent = atoi(splitstr[3].c_str());	
        	read_kmer_row.st_ev = atoi(splitstr[4].c_str());
        	read_kmer_row.end_ev = atoi(splitstr[5].c_str());
        }

       match.push_back(read_kmer_row);
}

void scan_seq_hash_map (string &seq_read, const string &id_read, const string & qual_read, std::unordered_multimap <string, string> &map, std::unordered_map <string, int> &count, std::vector<RetrReadRow> &match, const int kmer_size, string ind)
{
	// EVERY KMER IN THE READ IS OBTAINED (AND IN THE READ REVERSE COMPLEMENT) AND SEARCHED IN THE HASH TABLE

	int c_plus, c_minus, read_size;
	string plus_str, minus_str, rc_seq_read, ref_info, kmer;

	read_size = seq_read.length(); 

        rc_seq_read = reverse_complement(seq_read);

	// EVERY READ IS SPLITTED IN KMERS
	
	for ( c_plus = 0, c_minus = (read_size - kmer_size + 1); c_plus <= (read_size - kmer_size); c_plus++, c_minus--)
	{
		plus_str = seq_read.substr(c_plus, kmer_size);
		minus_str = rc_seq_read.substr(c_minus, kmer_size);

		std::unordered_map<std::string, std::string>::const_iterator got_plus = map.find(plus_str);
		std::unordered_map<std::string, std::string>::const_iterator got_minus = map.find(minus_str);
		
		// IF THE FORWARD KMER IS FOUND IN THE HASH TABLE, THE LOCALIZATION INFORMATION IS SAVED

		if (got_plus != map.end() ){
		    // IF THE CS BELONGS TO SEVERAL VSR'S
 		    // THE INFORMATION FOR THE ALIGNMENTS IS STORED FOR EVERY VSR (REF_INFO WILL BE DIFFERENT FOR EVERY VSR)
			if (count[plus_str] > 1){
				auto range = map.equal_range(plus_str);
				for (auto it = range.first; it != range.second; ++it){
					kmer = it ->first;
					ref_info = it -> second;
					add_read_kmer_row('+', c_plus+1, seq_read, id_read, qual_read, kmer, ref_info, match, ind);
				}
			}else{
				kmer = got_plus->first;
                                ref_info = got_plus->second;
				add_read_kmer_row('+', c_plus+1, seq_read, id_read, qual_read, kmer, ref_info, match, ind);
			}
		}
		// IF THE REVERSE KMER IS FOUND IN THE HASH TABLE, THE LOCALIZATION INFORMATION IS SAVED
	
		if (got_minus != map.end() ){
			// IF THE CS BELONGS TO SEVERAL VSR'S
		        // THE INFORMATION FOR THE ALIGNMENTS IS STORED FOR EVERY VSR (REF_INFO WILL BE DIFFERENT FOR EVERY VSR)

			if (count[minus_str] > 1){
				auto range = map.equal_range(minus_str);
				for (auto it = range.first; it != range.second; ++it){
						kmer = it -> first;
                                       		ref_info = it -> second;
                                	       	add_read_kmer_row('-', c_plus+1, seq_read, id_read, qual_read, kmer, ref_info, match, ind);
				}
			}else{
			 	kmer = got_minus->first;
                                ref_info = got_minus->second;
                                add_read_kmer_row('-', c_plus+1, seq_read, id_read, qual_read, kmer, ref_info, match, ind);
			}
		}
	}
 }

// MAIN FUNCTION
//
int main(int argc, char *argv[])
{

  int c, kmer_size;
  string cs_file, read_type, read_file_str, ind;

 while ((c = getopt (argc, argv, "hc:t:f:k:i:")) != -1)
  {

 switch(c)
 {
  case 'c':
   cs_file = std::string(optarg);
   break;
  case 't':
   read_type = std::string(optarg);
   break;
  case 'f':
   read_file_str = std::string(optarg);
   break;
  case 'k':
   kmer_size = atoi(optarg);
   break;
  case 'i':
   ind = std::string(optarg);
   break;
  case 'h':
   fprintf (stderr,  "USAGE:\n\tRetrieve_SignatureReads.out  -c SignatureCS -t ReadType -f ReadFileX -k kmer_size -i Ind > outX\n");
   fprintf (stderr, "PARAMETERS.\n\tc [SignatureCS]\t\tA multi-fasta file with the sequences of all Signature CSs");
   fprintf (stderr, "\n\tt [ReadType]\t\tFASTA or FASTQ formats are supported [FASTA|FASTQ]");
   fprintf (stderr, "\n\tf [ReadFileX]\t\tA fasta or fastq file with the reads of the sequencing project");
   fprintf (stderr, "\n\tk [kmer_size]\t\tCSs size = Jellyfish database kmer size");
   fprintf (stderr, "\n\ti [Ind ]\t\tSet to CHILD for SNV discovery in one individual\n\t\t\t\tSet to PARENT for SNV discovery in the parents");
   fprintf (stderr, "\n\toutX\t\t\tOutput file\n");
    return (0);
  }
 }

 std::unordered_multimap <string, string> cs_kmer;
 std::unordered_map <string, int> cs_count; 
 std::vector<std::string> kmer_info;

 // OPEN MULTI-FASTA FILE WITH SIGNATURE CS'S SEQUENCES
 //
 FILE * query_kmer = fopen(cs_file.c_str() , "r");
 FILE * read_file = fopen(read_file_str.c_str(), "r");
 std::vector<RetrReadRow> matching_reads;


// CREATING THE HASH WITH CS ENTRIES
 // CREATING A HASH, EVERY ENTRY IS A CS SEQUENCE 
 // HASH PAIR: CS SEQUENCE -  CS ID (>CHR:CS:PRECS:POSTCS)

 string seq, id, qual;

 while (Fasta_Read(query_kmer, seq, id))
 {
//	cout << id << "\n" << seq << endl;

	for(unsigned int k = 0; k < seq.length(); k++)
	    seq[k] = toupper(seq[k]);
	
	cs_kmer.insert (std::make_pair<string, string>(seq, id)); 

	// THE SAME CS COULD BE INVOLVED IN MORE THAN ONE VSR (SAME SEQUENCE WITH DIFFERENTS ID'S)
	// COUNT THE NUMBER OF VSR'S IN WHICH EVERY CS APPEARS
	// LOOK THE CS IN THE CS_COUNT STRUCTURE BY SEQUENCE

	std::unordered_map<string, int>::const_iterator got_count = cs_count.find(seq);
	if (got_count != cs_count.end()){
		int count = got_count -> second;
		count = count +1;
		cs_count[seq] = count;
	}else{
		cs_count[seq] = 1;
	}			
//	std::pair<string, std::vector<string> > new_kmer (seq, id);
//	cs_kmer.insert(new_kmer);
 }

// SCANNIG THE READ FILE -  EVERY FASTA/FASTQ RECORD IS ANALYZED

 if (read_type == "FASTQ"){
	while (Fastq_Read(read_file, seq, id, qual))
		scan_seq_hash_map(seq, id, qual, cs_kmer, cs_count, matching_reads, kmer_size, ind);	
 }else{
	while(Fasta_Read(read_file, seq, id))
		scan_seq_hash_map(seq, id, "NA", cs_kmer, cs_count, matching_reads, kmer_size, ind);
 }
 

// THE FINAL ENTRIES ARE SORTED ACCORDIN TO THEIR CHROMOSOME AND POSITION ALONG THE REFERENCE GENOME
 std::sort (matching_reads.begin(), matching_reads.end(), sort_RRR);
 
 if (ind.compare("CHILD") == 0){
 	for (unsigned i=0; i<matching_reads.size(); i++)
 	{
//	 fprintf ("%s\t%u\t%s\t%u\t%s\t%s\t%s\t%s\n",  matching_reads.at(i).chr, matching_reads.at(i).pos, matching_reads.at(i).kmer, matching_reads.at(i).pos_read, matching_reads.at(i).strand, matching_reads.at(i).id_read, matching_reads.at(i).seq_read, matching_reads.at(i).qual_read);
		cout << matching_reads.at(i).chr << "\t" << matching_reads.at(i).pos << "\t"  << matching_reads.at(i).st_ev << "\t"  << matching_reads.at(i).end_ev << "\t" << matching_reads.at(i).kmer << "\t" << matching_reads.at(i).pos_read << "\t" << matching_reads.at(i).strand << "\t" << matching_reads.at(i).id_read << "\t" << matching_reads.at(i).seq_read  << "\t" << matching_reads.at(i).qual_read << endl;

 	}
 }else{
 	for (unsigned i=0; i<matching_reads.size(); i++)
        {
	 cout << matching_reads.at(i).chr << "\t" << matching_reads.at(i).pos << "\t"  << matching_reads.at(i).st_ev_parent << "\t"  << matching_reads.at(i).end_ev_parent << "\t" << matching_reads.at(i).st_ev << "\t"  << matching_reads.at(i).end_ev << "\t" << matching_reads.at(i).kmer << "\t" << matching_reads.at(i).pos_read << "\t" << matching_reads.at(i).strand << "\t" << matching_reads.at(i).id_read << "\t" << matching_reads.at(i).seq_read  << "\t" << matching_reads.at(i).qual_read << endl;
	}
 }

return(0);
}  

