////////////////////
////////////   ALIGN READ-REFERENCE GENOME
////////////////////
//
//
// In this step every read is cut based on the positions recorded in the READALN files.
// This partial sequence read is aligned to the corresponding region in the RG.
// Incongruent read are filterred out, PCR duplicates are eliminated, variable regions are obtained.
// Alleles per base (reference  or alternative) are obtained, no gaps are allowed.
// FOR MORE INFORMATION SEE USE MANUAL
//
// COMMAND LINE
// 				Align_read_RG –a READALN-file –r REFDIR/ -s sufixREF –k k –t ReadRegion 
//				–p ReadPerbase –l LengthPartial –c ChildSNV –d RemoveDupFlag –i Ind –o out.total
//								 –q out.perbase
//
// ONE INDIVIDUAL SNV DISCOVERY PARAMETERS
//
//	a [READALN-FILE] 	File per chromosome that contains information about the 
//			alignment between the read and the respective Signature CS’s
//	r [REFDIR] 		RG directory path
//	s [sufixREF] 		Sufix of the fasta reference files
//	k [k] 			CSs size = Jellyfish database kmer size
//	t [ReadRegion] 	Each mismatch region between the reads and the RG must 
//					be supported by at least ReadRegion different reads
//	p [ReadPerbase] 	Each mismatch nucleotide between the reads and the RG 
//						must be supported by at least ReadPerbase different reads
//	l [LengthPartial]	Length of extensión of partial  alignment (set to 10)
//	c [ChildSNV]		Set to NA (Parameter used in family-framework)
//	d [RemoveDupFlag]	Set to TRUE to remove PCR duplicates (FALSE otherwise)
//	i [Ind]			For one individual SNV discovery set to CHILD
//	o [out.total]		Output file with regions of polymorphism data
//	q [out.perbase] 	Output file with single nucleotide polymorphism data
//
// FAMILY-BASED SNV DISCOVERY PARAMETERS
//
//	a [READALN-FILE ]	File per chromosome that contains information about the 
//			alignment between the read and the respective Signature CS’s
//	r [REFDIR]		RG directory path
//	a [sufixREF] 		Sufix of the fasta reference files
//	k [k ]			CSs size = Jellyfish database kmer size
//	t [ReadRegion]	Parameter used in One Individual SNV discovery(set to 1)
//	p [ReadPerbase ]	Each mismatch nucleotide between the reads and the RG 
//						must be supported by at least ReadPerbase different reads
//	l [LengthPartial]	Length of extensión of partial  alignment (set to 10)
//	c [ChildSNV]		Path to variable regions file for child individual (out.total)
//	d [RemoveDupFlag]	Set to TRUE to remove PCR duplicates (FALSE otherwise)
//	i [Ind]			Set to PARENT for SNV discovery in the parents
//	o [out.total] 		Set to NA for SNV discovery in the parents
//	q [out.perbase]	Output file with single nucleotide polymorphism data
//
//
// FOR MORE INFORMATION ABOUT OUTPUT FILES SEE USER MANUAL
//
////////////////////
////////////////////
////////////////////




#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "align.hh"
#include "util.cc"


using std::string;
using std::ifstream;


struct RetrReadRow{
	string chr, kmer, id_read, seq_read, qual_read;
	unsigned int pos, end_ev, st_ev, st_ev_child, end_ev_child, read_length;
	unsigned short pos_read;
	string strand;
};

struct AlnSum{
	string strand, event_ref, event_read, chr, type;
	long long unsigned int st_ref, end_ref;
	unsigned int  st_read, end_read, st_signature, end_signature, distance_ref, distance_read;
};

struct FinalEventsRegion {
	unsigned int st_ref, end_ref, count_total, count_partial, count;
	string event_ref, event_qry;
};

struct InfoAln{
        string read_id, type, ref_aln, read_aln;
        int st_read, end_read;
	unsigned int read_length;
};

struct AlnSumBase{
        string chr, type, strand;
        string event_ref, event_read;
        unsigned long long pos_ref;
        unsigned int  pos_read;
};

struct FinalEventsPerbase {
        unsigned long long ref;
        unsigned int count_total, count_partial, count;
        string event_ref, event_qry;
};

	
void read_subset(const string &file, std::vector<unsigned long long> &subset, std::vector<unsigned long long> &subset_ev_pos)
{
        std::vector<string> linesplit;
        string line;
        std::ifstream subset_file;
        unsigned long long st, ev_st, ev_end;
        subset_file.open(file.c_str());

	for (line; std::getline(subset_file, line);)
        {
            linesplit.clear();
       		split(line, '\t', linesplit);
       		st = atoi(linesplit[1].c_str());
	        subset.push_back(st);
            ev_st = atoi(linesplit[3].c_str());
            ev_end = atoi(linesplit[4].c_str());
       		for (unsigned long long pos = ev_st; pos <= ev_end; pos++){
                    subset_ev_pos.push_back(pos);
            }
        }
        std::sort(subset_ev_pos.begin(), subset_ev_pos.end());
}

void read_RRR_parent(const string &file, std::vector<RetrReadRow> &rrr_registry, std::vector<unsigned long long> &subset)
{
        std::vector<string> linesplit;
        string line;
        std::ifstream rrr_file;
        std::vector<unsigned long long>::iterator got_pos;
        unsigned int entry = 0;
        int st_event_child;
        rrr_file.open(file.c_str());

        for (line; std::getline(rrr_file, line);)
        {
                linesplit.clear();
                RetrReadRow record;
                split(line, '\t', linesplit);
                record.chr = linesplit[0];
                record.pos = atoi(linesplit[1].c_str());
                record.st_ev = atoi(linesplit[2].c_str());

                // THE INFORMATION TO COMPARE THE POSITIONS WUTH THE CHILD EVENTS SHOULD BE THE CHILD INFORMATION (MOD PARENT)
                st_event_child = atoi(linesplit[4].c_str());
                if (subset[entry] > st_event_child){    continue; }
                while (subset[entry] < st_event_child && entry < subset.size()) {
                        entry++;
                }
                if (subset[entry] == st_event_child){
                        // THE COLUMNS IN WHICH ALL THIS INFORMATION IS STORED ARE SPECIFIC FOR THE PARENTS FILES
                        record.end_ev = atoi(linesplit[3].c_str());
                        record.st_ev_child = st_event_child;
						record.end_ev_child = atoi(linesplit[5].c_str());
						record.kmer = linesplit[6];
                        record.pos_read = atoi(linesplit[7].c_str());
                        record.strand = linesplit[8];
                        record.id_read = linesplit[9];
                        record.seq_read = linesplit[10];
                        rrr_registry.push_back(record);
                }

        }
}

void read_RRR(const string &file, std::vector<RetrReadRow> &rrr_registry)
{
	std::vector<string> linesplit;
	string line;
	std::ifstream rrr_file;
	rrr_file.open(file.c_str());
	
	for (line; std::getline(rrr_file, line);)
	{
		
		linesplit.clear();
		RetrReadRow record;
		split(line, '\t', linesplit);

		record.chr = linesplit[0];
		record.pos = atoi(linesplit[1].c_str());
		record.st_ev = atoi(linesplit[2].c_str());
		record.end_ev = atoi(linesplit[3].c_str());
		record.kmer = linesplit[4];
		record.pos_read = atoi(linesplit[5].c_str());
		record.strand = linesplit[6];
		record.id_read = linesplit[7];
		record.seq_read = linesplit[8];	
		record.read_length = record.seq_read.length();
		rrr_registry.push_back(record);		
	}
}


void consensus_event_perbase(std::vector<AlnSumBase> &aln_summary_perbase, unsigned short thr_read, const int &end_end, const int &partial, RetrReadRow &rrr, const int kmer, std::ofstream &file_perbase)
{
	string event_data, all_data, final_event, partial_counts_string, total_counts_string, final_counts_string, ref_data, ref;
	unsigned long long pos;
    std::map <unsigned long long, std::vector<string> > nt_changes;
	std::map <unsigned long long, string> event_vector;
	std::map <string, unsigned long long> final_counts;
	std::vector <unsigned long long> event_pos;
    std::map <string, unsigned long long> final_total_counts, final_partial_counts;
    std::map<string, unsigned long long>::iterator got_count, got_total, got_partial;
    std::vector<string>::iterator it_nt;
    std::vector<unsigned long long>::iterator it;
    std::map<unsigned long long, string>::iterator it_event;

//	std::cout << "CONSENSUS-PERBASE " << std::endl;

	// COUNTING CASES PER EVENT
	 for (unsigned i = 0; i < aln_summary_perbase.size(); i++)
    {
        AlnSumBase aln_row = aln_summary_perbase.at(i);
        ref_data = std::to_string(aln_row.pos_ref);
        event_data =  aln_row.event_ref + ":" +  aln_row.event_read;
        all_data = ref_data + ":" + event_data;

        // STORING THE EVENTS WITH A NUCLEOTIDE IN THE REFERENCE
        // SAVE ALL POSSIBLE EVENTS IN THE READS FOR EACH REFERENCE POSITION
        // SAVE REFERENCE NUCLEOTIDE
        if (aln_row.event_ref.compare("-") != 0 ){
            it_nt = std::find(nt_changes[aln_row.pos_ref].begin(), nt_changes[aln_row.pos_ref].end(), aln_row.event_read);
            if  ( it_nt == nt_changes[aln_row.pos_ref].end() ){
                nt_changes[aln_row.pos_ref].push_back(aln_row.event_read);
            }
            event_vector[aln_row.pos_ref] = aln_row.event_ref;
        }

        // CREATING A VECTOR WITH ALL REFERENCE POSITION FOR VARIABLE NUCLEOTIDES
        it = std::find(event_pos.begin(), event_pos.end(), aln_row.pos_ref);
        if (it == event_pos.end()){
            event_pos.push_back(aln_row.pos_ref);
        }

        // COUNTING THE NUMBER OF READS THAT SUPPORT EVERY REFERENCE-READ EVENT
        got_count = final_counts.find(all_data);
        if (got_count != final_counts.end())
        {
            final_counts[all_data]++;
        }else{
            final_counts[all_data] = 1;
        }

		// COUNTING THE NUMBER OF READS WITH BOTH SIGNATURE CS'S THAT SUPPORT EVERY REFERENCE-READ EVENT
		if (aln_row.type.compare("TOTAL") == 0){
            got_total = final_total_counts.find(all_data);
            if (got_total != final_total_counts.end())
            {
                final_total_counts[all_data]++;
            }else{
                final_total_counts[all_data] = 1;
            }
        }else{ 
        	// COUNTING THE NUMBER OF READS WITH ONLY PRE-CS THAT SUPPORT EVERY REFERENCE-READ EVENT
			got_partial = final_partial_counts.find(all_data);
            if (got_partial != final_partial_counts.end())
            {
                final_partial_counts[all_data]++;
            }else{
                final_partial_counts[all_data] = 1;
            }
        }
	}
	
	// FOR EACH STORED REFERENCE POSITION
	for (unsigned int index = 0; index < event_pos.size(); index++)
    {
        pos = event_pos[index];
        it_event = event_vector.find(pos);
 
        if (it_event != event_vector.end()){
            ref = it_event -> second;
			final_event = "";
			partial_counts_string = "";
            total_counts_string = "";
            final_counts_string = "";
            unsigned long long partial_temp = 0;
			unsigned int gap_flag = 0;

			// FOR EACH REFERENCE POSITION
			// FOR EACH ALLELE AT THAT PARTICULAR POSITION
			// CONCATENATE INFORMATION FROM TOTAL-ALN, PARTIAL-ALN AND ALL-ALN
			for (unsigned j = 0; j < nt_changes[pos].size(); j++){
                string key = std::to_string(pos) + ":" + ref + ":" + nt_changes[pos].at(j);
                // ONLY CONCATENATE INFORMATION FOR EVENTS WITH MORE THAN THR_READ ALL-ALN
                // OTHER INFORMATION FOR OTHER EVENTS AT THAT REFERENCE POSITION COULD BE CONCATENATED IF IT PASSES THR_READ
				if (final_counts[ key ] >= thr_read){
                	if (final_event.compare("") != 0 ){
                    	final_event += "/";
                        partial_counts_string += "/";
                        total_counts_string += "/";
                        final_counts_string += "/";
                    }
                    final_event += nt_changes[pos].at(j);
                    partial_temp += final_partial_counts[key];
                    partial_counts_string += std::to_string(final_partial_counts[key]);
                    total_counts_string += std::to_string(final_total_counts[key]);
                    final_counts_string += std::to_string(final_counts[key]);
                    
                    // FILTER OUT THE EVENTS IN WHICH THERE IS A GAP AT THE QUERY GENOME
					if (nt_changes[pos].at(j).compare("-") == 0){
						gap_flag = 1;
					}
                }
            }
            // WRITE FINAL EVENTS
            if (final_event.compare("") != 0 && gap_flag == 0){
				file_perbase << rrr.chr << "\t" << rrr.st_ev << "\t" << rrr.end_ev << "\t" << pos << "\t" << ref << "\t" << final_event << "\t" << total_counts_string << "/" << end_end << "\t" << partial_counts_string << "/" << partial_temp  << "\t" << final_counts_string << "/" << partial_temp + end_end << std::endl;
 //           std::cout << rrr.chr << "\t" << rrr.st_ev << "\t" << rrr.end_ev << "\t" << pos << "\t" << ref << "\t" << final_event << "\t" << total_counts_string << "/" << end_end << "\t" << partial_counts_string << "/" << partial_temp  << "\t" << final_counts_string << "/" << partial_temp + end_end << std::endl;

		}
    	}
    }
}

// STORES ALL NO GAP EVENTS
void event_validation_perbase (const string &ref_aln, const string &read_aln, const string &chr, int pos_ref, int pos_read, int st_ref, int end_ref, string strand, int j, string &event_ref, string &event_read, std::vector<AlnSumBase> &aln_summary, const string &type, unsigned int belongs)
{
    AlnSumBase record;
	if (event_ref.compare("-") != 0 && event_read.compare("-") != 0){
        record.pos_ref = pos_ref;
        record.pos_read = pos_read;
        record.event_ref = event_ref;
        record.event_read = event_read;
        record.strand = strand;
        record.chr = chr;
        record.type = type;
        // FOR THE CHILD, BELOGS IS SET TO 1
        // FOR THE PARENTS, BELONGS=1, MEANS THAT THE RECORD BELONGS TO A VARIABLE POSITION IN THE CHILD
		if (belongs == 1){ 
            aln_summary.push_back(record);
//		std::cout << event_ref << "\t" << event_read << "\t" << type << std::endl;
        }
    }
}


// DO PARTIAL ALIGNMENTS
// STORES INFORMATION FOR NT IN VARIABLE REGION FROM TOTAL ALN IN THE CHILD
// STORES INFORMATION FOR NT IN VARIABLE REGION FROM THE CHILD IN THE PARENTS
// IT DOES NOT STORE GAPS
void aln_perbase( std::map<string, string> &reads_id, std::map<string, unsigned int> &st_cs, std::map<string, unsigned int> &end_cs, std::map<string, string> &strand_list, RetrReadRow rrr, string &chr_fa, string &chr_id, std::vector<AlnSumBase> &aln_summary_perbase, unsigned int &end_end, unsigned int &partial, int kmer, std::vector<unsigned long long> &subset_ev_pos, std::vector<InfoAln> &info_reads_aln, unsigned int max_ref, unsigned int max_read, string flag_ind){
        string ref_seq, read_id, read_seq, ref_seq_up, ref_aln, read_aln, type, cons_aln, rc,event_ref, event_read, strand;
        int st_ref, pos_read, st_read, end_read;
        unsigned int length, end_ref, pos_ref, st_event_ref, st_event_read, gap_read_end, gap_ref_end, belongs_child;
        int ENDFREE = 1;
        int VERBOSE = 0;
        int event = 0;
        std::map<string, string>::iterator it, got_read;
        std::vector<unsigned long long>::iterator it_aln;

		// FOR EVERY READ
		// IF THE READ CONTAIN BOTH SIGNATURE CS'S, LOAD ALN INFO
		// OTHERWISE, DO THE ALN
        for (unsigned int i=0; i < info_reads_aln.size(); i++){
			read_id = info_reads_aln.at(i).read_id;
			strand = strand_list[read_id];
            st_ref = rrr.st_ev;
			cons_aln.clear();
//			std::cout << read_id << std::endl;

			if(info_reads_aln.at(i).type.compare("TOTAL") == 0){
				// LOAD ALN INFO
            	ref_aln = info_reads_aln.at(i).ref_aln;
                read_aln = info_reads_aln.at(i).read_aln;
                end_end++;
            }
            else{
            	// DO THE ALN
				got_read = reads_id.find(read_id);
				read_seq = got_read->second;
                        
				// ESTABLISH NEW START OR END POSITIONS, BASED ON OPTIMAL DISTANCE CALCULATED FROM TOTAL ALN
				if( strand.compare("+") == 0){
                    st_read = info_reads_aln.at(i).st_read;
                    end_read = st_read + max_read -1;
				}else{
                    end_read = info_reads_aln.at(i).end_read;
                    st_read = end_read - max_read + 1;
				}
//				std::cout << "ALN-PB " << max_read << std::endl;
//				std::cout << st_read << "\t" << end_read << std::endl;
				// FILTER AMBIGUITIES
				if( st_read < 1 || end_read > info_reads_aln.at(i).read_length){ continue; }

				// UPPER AND FORWARD SEQUENCE
//				std::cout << st_ref <<"\t" << max_ref << "\t" << st_read << "\t" << max_read << std::endl;
				ref_seq = chr_fa.substr(st_ref - 1, max_ref + 1);
                read_seq = read_seq.substr(st_read - 1, max_read);

                ref_seq_up.clear();
                ref_seq_up = to_upper_seq(ref_seq);

                if (strand.compare("-") == 0)
                {
                    rc = reverse_complement(read_seq);
					read_seq = rc;
                }

                ref_aln.clear();
                read_aln.clear();
                global_align_aff(ref_seq_up, read_seq, ref_aln, read_aln, ENDFREE, VERBOSE);
                partial++;
            }
//        			std::cout << "MAX " << max_ref << "\t" << max_read << std::endl; 
//			std::cout << ref_aln << "\t" << read_aln << std::endl;
			// ANALYZE ALN INFORMATION BASE BY BASE
            
	if (strand.compare("+") == 0){
				pos_read = info_reads_aln.at(i).st_read;
            }else{
                pos_read = info_reads_aln.at(i).end_read;
            }
            pos_ref = st_ref;

            // IN THE CHILD ALL NT FOUND IN VARIABLE REGIONS RETRIEVED FROM TOTAL ALN ARE ANALYZED 
            // INFORMATION FOR ALN WITH NO VARIATION IS ALSO ANALYZED  
//			if (flag_ind.compare("CHILD") == 0){
				for (unsigned int il= 0; il < ref_aln.length(); il++)
	            {
	            	// IF THERE IS NO VARIATION
	            	// BUT THE POSITION IS FOUND IN A VARIABLE REGION (FROM TOTAL ALN), THE INFO IS STORED
  		            if (ref_aln[il] == read_aln[il]){
                		cons_aln = cons_aln + " ";
                        it_aln = std::find(subset_ev_pos.begin(), subset_ev_pos.end(), pos_ref);
                        if (it_aln != subset_ev_pos.end()){
                        	event_ref = ref_aln[il];
   	                     	event_read = read_aln[il];
	                        event_validation_perbase(ref_aln, read_aln, chr_id, pos_ref, pos_read, st_ref, end_ref, strand, il, event_ref, event_read, aln_summary_perbase, info_reads_aln.at(i).type, 1);
        	            }
                	}else{
                    	cons_aln = cons_aln + "*";
                        event_ref = ref_aln[il];
	                    event_read = read_aln[il];
				
						// COUNTER CHANGES IF EVENT REF IS A GAP
						if (event_ref.compare("-") == 0){ 
                	    	il++;
                        	while (ref_aln[il] == '-'){
                            	event_read += read_aln[il];
                                if (strand.compare("+") == 0){
                                    pos_read++;
          		              	}else{
                        	        pos_read--;
                                }
                            	il++;
         	        		}
                	    	il--;
						}
						// IF THERE IS VARIATION
	            		// AND THE POSITION IS FOUND IN A VARIABLE REGION (FROM TOTAL ALN), THE INFO IS STORED
	                    it_aln = std::find(subset_ev_pos.begin(), subset_ev_pos.end(), pos_ref);
        	            if (it_aln != subset_ev_pos.end()){
	                		event_validation_perbase(ref_aln, read_aln, chr_id, pos_ref, pos_read, st_ref, end_ref, strand, il, event_ref, event_read, aln_summary_perbase, info_reads_aln.at(i).type, 1);
                        }
					}
					
					// REF AND READ POSITION COUNTER 
					// IF THE REF IS A GAP IN THE LAST POSITION THE COUNTER HAS NOT BEEN INCREASED
                	if (ref_aln[il] != '-'){ 
                    	pos_ref++;
           	        }
                	if(read_aln[il] != '-'){
                    	if (strand.compare("+") == 0){
                        	pos_read++;
                	    }else{
                        	pos_read--;
                        }
                   	 }
        	}
 //               std::cout << read_id << "\n" << read_aln << "\n" << ref_aln << "\n" << cons_aln << std::endl;

	}
}


// CALCULATE OPTIMAL LENGTH FOR PARTIAL ALIGNMENTS
// ONLY KEEP REGIONS THAT ARE IN MORE THAN TOTAL-X NUMBER OF TOTAL ALN
// ALL THE EVENTS FOUND IN REGIONS WITH MORE THAN ONE DIFERENT START/END WOULD BE DISCARDED
void consensus_event_per_region(std::vector<AlnSum> &aln_summary, unsigned short thr_read, const int &end_end, RetrReadRow &rrr, const int kmer, std::ofstream &file_reg, std::vector<unsigned long long> &subset_ev_pos, unsigned int &max_ref, unsigned int &max_read, string flag_ind, unsigned int &info_events_size, unsigned int dist_aln)
{
	string ref_data, event_data, all_data, event;
	unsigned int count, st_local, end_local;
	std::map <unsigned int, unsigned int> st_pairs, end_pairs, dist_read, dist_ref;
	std::map <string, unsigned int> final_total_counts;
	std::map <string, string> final_pairs;
	std::vector<FinalEventsRegion> info_events;
	std::vector<string> split_ev;
	long long unsigned int end, st;
	std::map<string, unsigned int >::iterator got_total;
    std::map<unsigned int, unsigned int>::iterator got_dist;
    std::map<unsigned int, unsigned int>::iterator it_max;
 

	// COUNTING CASES PER EVENT, EVERY EVENT IS IDENTIFIED POR THE POSITIONS IN THE REFERENCE GENOME, AND THE SEQUENCE IN THE REFERENCE AND IN THE READ
	for (unsigned i = 0; i < aln_summary.size(); i++)
        {

		// TO CALCULATE THE OPTIMAL LENGTH OF PARTIAL ALIGNMENTS
		// DIST_REF WILL CALCULATE THE NUMBER OF TOTAL ALN THAT SUPPORT A EVENT AT X NT
		// DIST_REF[10] = 3, 3 TOTAL ALN SUPPORTED AN EVENT TO 10BP FROM REF START
		got_dist = dist_ref.find(aln_summary.at(i).distance_ref);
        if (got_dist != dist_ref.end()){
        	dist_ref[aln_summary.at(i).distance_ref]++;
        }else{
            dist_ref[aln_summary.at(i).distance_ref] = 1;
        }

        got_dist = dist_read.find(aln_summary.at(i).distance_read);
        if (got_dist != dist_read.end()){
        	dist_read[aln_summary.at(i).distance_read]++;
        }else{
            dist_read[aln_summary.at(i).distance_read] = 1;
        }

		// COUNT THE NUMBER OF TOTAL ALN SUPPORTING EVERY EVENT
		ref_data = std::to_string(aln_summary.at(i).st_ref) + ":" + std::to_string(aln_summary.at(i).end_ref);
		event_data =  aln_summary.at(i).event_ref + ":" +  aln_summary.at(i).event_read;
		all_data = ref_data + ":" + event_data;
		
		got_total = final_total_counts.find(all_data);
	    if (got_total != final_total_counts.end())
	    {
            final_total_counts[all_data]++;
        }else{
	        final_total_counts[all_data] = 1;
        }
	}

    // CALCULATE THE OPTIMAL LENGTH OF PARTIAL ALIGNMENTS
    // THE MAXIMUM DISTANCE TO COVER THE EVENT IN AT LEAST ONE TOTAL ALIGNMENT
    for (it_max = dist_ref.begin(); it_max != dist_ref.end(); it_max++){
		if (it_max -> first > max_ref){
            max_ref = it_max -> first;
        }
    }
    max_ref += dist_aln;

    for (it_max = dist_read.begin(); it_max != dist_read.end(); it_max++){
		if(it_max -> first > max_read){
            max_read = it_max -> first;
        }
    }
    max_read += dist_aln;        


	// ONLY FOR VSR REGIONS DISCOVERY IN CHILD (ONE INDIVIDUAL)
 //	if (flag_ind.compare("CHILD") == 0){                
		// ONLY KEEP REGIONS THAT ARE IN MORE THAN TOTAL-X NUMBER OF TOTAL ALN
		for (std::map<string, unsigned int>::iterator it = final_total_counts.begin(); it!=final_total_counts.end(); ++it)
		{
			event = it -> first;
			count = it -> second;
			if (count >= thr_read){
				FinalEventsRegion record;
				split_ev.clear();
				split(event, ':', split_ev);
				record.st_ref = atoi(split_ev[0].c_str()); 
				record.end_ref = atoi(split_ev[1].c_str());
				record.event_ref = split_ev[2];
				record.event_qry = split_ev[3];
				record.count_total = final_total_counts[event];
				info_events.push_back(record);
		
			}
		}	
	

		// ALL THE EVENTS FOUND IN REGIONS WITH MORE THAN ONE DIFERENT START/END WOULD BE DISCARDED
		// THESE FILTERS OUT:
		// TRIALLELIC SNPS
		// REGIONS WITH HETEROZYGOUS INDELS EACH VERSION OF A DIFFERENT SIZE
		// REGIONS WITH SNP'S AND INDELS AT THE SAME POSITION
		//
		for (unsigned i = 0; i <info_events.size(); i++){
			st_local = info_events.at(i).st_ref;
			end_local = info_events.at(i).end_ref;
			
			// To check start redundancy
        	std::map<unsigned int, unsigned int>::iterator got_st = st_pairs.find(st_local);
			if (got_st != st_pairs.end())
	        {
				st_pairs.erase(got_st); 
				st_pairs.insert( std::pair<unsigned int, unsigned int>(st_local, 0) );
			}else{
				st_pairs[st_local] = end_local;
			}
	
			// To check end redundancy
			std::map<unsigned int, unsigned int>::iterator got_end = end_pairs.find(end_local);
	        if (got_end != end_pairs.end())
        	{
	        	end_pairs.erase(got_end); 
        	    end_pairs.insert( std::pair<unsigned int, unsigned int>(end_local, 0) );
	        }else{
        	    end_pairs[end_local] = st_local;
         	}
		}
	
		for (unsigned i = 0; i <info_events.size(); i++){
			st_local = info_events.at(i).st_ref;
	        end_local = info_events.at(i).end_ref;
			if (st_pairs[st_local] == end_local && end_pairs[end_local] == st_local){
					// WRITE VARIATION REGIONS TO OUTPUT
					file_reg << rrr.chr << "\t" << rrr.st_ev << "\t" << rrr.end_ev << "\t" << info_events[i].st_ref << "\t" << info_events[i].end_ref << "\t" << info_events[i].event_ref << "\t" << info_events[i].event_qry << "\t" << info_events[i].count_total << "/" << end_end << std::endl;
					// STORE ALL VARIABLE POSITIONS
					for (unsigned long long pos = info_events[i].st_ref; pos <= info_events[i].end_ref; pos++){
						subset_ev_pos.push_back(pos);
  	        	  	}
				}
			}
//		}
		info_events_size = info_events.size();
}



// VALIDATE MISMATCH REGION
void event_validation (const string &ref_aln, const string &read_aln, const string &chr, int st_event_ref,  int st_event_read, int pos_ref, int pos_read, int st_ref, int end_ref, const string strand, int j, string &event_ref, string &event_read, std::vector<AlnSum> &aln_summary, const string &type, int pos_read_hit, string read_id)
{
	long long unsigned int end_event_ref;
	int end_event_read;
	AlnSum record;
	
	// CORRECT THE END POSITION, THE LAST NUCLEOTIDE HAS INCREASED THE COUNTER
	// THE EVENT END IS ONE POSITION BEFORE THE ACTUAL END POSITION
	// EXCEPT IN THOSE CASES IN WHICH THE EVENT ENDS IN GAP AND CONTAINS EXCLUSIVELY GAPS
	// DO FOR READ ADN R
	if (ref_aln[j-1] == '-'){
		// NOT CHANGE ONLY IF THE EVENT HAS EXCLUSIVELY CONTAINED GAPS
		if (end_event_ref == st_event_ref){  
			end_event_ref = pos_ref;
		}else{ 
			end_event_ref = pos_ref -1;
		}
	}
	else{ 
		end_event_ref = pos_ref - 1;
	}

	if (read_aln[j-1] == '-'){
		end_event_read = pos_read;
	}else{
		if (strand.compare("+") == 0){
			end_event_read = pos_read - 1;
		}else{
			end_event_read = pos_read + 1;
		}
		
	}

	// IF THE EVENT IS A GAP EITHER AT THE BEGINNING OR AT THE END OF THE REFERENCE OR READ SEQUENCE
	// THE EVENT IS NOT CONSIDERED 
	unsigned gap_read = 1;
	unsigned gap_ref = 1;
	if(st_event_ref == st_ref || end_event_ref == end_ref  || st_event_read == 1 || end_event_read == read_aln.length() || end_event_read <= 0 ){
		for (unsigned int g = 0; g < event_ref.length(); g++){
			if (event_ref[g] != '-') { gap_ref = 0; }
		}
		for (unsigned int g = 0 ; g < event_read.length(); g++){
			if (event_read[g] != '-') { gap_read = 0;}
		}
	}else{
		gap_read = 0;
		gap_ref = 0;
	}

	// ALL THE EVENTS ARE RECORDED IN aln_summary
	if (gap_read == 0 && gap_ref == 0 ){
		record.st_ref = st_event_ref;
		record.end_ref = end_event_ref;
		record.st_read = st_event_read;
		record.end_read = end_event_read;
		record.event_ref = event_ref;
		record.event_read = event_read;
		record.strand = strand;
		record.chr = chr;
		record.type = type;
	
		// CALCULATE THE DISTANCE FROM THE START OF REFERENCE VSR TO THE END OF THE EVENT
		// CALCULATE THE DISTANCE FROM THE START OF THE READ  TO THE END OF THE EVENT IN THE READ
		// IN ORDER TO RECORD THE OPTIMAL DISTANCE FOR REFERENCE AND READ PARTIAL ALIGNMENTS
	    record.distance_ref = abs(end_event_ref - st_ref);
	    record.distance_read = abs(end_event_read - pos_read_hit);
			
		aln_summary.push_back(record);
//		std::cout << "REGION\t" << event_ref << "\t" << event_read << "\t" << type << std::endl;
	}
}

// ALN FOR READS THAT CONTAIN BOTH PRE AND POST CS
void end_end_aln( std::map<string, string> &reads_id, std::map<string, unsigned int> st_cs, std::map<string, unsigned int> end_cs, std::map<string, string> strand_list, RetrReadRow rrr, string &chr_fa, string &chr_id, std::vector<AlnSum> &aln_summary, unsigned int &end_end, int kmer,  std::vector<InfoAln> &info_reads_aln, std::map<string, string> &read_type_vector, string flag_ind, unsigned int &i_pos){
	string ref_seq, read_id, read_seq, ref_seq_up, ref_aln, read_aln, type, cons_aln, read_type, strand;
	int st_read_aln, end_read_aln, st_ref, pos_read, pos_read_st_aln;
	unsigned int length_ref, length_read, end_ref, pos_ref, st_event_ref, st_event_read, read_length;
	std::map<string, string>::iterator it;
	InfoAln info_aln;
	
	int ENDFREE = 1;
	int VERBOSE = 0;
	int event = 0;
	string event_ref = "";
	string event_read = "";
  
	// THIS LOOPS ACROSS EVERY READ
	for (it=reads_id.begin(); it!=reads_id.end(); ++it){
		// GET READ ATTRIBUTES
		read_id = it->first;
		read_seq = it->second;
		read_length = read_seq.length();
		strand = strand_list[read_id];
		std::map<string, string >::iterator got_type = read_type_vector.find(read_id);
		read_type = got_type->second;
	
		// ONLY READS WITH BOTH COIN-STRINGS
		// DEFINING START AND END POSITIONS IN THE READ TO DO THE ALIGNMENT
		// DEFININD READ-ALN BASIC INFORMATION
		if (read_type.compare("BOTH") == 0) { 
			if (strand.compare("+") == 0){ 
				st_read_aln = st_cs[read_id];
				end_read_aln = end_cs[read_id] + kmer - 1;
			}else{
				st_read_aln = end_cs[read_id];
				end_read_aln = st_cs[read_id] + kmer - 1;
			}

			// FILTER AMBIGUOUS CASES
			if (end_read_aln <= st_read_aln){	
				continue;
			}
			
			// COUNT THE NUMBER OF READS WITH BOTH CS´S 
			type = "TOTAL";
			end_end++;
			
			info_aln.read_id = read_id;
			info_aln.st_read = st_read_aln;				
			info_aln.end_read = end_read_aln;
			info_aln.type = type;
			info_aln.read_length = read_length;
		}else{
			// STORE READ-ALN INFO AND SKIP DOING ALIGNMENT
			// ONLY THE START OR END POSITION WILL BE RECORDED (DEPENDING ON THE STRAND)
	        if (strand.compare("+") == 0){ 
	        	st_read_aln = st_cs[read_id];
                end_read_aln = 0;
             }else{
                st_read_aln = 0;
                end_read_aln = st_cs[read_id] + kmer - 1;
	        }
	        type = "PARTIAL";
			info_aln.read_id = read_id;
            info_aln.st_read = st_read_aln;          
            info_aln.end_read = end_read_aln;
			info_aln.type = type;
			info_aln.read_length = read_length;
			info_aln.read_aln = "NA";
			info_aln.ref_aln = "NA";
			info_reads_aln.push_back(info_aln);
			continue;
		}
			
		// DEFINING POSITIONS IN THE REFERENCE TO DO THE ALIGNMENT
		st_ref = rrr.st_ev;
		end_ref = rrr.end_ev + kmer - 1;

		// DEFINING POSITIONS IN THE READS TO DO THE ALIGNMENT
		length_read = end_read_aln - st_read_aln + 1;
		length_ref = end_ref - st_ref + 1;

         // IF THE DIFFERENCE BETWEEN THE READ REGION AND REFERENCE REGION IS BIG
         // THE ALIGNMENTS IN SUCH CASES WILL NOT BE CORRECT
		 if ( abs(length_read - length_ref) >= 20) { continue; } 


		// INCONGRUENCY OF READ START OR READ END
		if( st_read_aln < 1 || end_read_aln > read_length){ continue;}

		
		// OBTAINING THE SEQUENCES TO BE ALIGNNED, BOTH IN UPPERCASE AND IN THE SAME STRAND
		ref_seq = chr_fa.substr(st_ref - 1 , length_ref + 1);
		read_seq = read_seq.substr(st_read_aln - 1, length_read);
	
		ref_seq_up.clear();
		ref_seq_up = to_upper_seq(ref_seq);
		if (strand_list[read_id].compare("-") == 0)
		{
			string rc = reverse_complement(read_seq);
			read_seq = rc;
		}
	
		ref_aln.clear();
		read_aln.clear();
		cons_aln.clear();
		
		// DOING THE ALIGNMENT, FUNCTION OBTAINED FROM GIUSSEPE'S WORK
		global_align_aff(ref_seq_up, read_seq, ref_aln, read_aln, ENDFREE, VERBOSE);	

		info_aln.read_aln = read_aln;
        info_aln.ref_aln = ref_aln;
		info_reads_aln.push_back(info_aln);


		// LOOKING UP FOR THE DIFFERENCES BETWEEN THE REFERENCE AND THE READ
		// INITIALIZING START POSITION
		if (strand.compare("+")  == 0){
			pos_read = st_read_aln;
		}else{
			pos_read = end_read_aln;
		}
		pos_ref = st_ref;
		pos_read_st_aln = pos_read;

		// LOOKING UP FOR THE VARIABLE REGIONS
		for (unsigned int il= 0; il < ref_aln.length(); il++)
		{
			if (ref_aln[il] == read_aln[il]){
				cons_aln = cons_aln + " ";
				//IN THE NEXT CONSERVED POSITION THE EVENT IS VALIDATED
				if (event == 1){ 
					event = 0;
					event_validation(ref_aln, read_aln, chr_id, st_event_ref, st_event_read, pos_ref, pos_read,  st_ref, end_ref, strand_list[read_id], il, event_ref, event_read, aln_summary, type, pos_read_st_aln, read_id);
					event_ref.clear();
					event_read.clear();
				}
			}else{
				// THE MISMATCH REGION BETWEEN READ END REFERENCE IS SOTRED
				cons_aln = cons_aln + "*";
				event_ref += ref_aln[il];
				event_read += read_aln[il];
				// STORING A GAP
				if (event_ref.compare("-") == 0){
        	                        il++;
	                                while (ref_aln[il] == '-'){
                	                        event_read += read_aln[il];
                        	                if (strand.compare("+") == 0){
                                	                pos_read++;
                                        	}else{
                                                	pos_read--;
     	                                   }
        	                                il++;
                	                }
                        	        il--;
                   		}	
				// THE START OF THE EVENT IS RECORDED
				if (event == 0){ 
					st_event_ref = pos_ref;
					st_event_read = pos_read;
					event = 1;
				}
			}
			// THE MARKERS OF POSITION ARE CHANGED ACCORDINGLY TO THE STRAND AND IF THE POSITION IS EITHER A GAP OR NOT
			if (ref_aln[il] != '-'){
				pos_ref++;
			}
			if(read_aln[il] != '-'){ 
				if (strand_list[read_id].compare("+") == 0){
					pos_read++;
				}else{
					pos_read--;
				}
			}
			
			
		}
		// VALIDATE EVENT AT THE END OF THE SEQUENCE
		if (event == 1){
			event = 0;
			event_validation(ref_aln, read_aln, chr_id, st_event_ref, st_event_read, pos_ref, pos_read, st_ref, end_ref, strand_list[read_id], ref_aln.length(), event_ref, event_read, aln_summary, type, pos_read_st_aln, read_id);
			event_ref.clear();
            event_read.clear();
		}
	}
}


// ELIMINATE PCR DUPLICATES
void eliminate_duplicates( std::vector<RetrReadRow> &rrr_congruent, std::vector<unsigned int> &duplicate, std::map<string, unsigned int> st_cs, std::map<string, unsigned int> end_cs,  std::map<string, string> strand_list,  std::map<string, string> &reads_ids,  std::map<string, string> &reads_ids_prov, string flag_el_duplicate){
	unsigned int pos = 0;
	unsigned int read_length_ev, read_length_analyze;
	std::vector<string> duplicate_specific;
	string read_id;

	// IF ELIMINATE DUPLICATES FLAG IS SET TO FALSE
	if (flag_el_duplicate.compare("FALSE") == 0){
    	for (unsigned int ev = 0; ev < rrr_congruent.size(); ev++){
    		// TO SAVE ONLY ONCE EACH READ, ONLY THE HITS WITH THE PRE-CS ARE STORED
			if (rrr_congruent.at(ev).pos == rrr_congruent.at(ev).end_ev){ continue; }  
			read_id = rrr_congruent.at(ev).id_read;
            reads_ids.insert( std::pair <string, string>( read_id, reads_ids_prov.find(read_id)->second) );
        }
    }else{
    	// IF ELIMINATE DUPLICATES FLAG IS SET TO TRUE
		for (unsigned int ev = 0; ev < rrr_congruent.size(); ev++){
			read_length_ev = rrr_congruent.at(ev).read_length;
			read_id = rrr_congruent.at(ev).id_read;
			// THE DUPLICATES WILL BE CALCULATED BASED ONLY ON PRE-CS HITS
			// IF THE HIT CORRESPONDS TO THE PRE-CS THEN CONTINUE
			if (rrr_congruent.at(ev).pos == rrr_congruent.at(ev).end_ev){ continue; }
			
			// IF THE READ HAS BEEN MARKED AS A DUPLICATE THEN CONTINUE 
			if (std::find( duplicate.begin(), duplicate.end(), ev) != duplicate.end()) { continue; } 
			
			duplicate_specific.clear();
			duplicate_specific.push_back(read_id);

			// TRANSLATE MINUS STRAND POSITION TO FORWARD POSITIONS
			if(rrr_congruent.at(ev).strand.compare("+") == 0){
				pos = rrr_congruent.at(ev).pos_read;
			}else{
				pos = read_length_ev - (rrr_congruent.at(ev).pos_read + 25 -1) + 1;
			}

			// LOOK FOR DUPLICATED READ IN THE NEXT PART OF THE READ ARRAY
			for (unsigned analyze = ev + 1; analyze < rrr_congruent.size(); analyze++){
				// IF THE HIT CORRESPONDS TO THE PRE-CS THEN CONTINUE
				if (rrr_congruent.at(analyze).pos == rrr_congruent.at(analyze).end_ev){ continue; }  

				read_length_analyze = rrr_congruent.at(analyze).read_length;
				// FOR PLUS STRAND HITS COMPARE ALN CS-READ POSITION
				if (rrr_congruent.at(analyze).strand.compare("+") == 0 && rrr_congruent.at(analyze).pos_read == pos){
					duplicate.push_back(analyze);
					duplicate_specific.push_back(rrr_congruent.at(analyze).id_read);
				}
				// FOR MINUS STRAND HITS COMPARE ALN CS-READ POSITION
				if(rrr_congruent.at(analyze).strand.compare("-") == 0 && (read_length_analyze - (rrr_congruent.at(analyze).pos_read + 25 - 1) + 1) == pos){
					duplicate.push_back(analyze);
                    duplicate_specific.push_back(rrr_congruent.at(analyze).id_read);
				}
			}
				
			// GETTING THE FIRST OF THE READS MARKED AS DUPLICATES FOR THAT CS-READ POSITION
			// THE ID OF THE READ IS LINKED WITH ITS SEQUENCE
			string id = duplicate_specific.at(0);
			reads_ids.insert( std::pair <string, string>( id, reads_ids_prov.find(id)->second) ); 
		}
	}

}


// RETURN A std::vector<RetrReadRow> WITH ALL READ THAT PASS CONGRUENCY FILTERS:
// NO TWICE PRE-CS ALN, NO TWICE POST-CS ALN, NO INCONGRUENT STRAND
// IT ALSO RETURNS THE READ TYPE
// BOTH: READ WITH PRE AND POST-CS, PREVIOUS: READ WITH ONLY PRE-CS
void check_congruency( unsigned int event_start_i, unsigned int i, std::vector<RetrReadRow> &rrr_registry, std::map<string, unsigned int> st_cs, std::map<string, unsigned int> end_cs,  std::map<string, string> strand_list, std::vector<RetrReadRow> &rrr_congruent, std::map<string, string> &read_type){
        unsigned int pos = 0;
        string read_id, type;
        for (unsigned int ev = event_start_i; ev < i; ev++){
        	read_id = rrr_registry.at(ev).id_read;
            type = "NA";
			std::map<string, unsigned int >::iterator got_st = st_cs.find(read_id);
            std::map<string, unsigned int >::iterator got_end = end_cs.find(read_id);

            if (got_st != st_cs.end() && got_st -> second == 0){ continue; }
            if (got_end != end_cs.end() && got_end -> second == 0){ continue; }
            if (strand_list[read_id].compare("NA") == 0){ continue; } 

			if (got_st != st_cs.end() && got_end != end_cs.end()) { type = "BOTH"; }
			if (got_st != st_cs.end() && got_end == end_cs.end()) { type = "PREVIOUS"; }
			if (got_st == st_cs.end() && got_end != end_cs.end()) { continue; }

			rrr_congruent.push_back(rrr_registry.at(ev));
            read_type[read_id] = type;
		}
}
		
/////
// MAIN FUNCTION
/////


int main(int argc, char *argv[])
{

  int c, kmer;
  unsigned short thr_read, thr_read_perbase;
  unsigned int dist_aln;
  string read_info_file, chr_dir, chr_sufix, subset_event_file, flag_el_duplicate, flag_ind, out_total, out_perbase;

 while ((c = getopt (argc, argv, "ha:r:s:k:t:p:l:c:d:i:o:q:")) != -1)
  {
   switch(c)
 {
  case 'a':
   read_info_file = std::string(optarg);
   break;
  case 'r':
   chr_dir = std::string(optarg);
   break;
  case 's':
   chr_sufix = std::string(optarg);
   break;
  case 'k':
   kmer = atoi(optarg);
   break;
  case 't':
   thr_read = atoi(optarg);
   break;
  case 'p':
   thr_read_perbase = atoi(optarg);
   break;
  case 'l':
   dist_aln = atoi(optarg);
   break;
  case 'c':
   subset_event_file = std::string(optarg);
   break;
  case 'd':
   flag_el_duplicate = std::string(optarg);
   break;
  case 'i':
   flag_ind = std::string(optarg);
   break;
  case 'o':
   out_total = std::string(optarg);
   break;
  case 'q':
   out_perbase = std::string(optarg);
   break;
  case 'h':
    fprintf(stderr, "USAGE:\n\tAlign_read_RG.out -a READALN-file -r REFDIR/ -s sufixREF -k k -t ReadRegion -p ReadPerbase -l LengthPartial -c ChildSNV -d RemoveDupFlag -i Ind -o out.total -q out.perbase");
   fprintf (stderr, "\nPARAMETERS.\n\ta [READALN-FILE]\tFile per chromosome that contains information about the alignment between the read and the respective SignatureCS");
   fprintf (stderr, "\n\tr [REFDIR]\t\tRG directory path");
   fprintf (stderr, "\n\ts [sufixREF]\t\tSufix of the fasta reference files");
   fprintf (stderr, "\n\tk [k]\t\t\tCSs size = Jellyfish database kmer size");
   fprintf (stderr, "\n\tt [ReadRegion]\t\tEach mismatch region between the reads and the RG must be supported by at least ReadRegion different reads");
   fprintf (stderr, "\n\tp [ReadPerbase]\t\tEach mismatch nucleotide between the reads and the RG must be supported by at least ReadPerbase different reads");
   fprintf (stderr, "\n\tl [LengthPartial]\tLength of extensión of partial  alignment (set to 10)");
   fprintf (stderr, "\n\tc [ChildSNV]\t\tSet to NA (Parameter used in family-framework)");
   fprintf (stderr, "\n\td [RemoveDupFlag]\tSet to TRUE to remove PCR duplicates (FALSE otherwise)");
   fprintf (stderr, "\n\ti [Ind]\t\t\tFor one individual SNV discovery set to CHILD\n\t\t\t\tSet to PARENT for SNV discovery in the parents");
   fprintf (stderr, "\n\to [out.total]\t\tOutput file with regions of polymorphism data\n\t\t\t\tSet to NA for SNV discovery in the parents");
   fprintf (stderr, "\n\tq [out.perbase]\t\tOutput file with single nucleotide polymorphism data\n\n");
   return (0);
  }
 }

	std::ofstream file_reg;
	std::ofstream file_perbase;

	string chr = "";
	string chr_file, fa, hdr;
	
	std::vector<RetrReadRow> rrr_registry, rrr_congruent;
	std::vector<AlnSum> aln_summary;
	std::vector<AlnSumBase> aln_summary_perbase;
	std::vector<unsigned int> duplicate;
	std::vector<unsigned long long> subset, subset_ev_pos, subset_ev_pos_read, subset_ev_pos_child;
	std::vector<InfoAln> info_reads_aln;	

	std::map<string, string> reads_ids, reads_ids_prov;
	std::map<string, unsigned int> st_cs;
	std::map<string, unsigned int> end_cs;
	std::map<string, string> strand_list;
	std::map<string, string> read_type;

	unsigned int end_end = 0;
	unsigned int partial = 0;
	unsigned int event_start_i = 0;
	unsigned int i_pos = 0;
	unsigned int max_ref = 0;
    unsigned int max_read = 0;
	unsigned int info_events_size = 0;

	RetrReadRow rrr;

//	if (flag_ind.compare("CHILD") == 0){
       	file_reg.open(out_total);
//	}
    file_perbase.open(out_perbase);
 

	// STORE SIGNATURE-CS MATCH INFORMATION ON READS 
	if (flag_ind.compare("CHILD") == 0){
    	read_RRR(read_info_file, rrr_registry);
    }else{
    // STORE SIGNATURE-CS MATCH INFORMATION ON READS 
    // STORE SUBSET OF POSOTIONS TO BE ANALIZED FROM CHILD SNV FILE
    	read_subset(subset_event_file, subset, subset_ev_pos_child);
	    read_RRR_parent(read_info_file, rrr_registry, subset);
    }

//	std::cout << "REGISTRY\t" << rrr_registry.size() << std::endl;
	
	unsigned int event_base = (rrr_registry.at(0).st_ev);
	// LOAD EVERY CS-READ MATCH INFORMATION
	for (unsigned int i = 0; i < rrr_registry.size(); i++)
 	{
 		// UPLOAD NEXT CHROMOSOME
		if (rrr_registry.at(i).chr.compare(chr) != 0 ) 
		{
			chr = rrr_registry.at(i).chr;
			chr_file = chr_dir + chr + chr_sufix;
			FILE * chr_fasta = fopen(chr_file.c_str(), "r");
			if (!Fasta_Read(chr_fasta, fa, hdr))
			{
				std::cerr << "ERROR: Can't read fasta record from: " << chr_fasta << std::endl;
				exit(1);
			}
		}

		// THE ALIGNMENTS ARE DONE WHEN THE INFORMATION OF ALL THE READS SPANNING A SPECIFIC VSR HAS BEEN RECORDED		
		if(rrr_registry.at(i).st_ev != event_base){

//			std::cout << event_base << std::endl;
			rrr = rrr_registry.at(i-1);
			// KEEP ONLY CONGRUENT READS 
			check_congruency(event_start_i, i, rrr_registry, st_cs, end_cs,  strand_list, rrr_congruent, read_type);

			// TO ELIMINATE THE PCR DUPLICATES, IN ID_READS ONLY THE BEST READ OF ALL THE DUPLICATED WOULD BE KEPT
			eliminate_duplicates(rrr_congruent, duplicate, st_cs, end_cs,  strand_list, reads_ids, reads_ids_prov, flag_el_duplicate);

			// COMPUTE ALN WITH PRE AND POST SIGNATURE CS'S (TOTAL ALN)
			end_end_aln(reads_ids, st_cs, end_cs, strand_list, rrr, fa, rrr.chr, aln_summary, end_end, kmer, info_reads_aln, read_type, flag_ind, i_pos);
			
			// COMPUTE CONSENSUS VARITION FROM TOTAL ALN
			consensus_event_per_region(aln_summary, thr_read, end_end, rrr, kmer, file_reg, subset_ev_pos, max_ref, max_read, flag_ind,  info_events_size, dist_aln);
//			std::cout << "CONSENSUS\t" << aln_summary.size() << "\t" << info_events_size <<  "\t" << subset_ev_pos.size() << std::endl; 

			// ONLY FOR PARENTAL SNV DISCOVERY
			// GET ALL CHILD SNV'S THAT ARE CONTAINED IN THE VSR BEING ANALYZED
			if (flag_ind.compare("CHILD") != 0){
//	                      std::cout << "PARENT\t" << subset_ev_pos_read.size() << std::endl;
	//	    		// FOR THE PARENTS THE SIZE OF PARTIAL ALIGNMENTS
		    		//                         // WILL BE DETERMINED FOR THE MOST 3' VARIATON OF THE CHILD TOTAL ALN +5
				for (i_pos; i_pos < subset_ev_pos_child.size(); i_pos++){
        				if (subset_ev_pos_child[i_pos] >= rrr.st_ev_child && subset_ev_pos_child[i_pos] <= (rrr.end_ev_child + kmer)){
            					subset_ev_pos_read.push_back(subset_ev_pos_child[i_pos]);
	        			//	std::cout << subset_ev_pos_child[i_pos] << std::endl;
					}
        				if (subset_ev_pos_child[i_pos] > (rrr.end_ev_child + kmer)){
            					i_pos = i_pos-1;
			    			break;
	        			}
        			}
  //                             std::cout << "PARENT\t" << subset_ev_pos_read.size() << std::endl;
				
//				std::cout << subset_ev_pos_read.back() << "\t"<< event_base <<"\t" << dist_aln << std::endl;
				// MUST BE CHILD EVENTS IN THAT READ
				if (subset_ev_pos_read.size() > 0 ){
	                                // IF THERE ARE TOTAL ALN WITH VARIATION, NO REGION INSONSISTENCY MUST EXist
					if (subset_ev_pos.size() == 0){
        	        			if (event_base < subset_ev_pos_read.back() ){
	        	        	        	max_ref = subset_ev_pos_read.back() - event_base + dist_aln;
                                        		max_read = max_ref;
						}else{
							max_ref = event_base - subset_ev_pos_read.back() + dist_aln;
							max_read = max_ref;
						}
                                	}
                                	end_end = 0;
//                               	std::cout << "CICLO PARENT" << std::endl;
					if ( (info_events_size == 0) || ( info_events_size > 0 && subset_ev_pos.size() > 0) ){                                // COMPUTE VARIATION PER NUCLEOTIDE EXCLUDING GAPS
                                                aln_perbase(reads_ids, st_cs, end_cs, strand_list, rrr_registry.at(i-1), fa, rrr_registry.at(i-1).chr, aln_summary_perbase, end_end, partial, kmer, subset_ev_pos_read, info_reads_aln, max_ref, max_read, flag_ind);
                                                consensus_event_perbase(aln_summary_perbase, thr_read_perbase, end_end, partial, rrr_registry.at(i-1), kmer, file_perbase);

					}
				}

//			}

//			std::cout << "PARENT\t" << subset_ev_pos_read.size() << std::endl;
	
			// FOR THE PARENTS THE SIZE OF PARTIAL ALIGNMENTS
			// WILL BE DETERMINED FOR THE MOST 3' VARIATON OF THE CHILD TOTAL ALN +5
//			if (flag_ind.compare("CHILD") !=0){
//				if (subset_ev_pos.size() == 0){ 
//					max_ref = subset_ev_pos_read.back() - event_base + dist_aln;
//					max_read = max_ref;
//				}
//				end_end = 0;
//				std::cout << "CICLO PARENT" << std::endl;
				// IF THERE ARE TOTAL ALN WITH VARIATION, NO REGION INSONSISTENCY MUST EXIST
//				if (subset_ev_pos_read.size() > 0 ){
//					if ( (info_events_size == 0) || ( info_events_size > 0 && subset_ev_pos.size() > 0) ){ 
				// COMPUTE VARIATION PER NUCLEOTIDE EXCLUDING GAPS
//						aln_perbase(reads_ids, st_cs, end_cs, strand_list, rrr_registry.at(i-1), fa, rrr_registry.at(i-1).chr, aln_summary_perbase, end_end, partial, kmer, subset_ev_pos_read, info_reads_aln, max_ref, max_read, flag_ind);
                // 
  //              				consensus_event_perbase(aln_summary_perbase, thr_read_perbase, end_end, partial, rrr_registry.at(i-1), kmer, file_perbase);
//					}
//				}
			}else{
				if(subset_ev_pos.size() > 0){
//					 std::cout << "CICLO HIJO" << std::endl;
					end_end = 0;
					// COMPUTE VARIATION PER NUCLEOTIDE EXCLUDING GAPS
					aln_perbase(reads_ids, st_cs, end_cs, strand_list, rrr_registry.at(i-1), fa, rrr_registry.at(i-1).chr, aln_summary_perbase, end_end, partial, kmer, subset_ev_pos, info_reads_aln, max_ref, max_read, flag_ind);
					consensus_event_perbase(aln_summary_perbase, thr_read_perbase, end_end, partial, rrr_registry.at(i-1), kmer, file_perbase);
				}	
			}

			duplicate.clear();
			reads_ids.clear();
			reads_ids_prov.clear();
			st_cs.clear();
			end_cs.clear();
			strand_list.clear();
			aln_summary.clear();
			aln_summary_perbase.clear();	
			end_end = 0;
			partial = 0;
			subset_ev_pos.clear();
			subset_ev_pos_read.clear();
			info_reads_aln.clear();
			event_base = rrr_registry.at(i).st_ev;
			max_ref = 0; max_read = 0;
			rrr_congruent.clear();
            read_type.clear();
			info_events_size = 0;
		}
		
		// STORE THE INFORMATION OF ALL THE READS SPANNING A SPECIFIC VSR 
		
		//STORE VSR START
		if( st_cs.size() == 0){
			event_start_i = i;
		}

		// STORE ALN INFORMATION FOR THE PREV-CS
		if (rrr_registry.at(i).pos == rrr_registry.at(i).st_ev){
			std::map<string, unsigned int >::iterator got_st = st_cs.find(rrr_registry.at(i).id_read);
			// A READ WITH THE SAME PREV_CS IN MULTIPLE POSITIONS WOULD BE LABELED WITH 0 IN ST_CS 
			// OTHERWISE THE MAPPING POSITION WOULD BE RECORDED
        	if (got_st != st_cs.end()){
				st_cs[rrr_registry.at(i).id_read] = 0; 
			}else{
				st_cs[rrr_registry.at(i).id_read] = rrr_registry.at(i).pos_read; 
			}
		}

		// STORE ALN INFORMATION FOR THE POST-CS
		if(rrr_registry.at(i).pos == rrr_registry.at(i).end_ev){ 
			std::map<string, unsigned int >::iterator got_end = end_cs.find(rrr_registry.at(i).id_read);
			// A READ WITH THE SAME POST-CS IN MULTIPLE POSITIONS WOULD BE LABELED WITH 0 IN END_CS 
			// OTHERWISE THE MAPPING POSITION WOULD BE RECORDED
        	if (got_end != end_cs.end()){
            	end_cs[rrr_registry.at(i).id_read] = 0;
	        }else{
        	     end_cs[rrr_registry.at(i).id_read] = rrr_registry.at(i).pos_read;
            }
		}
		
		// STORE ALN INFORMATION FOR THE STRAND
		// ALL CS´S MAPPING TO THE SAME READ SHOULD BE IN THE SAME STRAND
		// OTHERWISE THE READ WOULD BE LABELED WITH NA
	 	std::map<string, string>::iterator got_strand = strand_list.find(rrr_registry.at(i).id_read);
        if (got_strand != strand_list.end()){ 
        	if (strand_list[rrr_registry.at(i).id_read] != rrr_registry.at(i).strand){
				strand_list[rrr_registry.at(i).id_read] = "NA"; 
			}
        }else{
        	strand_list[rrr_registry.at(i).id_read] = rrr_registry.at(i).strand;
        }

		// STORE ID-READ, SEQ-READ PAIR INFORMATION
		reads_ids_prov.insert( std::pair <string, string>(rrr_registry.at(i).id_read, rrr_registry.at(i).seq_read) ); //THE ID OF THE READS ARE LINKED WITH ITS SEQUENCE
	}

	// PROCESS LAST VSR
	unsigned int i;
	i = rrr_registry.size();
	rrr = rrr_registry.at(i-1);
	check_congruency(event_start_i, i, rrr_registry, st_cs, end_cs,  strand_list, rrr_congruent, read_type);
	eliminate_duplicates(rrr_congruent, duplicate, st_cs, end_cs,  strand_list, reads_ids, reads_ids_prov, flag_el_duplicate);
	end_end_aln(reads_ids, st_cs, end_cs, strand_list, rrr_registry.at(i-1), fa, rrr_registry.at(i-1).chr, aln_summary, end_end, kmer, info_reads_aln, read_type, flag_ind, i_pos);
    consensus_event_per_region(aln_summary, thr_read, end_end, rrr_registry.at(i-1), kmer, file_reg, subset_ev_pos, max_ref, max_read, flag_ind, info_events_size, dist_aln);

//	std::cout << "CONSENSUS\t" << aln_summary.size() << "\t" << info_events_size <<  "\t" << subset_ev_pos.size() << std::endl; 	
 
	// ONLY FOR PARENTAL SNV DISCOVERY
	// GET ALL CHILD SNV'S THAT ARE CONTAINED IN THE VSR BEING ANALYZED
//	                        std::cout << "PARENT\t" << subset_ev_pos_read.size() << std::endl;
  		// FOR THE PARENTS THE SIZE OF PARTIAL ALIGNMENTS
                // WILL BE DETERMINED FOR THE MOST 3' VARIATON OF THE CHILD TOTAL ALN +5
        if (flag_ind.compare("CHILD") != 0){

		for (i_pos; i_pos < subset_ev_pos_child.size(); i_pos++){
 			if (subset_ev_pos_child[i_pos] >= rrr.st_ev_child && subset_ev_pos_child[i_pos] <= (rrr.end_ev_child + kmer)){
            			subset_ev_pos_read.push_back(subset_ev_pos_child[i_pos]);
	       		 }
        		if (subset_ev_pos_child[i_pos] > (rrr.end_ev_child + kmer)){
            			i_pos = i_pos-1;
		  		break;
	  		}
        	}
  //                std::cout << "PARENT\t" << subset_ev_pos_read.size() << std::endl;

	// MUST BE CHILD EVENTS IN THAT READ
		if (subset_ev_pos_read.size() > 0 ){
	        	// IF THERE ARE TOTAL ALN WITH VARIATION, NO REGION INSONSISTENCY MUST EXist
			if (subset_ev_pos.size() == 0){
//				std::cout << subset_ev_pos_read.back() << "\t"<< event_base << std::endl;
				max_ref = subset_ev_pos_read.back()- event_base + dist_aln;
                        	max_read = max_ref;
  //             			std::cout << "MAX " << max_ref << "\t" << max_read << std::endl; 
	      		}
               		end_end = 0;
    //           		std::cout << "CICLO PARENT" << std::endl;
			if ( (info_events_size == 0) || ( info_events_size > 0 && subset_ev_pos.size() > 0) ){                                // COMPUTE VARIATION PER NUCLEOTIDE EXCLUDING GAPS
                                   aln_perbase(reads_ids, st_cs, end_cs, strand_list, rrr_registry.at(i-1), fa, rrr_registry.at(i-1).chr, aln_summary_perbase, end_end, partial, kmer, subset_ev_pos_read, info_reads_aln, max_ref, max_read, flag_ind);
                                    consensus_event_perbase(aln_summary_perbase, thr_read_perbase, end_end, partial, rrr_registry.at(i-1), kmer, file_perbase);
			}
		}


	}else{
		if(subset_ev_pos.size() > 0){
//			 std::cout << "CICLO HIJO" << std::endl;
			end_end = 0;
			// COMPUTE VARIATION PER NUCLEOTIDE EXCLUDING GAPS
			aln_perbase(reads_ids, st_cs, end_cs, strand_list, rrr_registry.at(i-1), fa, rrr_registry.at(i-1).chr, aln_summary_perbase, end_end, partial, kmer, subset_ev_pos, info_reads_aln, max_ref, max_read, flag_ind);
			consensus_event_perbase(aln_summary_perbase, thr_read_perbase, end_end, partial, rrr_registry.at(i-1), kmer, file_perbase);
		}	
	}
	

	file_reg.close();
	file_perbase.close();

	return 0;
}


