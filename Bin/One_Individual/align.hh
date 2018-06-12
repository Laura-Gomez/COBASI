#ifndef ALIGN_HH
#define ALIGN_HH 1


void global_align(const std::string & S, const std::string & T,
                  std::string & S_aln, std::string & T_aln,
                  int endfree, int verbose);

void global_align_aff(const std::string & S, const std::string & T,
                      std::string & S_aln, std::string & T_aln,
                      int endfree, int verbose);

#endif
