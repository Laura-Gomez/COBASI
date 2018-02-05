################
################ OBTAIN LANDSCAPE
################ 
#
# This script filters out the positions for every non-unique kmers. It outputs a list of start positions for every unique kmer (CSs) asociated wit its count.
#
#       COMMAND LINE.
#  perl Compute_Landscape.pl -cov OUTPUT.COVERAGE -unique DIR_CSs/ -suf_unique SUFIX -out_dir DIR_OUT /
#
# PARAMETERS.
# 	cov 			Coverage file
#	unique 			The directory containing the files with the start and end positions of the non-overlapping unique region
#	suf_unqiue 		Sufix of the unique regions files (SET to SUFIX =_cs_uniq_regions.tab)
#	out_dir			Output directory. 
#
#  The output files will be named {RG chromosome}{.land}
#
# OUTPUT.
# The Variation Landscape (VL) file per chromosome. It contains how many times each CS is found in the sequencing reads. The VL file contains two columns: 
#	column 1, RG position
#	column2, the number of ocurrences per each CS along the reads.
#
# COVERAGE FILE HEAD
#  >chrUn_KI270385v1
#  1	139309
#  2	139570
#  ...
################
################ 
################


use Getopt::Long;
use strict;

my (%opts, %st, %end);
GetOptions(\%opts, 'h', "cov=s", "unique=s", "suf_unique=s", "out_dir=s");


my $help = "USAGE:
perl Compute_Landscape.pl -cov OUTPUT.COVERAGE -unique DIR_CSs/ -suf_unique SUFIX -out_dir DIR_OUT /\n
PARAMETERS.
	cov 			Coverage file
	unique 		The directory containing the files with the start and end positions of the non-overlapping unique region
	suf_unqiue 		Sufix of the unique regions files(SET to SUFIX =_cs_uniq_regions.tab)		
	out_dir		Output directory. The output files will be named {RG chromosome}{.land}\n
OUTPUT.
	The Variation Landscape (VL) file per chromosome. It contains how many times each CS is found in the sequencing reads. The VL file contains two columns: 
	•	column 1, RG position
	•	column2, the number of ocurrences per each CS along the reads.\n";

die $help if $opts{h};



my ($line, $chr, $cs_file, $out_file, $pos_ref, $cov, $r);
my (@regions, @st, @end, @data);

open (COV, $opts{'cov'});
while(<COV>){
	chomp($_);
	$line = $_;
	
	#### IN THE FIRST LINE OF COVERAGE FILE IS THE CHROMOSOME INFORMATION
	#### RETRIEVE CHROMOSOME INFORMATION
	if ($line =~ m/>(\w+)/){
		$chr = $1;
		$cs_file = $opts{'unique'}.$chr.$opts{'suf_unique'};
		
		#### OPEN THE FILE WITH THE UNIQUE NON-OVERLAPPING REGIONS FOR THAT SPECIFIC CHROMOSOME
		#### CREATE ST POSITION AND END POSITIONS ARRAYS
		open(UN, $cs_file);
		@regions = <UN>;
		close(UN);
		my (@line_region);
		my ($st_r, $end_r);
		foreach $line (@regions){
			chomp($line);
			@line_region = split(/\t/, $line);
			$st_r = $line_region[0];
			$end_r = $line_region[1];
			chomp($end_r);
			chomp($st_r);
			push ( @st, $st_r);
			push ( @end, $end_r);
			
		}
		#### INITIALIZE COUNTER TO MOVE ON POSITIONS ARRAYS
		$r = 0;
		
		#### OPEN OUTPUT FILE
		$out_file = $opts{'out_dir'}.$chr.".land";
		open(OUT, ">$out_file");
		
	}else{
		#### COMPARE EVERY POSITION IN THE COVERAGE FILE AGAINST THE REGIONS ARRAYS
		@data = split(/\t/, $line);
		$pos_ref = $data[0];
		$cov = $data[1];
		
		#### INCREASE REFERENCE POSITION UNTIL A UNIQUE REGION IS FOUND
		if ($pos_ref < $st[$r]){
			next;
		}
		#### OUTPUT ALL POSITIONS IN THAT UNIQUE REGION
		if ($pos_ref >= $st[$r] && $pos_ref < $end[$r]){
			print OUT "$pos_ref\t$cov\n";
			next;
		}
		#### AT THE END OF THE REGION MOVE FORWARD TO THE NEXT UNIQUE REGION
		if($pos_ref == $end[$r]){
			print OUT "$pos_ref\t$cov\n";
			$r++;
		}
		
	}
}

close(OUT);
close(COV);

