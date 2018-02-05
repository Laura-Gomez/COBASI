################
################ OBTAIN NON-OVERLAPPING UNIQUE REGIONS
################ 
#
# This script concatenates all positions from adjacent unique kmers (found by bowtie). It generates a list of unique non-overlapping regions.
#
#                   COMMAND LINE.
#   perl Compute_CSDB_FromBowtie.pl –dir DIR_BOWTIE/ -sufix .BOWTIE.OUT –out DIR_OUT/
#
# PARAMETERS.
#	dir		Directory that contains the output from bowtie 
#	sufix	Bowtie output files sufix
#	out		Output directory
#
# OUTPUT.
# One file per chr with a list of start and end positions of unique non-overlapping regions 
# The ouput files will have the extension “_cs_uniq_regions.tab”
#
################
################ 
################

use Getopt::Long;
use strict;

my (%opts, %pos);
GetOptions(\%opts, 'h', "dir=s", "sufix=s", "out=s");

my $help = "USAGE:
	perl Compute_CSDB_FromBowtie.pl –dir DIR_BOWTIE/ -sufix .BOWTIE.OUT –out DIR_OUT/\n
PARAMETERS.
	dir		Directory that contains the output from bowtie 
	sufix		Bowtie output files sufix
	out		Output directory\n
OUTPUT.
	One file per chr with a list of start and end positions of unique non-overlapping regions.
	The ouput files will have the extension “_cs_uniq_regions.tab”\n\n";


die $help if $opts{h};

#### GET ALL FILE WITH SUFIX AT THE END OF THE NAME
opendir(DIR, "$opts{'dir'}");
my @bowtie_files = grep{/$opts{'sufix'}$/} readdir(DIR);
closedir(DIR);

my ($file, $chr, $file_name, $file_out);

foreach $file (@bowtie_files){
    $file =~ m/^(.+)\.$opts{'sufix'}/;
	$chr = $1;

	#### MERGE FILE PATH
	#### CREATE OUTPUT FILE PATH WITH EXTENXION [_cs_uniq_regions.tab]
	$file_name = $opts{'dir'}.$file;
	$file_out = $opts{'out'}.$chr."_cs_uniq_regions.tab";

	print "$file_name\n$file_out\n";
	
	open(IN, $file_name);

	my ($line, $pos, $st, $last);
	my (@data);
	my $cont = 0;	
	while(<IN>){
        chomp($_);
		$line = $_;
		@data = split(/\t/, $line);
		$pos = $data[0];
		
		#### WHEN THE FIRST LINE OF THE INPUT FILE IS READ:
		#### INITIALIZE THE START AND END OF THE FIRST NON-OVERLAPPING REGION TO POS
		#### OPEN OUTPUT FILE
		if ($cont == 0 ){
			open(OUT, ">$file_out");
			$st  = $pos;
			$last = $pos;
			$cont++;
			next;
		}
		
		#### WHEN A NON-ADJACENT POSITION IS FOUND:
		#### THE LAST CONTINUOUS REGION IS WRITTEN TO OUTPUT
		#### INITIALIZE THE START AND END OF THE NEXT NON-OVERLAPPING REGION TO POS
		if ($pos != $last + 1){
			print OUT "$st\t$last\n";
			$st = $pos;
			$last = $pos;
       	 	}
       	 	
       	#### THE END OF THE CURRENT REGION IS ASSUMED TO BE EVERY LINE
		$last = $pos;
  	}
	close(IN);
	#### THE LAST CONTINUOUS REGION IS WRITTEN TO OUTPUT
	print OUT "$st\t$last\n";
	close(OUT);
}	

