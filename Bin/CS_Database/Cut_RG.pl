################
################ CUT REFERENCE GENOME
################ 
#
#This script will cut one chromosome from the reference genome getting all kmers by sliding-windows of one base pair . 
#
#                            COMMAND LINE.
#           perl Cut_RG.pl -fna  FASTA -k K –out OUT.STR
#
# PARAMETERS.
#	fna		A fasta file [FASTA] with the sequence of one RG chromosome
#	k		kmer size (set to 30)
#	out		Output file [OUT.STR]
#
# OUTPUT.
# A multi-fasta file with the sequence for every kmer along each chromosome of the RG
# The ID for every multi-fasta sequence will correspond to the start position of each kmer
#
################
################ 
################

use Getopt::Long;
use strict;


my (%opts);
GetOptions(\%opts, 'h', "fna=s","k=i", "out=s");

my $help = "USAGE:
	perl Cut_RG.pl -fna  FASTA -k K –out OUT.STR\n
PARAMETERS:
	fna\tA fasta file [FASTA] with the sequence of one RG chromosome
	k\tkmer size (set to 30)
	out\tOutput file [OUT.STR]\n
OUTPUT:
	A multi-fasta file with the sequence for every kmer along each chromosome of the RG
	The ID for every multi-fasta sequence will correspond to the start position of each kmer\n\n";

die $help if $opts{h};


my ($cs, $cs2) = "";
my $pos_chr = 0;
my $k = $opts{'k'};
my $print_pos;

open(OUT, ">$opts{'out'}");

open(IN, $opts{'fna'});
while(<IN>){
	my $line = $_;
	chomp($line);
	if ($line =~ />/){ next; }

	my $pos_line = 0;
	my $extra = "";
	
	#### SLIDING-WINDOW OF ONE BASE PAIR
	for ($pos_line = 0; $pos_line < length($line); $pos_line++){
		if (length($cs) == $k){
			### REMOVE THE FIRST LETTER OF THE LAST RECORDED KMER
			$cs2 = substr($cs, 1);
			$cs = $cs2;
		}
		#### GET THE NEXT LETTER AFTER THE END OF THE LAST RECORDED KMER		
		$extra = substr($line, $pos_line, $k-length($cs));
		$cs = $cs.$extra;
		
		#### OUTPUT THE FINAL KMER
		$print_pos = $pos_chr+1;
		print OUT ">$print_pos\n$cs\n";
		
		#### INITIALIZE THE POINTER TO THE END OF THE FIRST RECORDED KMER
		if ($pos_chr == 0){
			$pos_line += ($k-1);
		}
		$pos_chr++;
	}
}

close(IN);
close(OUT);
