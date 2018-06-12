################
################ GET ONE-DROP REGIONS (INTERNAL PART OF VSR’S)
################ 
#
# This script obtain a list of partial VSR’s.
# It only looks for the internal drop and rise in coverage characteristics of any VSR. 
#
#         COMMAND LINE.
#   perl Get_Onedrop.pl -land VL.LAND -max M -fac RCI -rmin R -fst FST -nt NT –density DEN –var chr1.var
#
# PARAMETERS.
#	land		Input file. VL file
#	max		CS’s with a coverage count higher than M are skipped in the VSR identification process
#	fac		A Relative Coverage Index absolute value higher than RCI is considered as bona fide variation signal (SET RCI=0.3)
#	rmin		The PrevCS and PostCS coverage must be at least R (SET R=10)
#	fst		The median of the coverage values for InterCS’s must be lower than FST (SET FST=third quartile coverage+10(IQR))
#	nt		The length of the VSR must be longer than NT
#	density	DEN is the mínimum density of CSs required inside the VSR 
#	var		Output file
#
# OUTPUT.
# The output file contains either whole VSRs or partial VSRs (only the internal drop and rise). In this step, the PrevCS and PostCS will be the ones before and after this (possibly incomplete) VSR signal. These regions will be extended, if possible, in the next process. This output file is composed of 13 columns:
# 1.	PrevCS position
# 2.	PrevCS count
# 3.	Next CS after PrevCS position
# 4.	Next CS after PrevCS count
# 5.	Number of nucleotides in the internal region of the VSR
# 6.	Number of interCSs  (a VSR with a CS density of 1 will have the same value in columns 5 and 6)
# 7.	Last CS before PostCS position
# 8.	Last CS before PostCS count
# 9.	PostCS position
# 10.	PostCS count
# 11.	RCI value for the VSR start
# 12.	RCI value for the VSR end
# 13.	Median of the coverage values for the InterCS’s
#
################
################ 
################


use Getopt::Long;
use List::Util qw[min max];
use strict;

my (%opts);
my $region = 0;
my $cont = 0;
my (@before, @start, @end, @after, @line, @prev, @median, @median_sort, @nmedian);
my ($coin, $read, $max, $factor, $factor_abs, $cont_coin, $fac_st, $fac_end, $nt, $cs_ratio, $mediana, $size, $nmedian_ref, $class);

GetOptions(\%opts, 'h', "land=s", "fac=f", "rmin=i", "max=i", "fst=i", "out=s", "nt=s", "density=f" );

my $help = "USAGE:
	perl Get_Onedrop.pl -land FILE.LAND -max M -fac RCI -rmin R -fst FST -nt NT –density DEN –out FILE.ONEDROP\n
PARAMETERS:
	land		Input file. VL file [FILE.LAND]
	max		CS’s with a coverage count higher than M are skipped in the VSR identification process
	fac		A Relative Coverage Index absolute value higher than RCI is considered as bona fide variation signal (SET RCI=0.3)
	rmin		The PrevCS and PostCS coverage must be at least R (SET R=10)
	fst		The median of the coverage values for InterCS’s must be lower than FST (SET FST=third quartile coverage+10(IQR))
	nt		The length of the VSR must be longer than NT
	density	 	DEN is the mínimum density of CSs required inside the VSR 
	out		Output file [FILE.ONEDROP]\n
OUTPUT:
	The output file contains either whole VSRs or partial VSRs (only the internal drop and rise). In this step, the PrevCS and PostCS will be the ones before and after this (possibly incomplete) VSR signal. These regions will be extended, if possible, in the next process. This output file is composed of 13 columns:
	1.	PrevCS position
	2.	PrevCS count
	3.	Next CS after PrevCS position
	4.	Next CS after PrevCS count
	5.	Number of nucleotides in the internal region of the VSR
	6.	Number of interCSs  (a VSR with a CS density of 1 will have the same value in columns 5 and 6)
	7.	Last CS before PostCS position
	8.	Last CS before PostCS count
	9.	PostCS position
	10.	PostCS count
	11.	RCI value for the VSR start
	12.	RCI value for the VSR end
	13.	Median of the coverage values for the InterCS’s\n";

die $help if $opts{h};

print "$opts{'land'}\n";
open(LN, "$opts{'land'}");
open(VAR, ">$opts{'out'}");


#### OPEN LANDSCAPE FILE PER CHROMOSOME
while(<LN>){
	chomp($_);
	@line = split(/\s+/, $_);
	$coin = $line[0]; $read = $line[1];
	#print "$read\n";
	
	#### ALL CS'S WITH A COVERAGE ABOVE MAX THRESHOLD WILL NOT BE CONSIDERED
	if ($read >= $opts{'max'}){ next; }
	
	
	#### INTITIALIZING DATA FOR PREVIOUS LINES
	#### TO AVOID CONSECUTIVE PARTIAL VSR'S
	if ($cont == 0 ){
		@prev = @line;
		$cont++;	
		next;
	}
	
	#### RCI IS CALCULATED
	#### IF COVERGE FOR CS N AND CS N+1 EQUAL ZERO, RCI IS SET TO ZERO
	$max= max($read, $prev[1]);
	if ($max == 0){
		$factor = 0;
	}else{
		$factor = ($read - $prev[1])/$max;
	}
	
	
	#### START OR END OF VSR'S
	$factor_abs = abs($factor);
	if($factor_abs >= $opts{'fac'}){
		
		#### EVERY DOWN IN COVERAGE WILL BE CONSIDERED AS A POSSIBLE START OF A PARTIAL VSR
		#### DATA FOR PREV-CS AND START-CS ARE SOTRED
		#### A REGION FLAG, A CS COUNTER,  AND A VECTOR WITH THE COVERAGE FOR ALL INTER-CS'S ARE INITIALIZED
		if($factor < 0  && $prev[1] >= $opts{'rmin'}){
			@before = @prev;
			@start = @line;
			$fac_st = $factor;
	
			$region = 1;
			$cont_coin = 1;			
			@median = ();
			push (@median, $read);
		}

		#### AN UP IN COVERAGE AFTER A DOWN HAS BEEN FOUND WILL BE CONSIDERED AS A POSSIBLE END OF A PARTIAL VSR
		#### IF A PREVIOUS START OF REGIO
		if($factor > 0 && $region == 1 && $read >= $opts{'rmin'} ){	
			#### DATA FOR POST-CS AND END-CS ARE SOTRED		
			@end = @prev;
			@after = @line;
			$fac_end = $factor;

			#### STATISTICS FOR THE PARTIAL VSR ARE OBTAINED
			#### SIZE OF REGION (IN BASE PAIRS)
			#### THE MEDIAN OF THE COVERAGE FOR THE INTER-CS'S
			#### THE DENSITY OF CS'S IN THE REGION
			$nt = $after[0] - $before[0];

			($nmedian_ref, $class) = get_median_class_n($cont_coin);
			@nmedian = @$nmedian_ref;
			@median = sort {$a<=>$b} @median;
			$size = scalar(@median);
			if ($size == 1){
				$mediana = $median[0];
			}else{
				if ($class == 1){
					$mediana = $median[$nmedian[0]];
				}else{
					$mediana = ($median[$nmedian[0]] + $median[$nmedian[1]])/2;
				}
			}
			
			$cs_ratio = $cont_coin/$nt;
			if ($cs_ratio >= $opts{'density'} ){ 	
				if ($nt >= $opts{'nt'} && $mediana < $opts{'fst'}){
					print VAR "$before[0]\t$before[1]\t$start[0]\t$start[1]\t$cont_coin\t$nt\t$end[0]\t$end[1]\t$after[0]\t$after[1]\t$fac_st\t$mediana\t$fac_end\n";
				}
			}
			
			#### REGION FLAG INITIALIZED
			#### FLAG TO AVOID ADJACENT PARTIAL VSR'S
			$region = 0;
			$cont=0;			
		}
			
	}elsif ($region == 1){
		#### STORE COVERAGE DATA FOR INTER-CS'S
		$cont_coin++;
		push (@median, $read);
	}
	
	@prev = @line;	
}
close(LN);
close(VAR);


sub get_median_class_n {
	my @nmedian;
	my $class;
	my $n = $_[0];
	if ($n%2 == 1){
		@nmedian = (int(($n/2)+0.5));
		$class = 1;
	}
	if($n%2 == 0){
		@nmedian = ($n/2, ($n/2)+1);
		$class = 2;
	}
	return(\@nmedian, $class);
}
