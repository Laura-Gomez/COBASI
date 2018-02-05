##################
################## POST-PROCCESING
##################
#
#
# Once the candidate de novo SNVs have been identified, this script will analyze the regions corresponding to the child VSR for the three family individuals to identify undesired patterns.
# 1) regions with low  CS density; 
# 2) regions in which any CS has a coverage higher than expected;
# 3) for any individual regions with low coverage for the CSs corresponding to the child Signature CSs,
# 4) regions with additional peaks inside the region corresponding to the child VSR :
# in the case of the child if there is any additional drop or rise it should correspond to a region with almost no coverage;
# in the case of the parents there should not exist any drop or rise that indicates a posible heterozigosity for the child SNV position or there should not exist a drop and rise that correspond to the exact same child’s VSR boundaries; 
# and 5) for the child, regions with unequal coverage in both sides of the VSR
# 
# COMMAND LINE.
# perl Postprocessing.pl -novo Novo.out –land LAND-DIR/ -chr CHR –parent P -density DEN –max MAX -higher HIGH -cov COV -rci RCI –rmin RMIN -low LOW -peaks PEAK -out Novo.chr.genome.tp
#
#
# PARAMETERS
#	novo		File with de candidate de novo SNVs
#	land		Directory with the VLs (Variation landscapes)
#	chr		    The chromosome to be analyzed [CHR]
#	parent		If the landscape corresponds to the child set [P] to 0
#               If the landscape corresponds to the parent set [P] to 1
#	density   	Regions with a density of at least [DEN] will be kept
#	max/higher	Regions with at most [HIGH] CSs with a coverage higher than
#               [MAX] will be kept
#	cov		    Regions for which the CSs corresponding to the child’s Sginature
#               CSs have a coverage of at least COV will be kept 
#	rci		    A change higher than RCI in the RCI index will be considered 
#		        significative
#	rmin		A change in coverage will be considered significative only if the 
#		        coverage for any CS is no lower than RMIN
#	low	        For the child, a CS with a coverage lower than LOW will be 
#		        considered as a region of “almost no coverage”		
#	peaks		Regions with less than PEAK significative changes in coverage 
#		        inside the VSR will be kept
#	out		    Output file. This will contain all TP de novo SNVs
#
#
# OUTPUT
# Same columns as Novo.out
#
####



use Getopt::Long;
use List::Util qw[min max];
use strict;

my (%opts);
my $region = 0;
my $cont = 0;
my $temp_before;
my (@before, @start, @end, @after);
my (@line, @prev);

GetOptions(\%opts, 'h', "novo=s", "chr=s", "land=s", "rmin=i", "max=i", "density=f", "rci=f", "higher=i", "peaks=i", "parent=s", "out=s", "outf=s", "low=i", "cov=s" );

my $help = "USAGE:
	perl Postprocessing.pl -novo Novo.out –land LAND-DIR/ -chr CHR –parent P -density DEN –max MAX -higher HIGH -cov COV -rci RCI –rmin RMIN -low LOW -peaks PEAK -out Novo.chr.genome.tp\n
PARAMETERS
	novo		File with de candidate de novo SNVs
	land		Directory with the VLs (Variation landscapes)
	chr		The chromosome to be analyzed [CHR]
	parent		If the landscape corresponds to the child set [P] to 0. If the landscape corresponds to the parent set [P] to 1
	density		Regions with a density of at least [DEN] will be kept
	max/higher	Regions with at most [HIGH] CSs with a coverage higher than [MAX] will be kept
	cov		Regions for which the CSs corresponding to the child’s Sginature CSs have a coverage of at least COV will be kept 
	rci		A change higher than RCI in the RCI index will be considered significative
	rmin		A change in coverage will be considered significative only if the coverage for any CS is no lower than RMIN
	low		For the child, a CS with a coverage lower than LOW will be 	considered as a region of “almost no coverage”		
	peaks		Regions with less than PEAK significative changes in coverage 	inside the VSR will be kept
	out		Output file. This will contain all TP de novo SNVs\n
OUTPUT
	Same columns as the input file de de novo SNVs (Novo.out)\n";

die $help if $opts{h};

#########
my (%sr);
my (@line);
my ($info, $chromosome, $start, $end, $sr_name);
#########

%sr = {};
open(NOVO, "$opts{'novo'}");
## ONLY SRS WITH ONE DE NOVO SNV
while(<NOVO>){
	$info = $_;
        @line = split(/\t/, $info);
        $chromosome = $line[0];
        $start = $line[1];
        $end = $line[2];
        $sr_name = $start."-".$end;
	if (exists $sr{$sr_name}){
		$sr{$sr_name}++;
	}else{
		$sr{$sr_name} = 1;
	}
}

open(OUT, ">$opts{'out'}");
open(OUTF, ">$opts{'outf'}");


###########
my ($file_land, $pos_land, $event, $k);
my ($cs_no, $cs_higher, $pos_prev, $cov_prev, $st_parent, $end_parent, $cs_peaks, $cs_higher, $drop, $rise);
my ($pos, $cov, $cov_first, $cov_before, $cov_last, $cov_after, $pos_first, $pos_before, $pos_last, $pos_after, $max, $factor, $factor_abs);
my ($cs_sr_max, $density_cs, $factor1, $factor2, $min1, $min2, $max1, $max2);      

$file_land = $opts{'land'}.$opts{'chr'}.".land";
open(LAND, $file_land);
$pos_land = 0;

$k = $opts{'k'};

seek(NOVO,0,0);      
while(<NOVO>){
	$info = $_;
	print ($info);
	@line = split(/\t/, $info);
        $chromosome = $line[0]; 
	$start = $line[1];
	$end = $line[2];
	$event = $line[4];
	
	$sr_name = $start."-".$end;
	## ONLY SRS WITH ONE DE NOVO SNV
	if ($sr{$sr_name} > 1){
		next;
	}	

	## ONLY THE EVENTS FOR A CHROMOSOME
	if ($chromosome ne $opts{'chr'}){
		next;
	}

	$cs_no = 0;
	$cs_higher = 0;

	$cov_prev = 0;
	$pos_prev = 0;
	$st_parent = 0; $end_parent = 0;
	$cs_peaks = 0; $cs_higher = 0; $cs_no = 0;

	$drop = 0;
	$rise = 0;

	seek(LAND, $pos_land, 0);
	while(<LAND>){
		chomp($_);
		@line = split(/\t/, $_);
        	$pos = $line[0];
        	$cov = $line[1];
       		
		if ($pos == $start){
			$pos_land = tell(LAND);
		}

	
		if ($pos_prev == $start){
			$cov_first = $cov;
			$cov_before = $cov_prev;
			$pos_first = $pos;
			$pos_before = $pos_prev;
#			print "START\n$pos_first\t$cov_first\t$pos_before\t$cov_before\n";
		}
		
		if ($pos == $end){
			$cov_last = $cov_prev;
			$cov_after = $cov;
			$pos_last = $pos_prev;
			$pos_after = $pos;
#			print "END\n$pos_last\t$cov_last\t$pos_after\t$cov_after\n";
		}
	
		if ($pos >= $start && $pos <= $end){
#			print "HIGH $pos\t$cov\n";
			$cs_no++;
			if ($cov > $opts{'max'}){
				$cs_higher++;
			}
	
			$max= max($cov, $cov_prev);
			if ($max == 0){
               			$factor = 0;
        		}else{
               			$factor = ($cov - $cov_prev)/$max;
        		}

#			print "ROW\t$pos\t$cov\t$pos_prev\t$cov_prev\t$factor\n";

			if (($factor*(-1)) >= $opts{'rci'} && $cov_prev > $opts{'rmin'} && $cov > $opts{'low'} && $pos_prev != $start && $opts{'parent'} == 0){
				$cs_peaks++;
#				print "PEAK\t$factor\t$pos_prev\t$cov_prev\t$pos\t$cov\n";
			}
			if ($factor >= $opts{'rci'} && $cov > $opts{'rmin'} && $cov_prev > $opts{'low'}  && $pos != $end && $opts{'parent'} == 0){
                        	$cs_peaks++;
#				print "PEAK\t$factor\t$pos_prev\t$cov_prev\t$pos\t$cov\n";                
			}
			
			$factor_abs = abs($factor);
			if ($factor_abs >= $opts{'rci'} && $cov_prev > $opts{'rmin'} && $pos_prev == $start && $opts{'parent'} == 1){
				if ($event == ($start + 30) ){
                        		$st_parent = 1;
               	 		}
			}
			if ($factor_abs >= $opts{'rci'} && $cov > $opts{'rmin'} && $pos_prev == $event && $opts{'parent'} == 1){
				$end_parent = 1;
			}

			if (($factor*(-1)) >= $opts{'rci'} && $pos_prev == $start && $opts{'parent'} == 1){
				$drop = 1;
			}
			if ($factor >= $opts{'rci'} && $pos == $end && $opts{'parent'} == 1){
				$rise = 1;
			}

		}
		
		if ($pos > $end){
			last;
		}
		$pos_prev = $pos;
		$cov_prev = $cov;
	}

	### IF THE CS DENSITY IS LOWER THAN 'DENSITY'
	$cs_sr_max = $end - $start + 1;
	$density_cs = $cs_no/$cs_sr_max;
#	print "DENSITY\t$density_cs\n";
	if ($density_cs < $opts{'density'}){
		## IT DOES NOT APPLIES IF THE VSR HAS AT LEAST ONE INTERNAL ADJACENT CS IN EVERY SIDE
		print "$pos_first\t$pos_before\t$pos_after\t$pos_last\n";
		if ($pos_first != ($pos_before + 1) || $pos_after != ($pos_last + 1)){
			print OUTF "DENSITY\t$info";
			next;
		}
	}

	## IF THERE ARE MORE THAN 'PEAKS' DROPS OR RISES INSIDE THE VSR 
	if ($cs_peaks >= $opts{'peaks'}){
		print OUTF "PEAKS\t$info";
		next;
	}
	
	## IF THERE ARE MORE THAN 'HIGHER' CS'S WITH A COVERAHE HIGHER THAN 'MAX'
	if ($cs_higher >= $opts{'higher'}){
		print OUTF "COVERAGE\t$info";
		next;
	}

	## IF THERE IS A DROP OR RISE IN EITHER PARENT
	if ($st_parent == 1 || $end_parent == 1){
		print OUTF "PARENT-RISE-DROP\t$info";
		next;
	}

	## IF FOR ANY PARENT THE COVERAGE IN EITHER SIDE OF THE VSR IS LOWER THAN RMIN
	if ($cov_before < $opts{'cov'} || $cov_after < $opts{'cov'} ){
 		print OUTF "LOWCOV\t$info";
                next;
	}

	## THE SAME VSR IN EITHER PARENT
	if ($drop == 1 && $rise == 1){
	        print OUTF "VSR\t$info";
                next;
	}

	#FOR THE CHILD
	## THE RCI FOR ANY ANY COMBINATION OF FIRST INTERNAL CS AND SCS MUST BE HIGHER THAN RATIO
	
	$max1 = max($cov_before, $cov_last);
	$min1 = min($cov_before, $cov_last);
	$max2 = max($cov_after, $cov_first);
	$min2 = min($cov_after, $cov_first);

	$factor1 = ($max1-$min1)/$max1;
	$factor2 = ($max2-$min2)/$max2;
#	print "FACTOR\t$factor1\t$factor2\n";
	if (($factor1 < $opts{'rci'} || $factor2 < $opts{'rci'} ) && $opts{'parent'} == 0) { 
		print OUTF "NOISE\t$info";
		next;
	}

	print OUT "$info";

}

close(OUT);
close(OUTF);
close(LAND);
close(NOVO);


