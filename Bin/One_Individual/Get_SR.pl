################
################ GET SR
################ 
#
# In the previous step possibly incomplete VSRs are detected. 
# In this step, additional drops will be searched at the start of the previously identitfied partial VSRs and adittional rises will be looked at the end of such VSRs.
# Finally adjacent partial VSR’s will be concatenated
# 
#       COMMAND LINE.
#  perl Get_SR.pl -land FILE.LAND –var FILE.ONEDROP -fac RCI -rmin R -max M -cs_size K –ratio RAT –sr_max SRM –out FILE.SR
#
# PARAMETERS.
#	land 		Input file. VL file [FILA.LAND]
#	var 		File obtained in the Get One-Drop Regions process [FILE.ONEDROP]
#	fac 		A Relative Coverage Index absolute value higher than RCI is considered as bona fide variation signal
#  				(set to RCI=0.3)
#	rmin 		The PrevCS and PostCS coverage must be at least R (set to R=10)
#	max 		CS’s with a count higher than M are skipped in the VSR process
#	cs_size		An additional drop or rise are searched in the K nucleotides before and after the partial VSR
#				 that has been detected in the var file
#	ratio		The ratio of the coverage between the PrevCS and PostCS should be < than RAT (recommended set RAT= 1.8)
#	sr_max		The distance between PrevCS and PostCS should be < thanSRM (recommenden set SRM=100000)
#	out 		Output file with final VSR information [FILE.SR]
#
# OUTPUT.
# 
# The output file contains final VSRs information (VSR files). This output file is composed of 14 columns:
# 1.	PrevCS position
# 2.	PrevCS count
# 3.	Next CS after PrevCS position
# 4.	Next CS after PrevCS count
# 5.	Number of nucleotides in the internal region of the VSR
# 6.	Number of inter CSs  (a VSR with a CS density of 1 will have the same value in columns 5 and 6)
# 7.	Last CS before PostCS position
# 8.	Last CS before PostCS count
# 9.	PostCS position
# 10.	PostCS count
# 11.	RCI value for the VSR start
# 12.	RCI value for the VSR end
# 13.	Median of the coverage values for InterCS’s
# 14.	 Flag indicating in which end a ladder is found (UP, DOWN, BOTH, NA)
################
################ 
################


use strict;
use Getopt::Long;
use List::Util qw[min max];

my (%opts);
my ($cs_prev, $nprev, $cs_start, $nstart, $cs_end, $nend, $cs_after, $nafter, $ni_cs, $nt); 
my ($upper_limit, $lower_limit, $upper_ladder, $lower_ladder, $coin, $read);
my ($fac_st, $fac_end, $factor, $factor_abs, $max, $mediana, $region, $cs_no_end, $no_cs);
my ($upper_limit, $lower_limit, $upper_ladder, $lower_ladder);
my (@before, @start, @end, @after);
my (@line_var, @line, @prev);
#my ($up_read, $down_read);
#my ($ratio_out);

GetOptions(\%opts, 'h', "land=s", "var=s", "fac=f", "rmin=i", "max=i", "cs_size=i", "ratio=f", "sr_max=i",  "out=s");

my $help = "USAGE
	perl Get_SR.pl -land FILE.LAND –var FILE.ONEDROP -fac RCI -rmin R -max M -cs_size K –ratio RAT –sr_max SRM –out FILE.SR\n
PARAMETERS.
	land 		Input file. VL file [FILA.LAND]
	var 		File obtained in the Get One-Drop Regions process [FILE.ONEDROP]
	fac 		A Relative Coverage Index absolute value higher than RCI is considered as bona fide variation signal (set to RCI=0.3)
	rmin 		The PrevCS and PostCS coverage must be at least R (set to R=10)
	max 		CS’s with a count higher than M are skipped in the VSR process
	cs_size		An additional drop or rise are searched in the K nucleotides before and after the partial VSR that has been detected in the var file
	ratio		The ratio of the coverage between the PrevCS and PostCS should be < than RAT (recommended set RAT= 1.8)
	sr_max	The distance between PrevCS and PostCS should be < than	SRM (recommenden set SRM=100000)
	out 		Output file with final VSR information [FILE.SR]\n
OUTPUT.
	The output file contains final VSRs information (VSR files). This output file is composed of 14 columns:
	1.	PrevCS position
	2.	PrevCS count
	3.	Next CS after PrevCS position
	4.	Next CS after PrevCS count
	5.	Number of nucleotides in the internal region of the VSR
	6.	Number of inter CSs  (a VSR with a CS density of 1 will have the same value in columns 5 and 6)
	7.	Last CS before PostCS position
	8.	Last CS before PostCS count
	9.	PostCS position
	10.	PostCS count
	11.	RCI value for the VSR start
	12.	RCI value for the VSR end
	13.	Median of the coverage values for InterCS’s
	14.	 Flag indicating in which end a ladder is found (UP, DOWN, BOTH, NA)\n";

die $help if $opts{h};

print "$opts{'land'}\n";

open(LN, "$opts{'land'}");
open(VAR, "$opts{'var'}");

open(OUT, ">$opts{'out'}");

my $temp_out = $opts{'out'}.".temp.txt";
open(TEMP, ">$temp_out");

my $pos_land = 0;
my $cs_size = $opts{'cs_size'};
my $ratio_out_thr = $opts{'ratio'};

#### OPEN ONEDROP FILE (PARTIAL VSR'S)
while (<VAR>){
	chomp($_);
    @line_var = split(/\s+/, $_);
	$cs_prev = $line_var[0]; $nprev = $line_var[1];
	$cs_start = $line_var[2]; $nstart = $line_var[3];
	$cs_end = $line_var[6]; $nend = $line_var[7];
	$cs_after = $line_var[8]; $nafter = $line_var[9];
	$no_cs = $line_var[4]; $nt = $line_var[5];
	$fac_st = $line_var[10]; $mediana = $line_var[11]; $fac_end = $line_var[12];
	
	#### THE SECOND DROP OR RISE WILL BE SEARCHED IN K NT UPSTREAM AND DOWNSTREAM OF THE PARTIAL VSR BOUNDARIES
	$upper_limit = $cs_prev - $cs_size;
	$lower_limit = $cs_after + $cs_size;

	#### OPEN LADSCAPE FILE
	$upper_ladder = 0; $lower_ladder = 0; $region = 0; $cs_no_end = 0;
	seek(LN, $pos_land, 0);
	while(<LN>){
		chomp($_);
		@line = split(/\s+/, $_);
		$coin = $line[0]; $read = $line[1];
		
		#### ALL CS'S WITH A COVERAGE ABOVE MAX THRESHOLD WILL NOT BE CONSIDERED
	 	if ($read >= $opts{'max'}){ next; }

		#### THE SERACH FOR THE NEXT SR WILL START AT THE START OF THE PREVIOUS SR
		if ($coin == $cs_prev){
			$pos_land = tell(LN);
		}
		#### AT THE LOWER LIMIT THE SEARCH WILL DE ENDED
		if ($prev[0] > $lower_limit){
			last;
		}
		
		#### A DOUBLE DROP WILL BE SEARCHED BETWEEN UPPER_LIMIT AND THE PREV-CS
		#### THE PREV-CS AND START-CS CAN BE NOT ADJACENT  
		if($coin >= $upper_limit && $coin <= $cs_prev){	
			$factor = obtain_factor($read, $prev[1]);
			$factor_abs = abs($factor);
			
			#### IF A DOUBLE DROP HAS BEEN FOUND, COUNT THE NUMBER OF CS'S IN THE REGION
			if ($region == 1){
				$no_cs++;
			}

			#### A DOUBLE DROP IS FOUND
			#### A POSSIBLE START FOR THE WHOLE VSR IS SET
        	if($factor_abs >= $opts{'fac'} ){
				if($factor < 0  && $prev[1] >= $opts{'rmin'} ){
					$upper_ladder = 1;
					@before = @prev;
					@start = @line;
					$fac_st = $factor;
					$region = 1;
					$no_cs++;
					@prev = @line;
					next;
				}
			}
			#### IF AN ADDITIONAL RISE IS FOUND AFTER THE ADDITIONAL DROP, THAT RISE IS NOT TAKEN AS THE BOUNDARY OF THE VSR
			if ($upper_ladder == 1 &&  $factor > $opts{'fac'} && $prev[1] >= $opts{'rmin'}){
				$upper_ladder = 0;
			}
		}

		#### A DOUBLE RISE WILL BE SEARCHED BETWEEN POST-CS AND THE LOWER_LIMIT
		#### THE PREV-CS AND START-CS COULD BE NOT ADJACENT  
		if ($coin > $cs_after && $prev[0] <= $lower_limit){	
			$factor = obtain_factor($read, $prev[1]);
            $factor_abs = abs($factor);
			
			#### COUNT THE NUMBER OF CS'S IN THE REGION
			$cs_no_end++;
            if($factor_abs >= $opts{'fac'}){
            	
            	#### IF THERE IS A DROP AFTER THE RISE, THE BOUNDARY OF THE VSR IS KEPT AS THE INTERNAL REGION
				if ($factor < 0  && $prev[1] >= $opts{'rmin'}){
					last;
				}
				#### THE DATA FOR THE SECOND RISE IS STORED	 
				if($factor > 0  && $read >= $opts{'rmin'} ){
					$lower_ladder = 1;
					@end = @prev;
					@after = @line;
					$no_cs = $no_cs + $cs_no_end;
					$fac_end = $factor;
					$cs_end = $end[0];
				}
			}
		}
		@prev = @line;
	}
	
	#### WRITE NEW VSR BOUNDARIES
	if ($upper_ladder == 1){
		print TEMP  "$before[0]\t$before[1]\t$start[0]\t$start[1]\t";
		$nt = $cs_end - $start[0] + 1;
	}else{
		print TEMP "$cs_prev\t$nprev\t$cs_start\t$nstart\t";
		$nt = $cs_end - $cs_start + 1;
	}

	print TEMP "$no_cs\t$nt\t";

	if ($lower_ladder == 1){
		print TEMP "$end[0]\t$end[1]\t$after[0]\t$after[1]\t";
	}else{
		print TEMP "$cs_end\t$nend\t$cs_after\t$nafter\t";
	}

	print TEMP "$fac_st\t$mediana\t$fac_end\t";
	
	if($upper_ladder == 1){
		if ($lower_ladder == 1){
			print TEMP "BOTH\n";
		}else{
			print TEMP "UP\n";
		}
	}else{
		if ($lower_ladder == 1){
        	print TEMP "DOWN\n";
     	}else{
            print TEMP "NA\n";
        }
	} 

}
close(TEMP);


my (%analyzed);
my ($prev_cs, $prev_read, $start_cs, $start_read, $cs, $nt, $end_cs, $end_read, $after_cs, $after_read);
my ($fprev, $median, $fafter, $status);
my ($pos_for, $prev_cs_for, $ratio_out);

#### CONCATENATE ADJACENT VSR'S
$pos_for = 0;

#### OPEN TEMPORAL FILE WITH COMPLETED VSR'S
open(FOR, $temp_out);
open(IN, $temp_out);
while(<IN>){
	chomp($_);
    @line = split(/\s+/, $_);
	$prev_cs = $line[0]; $prev_read = $line[1];
	$start_cs = $line[2]; $start_read = $line[3];
	$cs = $line[4]; $nt = $line[5];
	$end_cs = $line[6]; $end_read = $line[7];
	$after_cs = $line[8]; $after_read = $line[9];
	$fprev = $line[10]; $median = $line[11]; $fafter = $line[12];
	$status = $line[13];

	#### IF THAT LINE HAS ALREADY BEEN CONCATENATED TO PREVIOUS LINES, SKIP THE LINE
	if (exists $analyzed{$prev_cs}){
		next;
	}

	#### TO COMPARE, OPEN TEMPORAL FILE WITH COMPLETED VSR'S AGAIN
	seek(FOR, 0, $pos_for);
	while(<FOR>){
		chomp($_);
        @line = split(/\s+/, $_);	
		$prev_cs_for = $line[0];
		#### IF THE PREV-CS OF THIS VSR EQUAL THE AFTER-CS OF THE VSR OF THE LINE BEING COMPARED
		#### CONCATENATE VSR'S
		#### ADD THIS VSR'S TO THE LIST OF ALREADY CONCATENATED VSR'S
		#### UPDATE STATUS TO MIXED
		#### UPDATE AFTER_CS AND REPEAT THE PROCESS
		if ($prev_cs_for == $after_cs){
			$end_cs = $line[6]; $end_read = $line[7];
        	$after_cs = $line[8]; $after_read = $line[9];
			$cs = $cs + $line[4] - 1; $nt = $nt + $line[5] - 1;
			$median = int(($median + $line[11])/2);
			$fafter = $line[12]; $status = "MIXED";
			$analyzed{$prev_cs_for} = 1;
		}
		#### REPEAT THE LAST PROCESS UNTIL THE PREV_CS IS LARGER THAN THE AFTER_CS OF THE VSR BEING COMPARED
		#### CHECK VSR LENGTH AND PREV-CS/POST-CS RATIO
		#### WRITE OUTPUT
		if ($prev_cs_for > $after_cs){
			if (($after_cs - $prev_cs) < $opts{'sr_max'}){
				$ratio_out = max($prev_read, $after_read)/min($prev_read, $after_read);
        			if ($ratio_out < $ratio_out_thr){			
					print OUT "$prev_cs\t$prev_read\t$start_cs\t$start_read\t$cs\t$nt\t$end_cs\t$end_read\t$after_cs\t$after_read\t$fprev\t$median\t$fafter\t$status\n";

				}
			}
			last;
		}
		$pos_for = tell(FOR);
	}	
}

#### WRITE OUPUT FOR LAST LINE
if (!exists $analyzed{$prev_cs}){
	if (min($prev_read, $after_read) == 0){
		$ratio_out = 0;
	}else{
		$ratio_out = max($prev_read, $after_read)/min($prev_read, $after_read);
        }
	if ($ratio_out < $ratio_out_thr){	
		if (($after_cs - $prev_cs) < $opts{'sr_max'}){
			print OUT "$prev_cs\t$prev_read\t$start_cs\t$start_read\t$cs\t$nt\t$end_cs\t$end_read\t$after_cs\t$after_read\t$fprev\t$median\t$fafter\t$status\n";
		}
	}
}
close(FOR);
close(IN);
close(OUT);

my $rm = "rm ".$temp_out;
system ($rm);

sub obtain_factor { 
	my $read = $_[0]; my $prev = $_[1];
	$max= max($read, $prev);
        if ($max == 0){
                $factor = 0;
        }else{
                $factor = ($read - $prev)/$max;
        }
}
