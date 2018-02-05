################
################ COVERAGE AND PURITY FILTER
################ 
#
#
#
# Several characteristics are required for a mendelian incongruent SNV to be classified as de novo. 
# 1) A minimum coverage is required for the variant region and for the variant allele for every individual; 
# 2) The difference in coverage between the different alleles(in the child) must be lower than a required threshold; 
# 3)No parent should habe both child alleles contained in more than one total alignment. 
#
#
# COMMAND LINE.
# Filter_Novo.pl  -dir child-GENOTYPE/ -sufix mend.incongruent –cov_child 5 –total_child 5 –cov_parent 5 –total_parent –min MIN -stat  Stat.out > result.out
#
# PARAMETERS.
#	dir		Directory with mendelian incongruent genotype
#	sufix		Sufix for the mendelian incongruent genotype files 
#	cov_child	Minimum number of reads [either partial or total] in the child
#	total_child	Minimum number of total alignments in the child
#	cov_parent	Minimum number of reads [either partial or total] in either parent
#	total_parent	Minimum number of total alignments in either parent
#	min		Events with a ratio higher than MIN between the read counts for the different alleles are filtered out
#	stat		Output file with statistics of filtered SNVs
#	result.out	Output with the SNVs that PASSED the filtering criteria
#
################
################ 
################


#!/usr/bin/perl
##

use strict;
use warnings;
use Getopt::Long;
use Switch;

################
################ 
################

sub check_end_end{
	my ($end_end, $total, $al_ref, $al_all_ref)  = @_;
	my @end_al = split("/", $end_end);
	my @al = @{$al_ref};
	my @al_all = @{$al_all_ref};

	for (my $i=0; $i < scalar(@al); $i++){
		for (my $al_cont = 0; $al_cont < scalar(@al_all); $al_cont++){
			if ($al_all[$al_cont] eq $al[$i]){
				if ($end_al[$al_cont] < $total){ return (0); }
			}
		}
	}
	return (1);
}
################
################ 
################



my (%opts);

GetOptions(\%opts, 'h', "dir=s", "sufix=s", "cov_child=i", "total_child=i", "cov_parent=i", "total_parent=i", "min=f",  "stat=s", "out=s");

my $help = "USAGE:
	Filter_Novo.pl  -dir child-GENOTYPE/ -sufix mend.incongruent –cov_child 5 –total_child 5 –cov_parent 5 –total_parent –min MIN -stat  Stat.out > result.out\n
PARAMETERS:
	dir		Directory with mendelian incongruent genotype
	sufix		Sufix for the mendelian incongruent genotype files 
	cov_child	Minimum number of reads [either partial or total] in the child
	total_child	Minimum number of total alignments in the child
	cov_parent	Minimum number of reads [either partial or total] in either parent
	total_parent	Minimum number of total alignments in either parent
	min		Events with a ratio higher than MIN between the read counts for the different alleles are filtered out
	stat		Output file with statistics of filtered SNVs
	result.out	Output with the SNVs that PASSED the filtering criteria\n
OUTPUT.
	The same information as the  SUMMARY file with mendelian incongruent genotypes per chromosome [mend.incongruent]\n";

die $help if $opts{h};

open(OUT, ">$opts{'out'}");

opendir(DIR, "$opts{'dir'}");
my @novo_files = grep{/$opts{'sufix'}$/} readdir(DIR);
closedir(DIR);

my ($file, $line);
my (@info, @total_al_child, @total_al_father, @total_al_mother, @al_child, @al_mother, @al_father, @al_all_child, @al_all_mother, @al_all_father,  @end_al_mother, @end_al_father, @end_al_child);
my ($proximal, $distal, $distance, $ref, $end_end_child, $end_end_father, $end_end_mother, $total_child, $total_father, $total_mother, $al_info_child, $al_info_father, $al_info_mother, $ge_child, $ge_father, $ge_mother, $count_child, $count_father, $count_mother, $end_check_child, $end_check_father, $end_check_mother, $pure_child, $pure_father, $pure_mother, $al_specific_child, $al_specific_father, $al_specific_mother, $found_hete, $parents_nochild, $mother_allalleleschild, $father_allalleleschild);
my $pos_var;

my $cov_child = $opts{'cov_child'};
my $total_cov_child = $opts{'total_child'};

my $cov_parent = $opts{'cov_parent'};
my $total_cov_parent= $opts{'total_parent'};


my ($fail_cov, $fail_total, $fail_alleles, $fail_ratio);
$fail_cov = 0; $fail_total = 0; $fail_alleles=0; $fail_ratio=0;
foreach my $novo_file (@novo_files){
        $file = $opts{'dir'}.$novo_file;

	open(IN, $file) or die "Can't open < $file: $!";
	while(<IN>){
		$line = $_;
		chomp($line);
		@info = split(/\t/, $line);
		$proximal = $info[1]; $distal = $info[2];
		$distance = $info[2] - $info[1];
		$pos_var = $info[3];
		$ref = $info[4];
		$end_end_child = $info[6]; $end_end_father = $info[12]; $end_end_mother = $info[18];
		$total_child = $info[8]; $total_father = $info[14]; $total_mother = $info[20];
		$al_info_child = $info[5]; $al_info_father = $info[11]; $al_info_mother = $info[17];
		$ge_child = $info[9]; $ge_father = $info[15]; $ge_mother = $info[21];
		$al_specific_child = $info[10]; $al_specific_father = $info[16]; $al_specific_mother = $info[22];

		@total_al_child = split("/", $total_child);
		@total_al_father = split("/", $total_father);
		@total_al_mother = split("/", $total_mother);
			
		@al_child = split("/", $al_specific_child);
		@al_father = split("/", $al_specific_father);
		@al_mother = split("/", $al_specific_mother);

		@al_all_child = split("/", $al_info_child);
		@al_all_mother = split("/", $al_info_mother);
		@al_all_father = split("/", $al_info_father);

		@end_al_child = split("/", $end_end_child);
        @end_al_father = split("/", $end_end_father);
        @end_al_mother = split("/", $end_end_mother);
        
        my %total_read_child = ();
		my %total_read_father = ();
		my %total_read_mother = ();


		#### REMOVE N ALLELES
		for(my $i = 0; $i < scalar(@al_all_child); $i++){
			if ($al_all_child[$i] eq "N"){
				splice @al_all_child, $i, 1;
				splice @total_al_child, $i, 1;
			}
		}

		for(my $i = 0; $i < scalar(@al_all_mother); $i++){
                        if ($al_all_mother[$i] eq "N"){
                                splice @al_all_mother, $i, 1;
                                splice @total_al_mother, $i, 1;
                        }
                }

		 for(my $i = 0; $i < scalar(@al_all_father); $i++){
                        if ($al_all_father[$i] eq "N"){
                                splice @al_all_father, $i, 1;
                                splice @total_al_father, $i, 1;
                        }
                }
		
		$count_child = scalar(@total_al_child);
		$count_father = scalar(@total_al_father);
		$count_mother = scalar(@total_al_mother);
	
		## NUMBER OF ALIGNMENTS (ALL ALIGNMENTS) AT ANY GENOME SHOULD BE HIGHER THAN SOM THRESHOLD ( 15-20 )
		## THIS THRESHOLD COULD BE DIFFERENT IN THE PARENTS COMPARED TO THE CHILD
		if ($total_al_child[$count_child-1] < $cov_child || $total_al_father[$count_father-1] < $cov_parent || $total_al_mother[$count_mother-1] < $cov_parent){
			$fail_cov++;
			next;
		}
		
		## EVERY ALLELE SHOULD BE SUPPORTED BY AT LEAT total_cov END-END ALN I ANY INDIVIDUAL
		## THIS THRESHOLD COULD BE DIFFERENT IN THE PARENTS COMPARED TO THE CHILD
		
		$end_check_child = check_end_end($end_end_child, $total_cov_child, \@al_child, \@al_all_child);
		$end_check_father = check_end_end($end_end_father, $total_cov_parent, \@al_father, \@al_all_father);
		$end_check_mother = check_end_end($end_end_mother, $total_cov_parent, \@al_mother, \@al_all_mother);
		
		if ($end_check_child == 0 || $end_check_mother == 0 || $end_check_father == 0){
			$fail_total++;
			next;
		}

		#### THE CHILD ALLELE WITH THE LOWEST COUNT SHOULD BE IN MORE THAN 1/4 (0.25) READS (PARTIAL + TOTAL ALIGNMENTS)
		
		my $count = 1000;
		for (my $i = 0; $i < scalar(@al_child); $i++){
			for (my $i_all=0; $i_all < scalar(@al_all_child); $i_all++){
				if ($al_child[$i] eq $al_all_child[$i_all]){
					if ($total_al_child[$i_all] < $count){
						$count = $total_al_child[$i_all];
					}
					last;
				}
			}
		}
		my $min_ratio = $count/$total_al_child[$count_child-1];
		if ($min_ratio < $opts{'min'}){
			$fail_ratio++;
			next;
		}

		####  IF ANY PARENT HAS ALL CHILD ALLELES CONTAINED IN TOTAL ALN THAT SNV IS FILTERED
		
        for (my $i = 0; $i < scalar(@al_all_child); $i++){
        	for (my $j = 0; $j < scalar(@al_child); $j++){
            	if ($al_all_child[$i] eq $al_child[$j]){
                	$total_read_child{$al_all_child[$i]} = $end_al_child[$i];
                    print "CHILD - $al_all_child[$i] - $end_al_child[$i]\n";
                }
            }
        }

        for (my $i = 0; $i < scalar(@al_all_mother); $i++){
        	if ($end_al_mother[$i] > 1){
            	$total_read_mother{$al_all_mother[$i]} = $end_al_mother[$i];
               	print "MOTHER - $al_all_mother[$i] - $end_al_mother[$i]\n";
            }
        }

        for (my $i = 0; $i < scalar(@al_all_father); $i++){
        	if ($end_al_father[$i] > 1){
            	$total_read_father{$al_all_father[$i]} = $end_al_father[$i];
                print "FATHER - $al_all_father[$i] - $end_al_father[$i]\n";
            }
        }

        foreach my $allele (keys %total_read_child){
	        print "AL- $allele\n";
            if (!exists $total_read_mother{$allele} && !exists $total_read_father{$allele}){
				print OUT "$info[0]\t$proximal\t$distal\t$distance";
				for (my $i = 3; $i < scalar(@info); $i++){
					print OUT "\t$info[$i]";
				}
				print OUT "\n";
				last;
			}
		}
	}
#	die;
}

open(STAT, ">$opts{'stat'}");
print STAT "NOT ENOUGH COVERAGE IN AT LEAST ONE INDIVIDUAL\t$fail_cov\n";
print STAT "NOT ENOUGH TOTAL ALIGNMENTS IN AT LEAST ONE INDIVIDUAL\t$fail_total\n";
print STAT "RATIO LESS ABUNDANT ALLELE/MOST ABUNDANT ALLELE LOWER THAN THRESHOLD\t$fail_ratio\n";
close(STAT);
close(OUT);

