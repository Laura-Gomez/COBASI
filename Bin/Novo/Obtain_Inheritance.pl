################
################ OBTAIN MENDELIAN INCONGRUENT SNV’S
################ 
#
#
# In this scrip the inheritance mode for each SNV in the child is obtained (either congruent with mendelian inheritance or incongruent). 
# To interrogate any genomic position, genotypes must be succesfuly assigned for any individual and total alignments should exist in the region.
# 
#  COMMAND LINE.
#  Obtain_Inheritance.pl –hg3 childGenotype –hg1 fatherGenotype -hg2 motherGenotype –chr  chr –out OUTDIR3
#
# PARAMETERS.
#	dir_hg3 	File with genotype information for the child
#	dir_hg1	File with genotype information for the father
#	dir_hg2	File with genotype information for the mother
#	chr		Chromosome
#	out		Output diectory
#
##### THIS PERL WILL GENERATE SEVERAL FILES:
# *.mend the one will have one label per line: CONGRUENT/INCONGRUENT (con respecto a mendelian inheritance)
# *.mend.congruent the information of the three individuals for every congruent event will be recorded
# *.mend.incongruent the information of the three individuals for every incongruent event will be recorded
# Error files with all genomic positions that:
#	Are not present in either the father or the mother SNV list
#	Have failled to succesfully assign a genotype in any individual
#	Don’t have total alignments in any individual
#
# TO SEE OUTPUT DESCRIPTION, SEE USER MANUAL
################
################ 
################



#!/usr/bin/perl
#
use strict;
use warnings;
use Getopt::Long;
use Switch;

my (%opts);

GetOptions(\%opts, 'h', "hg3=s", "hg2=s", "hg1=s", "out=s", "chr=s");

my $help = "USAGE:
Obtain_Inheritance.pl –hg3 childGenotype –hg1 fatherGenotype -hg2 motherGenotype –chr  chr –out OUTDIR\n
PARAMETERS.
	dir_hg3 	File with genotype information for the child
	dir_hg1	File with genotype information for the father
	dir_hg2	File with genotype information for the mother
	chr		Chromosome
	out		Output diectory\n
OUTPUT.
	•	A file with all SNVs with their inheritance mode: either mendelian congruent or incongruent. The information contained in this file is the same information of the genotype file for the child. In addition the last column contained  the inheritance mode [ CONGRUENT | INCONGRUENT ]
	•	A SUMMARY file with mendelian congruent genotypes per chromosome [mend.congruent]
	•	A SUMMARY file with mendelian incongruent genotypes per chromosome [mend.incongruent]
	•	Error files with all genomic positions that:
		o	Are not present in either the father or the mother SNV list
		o	Have failled to succesfully assign a genotype in any individual
		o	Don’t have total alignments in any individual\n
	FOR INFORMATION ABOUT THE CONTENT OF THE SUMMARY FILES SEE THE USER GUIDE
	
	These files will be written to:
		•	OUTDIR /{chr}.mend
		•	OUTDIR /{chr}.mend.congruent
		•	OUTDIR /{chr}.mend.incongruent\n";

die $help if $opts{h};

# Subroutine prototypes
#sub assign_alleles($;$);

my ($line_hg1, $line_hg2, $line_hg3);
my (@print_hg1, @print_hg2, @info_hg1, @info_hg2, @info_hg3, @hg1_opt, @hg2_opt, @count_hg1, @count_hg2);
my ($chr, $bef, $aft, $partial_hg3, $total_hg3, $partial_hg1, $total_hg1, $partial_hg2,$total_hg2);

my ($end_column, $ge_column, $alleles_column, $partial_column, $total_column);

$end_column = 6;
$partial_column = 7;
$total_column = 8;
$ge_column = 9;
$alleles_column = 10;

my ($cont_father, $cont_mother, $cont_both, $cont, $cont_not_found, $cont_hg3_NA, $cont_hg3_gap, $cont_st, $cont_hg12_NA, $cont_cov, $cont_ge_na);
$cont_father = 0; $cont_mother = 0; $cont_both = 0; $cont = 0; $cont_not_found = 0; $cont_hg3_NA =0; $cont_hg3_gap=0; $cont_st=0; $cont_hg12_NA=0; $cont_cov=0; $cont_ge_na=0;

$chr = $opts{'chr'};
my $file_na = $opts{'out'}.$chr."_mendelian.na";
my $file_gap = $opts{'out'}.$chr."_mendelian.gap";
#my $file_st = $opts{'out'}.$chr."_mendelian.st";
my $file_nahg = $opts{'out'}.$chr."_mendelian.nahg";
my $file_cov = $opts{'out'}.$chr."_mendelian.cov";
my $file_gena = $opts{'out'}.$chr."_mendelian.gena";
my $file_nf = $opts{'out'}.$chr."_mendelian.notfound";

open(NA, ">$file_na");
open(GAP, ">$file_gap");
#open(ST, ">$file_st");
open(NAHG, ">$file_nahg");
open(COV, ">$file_cov");
open(GENA, "$file_gena");
open(NF, ">$file_nf");

my $path_hg3 = $opts{'hg3'};
my $path_hg1 = $opts{'hg1'};
my $path_hg2 = $opts{'hg2'};

my $path_out_hg3 = $opts{'out'}.$chr.".mend";   
open(OUT, ">$path_out_hg3");

my $path_incongruent = $opts{'out'}.$chr.".mend.incongruent";
open(INC, ">$path_incongruent");
	
my $path_congruent = $opts{'out'}.$chr.".mend.congruent";
open(CON, ">$path_congruent");

print "$path_hg1\n$path_hg2\n$path_hg3\n$path_out_hg3\n";


my ($st_hg3, $ref, $hg3_possible, $hg3_alleles, $end_hg3, $ge_hg3);
my ($st_hg1, $hg1_possible, $hg1_alleles, $end_hg1, $ge_hg1);
my ($st_hg2, $hg2_possible, $hg2_alleles, $end_hg2, $ge_hg2);

my (@hg3_alleles_array, @count_hg3, @pos_array_hg3);

my $pos_hg1 = 0;
my $pos_hg2 = 0;
	
my ($cov_hg3, $hg3_al1, $hg3_al2);
my ($cov_hg1, $hg1_al1, $hg1_al2);
my ($cov_hg2, $hg2_al1, $hg2_al2);

open(HG1, $path_hg1) or die "Can't open < $path_hg1: $!";
open(HG2, $path_hg2) or die "Can't open < $path_hg2: $!";
open(HG3, $path_hg3) or die "Can't open < $path_hg3: $!";

while(<HG3>){
	my $line_hg3 = $_;
	chomp($line_hg3);
	@info_hg3 = split(/\t/, $line_hg3);
	$st_hg3 = $info_hg3[3]; $ref = $info_hg3[4]; $hg3_possible = $info_hg3[5]; 
	
	$hg3_alleles = $info_hg3[$alleles_column];
	$end_hg3 = $info_hg3[$end_column];
	$ge_hg3 = $info_hg3[$ge_column]; 
	$cont++;

	# IF THE ALLELES COULDN'T BE DETERMINED, THAT POSITION IS SKIPPED
	if ($hg3_alleles =~ m/NA/){
		$cont_hg3_NA++;
		print NA "$line_hg3\n";
		next;
	}
	
	# IF EITHER THE REFERENCE OR THE QUERY IS A GAP AT THAT POSITION, THE POSITION IS SKIPPED
	if ($ref =~ m/-/ || $hg3_possible =~ m/-/){
		$cont_hg3_gap++;
		print GAP "$line_hg3\n";
		next;
	}
		
	@hg3_alleles_array = split("/", $hg3_alleles);
	@count_hg3 = split("/", $end_hg3);
	@pos_array_hg3 = split("/", $hg3_possible);

	($cov_hg3, $hg3_al1, $hg3_al2) = assign_alleles($ge_hg3, \@hg3_alleles_array, \@count_hg3, \@pos_array_hg3);
			
	my (@hg1_alleles_array, @count_hg1, @pos_array_hg1);
	my (@hg2_alleles_array, @count_hg2, @pos_array_hg2);
	my ($A_1, $B_1, $A_2, $B_2);
	my %gen_possible;

	my $found_hg1 = 0;
	my $found_hg2 = 0;

	seek(HG1, $pos_hg1, 0);
	while(<HG1>){
		$line_hg1 = $_;
	    chomp($line_hg1);
     	@info_hg1 = split(/\t/, $line_hg1);
     	$st_hg1 = $info_hg1[3];
                
		if ($st_hg1 > $st_hg3){
			last;
		}
		if ($st_hg1 == $st_hg3){			
			$hg1_possible = $info_hg1[5]; @pos_array_hg1 = split('/', $hg1_possible); 
			$hg1_alleles = $info_hg1[$alleles_column]; @hg1_alleles_array = split('/', $hg1_alleles);
			$end_hg1 = $info_hg1[$end_column]; @count_hg1 = split('/', $end_hg1);
			$ge_hg1 = $info_hg1[$ge_column];
			@print_hg1 = @info_hg1;
			$found_hg1 = 1;
			$pos_hg1 = tell(HG1);
			last;
		}
	}

	seek(HG2, $pos_hg2, 0);
    while(<HG2>){
    	$line_hg2 = $_;
        chomp($line_hg2);
        @info_hg2 = split(/\t/, $line_hg2);
        $st_hg2 = $info_hg2[3]; 
                       
       	 if ($st_hg2 > $st_hg3){
            last;
        }
       	 if ($st_hg2 == $st_hg3){                 
      		$hg2_possible = $info_hg2[5]; @pos_array_hg2 = split('/', $hg2_possible); 
			$hg2_alleles = $info_hg2[$alleles_column]; @hg2_alleles_array = split('/', $hg2_alleles);
			$end_hg2 = $info_hg2[$end_column]; @count_hg2 = split('/', $end_hg2);
			$ge_hg2 = $info_hg2[$ge_column];
			@print_hg2 = @info_hg2;
	        $found_hg2 = 1; 
			$pos_hg2 = tell(HG2);
		    last;
        }
    }
	
	# IF THE POSITION IS NOT REPORTED FOR EITHER THE FATHER O THE MOTHER IN THE SNV INFO FILES
	if ($found_hg1 == 0){
		$cont_father++;
	}
	if ($found_hg2 == 0){
		$cont_mother++;
	}
	if ($found_hg1 == 0 && $found_hg2==0){
		$cont_both++;
	}

	if ($found_hg1 == 0 || $found_hg2 == 0){
		$cont_not_found++;
		print NF "HG3\t$line_hg3\nHG1\t$line_hg1\nHG2\t$line_hg2\n";
		next;
	}

	# IF THE POSITION IS REPORTED BUT NO ALLELE ASSIGNATION WAS POSSIBLE FOR EITHER THE FATHER O THE MOTHER IN THE SNV INFO FILES
	if ($hg1_alleles =~ m/NA/ || $hg2_alleles =~ m/NA/){
		$cont_hg12_NA++;
		print NAHG "$line_hg3\t$hg1_alleles\t$hg2_alleles\n";
		next;
	}

  	($cov_hg1, $A_1, $B_1) = assign_alleles($ge_hg1, \@hg1_alleles_array, \@count_hg1, \@pos_array_hg1);
	($cov_hg2, $A_2, $B_2) = assign_alleles($ge_hg2, \@hg2_alleles_array, \@count_hg2, \@pos_array_hg2);

	$gen_possible{$A_1}{$A_2} = 1;
	$gen_possible{$A_2}{$A_1} = 1;
	$gen_possible{$A_1}{$B_2} = 1;
	$gen_possible{$B_2}{$A_1} = 1;
	$gen_possible{$B_1}{$A_2} = 1;
	$gen_possible{$A_2}{$B_1} = 1;
	$gen_possible{$B_1}{$B_2} = 1;
	$gen_possible{$B_2}{$B_1} = 1;

	# IF THERE IS NO TOTAL ALIGNMENTS TO COVER SOME OF THE VARIABLE ALLELES, SO SKIP THAT EVENT
	if ($cov_hg1 == 0 || $cov_hg2 == 0 || $cov_hg3 == 0){
		$cont_cov++;
		print COV "HG3\t$line_hg3\nHG1\t$line_hg1\nHG2\t$line_hg2\n$cov_hg3\t$cov_hg1\t$cov_hg2\n";
		next;
	}

	# IF NO GENOTYPE ASSIGNATION WAS POSSIBLE FOR EITHER INDIVIDUAL

	if ( ($ge_hg3 ne 'NA') && ($ge_hg1 ne 'NA') && ($ge_hg2 ne 'NA') ){
		$chr = $info_hg3[0]; $bef = $info_hg3[1]; $aft = $info_hg3[2];
		$partial_hg3 = $info_hg3[$partial_column]; $total_hg3 =$info_hg3[$total_column];
		$partial_hg1 = $print_hg1[$partial_column]; $total_hg1 = $print_hg1[$total_column];
		$partial_hg2 = $print_hg2[$partial_column]; $total_hg2 = $print_hg2[$total_column];

		if (exists $gen_possible{$hg3_al1}{$hg3_al2} || exists $gen_possible{$hg3_al2}{$hg3_al1}){
			print OUT "$line_hg3\tCONGRUENT\n";
			print CON "$chr\t$bef\t$aft\t$st_hg3\t$ref\t$hg3_possible\t$end_hg3\t$partial_hg3\t$total_hg3\t$ge_hg3\t$hg3_al1/$hg3_al2\t$hg1_possible\t$end_hg1\t$partial_hg1\t$total_hg1\t$ge_hg1\t$A_1/$B_1\t$hg2_possible\t$end_hg2\t$partial_hg2\t$total_hg2\t$ge_hg2\t$A_2/$B_2\n";
				
		}else{
			print OUT "$line_hg3\tINCONGRUENT\n";
			print INC "$chr\t$bef\t$aft\t$st_hg3\t$ref\t$hg3_possible\t$end_hg3\t$partial_hg3\t$total_hg3\t$ge_hg3\t$hg3_al1/$hg3_al2\t$hg1_possible\t$end_hg1\t$partial_hg1\t$total_hg1\t$ge_hg1\t$A_1/$B_1\t$hg2_possible\t$end_hg2\t$partial_hg2\t$total_hg2\t$ge_hg2\t$A_2/$B_2\n";			
		}
	}else{
		$cont_ge_na++;
		print GENA "$line_hg3\t$ge_hg3\t$ge_hg1\t$ge_hg2\n";
	}
}

close(HG1);
close(HG2);
close(HG3);
close(OUT);
close(INC);
close(CON);
		

print "$cont\n";
print "NA HG3\t$cont_hg3_NA\n";
print "GAP HG3\t$cont_hg3_gap\n";
#print "DIF ST\t$cont_st\n";
print "NA HG1|HG2\t$cont_hg12_NA\n";
print "COV\t$cont_cov\n";
print "NA GE\t$cont_ge_na\n";
print "NOT FOUND\t$cont_not_found\n";
print "NOT FOUND FATHER\t$cont_father\n";
print "NOT FOUND MOTHER\t$cont_mother\n";
print "NOT FOUND BOTH\t$cont_both\n";

close(NA);
close(GAP);
#close(ST);
close(NAHG);
close(COV);
close(GENA);
close(NF);

sub assign_alleles{
	my ($ge, $alleles_s, $count_s, $pos_alleles_s) = @_;
        my @alleles = @{$alleles_s};
        my @count = @{$count_s};
	my @pos_alleles = @{$pos_alleles_s};
	my ($al, $al1, $al2, $cov);
	$cov = 1;

	foreach $al (@alleles){
		for (my $i=0; $i < scalar(@pos_alleles); $i++){
				if ($pos_alleles[$i] eq $al){
					if ($count[$i] == 0 ){
						$cov = 0;
					}
				}
		}
	}

	if ($ge == 1 || $ge == 4){
		$al1 = $alleles[0];
		$al2 = $alleles[1];
	}else{
		$al1 = $alleles[0];
		$al2 = $alleles[0];
	}
	
	return($cov, $al1, $al2);
}


