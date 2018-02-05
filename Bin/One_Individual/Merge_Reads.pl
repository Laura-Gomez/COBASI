################
################ MERGE SIGNATURE READS INFORMATION FILES
################ 
#
# If  the sequencing experiments generates multiple FASTQ read files,
# the script that gets the Signature Reads will be run multiple times, 
# and the hits for every SignatureCS will be distributed over multiple files. 
# This script will concatenate such information.
#
#                     COMMAND LINE.
# perl Merge_Reads.pl -dir_read ALNDIR/ -sufix_read sufixREAD -dir_var SR/ -prefix_var prefixSR -sufix_var sufixSR  -ind Individual -dir_out OUTDIR/ -genome_prefix GenomePrefix -sufix_out sufixREADALN &
# PARAMETERS.
#	dir_read 	The directory where the READALN-UNORDER files were written
#	sufix_read 	Sufix of the READALN-UNORDER files
#	dir_var 	The directory where the VSR files were written.
#	prefix_var	VSR file name must be: {prefixSR}chrName{sufixSR}.
#				If there is no prefix set to NA
#	sufix_var 	Sufix of the VSR files. All the files from SR/ ending in sufixSR will be analyzed.
#	ind			Set to CHILD for one individual SNV discovery
#				Set to PARENT for SNV discovery in the parents
#	dir_out 	Output directory.
#	genome_prefix Prefix for the output files. Set to NA if no prefix is desired.
#	sufix_out	Sufix for the output files
#
# OUTPUT.
# One output file per chromosome which contains information about the alignment between the read and the respective Signature CS’s (READALN files).
# The information contained in these files is the same information contained in the READALN-UNORDER files.
#
################
################ 
################


use Getopt::Long;

my %opts;
my @pos;
GetOptions(\%opts, 'h', "dir_read=s", "sufix_read=s", "dir_var=s", "prefix_var=s", "sufix_var=s", "ind=s", "dir_out=s", "prefix_out=s", "sufix_out=s");


my $help = "USAGE:
	perl Merge_Reads.pl -dir_read ALNDIR/ -sufix_read sufixREAD -dir_var SR/ -prefix_var prefixSR -sufix_var sufixSR  -ind CHILD -dir_out OUTDIR/ -prefix_out GenomePrefix -sufix_out sufixREADALN &\n
PARAMETERS:
	dir_read 	The directory where the READALN-UNORDER files were written
	sufix_read 	Sufix of the READALN-UNORDER files
	dir_var 	The directory where the VSR files were written.
	prefix_var	VSR file name must be: {prefixSR}chrName{sufixSR}.If there is no prefix set to NA
	sufix_var 	Sufix of the VSR files. All the files from SR/ ending in sufixSR will be analyzed.
	ind		Set to CHILD for one individual SNV discovery
	dir_out 	Output directory.
	prefix_out 	Prefix for the output files. Set to NA if no prefix is desired.
	sufix_out	Sufix for the output files [READALN files]\n
OUTPUT.
	One output file per chromosome.
	Each output file contains ordered information about the alignment between the read and the respective Signature CS’s (READALN files).
	The information contained in these files is the same unordered information contained in the READALN-UNORDER files.\n";

die $help if $opts{h};

#### IF THE SNV IDENTIFICATION IS FOR THE CHILD
if ($opts{'ind'} == "CHILD"){
	$st_index = 2;
	$end_index = 3;
}else{
#### IF THE SNV IDENTIFICATION IS FOR THE PARENTS
	$st_index = 4;
	$end_index = 5;
}


#### STORE ALL NAMES FOR THE VSR FILES
opendir(VAR, "$opts{'dir_var'}");
my @var_files = grep{/$opts{'sufix_var'}$/} readdir(VAR);
closedir(VAR);

#### STORE ALL NAMES FOR THE FILES WITH THE SIGNATURE READ INFORMATION
opendir(DIR, "$opts{'dir_read'}");
my @read_files = grep{/$opts{'sufix_read'}$/} readdir(DIR);
closedir(DIR);

#### READS IN SIGNATURE READS FILES ARE ORDERED BY CHR
#### SORT VSR FILES TO READ THEM IN THE SAME ORDER
my @sort_var_files = (sort { $a cmp $b; } @var_files);

#### ALL READ FILES WILL BE READ AT THE SAMEN TIME
#### OPEN MULTIPLE FILE HANDLES
for ($i = 0; $i < scalar(@read_files); $i++){
	$path_read = $opts{'dir_read'}.$read_files[$i];
	$handle = "R".$i;
	open($handle, $path_read);
	$read_files[$i] = $path_read;
	$pos[$i] = 0;
}

#### FOR EVERY VSR, THE SIGNATURE READ INFORMATION WILL BE SEARCHED IN ALL SIGNATURE READS FILES
foreach $var (@sort_var_files){
	#### RETRIEVING CHR INFORMATION
	#### VSR FILE NAME MUST BE: {prefixSR}chrName{sufixSR}.
	if ($opts{'prefix_var'} eq 'NA'){
		$var =~ m/^(.+)$opts{'sufix_var'}/;
		$chr = $1;
	}else{
		$var =~ m/$opts{'prefix_var'}(.+)$opts{'sufix_var'}/;
		$chr = $1;
	}
	print "$chr\n";
	
	#### OPEN VSR FILE
	#### OPEN OUT FILE
	$file_var = $opts{'dir_var'}. $var;
	if ($opts{'genome_prefix'} eq "NA"){
		$out  = $opts{'dir_out'}.$var.$opts{'sufix_out'};
	}else{
		$out  = $opts{'dir_out'}.$opts{'genome_prefix'}.$chr.$opts{'sufix_out'};
	}
	open(EV, "$file_var");
	open(OUT, ">$out");

	while(<EV>){
		@ev_info = split(/\t/, $_);
		$st_ev = $ev_info[0];
		$end_ev = $ev_info[8];
		
		#### FOR EVERY VSR, THE SIGNATURE READ INFORMATION WILL BE SEARCHED IN ALL SIGNATURE READS FILES
		#### FIRST ALL INFORMATION FOR THE PRE-CS WILL BE FOUND
		for ($i = 0; $i < scalar(@read_files); $i++){
		    #### THE Ri SIGNATURE READ FILE WILL BE OPEN FROM THE LAST READING POINT
			$handle = "R".$i;
		    seek($handle, $pos[$i], 0);
			while(<$handle>){
				$read_line = $_;
				@info_read = split(/\t/, $read_line);
				$chr_read = $info_read[0]; $ev_read = $info_read[$st_index];
				#### THE READING WILL BE ENDED WHEN CHR CHANGES
				#### OR WHEN CS-POS CHANGES
				if ($chr_read ne $chr || $ev_read > $st_ev){
					last;
				}
				print OUT "$read_line";
				$pos[$i] = tell($handle);
			}
		}
		
		#### SECOND ALL INFORMATION FOR THE POST-CS WILL BE FOUND
		for ($i = 0; $i < scalar(@read_files); $i++){
			$handle = "R".$i;
		    seek($handle, $pos[$i], 0);
			while(<$handle>){
				$read_line = $_;
				@info_read = split(/\t/, $read_line);
				$chr_read = $info_read[0]; $ev_read = $info_read[$end_index];
				#### THE READING WILL BE ENDED WHEN CHR CHANGES
				#### OR WHEN CS-POS CHANGES
				if ($chr_read ne $chr || $ev_read > $end_ev){
					last;
				}
				print OUT "$read_line";
				$pos[$i] = tell($handle);
			}
		}
	}
	close(EV);
	close(OUT);
}

#### CLOSING ALL FILES	
for ($i = 0; $i < scalar(@read_files); $i++){
	close($handle);
}

