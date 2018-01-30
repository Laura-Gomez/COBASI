#### MOVE TO PARENT DIRECTORY FOR COBASI DOWNLOADS
mv /home/COBASI/

##### OBTAIN SIGNATURE CS'S SEQUENCE
python Bin/Novo/SignatureCS_Parent.py --SignatureCSs Data/Analysis/Regions.fa --OUTPUT Data/Analysis-TRIO/Father/Regions_father.fa
python Bin/Novo/SignatureCS_Parent.py --SignatureCSs Data/Analysis/Regions.fa --OUTPUT Data/Analysis-TRIO/Mother/Regions_mother.fa

#### GET SIGNATURE READS
Bin/One_Individual/Retrieve_SignatureReads -c Data/Analysis-TRIO/Father/Regions_father.fa -t FASTA -f Data/Sequence/example.father.fasta -k 30 -i PARENT > Data/Analysis-TRIO/Father/father.fasta.reads
Bin/One_Individual/Retrieve_SignatureReads -c Data/Analysis-TRIO/Mother/Regions_mother.fa -t FASTA -f Data/Sequence/example.mother.fasta -k 30 -i PARENT > Data/Analysis-TRIO/Mother/mother.fasta.reads


#### MERGE READS
perl Bin/One_Individual/Merge_Reads.pl -dir_read Data/Analysis-TRIO/Father/ -sufix_read .fasta.reads -dir_var Data/Analysis/ -prefix_var NA -sufix_var .sr -ind PARENT -dir_out Data/Analysis-TRIO/Father/ -genome_prefix HG1_ -sufix_out .reads
perl Bin/One_Individual/Merge_Reads.pl -dir_read Data/Analysis-TRIO/Mother/ -sufix_read .fasta.reads -dir_var Data/Analysis/ -prefix_var NA -sufix_var .sr -ind PARENT -dir_out Data/Analysis-TRIO/Mother/ -genome_prefix HG2_ -sufix_out .reads


#### ALIGN READ-REFERENCE GENOME
Bin/One_Individual/Align_Read_RG -a Data/Analysis-TRIO/Father/HG1_chr11.reads -r Data/Reference/ -s .fa -k 30 -t 2 -p 1 -l 5 -c Data/Analysis/chr11.total -d FALSE -i PARENT -o Data/Analysis-TRIO/Father/chr11.total -q Data/Analysis-TRIO/Father/chr11.perbase
Bin/One_Individual/Align_Read_RG -a Data/Analysis-TRIO/Father/HG1_chrX.reads -r Data/Reference/ -s .fa -k 30 -t 2 -p 1 -l 5 -c Data/Analysis/chrX.total -d FALSE -i PARENT -o Data/Analysis-TRIO/Father/chrX.total -q Data/Analysis-TRIO/Father/chrX.perbase

Bin/One_Individual/Align_Read_RG -a Data/Analysis-TRIO/Mother/HG2_chr11.reads -r Data/Reference/ -s .fa -k 30 -t 2 -p 1 -l 5 -c Data/Analysis/chr11.total -d FALSE -i PARENT -o Data/Analysis-TRIO/Mother/chr11.total -q Data/Analysis-TRIO/Mother/chr11.perbase
Bin/One_Individual/Align_Read_RG -a Data/Analysis-TRIO/Mother/HG2_chrX.reads -r Data/Reference/ -s .fa -k 30 -t 2 -p 1 -l 5 -c Data/Analysis/chrX.total -d FALSE -i PARENT -o Data/Analysis-TRIO/Mother/chrX.total -q Data/Analysis-TRIO/Mother/chrX.perbase


#### GET SNVs GENOTYPE
Bin/One_Individual/Compute_Genotype -p Data/Analysis-TRIO/Father/chr11.perbase -t P -i PARENT -u Data/Analysis-TRIO/Father/chr11.undet.genotype > Data/Analysis-TRIO/Father/chr11.genotype
Bin/One_Individual/Compute_Genotype -p Data/Analysis-TRIO/Mother/chr11.perbase -t P -i PARENT -u Data/Analysis-TRIO/Mother/chr11.undet.genotype > Data/Analysis-TRIO/Mother/chr11.genotype

Bin/One_Individual/Compute_Genotype -p Data/Analysis-TRIO/Father/chrX.perbase -t P -i PARENT -u Data/Analysis-TRIO/Father/chrX.undet.genotype > Data/Analysis-TRIO/Father/chrX.genotype
Bin/One_Individual/Compute_Genotype -p Data/Analysis-TRIO/Mother/chrX.perbase -t P -i PARENT -u Data/Analysis-TRIO/Mother/chrX.undet.genotype > Data/Analysis-TRIO/Mother/chrX.genotype

#### OBTAIN MENDELIAN INCONGRUENT SNV's
perl Bin/Novo/Obtain_Inheritance.pl -hg3 Data/Analysis/chr11.genotype -hg1 Data/Analysis-TRIO/Father/chr11.genotype -hg2 Data/Analysis-TRIO/Mother/chr11.genotype -chr chr11 -out Data/Analysis-TRIO/Novo/

perl Bin/Novo/Obtain_Inheritance.pl -hg3 Data/Analysis/chrX.genotype -hg1 Data/Analysis-TRIO/Father/chrX.genotype -hg2 Data/Analysis-TRIO/Mother/chrX.genotype -chr chrX -out Data/Analysis-TRIO/Novo/


#### COVERAGE AND PURITY FILTER
perl Bin/Novo/Filter_Novo.pl -dir Data/Analysis-TRIO/Novo/ -sufix .mend.incongruent -cov_child 10 -total_child 2 -cov_parent 10 -total_parent 2 -min 0.25 -stat Data/Analysis-TRIO/Novo/novo.filter.stat -out Data/Analysis-TRIO/Novo/novo.cov10-10.filter

#### POSTPROCESSING
perl Bin/Extra/Postprocessing.pl -novo Data/Analysis-TRIO/Novo/novo.cov10-10.filter -land Data/Analysis-TRIO/Father/ -chr chr11 -parent 1 -density 0.6 -max 105 -higher 1 -cov 15 -rci 0.15 -rmin 10 -low 5 -peaks 1 -out Data/Analysis-TRIO/Novo/Novo.chr11.hg1.tp &
perl Bin/Extra/Postprocessing.pl -novo Data/Analysis-TRIO/Novo/novo.cov10-10.filter -land Data/Analysis-TRIO/Mother/ -chr chr11 -parent 1 -density 0.6 -max 105 -higher 1 -cov 15 -rci 0.15 -rmin 10 -low 5 -peaks 1 -out Data/Analysis-TRIO/Novo/Novo.chr11.hg2.tp &
perl Bin/Extra/Postprocessing.pl -novo Data/Analysis-TRIO/Novo/novo.cov10-10.filter -land Data/Analysis/ -chr chr11 -parent 0 -density 0.6 -max 300 -higher 1 -cov 15 -rci 0.2 -rmin 15 -low 5 -peaks 1 -out Data/Analysis-TRIO/Novo/Novo.chr11.hg3.tp &

#### MERGE POST-PROCESSING

cat Data/Analysis-TRIO/Novo/Novo.*hg3.tp  > Data/Analysis-TRIO/Novo/Novo_HG3_TP.tab 
cat Data/Analysis-TRIO/Novo/Novo.*hg1.tp  > Data/Analysis-TRIO/Novo/Novo_HG1_TP.tab 
cat Data/Analysis-TRIO/Novo/Novo.*hg2.tp  > Data/Analysis-TRIO/Novo/Novo_HG2_TP.tab 

python Bin/Extra/Compare_TP.py -f Data/Analysis-TRIO/Novo/Novo_HG1_TP.tab -m Data/Analysis-TRIO/Novo/Novo_HG2_TP.tab -c Data/Analysis-TRIO/Novo/Novo_HG3_TP.tab -o Data/Analysis-TRIO/Novo/Novo_TP.tab -p Data/Analysis-TRIO/Novo/Novo_TP.info.tab

######
###### LANDSCAPE PARENTS
######

## COUNT KMERS JELLYFISH
jellyfish count -o Data/Analysis-TRIO/Father/example_HG1 -m 30 --both-strands -s 5G -t 36 -L 2 Data/Sequence/example.father.fasta &
jellyfish count -o Data/Analysis-TRIO/Mother/example_HG2 -m 30 --both-strands -s 5G -t 36 -L 2 Data/Sequence/example.mother.fasta &

## COVERAGE ALL KMERS
kmer-cov-plot --jellyfish -s Data/Analysis-TRIO/Father/example_HG1_0 < Data/Reference/chr11.fa > Data/Analysis-TRIO/Father/chr11.coverage &
kmer-cov-plot --jellyfish -s Data/Analysis-TRIO/Father/example_HG1_0 < Data/Reference/chrX.fa > Data/Analysis-TRIO/Father/chrX.coverage &

kmer-cov-plot --jellyfish -s Data/Analysis-TRIO/Mother/example_HG2_0 < Data/Reference/chr11.fa > Data/Analysis-TRIO/Mother/chr11.coverage &
kmer-cov-plot --jellyfish -s Data/Analysis-TRIO/Mother/example_HG2_0 < Data/Reference/chrX.fa > Data/Analysis-TRIO/Mother/chrX.coverage &


## LANDSCAPE
perl Bin/One_Individual/Compute_Landscape.pl -cov Data/Analysis-TRIO/Father/chr11.coverage -unique Data/CS-Ref/ -suf_unique .str_cs_uniq_regions.tab  -out_dir Data/Analysis-TRIO/Father/ &
perl Bin/One_Individual/Compute_Landscape.pl -cov Data/Analysis-TRIO/Father/chrX.coverage -unique Data/CS-Ref/ -suf_unique .str_cs_uniq_regions.tab  -out_dir Data/Analysis-TRIO/Father/ &

perl Bin/One_Individual/Compute_Landscape.pl -cov Data/Analysis-TRIO/Mother/chr11.coverage -unique Data/CS-Ref/ -suf_unique .str_cs_uniq_regions.tab  -out_dir Data/Analysis-TRIO/Mother/ &
perl Bin/One_Individual/Compute_Landscape.pl -cov Data/Analysis-TRIO/Mother/chrX.coverage -unique Data/CS-Ref/ -suf_unique .str_cs_uniq_regions.tab  -out_dir Data/Analysis-TRIO/Mother/ &


