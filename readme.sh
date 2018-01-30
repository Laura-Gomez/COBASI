#### MOVE TO PARENT DIRECTORY FOR COBASI DOWNLOADS
mv /home/COBASI/

### COUNT KMERS
jellyfish count -o example_DB -m 30 --both-strands -s 5G -t 36 -L 2 Data/Sequence/example.fastq Data/Sequence/example.fasta 

### OBTAIN WHOLE-GENOME COVERAGE
/home/lgomezro/bin/AMOS/bin/kmer-cov-plot --jellyfish -s example_DB_0 < Data/Reference/chr11.fa > Data/Analysis/chr11.coverage
/home/lgomezro/bin/AMOS/bin/kmer-cov-plot --jellyfish -s example_DB_0 < Data/Reference/chrX.fa > Data/Analysis/chrX.coverage

### OBTAIN LANDSCAPE
perl Bin/One_Individual/Compute_Landscape.pl -cov Data/Analysis/chrX.coverage -unique Data/CS-Ref/ -suf_unique .str_cs_uniq_regions.tab  -out_dir Data/Analysis/
perl Bin/One_Individual/Compute_Landscape.pl -cov Data/Analysis/chr11.coverage -unique Data/CS-Ref/ -suf_unique .str_cs_uniq_regions.tab  -out_dir Data/Analysis/

### GET ONE-DROP REGIONS (INTERNAL PART OF VSR'S)
perl Bin/One_Individual/Get_Onedrop.pl -land chr11.land -max 540 -fac 0.2 -rmin 10 -fst 540 -nt 0 -density 0 -out Data/Analysis/chr11.onedrop
perl Bin/One_Individual/Get_Onedrop.pl -land chrX.land -max 540 -fac 0.2 -rmin 10 -fst 540 -nt 0 -density 0 -out Data/Analysis/chrX.onedrop

### GET VARIATION SIGNATURE REGIONS
perl Bin/One_Individual/Get_SR.pl -land Data/Analysis/chr11.land -var Data/Analysis/chr11.onedrop -fac 0.20 -rmin 10 -max 540 -cs_size 30 -ratio 2.0 -sr_max 100000 -out Data/Analysis/chr11.sr
perl Bin/One_Individual/Get_SR.pl -land Data/Analysis/chrX.land -var Data/Analysis/chrX.onedrop -fac 0.20 -rmin 10 -max 540 -cs_size 30 -ratio 2.0 -sr_max 100000 -out Data/Analysis/chrX.sr

### OBTAIN SIGNATURE CS'S SEQUENCE
python Bin/One_Individual/Cut_SignatureCSs.py --vsr=Data/Analysis/ --REFDIR=Data/Reference/ --sufixREF=fa --sufixVSR=.sr --prefixVSR=NA --kSIZE=30 --output=Data/Analysis/Regions.fa

### GET SIGNATURE READS
Bin/One_Individual/Retrieve_SignatureReads -c Data/Analysis/Regions.fa -t FASTA -f Data/Sequence/example.fasta -k 30 -i CHILD > Data/Analysis/example.fasta.reads
Bin/One_Individual/Retrieve_SignatureReads -c Data/Analysis/Regions.fa -t FASTQ -f Data/Sequence/example.fastq -k 30 -i CHILD > Data/Analysis/example.fastq.reads

### MERGE READS
perl Bin/One_Individual/Merge_Reads.pl -dir_read Data/Analysis/ -sufix_read .reads -dir_var Data/Analysis/ -prefix_var NA -sufix_var .sr -ind CHILD -dir_out Data/Analysis/ -prefix_out NA -sufix_out .reads

### ALIGN READ-REFERENCE GENOME
Bin/One_Individual/Align_Read_RG -a Data/Analysis/chr11.sr.reads -r Data/Reference/ -s .fa -k 30 -t 3 -p 1 -l 10 -c NA -d FALSE -i CHILD -o Data/Analysis/chr11.total -q Data/Analysis/chr11.perbase
Bin/One_Individual/Align_Read_RG -a Data/Analysis/chrX.sr.reads -r Data/Reference/ -s .fa -k 30 -t 3 -p 1 -l 10 -c NA -d FALSE -i CHILD -o Data/Analysis/chrX.total -q Data/Analysis/chrX.perbase

### COMPUTE GENOTYPE
Bin/One_Individual/Compute_Genotype -p Data/Analysis/chr11.perbase -t P -i CHILD -u Data/Analysis/chr11.undet.genotype > Data/Analysis/chr11.genotype
Bin/One_Individual/Compute_Genotype -p Data/Analysis/chrX.perbase -t P -i CHILD -u Data/Analysis/chrX.undet.genotype > Data/Analysis/chrX.genotype



