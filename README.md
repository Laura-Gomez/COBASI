# COBASI

### COBASI allows the accurate identification of *de novo* Single Nucleotide Variants from sequencing data. 

#### REPOSITORY CONTENT:

- **UserManual_v5.1.pdf**	COBASI User Manual
- **readme.sh** Bash script to run the One-Individual discovery pipeline
- **readme-trio.sh** Bash script to run the De novo discovery pipeline (Family-Based Framework)
- **Bin**:			         All scripts required for the COBASI pipeline
    * Bin/One_Individual:	All scripts required for COBASI One Individual Framework
    * Bin/Novo:		         All scripts required for COBASI Family-based framework
    * Bin/CS_Database:	   All scripts required to generate the CS regions (see the the User Manual, no bash script to complete this step is provided)
    * Bin/Extra:		      Post and pre-processing scripts\n

- **Data**:			      All data required for the COBASI pipeline and all data generated for the COBASI pipeline
    * Data/Sequence:		Fastq or fasta file for the sequencing raw reads
    * Data/Reference:		The reference sequences can be downloaded from: http://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/
    * Data/CS-Ref:		   CS regions obtained from the RG
    * Data/Analysis:		Results for the COBASI One Individual Framework and the bash script
    * Data/Analysis-TRIO:	Results for the COBASI One Family-based Framework and the bash script

#### ANALYSIS
To re-run the example data:

- Decompress data in CS-Ref and Sequence directories
- Download reference data (one reference chromosome per fasta file)
- Empty the Analysis and Analysis-TRIO folders
- For each bash script, modify the paths 
- Run the analysis in the Analysis directory (follow the instructions in readme.sh)
- In the Analysis-TRIO directory create the directories: Father, Mother, Novo
- Run the analysis in the Analysis-TRIO directory (follow the instructions in readme-trio.sh)

Each bash script should not be run automatically, there are some steps that can be run in parallel, read user manual
Read User Manual to understand the parameters for every script
Some intermediary files were not uploaded into the repository because they are too big (coverage and landscape files)

