# COBASI

### COBASI allows the accurate identification of Single Nucleotide Variants from sequencing data. It also allows to precisely identify de novo SNVs.

#### REPOSITORY CONTENT:

- **UserManual_v4.1.docx**	COBAIS User Manual
- **Bin**:			         All scripts required for the COBASI pipeline
    * Bin/One_Individual:	All scripts required for COBASI One Individual Framework
    * Bin/Novo:		         All scripts required for COBASI Family-based framework
    * Bin/CS_Database:	   All scripts required to generate the CS regions. The usage of this scipts is describe in the User Manual. No bash bash scripts to complete this step are provided.
    * Bin/Extra:		      Post and pre-processing scripts\n

- **Data**:			      All data required for the COBASI pipeline and all data generated for the COBASI pipeline
    * Data/Sequence:		Fastq or fasta file for the sequencing raw reads
    * Data/Reference:		Fasta files for Reference Genome (RG)
    * Data/CS-Ref:		   CS regions obtained from the RG
    * Data/Analysis:		Results for the COBASI One Individual Framework and the bash script
    * Data/Analysis-TRIO:	Results for the COBASI One Family-based Framework and the bash script

#### ANALYSIS
To re-run the example data:

- Decompress data in CS-Ref, Reference, and Sequence directories
- Empty the Analysis and Analysis-TRIO folders (you should keep the readme.sh files that are contained in each directory)
- For each bash script, modify the paths 
- Run the analysis in the Analysis directory
- In the Analysis-TRIO directory create the directories: Father, Mother, Novo
- Run the analysis in the Analysis-TRIO directory

Each bash script should not be run automatically, there are some steps that can be run in parallel, read user manual
Read User Manual to understand the parameters
