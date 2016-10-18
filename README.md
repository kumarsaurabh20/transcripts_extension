#Partial trainscripts extension
This ruby based script is a  helper program (only for linux and mac based OS) to aid in partial transcripts extension workflow. It creates two read databases, database A (for sequence reads searcing) and database B (fetching mates of paired-end reads). To learn more about the concept behind reads databases, please go through the manuscript.

##Reads Database A
Reads database A is just a reads repository as created by blast tool. If you are using nhmmer for very sensitive searches, just merge all the reads, form fastq files, in to a single fasta formatted file.

##Reads Database B
Reads Database B is quality and length trimmed fastq files (R1 and R2). We keep the R1 and R2 file separate to effectively retrieve good quality reads using the hits obtained from BLAST/NHMMER searches. 

##Dependencies
this program is dependent on 3 programs namely, blast command line tools, seqtk, sickle and an optional nhmmer tool.The whole workflow mainly depends on Geneious software for the reads mapping step.

Keep your gzipped fastq files (RNASeq datasets: merge all samples together in to R1 and R2 files) in Data folder. (Remove the test files before running your datasets). Create a multi-fasta file with all partial transcripts from de-novo transcriptome assembly and place it in program's root folder. (You can replace the partials.fasta file with your partial transcripts file)

##Output
The program will create a new DB folder with blast/nhmmer databases. The partial transcripts fasta is first moved to DB and splitted in to individual sequence files. Next, database searches are performed and later all the files are transferred to Data folder. A new folder is created, for individual partial transcripts, with the same name as your partial transcripts header sequences. These folder have partial transcript sequence, blast/nhmmer results (just the hits name) and filtered/targeted subset of reads.

##Usage

Usage: ruby init.rb [options]

Required options:
-----------------
    -s, --transcripts <string>       incompletely assembled or partial transcripts. Accepts a fasta file with single partial transcript or a consensus sequence.
    -a, --algorithm <blast|nhmmer>   Search algorithm (blast for faster searches and hmmer for a sensitive search).
    -p, --prefix <String>            Database Prefix.

Common options:
---------------
    -k, --precheck                   Use this option to check the availability of required tools in the system
    -d, --database <database_prefix> Use this option to create blast and hmmer database
    -v, --version                    Show version
    -i, --info                       Display tool information
    -h, --help                       Displays all options


Follow the steps:

1> Check if you have all the supporting tools in your path:

`ruby init.rb --precheck`

2> Create the databases:

`ruby init.rb --database <databases_name>`

3> Create a targeted subset of reads to further use them in the mapping step in Geneious:

`ruby init.rb -s partials.fasta -a blast -p <database_name>`

NOTE: The database name should be same as in step 2.

4> Import the filtered R1 and R2 fastq/fasta files from the Data folder in to Geneious along with related partial transcript sequence.

5> Perform the mapping using Geneious mapper (Mapping using reference. Here reference sequence should be your partial transcript. In Geneious, use the sensitivity as recommended by the program and check the option "do not trim reads").

6> Once the mapping is done, export the consensus sequence of the contig and repeat steps 3-6 multiple times until a desired gene length is obtained.  

##Test Data
Sample data has been included with the program just to understand the workflow and file locations. Before running your own datasets, it is recommended to run the program using sample data and using the above steps and later remove the sample files with your datasets.
