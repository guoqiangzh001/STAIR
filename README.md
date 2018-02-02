# STAIR
DESCRIPTION
----------------------------
STAIR is a computational pipeline identifying and validating the splicing isoforms collapsed from the full-length transcripts (such as produced by PacBio RSII or Sequel System) after being mapped against a reference genome. It classifies the Circular Consensus Sequences (CCSs) into Full-Length Non-Chimeric (FLNC) and Non-Full Length (NFL) group, applies the high confidence short reads to error-correct the FLNC CCSs and validate the collapsed splicing isoforms by checking their chained splicing junctions.

REQUIREMENTS
------------
* python 3.X

* hmmer-3.1 or higer.

* proovread, a hybrid correction pipeline for SMRT reads through iterative short read consensus.

* gmap, a genomic mapping and alignment program for mRNA sequences.

* STAR, an ultrafast universal RNA-seq aligner.

* pbtranscript-tofu, a scipt (Pacific Biosciences) classifying CCSs into FLNC and NFL groups and collapsing the FLNC CCSs into splicing isoforms.

GETTING STARTED
---------------
* Downloading STAIR pipeline.

* Adding the environmental variables of proovread, gmap and STAR to current environmental environmental variable, for example:

    export PERL5LIB=$proovread_Path/lib

    export PATH=$blast_Path/bin:$samtools_Path/bin:$proovread_Path/bin:$STAR_Path/bin/Linux_x86_64_static:$gmap_Path/bin:$hmmer_Path/binaries:$PATH

* Activating the VENV_TOFU virtual environment:

    source $VENV_TOFU_Path/bin/activate

* Running:

    You can take a look at the command line parameters using the following command line:

    python STAIR_test_full_version.py -h

    Identifying and validating: to error-correct FLNC CCSs with short reads, identify splicing isoforms, detect splicing junctions and validate splicing isoforms, run STAIR as follows:

    python STAIR.py --input ccs.fasta -read1 shortReads_1.fastq -read2 shortReads_2.fastq -subreads merged.subreads.fasta -corr Y -classify Y -ref reference_genome.fa -o Output #in the case of the pair-end reads.

    python STAIR.py --input flnc_ccs.fasta -read1 shortReads_1.fastq -read2 shortReads_2.fastq -corr Y -ref reference_genome.fa -o Output #in the case of FLNC CCSs without error-correction.

    python STAIR.py --input flnc_ccs.fasta -read1 shortReads_1.fastq -read2 shortReads_2.fastq -stardir genome_directory_STAR -gmapdir genome_directory_gmap -gmapdb genome_database_gmap -corr Y -ref reference_genome.fa -o Output #in the case of the genome database being provided.

    python STAIR.py --input flnc_ccs.fasta -reads reads.fastq -stardir genome_directory_STAR -gmapdir genome_directory_gmap -gmapdb genome_database_gmap -corr Y -ref reference_genome.fa -o Output ##in the case of the of single-end reads.

