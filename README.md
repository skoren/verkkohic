This contains scripts for HiC integration in verkko.

Main script to run is all_wrapper.sh
It runs hi-c scaffolding pipeline over the preexisting previous verkko (either hifi+ont or hifi+ont+trio) run.

As a prerequirements you'll need to be installed and in $PATH:
* mashmap
* python3 with networkx package

Also, you'll need 
* pstools (for fast alignment of hi-c reads to the assembly graph). pstools release on github is too old, so we downloaded binaries here https://pstools.s3.us-east-2.amazonaws.com/pstools_1 
(recommendation by pstools' author, https://github.com/shilpagarg/DipAsm/issues/16 )
* verkko assembler (for running consensus on resulting paths and for some auxilary scripts)
 
$PSTOOLS and $VERKKO environment variables should be set to the locations of corresponding tools.

Usage:  ./all_wrapper.sh <previous_verkko_run> <output_folder> <reads_folder>

For interface simplicity script expects that all input reads are located in <reads_folder>, with hifi/ont/hic subfolders for each data type.
Reads can be in {fastq/fasta}.gz format; hic paired reads should be in separate files with names like smthR1.fastq.gz and smthR2.fastq.gz

Final assembly will be in <output_folder>/final_consensus/assembly.fasta
