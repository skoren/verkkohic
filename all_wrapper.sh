#!/bin/bash
module load gcc/9.2.0 
module load mashmap
PATH=$PATH:/data/korens/devel/verkko-tip/lib/verkko/scripts/:/gpfs/gsfs11/users/antipovd2/devel/tmp_pstools/
PSTOOLS=/gpfs/gsfs11/users/antipovd2/devel/tmp_pstools/


HIC1=`ls $3/*R1*fastq.gz`
HIC2=`ls $3/*R2*fastq.gz`
#verkko_output_folder script_output hic-reads

#echo $HIC1 
#bwa index assembly_graph.fasta
mkdir -p $2
cat $1/6-layoutContigs/unitig-popped.layout.scfmap | awk '{if (match($1, "path")) print $2"\t"$3}' > $2/contigs_rename.map 
process_reads.py rename $2/assembly.fasta $2/contigs_rename.map $1/assembly.fasta


#homopolymer compressed assembly for mashmap
awk '/^S/{print ">"$2"\n"$3}' $1/assembly.homopolymer-compressed.gfa | fold > $2/assembly.hpc.fasta
mashmap -r $2/assembly.hpc.fasta -q $2/assembly.hpc.fasta -t 8 -f none --pi 95 -s 10000 -o $2/mashmap.out
cat $2/mashmap.out |awk '{if ($NF > 99 && $4-$3 > 500000 && $1 != $6) print $1"\t"$6}'|sort |uniq > $2/assembly.matches


#hicmapping
echo "Mapping $HIC1 and $HIC2 to $2/assembly.fasta"
echo "pstools will write BIG temporary files to current directory $PWD"
if [ ! -e $2/map_uncompressed.out ]; then
   echo "$PSTOOLS/pstools"
   $PSTOOLS/pstools hic_mapping_unitig -k19 -t60 -o $2/map_uncompressed.out $2/assembly.fasta $HIC1 $HIC2
fi
mv hic_name_connection.output $2/hic_mapping.byread.output
cp $1/assembly.homopolymer-compressed.noseq.gfa $2/assembly.hpc.gfa

echo "$( dirname -- "$0"; )"
python3 $( dirname -- "$0"; )/hicverkko.py $1 $2 
