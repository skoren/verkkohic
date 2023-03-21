#!/bin/bash
#module load gcc/9.2.0 
#module load python/3.8
#module load mashmap
#module load samtools
#module load snakemake
#PSTOOLS=/gpfs/gsfs11/users/antipovd2/devel/tmp_pstools/
#VERKKO=/data/korens/devel/verkko-tip/

#HIC1=$3/*hic*/*R1*fastq.gz
#HIC2=$3/*hic*/*R2*fastq.gz
#verkko_output_folder script_output hic-reads

set -e -o pipefail

display_usage() { 
    echo -e "\nUsage: $0 [previous_verkko_run] [output_folder] [reads folder] [aligned_sorted_bam]\n" 
    } 
# if less than two arguments supplied, display usage 
if [  $# -le 1 ] 
   then 
   display_usage
   exit 1
fi 
echo "Starting with previous run $1, output will go to $2 using reads from $3"
echo "I have slurm is $slurm and cores is $cores"

mkdir -p $2

if [[ "$#" -ge 4 ]]; then
   in_mapping=$(realpath $4)
   echo "Re-using provided bam file of $in_mapping"
   ln -sf $in_mapping $2/hic_to_assembly.sorted_by_read.bam
fi

THEPATH=$(realpath $0)
cores=$SLURM_CPUS_PER_TASK
slurm="--slurm"
if [ "x$cores" == "x" ]; then
   cores=`cat /proc/cpuinfo |grep -c processor`
   slurm=""
fi
echo "Path of script detected $THEPATH"
SCRIPT_DIR=$(dirname $THEPATH)
SCRIPT_DIR=$(realpath $SCRIPT_DIR)
SCRIPT_DIR=`echo "${SCRIPT_DIR}" | head -1`
echo "---Running main script"
echo "Directory of the script $SCRIPT_DIR "
echo "Path to script..."
echo "$SCRIPT_DIR/hicverkko.py"

#echo "---Running consensus on graph edges to get homopolymer uncompressed seqs"
#sh $VERKKO/bin/verkko  --paths $1/6-layoutContigs/consensus_paths.txt --assembly $1 -d $2/consensus_unitigs/ --hifi $3/hifi/*fast*.gz --nano $3/ont/*fast*.gz

echo "---Preprocessing graph files"
mkdir -p $2

awk '/^S/{print ">"$2"\n"$3}' $1/assembly.homopolymer-compressed.gfa | fold > $2/unitigs.hpc.fasta
cp $1/assembly.homopolymer-compressed.noseq.gfa $2/unitigs.hpc.noseq.gfa


echo "---Running mashmap"
#homopolymer compressed unitigs for mashmap
mashmap -r $2/unitigs.hpc.fasta -q $2/unitigs.hpc.fasta -t $cores -f none --pi 95 -s 10000 -o $2/mashmap.out
cat $2/mashmap.out |awk '{if ($NF > 99 && $4-$3 > 500000 && $1 != $6) print $1"\t"$6}'|sort |uniq > $2/unitigs.matches

bash $SCRIPT_DIR/uncompress.sh $1 $2 $cores

bash $SCRIPT_DIR/align.sh $1 $2 $3 $cores

BWABAM="$2/hic_to_assembly.sorted_by_read.bam"

echo "---Parsing alignments"

#samtools view -F 0x800 -q 1 $BWABAM > $2/unique.sam
samtools view -F 0x800 -q 1 $BWABAM | python3 $SCRIPT_DIR/parse_sam_pairs.py  > $2/hic_mapping.byread.output
python3 $SCRIPT_DIR/hicverkko.py $1 $2 

mkdir -p $2/6-gfase_rukki
ln -sf $(realpath $2/hicverkko.colors.tsv) $2/6-gfase_rukki/unitig-popped-unitig-normal-connected-tip.colors.csv

bash $SCRIPT_DIR/rukki.sh $1 $2 $3 $slurm
