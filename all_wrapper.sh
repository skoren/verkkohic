#!/bin/bash
module load gcc/9.2.0 
module load mashmap
#module load snakemake
#PATH=$PATH:/data/korens/devel/verkko-tip/:/gpfs/gsfs11/users/antipovd2/devel/tmp_pstools/
#PSTOOLS=/gpfs/gsfs11/users/antipovd2/devel/tmp_pstools/
#VERKKO =/data/korens/devel/verkko-tip/
HIC1=`ls $3/*hic*/*R1*fastq.gz`
HIC2=`ls $3/*hic*/*R2*fastq.gz`
#verkko_output_folder script_output hic-reads
echo $VERKKO
echo $SLURM_JOB_ID


#echo $HIC1 
#bwa index assembly_graph.fasta

echo "---Preprocessing graph files"
mkdir -p $2
cat $1/6-layoutContigs/unitig-popped.layout.scfmap | awk '{if (match($1, "path")) print $2"\t"$3}' > $2/contigs_rename.map 
python3 $VERKKO/lib/verkko/scripts/process_reads.py rename $2/assembly.fasta $2/contigs_rename.map $1/assembly.fasta
awk '/^S/{print ">"$2"\n"$3}' $1/assembly.homopolymer-compressed.gfa | fold > $2/assembly.hpc.fasta
cp $1/assembly.homopolymer-compressed.noseq.gfa $2/assembly.hpc.noseq.gfa


echo "---Running msahmap"
#homopolymer compressed assembly for mashmap
mashmap -r $2/assembly.hpc.fasta -q $2/assembly.hpc.fasta -t 8 -f none --pi 95 -s 10000 -o $2/mashmap.out
cat $2/mashmap.out |awk '{if ($NF > 99 && $4-$3 > 500000 && $1 != $6) print $1"\t"$6}'|sort |uniq > $2/assembly.matches

echo "---Mapping reads with pstools"
#hicmapping
echo "Mapping $HIC1 and $HIC2 to $2/assembly.fasta"
echo "pstools will write BIG temporary files to current directory $PWD"
if [ ! -e $2/map_uncompressed.out ]; then
   echo "$PSTOOLS/pstools"
   $PSTOOLS/pstools hic_mapping_unitig -k19 -t60 -o $2/map_uncompressed.out $2/assembly.fasta $HIC1 $HIC2
   mv hic_name_connection.output $2/hic_mapping.byread.output
fi


if [ -n "$SLURM_JOB_ID" ] ; then
THEPATH=$(scontrol show job $SLURM_JOBID | awk -F= '/Command=/{print $2}')
echo "Detected SLURM"
else
THEPATH=$(realpath $0)
fi

echo "Path of script detected $THEPATH"
SCRIPT_DIR=$(dirname $THEPATH)
echo "---Running main script"
python3 $SCRIPT_DIR/hicverkko.py $1 $2 


echo "---Running rukki on the resulting clustering"
params=""
params="$params --init-assign out_init_ann.csv"
params="$params --refined-assign out_refined_ann.csv"
params="$params --final-assign out_final_ann.csv"
params="$params --marker-sparsity 5000"
params="$params --issue-sparsity 1000"
params="$params --try-fill-bubbles"

if [ xtrio = xtrio ]; then
   params="$params --issue-len 200000  --marker-ratio 5. --issue-ratio 3. --issue-cnt 100"
else
   params="$params --issue-len 2000000 --marker-ratio 3. --issue-ratio 2. --issue-cnt 1000"
fi

$VERKKO/lib/verkko/bin/rukki trio -g assembly.hpc.noseq.gfa -m hicverkko.colors.csv              -p rukki.paths.tsv $params
$VERKKO/lib/verkko/bin/rukki trio -g assembly.hpc.noseq.gfa -m hicverkko.colors.csv --gaf-format -p rukki.paths.gaf $params

echo "---final verkko consensus on paths"
sh $VERKKO/bin/verkko --slurm --paths $2/rukki.paths.gaf --assembly $1 -d $2/consensus/ --hifi $3/hifi/*fastq.gz --nano $3/ont/*fastq.gz
