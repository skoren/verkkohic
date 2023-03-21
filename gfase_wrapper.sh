#!/bin/bash
set -e -o pipefail

display_usage() {
    echo -e "\nUsage: $0 [previous_verkko_run] [output_folder] [reads folder]\n"
    }
# if less than two arguments supplied, display usage
if [  $# -le 1 ]
   then
   display_usage
   exit 1
fi

s=`ps $$|grep -c bash`
if [ $s -le 0 ]; then
   echo "Error: please run this script under bash interpreter not sh or others"
   display_usage
   exit 1
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
echo "Starting with previous run $1, output will go to $2 using reads from $3"
echo "I have slurm is $slurm and cores is $cores"

mkdir -p $2

if [[ "$#" -ge 4 ]]; then
   in_mapping=$(realpath $4)
   echo "Re-using provided bam file of $in_mapping"
   ln -sf $in_mapping $2/hic_to_assembly.sorted_by_read.bam
fi

bash $SCRIPT_DIR/uncompress.sh $1 $2 $cores

bash $SCRIPT_DIR/align.sh $1 $2 $3 $cores

echo "GFASing..."
if [ ! -e $2/6-gfase_output/phases.csv ]; then
   rm -rf $2/6-gfase_output

   $GFASE/phase_contacts_with_monte_carlo \
      -i $2/hic_to_assembly.sorted_by_read.bam \
      -g $2/unitigs.gfa \
      -o $2/6-gfase_output \
      --use_homology \
      --skip_unzip \
      -m 3 \
      -t $cores
fi

mkdir -p $2/6-gfase_rukki

echo -e "node\tmat\tpat\tmat:pat\tcolor" > $2/6-gfase_rukki/unitig-popped-unitig-normal-connected-tip.colors.csv
cat $2/6-gfase_output/phases.csv |awk -F "," '{if ($2 == -1) { print $1"\t0\t100000\t0:100000\t#8888FF"} else if ($2 == 1) { print $1"\t100000\t0\t100000:0\t#FF8888"}}' >> $2/6-gfase_rukki/unitig-popped-unitig-normal-connected-tip.colors.csv

bash $SCRIPT_DIR/rukki.sh $1 $2 $3 $slurm
