#!/bin/bash

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

cores=$SLURM_CPUS_PER_TASK
slurm="--slurm"
if [ "x$cores" == "x" ]; then
   cores=`cat /proc/cpuinfo |grep -c processor`
   slurm=""
fi
echo "Starting with previous run $1, output will go to $2 using reads from $3"
echo "I have slurm is $slurm and cores is $cores"

mkdir -p $2

if [[ "$#" -ge 4 ]]; then
   in_mapping=$(realpath $4)
   echo "Re-using provided bam file of $in_mapping"
   ln -s $in_mapping $2/hic_to_assembly.sorted_by_read.bam
fi

if [ ! -s $2/assembly.gfa ]; then
   echo "Generating uncompressed gfa..."

   cp $1/assembly.homopolymer-compressed.noseq.gfa $2/

   # create the appropriate input for mapping
   ## create renaming map...
   cat $1/6-*/*scfmap \
       | grep utig4 \
       | awk '{print $2"\t"$NF}' \
       > $2/rename.map

   ## rename fasta from ">unassigned-0000001" to names like ">utig4-0" (etc)
   python $VERKKO/lib/verkko/scripts/fasta_combine.py rename $2/assembly.fasta $2/rename.map $1/assembly.fasta

   ## align to get homopolymer uncompressed GFA?
   $VERKKO/lib/verkko/bin/alignGFA \
       -V -e 0.30 \
       -gfa \
       -i $1/assembly.homopolymer-compressed.gfa \
       -T $2/assembly.fasta 0 \
       -t $cores \
       -o $2/assembly.gfa
fi

echo "Aligning reads from $3..."
if [ -e $3/hic ]; then
   HIC1=`ls $3/hic/*R1_001.fast[aq].gz|tr '\n' ' ' |awk '{print substr($0, 1, length($0)-1)}'`
   HIC2=`ls $3/hic/*R2_001.fast[aq].gz|tr '\n' ' ' |awk '{print substr($0, 1, length($0)-1)}'`
   echo "Hic pairs 1: '$HIC1'"
   echo "Hic pairs 2: '$HIC2'"

   if [ ! -s $2/hic_to_assembly.sorted_by_read.bam ]; then
      echo "Mapping Illumina HiC w/BWA"
      bwa index $2/assembly.fasta && bwa mem -t $cores -5 -S -P $2/assembly.fasta <(zcat $HIC1) <(zcat $HIC2) | samtools view -bh -@ $cores -q 1 - | samtools sort -n -@ $cores - -o $2/hic_to_assembly.sorted_by_read.bam
   fi
elif [ -e $3/porec ]; then
   POREC=`ls $3/porec/*.fast[aq].gz |tr '\n' ' ' |awk '{print substr($0, 1, length($0)-1)}'`
   echo "Porec file: $POREC"

   if [ ! -s $2/hic_to_assembly.sorted_by_read.bam ]; then
      echo "Mapping PoreC w/Minimap"
      minimap2 -a -x map-ont -k 17 -t $cores -K 10g -I 8g $2/assembly.fasta <(zcat $POREC) | samtools view -bh -@ $cores -q 1 - > $2/hic_to_assembly.sorted_by_read.bam
   fi
else
   echo "Error no reads found, expected hic or porec folder"
   exit -1  
fi

echo "GFASing..."
rm -rf $2/6-gfase_output

$GFASE/phase_contacts_with_monte_carlo \
-i $2/hic_to_assembly.sorted_by_read.bam \
-g $2/assembly.gfa \
-o $2/6-gfase_output \
--use_homology \
--skip_unzip \
-m 3 \
-t $cores

mkdir -p $2/6-gfase_rukki

echo -e "node\tmat\tpat\tmat:pat\tcolor" > $2/6-gfase_rukki/unitig-popped-unitig-normal-connected-tip.colors.csv
cat $2/6-gfase_output/phases.csv |awk -F "," '{if ($2 == -1) { print $1"\t0\t100000\t0:100000\t#8888FF"} else if ($2 == 1) { print $1"\t100000\t0\t100000:0\t#FF8888"}}' >> $2/6-gfase_rukki/unitig-popped-unitig-normal-connected-tip.colors.csv

echo "---Running rukki on the resulting clustering"
params=""
params="$params --init-assign $2/6-gfaserukki/out_init_ann.csv"
params="$params --refined-assign $2/6-gfase_rukki/out_refined_ann.csv"
params="$params --final-assign $2/6-gfase_rukki/out_final_ann.csv"
params="$params --marker-sparsity 5000"
params="$params --issue-sparsity 1000"
params="$params --try-fill-bubbles"
params="$params --fillable-bubble-diff 1000"
params="$params --fillable-bubble-len 500000"
params="$params --assign-tangles --tangle-allow-deadend"
params="$params --issue-ratio 1."
params="$params --solid-homozygous-cov-coeff 1.1"
params="$params --solid-ratio 1.5"
params="$params --hap-names haplotype1,haplotype2"

if [ xtrio = xtrio ]; then
   params="$params --marker-ratio 5."
else
   params="$params --marker-ratio 3."
fi

$VERKKO/lib/verkko/bin/rukki trio -g $1/assembly.homopolymer-compressed.noseq.gfa -m $2/6-gfase_rukki/unitig-popped-unitig-normal-connected-tip.colors.csv              -p $2/6-gfase_rukki/rukki.paths.tsv $params
$VERKKO/lib/verkko/bin/rukki trio -g $1/assembly.homopolymer-compressed.noseq.gfa -m $2/6-gfase_rukki/unitig-popped-unitig-normal-connected-tip.colors.csv --gaf-format -p $2/6-gfase_rukki/rukki.paths.gaf $params

sh $VERKKO/bin/verkko $slurm --screen human --paths $2/6-gfase_rukki/rukki.paths.gaf --assembly $1 -d $2/7-final_consensus/ --hifi $3/hifi/*fast*.gz --nano $3/ont/*fast*.gz
