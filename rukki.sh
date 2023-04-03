#!/bin/bash

set -e -o pipefail

mkdir -p $2/6-gfase_rukki

echo "---Running rukki on the resulting clustering"
params=""
params="$params --init-assign $2/6-gfase_rukki/out_init_ann.csv"
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

$VERKKO/lib/verkko/bin/rukki trio -g $2/unitigs.homopolymer-compressed.noseq.gfa -m $2/6-gfase_rukki/unitig-popped-unitig-normal-connected-tip.colors.csv              -p $2/6-gfase_rukki/rukki.paths.tsv $params
$VERKKO/lib/verkko/bin/rukki trio -g $2/unitigs.homopolymer-compressed.noseq.gfa -m $2/6-gfase_rukki/unitig-popped-unitig-normal-connected-tip.colors.csv --gaf-format -p $2/6-gfase_rukki/rukki.paths.gaf $params

sh $VERKKO/bin/verkko $4 --screen human --paths $2/6-gfase_rukki/rukki.paths.gaf --assembly $1 -d $2/7-final_consensus/ --hifi $3/hifi/*fast*.gz --nano $3/ont/*fast*.gz

if [ ! -e $2/7-final_consensus/assembly.haplotype1.fasta && -e $2/7-final_consensus/assembly.fasta ]; then
   echo "--------------------"
   echo "Finding haplotypes for label1='haplotype1' and label2='haplotype2'."
   echo ""

   touch $2/7-final_consensus/7-consensus/paths-haplotype1    #  Because if we don't, the awk below isn't guaranteed
   touch $2/7-final_consensus/7-consensus/paths-haplotype2    #  to make them, and then the extract fails.
   touch $2/7-final_consensus/7-consensus/paths-unassigned

   awk < $2/7-final_consensus/6-layoutContigs/unitig-popped.layout.scfmap \
     -v d=$2 'BEGIN {
        FS="[ \t]+"; OFS="\t";
      }
      ($1 = "path") && ($2 ~ /^haplotype1-[0-9]+$/)    { print $2 > d"/7-final_consensus/paths-haplotype1"; }
      ($1 = "path") && ($2 ~ /^haplotype2-[0-9]+$/)    { print $2 > d"/7-final_consensus/paths-haplotype2"; }
      ($1 = "path") && ($2 ~ /^unassigned-[0-9]+$/)    { print $2 > d"/7-final_consensus/paths-unassigned"; }'

   echo "--------------------"
   echo "Extracting haplotypes."
   echo ""

   if [ -e $2/7-final_consensus/paths-haplotype1 ] ; then
      $VERKKO/lib/verkko/scripts/fasta_extract.py extract $2/7-final_consensus/assembly.haplotype1.fasta $2/7-final_consensus/paths-haplotype1 $2/7-final_consensus/assembly.fasta
   fi
   if [ -e $2/7-final_consensus/paths-haplotype2 ] ; then
      $VERKKO/lib/verkko/scripts/fasta_extract.py extract $2/7-final_consensus/assembly.haplotype2.fasta $2/7-final_consensus/paths-haplotype2 $2/7-final_consensus/assembly.fasta
   fi
   if [ -e $2/7-final_consensus/paths-unassigned ] ; then
      $VERKKO/lib/verkko/scripts/fasta_extract.py extract $2/7-final_consensus/assembly.unassigned.fasta $2/7-final_consensus/paths-unassigned $2/7-final_consensus/assembly.fasta
   fi
fi

cp $2/7-final_consensus/assembly.*fasta $2/
