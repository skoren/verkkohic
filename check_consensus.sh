#!/bin/bash

#opts=$1
module load gcc/9.2.0
module load python/3.7
module load snakemake/7.3.7

#sh /data/korens/devel/verkko-tip/bin/verkko --slurm -paths $1 -assembly $2 
sh /data/korens/devel/verkko-tip/bin/verkko --slurm --paths ../../data/HG00544/unitig-popped-unitig-normal-connected-tip.paths.gaf --assembly /data/Phillippy/projects/verkko/hg00544/beta_asm_faster_rukki/ -d ../../data/HG00544/verkko_rerun_test3/ --hifi /data/Phillippy/projects/verkko/hg00544/hifi/*fastq.gz
#$opts --ovb-run 8 32 96 --lay-run 1 250 48 --slurm -d beta_asm_faster_rukki --hifi `pwd`/hifi/*.fastq.gz --nano `pwd`/ont/*.fast[aq].gz --hap-kmers beta_asm/5-untip/maternal.k30.hapmer.meryl beta_asm/5-untip/paternal.k30.hapmer.meryl trio
