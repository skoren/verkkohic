R1=`ls *hic*/*R1*fastq.gz`
R2=`ls *hic*/*R2*fastq.gz`
echo "Mapping $R1 and $R2 to unitig-popped-unitig-normal-connected-tip.uncompressed.fasta"
if [ ! -e map_uncompressed.out ]; then
   $BIN/pstools hic_mapping_unitig -k19 -t60 -o map_uncompressed.out unitig-popped-unitig-normal-connected-tip.uncompressed.fasta $R1 $R2
fi
