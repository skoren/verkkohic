#!/usr/bin/bash
PSTOOLS="/gpfs/gsfs11/users/antipovd2/devel/source_pstools/pstools"
HIC1=$1/*hic*/*R1*fastq.gz
HIC2=$1/*hic*/*R2*fastq.gz
echo "pstools will write BIG temporary files to current directory $PWD"
if [ ! -e $2/map_uncompressed.out ]; then
   echo "$PSTOOLS/pstools hic_mapping_unitig -k19 -t60 -o $2/map_uncompressed.out $2/unitigs.fasta <(zcat $HIC1) <(zcat $HIC2)"
   $PSTOOLS/pstools hic_mapping_unitig -k19 -t60 -o $2/map_uncompressed.out $2/unitigs.fasta <(zcat $HIC1) <(zcat $HIC2)
   mv hic_name_connection.output $2/hic_mapping.byread.output
fi
