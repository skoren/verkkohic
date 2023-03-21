COV=`cat unitig-popped-unitig-normal-connected-tip.noseq.gfa |grep "^S" |sed s/LN:i://g |sed s/ll:f://g |awk '{if ($4 > 500000) { SUM+=$NF*$4; TOTAL+=$4; }} END { print SUM/TOTAL}'`
echo $COV
cat unitig-popped-unitig-normal-connected-tip.noseq.gfa |grep "^S" |sed s/LN:i://g |sed s/ll:f://g |awk -v C=$COV '{if ($NF > C*1.5 || $4 < 50000) print $2}' > unassigned
wc -l unassigned

echo -e "node\tmat\tpat\tmat:pat\tcolor" > unitig-popped-unitig-normal-connected-tip.colors.csv 
cat cluster.out |grep -A 2 Connected|grep -v Initial | grep -v Connected |awk -F "}," '{alen=split($1, a, ","); blen=split($2, b, ","); for (i = 1; i<=alen; i++) { print a[i]"\t0\t100000\t0:100000\t#8888FF"} for (i = 1; i<= blen; i++) {print b[i]"\t100000\t0\t100000:0\t#FF8888"} }'|sed 's/({//g' |sed 's/})//g' |sed s/\'//g|sed s/\ //g |sed 's/{//g' | sort |uniq |grep -w -v -f unassigned >> unitig-popped-unitig-normal-connected-tip.colors.csv
