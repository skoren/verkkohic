#!/bin/sh
set -e

#
#  Run Rukki.  Once for Bandage input, once for Consensus.
#

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

rukki trio -g unitig-popped-unitig-normal-connected-tip.noseq.gfa -m unitig-popped-unitig-normal-connected-tip.colors.csv              -p unitig-popped-unitig-normal-connected-tip.paths.tsv $params
rukki trio -g unitig-popped-unitig-normal-connected-tip.noseq.gfa -m unitig-popped-unitig-normal-connected-tip.colors.csv --gaf-format -p unitig-popped-unitig-normal-connected-tip.paths.gaf $params
