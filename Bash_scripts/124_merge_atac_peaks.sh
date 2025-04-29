#!/usr/bin/env bash
#
# =============================================================================
# Job Script
# =============================================================================
#

#SBATCH --job-name=merge_frag
#SBATCH --mail-type=ALL
#SBATCH --mail-user=manuel.tardaguila@fht.org
#SBATCH --partition=cpuq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --output=merge_frag_%j.log
#SBATCH --mem=96G


echo "Running on $(hostname)"
echo "Started: $(date)"
echo "========================"

MASTER_ROUTE=$1

module load htslib-tools

#cd $MASTER_ROUTE

for SAMPLE in MCO_01326 MCO_01327 MCO_01328 MCO_01329; do zcat $MASTER_ROUTE/${SAMPLE}/outs/atac_fragments.tsv.gz | awk -v SAMPLE=$SAMPLE \ 'BEGIN{FS=OFS="\t"} {print $1,$2,$3,SAMPLE"_"$4,$5}\'; done | sort -k1,1 -k2,2n -S 72G | bgzip -c -@ 16 > $MASTER_ROUTE/processing_outputs/merged.atac_fragments.tsv.gz 
tabix -p bed $MASTER_ROUTE/processing_outputs/merged.atac_fragments.tsv.gz



echo "========================"
echo "Completed: $(date)"
