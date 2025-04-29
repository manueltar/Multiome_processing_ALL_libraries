#!/bin/bash
 
eval "$(conda shell.bash hook)"
 
 
Rscripts_path=$(echo "/home/manuel.tardaguila/Scripts/R/")

MASTER_ROUTE=$1
analysis=$2

output_dir=$(echo "$MASTER_ROUTE""$analysis""/")

Log_files=$(echo "$output_dir""/""Log_files/")
 
conda activate multiome_QC_DEF

### POST_G_RPCA_integration_and_clustering

type=$(echo "POST_G_RPCA_integration_and_clustering")
outfile_POST_G_RPCA_integration_and_clustering=$(echo "$Log_files""outfile_11_""$type"".log")
touch $outfile_POST_G_RPCA_integration_and_clustering
echo -n "" > $outfile_POST_G_RPCA_integration_and_clustering
name_POST_G_RPCA_integration_and_clustering=$(echo "$type""_job")


Rscript_POST_G_RPCA_integration_and_clustering=$(echo "$Rscripts_path""465_rpca_and_clustering_post_genotyping.R")

object_post_genotyping=$(echo "/group/soranzo/manuel.tardaguila/2025_hESC_lymph_multiome/Multiome/processing_outputs/merged_unprocessed_db_db_filt_clustered_QCed_cell_annotated_rpca_integrate_rpca_integrate_clustered_only_genotyped.rds")

# --ntasks-per-node=30 --mem-per-cpu=8192

myjobid_POST_G_RPCA_integration_and_clustering=$(sbatch --job-name $name_POST_G_RPCA_integration_and_clustering --output=$outfile_POST_G_RPCA_integration_and_clustering --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=20 --mem-per-cpu=8192 --parsable --wrap="Rscript $Rscript_POST_G_RPCA_integration_and_clustering --object_post_genotyping $object_post_genotyping --type $type --out $output_dir")
myjobid_seff_POST_G_RPCA_integration_and_clustering=$(sbatch --dependency=afterany:$myjobid_POST_G_RPCA_integration_and_clustering --open-mode=append --output=$outfile_POST_G_RPCA_integration_and_clustering --job-name="seff" --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_POST_G_RPCA_integration_and_clustering >> $outfile_POST_G_RPCA_integration_and_clustering")


 
conda deactivate
