#!/bin/bash

eval "$(conda shell.bash hook)"
 
  
Rscripts_path=$(echo "/home/manuel.tardaguila/Scripts/R/")

MASTER_ROUTE=$1
analysis=$2


##########################################################################

output_dir=$(echo "$MASTER_ROUTE""$analysis""/")

Log_files=$(echo "$output_dir""/""Log_files/")
 
conda activate multiome_QC_DEF
 
### Export_h5ad_RNA_and_ATAC

type=$(echo "Export_h5ad_RNA_and_ATAC")
outfile_Export_h5ad_RNA_and_ATAC=$(echo "$Log_files""outfile_8_""$type"".log")
touch $outfile_Export_h5ad_RNA_and_ATAC
echo -n "" > $outfile_Export_h5ad_RNA_and_ATAC
name_Export_h5ad_RNA_and_ATAC=$(echo "$type""_job")


Rscript_Export_h5ad_RNA_and_ATAC=$(echo "$Rscripts_path""474_export_h5ad_RNA_and_ATAC_CUX1_multiome.R")


Seurat_object=$(echo "/group/soranzo/manuel.tardaguila/SC_RNA_seq/k562_multiome/NEW_object_output/merged_downstream_analysis_FOR_DANIELA.rds")

myjobid_Export_h5ad_RNA_and_ATAC=$(sbatch --job-name $name_Export_h5ad_RNA_and_ATAC --output=$outfile_Export_h5ad_RNA_and_ATAC --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=8 --mem-per-cpu=8192 --parsable --wrap="Rscript $Rscript_Export_h5ad_RNA_and_ATAC --Seurat_object $Seurat_object --type $type --out $output_dir")
myjobid_seff_Export_h5ad_RNA_and_ATAC=$(sbatch --dependency=afterany:$myjobid_Export_h5ad_RNA_and_ATAC --open-mode=append --output=$outfile_Export_h5ad_RNA_and_ATAC --job-name="seff" --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_Export_h5ad_RNA_and_ATAC >> $outfile_Export_h5ad_RNA_and_ATAC")

conda deactivate
