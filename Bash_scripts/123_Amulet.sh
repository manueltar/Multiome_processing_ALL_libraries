#!/bin/bash

eval "$(conda shell.bash hook)"
 

Rscripts_path=$(echo "/home/manuel.tardaguila/Scripts/R/")

MASTER_ROUTE=$1
analysis=$2

output_dir=$(echo "$MASTER_ROUTE""$analysis""/")

Log_files=$(echo "$output_dir""/""Log_files/")
 
conda activate multiome_QC_DEF

#sample_array=$(echo 'MCO_01326,MCO_01327,MCO_01328,MCO_01329')
sample_array=$(echo 'MCO_01330,MCO_01331,MCO_01332,MCO_01333')

a=($(echo "$sample_array" | tr "," '\n'))

declare -a arr

for i  in "${a[@]}"
do

    sample_array_sel=$i
    echo "$sample_array_sel"


    ### Amulet_run

    type=$(echo "$sample_array_sel""_""Amulet_run")
    outfile_Amulet_run=$(echo "$Log_files""outfile_4_""$type"".log")
    touch $outfile_Amulet_run
    echo -n "" > $outfile_Amulet_run
    name_Amulet_run=$(echo "$type""_job")

 
    Rscript_Amulet_run=$(echo "$Rscripts_path""404_Amulet.R")

    sample_name=$sample_array_sel
    frag_file=$(echo "$MASTER_ROUTE""$sample_array_sel""/""outs""/""atac_fragments.tsv.gz")
    barcode_file=$(echo "$output_dir""$sample_array_sel""/""intermediate/""keep_barcodes_step1.txt")
    repeats=$(echo "/group/soranzo/paola.benaglio/references/blacklist_repeats_segdups_rmsk_hg38.bed")

    echo "$frag_file"
    echo "$barcode_file"
    echo "$repeats"

    myjobid_Amulet_run=$(sbatch --job-name $name_Amulet_run --output=$outfile_Amulet_run --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=8 --mem-per-cpu=8192 --parsable --wrap="Rscript $Rscript_Amulet_run --sample_name $sample_name --frag_file $frag_file --barcode_file $barcode_file --repeats $repeats --type $type --out $output_dir")
    myjobid_seff_Amulet_run=$(sbatch --dependency=afterany:$myjobid_Amulet_run --open-mode=append --output=$outfile_Amulet_run --job-name="seff" --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_Amulet_run >> $outfile_Amulet_run")


done

conda deactivate
