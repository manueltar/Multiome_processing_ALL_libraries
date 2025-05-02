#!/bin/bash

eval "$(conda shell.bash hook)"
  

Rscripts_path=$(echo "/home/manuel.tardaguila/Scripts/R/")

MASTER_ROUTE=$1
analysis=$2
annotation_file=$3
selected_annotations=$(echo "Dorothea_ABCD_GATA6_targets,Dorothea_ABCD_BCL11A_targets,Dorothea_AB_GATA2_targets")
selected_annotations_other=$(echo "G2M,In-house")

output_dir=$(echo "$MASTER_ROUTE""$analysis""/")

Log_files=$(echo "$output_dir""/""Log_files/")


Diff_array=$(echo 'Diff_lymph,Diff_MK')

a=($(echo "$Diff_array" | tr "," '\n'))


declare -a array_2_length

array_2_length=${#a[@]}

for (( i=0; i<${array_2_length}; i=i+1 ));
do

    Diff_array_sel=${a[$i]}
    echo "$Diff_array_sel"

    conda activate /home/manuel.tardaguila/conda_envs/multiome_NEW_downstream_analysis

    type=$(echo "$Diff_array_sel""_""$analysis""_""bespoke_heatmap_function")
    outfile_bespoke_heatmap_function=$(echo "$Log_files""outfile_6_""$type"".out")
    touch $outfile_bespoke_heatmap_function
    echo -n "" > $outfile_bespoke_heatmap_function
    name_bespoke_heatmap_function=$(echo "$type""_job")
    seff_name=$(echo "seff""_""$type")

    Rscript_bespoke_heatmap_function=$(echo "$Rscripts_path""471_bespoke_DE_heatmap_ALL_DIFFs.R")

    DE_results=$(echo "$output_dir""DE_results_""$Diff_array_sel"".rds")
    normalised_counts=$(echo "$output_dir""norcounts_FINAL_""$Diff_array_sel"".rds")
    selected_clone_lines=$(echo "wt_1,wt_2,wt_3,rs62237617_1,rs62237617_2,rs62237617_3,DNMT3A_1,DNMT3A_2,DNMT3A_3,rs62237617_DNMT3A_1,rs62237617_DNMT3A_2,rs62237617_DNMT3A_3")
    selected_contrasts=$(echo "Genotype_DNMT3A_vs_wt,Genotype_rs62237617_DNMT3A_vs_wt,Genotype_rs62237617_vs_wt")


    myjobid_bespoke_heatmap_function=$(sbatch --job-name=$name_bespoke_heatmap_function --output=$outfile_bespoke_heatmap_function --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=4096 --parsable --wrap="Rscript $Rscript_bespoke_heatmap_function --DE_results $DE_results --normalised_counts $normalised_counts --selected_annotations $selected_annotations --annotation_file $annotation_file --selected_annotations_other $selected_annotations_other --selected_clone_lines $selected_clone_lines --selected_contrasts $selected_contrasts --Diff_sel $Diff_array_sel --type $type --out $output_dir")
    myjobid_seff_bespoke_heatmap_function=$(sbatch --dependency=afterany:$myjobid_bespoke_heatmap_function --open-mode=append --output=$outfile_MSigDB_ORA --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_bespoke_heatmap_function >> $outfile_bespoke_heatmap_function")

    conda deactivate


done


