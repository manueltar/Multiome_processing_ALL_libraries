### Block 1

$ sbatch ~/Scripts/sbatch/5_Cell_ranger_arccount.sh /group/soranzo/manuel.tardaguila/2025_hESC_lymph_multiome/Multiome/ MCO_01330

### Block 2

$ bash ~/Scripts/Wraper_scripts/120_Seurat_first_v2.sh /group/soranzo/manuel.tardaguila/2025_hESC_lymph_multiome/Multiome/ processing_outputs

## have to be launch from /group/soranzo/manuel.tardaguila/2025_hESC_lymph_multiome/Multiome/processing_outputs/MCO_01331

$ sbatch ~/Scripts/sbatch/7_CellBender.sh /group/soranzo/manuel.tardaguila/2025_hESC_lymph_multiome/Multiome/ MCO_01330

$ bash ~/Scripts/Wraper_scripts/124_merge_atac_peaks.sh /group/soranzo/manuel.tardaguila/2025_hESC_lymph_multiome/Multiome/processing_outputs/

$ bash ~/Scripts/Wraper_scripts/122_snATAC_pipeline.sh /group/soranzo/manuel.tardaguila/2025_hESC_lymph_multiome/Multiome/ processing_outputs ########### There is an error in this code, it adds the sample names as the last lines of the merge file and it crashes later when using macs2 to do the peak calling

$ sbatch ~/Scripts/sbatch/8_merge_atac_peaks_v2.sh /group/soranzo/manuel.tardaguila/2025_hESC_lymph_multiome/Multiome/

$ sbatch ~/Scripts/sbatch/6_align_to_barcodes_v2.sh /group/soranzo/manuel.tardaguila/2025_hESC_lymph_multiome/Multiome/ /group/soranzo/manuel.tardaguila/Multiome/RITM0023280/special_reference_files/GFP_transgene_vCHEK2_and_DNMT3A.fa MCO_01330

$ bash ~/Scripts/Wraper_scripts/119_Filter_Larry_and_graphs_v2.sh /group/soranzo/manuel.tardaguila/2025_hESC_lymph_multiome/Multiome/deconvolute_LARRY/ count_and_filter /group/soranzo/manuel.tardaguila/2025_hESC_lymph_multiome/Multiome/deconvolute_LARRY/

### Block 3

$ bash ~/Scripts/Wraper_scripts/123_Amulet.sh /group/soranzo/manuel.tardaguila/2025_hESC_lymph_multiome/Multiome/ processing_outputs

### Block 4

$ bash ~/Scripts/Wraper_scripts/125_Seurat_second_pass_v2.sh


### Block 5

$ bash ~/Scripts/Wraper_scripts/126_merge_pre_merged_per_sample_v2.sh /group/soranzo/manuel.tardaguila/2025_hESC_lymph_multiome/Multiome/ processing_outputs
 /group/soranzo/manuel.tardaguila/2025_hESC_lymph_multiome/Multiome/processing_outputs/merged.atac_fragments.tsv.gz

### Block 6: Final QC

----> Jupyter notebook: Final_QC_in_the_merged_object.ipynb

### Recluster and export h5ad for rpca

$ bash ~/Scripts/Wraper_scripts/153_Recluster_and_export_h5ad.sh /group/soranzo/manuel.tardaguila/2025_hESC_lymph_multiome/Multiome/ processing_outputs

### Block 7 Cell typist prediction

----> Jupyter notebook: Cell_Typist_triple_prediction_cell_identity.ipynb
----> Jupyter notebook: mapping_cell_types.ipynb

### Block 8 Rpca #1

$ bash ~/Scripts/Wraper_scripts/163_rpca_integration.sh /group/soranzo/manuel.tardaguila/2025_hESC_lymph_multiome/Multiome/ processing_outputs

----> Jupyter notebook:RPCA_graphs_ALL_SAMPLES.ipynb
----> Jupyter notebook:RPCA_explore_ALL_SAMPLES.ipynb

### Block 9 notebook to assign barcodes

----> Jupyter notebook: notebook_to_assign_barcodes.ipynb

### Block 10 Rpca #2

$ bash /home/manuel.tardaguila/Scripts/Wraper_scripts/172_rpca_post_genotypiong.sh /group/soranzo/manuel.tardaguila/2025_hESC_lymph_multiome/Multiome/ processing_outputs

### Block 11 Final cell annotation

NOTEBOOK: /group/soranzo/manuel.tardaguila/2025_hESC_lymph_multiome/Multiome/processing_outputs/Post_G_final_cell_annotation.ipynb

### Block 12 Call ATAC peaks by

$ bash ~/Scripts/Wraper_scripts/138_MACS2_recall_peaks_by_cell_type_integrated_annotation.sh /group/soranzo/manuel.tardaguila/2025_hESC_lymph_multiome/Multiome/ processing_outputs

### Block 13 DE

$ bash ~/Scripts/Wraper_scripts/173_Multiome_DE_per_identity_both_Diffs.sh /group/soranzo/manuel.tardaguila/2025_hESC_lymph_multiome/Multiome/processing_outputs/Downstream_analysis/ DE_per_identity


<Bespoke Heatmaps>

$ bash ~/Scripts/Wraper_scripts/174_Multiome_bespoke_heatmaps_ALL_Diffs.sh /group/soranzo/manuel.tardaguila/2025_hESC_lymph_multiome/Multiome/processing_outputs/Downstream_analysis/ DE_per_identity /group/soranzo/manuel.tardaguila/2025_hESC_lymph_multiome/Multiome/processing_outputs/Downstream_analysis/DE_per_identity/genes_ORA_annotated_Diff_lymph.tsv

### Block 14 DA

$ bash ~/Scripts/Wraper_scripts/138_MACS2_recall_peaks_by_cell_type_integrated_annotation.sh /group/soranzo/manuel.tardaguila/2025_hESC_lymph_multiome/Multiome/ processing_outputs