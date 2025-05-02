### Block 1

$ sbatch ~/Scripts/sbatch/5_Cell_ranger_arccount.sh /group/soranzo/manuel.tardaguila/2025_hESC_lymph_multiome/Multiome/ MCO_01330

### Block 2

$ bash ~/Scripts/Wraper_scripts/120_Seurat_first_v2.sh /group/soranzo/manuel.tardaguila/2025_hESC_lymph_multiome/Multiome/ processing_outputs

## have to be launch from /group/soranzo/manuel.tardaguila/2025_hESC_lymph_multiome/Multiome/processing_outputs/MCO_01331

$ sbatch ~/Scripts/sbatch/7_CellBender.sh /group/soranzo/manuel.tardaguila/2025_hESC_lymph_multiome/Multiome/ MCO_01330

$ sbatch ~/Scripts/sbatch/8_merge_atac_peaks_v2.sh /group/soranzo/manuel.tardaguila/2025_hESC_lymph_multiome/Multiome/

$ bash ~/Scripts/Wraper_scripts/122_snATAC_pipeline.sh /group/soranzo/manuel.tardaguila/2025_hESC_lymph_multiome/Multiome/ processing_outputs 

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

#### Index the barcode genome with cellranger

$ cellranger mkref --fasta /group/soranzo/manuel.tardaguila/Multiome/RITM0023280/special_reference_files/GFP_transgene_vCHEK2_and_DNMT3A.fa --genes /group/soranzo/manuel.tardaguila/Multiome/RITM0023280/special_reference_files/STAR.gtf --genome GFP_transgene_vCHEK2_and_DNMT3A_cellranger


#### Modify targeted amplification of GEX to make it pass as CellRanger input

$ bash ~/Scripts/Wraper_scripts/127_cellranger_alignment_of_targeted_amp_GEX.sh /group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/GEX_reseq/ alignment /group/soranzo/manuel.tardaguila/Multiome/MCO_20250123/250124_A02059_0109_AHWTHYDSXC/adapter_trimmed_fastq/

$ bash ~/Scripts/Wraper_scripts/146_cellranger_alignment_of_targeted_amp_GEX_lymph.sh /group/soranzo/manuel.tardaguila/2025_hESC_lymph_multiome/ targeted_amplicon_GEX /group/soranzo/manuel.tardaguila/Multiome/QC_20250310/250310_A02059_0119_AH2N2TDMX2/adapter_trimmed_fastq/

$ bash ~/Scripts/Wraper_scripts/171_cellranger_alignment_of_targeted_amp_GEX_lymph.sh /group/soranzo/manuel.tardaguila/2025_hESC_lymph_multiome/ targeted_amplicon_GEX /group/soranzo/manuel.tardaguila/Multiome/RITM0029357/250414_A01481_0283_BHMKVYDRX5/fastq_raw/

#### Run cellranger on the adapted reads

$ sbatch ~/Scripts/Wraper_scripts/128_cellranger_MCO_1326.sh /group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/GEX_reseq/alignment/cellranger/
$ sbatch ~/Scripts/Wraper_scripts/128_cellranger_MCO_1327.sh /group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/GEX_reseq/alignment/cellranger/
$ sbatch ~/Scripts/Wraper_scripts/128_cellranger_MCO_1328.sh /group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/GEX_reseq/alignment/cellranger/
$ sbatch ~/Scripts/Wraper_scripts/128_cellranger_MCO_1329.sh /group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/GEX_reseq/alignment/cellranger/

#### Run CellBender to correct background of empty beads

$ sbatch ~/Scripts/Wraper_scripts/131_Cell_Bender_for_targeted_amplified_libraries.sh /group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/GEX_reseq/alignment/cellranger/ MCO_01326
$ sbatch ~/Scripts/Wraper_scripts/131_Cell_Bender_for_targeted_amplified_libraries.sh /group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/GEX_reseq/alignment/cellranger/ MCO_01327
$ sbatch ~/Scripts/Wraper_scripts/131_Cell_Bender_for_targeted_amplified_libraries.sh /group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/GEX_reseq/alignment/cellranger/ MCO_01328
$ sbatch ~/Scripts/Wraper_scripts/131_Cell_Bender_for_targeted_amplified_libraries.sh /group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/GEX_reseq/alignment/cellranger/ MCO_01329

#### ----> Jupyter notebook: notebook_to_assign_barcodes.ipynb

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