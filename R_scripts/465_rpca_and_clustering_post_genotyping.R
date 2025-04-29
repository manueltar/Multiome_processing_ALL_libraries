
.libPaths()
.libPaths(new = c("/home/manuel.tardaguila/conda_envs/multiome_QC_DEF/lib/R/library"))
.libPaths()
# sessionInfo()

Sys.setenv(RETICULATE_PYTHON="/home/manuel.tardaguila/conda_envs/multiome_QC_DEF/bin/python")
library(reticulate)
reticulate::use_python("/home/manuel.tardaguila/conda_envs/multiome_QC_DEF/bin/python")
reticulate::use_condaenv("/home/manuel.tardaguila/conda_envs/multiome_QC_DEF")
reticulate::py_module_available(module='leidenalg')
reticulate::import('leidenalg')
suppressMessages(library("optparse"))
suppressMessages(library(hdf5r))
suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(EnsDb.Hsapiens.v86))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(Matrix))
suppressMessages(library(data.table))
suppressMessages(library(ggpubr))
suppressMessages(library(ggplot2))
suppressMessages(library(scDblFinder))
suppressMessages(library("tidyr"))
suppressMessages(library("tibble"))
suppressMessages(library("biovizBase"))
suppressMessages(library("patchwork"))
suppressMessages(library(glmGamPoi))
suppressMessages(library(SeuratData))
suppressMessages(library(SeuratDisk))



opt = NULL

options(warn = -1)

rpca_integration = function(option_list)
{
  
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  
 
  #### READ and transform type ----
  
  type = opt$type
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type)))
  cat("\n")
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("out_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  #### Read filtered object by doublets -----
  
  
  adata<-readRDS(file=opt$object_post_genotyping)
  
  cat("adata_0\n")
  # cat(str(adata))
  cat("\n")
  
 
  
  
  
  #### Perform SCT normalization ----------------------
  
  
  DefaultAssay(adata)<-'RNA'
  
  adata.list <- SplitObject(adata, split.by = "orig.ident")
  
  rm(adata)
  
  adata.list <- lapply(X = adata.list, FUN = SCTransform, method = "glmGamPoi")
  
  
  ## select features that are repeatedly variable across datasets for integration run PCA on each dataset using these features --------------------
  
  
  features <- SelectIntegrationFeatures(object.list = adata.list, nfeatures = 3000)
  
  cat("features_0\n")
  # cat(str(features))
  cat("\n")
  
  adata.list <- PrepSCTIntegration(object.list = adata.list, anchor.features = features)
  
  ## Run PCA on features ----------------------------------------------------------
  
  adata.list <- lapply(X = adata.list, FUN = RunPCA, features = features)
  
  
  ## Find integrator anchors RPCA ----------------------------------------------------------
  
  
  anchors <- FindIntegrationAnchors(object.list = adata.list, normalization.method = "SCT",
                                    anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 5)
  
  
  cat("anchors_0\n")
  # cat(str(anchors))
  cat("\n")
  
  ## IntegrateData ---------------------------------------------------------------------------
  
  cat("Run IntegrateData\n")
  
  
  combined.sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = 1:30)
  
  
  ## Run PCA on integrated -----------------------------------------------------------------------
  
  cat("Run PCA\n")
  
  combined.sct <- RunPCA(combined.sct, verbose = FALSE, reduction.name= "PCA_POST_G")
   
  
  ## Run UMAP on integrated -----------------------------------------


  cat("Run UMAP\n")

  combined.sct <- RunUMAP(combined.sct, reduction = "PCA_POST_G", dims = 1:30, reduction.name="UMAP_POST_G")
  
  cat("Saving intermediate\n")
  
  setwd(out)
  
  saveRDS(combined.sct, file="merged_post_rpca_intermediate.rds")
  
  
  #### Name the new graph -----
  
  cat("FindNeighbors\n")
  
  combined.sct <- FindNeighbors(combined.sct, reduction = "PCA_POST_G", graph.name = "integrated_nn_POST_G")

  #### Find new clusters at resolution 0.5 ----
  
  cat("FindClusters\n")
  
  
  combined.sct <- FindClusters(combined.sct, graph.name='integrated_nn_POST_G', algorithm=4, resolution = 0.5, verbose=FALSE, method = "igraph")
  
  #### Prepare SCT for Find Markers see https://satijalab.org/seurat/reference/prepsctfindmarkers ----
  
  cat("PrepSCTFindMarkers\n")
  
  
  combined.sct<-PrepSCTFindMarkers(combined.sct, assay = "SCT", verbose = TRUE)
  
  #### Find marker genes ----
  
  cat("FindAllMarkers\n")
  
  
  combined.sct.markers <- FindAllMarkers(combined.sct, logfc.threshold = 0.25, test.use = "roc")
  
  cat("combined.sct.markers_0\n")
  cat(str(combined.sct.markers))
  cat("\n")
  
  
  
  ##### Saving ---------------
  
  cat("Saving\n")
  
  setwd(out)
  
  saveRDS(combined.sct, file="merged_post_rpca_reclustered.rds")
  
  saveRDS(combined.sct.markers, file="merged_post_rpca_reclustered_MARKER_GENES.rds")
  
  
  
  
 
}


printList = function(l, prefix = "    ") {
  list.df = data.frame(val_name = names(l), value = as.character(l))
  list_strs = apply(list.df, MARGIN = 1, FUN = function(x) { paste(x, collapse = " = ")})
  cat(paste(paste(paste0(prefix, list_strs), collapse = "\n"), "\n"))
}


#### main script ----

main = function() {
  cmd_line = commandArgs()
  cat("Command line:\n")
  cat(paste(gsub("--file=", "", cmd_line[4], fixed=T),
            paste(cmd_line[6:length(cmd_line)], collapse = " "),
            "\n\n"))
  option_list <- list(
    make_option(c("--object_post_genotyping"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--type"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--out"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required.")
  )
  parser = OptionParser(usage = "140__Rscript_v106.R
                        --subset type
                        --TranscriptEXP FILE.txt
                        --cadd FILE.txt
                        --ncboost FILE.txt
                        --type type
                        --out filename",
                        option_list = option_list)
  opt <<- parse_args(parser)
  
  rpca_integration(opt)
 

}


###########################################################################

system.time( main() )
  