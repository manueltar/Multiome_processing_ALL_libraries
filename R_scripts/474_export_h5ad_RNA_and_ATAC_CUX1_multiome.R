


.libPaths()

assign(".lib.loc", "/home/manuel.tardaguila/conda_envs/multiome_QC_DEF/lib/R/library", envir = environment(.libPaths))

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
suppressMessages(library(plyr))
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
suppressMessages(library(SeuratData))
suppressMessages(library(SeuratDisk))



opt = NULL

options(warn = -1)

cluster_at_low_res_and_exporth5ad = function(option_list)
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
  
 
  #### Read seurat object -----
  
  
  adata<-readRDS(file=opt$Seurat_object)
  
  cat("adata_0\n")
  # cat(str(adata))
  cat("\n")
  
  adata@meta.data$clone_line<-revalue(adata@meta.data$Assigned_GFPbc,
                                      c('chrGFP_WTA' = 'wt_1',
                                        'chrGFP_WTB' = 'wt_2',
                                        'chrGFP_WTC' = 'wt_3',
                                        'chrGFP_KI_13' = 'rs139141690_1',
                                        'chrGFP_KI_27' = 'rs139141690_2',
                                        'chrGFP_KI_29' = 'rs139141690_3',
                                        'chrGFP_HET' = 'rs139141690_HET_1',
                                        'chrGFP_Del_16bp' = 'Del_16bp_1',
                                        'chrGFP_Del_233' = 'Del_80bp_1',
                                        'chrGFP_Del_235' = 'Del_80bp_2',
                                        'chrGFP_Del_287' = 'Del_80bp_3'))
  
  
  adata@meta.data$Genotype<-NA
  
  
  adata@meta.data$Genotype[which(adata@meta.data$clone_line%in%c('wt_1','wt_2','wt_3'))]<-"wt"
  adata@meta.data$Genotype[which(adata@meta.data$clone_line%in%c('rs139141690_1','rs139141690_2','rs139141690_3'))]<-"rs139141690"
  adata@meta.data$Genotype[which(adata@meta.data$clone_line%in%c('rs139141690_HET_1'))]<-"rs139141690_HET"
  adata@meta.data$Genotype[which(adata@meta.data$clone_line%in%c('Del_16bp_1'))]<-"Del_16bp"
  adata@meta.data$Genotype[which(adata@meta.data$clone_line%in%c('Del_80bp_1','Del_80bp_2','Del_80bp_3'))]<-"Del_80bp"
  
  
  
  
   
   
   
   adata@meta.data$seurat_clusters<-as.character(adata@meta.data$seurat_clusters)
   adata@meta.data$clone_line<-as.character(adata@meta.data$clone_line)
   adata@meta.data$Genotype<-as.character(adata@meta.data$Genotype)
   adata@meta.data$time_point<-as.character(adata@meta.data$time_point)
   adata@meta.data$diff_groups<-as.character(adata@meta.data$diff_groups)
   

   met<-adata[[]]
   
   cat("met_0\n")
   cat(str(met))
   cat("\n")

  ##### Reduce the Seurat object to h5ad with RNA counts not corrected by CellBender It doesn't work with CellBender corrected counts------

  DefaultAssay(adata)<-'RNA_raw'
  RNA_only<-DietSeurat(adata, assays = "RNA_raw")


  setwd(out)
  
  unlink(c("RNA.h5Seurat","RNA.h5ad"))

  SaveH5Seurat(RNA_only, filename = "RNA.h5Seurat")
  Convert("RNA.h5Seurat", dest = "h5ad")


  ##### Reduce the Seurat object to h5ad with ATAC counts not corrected by CellBender It doesn't work with CellBender corrected counts------

  DefaultAssay(adata)<-'ATAC'
  
  
  
  
  Peaks<-Features(adata)
  
  cat("Peaks_0\n")
  cat(str(Peaks))
  cat("\n")
  
  Peaks<-unique(Peaks)
  
  cat("Peaks_1\n")
  cat(str(Peaks))
  cat("\n")
  
  
  tmp.gather<- data.frame(matrix(vector(), length(Peaks), 4,
                                 dimnames=list(c(),
                                               c("chr","start","end","name"))),
                          stringsAsFactors=F)
  

  tmp.gather$name<-Peaks
  tmp.gather$chr<-gsub("-.+$","",Peaks)
  
  tmp.gather$start<-gsub("^[^-]+-","",Peaks)
  tmp.gather$start<-as.integer(gsub("-.+$","",tmp.gather$start))
  tmp.gather$end<-as.integer(gsub("^[^-]+-[^-]+-","",Peaks))
  
  
  cat("tmp.gather_0\n")
  cat(str(tmp.gather))
  cat("\n")
  
  tmp.gather$chr<-factor(tmp.gather$chr,
                           levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11",
                                    "chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21",
                                    "chr22","chr23","chrX","chrY"), ordered=T)
  
  tmp.gather<-tmp.gather[order(tmp.gather$chr, tmp.gather$start, decreasing = F),]
  
  cat("tmp.gather_1\n")
  cat(str(tmp.gather))
  cat("\n")
  

  # Peaks$chr<-gsub("-.+$","",)


  
  #### remove motifs before writing h5ad ----
  
  if (!is.null(adata[["ATAC"]]@motifs)) {
    cat("Stripping motif information from ATAC assay to prevent HDF5 write errors\n")
    adata[["ATAC"]]@motifs <- NULL
  }
  
  ATAC_only <- DietSeurat(adata, assays = "ATAC", counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, dimreducs = NULL, graphs = NULL, misc = FALSE)
  
  setwd(out)
  
  unlink(c("ATAC.h5Seurat","ATAC.h5ad"))

  SaveH5Seurat(ATAC_only, filename = "ATAC.h5Seurat")
  Convert("ATAC.h5Seurat", dest = "h5ad")

  
  write.table(tmp.gather, file="Peaks.bed", sep="\t", row.names = F,quote = F,col.names = F)
  
  
  
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
    make_option(c("--Seurat_object"), type="character", default=NULL, 
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
  
  cluster_at_low_res_and_exporth5ad(opt)
 

}


###########################################################################

system.time( main() )
  