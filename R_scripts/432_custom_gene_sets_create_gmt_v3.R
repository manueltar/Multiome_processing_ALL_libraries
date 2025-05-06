
suppressMessages(library("plyr", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("data.table", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("crayon", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("withr", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("ggplot2", lib.loc = "/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("farver", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("labeling", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("optparse", lib.loc = "/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("dplyr", lib.loc = "/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("withr", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("backports", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("broom", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("rstudioapi", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("cli", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("tzdb", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("svglite", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("ggeasy", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("sandwich", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("digest", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("tidyverse", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("RColorBrewer", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("svglite", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("cowplot", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("org.Hs.eg.db", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("ActivePathways", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("splitstackshape", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("dorothea", lib="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("decoupleR", lib="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("rWikiPathways", lib="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))

suppressMessages(library("splitstackshape", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))


opt = NULL

options(warn = 1)

create_gmt = function(option_list)
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
  
  cat("OUT_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  
  #### Table_of_gene_sets_Manual_curation ----
  
  Table_of_gene_sets<-read.table(opt$Table_of_gene_sets, sep="\t", header=T)
  

  cat("Table_of_gene_sets_0\n")
  cat(str(Table_of_gene_sets))
  cat("\n")
  cat(str(unique(Table_of_gene_sets$GeneSet)))
  cat("\n")
  
  Table_of_gene_sets_long<-unique(as.data.frame(cSplit(Table_of_gene_sets,sep = ',', direction = "long",
                                               splitCols = "Genes"),stringsAsFactors=F))
  
  cat("Table_of_gene_sets_long_0\n")
  cat(str(Table_of_gene_sets_long))
  cat("\n")
  cat(str(unique(Table_of_gene_sets_long$GeneSet)))
  cat("\n")
  
  Table_of_gene_sets_long$ENTREZID<-mapIds(org.Hs.eg.db, keys=Table_of_gene_sets_long$Genes, keytype="SYMBOL",column="ENTREZID")
  Table_of_gene_sets_long$ensembl_gene_id<-mapIds(org.Hs.eg.db, keys=Table_of_gene_sets_long$Genes, keytype="SYMBOL",column="ENSEMBL")
  
  cat("Table_of_gene_sets_long_1\n")
  cat(str(Table_of_gene_sets_long))
  cat("\n")
  
  Table_of_gene_sets_long_NO_NA<-Table_of_gene_sets_long[!is.na(Table_of_gene_sets_long$ENTREZID),]
  
  cat("Table_of_gene_sets_long_NO_NA_0\n")
  cat(str(Table_of_gene_sets_long_NO_NA))
  cat("\n")
  
  #### Prepare the file for gmt export ----
  
  for_gmt_Table_of_gene_sets_long_NO_NA<-unique(Table_of_gene_sets_long_NO_NA[,c(which(colnames(Table_of_gene_sets_long_NO_NA) == 'GeneSet'),which(colnames(Table_of_gene_sets_long_NO_NA) == 'ENTREZID'))])
  
  colnames(for_gmt_Table_of_gene_sets_long_NO_NA)[which(colnames(for_gmt_Table_of_gene_sets_long_NO_NA) == 'GeneSet')]<-'id'
  colnames(for_gmt_Table_of_gene_sets_long_NO_NA)[which(colnames(for_gmt_Table_of_gene_sets_long_NO_NA) == 'ENTREZID')]<-'gene'
  
  cat("for_gmt_Table_of_gene_sets_long_NO_NA_0\n")
  cat(str(for_gmt_Table_of_gene_sets_long_NO_NA))
  cat("\n")
  
  
  for_gmt_Table_of_gene_sets_long_NO_NA$name<-for_gmt_Table_of_gene_sets_long_NO_NA$id
  
  
  
  cat("for_gmt_Table_of_gene_sets_long_NO_NA_1\n")
  cat(str(for_gmt_Table_of_gene_sets_long_NO_NA))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(for_gmt_Table_of_gene_sets_long_NO_NA$name))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(for_gmt_Table_of_gene_sets_long_NO_NA$name)))))
  cat("\n")
  
  #### Reorder
  
  indx.reorder<-c(which(colnames(for_gmt_Table_of_gene_sets_long_NO_NA) == 'id'),which(colnames(for_gmt_Table_of_gene_sets_long_NO_NA) == 'name'),which(colnames(for_gmt_Table_of_gene_sets_long_NO_NA) == 'gene'))
  
  
  for_gmt_Table_of_gene_sets_long_NO_NA_reordered<-for_gmt_Table_of_gene_sets_long_NO_NA[,indx.reorder]
  
  cat("for_gmt_Table_of_gene_sets_long_NO_NA_reordered_0\n")
  cat(str(for_gmt_Table_of_gene_sets_long_NO_NA_reordered))
  cat("\n")
  
  #### miRNA_table_Manual_curation ----
  
  miRNA_table<-read.table(opt$miRNA_table, sep=",", header=T)
  
  
  cat("miRNA_table_0\n")
  cat(str(miRNA_table))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(miRNA_table$Evidence))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(miRNA_table$Evidence)))))
  cat("\n")
  
  miRNA_table$ENTREZID<-mapIds(org.Hs.eg.db, keys=miRNA_table$Target, keytype="SYMBOL",column="ENTREZID")
  miRNA_table$ensembl_gene_id<-mapIds(org.Hs.eg.db, keys=miRNA_table$Target, keytype="SYMBOL",column="ENSEMBL")
  
  cat("miRNA_table_1\n")
  cat(str(miRNA_table))
  cat("\n")
  
  miRNA_table_NO_NA<-miRNA_table[!is.na(miRNA_table$ENTREZID),]
  
  cat("miRNA_table_NO_NA_0\n")
  cat(str(miRNA_table_NO_NA))
  cat("\n")
  
  indx.experimental<-which(miRNA_table_NO_NA$Evidence%in%c('experimental (all)',
                                                           'experimental (all) + predicted (intersection) + predicted (union)',
                                                           'experimental (all) + predicted (union)'))
  
  
  
  
  
  miRNA_table_experimental<- data.frame(matrix(vector(), length(indx.experimental), 2,
                             dimnames=list(c(),
                                           c("GeneSet","ENTREZID"))),stringsAsFactors=F)
  
  
  
  miRNA_table_experimental$ENTREZID<-miRNA_table_NO_NA$ENTREZID[indx.experimental]
  miRNA_table_experimental$GeneSet<-rep("mir5739_experimental_targets",length(indx.experimental))                                    
  

  indx.intersection<-c(indx.experimental, which(miRNA_table_NO_NA$Evidence%in%c('predicted (intersection) + predicted (union)')))
  
  
  
  miRNA_table_intersection<- data.frame(matrix(vector(), length(indx.intersection), 2,
                                               dimnames=list(c(),
                                                             c("GeneSet","ENTREZID"))),stringsAsFactors=F)
  
  
  
  miRNA_table_intersection$ENTREZID<-miRNA_table_NO_NA$ENTREZID[indx.intersection]
  miRNA_table_intersection$GeneSet<-rep("mir5739_experimental_and_intersection_targets",length(indx.intersection))
  

  
  indx.union<-c(indx.experimental,indx.intersection, which(miRNA_table_NO_NA$Evidence%in%c('predicted (union)')))
  
  
  
  miRNA_table_union<- data.frame(matrix(vector(), length(indx.union), 2,
                                               dimnames=list(c(),
                                                             c("GeneSet","ENTREZID"))),stringsAsFactors=F)
  
  
  
  miRNA_table_union$ENTREZID<-miRNA_table_NO_NA$ENTREZID[indx.union]
  miRNA_table_union$GeneSet<-rep("mir5739_experimental_intersection_and_union_targets",length(indx.union))            
  
  
  #### Prepare the file for gmt export ----
  
  for_gmt_mir5739_targets<-rbind(miRNA_table_experimental,miRNA_table_intersection,miRNA_table_union)
  
  colnames(for_gmt_mir5739_targets)[which(colnames(for_gmt_mir5739_targets) == 'GeneSet')]<-'id'
  colnames(for_gmt_mir5739_targets)[which(colnames(for_gmt_mir5739_targets) == 'ENTREZID')]<-'gene'
  
  cat("for_gmt_mir5739_targets_0\n")
  cat(str(for_gmt_mir5739_targets))
  cat("\n")
  
  
  for_gmt_mir5739_targets$name<-for_gmt_mir5739_targets$id
  
  
  
  cat("for_gmt_mir5739_targets_1\n")
  cat(str(for_gmt_mir5739_targets))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(for_gmt_mir5739_targets$name))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(for_gmt_mir5739_targets$name)))))
  cat("\n")
  
  
  #### Reorder
  
  indx.reorder<-c(which(colnames(for_gmt_mir5739_targets) == 'id'),which(colnames(for_gmt_mir5739_targets) == 'name'),which(colnames(for_gmt_mir5739_targets) == 'gene'))
  
  
  for_gmt_mir5739_targets_reordered<-for_gmt_mir5739_targets[,indx.reorder]
  
  cat("for_gmt_mir5739_targets_reordered_0\n")
  cat(str(for_gmt_mir5739_targets_reordered))
  cat("\n")
  
  
 #### put together everything ----
  
  DEF<-rbind(for_gmt_Table_of_gene_sets_long_NO_NA_reordered,
             for_gmt_mir5739_targets_reordered)
  
  cat("DEF_0\n")
  cat(str(DEF))
  cat("\n")
 
 ### save as gmt ----
  
  setwd(out)
  
  writeGMT(DEF, paste("Custom_Soranzo",'_Hs.entrez.gmt',sep=''))
  
  
  
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
    make_option(c("--Table_of_gene_sets"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--miRNA_table"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--type"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--out"), type="character", default=NULL, 
                metavar="filename", 
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
  
  create_gmt(opt)

  
}


###########################################################################

system.time( main() )