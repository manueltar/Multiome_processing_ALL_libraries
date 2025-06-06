# This script scans specified regions for kmers or/and motifs using JASPAR2020 database.
# It outputs regions-by-kmers/motifs frequency matrix in .h5 format

# Author: Huidong Chen
# Contact information: hd7chen AT gmail DOT com

suppressMessages(library(optparse,quietly = TRUE))

main <- function(){
  option_list = list(
    make_option(c("-i", "--input"), type="character", default=NULL, 
                help="input region file name in .bed format", metavar="character"),
    make_option(c("-g", "--genome"), type="character", default=NULL, 
                help="Path to reference genome", metavar="character"),  
    make_option(c("--no_kmer"), action = "store_true",default=FALSE,
                help="disable scanning for kmers"),  
    make_option(c("--no_motif"), action = "store_true",default=FALSE,
                help="disable scanning for motifs"),          
    make_option(c("-k","--k_kmer"), type="integer", default=6, 
                help="k-mer length [default = %default].", metavar="integer"),       
    make_option(c("-s","--species"), type="character", default=NULL, 
                help="Species of motifs in the JASPAR database.
                Choose from 'Homo sapiens','Mus musculus'. Only valid when motif is used", 
                metavar="character"),
    make_option(c("-o", "--output"), type="character", default='output_kmers_motifs', 
                help="Output folder [default = %default]", metavar="character")    
  )
  
  opt_parser = OptionParser(option_list=option_list)
  opt = parse_args(opt_parser)
  
  if(is.null(opt$input)){
    print_help(opt_parser)
    stop("input region file must be specified", call.=FALSE)
  }
  if(!opt$no_motif){
    if(any(is.null(opt$genome),is.null(opt$species))){
      print_help(opt_parser)
      stop("reference genome and species must be both specified", call.=FALSE)
    }
  }
  
  file.input = opt$input
  genome = opt$genome
  no_kmer = opt$no_kmer
  no_motif = opt$no_motif
  k = opt$k_kmer
  species = gsub("_"," ",opt$species)
  dir.output = opt$output
  
  .libPaths()
  
  assign(".lib.loc", "/home/manuel.tardaguila/conda_envs/SIMBA_R/lib/R/library", envir = environment(.libPaths))
  
  .libPaths()
  
  suppressMessages(library(rhdf5))
  suppressMessages(library(HDF5Array))  # used for saving sparse matrix
  suppressMessages(library(Biostrings))
  suppressMessages(library(Matrix))
  suppressMessages(library(TFBSTools))
  suppressMessages(library(JASPAR2020))
  suppressMessages(library(motifmatchr))
  suppressMessages(library(SummarizedExperiment))
  suppressMessages(library(doParallel))
  suppressMessages(library(monaLisa))
  
  
  set.seed(2020)
  
  system(paste0('mkdir -p ',dir.output))

  print('Converting .bed to .fasta ...')
  ### convert peaks bed file to fasta file
  file.input.fa = paste0(basename(file.input),'.fa')
  system(paste("bedtools getfasta -fi",genome,
               "-bed",file.input,
               "-fo",file.path(dir.output,file.input.fa)))

  peaks_seq <- readDNAStringSet(file.path(dir.output,file.input.fa), "fasta")
  peaks_name = gsub(":|-",'_',names(peaks_seq))

  ### count kmers
  if(!no_kmer){
    print('Scanning for kmers ...')
    freq_k = oligonucleotideFrequency(peaks_seq, k)
    rownames(freq_k) = peaks_name
    freq_k = as(freq_k, "sparseMatrix")
  }

  ### scan for TF motifs
  if(!no_motif){
    print('Scanning for TF motifs JASPAR2020')
    opts <- list()
    opts["species"] <- species
    opts["collection"] <- "CORE"
    PFMatrixList = TFBSTools::getMatrixSet(JASPAR2020::JASPAR2020,opts = opts)
    motif_ix_scores <- motifmatchr::matchMotifs(PFMatrixList,peaks_seq, out = "scores")
    freq_motif = motifCounts(motif_ix_scores)
    motif_names = c()
    for (x in names(PFMatrixList)){
        motif_names = c(motif_names,PFMatrixList[[x]]@name)
    }
    colnames(freq_motif) = gsub("::",'_',motif_names)
    rownames(freq_motif) = peaks_name
    
    
    print('Scanning for TF motifs HOMER')
    
    homer_motif_file <- file.path("/home/manuel.tardaguila/HOMER_db/HOMER_motifs.txt")
    
    ## from monaLisa package
    
    PFMatrixList_homer <- homerToPFMatrixList(homer_motif_file)
    
    ## motif match against the peaks
    
    motif_ix_scores_homer <- motifmatchr::matchMotifs(PFMatrixList_homer, peaks_seq, out = "scores")
    
    ## access counts from Homer
    
    freq_motif_homer <- assays(motif_ix_scores_homer)[["motifCounts"]]
    rownames(freq_motif_homer) <- peaks_name
    colnames(freq_motif_homer) <- motif_ix_scores_homer$name
    
    ### Merge JASPAR and HOMER
    
    print('Merging TF motifs from JASPAR and HOMER')
    
    
    if(length(freq_motif_homer) > 0){
      
      # Check if the row names are identical and in the same order
      identical_rownames <- identical(rownames(freq_motif), rownames(freq_motif_homer))
      
      if (identical_rownames) {
        # Merge the two sparse matrices by columns
        merged_motif <- cbind(freq_motif, freq_motif_homer)
        
        # The column names will be a concatenation of the column names
        # from both original matrices. You might want to inspect or modify them.
        print("Dimensions of the merged motif matrix:")
        print(dim(merged_motif))
        print("First few column names of the merged matrix:")
        print(head(colnames(merged_motif)))
        
      } else {
        print("Warning: Row names of freq_motif and freq_motif_homer are not identical. Cannot merge directly by columns.")
        # If row names are not identical, you would need to align them first
        # before merging, which could involve subsetting or reordering rows.
      }
      
    }else{
      
      merged_motif<-freq_motif
      
    }#length(freq_motif_homer) > 0
    
    
    
    
  }

  ### save results
  ### save kmers
  if(!no_kmer){
    print('Saving kmer matrix ...')
    
   
    filename = 'freq_kmer.h5'
    filepath = file.path(dir.output, filename)
    
    if (file.exists(filepath)) {
      print(paste("Warning:", filepath, "already exists. Deleting it to handle existing files."))
      # You might want to add code here to delete the file:
      file.remove(filepath)
    }
    
    # writeHDF5Array internally transposes the matrix so `t()` is used to counteract this operation
    writeHDF5Array(t(freq_k), file.path(dir.output,filename), name="mat", with.dimnames=FALSE, verbose=FALSE)
    # using this structure in order for  anndata 'read_hdf' to recognize row names and column names
    h5write(rownames(freq_k), file.path(dir.output,filename), "row_names")
    h5write(colnames(freq_k), file.path(dir.output,filename), "col_names")
  }

  ### save motifs
  if(!no_motif){
    print('Saving motif matrix ...')

    filename = 'freq_motif.h5'
    filepath = file.path(dir.output, filename)
    
    if (file.exists(filepath)) {
      print(paste("Warning:", filepath, "already exists. Deleting it to handle existing files."))
      # You might want to add code here to delete the file:
      file.remove(filepath)
    }
    

    # writeHDF5Array internally transposes the matrix so `t()` is used to counteract this operation
    writeHDF5Array(t(merged_motif), file.path(dir.output,filename), name="mat", with.dimnames=FALSE, verbose=FALSE)
    # using this structure in order for  anndata 'read_hdf' to recognize row names and column names
    h5write(rownames(merged_motif), file.path(dir.output,filename), "row_names")
    h5write(colnames(merged_motif), file.path(dir.output,filename), "col_names")
  }
  
  print('Finished.')
}

main()