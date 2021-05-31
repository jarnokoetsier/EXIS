

#############   ####        ####   ####      #########
#############    ####      ####    ####    ##########
#####             ###########      ####   #####
#####               #######        ####   #####
#############     ###########      ####    ######## 
#############    ####      ####    ####      ########
#####           ####        ####   ####         ######
#####                                            #####
#####################################################
#       Isoform-specific expression analysis        #
#####################################################
#CREATED BY JARNO KOETSIER


###############################################################################################################################

#INSTRUCTIONS GLOBAL SCRIPT

###############################################################################################################################



#1. Source run this script before using the EXIS app.

#   IMPORTANT: The current analysis depends on the modified affy package from Brainarray.
#   If you have installed the original affy package previously, 
#   please remove this package from your package library before source running this script.


#2. Make sure that the exon/transcript-specific CDF file is in your working directory 
#   OR has been installed as a CDF package previously.


#3. When you would like to use the exon-specific probe set definition, 
#   also make sure to place exon annotation file in your working directory





#################################################################################################################################

#Install and load packages

#################################################################################################################################

pkg <- installed.packages()[, "Package"]


#Modified affy package from Brainarray

if(!('affy' %in% pkg)) {
  install.packages("http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/25.0.0/affy_1.68.0.tar.gz", type = "source", repos = NULL)
}


#Bioconductor packages

s_requiredpackages =
  c(
    "biomaRt",
    "GEOquery",
    "limma"
  )

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", ask = F)

for (i in s_requiredpackages) {
  if (!requireNamespace(i, quietly = TRUE))
    BiocManager::install(i, ask = F)
  require(as.character(i), character.only = TRUE)
}  


#Other packages

if(!('tidyverse' %in% pkg)){
  install.packages("tidyverse")
}

if(!('fuzzyjoin' %in% pkg)){
  install.packages("fuzzyjoin")
}

if(!("plotly" %in% pkg)){
  install.packages("plotly")
}

if(!("makecdfenv" %in% pkg)){
  install.packages("makecdfenv")
}


#Load packages
library(affy)
library(GEOquery)
library(biomaRt)
library(limma)
library(fuzzyjoin)
library(tidyverse)
library(plotly)
library(makecdfenv)


###################################################################################################################################

#get_grouping

###################################################################################################################################

get_grouping <- function(gset) {
  
  GEO_accession <- str_remove(names(gset), "_series_matrix.txt.gz")
  
  grouping <- as.data.frame(gset[[paste0(GEO_accession,"_series_matrix.txt.gz")]]@phenoData@data)[1]
  
  for (l in 1:length(gset[[paste0(GEO_accession,"_series_matrix.txt.gz")]]@phenoData@data)){
    if (((length(unique(gset[[paste0(GEO_accession,"_series_matrix.txt.gz")]]@phenoData@data[[l]]))) > 1) & 
        ((length(unique(gset[[paste0(GEO_accession,"_series_matrix.txt.gz")]]@phenoData@data[[l]]))) < length(gset[[paste0(GEO_accession,"_series_matrix.txt.gz")]]@phenoData@data[[l]]))){
      
      grouping <- cbind(grouping, as.data.frame(gset[[paste0(GEO_accession,"_series_matrix.txt.gz")]]@phenoData@data)[l])
    }
    
  }
  
  grouping <- as.data.frame(grouping[,-1])
  
  #remove columns with the same information
  uni = c(1)
  
  suppressWarnings(
    
    for (i in 1:(ncol(grouping)-1)) {
      for (j in (i+1):ncol(grouping)) {
        if (all(str_detect(grouping[,i], grouping[,j]))) {
          if (length(unique(str_remove(grouping[,i], grouping[,j]))) == 1){
            uni <- c(uni, i)
          }
        }
        if (all(str_detect(grouping[,j], grouping[,i]))) {
          if (length(unique(str_remove(grouping[,j], grouping[,i]))) == 1) {
            uni <- c(uni, j)
          }
        }
      }
    }
  )
  
  uni <- unique(uni[-1])
  grouping <- grouping[, -uni]
  rownames(grouping) <- NULL
  
  return(grouping)
  
}


###################################################################################################################################

#auto_group

###################################################################################################################################

auto_group <- function(gset, attempt = 1){
  
  GEO_accession <- str_remove(names(gset), "_series_matrix.txt.gz")
  
  participants <- gset[[paste0(GEO_accession,"_series_matrix.txt.gz")]]@phenoData@data[["geo_accession"]]
  
  participant_group <- NULL
  
  n_ch <- c(100,100)
  for (l in 1:length(gset[[paste0(GEO_accession,"_series_matrix.txt.gz")]]@phenoData@data)){
    if ((length(grep("control|non|healthy|treat", gset[[paste0(GEO_accession,"_series_matrix.txt.gz")]]@phenoData@data[[l]], ignore.case = TRUE)) > 0) &
        (length(unique(gset[[paste0(GEO_accession,"_series_matrix.txt.gz")]]@phenoData@data[[l]]))) < length(gset[[paste0(GEO_accession,"_series_matrix.txt.gz")]]@phenoData@data[[l]])){
      
      n_ch <- rbind(n_ch, c(l,nchar(gset[[paste0(GEO_accession,"_series_matrix.txt.gz")]]@phenoData@data[[l]][1])))
    }
  }
  
  if (attempt < nrow(n_ch)) {
    y = 1
    repeat {
      
      if (y == attempt) {
        break
      }
      y = y + 1
      n_ch = n_ch[-which.min(n_ch[,2]),]
    }
    
    int1 <- n_ch[which.min(n_ch[,2]),1]
    participant_group <- as.data.frame(gset[[paste0(GEO_accession,"_series_matrix.txt.gz")]]@phenoData@data)[int1]
    print(noquote("Grouping has been done automatically by default. Please check whether grouping has occured correctly."))
    
  }
  
  if (attempt >= nrow(n_ch)) {
    print(noquote("grouping information not found. Upload grouping data manually at grouping_column"))
  }
  
  return(colnames(participant_group))
}


###################################################################################################################################

#get_meta

###################################################################################################################################


get_meta <- function(gset, attempt = 1, grouping_column = NULL, pairing_column = NULL, file_name = NULL) {
  
  GEO_accession <- str_remove(names(gset), "_series_matrix.txt.gz")
  
  participants <- gset[[paste0(GEO_accession,"_series_matrix.txt.gz")]]@phenoData@data[["geo_accession"]]
  
  
  #automatic selection
  if (is.null(grouping_column)){
    
    n_ch <- c(100,100)
    for (l in 1:length(gset[[paste0(GEO_accession,"_series_matrix.txt.gz")]]@phenoData@data)){
      if ((length(grep("control|non|healthy|treat", gset[[paste0(GEO_accession,"_series_matrix.txt.gz")]]@phenoData@data[[l]], ignore.case = TRUE)) > 0) &
          (length(unique(gset[[paste0(GEO_accession,"_series_matrix.txt.gz")]]@phenoData@data[[l]]))) < length(gset[[paste0(GEO_accession,"_series_matrix.txt.gz")]]@phenoData@data[[l]])){
        
        n_ch <- rbind(n_ch, c(l,nchar(gset[[paste0(GEO_accession,"_series_matrix.txt.gz")]]@phenoData@data[[l]][1])))
      }
    }
    
    if (attempt < nrow(n_ch)) {
      y = 1
      repeat {
        
        if (y == attempt) {
          break
        }
        y = y + 1
        n_ch = n_ch[-which.min(n_ch[,2]),]
      }
      
      int1 <- n_ch[which.min(n_ch[,2]),1]
      participant_group <- gset[[paste0(GEO_accession,"_series_matrix.txt.gz")]]@phenoData@data[[int1]]
      print(noquote("Grouping has been done automatically by default. Please check whether grouping has occured correctly."))
      
    }
    
    if (attempt >= nrow(n_ch)) {
      print(noquote("grouping information not found. Upload grouping data manually at grouping_column"))
    }
  }
  
  
  #use grouping column
  if (!is.null(grouping_column)){
    if (!is.vector(grouping_column)){
      if (ncol(grouping_column) > 1) {
        participant_group <- unite(grouping_column, col = grouping_column, sep = ".")
      }
      if (ncol(grouping_column) == 1){
        participant_group <- grouping_column
      } 
    }
    if (is.vector(grouping_column)){
      participant_group <- grouping_column
    }
    
  }
  
  
  
  #independ samples
  if (is.null(pairing_column) | length(pairing_column) < 1) {
    if (exists("participant_group")){
      meta <- as.data.frame(cbind(participants, participant_group))
      rownames(meta) <- NULL
      colnames(meta) <- c("GEO ID", "Grouping")
      return(meta)
    }
  }
  
  
  
  #dependent samples
  if (!is.null(pairing_column) & length(pairing_column) >= 1){
    if (!is.vector(pairing_column)){
      if (ncol(pairing_column) > 1) {
        pairs <- unite(pairing_column, col = pairing_column, sep = ".")
      }
      if (ncol(pairing_column) == 1){
        pairs <- pairing_column
      } 
    }
    if (is.vector(pairing_column)){
      pairs <- pairing_column
    }
    
    meta <- as.data.frame(cbind(participants, participant_group, pairs))
    rownames(meta) <- NULL
    colnames(meta) <- c("GEO ID", "Grouping", "Pairing")
    return(meta)
  }
  
  
  
  #file upload
  if (!is.null(file_name)){
    meta <- read_excel(file_name)
    return(meta)
  }
}



###################################################################################################################################

#get_chiptype

###################################################################################################################################

get_chiptype <- function(gset){
  GEO_accession <- str_remove(names(gset), "_series_matrix.txt.gz")
  
  chiptype = "unknown"
  organism = "unknown"
  
  for (j in 1:length(gset[[paste0(GEO_accession,"_series_matrix.txt.gz")]]@phenoData@data)) {
    if (length(grep("huex", gset[[paste0(GEO_accession,"_series_matrix.txt.gz")]]@phenoData@data[[j]], ignore.case = TRUE)) > 0) {
      chiptype = "huex10st"
      organism = "hs"
    }
  }
  
  
  for (j in 1:length(gset[[paste0(GEO_accession,"_series_matrix.txt.gz")]]@phenoData@data)) {
    if (length(grep("hugene-1_0 | human gene 1.0", gset[[paste0(GEO_accession,"_series_matrix.txt.gz")]]@phenoData@data[[j]], ignore.case = TRUE)) > 0) {
      chiptype = "hugene10st"
      organism = "hs"
    }
  }
  
  
  for (j in 1:length(gset[[paste0(GEO_accession,"_series_matrix.txt.gz")]]@phenoData@data)) {
    if (length(grep("hugene-2_0 | human gene 2.0", gset[[paste0(GEO_accession,"_series_matrix.txt.gz")]]@phenoData@data[[j]], ignore.case = TRUE)) > 0) {
      chiptype = "hugene20st"
      organism = "hs"
    }
  }
  
  for (j in 1:length(gset[[paste0(GEO_accession,"_series_matrix.txt.gz")]]@phenoData@data)) {
    if (length(grep("HG-U133_Plus_2 | Human Genome U133 Plus 2.0", gset[[paste0(GEO_accession,"_series_matrix.txt.gz")]]@phenoData@data[[j]], ignore.case = TRUE)) > 0) {
      chiptype = "hgu133plus2"
      organism = "hs"
    }
  }
  
  
  if (chiptype == "unknown") {
    print(noquote("Chip type could not be found. Please select chip type and organism manually"))
  }
  
  
  if (chiptype != "unknown") {
    print(noquote(paste(chiptype, "chip type was selected.", 
                        "If this is the incorrect chip type, please select the chip type manually", sep = " ")))
  }
  
  return(chiptype)
}




###################################################################################################################################

#get_organism

###################################################################################################################################

get_organism <- function(gset){
  GEO_accession <- str_remove(names(gset), "_series_matrix.txt.gz")
  
  chiptype = "unknown"
  organism = "unknown"
  
  for (j in 1:length(gset[[paste0(GEO_accession,"_series_matrix.txt.gz")]]@phenoData@data)) {
    if (length(grep("huex", gset[[paste0(GEO_accession,"_series_matrix.txt.gz")]]@phenoData@data[[j]], ignore.case = TRUE)) > 0) {
      chiptype = "huex10st"
      organism = "hs"
    }
  }
  
  
  for (j in 1:length(gset[[paste0(GEO_accession,"_series_matrix.txt.gz")]]@phenoData@data)) {
    if (length(grep("hugene-1_0 | human gene 1.0", gset[[paste0(GEO_accession,"_series_matrix.txt.gz")]]@phenoData@data[[j]], ignore.case = TRUE)) > 0) {
      chiptype = "hugene10st"
      organism = "hs"
    }
  }
  
  
  for (j in 1:length(gset[[paste0(GEO_accession,"_series_matrix.txt.gz")]]@phenoData@data)) {
    if (length(grep("hugene-2_0 | human gene 2.0", gset[[paste0(GEO_accession,"_series_matrix.txt.gz")]]@phenoData@data[[j]], ignore.case = TRUE)) > 0) {
      chiptype = "hugene20st"
      organism = "hs"
    }
  }
  
  for (j in 1:length(gset[[paste0(GEO_accession,"_series_matrix.txt.gz")]]@phenoData@data)) {
    if (length(grep("HG-U133_Plus_2 | Human Genome U133 Plus 2.0", gset[[paste0(GEO_accession,"_series_matrix.txt.gz")]]@phenoData@data[[j]], ignore.case = TRUE)) > 0) {
      chiptype = "hgu133plus2"
      organism = "hs"
    }
  }
  
  
  if (chiptype == "unknown") {
    print(noquote("Chip type could not be found. Please select chip type and organism manually"))
  }
  
  
  if (chiptype != "unknown") {
    print(noquote(paste(organism, " was selected.", 
                        "If this is the incorrect organism, please select the organism manually", sep = " ")))
  }
  
  return(organism)
}



###################################################################################################################################

#readcels

###################################################################################################################################

readcels <- function(gset, chiptype, organism , version = 25, annotation = "enst", robust = TRUE, outliers = NULL) {
  
  pkg <- installed.packages()[, "Package"]
  GEO_accession <- str_remove(names(gset), "_series_matrix.txt.gz")
  exdir = paste0("data_", GEO_accession)
  
  #Download data
  if (!file.exists(exdir)){
    
    getGEOSuppFiles(GEO_accession)
    tarfile = paste0(GEO_accession, "/", GEO_accession, "_RAW.tar")
    untar(tarfile, exdir = exdir)
    
    tarfile = paste0(wd, "/", GEO_accession, "/", GEO_accession, "_RAW.tar")
    
    untar(tarfile, exdir = exdir)
    
  }
  
  #Get CEL files
  celfiles = list.files(paste0(exdir, "/"), pattern = "CEL", full.names = TRUE)
  
  if (!is.null(outliers)){
    
    outliers_select <- outliers[1]
    
    if (length(outliers) > 1){
      for (i in 1:(length(outliers)-1)){
        outliers_select <- paste(outliers_select, outliers[1+i], sep = "|")
      }
      
    }
    
    celfiles = celfiles[-grep(outliers_select, celfiles)]
  }
  
  
  #Get CDF 
  
  if (annotation == "ense") {
    cdf <- paste0(chiptype,"ense", version, "cdf")
    
  }
  
  if (annotation == "enst") {
    if (robust == TRUE){
      cdf <- paste0(chiptype, "enst", "robust", version, "cdf")
    }
    if (robust == FALSE){
      cdf <- paste0(chiptype, "enst", "all", version, "cdf")
    }
  }
  
  
  if(!(cdf %in% pkg)) {
    
    if (annotation == "ense") {
      cdf.file = paste0(chiptype, "ense", version)
    }
    
    if (annotation == "enst") {
      
      if (robust == FALSE) {
        cdf.file = paste0(chiptype, "enst", "all", version)
      }
      
      if (robust == TRUE) {
        cdf.file = paste0(chiptype, "enst", "robust", version)
      }
      
    }
    
    
    if (exists(cdf.file)) {
      make.cdf.package(paste0(cdf.file, ".cdf"), species = "Homo_sapiens")
      install.packages(paste0(cdf.file, "cdf/"), repos = NULL, type = "source")
      
    }
    
    if (!exists(cdf.file)) {
      print("Could not find correct CDF package or file. Please use the makeCDFfile function.")
      cdf = NULL
    }
    
  }
  
  
  
  #read microarray data
  if (!is.null(celfiles) & !is.null(cdf)) {
    data = ReadAffy(filenames=celfiles, cdfname = cdf)
    
    return(data)
  }
}


###################################################################################################################################

#pca.plot

###################################################################################################################################


pca.plot <- function(data.PC, meta, hpc = 1, vpc = 2) {
  
  
  samples = as.data.frame(rownames(data.PC$x))
  colnames(samples) <- "cel_names"
  
  samples_group <- samples %>% fuzzyjoin::fuzzy_inner_join(meta, by = c("cel_names" = "GEO ID"), match_fun = str_detect)
  
  Groups <- as.factor(samples_group[,3])
  Sample <- rownames(data.PC$x)
  
  
  X = as.data.frame(data.PC$x)[,hpc]
  Y = as.data.frame(data.PC$x)[,vpc]
  
  pca <- ggplot(data = as.data.frame(data.PC$x), aes(x = X, y = Y, col = Groups, key = Sample)) +
    geom_point() +
    ggtitle("PCA plot") +
    xlab(paste0("PC", hpc)) +
    ylab(paste0("PC", vpc)) +
    theme_light()
  
  return(ggplotly(pca))
  
}



###################################################################################################################################

#get_contrasts

###################################################################################################################################

get_contrasts <- function(meta) {
  
  levels1 <- unique(meta[,2])
  
  contrast <- paste(make.names(levels1)[2], make.names(levels1)[1], sep = " - ")
  
  for (m in 1:(length(levels1)-1)){
    for (n in (m+1):length(levels1)) {
      
      if (grepl("non|healthy|control", levels1[m])) {
        contrast <- c(contrast, paste(make.names(levels1)[n], make.names(levels1)[m], sep = " - "))
      }
      
      if (!grepl("non|healthy|control", levels1[m])){
        contrast <- c(contrast, paste(make.names(levels1)[m], make.names(levels1)[n], sep = " - "))
      }
    }
  }
  
  return(sort(contrast[-1]))
  
}

###################################################################################################################################

#auto_contrast

###################################################################################################################################

auto_contrasts <- function(meta) {
  
  levels1 <- unique(meta[,2])
  
  contrast <- paste(make.names(levels1)[2], make.names(levels1)[1], sep = " - ")
  
  for (m in 1:(length(levels1)-1)){
    for (n in (m+1):length(levels1)) {
      
      if (grepl("non|healthy|control", levels1[m])) {
        contrast <- c(contrast, paste(make.names(levels1)[n], make.names(levels1)[m], sep = " - "))
      }
      
      if (!grepl("non|healthy|control", levels1[m])){
        contrast <- c(contrast, paste(make.names(levels1)[m], make.names(levels1)[n], sep = " - "))
      }
    }
  }
  
  
  contrast <- contrast[-1]
  
  contrast <- str_replace_all(contrast, "_", "\\.") 
  
  contrast1 <- contrast
  
  
  for (k in 1:length(contrast)) {
    
    
    a <- strsplit(strsplit(contrast[k], split= " - ")[[1]][1], split = "\\.")
    a <- a[[1]]
    b <- strsplit(strsplit(contrast[k], split= " - ")[[1]][2], split = "\\.")
    b <- b[[1]]
    
    c <- setdiff(a,b)
    d <- setdiff(b,a)
    
    length_c <- length(c)
    length_d <- length(d)
    
    if (length_c > 1) {
      
      for (i in 1:(length_c-1)) {
        for (j in (i+1):(length_c)) {
          
          e <- paste(c[i], c[j], sep = ".")
          
          if (length(grep(e, make.names(unique(meta[,2])))) == length(match(grep(c[i], make.names(unique(meta[,2]))), grep(c[j], contrast)))) {
            length_c = length_c - 1
          }
        }
      }
    }
    
    
    if (length_d > 1) {
      
      for (i in 1:(length(d)-1)) {
        for (j in (i+1):(length(d))) {
          
          f <- paste(d[i], d[j], sep = ".")
          
          if (length(grep(f, make.names(unique(meta[,2])))) == length(match(grep(d[i], make.names(unique(meta[,2]))), grep(d[j], contrast)))) {
            length_d = length_d - 1
          }
          
          
        }
      }
    }
    
    if (length_c <= 1 & length_d <= 1) {
    }
    
    if (length_c > 1 | length_d > 1) {
      contrast1 <- contrast1[-which(contrast1 == contrast[k])]
    }
    
  }
  
  return(contrast1)
  
}


###################################################################################################################################

#diff_expr

###################################################################################################################################


diff_expr <- function(data.expr, meta, comparisons) {
  
  
  #get grouping variable in correct order
  
  ph = as.data.frame(colnames(data.expr))
  colnames(ph) <- "cel_names"
  
  ph1 <- ph %>% fuzzyjoin::fuzzy_inner_join(meta, by = c("cel_names" = "GEO ID"), match_fun = str_detect)
  
  
  #model design
  
  source <- ph1[,3]
  levels <- unique(ph1[,3])
  
  f <- factor(source,levels=levels)
  
  design <- model.matrix(~ 0 + f)
  colnames(design) <- make.names(levels)
  
  
  #fit linear model
  
  data.fit <- lmFit(data.expr,design)
  
  
  #get contrasts
  contrast <- comparisons
  
  
  
  #make top table for each comparison
  top.table <- list()
  
  for (i in contrast) {
    contrast.matrix <- makeContrasts(contrasts = i,levels=design)
    
    data.fit.con <- contrasts.fit(data.fit,contrast.matrix)
    
    data.fit.eb <- eBayes(data.fit.con)
    
    top.table[[i]] <- topTable(data.fit.eb, sort.by = "P", n = Inf)
    
  }
  
  return(top.table)
}

###################################################################################################################################

#transcript_selection

###################################################################################################################################


transcript_selection <- function(top.table, genes = NULL, transcripts = NULL, genome_wide = TRUE) {
  mart <- useMart("ENSEMBL_MART_ENSEMBL")
  mart <- useDataset("hsapiens_gene_ensembl", mart)
  
  if (genome_wide == FALSE){
    if (!is.null(genes)){
      gene_info <- getBM(
        mart=mart,
        attributes=c(
          "ensembl_transcript_id",
          "ensembl_gene_id"
        ),
        filter = "ensembl_gene_id",
        values = genes, uniqueRows=TRUE)
      
      select.top.table <- list()
      
      for (i in names(top.table)){
        transcript.top.table <- top.table[[i]]
        transcript.top.table <- transcript.top.table[,1:4]
        transcript.top.table <- cbind(rownames(transcript.top.table), transcript.top.table)
        rownames(transcript.top.table) <- NULL
        colnames(transcript.top.table) <- c("Ensembl.ID", "logFC", "Average expression", "t", "P.value")
        
        transcript.top.table <- inner_join(gene_info, transcript.top.table, by = c("ensembl_transcript_id" = "Ensembl.ID"))
        
        if (nrow(transcript.top.table) > 0) {
          transcript.top.table <- cbind(transcript.top.table, p.adjust(transcript.top.table[,"P.value"], method = "fdr"))
          colnames(transcript.top.table) <- c("Ensembl.transcript.ID", "Ensembl.gene.ID", "logFC", "Average expression", "t", "P.value", "FDR")
          transcript.top.table <- arrange(transcript.top.table, P.value)
        }
        
        select.top.table[[i]] <- transcript.top.table
        
        
      }
    }
    
    
    #based on list of transcripts
    
    if (!is.null(transcripts)) {
      
      gene_info <- getBM(
        mart=mart,
        attributes=c(
          "ensembl_transcript_id",
          "ensembl_gene_id"
        ),
        filter = "ensembl_transcript_id",
        values = transcripts, uniqueRows=TRUE)
      
      
      select.top.table <- list()
      
      for (i in names(top.table)){
        transcript.top.table <- top.table[[i]]
        transcript.top.table <- transcript.top.table[,1:4]
        transcript.top.table <- cbind(rownames(transcript.top.table), transcript.top.table)
        rownames(transcript.top.table) <- NULL
        colnames(transcript.top.table) <- c("Ensembl.ID", "logFC", "Average expression", "t", "P.value")
        
        transcript.top.table <- inner_join(gene_info, transcript.top.table, by = c("ensembl_transcript_id" = "Ensembl.ID"))
        transcript.top.table <- cbind(transcript.top.table, p.adjust(transcript.top.table[,"P.value"], method = "fdr"))
        colnames(transcript.top.table) <- c("Ensembl.transcript.ID", "Ensembl.gene.ID", "logFC", "Average expression", "t", "P.value", "FDR")
        
        transcript.top.table <- arrange(transcript.top.table, P.value)
        select.top.table[[i]] <- transcript.top.table
      }
    }
  }
  
  
  
  if (genome_wide == TRUE){
    
    gene_info <- getBM(
      mart=mart,
      attributes=c(
        "ensembl_transcript_id",
        "ensembl_gene_id"
      ),
      filter = "ensembl_transcript_id",
      values = rownames(top.table[[1]]), uniqueRows=TRUE)
    
    select.top.table <- list()
    
    for (i in names(top.table)){
      transcript.top.table <- top.table[[i]]
      transcript.top.table <- transcript.top.table[,1:4]
      transcript.top.table <- cbind(rownames(transcript.top.table), transcript.top.table)
      rownames(transcript.top.table) <- NULL
      colnames(transcript.top.table) <- c("Ensembl.ID", "logFC", "Average expression", "t", "P.value")
      
      transcript.top.table <- inner_join(gene_info, transcript.top.table, by = c("ensembl_transcript_id" = "Ensembl.ID"))
      transcript.top.table <- cbind(transcript.top.table, p.adjust(transcript.top.table[,"P.value"], method = "fdr"))
      colnames(transcript.top.table) <- c("Ensembl.transcript.ID", "Ensembl.gene.ID", "logFC", "Average expression", "t", "P.value", "FDR")
      
      transcript.top.table <- arrange(transcript.top.table, P.value)
      select.top.table[[i]] <- transcript.top.table
    }
  }
  
  
  return(select.top.table)
}



###################################################################################################################################

#exon_selection

###################################################################################################################################


exon_selection <- function(top.table, annotated, genes = NULL, transcripts = NULL, genome_wide = TRUE, unique_exons = TRUE) {
  
  annotation1 <- annotated[complete.cases(annotated),]
  
  if (genome_wide == FALSE) {
    if (!is.null(genes)){
      genes <- as.data.frame(genes)
      annotation1 <- inner_join(annotation1, genes, by = c("Gene" = "genes"))
    }
    
    if (!is.null(transcripts)){
      transcripts <- as.data.frame(transcripts)
      annotation1 <- inner_join(annotation1, transcripts, by = c("Transcript" = "transcripts"))
    }
  }
  
  
  
  if (unique_exons == TRUE){
    freq_table <- as.data.frame(table(unlist(annotation1[,"Exon.group"])))
    
    exons_unique_for_one_transcripts <- freq_table %>%
      filter(Freq == 1)
    
    unique_exons <- annotation1 %>%
      inner_join(exons_unique_for_one_transcripts, by = c("Exon.group" = "Var1"))
    
    annotation1 <- unique_exons[,-6]
    
  }
  
  
  
  select.top.table <- list()
  
  for (i in names(top.table)){
    transcript.top.table <- top.table[[i]]
    transcript.top.table <- transcript.top.table[,1:4]
    
    transcript.top.table <- cbind(rownames(transcript.top.table), transcript.top.table)
    rownames(transcript.top.table) <- NULL
    colnames(transcript.top.table) <- c("Ensembl.ID", "logFC", "Average expression", "t", "P.value")
    
    transcript.top.table <- inner_join(annotation1, transcript.top.table, by = c("Probe.set" = "Ensembl.ID"))
    transcript.top.table <- transcript.top.table[, c("Exon.group", "Transcript", "Gene", "logFC", "Average expression", "t", "P.value")]
    
    transcript.top.table <- cbind(transcript.top.table, p.adjust(transcript.top.table[,"P.value"], method = "fdr"))
    colnames(transcript.top.table) <- c("Ensembl.exon.ID", "Ensembl.transcript.ID", "Ensembl.gene.ID", "logFC", "Average expression", "t", "P.value", "FDR")
    transcript.top.table <- arrange(transcript.top.table, P.value)
    
    select.top.table[[i]] <- transcript.top.table
    
  }
  
  return(select.top.table)
  
}

###################################################################################################################################

#makeBoxplots

###################################################################################################################################

makeBoxplots <- function(contrast, meta, data.expr, annotated, gene = NULL, transcripts = NULL, unique_exons = TRUE, genome_wide = FALSE) {
  
  annotation1 <- annotated[complete.cases(annotated),]
  
  if (genome_wide == FALSE) {
    if (!is.null(gene)){
      gene <- as.data.frame(gene)
      annotation1 <- inner_join(annotation1, gene, by = c("Gene" = "gene"))
    }
    
    if (!is.null(transcripts)){
      transcripts <- as.data.frame(transcripts)
      annotation1 <- inner_join(annotation1, transcripts, by = c("Transcript" = "transcripts"))
    }
  }
     
  if (unique_exons == TRUE){
    freq_table <- as.data.frame(table(unlist(annotation1[,"Exon.group"])))
    
    exons_unique_for_one_transcripts <- freq_table %>%
      filter(Freq == 1)
    
    unique_exons <- annotation1 %>%
      inner_join(exons_unique_for_one_transcripts, by = c("Exon.group" = "Var1"))
    
    annotation1 <- unique_exons[,-6]
    
  }
    
  #Get samples from contrasts
  contrast <- strsplit(contrast, split = " - ")
  contrast <- rev(contrast[[1]])
  
  meta1 <- meta
  meta1$Grouping <- make.names(meta1$Grouping)
  meta1 <- inner_join(meta1, as.data.frame(contrast), by = c("Grouping" = "contrast"))
  
  
  #select samples of contrast from the expression matrix
  samples <- fuzzyjoin::fuzzy_inner_join(as.data.frame(colnames(data.expr)), meta1, by = c("colnames(data.expr)" = "GEO ID"), match_fun = str_detect)
  samples$Grouping <- factor(samples$Grouping, levels = contrast)
  
  data.expr1 <- data.expr[,samples[,1]]
  data.expr1 <- as.data.frame(data.expr1)
  data.expr1 <- cbind(rownames(data.expr1), data.expr1)
  rownames(data.expr1) <- NULL
  colnames(data.expr1) <- c("Probe.set", colnames(data.expr1)[-1])
  data.expr1 <- inner_join(annotation1[,1:2], data.expr1, by = c("Probe.set" = "Probe.set"))
  data.expr1 <- data.expr1[,-1]
  
  
  #arrange data in correct order for making the figures
  
  
  datafigure1 <- as.data.frame(t(data.expr1[,-1]))
  colnames(datafigure1) <- data.expr1[,1]
  datafigure1 <- cbind(rownames(datafigure1), datafigure1)
  rownames(datafigure1) <- NULL
  colnames(datafigure1) <- c("Samples", colnames(datafigure1)[-1])
  datafigure1 <- datafigure1[,colnames(datafigure1)[!duplicated(colnames(datafigure1))]]
  datafigure1 <- arrange(datafigure1, Samples)
  
  final <- pivot_longer(datafigure1, data.expr1[,1])
  colnames(final) <- c("Samples", "Exons", "logExpr")
  
  
  final <- inner_join(final, samples, by = c("Samples" = "colnames(data.expr)"))
  
  ggplot(data = final, aes(x = Grouping, y = logExpr, fill = Exons)) +
    geom_boxplot() +
    theme_minimal() +
    theme(panel.grid.major.x = element_blank()) +
    xlab("") +
    ylab("logExpr")
}



###################################################################################################################################

#makeBoxplots1

###################################################################################################################################

makeBoxplots1 <- function(contrast, meta, data.expr, select.top.table) {
  #Get samples from contrasts
  contrast <- strsplit(contrast, split = " - ")
  contrast <- rev(contrast[[1]])
  
  meta1 <- meta
  meta1$Grouping <- make.names(meta1$Grouping)
  meta1 <- inner_join(meta1, as.data.frame(contrast), by = c("Grouping" = "contrast"))
  
  
  #select samples of contrast from the expression matrix
  samples <- fuzzyjoin::fuzzy_inner_join(as.data.frame(colnames(data.expr)), meta1, by = c("colnames(data.expr)" = "GEO ID"), match_fun = str_detect)
  samples$Grouping <- factor(samples$Grouping, levels = contrast)
  
  data.expr1 <- data.expr[,samples[,1]]
  data.expr1 <- as.data.frame(data.expr1)
  data.expr1 <- cbind(rownames(data.expr1), data.expr1)
  rownames(data.expr1) <- NULL
  colnames(data.expr1) <- c("Transcript", colnames(data.expr1)[-1])
  data.expr1 <- inner_join(data.expr1, as.data.frame(select.top.table[[1]]$Ensembl.transcript.ID), by = c("Transcript" = "select.top.table[[1]]$Ensembl.transcript.ID"))
  
  
  #arrange data in correct order for making the figures
  
  
  datafigure1 <- as.data.frame(t(data.expr1[,-1]))
  colnames(datafigure1) <- data.expr1[,1]
  datafigure1 <- cbind(rownames(datafigure1), datafigure1)
  rownames(datafigure1) <- NULL
  colnames(datafigure1) <- c("Samples", colnames(datafigure1)[-1])
  datafigure1 <- arrange(datafigure1, Samples)
  
  final <- pivot_longer(datafigure1, data.expr1[,1])
  colnames(final) <- c("Samples", "Transcripts", "logExpr")
  
  
  final <- inner_join(final, samples, by = c("Samples" = "colnames(data.expr)"))
  
  ggplot(data = final, aes(x = Grouping, y = logExpr, fill = Transcripts)) +
    geom_boxplot() +
    theme_minimal() +
    theme(panel.grid.major.x = element_blank()) +
    xlab("") +
    ylab("logExpr")
}


###################################################################################################################################

#probe mapping enst

###################################################################################################################################

probemapping.enst <- function(chiptype, organism, version, gene){
  
  #Get mapping file
  mapping_file <- paste0(chiptype, "_", toupper(substring(organism, 1,1)), substring(organism, 2), "_ENST_mapping.txt")
  
  
  if (!file.exists(mapping_file)){
    file_download <- paste0("http://mbni.org/customcdf/", version, ".0.0/enst.download/", chiptype, "_",
                            toupper(substring(organism, 1,1)), substring(organism, 2), "_ENST_", version, ".0.0.zip")
    
    download.file(file_download, "probe_information_file.zip")
    unzip("probe_information_file.zip", mapping_file)
    unlink("probe_information_file.zip")
    
  }
  
  mapping <- read.table(mapping_file, sep = "\t", header = TRUE)
  mapping$probe <- paste(mapping$Probe.X, mapping$Probe.Y, sep = "_")
  mapping$Probe.Set.Name <- str_replace_all(mapping$Probe.Set.Name, "\\.._at", "")
  mapping$Probe.Set.Name <- str_replace_all(mapping$Probe.Set.Name, "\\..._at", "")
  
  
  #Get transcripts of gene
  mart <- useMart("ENSEMBL_MART_ENSEMBL")
  mart <- useDataset("hsapiens_gene_ensembl", mart)
  
  gene_info <- getBM(
    mart=mart,
    attributes=c(
      "ensembl_gene_id",
      "ensembl_transcript_id"
    ),
    filter = "ensembl_gene_id",
    values = gene, uniqueRows=TRUE)
  
  
  
  #Link transcripts to probe
  mapping <- arrange(mapping, Chr.From)
  mapping$probe <- factor(mapping$probe, levels = unique(mapping$probe))
  transcripts_probe_all <- mapping[, c("Probe.Set.Name", "probe", "Chr", "Chr.Strand", "Chr.From")]
  gene_probes_all <- inner_join(transcripts_probe_all, gene_info, by = c("Probe.Set.Name" = "ensembl_transcript_id"))
  
  #Filter for probes that are unique for 1 transcript
  probefreq <- as.data.frame(table(unlist(gene_probes_all$probe)))
  
  probefreq_filtered <- filter(probefreq, Freq == 1)
  probes_filtered <- inner_join(gene_probes_all, probefreq_filtered, by = c("probe" = "Var1"))
  probes_filtered <- cbind(probes_filtered, rep("Yes", nrow(probes_filtered)))
  colnames(probes_filtered) <- c("Transcript", "Probe", "Chr", "Chr.Strand", "Chr.From", "Gene", "Freq", "Unique.probe")
  
  probefreq_remainder <- filter(probefreq, Freq > 1)
  probes_remainder <- inner_join(gene_probes_all, probefreq_remainder, by = c("probe" = "Var1"))
  probes_remainder <- cbind(probes_remainder, rep("No", nrow(probes_remainder)))
  colnames(probes_remainder) <- c("Transcript", "Probe", "Chr", "Chr.Strand", "Chr.From", "Gene", "Freq", "Unique.probe")
  
  total <- rbind(probes_filtered, probes_remainder)
  number <- (4*length(unique(total$Transcript)) + length(unique(total$Probe)))/5
  
  #Make plot
  gg <- ggplot(total, aes(y = Transcript, x = Probe, fill = Unique.probe, text = paste(" Chr:", Chr, "<br>","Chr.Strand:", Chr.Strand, "<br>", "Chr.From:", Chr.From))) +
    geom_point(size = ((5*14)/number), shape = 22, color = "black") +
    theme_light() +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
    ggtitle("Probe mapping") +
    ylab("Transcripts") +
    xlab("Probes") +
    labs(fill = "Unique?")
  
  return(gg)
}

###################################################################################################################################

#probe mapping ense

###################################################################################################################################

probemapping.ense <- function(chiptype, organism, version, gene, annotated){
  
  #Get mapping file
  mapping_file <- paste0(chiptype, "_", toupper(substring(organism, 1,1)), substring(organism, 2), "_ENSE_mapping.txt")
  
  
  if (!file.exists(mapping_file)){
    file_download <- paste0("http://mbni.org/customcdf/", version, ".0.0/ense.download/", chiptype, "_",
                            toupper(substring(organism, 1,1)), substring(organism, 2), "_ENSE_", version, ".0.0.zip")
    
    download.file(file_download, "probe_information_file.zip")
    unzip("probe_information_file.zip", mapping_file)
    unlink("probe_information_file.zip")
    
  }
  
  mapping <- read.table(mapping_file, sep = "\t", header = TRUE)
  mapping$probe <- paste(mapping$Probe.X, mapping$Probe.Y, sep = "_")
  mapping <- arrange(mapping, Chr.From)
  mapping$probe <- factor(mapping$probe, levels = unique(mapping$probe))
  mapping <- mapping[, c("Probe.Set.Name", "probe", "Chr", "Chr.Strand", "Chr.From")]
  
  #Get exons of gene
  annotation1 <- annotated
  gene <- as.data.frame(gene)
  annotation1 <- inner_join(annotation1, gene, by = c("Gene" = "gene"))
  
  gene_probes_all <- inner_join(annotation1, mapping, by = c("Probe.set" = "Probe.Set.Name"))
  
  number <- (4*length(unique(gene_probes_all$Exon.group)) + length(unique(gene_probes_all$probe)))/5
  
  #Make plot
  gg <- ggplot(gene_probes_all, aes(y = Exon.group, x = probe, text = paste(" Chr:", Chr, "<br>","Chr.Strand:", Chr.Strand, "<br>", "Chr.From:", Chr.From))) +
    geom_point(size = ((5*14)/number), shape = 22, fill = "blue", color = "black") +
    theme_light() +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
    ggtitle("Probe mapping") +
    ylab("") +
    xlab("Probes") +
    labs(fill = "Unique?")
  
  return(gg)
}

###################################################################################################################################

#Flat2CDF

###################################################################################################################################

flat2Cdf<-function(file,chipType,tags=NULL,rows=2560,cols=2560,verbose=10,xynames=c("X","Y"),
                   gcol=5,ucol=6,splitn=4,col.class=c("integer","character")[c(1,1,1,2,2,2)],...) {
  split.quick<-
    function(r,ucol,splitn=3,verbose=TRUE) {
      rn3<-substr(r[,ucol],1,splitn)
      split.matrix<-split.data.frame
      rr<-split(r,factor(rn3))
      if (verbose) cat(" split into",length(rr),"initial chunks ...")
      rr<-unlist(lapply(rr,FUN=function(u) split(u,u[,ucol])),recursive=FALSE)
      if (verbose) cat(" unwrapped into",length(rr),"chunks ...")
      names(rr)<-substr(names(rr),splitn+2,nchar(rr))
      rr
    }
  
  if (verbose) cat("Reading TXT file ...")
  file<-read.table(file,header=TRUE,colClasses=col.class,stringsAsFactors=FALSE,comment.char="",...)
  if (verbose) cat(" Done.\n")
  
  if (verbose) cat("Splitting TXT file indices into units ...")
  #gxys<-split.quick(file,ucol,splitn)
  gxysInd<-split(seq_len(nrow(file)),file[,ucol])
  if (verbose) cat(" Done.\n")
  
  l<-vector("list",length(gxysInd))
  if (verbose) cat("Creating structure for",length(gxysInd),"units (dot=250):\n")
  for(i in  1:length(gxysInd)) {
    thisTab <- file[ gxysInd[[i]], ]
    sp<-split(thisTab, factor(thisTab[,gcol]))
    #sp<-split(gxys[[i]],factor(gxys[[i]][,gcol]))
    e<-vector("list",length(sp))
    for(j in 1:length(sp)) {
      np<-nrow(sp[[j]])
      e[[j]]<-list(x=sp[[j]][,xynames[1]],y=sp[[j]][,xynames[2]],pbase=rep("A",np),tbase=rep("T",np),atom=0:(np-1),indexpos=0:(np-1),
                   groupdirection="sense",natoms=np,ncellsperatom=1)
    }
    names(e)<-names(sp)
    #l[[i]]<-list(unittype=1,unitdirection=1,groups=e,natoms=nrow(gxys[[i]]),ncells=nrow(gxys[[i]]),ncellsperatom=1,unitnumber=i)
    l[[i]]<-list(unittype=1,unitdirection=1,groups=e,natoms=nrow(thisTab),ncells=nrow(thisTab),ncellsperatom=1,unitnumber=i)
    if (verbose) { if(i %% 250==0) cat("."); if(i %% 5000==0) cat("(",i,")\n",sep="") }
  }
  rm(file,e,sp,thisTab); gc()
  cat("\n")
  #names(l)<-names(gxys)
  names(l)<-names(gxysInd)  
  rm(gxysInd); gc()
  if(!is.null(tags) && tags!="") filename<-paste(chipType,tags,sep=",")
  else filename<-chipType
  filename<-paste(filename,"cdf",sep=".")
  hdr<-list(probesets=length(l),qcprobesets=0,reference="",chiptype=chipType,filename=filename,
            nqcunits=0,nunits=length(l),rows=rows,cols=cols,refseq="",nrows=rows,ncols=cols)
  require(affxparser)
  writeCdf(hdr$filename, cdfheader=hdr, cdf=l, cdfqc=NULL, overwrite=TRUE, verbose=verbose)
  invisible(list(cdfList=l,cdfHeader=hdr))
}



###################################################################################################################################

#MakeCDFfile

###################################################################################################################################
MakeCDFfile <- function(chiptype, organism, version, transcript = TRUE, robust = TRUE){
  
  if (transcript == FALSE) {
    ##############
    #Exon
    ##############
    mapping_file <- paste0(chiptype, "_", toupper(substring(organism, 1,1)), substring(organism, 2), "_ENSE_mapping.txt")
    mapping <- read.table(mapping_file, sep = "\t", header = TRUE)
    mapping$probe <- paste(mapping$Probe.X, mapping$Probe.Y, sep = "_")
    
    #Exon groups
    annotated <- read.table(file = paste(chiptype, organism, version, "OfficialNewProbeSets2.txt", sep = "_"), header = TRUE)
    probesets <- as.data.frame(unique(annotated$Probe.set))
    colnames(probesets) <- "Probe.Set.Name"
    
    mapping_ense <- inner_join(probesets, mapping, by = c("Probe.Set.Name" = "Probe.Set.Name"))
    
  }
  
  if (transcript == TRUE) {
    ##############
    #Transcripts
    ##############
    mapping_file <- paste0(chiptype, "_", toupper(substring(organism, 1,1)), substring(organism, 2), "_ENST_mapping.txt")
    mapping <- read.table(mapping_file, sep = "\t", header = TRUE)
    mapping$probe <- paste(mapping$Probe.X, mapping$Probe.Y, sep = "_")
    mapping$Probe.Set.Name <- str_replace_all(mapping$Probe.Set.Name, "\\.._at", "")
    mapping$Probe.Set.Name <- str_replace_all(mapping$Probe.Set.Name, "\\..._at", "")
    
    #Filter for probes that are unique for 1 transcript
    probefreq <- as.data.frame(table(unlist(mapping$probe)))
    probefreq <- filter(probefreq, Freq == 1)
    mapping_enst_all <- inner_join(mapping, probefreq, by = c("probe" = "Var1"))
    
    #Filter for probe sets with at least 3 probes
    transcriptfreq <- as.data.frame(table(unlist(mapping_enst_all$Probe.Set.Name)))
    transcriptfreq <- filter(transcriptfreq, Freq > 2)
    mapping_enst_robust <- inner_join(mapping_enst_all, transcriptfreq, by = c("Probe.Set.Name" = "Var1"))
    
    mapping_enst_robust <- mapping_enst_robust[,1:8]
    mapping_enst_all <- mapping_enst_all[,1:8]
    
  }
  
  #Get mapping file
  mapping_file_gene <- paste0(chiptype, "_", toupper(substring(organism, 1,1)), substring(organism, 2), "_ENSG_mapping.txt")
  mapping_file_transcript <- paste0(chiptype, "_", toupper(substring(organism, 1,1)), substring(organism, 2), "_ENST_mapping.txt")
  mapping_file_exon <- paste0(chiptype, "_", toupper(substring(organism, 1,1)), substring(organism, 2), "_ENSE_mapping.txt")
  
  if (!file.exists(mapping_file_gene)){
    file_download <- paste0("http://mbni.org/customcdf/", version, ".0.0/ensg.download/", chiptype, "_",
                            toupper(substring(organism, 1,1)), substring(organism, 2), "_ENSG_", version, ".0.0.zip")
    
    download.file(file_download, "probe_information_file.zip")
    unzip("probe_information_file.zip", mapping_file_gene)
    unlink("probe_information_file.zip")
    
  }
  
  if (!file.exists(mapping_file_transcript)){
    file_download <- paste0("http://mbni.org/customcdf/", version, ".0.0/enst.download/", chiptype, "_",
                            toupper(substring(organism, 1,1)), substring(organism, 2), "_ENST_", version, ".0.0.zip")
    
    download.file(file_download, "probe_information_file.zip")
    unzip("probe_information_file.zip", mapping_file_transcript)
    unlink("probe_information_file.zip")
    
  }
  
  if (!file.exists(mapping_file_exon)){
    file_download <- paste0("http://mbni.org/customcdf/", version, ".0.0/ense.download/", chiptype, "_",
                            toupper(substring(organism, 1,1)), substring(organism, 2), "_ENSE_", version, ".0.0.zip")
    
    download.file(file_download, "probe_information_file.zip")
    unzip("probe_information_file.zip", mapping_file_exon)
    unlink("probe_information_file.zip")
    
  }
  
  
  #Read mapping file
  mapping_gene <- read.table(mapping_file_gene, sep = "\t", header = TRUE)
  mapping_transcript <- read.table(mapping_file_transcript, sep = "\t", header = TRUE)
  mapping_exon <- read.table(mapping_file_exon, sep = "\t", header = TRUE)
  
  
  #Combine mapping files
  mapping_total <- rbind(mapping_gene, mapping_transcript, mapping_exon)
  mapping_total$probe <- paste(mapping_total$Probe.X, mapping_total$Probe.Y, sep = "_")
  mapping_total <- mapping_total[!duplicated(mapping_total$probe),]
  
  
  if (transcript == FALSE) {
    #Ense
    nonsense.probes <- anti_join(mapping_total, mapping_ense, by = c("probe" = "probe"))
    nonsense.probes$Probe.Set.Name <- rep("nonsense", nrow(nonsense.probes))
    
    mapping <- rbind(mapping_ense, nonsense.probes)
    
    cdf = paste0(chiptype, "ense", "", version)
    
    
  }
  
  if (transcript == TRUE) {
    
    if (robust == FALSE) {
      #Enst_all
      nonsense.probes <- anti_join(mapping_total, mapping_enst_all, by = c("probe" = "probe"))
      nonsense.probes$Probe.Set.Name <- rep("nonsense", nrow(nonsense.probes))
      
      mapping <- rbind(mapping_enst_all, nonsense.probes)
      
      cdf = paste0(chiptype, "ense", "all", version)
      
    }
    
    if (robust == TRUE){
      #Enst_robust
      nonsense.probes <- anti_join(mapping_total, mapping_enst_robust, by = c("probe" = "probe"))
      nonsense.probes$Probe.Set.Name <- rep("nonsense", nrow(nonsense.probes))
      
      mapping <- rbind(mapping_enst_robust, nonsense.probes)
      
      cdf = paste0(chiptype, "ense", "robust", version)
    }
    
  }
  
  final_all <- mapping[, c(5,6,1)]
  colnames(final_all) <- c("X", "Y", "Group_ID")
  final_all$Unit_ID <- final_all$Group_ID
  final_all$Probe_ID <- seq(1:nrow(final_all))
  final_all$Probe_Sequence <- mapping$probe
  final_all <- final_all[, c("Probe_ID", "X", "Y", "Probe_Sequence", "Group_ID", "Unit_ID")]
  
  write.table(final_all, file = paste(chiptype, organism, version, "final_all.txt", sep = "_"), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
  
  if (chiptype == "hugene10st"){
    rows = 1050
  }
  
  
  if (chiptype == "huex10st"){
    rows = 2560
  }
  
  flat2Cdf(file = paste(chiptype, organism, version, "final_all.txt", sep = "_"), 
           chipType = cdf, 
           rows = rows, 
           cols = rows)
  
  
}

