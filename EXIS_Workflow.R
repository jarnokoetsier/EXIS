


################################################################################################################################
#get GEO
################################################################################################################################



#Enter GEO accession number
GEO = "GSE36980"

#Download data set from GEO
gset <- getGEO(GEO, GSEMatrix =TRUE, getGPL = FALSE)




################################################################################################################################
#get meta
################################################################################################################################



#Get grouping variables from the data set
groups <- get_grouping(gset)

#Select grouping variable
groupingvar <- groups[,auto_group(gset)]

#NOTE: You can also select a pairing variable in case of a dependent study design.

#Make meta table
meta <- get_meta(gset = gset, grouping_column = groupingvar, pairing_column = NULL)




################################################################################################################################
#get data
################################################################################################################################

#Get chiptype
chiptype = get_chiptype(gset)

#Get organism
organism = get_organism(gset)

#NOTE: manual chiptype and/or organism selection is also possible.


#Annotation: 
#  "ense" for exon-specific probe sets
#  "enst" for transcript-specific probe sets
annotation = "ense"


#Make affyBatch
data1 <- readcels(gset = gset, 
                  chiptype = chiptype, 
                  organism = organism, 
                  version = 25,
                  annotation = annotation,
                  robust = TRUE,
                  outliers = "GSM4764672")



################################################################################################################################
#Quality control
################################################################################################################################


#Boxplots


#Density plots


#PCA plot




################################################################################################################################
#Normalization
################################################################################################################################

#Perform RMA normalization
data.rma <- affy::rma(data1, normalize = TRUE, background = TRUE)

#Retrieve probe set expression
data.expr <- exprs(data.rma)

#Remove the "nonsense" probe set. This probe set is only included for normalization purposes.
data.expr <- data.expr[rownames(data.expr) != "nonsense",]



################################################################################################################################
#Differential expression analysis
################################################################################################################################

#Get comparisons
comparisons <- get_contrasts(meta)

#NOTE: the get_contrast function retrieves all possible comparisons that are possible according to the meta table.

#Make top table
top.table <- diff_expr(data.expr = data.expr, 
                       meta = meta, 
                       comparisons = comparisons)



if (annotation == "ense"){
  #IF YOU WANT TO USE THE EXON-SPECIFIC PROBESET DEFINITION:
  
  #Get annotation file from the working directory
  annotated <- read.table(file = paste(chiptype, organism, 25, "OfficialNewProbeSets2.txt", sep = "_"), header = TRUE)
  
  
  #Select exons of interest AND annotate the exons with their corresponding transcript(s) and gene
  select.top.table <- exon_selection(top.table = top.table, 
                                     annotated = annotated,
                                     genes = NULL,
                                     transcripts = NULL,
                                     genome_wide = TRUE,
                                     unique_exons = FALSE)
  
}


if (annotation == "enst"){
  #IF YOU WANT TO USE THE TRANSCRIPT-SPECIFIC PROBESET DEFINITION:
  
  #Select transcripts of interest AND annotate the transcripts with their corresponding gene
  select.top.table <- transcript_selection(top.table = top.table,
                                           genes = "ENSG00000065989",
                                           transcripts = NULL,
                                           genome_wide = FALSE)
}


    

