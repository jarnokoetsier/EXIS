

#############   ####        ####   ####      #########
#############    ####      ####    ####    ##########
#####             ###########      ####   #####
#####               #######        ####   #####
#############     ###########      ####    ######## 
#############    ####      ####    ####      ########
#####           ####        ####   ####         ######
#####                                            #####
#####################################################
#####################################################
#Isoform expression analysis


#CREATED BY JARNO KOETSIER

###############################################################################################################################

#INSTRUCTIONS WORKFLOW SCRIPT

###############################################################################################################################



#1. Source run the global R script every time you use the workflow.

#2. Make sure that the exon/transcript-specific CDF file is in your working directory 
#   OR has been installed as a CDF package previously.

#3. When you would like to use the exon-specific probe set definition, 
#   also make sure to place the exon annotation file in your working directory



################################################################################################################################
#get GEO
################################################################################################################################



#Enter GEO accession number
GEO = "GSE36980"
GEO = "GSE37264"

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

#Get Brainarray version
version = 25

#Annotation: 
#  "ense" for exon-specific probe sets
#  "enst" for transcript-specific probe sets
annotation = "ense"


#Make affyBatch
data1 <- readcels(gset = gset, 
                  chiptype = chiptype, 
                  organism = organism, 
                  version = version,
                  annotation = annotation,
                  robust = FALSE,
                  outliers = "GSM4764672")

data1 <- readcels(gset = gset, 
                  chiptype = chiptype, 
                  organism = organism, 
                  version = version,
                  annotation = annotation,
                  robust = FALSE,
                  outliers = NULL)


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
#Quality control
################################################################################################################################

#Boxplots

##Raw data
par(mar=c(8,4,2,1))
boxplot(data1,which='pm', col = "red", las = 2, main = "Boxplot of raw data", ylab = "log2 intensity")

##Normalized data
par(mar=c(10,2,1,1))
boxplot(data.expr, col = "blue", las = 2, main = "Boxplot of normalized data", ylab = "log2 intensity")




#Density plots

##Raw data
par(mar=c(5,4,2,1))
hist(data1,lwd=2,which='pm',ylab='Density',xlab='log2 intensity',main='Density plot of raw data')


##Normalized data
par(mar=c(5,4,2,1))
plotDensity(data.expr,lwd=2,ylab='Density',xlab='log2 intensity',main='Density plot of normalized data')



#PCA plot
data.PC <- prcomp(t(data.expr,scale.=TRUE))

pca.plot(data.PC = data.PC,
         meta = meta,
         hpc = 1,           #horizontal axis (1 = PC1 on horizontal axis, 2 = PC2 on horizontal axis, etc.)
         vpc = 2)           #vertical exis (2 = PC2 on vertical axis, 3 = PC3 on vertical axis, etc.)




################################################################################################################################
#Differential expression analysis
################################################################################################################################

#Get comparisons
comparisons <- get_contrasts(meta)

#NOTE: the get_contrast function retrieves all possible comparisons that are possible according to the meta table.

#Make top table using the meta data
#This function performs DE analysis for all possible comparisons and returns a list object.
top.table <- diff_expr(data.expr = data.expr, 
                       meta = meta, 
                       comparisons = comparisons)





#IF YOU WANT TO USE THE EXON-SPECIFIC PROBESET DEFINITION:
#IMPORTANT: you can only use run these codes if you used the "ense" annotation in the readcels function

################################################################################################################################

#Get annotation file from the working directory
annotated <- read.table(file = paste(chiptype, organism, version, "ExonAnnotation.txt", sep = "_"), header = TRUE)


#Select exons of interest AND annotate the exons with their corresponding transcript(s) and gene
select.top.table <- exon_selection(top.table = top.table, 
                                   annotated = annotated,
                                   genes = c("ENSG00000065989", 
                                             "ENSG00000184588", 
                                             "ENSG00000105650",
                                             "ENSG00000113448"),
                                   transcripts = NULL,
                                   genome_wide = FALSE,
                                   unique_exons = FALSE)
  




#IF YOU WANT TO USE THE TRANSCRIPT-SPECIFIC PROBESET DEFINITION:
#IMPORTANT: you can only use run these codes if you used the "enst" annotation in the readcels function

################################################################################################################################

#Select transcripts of interest AND annotate the transcripts with their corresponding gene
select.top.table <- transcript_selection(top.table = top.table,
                                         genes = c("ENSG00000065989", 
                                                   "ENSG00000184588", 
                                                   "ENSG00000105650",
                                                   "ENSG00000113448"),
                                         transcripts = NULL,
                                         genome_wide = FALSE)






################################################################################################################################
#Output plots
################################################################################################################################


#Boxplots

makeBoxplots(contrast = comparisons[1],
             meta = meta,
             data.expr = data.expr,
             annotated = annotated,
             genes = "ENSG00000113448",
             transcripts = NULL,
             exons = NULL,
             unique_exons = FALSE, 
             genome_wide = FALSE)
             
             
makeBoxplots1(contrast = comparisons[1],
              meta = meta,
              data.expr = data.expr,
              select.top.table = select.top.table)
           




#Probe mapping plot

#Transcript
probemapping <- probemapping.enst(chiptype = chiptype, 
                                  organism = organism,
                                  version = version, 
                                  gene = "ENSG00000113448")
ggplotly(probemapping)


#Exon
probemapping <- probemapping.ense(chiptype = chiptype, 
                                  organism = organism,
                                  version = version, 
                                  gene = "ENSG00000113448",
                                  annotated = annotated)
ggplotly(probemapping)



test <- select.top.table[[7]][,1:3]
PDE4A <- test[test$Ensembl.gene.ID == "ENSG00000065989",]
PDE4B <- test[test$Ensembl.gene.ID == "ENSG00000184588",]
PDE4C <- test[test$Ensembl.gene.ID == "ENSG00000105650",]
PDE4D <- test[test$Ensembl.gene.ID == "ENSG00000113448",]
test1 <- rbind(PDE4A, PDE4B, PDE4C, PDE4D)


final$Transcripts <- factor(final$Transcripts, levels = test1$Ensembl.transcript.ID)
final$Exons <- factor(final$Exons, levels = unique(test1$Ensembl.exon.ID))




library(writexl)
write_xlsx(select.top.table[[7]], "hippocampus.xlsx")
write_xlsx(select.top.table[[1]], "frontalcortex.xlsx")
write_xlsx(select.top.table[[14]], "temporalcortex.xlsx")

write_xlsx(select.top.table[[1]], "temporalcortex.xlsx")
