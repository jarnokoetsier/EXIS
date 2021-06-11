

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

#EXIS 1.0.0 (June, 2021)

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
#In this step, the dataset is loaded from the Gene Expression Omnibus (GEO)



#Enter GEO accession number
#You can run the code with "GSE36980" for an example analysis.
GEO = "GSE36980"


#Download data set from GEO
gset <- getGEO(GEO, GSEMatrix =TRUE, getGPL = FALSE)




################################################################################################################################
#get meta
################################################################################################################################
#In this step, the grouping variable is selected from the meta data.
#This grouping variable is later used in the differentially expression analysis.



#Get grouping variables from the data set
groups <- get_grouping(gset)


#Select grouping variable
#auto-group automatically selects the relevant grouping variable.
#However, you can also select your grouping variable manually instead.
groupingvar <- groups[,auto_group(gset)]


#NOTE: You can also select a pairing variable in case of a dependent study design.


#Make meta table
#This table indicates to which group each sample belongs.
meta <- get_meta(gset = gset, grouping_column = groupingvar, pairing_column = NULL)



################################################################################################################################
#get data
################################################################################################################################



#GET CHIPTYPE

#The get_chiptype function automatically searches for this chiptype in the dataset.
#You can also select the chipytpe manually.
#IMPORTANT: the chiptype should be written without capitals, spaces, dots, comma's, etc.
#e.g. "hugene10st" or "huex10st"
chiptype = get_chiptype(gset)


#GET ORGANISM

#The get_organism function automatically searches for the organism in the dataset.
#You can also select the organism manually.
#IMPORTANT: the organism should be written in a two-letter code with lowercase letters only.
#e.g. "hs" or "mm"
organism = get_organism(gset)


#GET BRAINARRAY VERSION

#Currently, 25 is the latest version.
version = 25


#GET PROBESET ANNOTATION

#  "ense" for exon-specific probsets
#  "enst" for transcript-specific probesets
annotation = "ense"


#MAKE AFFYBATCH

#The readcels function requires the CDF to be in the working directory 
#unless the CDF has been installed as a package previously

data1 <- readcels(gset = gset, 
                  chiptype = chiptype, 
                  organism = organism, 
                  version = version,
                  annotation = annotation,
                  robust = FALSE,
                  outliers = "GSM4764672")



################################################################################################################################
#Normalization
################################################################################################################################



#Perform RMA normalization
data.rma <- affy::rma(data1, normalize = TRUE, background = TRUE)


#Retrieve probe set expression
data.expr <- exprs(data.rma)


#Remove the "nonsense" probeset. 
#This probeset is only included for normalization purposes.
data.expr <- data.expr[rownames(data.expr) != "nonsense",]



################################################################################################################################
#Quality control
################################################################################################################################



#BOXPLOTS OF LOG2 INTENSITIES

##Raw data
par(mar=c(8,4,2,1))
boxplot(data1,which='pm', col = "red", las = 2, main = "Boxplot of raw data", ylab = "log2 intensity")


##Normalized data
par(mar=c(8,4,2,1))
boxplot(data.expr, col = "blue", las = 2, main = "Boxplot of normalized data", ylab = "log2 intensity")




#DENSITY PLOTS OF LOG2 INTENSITIES

##Raw data
par(mar=c(5,4,2,1))
hist(data1,lwd=2,which='pm',ylab='Density',xlab='log2 intensity',main='Density plot of raw data')


##Normalized data
par(mar=c(5,4,2,1))
plotDensity(data.expr,lwd=2,ylab='Density',xlab='log2 intensity',main='Density plot of normalized data')



#PCA PLOT

data.PC <- prcomp(t(data.expr),scale.=TRUE)

pca.plot(data.PC = data.PC,
         meta = meta,
         hpc = 1,            #horizontal axis (1 = PC1 on horizontal axis, 2 = PC2 on horizontal axis, etc.)
         vpc = 2)            #vertical axis (2 = PC2 on vertical axis, 3 = PC3 on vertical axis, etc.)



################################################################################################################################
#Differential expression analysis
################################################################################################################################



#GET COMPARISONS

#Get comparisons for the differential expression analysis
#The get_contrast function retrieves all possible comparisons that are possible according to the meta table.
comparisons <- get_contrasts(meta)


#GET TOP TABLE

#Make top table using the meta data
#This function performs DE analysis for all possible comparisons and returns a list object.
top.table <- diff_expr(data.expr = data.expr, 
                       meta = meta, 
                       comparisons = comparisons)



################################################################################################################################
#SELECT AND ANNOTATE THE EXON-SPECIFIC PROBESETS
#IMPORTANT: you can only use run these codes if you used the "ense" annotation in the readcels function


#Get annotation file from the working directory
exonannotation <- read.table(file = paste(chiptype, organism, version, "ExonAnnotation.txt", sep = "_"), header = TRUE)


#Select exons of interest AND annotate the exons with their corresponding transcript(s) and gene
select.top.table <- exon_selection(top.table = top.table, 
                                   annotated = exonannotation,
                                   genes = c("ENSG00000065989", 
                                             "ENSG00000184588", 
                                             "ENSG00000105650",
                                             "ENSG00000113448"),
                                   transcripts = NULL,
                                   genome_wide = TRUE,
                                   unique_exons = FALSE)
  


################################################################################################################################
#SELECT AND ANNOTATE THE TRANSCRIPT-SPECIFIC PROBESETS
#IMPORTANT: you can only use run these codes if you used the "enst" annotation in the readcels function


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



#BOXPLOTS

#IMPORTANT: this function is not suitable for genome-wide analyses
#Thus, only run this function if you selected a hand-full of transcripts/exons in the select.top.table.
makeBoxplots1(contrast = comparisons[7],
              meta = meta,
              data.expr = data.expr,
              select.top.table = select.top.table,
              annotated = exonannotation)                    #In case of exon-specific probesets:
                                                   #annotated = exonannotation
           


#PROBE MAPPING PLOT

#These plots show the probes included in each transcript/exon-specific probesets
#These plots can not be used in genome-wide analyses

#Transcript-specific probe sets
probemapping.enst(chiptype = chiptype, 
                  organism = organism,
                  version = version, 
                  gene = "ENSG00000113448")



#Exon-specific probe sets
probemapping.ense(chiptype = chiptype, 
                  organism = organism,
                  version = version, 
                  gene = "ENSG00000113448",
                  annotated = exonannotation)



#VOLCANO PLOT

plotvolcano(select.top.table = select.top.table, 
            contrast = comparisons[1], 
            p = "raw",                           #choose from "raw" or "FDR"
            p.threshold = 0.05,
            logFC.threshold = 1)                            
  

