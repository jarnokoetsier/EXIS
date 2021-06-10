

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

#INSTRUCTIONS MAKE EXON ANNOTATION FILE SCRIPT

###############################################################################################################################



#1. Fill-in the chiptype, organism, and Brainarray version

#2. Run the below-written code to generate the exon-annotation file

#IMPORTANT: running the code can take a (really!) long time.
#So, only run this code if the exon annotation file of your chip type is not available.
#Furthermore, consider running the code step by step and saving intermediate variables.




###############################################################################################################################

#FILL-IN

###############################################################################################################################


chiptype = "hugene10st"
organism = "hs"
version = 25



###############################################################################################################################

#RUN THE BELOW-WRITTEN CODE

###############################################################################################################################


#Download ENSE probe mapping file from Brainarray
mapping_file <- paste0(chiptype, "_", toupper(substring(organism, 1,1)), substring(organism, 2), "_ENSE_mapping.txt")


if (!file.exists(mapping_file)){
  file_download <- paste0("http://mbni.org/customcdf/", version, ".0.0/ense.download/", chiptype, "_",
                          toupper(substring(organism, 1,1)), substring(organism, 2), "_ENSE_", version, ".0.0.zip")
  
  download.file(file_download, "probe_information_file.zip")
  unzip("probe_information_file.zip", mapping_file)
  unlink("probe_information_file.zip")
  
}

probe_info <- read.table(mapping_file, sep = "\t", header = TRUE)


#Combine X and Y coordinate of probe to create unique probe ID's
probe_info$probe <- paste(probe_info$Probe.X, probe_info$Probe.Y, sep = "_")


#Remove "_at" from probe set name
probe_info$Probe.Set.Name <- str_replace_all(probe_info$Probe.Set.Name, "_at", "")


#Get probes that are present in multiple exons
probelist <- as.data.frame(table(unlist(probe_info$probe)))
probelist <- filter(y, Freq > 1)

redundantprobes <- inner_join(probe_info, probelist, by = c("probe" = "Var1"))



####################################################################
#Get exons that contain unique probes only
####################################################################

#Remove exons that contain probes that are present in multiple exons
removed_exons <- as.data.frame(unique(redundantprobes$Probe.Set.Name))
colnames(removed_exons) <- "removed_exons"
unique_exons <- anti_join(probe_info, removed_exons, by = c("Probe.Set.Name" = "removed_exons")) 
unique_exons <- unique(unique_exons$Probe.Set.Name)
Annotation <- cbind(paste0(unique_exons, "_at"), unique_exons, unique_exons)
colnames(Annotation) <- c("Probe.set", "Exon.group", "Exon")


####################################################################
#Group exons together of they contain the same probes
####################################################################

#From all the probes with a frequency > 2, arrange the frequency.
redundantprobes <- arrange(redundantprobes, desc(Freq))

#Get the frequency of the amount of probes per exon
exonfreq <- as.data.frame(table(unlist(probe_info$Probe.Set.Name)))
colnames(exonfreq) <- c("Var1", "exonfreq")

redundantprobes <- inner_join(redundantprobes, exonfreq, by = c("Probe.Set.Name" = "Var1"))


#Get all exons that contain the same probe.
#Afterwards, from these exons select the exon with the lowest number of probes
#Check whether all probes of this exon are present in the same exons
#If so, take this probe set and link all exons to this single probe set
Annotation2 <- NULL
z1 <- redundantprobes

repeat{
  
  int1 <- filter(z1, probe == z1$probe[1])
  int1 <- arrange(int1, exonfreq)
  
  test <- filter(z1, Probe.Set.Name == int1$Probe.Set.Name[1])
  
  if ((length(unique(test$Freq)) == 1) & (nrow(int1) > 1)){
    
    Exon.group <- rep(paste0(int1$Probe.Set.Name[1], "g"), nrow(int1))
    Probe.set <- rep(paste0(int1$Probe.Set.Name[1], "_at"), nrow(int1))
    Exon <- int1$Probe.Set.Name
    
    Annotation1 <- cbind(Probe.set, Exon.group, Exon)
    Annotation2 <- as.data.frame(rbind(Annotation2, Annotation1))
  }
  
  z1 <- anti_join(z1, as.data.frame(int1$Probe.Set.Name), by = c("Probe.Set.Name" = "int1$Probe.Set.Name"))
  
  if (nrow(z1) == 0){
    break
  }
  
}

all <- as.data.frame(rbind(Annotation, Annotation2))


#Get transcripts and gene for each exon
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)

exontranscriptgene <- getBM(
  mart=mart,
  attributes=c(
    "ensembl_exon_id",
    "ensembl_transcript_id",
    "ensembl_gene_id"
  ),
  filter = "ensembl_exon_id",
  values = all$Exon, uniqueRows=TRUE)


all1 <- left_join(all, exontranscriptgene, by = c("Exon" = "ensembl_exon_id"))
colnames(all1) <- c("Probe.set", "Exon.group", "Exon", "Transcript", "Gene")




############################################################################
#Remove probesets that contain probes that are present in the absent exons
############################################################################

annotated <- all1[complete.cases(all1),]

probe_info <- read.table(paste0(chiptype, "_", toupper(substring(organism, 1,1)), substring(organism, 2), "_ENSE_mapping.txt"), sep = "\t", header = TRUE)
probe_info$probe <- paste(probe_info$Probe.X, probe_info$Probe.Y, sep = "_")
probe_info1 <- probe_info[,c("Probe.Set.Name", "Chr.From", "probe")]

annotated1 <- inner_join(annotated, probe_info1, by = c("Probe.set" = "Probe.Set.Name"))


#Get chromosome location for each exon (also absent exons)
exonposition <- getBM(
  mart=mart,
  attributes=c(
    "ensembl_gene_id",
    "ensembl_exon_id",
    "exon_chrom_start",
    "exon_chrom_end"),
  filter = "ensembl_gene_id",
  values = unique(annotated$Gene), uniqueRows=TRUE)


#Get the absent exons in Brainarray's probe mapping
exon_without_probeset <- anti_join(exonposition, annotated1, by = c("ensembl_exon_id" = "Exon"))


#Test whether probes of each gene link to absent exons of that gene
annotated_filtered <- annotated

for (i in unique(annotated$Gene)[42877:46250]){
  gene_probes <- annotated1[annotated1$Gene == i,][,5:7] #Gene, Chr.From, probe
  
  test1 <- inner_join(exon_without_probeset, gene_probes, by = c("ensembl_gene_id" = "Gene"))
  
  #Test whether probes of a gene are present in exons without probe set
  test1$start <- test1$Chr.From - test1$exon_chrom_start
  test1$end <- test1$exon_chrom_end - test1$Chr.From
  
  final_test <- test1 %>%
    filter(start > 0) %>%
    filter(end > 0)
  
  final_test <- final_test[!duplicated(final_test),]
  
  #Remove exon probe sets which contain probes that are also present in exons without probe set
  exon_removed <- inner_join(annotated1[,c("Exon.group", "probe")], as.data.frame(final_test$probe), by = c("probe" = "final_test$probe"))
  exon_removed <- unique(exon_removed$Exon.group)
  
  annotated_filtered <- anti_join(annotated_filtered, as.data.frame(exon_removed), by = c("Exon.group" = "exon_removed"))
  
}



####################################################################
#Remove exons that map the less transcripts than its probes
####################################################################

#Get transcript mapping from brainarray
mapping_file <- paste0(chiptype, "_", toupper(substring(organism, 1,1)), substring(organism, 2), "_ENST_mapping.txt")


if (!file.exists(mapping_file)){
  file_download <- paste0("http://mbni.org/customcdf/", version, ".0.0/enst.download/", chiptype, "_",
                          toupper(substring(organism, 1,1)), substring(organism, 2), "_ENSt_", version, ".0.0.zip")
  
  download.file(file_download, "probe_information_file.zip")
  unzip("probe_information_file.zip", mapping_file)
  unlink("probe_information_file.zip")
  
}

#####TRANSCRIPT mapping
#Use transcript mapping to get the probes linked to each transcript
tprobe_info <- read.table(mapping_file, sep = "\t", header = TRUE)
tprobe_info$probe <- paste(tprobe_info$Probe.X, tprobe_info$Probe.Y, sep = "_")
tprobe_info$Probe.Set.Name <- str_replace_all(tprobe_info$Probe.Set.Name, "\\.._at", "")
tprobe_info$Probe.Set.Name <- str_replace_all(tprobe_info$Probe.Set.Name, "\\..._at", "")

#get probe frequency
probe_freq_transcript <- as.data.frame(table(unlist(tprobe_info$probe)))


#####EXON mapping
#Use exon mapping to get the probes linked to each transcript
annotated <- annotated_filtered
probes <- probe_info[,c("Probe.Set.Name", "probe")]

annotated1 <- inner_join(annotated, probes, by = c("Probe.set" = "Probe.Set.Name"))
annotated1 <- annotated1[, c("Transcript", "probe")]
annotated1 <- annotated1[!duplicated(annotated1),]


#get probe frequency
probe_freq_exon <- as.data.frame(table(unlist(annotated1$probe)))




#Number transcripts to which a probe links, should be the same for exon and transcript mapping

test <- inner_join(probe_freq_transcript, probe_freq_exon, by = c("Var1" = "Var1"))

test$diff <- abs(test$Freq.x - test$Freq.y)

probes_removed <- filter(test, diff > 0)


#Remove exon probe sets that contain these probes

annotated2 <- inner_join(annotated, probes, by = c("Probe.set" = "Probe.Set.Name"))
removed_exons <- inner_join(annotated2, probes_removed, by = c("probe" = "Var1"))
removed_exons <- as.data.frame(removed_exons$Exon.group)
colnames(removed_exons) <- "Exon.group"

filtered_annotated <- anti_join(annotated, removed_exons, by = c("Exon.group" = "Exon.group"))



###########################################################################################
#Remove exon probe sets that link to transcripts that are not present in transcript data
##########################################################################################


transcripts_exonsets <- as.data.frame(filtered_annotated$Transcript)
transcripts_transcriptset <- as.data.frame(tprobe_info$Probe.Set.Name)

transcripts_absent <- anti_join(transcripts_exonsets, 
                                transcripts_transcriptset, 
                                by = c("filtered_annotated$Transcript" = "tprobe_info$Probe.Set.Name"))

transcripts_absent <- unique(transcripts_absent)
removed_exons1 <- inner_join(filtered_annotated, transcripts_absent, by = c("Transcript" = "filtered_annotated$Transcript"))
removed_exons1 <- as.data.frame(unique(removed_exons1$Exon.group))
colnames(removed_exons1) <- "Exon.group"

filtered_annotated1 <- anti_join(filtered_annotated, removed_exons1, by = c("Exon.group" = "Exon.group"))


write.table(filtered_annotated1, file = paste(chiptype, organism, version, "ExonAnnotation.txt", sep = "_"), row.names = FALSE)




