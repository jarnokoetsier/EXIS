

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

#INSTRUCTIONS APP SCRIPT

###############################################################################################################################



#1. Source run the global R script.

#2. If you would like to generate an exon-specific CDF, 
#   you need to place the exon annotation file in your working directory.

#4. Use the MakeCDFfile function as decribed below to generate your exon- and/or transcript-specific CDF(s).



#############################################################################################################

#Make your own exon- and transcript-specific CDF

#############################################################################################################


MakeCDFfile(chiptype = "hugene10st",
            organism = "hs",
            version = 25,
            annotation = "ense",
            robust = TRUE,
            rows = 1050)



#CHIPTYPE: the chip type for which you would like to generate the CDF.
#Make sure to write the chip type without capital letters, dots, spaces, etc.
#e.g. "hugene10st", "huex10st", etc.


#ORGANISM: the organism for which the chip has been designed.
#Write the organism in a two-letter code with lowercase letters only.
#e.g. "hs", "mm", etc. 


#VERSION: Brainarray version used for the generation of the exon- or transcript-specific CDF

#ANNOTATION: "ense" if you want an exon-specific CDF and "enst" for a transcript-specifc CDF 

#ROBUST: "Set TRUE if you only want to include probesets with at least 3 probes. 
#Only applicable for transcript-specific probesets

#ROWS: number of rows of the chip.
#1050 for "hugene10st" and 2560 for "huex10st"



