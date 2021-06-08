![image](https://user-images.githubusercontent.com/79576459/120113246-e3ae8100-c179-11eb-9ca3-7271d1fbe638.png)

# EXIS 1.0.0
EXIS is an interactive R shiny application for isoform-specific expression analysis using Affymetrix ST arrays.
The isoform-specific expresssion analysis can also be performed outside R shiny using the EXIS workflow.
The current version only supports the analysis of the Human Exon 1.0 ST array and the Human Gene 1.0. ST array.

## Instructions
1) Before using the app or workflow, source run the global R script (EXIS_Global.R).

IMPORTANT: The current analysis depends on the modified affy package from Brainarray.
If you have installed the original affy package previously, please remove this package from your package library before source running the global R script.

2) Make sure that the exon/transcript-specific CDF is in your working directory or has been installed as a CDF package previously.
Furthermore, when you would like to use the exon-specific probe set definition, also make sure to place the exon annotation file in your working directory.

3) If you would like to use the EXIS app, run the app R script (EXIS_App.R). If you open this script in RStudio, you can simply press "Run App" in the top-right corner. 
If you prefer using the functions outside the R shiny environment, you can use the workflow R script (EXIS_Workflow.R) to guide you through the analysis.
