![image](https://user-images.githubusercontent.com/79576459/120113246-e3ae8100-c179-11eb-9ca3-7271d1fbe638.png)

# EXIS 1.0.0
EXIS is an interactive R shiny application for isoform-specific expression analysis using Affymetrix ST arrays.
The isoform-specific expresssion analysis can also be performed outside R shiny using the EXIS workflow.
The current version only supports the analysis of the Human Exon 1.0 ST array and the Human Gene 1.0. ST array.

## Instructions
1) Before using the app or workflow, source run the global R script (EXIS_Global.R).

IMPORTANT: The current depends on the modified affy package from Brainarray.
If you have installed the original affy package previously, please remove this package from your package library before source running the global R script.

