# survival_analysis_indian_adult_sepsis_cohort
This repo contains code and data to reproduce the results of the manuscript "Immunosuppression, rather than inflammation, is a salient feature of sepsis in an Indian cohort".

# Immunosuppression, rather than inflammation, is a salient feature of sepsis in an Indian cohort
Samanwoy Mukhopadhyay, Pravat K Thatoi, Bidyut K Das, Saroj K Mohapatra1*

## Availability of data
The dataset(s) and the R code supporting the conclusions of this article are available under the project ssnibmg in figshare https://figshare.com/projects/ssnibmgsurv/67721.

### Step 1: Installation of the Data Package ssnibmg
1. Dowload the files ssgeosurv_1.0.tar.gz, ssnibmgsurv_1.0.tar.gz from
https://figshare.com/ (search for ssnibmgsurv)
2. Change the directory to where you saved the file. Start R.
3. At the R prompt, issue the following command:
> install.packages(pkgs=``ssnibmgsurv_1.0.tar.gz'', repos=NULL)
4. Now the data package ssnibmg is installed on your computer.
5. Check with the following command:

> library(“ssnibmgsurv”)

> library(“ssgeosurv”)

### Step 2: Running the analysis code
6. Dowload the file ssnibmgsurvdoc.zip from
https://figshare.com/ (search for ssnibmgsurv)
7. Save the file at a suitable location and extract the contents. Change into the newly created directory ssnibmgsurvdoc.
8. The code for data analysis are listed in the vignette ssnibmgsurvdoc.pdf.
9. The Rcode and Metadata are subdirectories under ssnibmgsurvdoc (your current working directory in R).
10. The vignette is a self-computable document. Run the R code in steps to reproduce the results described in the vignette.

