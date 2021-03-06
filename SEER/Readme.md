This folder contains code for processing the SEER data (Surveillance epidemiology and end results program, 1973-2013 Data, Nov 2015 Submission).

To run the code, you will need CancersSel.RData file, which contains indivudual SEER records. 
Our Data-Use Agreement for the SEER 1973-2013 Research Data Files 
constrains that we can only release this data to a person with the written approval of the SEER Program. 
In particular, all members of a research team who will have access to the data must sign this data-use
agreement.
To see exactly which SEER records were preselected for our analysis, please refer to file STable_1.xls, which contains Histology, Histology Description, Histology/Behavior, Histology/Behavior Description from the SEER database (first four columns), for the cancer types we chose to analyze (last column).

The code generates a subfolder Cancers with several summary plots, which are not directly relevant to the project.

The main output file is called Cancers.RData and contains summary statisics for clinical variables in 14 cancers 
and records their dependence on primary tumor diameter. This file is the input for the part of the project which contains
our metastasis formation model.

This is a piece of R code, which depends on the following libraries:
ggplot2, calibrate, gridExtra, survival.

Running the code:

source("Metastasis_Main.R") 

in R.
