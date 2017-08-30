This folder contains code for modeling metastasis formation and for reporting the results.

In order to fit the models to the data, run 

source("BCD_FitCancers.R") 

in R in a folder that contains the Cancers.RData file. 
Cancers.RData contains all data extracted from the SEER database, using the code in the SEER folder. We provide a copy of this file here.

Running this code produces a folder called Fit_BCDf_largeInit, which is also copied here for your convenience. This folder contains fitted parameters.

In order to reproduce all result figures & tables, please run 

source("BCD_SummarizeCancers.R")

in a folder that contains the Fit_BCDf_largeInit folder and HaenoData.txt file (also provided).

