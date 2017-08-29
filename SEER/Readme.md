This folder contains code for processing the SEER data.

To run the code, you will need CancersSel.RData file, which contains indivudual SEER records. 
Our Data-Use Agreement for the SEER 1973-2013 Research Data Files 
constrains that we can only release this data to a person with the written approval of the SEER Program. 
In particular, all members of a research team who will have access to the data must sign this data-use
agreement.

The code generates a subfolder Cancers with several summary plots, which are not directly relevant to the project.

The main output file is called Cancers.RData and contains summary statisics for clinical variables in 14 cancers 
and records their dependence on primary tumor diameter. This file is the input for the part of the project which contains
our metastasis formation project.
