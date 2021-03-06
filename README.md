# Rare-OTUs-ACE
Python script to estimate ACE (Abundance-based Coverage Estimator) values per sample from OTU tables specifying n. of samples for OTUs to be considered rare

### Definitions
- ACE: abundance-based coverage estimator. One of the alpha diversity indexes used to measure diversity in ecological communities. Specifically, the ACE measures species richness.
- OTU: operational taxonomic unit. Distinct taxonomic elements (at the desired phylogenetic level) that are present in a sampled community.

### Description
This script estimates per-sample ACE values using the *diversity.alpha.ace* function from the **scikit-bio** Python package. Compared to what is implemented in the QIIME pipeline (Caporaso et al. 2010), users can here specify the threshold (n. of samples) for OTUs to be considered rare.

The script accepts as parameters:

1. a file with OTU counts per sample: *m* OTUs (rows) x *n* samples (columns), plus columns for OTU IDs and taxonomy (as generated by the closed OTU picking step in QIIME, with taxonomy specified
2. the threshold (n. of sample) under which an OTU is considered rare (relevant for the calculation of the ACE indicator)

Launch as:

   `python ace.py --file otu_table.csv -rt number_of_samples`
   
As output, the scripts produces a **csv file** (*ace.csv*) with two columns: sample ID and estimated ACE value.




