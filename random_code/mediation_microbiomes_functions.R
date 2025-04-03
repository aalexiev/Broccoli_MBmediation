## Mediation functions
# this script is for running mediation analysis using one ASV or taxon 
# in microbiome table with two variables
# in our case, these were dietary components from diet journals leading up to
# the experiment, and different types of SFN during the hours after broccoli
# consumption

# minimum input will be an ASV table (asv_tab), table of variable X (x_tab), 
# and table of variable Y (y_tab), with a column named "subject_id" to denote 
# the sample IDs that should be used to combine data in the script

# extra features:
# - taxa table (tax_tab) in case there needs to be changes in phylogenetic
# level of analysis (phylo_level)

MB_asvmed <- function(asv_tab, x_tab, y_tab, tax_tab, phylo_level) {
  
}