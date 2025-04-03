## what is the sequence for ASV_516?

# this is the sequence table with ASV sequences instead of ASV number
seq_tab <- readRDS("/Users/alexieva/Documents/Projects/Analysis/broccoli project/00_data/Microbiome/seqtab_nochim.rds")

seq_tab <- as.data.frame(seq_tab) %>%
  rownames_to_column("subject_id") %>%
  dplyr::filter(grepl("_0h_", subject_id, ignore.case = TRUE))
seq_tab$subject_id <- strtrim(seq_tab$subject_id, 6)
subs <- med_input$subject_id
seq_tab <- seq_tab %>%
  dplyr::filter(subject_id %in% subs) %>%
  arrange(factor(subject_id, levels = med_input$subject_id))

# this is the asv table with ASV_###, filter for 516 only
asv_tab <- read.csv("/Users/alexieva/Documents/Projects/Analysis/broccoli project/00_data/Microbiome_Flat/ASVTable_BSS.csv") %>%
  dplyr::filter(grepl("_0h_", sample, ignore.case = TRUE)) %>%
  rename(subject_id = "sample") %>%
  dplyr::select(c(subject_id, ASV516))
asv_tab$subject_id <- strtrim(asv_tab$subject_id, 6)
asv_tab <- asv_tab %>%
  dplyr::filter(subject_id %in% subs) %>%
  arrange(factor(subject_id, levels = med_input$subject_id))

# are they in the same subject ID order? Then remove that column for next steps to work
seq_tab$subject_id == asv_tab$subject_id # same order
seq_tab <- column_to_rownames(seq_tab, "subject_id")
asv_tab <- column_to_rownames(asv_tab, "subject_id")

# match which sequence is the same as the counts in ASV_516 from the ASV table
out <- plyr::match_df(data.frame(t(seq_tab)), data.frame(t(asv_tab)), on = NULL)
rownames(out)
# TACGTAGGATGCAAGCGTTATCCGGAATGACTGGGCGTAAAGGGTGCGTAGGTGGTTTGTCAAGTTGGCAGCGTAATTCCGTGGCTTAACCGCGGAACTACTGCCAAAACTGATAGGCTTGAGTGCGGCAGGGGTATGTGGAATTCCTAGTGTAGCGGTGGAATGCGTAGATATTAGGAGGAACACCGGTGGCGAAAGCGACATACTGGGCCGTAACTGACACTGAAGCACGAAAGCGTGGGGAGCAAACAGG

## closest match in BLAST with 16S is Christensenella intestinihominis, then Desulfolucanica intricata
# Christensenella is in order Clostridiales, which matches taxa table



