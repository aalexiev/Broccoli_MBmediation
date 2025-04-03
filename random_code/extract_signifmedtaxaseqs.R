## read inputs
set.seed(3)
setwd("~/Documents/Projects/Analysis/broccoli project/00_data/Microbiome/")
asvtab <- as.data.frame(readRDS("seqtab_nochim.rds"))
library("dplyr")
# read in list of significant taxa
tax_list <- c("Collinsella", "Bifidobacterium",
              "Subdoligranulum", "Faecalibacterium",
              "[Ruminococcus] gauvreauii group",
              "[Ruminococcus] gnavus group",
              "Dialister")
tax <- readRDS("tax.rds") %>%
  as.data.frame() %>%
  tibble::rownames_to_column("seq") %>%
  dplyr::filter(Genus %in% tax_list)

asv516 <- c("TACGTAGGGGGCGAGCGTTGTCCGGAATTACTGGGCGTAAAGGGAGCGTAGGCGGTTAATTAAGTTAGATGTGAAATACCCGGGCTTAACTTGGGGGGTGCATCTAATACTGGTAAACTAGAGTACAGGAGAGGAAAGCGGAATTCCTAGTGTAGCGGTGAAATGCATAGATATTAGGAGGAACATCGGTGGCGAAGGCGGCTTTCTGGACTGACACTGACGCTGAGGCTCGAAAGCGTGGGGAGCAAACAGG",
"Bacteria", "Firmicutes", "Clostridia", NA, NA, NA, NA)

tax_2 <- rbind(tax, asv516)

seqs <- tax_2$seq
# seqs2 <- apply(array(seqs), 1, function (x) paste(">", x))

# write.table(tax_2,
#             file = "~/Documents/Projects/Analysis/broccoli project/01_explAnalysis/figstabs/signifmedtaxtable.txt",
#             sep = "\t", quote = F, row.names = F)

# write.table(seqs,
#             file = "~/Documents/Projects/Analysis/broccoli project/01_explAnalysis/figstabs/signifmedtaxSEQS.txt",
#             sep = "\t", quote = F, col.names = F, row.names = F)

asvtab_filt <- asvtab %>%
  dplyr::select(tax$seq)



