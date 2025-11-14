## libraries and directories

    library(tibble)
    library(dplyr)
    library(vegan)
    library(ggplot2)
    set.seed(3)

    setwd("/Users/alexieva/Documents/Projects/Analysis/broccoli project/02_pubAnalysis/Enterotyping/")

    load("/Users/alexieva/Documents/Projects/Analysis/broccoli project/02_pubAnalysis/Enterotyping//broccolipreMBpermanova.RData")

## Input data

    ## microbiome data
    # taxa table
    MB_taxtab <- read.csv("/Users/alexieva/Documents/Projects/Analysis/broccoli project/00_data/Microbiome_Flat/TaxTable_BSS.csv",
                          header = TRUE) 

    # metadata for microbiome stuff
    metadata <- read.csv("/Users/alexieva/Documents/Projects/Analysis/broccoli project/00_data/Microbiome_Flat/SampleData_BSS.csv",
                          header = TRUE) %>%
      dplyr::filter(time == 0) %>% # for now just Time 0 data
      mutate(broc_consum = case_when(veg == "broc" ~ "1",
                                     veg == "alf" ~ "0"))

    # ASV table (ASV reads per sample)
    MB_asvtab <- read.csv("/Users/alexieva/Documents/Projects/Analysis/broccoli project/00_data/Microbiome_Flat/ASVTable_BSS.csv",
                          header = TRUE) %>%
      dplyr::filter(sample %in% metadata$sample) # just time 0 data

    # clean up subject ID names so they will match the metadata
    MB_asvtab$sample <- sub("_[^_]+_[^_]+$", "", MB_asvtab$sample)

    # rename the ASVs by their taxonomy (genus level)
    MB_asvtab <-as.data.frame(t(column_to_rownames(MB_asvtab, "sample")))
    MB_asvtab$Genus[which(rownames(MB_asvtab) == MB_taxtab$ASV)] <- MB_taxtab$Genus
    MB_asvtab <- MB_asvtab %>%
      group_by(Genus) %>%
      summarise(across(where(is.numeric), sum)) %>%
      column_to_rownames("Genus") %>%
      t() %>%
      as.data.frame()

    ## SFN data
    # cumulative across all timepoints
    # include all SFN types
    ## SFN data

    Urine_SFN <- read.csv("/Users/alexieva/Documents/Projects/Analysis/broccoli project/00_data/SFN_Data/Urine_R_Format.csv",
                          header = TRUE) %>%
        group_by(subject_id) %>%
      summarise(sum_SFN = sum(SFN),
                sum_SFN_Cys = sum(SFN_Cys),
                sum_SFN_NAC = sum(SFN_NAC),
                sum_SFN_CG = sum(SFN_CG),
                sum_SFN_GSH = sum(SFN_GSH),
                sum_SFN_NIT = sum(SFN_NIT),
                sum_SFN_Tot = sum(SFN_Tot))

    # make a mega metadata file
    mega_meta <- inner_join(metadata, Urine_SFN, by = "subject_id")
    mega_meta <- mega_meta[order(mega_meta$subject_id),] # this also takes out the alfalfa eaters

    ## diet data
    # make a list of the "significant" diet factors we are interested in 
    diet_compnts <- c("diet_Chol_mg", "diet_SatFat_g", "diet_Alc_g", "diet_TotSolFib_g", 
                      "diet_MonSac_g", "diet_Disacc_g", "diet_BetaCaro_mcg", "diet_Vit_B1_mg",
                      "diet_Vit_B6_mg", "diet_Fol.DFE_mcg_DFE", "diet_Mang_mg", "diet_Copp_mg", 
                      "diet_TotFib_g", "diet_OCarb_g", "diet_Iron_mg", "diet_TotInsolFib_g")
    diet <- read.csv("/Users/alexieva/Documents/Projects/Analysis/broccoli project/00_data/Diet/New_Diet_Small.csv", header = TRUE) %>%
      dplyr::select(c("subject_id", all_of(diet_compnts)))

    library(phyloseq)

    # make per sample taxa plots (stacked bar plots)
    # ASV table (ASV reads per sample)
    MB_asvtab_2 <- read.csv("/Users/alexieva/Documents/Projects/Analysis/broccoli project/00_data/Microbiome_Flat/ASVTable_BSS.csv",
                          header = TRUE) %>%
      dplyr::filter(sample %in% metadata$sample) # just time 0 data
    # clean up subject ID names so they will match the metadata
    MB_asvtab_2$sample <- sub("_[^_]+_[^_]+$", "", MB_asvtab_2$sample)
    MB_asvtab_2 <- column_to_rownames(MB_asvtab_2, "sample")
    # transpose
    trans_mbasvtab <- as.data.frame(t(MB_asvtab_2))
    # make col into rownames in taxa table
    MB_taxtab2 <- MB_taxtab %>%
      column_to_rownames("ASV") %>%
      rename_at('Kingdom', ~ 'Domain') %>%
      as.matrix()

    # make phyloseq obj and plot
    OTU <- otu_table(trans_mbasvtab, taxa_are_rows = TRUE)
    TAX <- tax_table(MB_taxtab2)
    physeq <- phyloseq(OTU, TAX)
    TopNOTUs <- names(which(MB_taxtab2[,6] == "Prevotella_9" | MB_taxtab2[,6] == "Bacteroides" | MB_taxtab2[,6] == "Ruminococcus"))
    ent10 <- prune_taxa(TopNOTUs, physeq)
    plot_bar(ent10, fill = "Genus")

<img src="broccolipreMBpermanova_files/figure-markdown_strict/determine and assign enterotype-1.png" width="98%" height="98%" />


    # there does seem to be basis for 2-3 enterotypes here

## permanova tests on whole data set (with alfalfa eaters)

-   wanted to test permanova models with alfalfa eaters and no SFN so
    there is variation in dietary components like fiber

<!-- -->

    # with alfalfa eaters, there are 2 groups here, add diet
    meta_groups_alfbrocc <- read.table("/Users/alexieva/Documents/Projects/Analysis/broccoli project/01_explAnalysis/brocc_mediation_Rfiles/inputs/meta_groups_alfbrocc.txt",
                                       header = T, sep = "\t") %>%
      inner_join(diet, by = "subject_id")

    # without alfalfa eaters (just brocc), there are 3 groups here, add diet
    meta_groups_brocc <- read.table("/Users/alexieva/Documents/Projects/Analysis/broccoli project/01_explAnalysis/brocc_mediation_Rfiles/inputs/meta_groups_brocc.txt",
                                       header = T, sep = "\t") %>%
      inner_join(diet, by = "subject_id")

-   Ordistep used to choose best model for dbRDA.
-   This best model was then used in the PERMANOVA to get significance.

<!-- -->

    ALLmod0 <- capscale(MB_asvtab ~ 1, meta_groups_alfbrocc, distance = "bray")  # Model with intercept only
    ALLmod1 <- capscale(MB_asvtab ~ . + .*., meta_groups_alfbrocc, distance = "bray")  # Model with all explanatory variables
    ## problem: model is overfitted

    ordiALL <- ordistep(ALLmod0, scope = formula(ALLmod1)) # this determines what the best model is to run RDA on
    ## 
    ## Start: MB_asvtab ~ 1
    ## 
    ##                        Df  AIC       F Pr(>F)   
    ## + groups_alf.cluster    1  152 13.4221  0.005 **
    ## + group                 1  152 13.4221  0.005 **
    ## + diet_TotInsolFib_g    1  161  3.7486  0.005 **
    ## + diet_TotSolFib_g      1  162  2.7201  0.010 **
    ## + weight_lb             1  162  2.2776  0.020 * 
    ## + age                   1  162  2.3086  0.025 * 
    ## + diet_OCarb_g          1  162  2.1973  0.030 * 
    ## + height_in             1  162  1.9050  0.040 * 
    ## + sex                   1  163  1.7938  0.060 . 
    ## + sample_issue          1  162  1.9281  0.085 . 
    ## + diet_condensed        5  166  1.2988  0.100 . 
    ## + ethnicity             2  163  1.4760  0.105   
    ## + race                  9  168  1.2617  0.135   
    ## + diet_MonSac_g         1  163  1.3167  0.150   
    ## + diet_SatFat_g         1  163  1.4536  0.155   
    ## + diet_Vit_B6_mg        1  163  1.4415  0.160   
    ## + diet_TotFib_g         1  163  1.3853  0.180   
    ## + relation              6  167  1.2102  0.195   
    ## + diet_Copp_mg          1  163  1.2044  0.260   
    ## + diet_Alc_g            1  163  1.1021  0.360   
    ## + diet_Mang_mg          1  163  0.9984  0.385   
    ## + diet_Chol_mg          1  163  0.9620  0.465   
    ## + diet_Fol.DFE_mcg_DFE  1  163  0.9590  0.475   
    ## + cohort                1  163  0.8509  0.640   
    ## + bmi                   1  164  0.7307  0.665   
    ## + diet_Iron_mg          1  164  0.6967  0.785   
    ## + diet                 51  176  0.9011  0.790   
    ## + diet_Vit_B1_mg        1  164  0.6963  0.800   
    ## + treatment             3  166  0.8001  0.820   
    ## + diet_Disacc_g         1  164  0.6092  0.830   
    ## + diet_BetaCaro_mcg     1  164  0.5849  0.905   
    ## + broc_consum           1  164  0.4963  0.955   
    ## + veg                   1  164  0.4963  0.965   
    ## + sample               69 -Inf                  
    ## + subject_id           69 -Inf                  
    ## + time                  0  162                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Step: MB_asvtab ~ groups_alf.cluster 
    ## 
    ##                      Df    AIC      F Pr(>F)   
    ## - groups_alf.cluster  1 162.35 13.422  0.005 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##                        Df  AIC      F Pr(>F)   
    ## + diet_TotInsolFib_g    1  150 3.4861  0.005 **
    ## + age                   1  151 2.8820  0.005 **
    ## + diet_OCarb_g          1  151 2.8329  0.010 **
    ## + diet_TotSolFib_g      1  151 2.8323  0.010 **
    ## + race                  9  156 1.4264  0.010 **
    ## + sample_issue          1  151 2.2735  0.030 * 
    ## + diet_TotFib_g         1  152 1.6714  0.040 * 
    ## + weight_lb             1  152 1.7022  0.045 * 
    ## + ethnicity             2  152 1.7002  0.050 * 
    ## + diet_SatFat_g         1  152 1.7622  0.065 . 
    ## + diet_condensed        5  155 1.2670  0.075 . 
    ## + diet_Vit_B6_mg        1  152 1.5163  0.090 . 
    ## + diet_MonSac_g         1  152 1.5966  0.105   
    ## + height_in             1  152 1.3861  0.125   
    ## + diet_Alc_g            1  152 1.3998  0.140   
    ## + relation              6  156 1.2344  0.195   
    ## + diet_Copp_mg          1  152 1.2918  0.200   
    ## + sex                   1  153 1.1307  0.310   
    ## + diet_Fol.DFE_mcg_DFE  1  153 1.0934  0.320   
    ## + diet_Chol_mg          1  153 1.1135  0.385   
    ## + diet                 51  155 1.0372  0.435   
    ## + cohort                1  153 0.9354  0.500   
    ## + diet_Vit_B1_mg        1  153 0.8406  0.615   
    ## + diet_Iron_mg          1  153 0.8483  0.640   
    ## + treatment             3  155 0.8620  0.690   
    ## + bmi                   1  153 0.7443  0.740   
    ## + broc_consum           1  153 0.6838  0.780   
    ## + diet_Disacc_g         1  153 0.6863  0.790   
    ## + veg                   1  153 0.6838  0.815   
    ## + diet_Mang_mg          1  153 0.6144  0.870   
    ## + diet_BetaCaro_mcg     1  153 0.5416  0.925   
    ## + sample               68 -Inf                 
    ## + subject_id           68 -Inf                 
    ## + time                  0  152                 
    ## + group                 0  152                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Step: MB_asvtab ~ groups_alf.cluster + diet_TotInsolFib_g 
    ## 
    ##                      Df    AIC       F Pr(>F)   
    ## - diet_TotInsolFib_g  1 151.74  3.4861  0.005 **
    ## - groups_alf.cluster  1 160.59 12.9893  0.005 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##                                         Df  AIC      F Pr(>F)   
    ## + diet_OCarb_g                           1  148 3.5972  0.005 **
    ## + age                                    1  150 2.3454  0.015 * 
    ## + sample_issue                           1  150 2.4931  0.030 * 
    ## + diet_SatFat_g                          1  150 1.9419  0.030 * 
    ## + weight_lb                              1  150 1.7642  0.050 * 
    ## + diet_Alc_g                             1  150 1.6714  0.050 * 
    ## + race                                   9  155 1.3173  0.050 * 
    ## + height_in                              1  151 1.5013  0.080 . 
    ## + diet_condensed                         5  153 1.2477  0.085 . 
    ## + diet_Vit_B6_mg                         1  151 1.4042  0.115   
    ## + groups_alf.cluster:diet_TotInsolFib_g  1  151 1.3125  0.190   
    ## + diet_Vit_B1_mg                         1  151 1.1494  0.245   
    ## + diet_Chol_mg                           1  151 1.1823  0.250   
    ## + sex                                    1  151 1.1688  0.285   
    ## + relation                               6  155 1.0586  0.315   
    ## + diet_TotSolFib_g                       1  151 1.0760  0.320   
    ## + ethnicity                              2  152 1.1215  0.320   
    ## + diet                                  51  148 1.0697  0.325   
    ## + diet_Mang_mg                           1  151 0.9866  0.400   
    ## + treatment                              3  153 0.9539  0.530   
    ## + cohort                                 1  151 0.9555  0.535   
    ## + diet_TotFib_g                          1  151 0.9218  0.555   
    ## + diet_Fol.DFE_mcg_DFE                   1  151 0.8256  0.570   
    ## + diet_Iron_mg                           1  151 0.8594  0.585   
    ## + diet_Copp_mg                           1  151 0.8536  0.640   
    ## + veg                                    1  151 0.7965  0.660   
    ## + bmi                                    1  151 0.8331  0.665   
    ## + broc_consum                            1  151 0.7965  0.680   
    ## + diet_Disacc_g                          1  151 0.7881  0.680   
    ## + diet_MonSac_g                          1  152 0.6380  0.855   
    ## + diet_BetaCaro_mcg                      1  152 0.5585  0.930   
    ## + sample                                67 -Inf                 
    ## + subject_id                            67 -Inf                 
    ## + time                                   0  150                 
    ## + group                                  0  150                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Step: MB_asvtab ~ groups_alf.cluster + diet_TotInsolFib_g + diet_OCarb_g 
    ## 
    ##                      Df    AIC       F Pr(>F)   
    ## - diet_OCarb_g        1 150.19  3.5972  0.005 **
    ## - diet_TotInsolFib_g  1 150.84  4.2481  0.005 **
    ## - groups_alf.cluster  1 159.66 13.6781  0.005 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##                                         Df  AIC      F Pr(>F)  
    ## + weight_lb                              1  149 1.8011  0.030 *
    ## + sample_issue                           1  148 2.5256  0.040 *
    ## + age                                    1  149 1.7287  0.045 *
    ## + diet_OCarb_g:diet_TotInsolFib_g        1  149 1.7239  0.065 .
    ## + height_in                              1  149 1.5210  0.090 .
    ## + diet_Vit_B6_mg                         1  149 1.4456  0.100 .
    ## + race                                   9  154 1.1828  0.125  
    ## + groups_alf.cluster:diet_TotInsolFib_g  1  149 1.3934  0.135  
    ## + sex                                    1  149 1.2732  0.135  
    ## + diet_Alc_g                             1  149 1.3604  0.145  
    ## + groups_alf.cluster:diet_OCarb_g        1  149 1.2802  0.220  
    ## + diet_TotSolFib_g                       1  149 1.1677  0.235  
    ## + diet_TotFib_g                          1  149 1.1224  0.255  
    ## + ethnicity                              2  150 1.1130  0.275  
    ## + cohort                                 1  149 1.1090  0.295  
    ## + diet_condensed                         5  152 1.1078  0.295  
    ## + diet                                  51  144 1.0535  0.305  
    ## + relation                               6  153 1.0962  0.320  
    ## + diet_Chol_mg                           1  149 0.9250  0.530  
    ## + diet_SatFat_g                          1  150 0.8727  0.575  
    ## + diet_Fol.DFE_mcg_DFE                   1  150 0.8673  0.585  
    ## + bmi                                    1  150 0.8660  0.600  
    ## + diet_Copp_mg                           1  150 0.8732  0.610  
    ## + diet_Mang_mg                           1  150 0.8726  0.645  
    ## + diet_Disacc_g                          1  150 0.8102  0.650  
    ## + diet_BetaCaro_mcg                      1  150 0.8060  0.670  
    ## + diet_MonSac_g                          1  150 0.6512  0.825  
    ## + treatment                              3  152 0.8168  0.825  
    ## + diet_Iron_mg                           1  150 0.6729  0.830  
    ## + diet_Vit_B1_mg                         1  150 0.6882  0.870  
    ## + veg                                    1  150 0.6308  0.885  
    ## + broc_consum                            1  150 0.6308  0.890  
    ## + sample                                66 -Inf                
    ## + subject_id                            66 -Inf                
    ## + time                                   0  148                
    ## + group                                  0  148                
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Step: MB_asvtab ~ groups_alf.cluster + diet_TotInsolFib_g + diet_OCarb_g +      weight_lb 
    ## 
    ##                      Df    AIC       F Pr(>F)   
    ## - weight_lb           1 148.48  1.8011  0.055 . 
    ## - diet_OCarb_g        1 150.34  3.6080  0.005 **
    ## - diet_TotInsolFib_g  1 151.04  4.2931  0.005 **
    ## - groups_alf.cluster  1 159.24 12.9066  0.005 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##                                         Df  AIC      F Pr(>F)  
    ## + weight_lb:diet_OCarb_g                 1  149 1.8407  0.025 *
    ## + sample_issue                           1  148 2.3160  0.050 *
    ## + age                                    1  149 1.7764  0.055 .
    ## + diet_OCarb_g:diet_TotInsolFib_g        1  149 1.7551  0.055 .
    ## + weight_lb:groups_alf.cluster           1  149 1.4751  0.105  
    ## + groups_alf.cluster:diet_TotInsolFib_g  1  149 1.4047  0.105  
    ## + race                                   9  154 1.1823  0.125  
    ## + relation                               6  153 1.1577  0.155  
    ## + diet_Alc_g                             1  149 1.3307  0.170  
    ## + diet_Vit_B6_mg                         1  149 1.3340  0.175  
    ## + groups_alf.cluster:diet_OCarb_g        1  149 1.3040  0.185  
    ## + sex                                    1  149 1.1541  0.185  
    ## + diet_condensed                         5  152 1.1389  0.195  
    ## + diet_TotFib_g                          1  149 1.1773  0.245  
    ## + cohort                                 1  149 1.1195  0.250  
    ## + diet_TotSolFib_g                       1  149 1.1555  0.285  
    ## + diet                                  51  142 1.0288  0.380  
    ## + ethnicity                              2  150 1.0604  0.380  
    ## + diet_Fol.DFE_mcg_DFE                   1  150 0.9074  0.520  
    ## + diet_SatFat_g                          1  150 0.8786  0.530  
    ## + diet_Mang_mg                           1  150 0.9285  0.550  
    ## + diet_Chol_mg                           1  150 0.8504  0.595  
    ## + diet_Copp_mg                           1  150 0.8510  0.625  
    ## + diet_Disacc_g                          1  150 0.8079  0.625  
    ## + diet_BetaCaro_mcg                      1  150 0.8129  0.660  
    ## + weight_lb:diet_TotInsolFib_g           1  150 0.7717  0.720  
    ## + diet_Iron_mg                           1  150 0.6821  0.845  
    ## + diet_Vit_B1_mg                         1  150 0.6759  0.855  
    ## + diet_MonSac_g                          1  150 0.6137  0.855  
    ## + veg                                    1  150 0.6179  0.900  
    ## + bmi                                    1  150 0.6262  0.910  
    ## + broc_consum                            1  150 0.6179  0.925  
    ## + treatment                              3  152 0.7295  0.950  
    ## + height_in                              1  150 0.5478  0.960  
    ## + sample                                65 -Inf                
    ## + subject_id                            65 -Inf                
    ## + time                                   0  149                
    ## + group                                  0  149                
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Step: MB_asvtab ~ groups_alf.cluster + diet_TotInsolFib_g + diet_OCarb_g +      weight_lb + diet_OCarb_g:weight_lb 
    ## 
    ##                          Df    AIC       F Pr(>F)   
    ## - diet_OCarb_g:weight_lb  1 148.56  1.8407  0.050 * 
    ## - diet_TotInsolFib_g      1 150.54  3.7303  0.005 **
    ## - groups_alf.cluster      1 159.69 13.1900  0.005 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##                                         Df  AIC      F Pr(>F)  
    ## + sample_issue                           1  148 2.3993  0.030 *
    ## + age                                    1  149 1.6330  0.070 .
    ## + diet_OCarb_g:diet_TotInsolFib_g        1  149 1.5597  0.080 .
    ## + weight_lb:groups_alf.cluster           1  149 1.5018  0.090 .
    ## + diet_Alc_g                             1  149 1.4128  0.095 .
    ## + race                                   9  154 1.2189  0.110  
    ## + groups_alf.cluster:diet_TotInsolFib_g  1  149 1.4549  0.120  
    ## + groups_alf.cluster:diet_OCarb_g        1  149 1.3836  0.125  
    ## + diet_condensed                         5  152 1.1268  0.240  
    ## + relation                               6  153 1.0838  0.250  
    ## + diet_TotSolFib_g                       1  149 1.1834  0.265  
    ## + sex                                    1  149 1.1724  0.265  
    ## + diet_Vit_B6_mg                         1  149 1.0563  0.300  
    ## + cohort                                 1  149 1.1375  0.315  
    ## + ethnicity                              2  150 1.0628  0.320  
    ## + diet_BetaCaro_mcg                      1  150 0.9501  0.475  
    ## + diet_Mang_mg                           1  150 0.9625  0.480  
    ## + diet_SatFat_g                          1  150 0.9371  0.490  
    ## + diet_TotFib_g                          1  149 0.9829  0.535  
    ## + diet                                  51  139 0.9928  0.540  
    ## + diet_Chol_mg                           1  150 0.8770  0.580  
    ## + diet_Disacc_g                          1  150 0.8419  0.645  
    ## + diet_Fol.DFE_mcg_DFE                   1  150 0.7579  0.675  
    ## + diet_Copp_mg                           1  150 0.7919  0.690  
    ## + weight_lb:diet_TotInsolFib_g           1  150 0.6763  0.785  
    ## + treatment                              3  152 0.8193  0.830  
    ## + diet_Vit_B1_mg                         1  150 0.6506  0.845  
    ## + bmi                                    1  150 0.6872  0.870  
    ## + diet_Iron_mg                           1  150 0.6490  0.875  
    ## + diet_MonSac_g                          1  150 0.6140  0.885  
    ## + veg                                    1  150 0.6204  0.905  
    ## + broc_consum                            1  150 0.6204  0.910  
    ## + height_in                              1  150 0.5596  0.960  
    ## + sample                                64 -Inf                
    ## + subject_id                            64 -Inf                
    ## + time                                   0  149                
    ## + group                                  0  149                
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Step: MB_asvtab ~ groups_alf.cluster + diet_TotInsolFib_g + diet_OCarb_g +      weight_lb + sample_issue + diet_OCarb_g:weight_lb 
    ## 
    ##                          Df    AIC       F Pr(>F)   
    ## - sample_issue            1 148.58  2.3993  0.035 * 
    ## - diet_OCarb_g:weight_lb  1 148.07  1.9306  0.030 * 
    ## - diet_TotInsolFib_g      1 150.23  3.9593  0.005 **
    ## - groups_alf.cluster      1 159.52 13.4673  0.005 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##                                         Df  AIC      F Pr(>F)  
    ## + diet_OCarb_g:diet_TotInsolFib_g        1  148 1.5590  0.060 .
    ## + weight_lb:groups_alf.cluster           1  148 1.4868  0.090 .
    ## + groups_alf.cluster:diet_TotInsolFib_g  1  148 1.4779  0.095 .
    ## + age                                    1  148 1.3626  0.100 .
    ## + diet_Alc_g                             1  148 1.4953  0.135  
    ## + race                                   9  154 1.1373  0.145  
    ## + groups_alf.cluster:diet_OCarb_g        1  148 1.3839  0.155  
    ## + diet_TotSolFib_g                       1  149 1.2271  0.225  
    ## + diet_condensed                         5  152 1.1124  0.255  
    ## + cohort                                 1  149 1.1309  0.265  
    ## + sex                                    1  149 1.1780  0.275  
    ## + relation                               6  152 1.1047  0.320  
    ## + diet_TotFib_g                          1  149 1.0036  0.425  
    ## + diet_SatFat_g                          1  149 0.9681  0.465  
    ## + diet_Vit_B6_mg                         1  149 0.9237  0.495  
    ## + diet_Copp_mg                           1  149 0.9025  0.530  
    ## + diet_Chol_mg                           1  149 0.8886  0.565  
    ## + diet_BetaCaro_mcg                      1  149 0.8525  0.570  
    ## + diet_Mang_mg                           1  149 0.8933  0.580  
    ## + diet                                  50  139 0.9660  0.615  
    ## + veg                                    1  149 0.7723  0.685  
    ## + broc_consum                            1  149 0.7723  0.710  
    ## + diet_Disacc_g                          1  149 0.7951  0.720  
    ## + weight_lb:diet_TotInsolFib_g           1  149 0.6766  0.775  
    ## + diet_Iron_mg                           1  149 0.6975  0.780  
    ## + treatment                              3  151 0.8691  0.785  
    ## + diet_Vit_B1_mg                         1  149 0.6740  0.830  
    ## + diet_Fol.DFE_mcg_DFE                   1  149 0.6604  0.840  
    ## + ethnicity                              2  150 0.7040  0.850  
    ## + diet_MonSac_g                          1  149 0.6267  0.875  
    ## + bmi                                    1  149 0.6937  0.885  
    ## + height_in                              1  149 0.5795  0.950  
    ## + sample                                63 -Inf                
    ## + subject_id                            63 -Inf                
    ## + time                                   0  148                
    ## + group                                  0  148                
    ## + weight_lb:sample_issue                 0  148                
    ## + sample_issue:groups_alf.cluster        0  148                
    ## + sample_issue:diet_OCarb_g              0  148                
    ## + sample_issue:diet_TotInsolFib_g        0  148                
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ordiALL
    ## Call: capscale(formula = MB_asvtab ~ groups_alf.cluster +
    ## diet_TotInsolFib_g + diet_OCarb_g + weight_lb + sample_issue +
    ## diet_OCarb_g:weight_lb, data = meta_groups_alfbrocc, distance = "bray")
    ## 
    ## -- Model Summary --
    ## 
    ##               Inertia Proportion Rank
    ## Total          9.0628                
    ## RealTotal     10.0254     1.0000     
    ## Constrained    3.1488     0.3141    6
    ## Unconstrained  6.8766     0.6859   42
    ## Imaginary     -0.9627                
    ## 
    ## Inertia is squared Bray distance
    ## 
    ## -- Note --
    ## 
    ## Species scores projected from 'MB_asvtab'
    ## 
    ## -- Eigenvalues --
    ## 
    ## Eigenvalues for constrained axes:
    ##   CAP1   CAP2   CAP3   CAP4   CAP5   CAP6 
    ## 1.7558 0.6052 0.3833 0.2091 0.1216 0.0738 
    ## 
    ## Eigenvalues for unconstrained axes:
    ##   MDS1   MDS2   MDS3   MDS4   MDS5   MDS6   MDS7   MDS8 
    ## 1.2571 0.7334 0.5614 0.4406 0.4113 0.3895 0.3021 0.2941 
    ## (Showing 8 of 42 unconstrained eigenvalues)
    # best model is MB_asvtab ~ groups_alf.cluster + diet_TotInsolFib_g + diet_OCarb_g + age
    # + sample_issue + weight_lb + diet_TotInsolFib_g:age
    # model inertia was 9.0628; Inertia is squared Bray distance

    # calculate RDA (distance based)
    rda <- capscale(formula = MB_asvtab ~ groups_alf.cluster + diet_TotInsolFib_g + diet_OCarb_g + age
    + sample_issue + weight_lb + diet_TotInsolFib_g:age, 
                       data = meta_groups_alfbrocc, distance = "bray")
    RsquareAdj(rda)
    ## $r.squared
    ## [1] 0.3316808
    ## 
    ## $adj.r.squared
    ## [1] 0.2562254
    # 0.256 of variation explained by this model
    summary(rda)
    ## 
    ## Call:
    ## capscale(formula = MB_asvtab ~ groups_alf.cluster + diet_TotInsolFib_g +      diet_OCarb_g + age + sample_issue + weight_lb + diet_TotInsolFib_g:age,      data = meta_groups_alfbrocc, distance = "bray") 
    ## 
    ## Partitioning of squared Bray distance:
    ##               Inertia Proportion
    ## Total          10.025     1.0000
    ## Constrained     3.325     0.3317
    ## Unconstrained   6.700     0.6683
    ## 
    ## Eigenvalues, and their contribution to the squared Bray distance 
    ## 
    ## Importance of components:
    ##                         CAP1    CAP2    CAP3    CAP4    CAP5    CAP6     CAP7
    ## Eigenvalue            1.7484 0.64697 0.33609 0.21510 0.19026 0.12519 0.063182
    ## Proportion Explained  0.1744 0.06453 0.03352 0.02146 0.01898 0.01249 0.006302
    ## Cumulative Proportion 0.1744 0.23893 0.27246 0.29391 0.31289 0.32538 0.331681
    ##                         MDS1    MDS2    MDS3    MDS4    MDS5   MDS6   MDS7
    ## Eigenvalue            1.2169 0.78674 0.54165 0.42199 0.40157 0.3589 0.2917
    ## Proportion Explained  0.1214 0.07847 0.05403 0.04209 0.04005 0.0358 0.0291
    ## Cumulative Proportion 0.4531 0.53153 0.58556 0.62765 0.66771 0.7035 0.7326
    ##                          MDS8    MDS9   MDS10   MDS11   MDS12   MDS13   MDS14
    ## Eigenvalue            0.26057 0.25100 0.24482 0.21432 0.17379 0.14046 0.13853
    ## Proportion Explained  0.02599 0.02504 0.02442 0.02138 0.01734 0.01401 0.01382
    ## Cumulative Proportion 0.75860 0.78363 0.80805 0.82943 0.84677 0.86078 0.87459
    ##                         MDS15   MDS16   MDS17    MDS18    MDS19    MDS20
    ## Eigenvalue            0.13068 0.11491 0.10492 0.093565 0.088774 0.083472
    ## Proportion Explained  0.01304 0.01146 0.01047 0.009333 0.008855 0.008326
    ## Cumulative Proportion 0.88763 0.89909 0.90956 0.918890 0.927744 0.936071
    ##                          MDS21    MDS22    MDS23    MDS24    MDS25    MDS26
    ## Eigenvalue            0.081648 0.066187 0.060222 0.057592 0.051457 0.046033
    ## Proportion Explained  0.008144 0.006602 0.006007 0.005745 0.005133 0.004592
    ## Cumulative Proportion 0.944215 0.950817 0.956824 0.962568 0.967701 0.972292
    ##                          MDS27    MDS28    MDS29    MDS30    MDS31    MDS32
    ## Eigenvalue            0.045408 0.041475 0.029398 0.027455 0.026013 0.018506
    ## Proportion Explained  0.004529 0.004137 0.002932 0.002739 0.002595 0.001846
    ## Cumulative Proportion 0.976822 0.980959 0.983891 0.986630 0.989224 0.991070
    ##                          MDS33    MDS34    MDS35    MDS36     MDS37    MDS38
    ## Eigenvalue            0.017128 0.016087 0.011193 0.010445 0.0091846 0.008502
    ## Proportion Explained  0.001708 0.001605 0.001116 0.001042 0.0009161 0.000848
    ## Cumulative Proportion 0.992779 0.994383 0.995500 0.996542 0.9974577 0.998306
    ##                           MDS39     MDS40     MDS41     MDS42
    ## Eigenvalue            0.0064944 0.0052861 0.0032951 0.0019107
    ## Proportion Explained  0.0006478 0.0005273 0.0003287 0.0001906
    ## Cumulative Proportion 0.9989535 0.9994807 0.9998094 1.0000000
    ## 
    ## Accumulated constrained eigenvalues
    ## Importance of components:
    ##                         CAP1   CAP2   CAP3    CAP4    CAP5    CAP6    CAP7
    ## Eigenvalue            1.7484 0.6470 0.3361 0.21510 0.19026 0.12519 0.06318
    ## Proportion Explained  0.5258 0.1946 0.1011 0.06469 0.05722 0.03765 0.01900
    ## Cumulative Proportion 0.5258 0.7204 0.8214 0.88613 0.94335 0.98100 1.00000

    plot(rda)

<img src="broccolipreMBpermanova_files/figure-markdown_strict/dbRDA run with best model from ordistep-1.png" width="98%" height="98%" />

    dm_walf <- vegdist(MB_asvtab, method = "bray")
    Perm <- adonis2(dm_walf ~ groups_alf.cluster + diet_TotInsolFib_g + diet_OCarb_g + age
    + sample_issue + weight_lb + diet_TotInsolFib_g:age, 
                                data = meta_groups_alfbrocc)
    Perm
    ## Permutation test for adonis under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = dm_walf ~ groups_alf.cluster + diet_TotInsolFib_g + diet_OCarb_g + age + sample_issue + weight_lb + diet_TotInsolFib_g:age, data = meta_groups_alfbrocc)
    ##          Df SumOfSqs      R2      F Pr(>F)    
    ## Model     7   3.2647 0.36023 4.9871  0.001 ***
    ## Residual 62   5.7981 0.63977                  
    ## Total    69   9.0628 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    # write.table(Perm,
    #             "Perm_w_alf.txt",
    #             sep = "\t")

    # create files for graphing
    smry_rda <- summary(rda)

    # this is the coordinates for the points and metadata
    meta_groups_filt <- meta_groups_alfbrocc %>%
      dplyr::select(c(subject_id, group, diet_TotInsolFib_g, diet_OCarb_g, veg, age, weight_lb, sample_issue))
    pw_PC1  <- data.frame(smry_rda$sites[,1:2]) %>%  # these are the x, y coordinates for the sample points
      rownames_to_column("subject_id") %>%
      inner_join(meta_groups_filt, by = "subject_id") %>%
      column_to_rownames("subject_id")
    mylims <- range(with(pw_PC1, c(CAP1, CAP2)))
    pw_PC1$group <- as.factor(pw_PC1$group)

    # the biplot scores are the correlations between your environmental variables and axes
    pw_PC2  <- data.frame(smry_rda$biplot[1:6,])

    # put plot together
    rda_plot <- ggplot(pw_PC1, aes(x = CAP1, y = CAP2)) + 
      geom_point(aes(color = group), size = 2) +
      theme_classic() +
      geom_hline(yintercept = 0, linetype = "dotted") +
      geom_vline(xintercept = 0, linetype = "dotted") +
      stat_ellipse(aes(group = veg, linetype = veg), show.legend = T) +
      geom_segment(data = pw_PC2, aes(x = 0, xend = CAP1, y = 0, yend = CAP2), 
                   color = "black", arrow = arrow(length = unit(0.01, "npc"))) +
      geom_text(data = pw_PC2, 
                aes(x = CAP1, y = CAP2), label = row.names(pw_PC2), 
                color = "black", size = 4, hjust = 0, vjust = 1) +
      labs(x = paste0("CAP1 (",round(100*smry_rda$cont$importance[2, "CAP1"], digits = 2),"%)"),
           y = paste0("CAP2 (",round(100*smry_rda$cont$importance[2, "CAP2"], digits = 2),"%)")) +  
      scale_color_brewer(palette = "Set1",
                        name = "pre-treatment microbiome group") +
      theme(text = element_text(size = 20))
    rda_plot

<img src="broccolipreMBpermanova_files/figure-markdown_strict/rda graph-1.png" width="98%" height="98%" />


    # ggsave("wALFpwrda_plot.tiff",
    #        width = 12,
    #        height = 7,
    #        dpi = 300)

## permanova tests on only broccoli consumers

-   removed alfalfa eaters and add sfn data

<!-- -->

    MB_asvtab_brocc <- MB_asvtab %>%
      rownames_to_column("subject_id") %>%
      dplyr::filter(subject_id %in% meta_groups_brocc$subject_id) %>%
      column_to_rownames("subject_id")

-   Ordistep used to choose best model for dbRDA.
-   This best model was then used in the PERMANOVA to get significance.

<!-- -->

    meta_groups_brocc <- meta_groups_brocc %>%
      dplyr::select(-c(sample, treatment, sample_issue, time, veg, broc_consum, diet))
    broccmod0 <- capscale(MB_asvtab_brocc ~ 1, meta_groups_brocc, distance = "bray")  # Model with intercept only
    broccmod1 <- capscale(MB_asvtab_brocc ~ . + .*., meta_groups_brocc, distance = "bray")  # Model with all explanatory variables
    ordibrocc <- ordistep(broccmod0, scope = formula(broccmod1)) # this determines what the best model is to run RDA on
    ## 
    ## Start: MB_asvtab_brocc ~ 1
    ## 
    ##                        Df  AIC      F Pr(>F)   
    ## + group                 2   50 8.8017  0.005 **
    ## + diet_TotInsolFib_g    1   61 2.5693  0.010 **
    ## + age                   1   61 2.2433  0.025 * 
    ## + diet_condensed        4   63 1.5741  0.025 * 
    ## + diet_OCarb_g          1   61 2.0744  0.030 * 
    ## + diet_SatFat_g         1   61 2.0721  0.045 * 
    ## + sex                   1   61 1.8656  0.045 * 
    ## + diet_TotSolFib_g      1   62 1.5053  0.090 . 
    ## + height_in             1   62 1.5384  0.110   
    ## + race                  6   65 1.2796  0.115   
    ## + weight_lb             1   62 1.4719  0.150   
    ## + diet_Fol.DFE_mcg_DFE  1   62 1.3890  0.150   
    ## + relation              3   63 1.2863  0.165   
    ## + diet_MonSac_g         1   62 1.2712  0.205   
    ## + diet_Chol_mg          1   62 1.2888  0.230   
    ## + diet_Mang_mg          1   62 1.0341  0.280   
    ## + diet_BetaCaro_mcg     1   62 1.1219  0.340   
    ## + ethnicity             2   63 1.0987  0.345   
    ## + diet_Copp_mg          1   62 1.1123  0.385   
    ## + diet_TotFib_g         1   62 0.9651  0.445   
    ## + diet_Iron_mg          1   62 0.8929  0.535   
    ## + diet_Alc_g            1   62 0.8545  0.545   
    ## + diet_Vit_B6_mg        1   62 0.6933  0.735   
    ## + diet_Vit_B1_mg        1   62 0.6390  0.760   
    ## + bmi                   1   62 0.7046  0.790   
    ## + diet_Disacc_g         1   63 0.5700  0.890   
    ## + cohort                1   63 0.5025  0.935   
    ## + subject_id           37 -Inf                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Step: MB_asvtab_brocc ~ group 
    ## 
    ##         Df    AIC      F Pr(>F)   
    ## - group  2 61.142 8.8017  0.005 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##                        Df  AIC      F Pr(>F)   
    ## + diet_TotSolFib_g      1   49 2.0965  0.010 **
    ## + diet_OCarb_g          1   50 1.8955  0.020 * 
    ## + diet_condensed        4   51 1.5002  0.020 * 
    ## + diet_SatFat_g         1   50 1.8808  0.030 * 
    ## + diet_TotInsolFib_g    1   50 1.9506  0.040 * 
    ## + weight_lb             1   50 1.7457  0.040 * 
    ## + diet_Chol_mg          1   50 1.6434  0.055 . 
    ## + race                  6   52 1.3838  0.055 . 
    ## + sex                   1   50 1.6218  0.060 . 
    ## + height_in             1   50 1.5661  0.075 . 
    ## + ethnicity             2   50 1.4847  0.075 . 
    ## + diet_TotFib_g         1   50 1.2106  0.200   
    ## + diet_MonSac_g         1   50 1.2001  0.240   
    ## + diet_Alc_g            1   50 1.1696  0.240   
    ## + relation              3   52 1.1208  0.265   
    ## + diet_Iron_mg          1   50 1.2390  0.285   
    ## + diet_Copp_mg          1   50 1.1097  0.345   
    ## + diet_Fol.DFE_mcg_DFE  1   50 1.0550  0.385   
    ## + diet_BetaCaro_mcg     1   50 1.0870  0.395   
    ## + diet_Vit_B6_mg        1   50 1.0688  0.405   
    ## + bmi                   1   51 0.9973  0.460   
    ## + age                   1   51 0.9032  0.580   
    ## + diet_Disacc_g         1   51 0.6028  0.880   
    ## + diet_Mang_mg          1   51 0.5408  0.945   
    ## + diet_Vit_B1_mg        1   51 0.5112  0.950   
    ## + cohort                1   51 0.4439  0.980   
    ## + subject_id           35 -Inf                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Step: MB_asvtab_brocc ~ group + diet_TotSolFib_g 
    ## 
    ##                    Df    AIC      F Pr(>F)   
    ## - diet_TotSolFib_g  1 49.659 2.0965  0.035 * 
    ## - group             2 61.585 9.0370  0.005 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##                          Df  AIC      F Pr(>F)  
    ## + diet_SatFat_g           1   49 1.8331  0.025 *
    ## + diet_OCarb_g            1   49 1.8757  0.040 *
    ## + ethnicity               2   50 1.4946  0.050 *
    ## + sex                     1   50 1.6509  0.055 .
    ## + weight_lb               1   49 1.7374  0.060 .
    ## + race                    6   52 1.3075  0.070 .
    ## + height_in               1   50 1.4035  0.090 .
    ## + diet_Chol_mg            1   50 1.3722  0.110  
    ## + diet_condensed          4   52 1.2436  0.115  
    ## + diet_TotInsolFib_g      1   50 1.3854  0.185  
    ## + relation                3   51 1.1210  0.280  
    ## + diet_Vit_B6_mg          1   50 1.1248  0.345  
    ## + diet_BetaCaro_mcg       1   50 1.0986  0.345  
    ## + bmi                     1   50 1.0245  0.455  
    ## + diet_Iron_mg            1   50 0.9072  0.525  
    ## + diet_Alc_g              1   50 0.9641  0.540  
    ## + diet_MonSac_g           1   50 0.9331  0.555  
    ## + diet_Fol.DFE_mcg_DFE    1   50 0.8356  0.600  
    ## + age                     1   50 0.8514  0.660  
    ## + diet_TotFib_g           1   51 0.7662  0.720  
    ## + diet_Copp_mg            1   51 0.6830  0.820  
    ## + diet_Mang_mg            1   51 0.6584  0.865  
    ## + diet_Disacc_g           1   51 0.6045  0.905  
    ## + cohort                  1   51 0.5523  0.905  
    ## + diet_Vit_B1_mg          1   51 0.5272  0.955  
    ## + group:diet_TotSolFib_g  2   52 0.5973  0.955  
    ## + subject_id             34 -Inf                
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Step: MB_asvtab_brocc ~ group + diet_TotSolFib_g + diet_SatFat_g 
    ## 
    ##                    Df    AIC      F Pr(>F)   
    ## - diet_SatFat_g     1 49.386 1.8331  0.045 * 
    ## - diet_TotSolFib_g  1 49.614 2.0426  0.010 **
    ## - group             2 61.469 8.7300  0.005 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##                                  Df  AIC      F Pr(>F)  
    ## + weight_lb                       1   49 1.7870  0.015 *
    ## + diet_SatFat_g:diet_TotSolFib_g  1   49 1.7423  0.035 *
    ## + sex                             1   49 1.6796  0.035 *
    ## + ethnicity                       2   50 1.5422  0.050 *
    ## + race                            6   51 1.3397  0.065 .
    ## + height_in                       1   50 1.4426  0.100 .
    ## + diet_TotInsolFib_g              1   50 1.5388  0.105  
    ## + diet_condensed                  4   52 1.1958  0.140  
    ## + diet_Vit_B6_mg                  1   50 1.3569  0.145  
    ## + diet_Chol_mg                    1   50 1.3550  0.165  
    ## + diet_TotFib_g                   1   50 1.3275  0.175  
    ## + relation                        3   51 1.1473  0.220  
    ## + diet_OCarb_g                    1   50 1.2055  0.265  
    ## + diet_MonSac_g                   1   50 1.2612  0.275  
    ## + diet_BetaCaro_mcg               1   50 1.1673  0.290  
    ## + diet_Fol.DFE_mcg_DFE            1   50 1.0638  0.370  
    ## + bmi                             1   50 1.0668  0.410  
    ## + diet_Copp_mg                    1   50 1.0326  0.440  
    ## + age                             1   50 1.0006  0.450  
    ## + diet_Alc_g                      1   50 0.9661  0.500  
    ## + diet_Mang_mg                    1   50 0.8667  0.550  
    ## + diet_Iron_mg                    1   50 0.8646  0.555  
    ## + group:diet_SatFat_g             2   51 0.8856  0.640  
    ## + diet_Disacc_g                   1   51 0.6411  0.810  
    ## + group:diet_TotSolFib_g          2   52 0.6786  0.905  
    ## + diet_Vit_B1_mg                  1   51 0.5897  0.910  
    ## + cohort                          1   51 0.5501  0.990  
    ## + subject_id                     33 -Inf                
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Step: MB_asvtab_brocc ~ group + diet_TotSolFib_g + diet_SatFat_g +      weight_lb 
    ## 
    ##                    Df    AIC      F Pr(>F)   
    ## - weight_lb         1 49.331 1.7870  0.085 . 
    ## - diet_SatFat_g     1 49.436 1.8800  0.035 * 
    ## - diet_TotSolFib_g  1 49.612 2.0377  0.035 * 
    ## - group             2 61.768 8.7008  0.005 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##                                  Df  AIC      F Pr(>F)  
    ## + diet_SatFat_g:diet_TotSolFib_g  1   49 1.7711  0.040 *
    ## + sex                             1   49 1.5651  0.070 .
    ## + weight_lb:group                 2   50 1.3377  0.075 .
    ## + race                            6   51 1.3213  0.080 .
    ## + diet_Chol_mg                    1   50 1.3909  0.095 .
    ## + ethnicity                       2   50 1.4585  0.100 .
    ## + diet_condensed                  4   51 1.2724  0.110  
    ## + diet_TotInsolFib_g              1   49 1.5365  0.120  
    ## + diet_TotFib_g                   1   50 1.3670  0.140  
    ## + diet_Vit_B6_mg                  1   50 1.3353  0.165  
    ## + diet_OCarb_g                    1   50 1.2669  0.200  
    ## + relation                        3   51 1.1940  0.205  
    ## + diet_MonSac_g                   1   50 1.2289  0.270  
    ## + group:diet_SatFat_g             2   51 1.0581  0.305  
    ## + diet_BetaCaro_mcg               1   50 1.0910  0.350  
    ## + diet_Fol.DFE_mcg_DFE            1   50 1.0301  0.360  
    ## + age                             1   50 1.0139  0.365  
    ## + weight_lb:diet_SatFat_g         1   50 1.0662  0.380  
    ## + diet_Copp_mg                    1   50 1.0479  0.415  
    ## + diet_Alc_g                      1   50 0.9780  0.475  
    ## + weight_lb:diet_TotSolFib_g      1   50 0.8617  0.550  
    ## + diet_Mang_mg                    1   50 0.8696  0.555  
    ## + diet_Iron_mg                    1   50 0.8615  0.635  
    ## + height_in                       1   50 0.6966  0.840  
    ## + diet_Disacc_g                   1   51 0.6316  0.860  
    ## + bmi                             1   50 0.6534  0.865  
    ## + group:diet_TotSolFib_g          2   52 0.6961  0.880  
    ## + cohort                          1   51 0.5666  0.930  
    ## + diet_Vit_B1_mg                  1   51 0.4653  0.980  
    ## + subject_id                     32 -Inf                
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Step: MB_asvtab_brocc ~ group + diet_TotSolFib_g + diet_SatFat_g +      weight_lb + diet_TotSolFib_g:diet_SatFat_g 
    ## 
    ##                                  Df    AIC      F Pr(>F)   
    ## - diet_TotSolFib_g:diet_SatFat_g  1 49.267 1.7711  0.070 . 
    ## - weight_lb                       1 49.317 1.8145  0.055 . 
    ## - group                           2 62.010 8.6525  0.005 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##                              Df  AIC      F Pr(>F)  
    ## + race                        6   50 1.3967  0.045 *
    ## + diet_TotInsolFib_g          1   49 1.6146  0.085 .
    ## + diet_Vit_B6_mg              1   49 1.3966  0.095 .
    ## + diet_condensed              4   51 1.2585  0.095 .
    ## + weight_lb:group             2   50 1.3812  0.100 .
    ## + sex                         1   49 1.5684  0.115  
    ## + ethnicity                   2   49 1.5077  0.120  
    ## + diet_Copp_mg                1   50 1.2725  0.170  
    ## + diet_MonSac_g               1   50 1.2962  0.185  
    ## + diet_TotFib_g               1   50 1.3017  0.190  
    ## + diet_OCarb_g                1   50 1.2920  0.190  
    ## + diet_Chol_mg                1   50 1.2602  0.240  
    ## + weight_lb:diet_SatFat_g     1   50 1.1750  0.255  
    ## + group:diet_SatFat_g         2   50 1.1157  0.285  
    ## + relation                    3   51 1.2137  0.290  
    ## + diet_Fol.DFE_mcg_DFE        1   50 1.2033  0.300  
    ## + diet_BetaCaro_mcg           1   50 0.9112  0.440  
    ## + diet_Alc_g                  1   50 0.9827  0.460  
    ## + diet_Iron_mg                1   50 0.9179  0.505  
    ## + diet_Mang_mg                1   50 0.9285  0.515  
    ## + weight_lb:diet_TotSolFib_g  1   50 0.9121  0.555  
    ## + age                         1   50 0.9178  0.565  
    ## + height_in                   1   50 0.6874  0.795  
    ## + diet_Disacc_g               1   50 0.6323  0.835  
    ## + group:diet_TotSolFib_g      2   51 0.7747  0.860  
    ## + cohort                      1   50 0.6558  0.875  
    ## + bmi                         1   50 0.6399  0.880  
    ## + diet_Vit_B1_mg              1   51 0.4783  0.975  
    ## + subject_id                 31 -Inf                
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Step: MB_asvtab_brocc ~ group + diet_TotSolFib_g + diet_SatFat_g +      weight_lb + race + diet_TotSolFib_g:diet_SatFat_g 
    ## 
    ##                                  Df    AIC      F Pr(>F)   
    ## - weight_lb                       1 50.602 1.6523  0.080 . 
    ## - race                            6 49.155 1.3967  0.045 * 
    ## - diet_TotSolFib_g:diet_SatFat_g  1 51.154 2.0419  0.045 * 
    ## - group                           2 65.632 8.3611  0.005 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##                              Df  AIC      F Pr(>F)  
    ## + diet_Vit_B6_mg              1   49 1.9020  0.025 *
    ## + diet_TotInsolFib_g          1   50 1.6588  0.055 .
    ## + race:diet_SatFat_g          3   50 1.3420  0.100 .
    ## + diet_condensed              4   50 1.2754  0.105  
    ## + weight_lb:group             2   50 1.3570  0.105  
    ## + diet_Iron_mg                1   50 1.5506  0.115  
    ## + relation                    3   50 1.2662  0.120  
    ## + diet_Vit_B1_mg              1   50 1.3923  0.140  
    ## + race:weight_lb              3   50 1.3331  0.150  
    ## + sex                         1   50 1.4279  0.150  
    ## + race:diet_TotSolFib_g       3   50 1.2025  0.185  
    ## + diet_Fol.DFE_mcg_DFE        1   50 1.2648  0.230  
    ## + diet_Copp_mg                1   50 1.1643  0.280  
    ## + group:diet_SatFat_g         2   51 1.1193  0.280  
    ## + diet_Mang_mg                1   50 1.2006  0.285  
    ## + diet_MonSac_g               1   51 1.0687  0.285  
    ## + diet_Chol_mg                1   50 1.1593  0.305  
    ## + diet_Disacc_g               1   50 1.0994  0.345  
    ## + weight_lb:diet_SatFat_g     1   51 1.0308  0.360  
    ## + race:group                  3   51 1.0807  0.360  
    ## + weight_lb:diet_TotSolFib_g  1   50 1.0986  0.390  
    ## + diet_OCarb_g                1   51 1.0186  0.400  
    ## + diet_TotFib_g               1   50 1.0949  0.410  
    ## + ethnicity                   1   51 0.7745  0.645  
    ## + age                         1   51 0.7235  0.745  
    ## + bmi                         1   51 0.7522  0.770  
    ## + diet_Alc_g                  1   51 0.6902  0.785  
    ## + height_in                   1   51 0.6720  0.840  
    ## + group:diet_TotSolFib_g      2   52 0.6629  0.890  
    ## + diet_BetaCaro_mcg           1   51 0.5518  0.910  
    ## + cohort                      1   51 0.5120  0.945  
    ## + subject_id                 25 -Inf                
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Step: MB_asvtab_brocc ~ group + diet_TotSolFib_g + diet_SatFat_g +      weight_lb + race + diet_Vit_B6_mg + diet_TotSolFib_g:diet_SatFat_g 
    ## 
    ##                                  Df    AIC      F Pr(>F)   
    ## - weight_lb                       1 49.784 1.6403  0.105   
    ## - diet_TotSolFib_g:diet_SatFat_g  1 50.498 2.1267  0.050 * 
    ## - diet_Vit_B6_mg                  1 50.170 1.9020  0.045 * 
    ## - race                            6 49.426 1.5076  0.020 * 
    ## - group                           2 65.272 8.3121  0.005 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Step: MB_asvtab_brocc ~ group + diet_TotSolFib_g + diet_SatFat_g +      race + diet_Vit_B6_mg + diet_TotSolFib_g:diet_SatFat_g
    ## 
    ##                                   Df  AIC      F Pr(>F)  
    ## + race:diet_Vit_B6_mg              3   48 1.5847  0.020 *
    ## + race:diet_SatFat_g               3   49 1.4777  0.040 *
    ## + diet_Iron_mg                     1   49 1.6670  0.050 *
    ## + sex                              1   49 1.6995  0.060 .
    ## + diet_TotInsolFib_g               1   49 1.7021  0.075 .
    ## + weight_lb                        1   49 1.6403  0.085 .
    ## + diet_Vit_B1_mg                   1   50 1.4410  0.105  
    ## + diet_Chol_mg                     1   50 1.3520  0.105  
    ## + race:diet_TotSolFib_g            3   50 1.2852  0.110  
    ## + diet_Mang_mg                     1   50 1.4092  0.145  
    ## + diet_condensed                   4   50 1.2134  0.170  
    ## + diet_SatFat_g:diet_Vit_B6_mg     1   50 1.2828  0.180  
    ## + relation                         3   50 1.2338  0.185  
    ## + height_in                        1   50 1.3070  0.200  
    ## + diet_Disacc_g                    1   50 1.3117  0.210  
    ## + diet_TotSolFib_g:diet_Vit_B6_mg  1   50 1.2410  0.235  
    ## + bmi                              1   50 1.1759  0.235  
    ## + diet_MonSac_g                    1   50 1.2183  0.245  
    ## + race:group                       3   50 1.1656  0.245  
    ## + diet_Copp_mg                     1   50 1.1523  0.305  
    ## + diet_OCarb_g                     1   50 0.9823  0.400  
    ## + group:diet_Vit_B6_mg             2   51 1.0004  0.455  
    ## + diet_TotFib_g                    1   50 0.8939  0.520  
    ## + group:diet_SatFat_g              2   51 0.9449  0.540  
    ## + diet_Fol.DFE_mcg_DFE             1   51 0.7885  0.680  
    ## + ethnicity                        1   50 0.8392  0.720  
    ## + age                              1   51 0.7561  0.720  
    ## + group:diet_TotSolFib_g           2   51 0.7355  0.860  
    ## + diet_Alc_g                       1   51 0.5646  0.930  
    ## + diet_BetaCaro_mcg                1   51 0.5413  0.965  
    ## + cohort                           1   51 0.4516  0.985  
    ## + subject_id                      25 -Inf                
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Step: MB_asvtab_brocc ~ group + diet_TotSolFib_g + diet_SatFat_g +      race + diet_Vit_B6_mg + diet_TotSolFib_g:diet_SatFat_g +      race:diet_Vit_B6_mg 
    ## 
    ##                                  Df    AIC      F Pr(>F)   
    ## - race:diet_Vit_B6_mg             3 49.784 1.5847  0.050 * 
    ## - diet_TotSolFib_g:diet_SatFat_g  1 49.633 1.9854  0.035 * 
    ## - group                           2 63.646 7.2781  0.005 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##                                   Df  AIC      F Pr(>F)  
    ## + diet_Chol_mg                     1   47 1.6951  0.095 .
    ## + weight_lb                        1   48 1.4091  0.150  
    ## + diet_Vit_B1_mg                   1   48 1.3126  0.210  
    ## + diet_Mang_mg                     1   48 1.2931  0.245  
    ## + diet_condensed                   4   48 1.1321  0.255  
    ## + diet_SatFat_g:diet_Vit_B6_mg     1   48 1.2312  0.265  
    ## + height_in                        1   48 1.1517  0.310  
    ## + diet_Iron_mg                     1   48 1.1687  0.315  
    ## + sex                              1   48 1.0720  0.360  
    ## + bmi                              1   48 1.0837  0.390  
    ## + diet_TotInsolFib_g               1   49 0.9760  0.420  
    ## + diet_TotSolFib_g:diet_Vit_B6_mg  1   49 1.0065  0.440  
    ## + group:diet_Vit_B6_mg             2   49 1.0096  0.455  
    ## + diet_Copp_mg                     1   49 1.0009  0.470  
    ## + diet_Disacc_g                    1   49 0.9802  0.490  
    ## + race:diet_SatFat_g               1   49 0.9065  0.530  
    ## + ethnicity                        1   49 0.8974  0.530  
    ## + relation                         2   49 0.9131  0.545  
    ## + group:diet_SatFat_g              2   49 0.8740  0.605  
    ## + age                              1   49 0.8046  0.645  
    ## + diet_Fol.DFE_mcg_DFE             1   49 0.7593  0.680  
    ## + diet_BetaCaro_mcg                1   49 0.7937  0.695  
    ## + group:diet_TotSolFib_g           2   49 0.8132  0.730  
    ## + diet_TotFib_g                    1   49 0.7493  0.735  
    ## + race:group                       2   49 0.8068  0.760  
    ## + diet_OCarb_g                     1   49 0.6298  0.860  
    ## + diet_MonSac_g                    1   49 0.5926  0.865  
    ## + race:diet_TotSolFib_g            1   50 0.4556  0.965  
    ## + cohort                           1   49 0.5093  0.985  
    ## + diet_Alc_g                       1   49 0.4794  0.990  
    ## + subject_id                      22 -Inf                
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ordibrocc
    ## Call: capscale(formula = MB_asvtab_brocc ~ group + diet_TotSolFib_g +
    ## diet_SatFat_g + race + diet_Vit_B6_mg + diet_TotSolFib_g:diet_SatFat_g +
    ## race:diet_Vit_B6_mg, data = meta_groups_brocc, distance = "bray")
    ## 
    ## -- Model Summary --
    ## 
    ##               Inertia Proportion Rank
    ## Total          4.6352                
    ## RealTotal      4.8697     1.0000     
    ## Constrained    3.2905     0.6757   15
    ## Unconstrained  1.5792     0.3243   22
    ## Imaginary     -0.2345                
    ## 
    ## Inertia is squared Bray distance
    ## 
    ## -- Note --
    ## 
    ## Species scores projected from 'MB_asvtab_brocc'
    ## 
    ## Some constraints or conditions were aliased because they were redundant.
    ## This can happen if terms are linearly dependent (collinear): 'raceMore than
    ## one race:diet_Vit_B6_mg', 'raceOther:diet_Vit_B6_mg',
    ## 'raceWhite:diet_Vit_B6_mg'
    ## 
    ## -- Eigenvalues --
    ## 
    ## Eigenvalues for constrained axes:
    ##   CAP1   CAP2   CAP3   CAP4   CAP5   CAP6   CAP7   CAP8   CAP9  CAP10  CAP11 
    ## 1.1552 0.6415 0.4029 0.2485 0.1939 0.1580 0.1261 0.0892 0.0749 0.0592 0.0526 
    ##  CAP12  CAP13  CAP14  CAP15 
    ## 0.0375 0.0309 0.0135 0.0067 
    ## 
    ## Eigenvalues for unconstrained axes:
    ##    MDS1    MDS2    MDS3    MDS4    MDS5    MDS6    MDS7    MDS8 
    ## 0.26215 0.21118 0.20354 0.13610 0.11380 0.10303 0.09044 0.08187 
    ## (Showing 8 of 22 unconstrained eigenvalues)
    # best model is MB_asvtab_brocc ~ group + diet_TotSolFib_g + diet_OCarb_g + weight_lb +
    # diet_TotSolFib_g:diet_OCarb_g
    # model inertia was 4.6352; Inertia is squared Bray distance

    # calculate RDA (distance based)
    rda_brocc <- capscale(formula = MB_asvtab_brocc ~ group + diet_TotSolFib_g + diet_OCarb_g + weight_lb +
    diet_TotSolFib_g:diet_OCarb_g, data = meta_groups_brocc, 
                          distance = "bray")
    RsquareAdj(rda_brocc)
    ## $r.squared
    ## [1] 0.4757074
    ## 
    ## $adj.r.squared
    ## [1] 0.3742315
    # 0.374 of variation explained by this model
    summary(rda_brocc)
    ## 
    ## Call:
    ## capscale(formula = MB_asvtab_brocc ~ group + diet_TotSolFib_g +      diet_OCarb_g + weight_lb + diet_TotSolFib_g:diet_OCarb_g,      data = meta_groups_brocc, distance = "bray") 
    ## 
    ## Partitioning of squared Bray distance:
    ##               Inertia Proportion
    ## Total           4.870     1.0000
    ## Constrained     2.317     0.4757
    ## Unconstrained   2.553     0.5243
    ## 
    ## Eigenvalues, and their contribution to the squared Bray distance 
    ## 
    ## Importance of components:
    ##                         CAP1   CAP2    CAP3    CAP4    CAP5    CAP6    MDS1
    ## Eigenvalue            1.0950 0.6021 0.26424 0.15266 0.12610 0.07648 0.41759
    ## Proportion Explained  0.2249 0.1236 0.05426 0.03135 0.02589 0.01571 0.08575
    ## Cumulative Proportion 0.2249 0.3485 0.40276 0.43411 0.46000 0.47571 0.56146
    ##                          MDS2    MDS3    MDS4   MDS5   MDS6    MDS7    MDS8
    ## Eigenvalue            0.34611 0.27766 0.24558 0.1607 0.1578 0.13238 0.11850
    ## Proportion Explained  0.07107 0.05702 0.05043 0.0330 0.0324 0.02718 0.02433
    ## Cumulative Proportion 0.63254 0.68955 0.73998 0.7730 0.8054 0.83257 0.85690
    ##                          MDS9   MDS10   MDS11   MDS12   MDS13   MDS14    MDS15
    ## Eigenvalue            0.10996 0.10329 0.08300 0.07310 0.05578 0.05167 0.043990
    ## Proportion Explained  0.02258 0.02121 0.01704 0.01501 0.01146 0.01061 0.009033
    ## Cumulative Proportion 0.87948 0.90069 0.91774 0.93275 0.94420 0.95482 0.963849
    ##                          MDS16    MDS17   MDS18    MDS19    MDS20    MDS21
    ## Eigenvalue            0.034985 0.034376 0.02264 0.020635 0.016427 0.015686
    ## Proportion Explained  0.007184 0.007059 0.00465 0.004237 0.003373 0.003221
    ## Cumulative Proportion 0.971033 0.978092 0.98274 0.986979 0.990352 0.993574
    ##                          MDS22    MDS23    MDS24     MDS25     MDS26     MDS27
    ## Eigenvalue            0.012577 0.007638 0.005103 0.0035069 0.0021824 2.871e-04
    ## Proportion Explained  0.002583 0.001568 0.001048 0.0007202 0.0004481 5.896e-05
    ## Cumulative Proportion 0.996156 0.997725 0.998773 0.9994929 0.9999410 1.000e+00
    ## 
    ## Accumulated constrained eigenvalues
    ## Importance of components:
    ##                         CAP1   CAP2   CAP3   CAP4    CAP5    CAP6
    ## Eigenvalue            1.0950 0.6021 0.2642 0.1527 0.12610 0.07648
    ## Proportion Explained  0.4727 0.2599 0.1141 0.0659 0.05443 0.03301
    ## Cumulative Proportion 0.4727 0.7326 0.8466 0.9126 0.96699 1.00000

    plot(rda_brocc)

<img src="broccolipreMBpermanova_files/figure-markdown_strict/broccoli dbRDA run with best model from ordistep-1.png" width="98%" height="98%" />

    dm_brocc <- vegdist(MB_asvtab_brocc, method = "bray")
    Perm_brocc <- adonis2(dm_brocc ~ group + diet_TotSolFib_g + diet_OCarb_g + weight_lb +
    diet_TotSolFib_g:diet_OCarb_g, 
                                data = meta_groups_brocc)
    Perm_brocc
    ## Permutation test for adonis under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = dm_brocc ~ group + diet_TotSolFib_g + diet_OCarb_g + weight_lb + diet_TotSolFib_g:diet_OCarb_g, data = meta_groups_brocc)
    ##          Df SumOfSqs      R2      F Pr(>F)    
    ## Model     6   2.3012 0.49647 5.0941  0.001 ***
    ## Residual 31   2.3340 0.50353                  
    ## Total    37   4.6352 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    p.adjust(Perm_brocc$`Pr(>F)`, "fdr")
    ## [1] 0.001    NA    NA

    # write.table(Perm_brocc,
    #             "../outputs/Perm_brocc.txt",
    #             sep = "\t")

    # create files for graphing
    smry_rda_brocc <- summary(rda_brocc)

    # this is the coordinates for the points and metadata
    meta_groups_filt_brocc <- meta_groups_brocc %>%
      dplyr::select(c(subject_id, group, diet_OCarb_g, diet_TotSolFib_g, weight_lb))
    brocc_PC1  <- data.frame(smry_rda_brocc$sites[,1:2]) %>%  # these are the x, y coordinates for the sample points
      rownames_to_column("subject_id") %>%
      inner_join(meta_groups_filt_brocc, by = "subject_id") %>%
      column_to_rownames("subject_id")
    mylims <- range(with(brocc_PC1, c(CAP1, CAP2)))
    brocc_PC1$group <- as.factor(brocc_PC1$group)

    # the biplot scores are the correlations between your environmental variables and axes
    brocc_PC2  <- data.frame(smry_rda_brocc$biplot[1:5,])

    # put plot together
    rda_plot_brocc <- ggplot(brocc_PC1, aes(x = CAP1, y = CAP2)) + 
      geom_point(aes(color = group), size = 2) +
      theme_classic() +
      geom_hline(yintercept = 0, linetype = "dotted") +
      geom_vline(xintercept = 0, linetype = "dotted") +
      geom_segment(data = brocc_PC2[c(3:5),], aes(x = 0, xend = CAP1, y = 0, yend = CAP2), 
                   color = "black", arrow = arrow(length = unit(0.01, "npc"))) +
      geom_text(data = brocc_PC2[c(3:5),], 
                aes(x = CAP1, y = CAP2), label = c("Total Soluble Fiber", "Other Carbs", "Weight (lb)"), 
                color = "black", size = 4, hjust = 1, vjust = 0.5) +
      labs(x = paste0("CAP1 (",round(100*smry_rda$cont$importance[2, "CAP1"], digits = 2),"%)"),
           y = paste0("CAP2 (",round(100*smry_rda$cont$importance[2, "CAP2"], digits = 2),"%)")) +  
      scale_color_brewer(palette = "Set1",
                        name = "pre-treatment\nmicrobiome group") +
      theme(text = element_text(size = 20))
    rda_plot_brocc

<img src="broccolipreMBpermanova_files/figure-markdown_strict/broccoli rda graph-1.png" width="98%" height="98%" />


    # ggsave("../outputs/rda_plot_brocc_groupsinformula.png",
    #        width = 10,
    #        height = 7,
    #        dpi = 300)


    ## simple graph
    rda_plot_plain <- ggplot(brocc_PC1, aes(x = CAP1, y = CAP2)) + 
      geom_point(aes(color = group), size = 2) +
      theme_classic() +
      geom_hline(yintercept = 0, linetype = "dotted") +
      geom_vline(xintercept = 0, linetype = "dotted") +
      scale_color_brewer(palette = "Set1",
                        name = "pre-treatment\nmicrobiome group") +
      theme(text = element_text(size = 20))
    rda_plot_plain

<img src="broccolipreMBpermanova_files/figure-markdown_strict/broccoli rda graph-2.png" width="98%" height="98%" />


    # ggsave("../outputs/rda_plot_plain.png",
    #        width = 10,
    #        height = 7,
    #        dpi = 300)

## redo above but without group - this is for the main manuscript

-   broccoli only
-   remove group since it might be “stealing” some of the variation in
    BC ordination space
-   Sample issue removed; subject 40 is on the edge of the ordination
    space and was artificially stealing some variation as such but did
    not have any reason to be noted and isn’t truly an outlier

<!-- -->

    meta_nogroups_brocc <- meta_groups_brocc %>%
      dplyr::select(-c(group, ethnicity, relation, height_in, weight_lb)) %>% # removed colinear ones as well because they are creating overfitting in model
      column_to_rownames("subject_id")
    broccmod0 <- capscale(MB_asvtab_brocc ~ 1, meta_nogroups_brocc, distance = "bray")  # Model with intercept only
    broccmod1 <- capscale(MB_asvtab_brocc ~ . , meta_nogroups_brocc, distance = "bray")  # Model with all explanatory variables
    ordibrocc <- ordistep(broccmod0, scope = formula(broccmod1), permutations = how(nperm = 499)) # this determines what the best model is to run RDA on
    ## 
    ## Start: MB_asvtab_brocc ~ 1 
    ## 
    ##                        Df    AIC      F Pr(>F)   
    ## + diet_TotInsolFib_g    1 60.522 2.5693  0.008 **
    ## + diet_OCarb_g          1 61.013 2.0744  0.024 * 
    ## + diet_condensed        4 62.506 1.5741  0.024 * 
    ## + age                   1 60.845 2.2433  0.026 * 
    ## + diet_SatFat_g         1 61.015 2.0721  0.036 * 
    ## + sex                   1 61.222 1.8656  0.054 . 
    ## + race                  6 64.734 1.2796  0.110   
    ## + diet_TotSolFib_g      1 61.585 1.5053  0.126   
    ## + diet_Fol.DFE_mcg_DFE  1 61.703 1.3890  0.170   
    ## + diet_MonSac_g         1 61.823 1.2712  0.196   
    ## + diet_Chol_mg          1 61.805 1.2888  0.238   
    ## + diet_Copp_mg          1 61.986 1.1123  0.284   
    ## + diet_BetaCaro_mcg     1 61.976 1.1219  0.302   
    ## + diet_Mang_mg          1 62.066 1.0341  0.378   
    ## + diet_TotFib_g         1 62.137 0.9651  0.434   
    ## + diet_Iron_mg          1 62.211 0.8929  0.540   
    ## + diet_Alc_g            1 62.251 0.8545  0.558   
    ## + bmi                   1 62.405 0.7046  0.752   
    ## + diet_Vit_B6_mg        1 62.417 0.6933  0.772   
    ## + diet_Vit_B1_mg        1 62.473 0.6390  0.790   
    ## + diet_Disacc_g         1 62.545 0.5700  0.882   
    ## + cohort                1 62.615 0.5025  0.952   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Step: MB_asvtab_brocc ~ diet_TotInsolFib_g 
    ## 
    ##                      Df    AIC      F Pr(>F)   
    ## - diet_TotInsolFib_g  1 61.142 2.5693   0.01 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##                        Df    AIC      F Pr(>F)  
    ## + diet_OCarb_g          1 60.009 2.3927  0.016 *
    ## + diet_SatFat_g         1 60.041 2.3615  0.020 *
    ## + age                   1 60.289 2.1189  0.020 *
    ## + diet_condensed        4 62.424 1.3927  0.056 .
    ## + sex                   1 60.737 1.6833  0.094 .
    ## + race                  6 63.792 1.2914  0.104  
    ## + diet_TotSolFib_g      1 61.019 1.4127  0.142  
    ## + diet_Chol_mg          1 61.203 1.2362  0.262  
    ## + diet_BetaCaro_mcg     1 61.331 1.1145  0.316  
    ## + diet_Mang_mg          1 61.344 1.1020  0.328  
    ## + bmi                   1 61.615 0.8453  0.590  
    ## + diet_Alc_g            1 61.596 0.8639  0.614  
    ## + cohort                1 61.629 0.8327  0.652  
    ## + diet_Iron_mg          1 61.700 0.7660  0.658  
    ## + diet_Vit_B1_mg        1 61.682 0.7829  0.660  
    ## + diet_Fol.DFE_mcg_DFE  1 61.684 0.7803  0.664  
    ## + diet_Vit_B6_mg        1 61.720 0.7471  0.704  
    ## + diet_Copp_mg          1 61.817 0.6553  0.798  
    ## + diet_MonSac_g         1 61.858 0.6169  0.840  
    ## + diet_TotFib_g         1 61.897 0.5803  0.872  
    ## + diet_Disacc_g         1 61.937 0.5434  0.908  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Step: MB_asvtab_brocc ~ diet_TotInsolFib_g + diet_OCarb_g 
    ## 
    ##                      Df    AIC      F Pr(>F)   
    ## - diet_OCarb_g        1 60.522 2.3927  0.020 * 
    ## - diet_TotInsolFib_g  1 61.013 2.8788  0.002 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##                        Df    AIC      F Pr(>F)  
    ## + sex                   1 60.172 1.6848  0.078 .
    ## + age                   1 60.252 1.6095  0.078 .
    ## + diet_condensed        4 62.256 1.2669  0.104  
    ## + diet_TotSolFib_g      1 60.383 1.4872  0.122  
    ## + race                  6 63.923 1.1462  0.220  
    ## + diet_BetaCaro_mcg     1 60.699 1.1927  0.290  
    ## + cohort                1 60.828 1.0735  0.372  
    ## + diet_SatFat_g         1 60.859 1.0452  0.400  
    ## + diet_Chol_mg          1 60.861 1.0436  0.404  
    ## + diet_Mang_mg          1 60.923 0.9861  0.430  
    ## + bmi                   1 60.964 0.9486  0.506  
    ## + diet_Fol.DFE_mcg_DFE  1 61.111 0.8136  0.580  
    ## + diet_Vit_B6_mg        1 61.134 0.7919  0.642  
    ## + diet_Alc_g            1 61.173 0.7563  0.718  
    ## + diet_TotFib_g         1 61.255 0.6822  0.770  
    ## + diet_MonSac_g         1 61.281 0.6578  0.800  
    ## + diet_Copp_mg          1 61.301 0.6401  0.820  
    ## + diet_Iron_mg          1 61.320 0.6228  0.842  
    ## + diet_Disacc_g         1 61.400 0.5494  0.908  
    ## + diet_Vit_B1_mg        1 61.620 0.3503  0.992  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ordibrocc
    ## Call: capscale(formula = MB_asvtab_brocc ~ diet_TotInsolFib_g +
    ## diet_OCarb_g, data = meta_nogroups_brocc, distance = "bray")
    ## 
    ## -- Model Summary --
    ## 
    ##               Inertia Proportion Rank
    ## Total          4.6352                
    ## RealTotal      4.8697     1.0000     
    ## Constrained    0.6152     0.1263    2
    ## Unconstrained  4.2545     0.8737   27
    ## Imaginary     -0.2345                
    ## 
    ## Inertia is squared Bray distance
    ## 
    ## -- Note --
    ## 
    ## Species scores projected from 'MB_asvtab_brocc'
    ## 
    ## -- Eigenvalues --
    ## 
    ## Eigenvalues for constrained axes:
    ##   CAP1   CAP2 
    ## 0.3792 0.2360 
    ## 
    ## Eigenvalues for unconstrained axes:
    ##   MDS1   MDS2   MDS3   MDS4   MDS5   MDS6   MDS7   MDS8 
    ## 1.0487 0.6276 0.4197 0.3263 0.2557 0.2454 0.1949 0.1607 
    ## (Showing 8 of 27 unconstrained eigenvalues)
    # best model is diet_TotInsolFib_g + diet_OCarb_g
    # model inertia was 4.6352; Inertia is squared Bray distance

    # calculate RDA (distance based) of first model
    rda_brocc <- capscale(formula = MB_asvtab_brocc ~ diet_TotInsolFib_g + diet_OCarb_g, data = meta_nogroups_brocc, 
                          distance = "bray")
    RsquareAdj(rda_brocc)
    ## $r.squared
    ## [1] 0.1263413
    ## 
    ## $adj.r.squared
    ## [1] 0.07641794
    # 0.08 of variation explained by this model

    Perm_brocc <- adonis2(formula = MB_asvtab_brocc ~ diet_TotInsolFib_g + diet_OCarb_g, 
                          data = meta_nogroups_brocc,
                          method = "bray",
                          by = "terms")
    Perm_brocc
    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = MB_asvtab_brocc ~ diet_TotInsolFib_g + diet_OCarb_g, data = meta_nogroups_brocc, method = "bray", by = "terms")
    ##                    Df SumOfSqs      R2      F Pr(>F)   
    ## diet_TotInsolFib_g  1   0.3226 0.06959 2.8057  0.004 **
    ## diet_OCarb_g        1   0.2885 0.06225 2.5096  0.010 **
    ## Residual           35   4.0241 0.86816                 
    ## Total              37   4.6352 1.00000                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    p.adjust(Perm_brocc$`Pr(>F)`, "fdr")
    ## [1] 0.008 0.010    NA    NA

    # write.table(Perm_brocc,
    #             "../outputs/Perm_brocc_nogroups.txt",
    #             sep = "\t")

    # create files for graphing
    smry_rda_brocc <- summary(rda_brocc)

    # this is the coordinates for the points and metadata
    meta_groups_filt_brocc <- meta_groups_brocc %>%
      dplyr::select(c(subject_id, diet_TotInsolFib_g, diet_OCarb_g, group))
    brocc_PC1  <- data.frame(smry_rda_brocc$sites[,1:2]) %>%  # these are the x, y coordinates for the sample points
      rownames_to_column("subject_id") %>%
      inner_join(meta_groups_filt_brocc, by = "subject_id") %>%
      column_to_rownames("subject_id")
    mylims <- range(with(brocc_PC1, c(CAP1, CAP2)))
    brocc_PC1$group <- as.factor(brocc_PC1$group)

    # the biplot scores are the correlations between your environmental variables and axes
    brocc_PC2  <- data.frame(smry_rda_brocc$biplot)

    # create hulls
    # hull_data <- brocc_PC1 %>%
    #   group_by(group) %>%
    #   slice(chull(CAP1, CAP2))

    # put plot together
    rda_plot_brocc <- ggplot(brocc_PC1, aes(x = CAP1, y = CAP2)) + 
      geom_point(aes(color = group), size = 2) +
      theme_classic() +
      geom_hline(yintercept = 0, linetype = "dotted") +
      geom_vline(xintercept = 0, linetype = "dotted") +
      geom_segment(data = brocc_PC2, aes(x = 0, xend = CAP1, y = 0, yend = CAP2), 
                   color = "black", arrow = arrow(length = unit(0.01, "npc"))) +
      geom_text(data = brocc_PC2, 
                aes(x = CAP1, y = CAP2), label = c("Total Insoluble Fiber (g)",
                                                   "Other Carbs (g)"), 
                color = "black", size = 4, hjust = c(0,1), vjust = c(0, 0.5)) +
      # stat_ellipse(aes(linetype = group), position = "identity") +
      labs(x = paste0("CAP1 (",round(100*smry_rda$cont$importance[2, "CAP1"], digits = 2),"%)"),
           y = paste0("CAP2 (",round(100*smry_rda$cont$importance[2, "CAP2"], digits = 2),"%)")) +  
      scale_color_brewer(palette = "Set1",
                        name = "pre-treatment\nmicrobiome group") +
      theme(text = element_text(size = 20))
    rda_plot_brocc

<img src="broccolipreMBpermanova_files/figure-markdown_strict/broccoli rda graph without group-1.png" width="98%" height="98%" />


    # ggsave("../outputs/rda_plot_brocc_V2.png",
    #        width = 9,
    #        height = 6,
    #        dpi = 300)

## what could groups be correlated with?

    # alfalfa and broccoli
    meta_groups_alfbrocc$group_numb[meta_groups_alfbrocc$group == "A"] <- 1
    meta_groups_alfbrocc$group_numb[meta_groups_alfbrocc$group == "B"] <- 2
    names(meta_groups_alfbrocc)
    ##  [1] "sample"               "subject_id"           "time"                
    ##  [4] "treatment"            "cohort"               "sex"                 
    ##  [7] "race"                 "ethnicity"            "age"                 
    ## [10] "height_in"            "weight_lb"            "bmi"                 
    ## [13] "diet"                 "relation"             "sample_issue"        
    ## [16] "diet_condensed"       "veg"                  "broc_consum"         
    ## [19] "groups_alf.cluster"   "group"                "diet_Chol_mg"        
    ## [22] "diet_SatFat_g"        "diet_Alc_g"           "diet_TotSolFib_g"    
    ## [25] "diet_MonSac_g"        "diet_Disacc_g"        "diet_BetaCaro_mcg"   
    ## [28] "diet_Vit_B1_mg"       "diet_Vit_B6_mg"       "diet_Fol.DFE_mcg_DFE"
    ## [31] "diet_Mang_mg"         "diet_Copp_mg"         "diet_TotFib_g"       
    ## [34] "diet_OCarb_g"         "diet_Iron_mg"         "diet_TotInsolFib_g"  
    ## [37] "group_numb"

    lm_alfs <- lm(group_numb ~ sex + age + race + weight_lb + diet_condensed + cohort + ethnicity + height_in + sample_issue + diet_Chol_mg + diet_SatFat_g + diet_Alc_g + diet_TotSolFib_g + diet_MonSac_g + diet_Disacc_g + diet_BetaCaro_mcg + diet_Vit_B1_mg + diet_Vit_B6_mg + diet_Fol.DFE_mcg_DFE + diet_Mang_mg + diet_Copp_mg + diet_TotFib_g + diet_OCarb_g + diet_Iron_mg, 
                  data = meta_groups_alfbrocc)
    summary(lm_alfs)
    ## 
    ## Call:
    ## lm(formula = group_numb ~ sex + age + race + weight_lb + diet_condensed + 
    ##     cohort + ethnicity + height_in + sample_issue + diet_Chol_mg + 
    ##     diet_SatFat_g + diet_Alc_g + diet_TotSolFib_g + diet_MonSac_g + 
    ##     diet_Disacc_g + diet_BetaCaro_mcg + diet_Vit_B1_mg + diet_Vit_B6_mg + 
    ##     diet_Fol.DFE_mcg_DFE + diet_Mang_mg + diet_Copp_mg + diet_TotFib_g + 
    ##     diet_OCarb_g + diet_Iron_mg, data = meta_groups_alfbrocc)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.49939 -0.17514 -0.01504  0.10269  0.71340 
    ## 
    ## Coefficients: (1 not defined because of singularities)
    ##                              Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)                -1.641e+00  1.774e+00  -0.925   0.3615  
    ## sexM                        9.924e-02  1.581e-01   0.628   0.5346  
    ## age                         3.304e-06  6.359e-03   0.001   0.9996  
    ## raceAsian                   1.587e+00  8.646e-01   1.835   0.0755 .
    ## raceBlack/African American  1.354e+00  9.826e-01   1.378   0.1776  
    ## raceDecline                 8.846e-01  9.698e-01   0.912   0.3683  
    ## raceDecline to state        2.244e+00  1.154e+00   1.944   0.0605 .
    ## raceMore than one race      5.830e-01  1.078e+00   0.541   0.5925  
    ## raceMTOR                    1.823e+00  9.864e-01   1.848   0.0736 .
    ## raceOther                   2.438e+00  1.038e+00   2.349   0.0250 *
    ## raceWhite                   1.465e+00  8.537e-01   1.716   0.0956 .
    ## raceWhite (Ashkanazi Jew)   1.243e+00  9.400e-01   1.323   0.1950  
    ## weight_lb                  -2.917e-03  3.594e-03  -0.812   0.4228  
    ## diet_condensedGluten_Free  -2.732e-01  3.110e-01  -0.878   0.3861  
    ## diet_condensedHigh_Protein -1.291e-01  3.573e-01  -0.361   0.7201  
    ## diet_condensedOmnivore     -2.113e-02  1.713e-01  -0.123   0.9026  
    ## diet_condensedPescetarian  -8.320e-02  3.084e-01  -0.270   0.7890  
    ## diet_condensedVegetarian    3.713e-01  2.459e-01   1.510   0.1406  
    ## cohort                     -1.209e-02  2.818e-02  -0.429   0.6706  
    ## ethnicityHL                        NA         NA      NA       NA  
    ## ethnicityNHL                9.577e-01  4.885e-01   1.961   0.0584 .
    ## height_in                   1.025e-02  2.364e-02   0.434   0.6673  
    ## sample_issueTRUE            6.323e-01  8.024e-01   0.788   0.4363  
    ## diet_Chol_mg               -1.278e+00  1.195e+00  -1.069   0.2929  
    ## diet_SatFat_g               1.996e+01  2.092e+01   0.954   0.3469  
    ## diet_Alc_g                 -5.774e+00  1.846e+01  -0.313   0.7564  
    ## diet_TotSolFib_g            6.692e+01  2.067e+02   0.324   0.7482  
    ## diet_MonSac_g              -9.376e+00  1.629e+01  -0.576   0.5688  
    ## diet_Disacc_g               9.206e+01  4.132e+01   2.228   0.0328 *
    ## diet_BetaCaro_mcg           6.519e-02  4.092e-02   1.593   0.1206  
    ## diet_Vit_B1_mg              8.290e+02  5.304e+02   1.563   0.1276  
    ## diet_Vit_B6_mg              6.071e+02  3.874e+02   1.567   0.1267  
    ## diet_Fol.DFE_mcg_DFE       -1.186e+00  1.496e+00  -0.793   0.4334  
    ## diet_Mang_mg                3.301e+02  1.492e+02   2.212   0.0340 *
    ## diet_Copp_mg               -4.211e+02  9.185e+02  -0.459   0.6496  
    ## diet_TotFib_g              -4.563e+01  2.863e+01  -1.594   0.1206  
    ## diet_OCarb_g               -7.609e-01  5.049e+00  -0.151   0.8811  
    ## diet_Iron_mg               -4.171e+01  5.302e+01  -0.787   0.4371  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.3578 on 33 degrees of freedom
    ## Multiple R-squared:  0.5072, Adjusted R-squared:  -0.0304 
    ## F-statistic: 0.9435 on 36 and 33 DF,  p-value: 0.5693

    # raceAsian                   1.587e+00  8.646e-01   1.835   0.0755 .
    # raceDecline to state        2.244e+00  1.154e+00   1.944   0.0605 .
    # raceMTOR                    1.823e+00  9.864e-01   1.848   0.0736 .
    # raceOther                   2.438e+00  1.038e+00   2.349   0.0250 *
    # raceWhite                   1.465e+00  8.537e-01   1.716   0.0956 .
    # ethnicityNHL                9.577e-01  4.885e-01   1.961   0.0584 .
    # diet_Disacc_g               9.206e+01  4.132e+01   2.228   0.0328 *
    # diet_Mang_mg                3.301e+02  1.492e+02   2.212   0.0340 *

    # broccoli only
    meta_groups_brocc$group_numb[meta_groups_brocc$group == "A"] <- 1
    meta_groups_brocc$group_numb[meta_groups_brocc$group == "B"] <- 2
    meta_groups_brocc$group_numb[meta_groups_brocc$group == "C"] <- 3

    lm_brocc <- lm(group_numb ~ sex + age + race + weight_lb + diet_condensed + cohort + ethnicity + height_in + diet_Chol_mg + diet_SatFat_g + diet_Alc_g + diet_TotSolFib_g + diet_MonSac_g + diet_Disacc_g + diet_BetaCaro_mcg + diet_Vit_B1_mg + diet_Vit_B6_mg + diet_Fol.DFE_mcg_DFE + diet_Mang_mg + diet_Copp_mg + diet_TotFib_g + diet_OCarb_g + diet_Iron_mg, 
                  data = meta_groups_brocc)
    sum <- summary(lm_brocc)

    # age                         5.258e-02  1.609e-02   3.267  0.01710 *
    # raceAsian                   4.083e+00  1.782e+00   2.292  0.06180 .
    # raceMTOR                    4.792e+00  1.860e+00   2.577  0.04196 *
    # raceWhite                   4.484e+00  1.752e+00   2.560  0.04293 *
    # weight_lb                  -1.590e-02  7.742e-03  -2.053  0.08583 .
    # diet_condensedPescetarian  -2.774e+00  7.440e-01  -3.729  0.00975 **
    # diet_Alc_g                 -1.189e+02  5.092e+01  -2.336  0.05818 .
    # diet_TotSolFib_g           -2.037e+03  6.975e+02  -2.920  0.02663 *
    # diet_Disacc_g               2.432e+02  8.086e+01   3.008  0.02377 *
    # diet_TotFib_g              -1.486e+02  6.761e+01  -2.198  0.07028 .
    write.table(sum$coefficients, "../outputs/Pred2_grpLMs.txt")

## Permanova on each group

# group A

    meta_groupA_brocc <- meta_groups_brocc %>%
      dplyr::filter(group == "A") %>%
      dplyr::select(-c(group, group_numb, ethnicity, relation, height_in, weight_lb)) %>%
      column_to_rownames("subject_id")
    MB_asvtab_broccA <- MB_asvtab_brocc %>%
      rownames_to_column("subject_id") %>%
      dplyr::filter(subject_id %in% rownames(meta_groupA_brocc)) %>%
      column_to_rownames("subject_id")
    Abroccmod0 <- capscale(MB_asvtab_broccA ~ 1, meta_groupA_brocc, distance = "bray")  # Model with intercept only
    Abroccmod1 <- capscale(MB_asvtab_broccA ~ . + .*., meta_groupA_brocc, distance = "bray")  # Model with all explanatory variables
    ordibroccA <- ordistep(Abroccmod0, scope = formula(Abroccmod1)) # this determines what the best model is to run RDA on
    ## 
    ## Start: MB_asvtab_broccA ~ 1 
    ## 
    ##                        Df     AIC      F Pr(>F)   
    ## + diet_Alc_g            1 -4.1358 1.9915  0.010 **
    ## + bmi                   1 -4.1330 1.9888  0.020 * 
    ## + diet_Chol_mg          1 -3.6197 1.4860  0.105   
    ## + diet_Iron_mg          1 -3.5562 1.4251  0.135   
    ## + diet_OCarb_g          1 -3.4778 1.3503  0.230   
    ## + sex                   1 -3.3616 1.2404  0.240   
    ## + diet_BetaCaro_mcg     1 -3.3684 1.2469  0.275   
    ## + diet_MonSac_g         1 -3.3247 1.2057  0.300   
    ## + diet_condensed        4 -1.5627 1.0745  0.365   
    ## + diet_Disacc_g         1 -3.3591 1.2381  0.410   
    ## + diet_Fol.DFE_mcg_DFE  1 -3.1318 1.0260  0.420   
    ## + race                  4 -1.4250 1.0422  0.420   
    ## + diet_SatFat_g         1 -2.9318 0.8424  0.600   
    ## + diet_Copp_mg          1 -2.9290 0.8398  0.610   
    ## + diet_TotFib_g         1 -2.9682 0.8756  0.650   
    ## + diet_TotSolFib_g      1 -2.8839 0.7989  0.675   
    ## + diet_TotInsolFib_g    1 -2.8771 0.7927  0.690   
    ## + diet_Mang_mg          1 -2.8517 0.7697  0.735   
    ## + age                   1 -2.8601 0.7772  0.750   
    ## + cohort                1 -2.8012 0.7240  0.780   
    ## + diet_Vit_B1_mg        1 -2.6157 0.5579  0.870   
    ## + diet_Vit_B6_mg        1 -2.4838 0.4412  0.960   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Step: MB_asvtab_broccA ~ diet_Alc_g 
    ## 
    ##              Df     AIC      F Pr(>F)  
    ## - diet_Alc_g  1 -3.9725 1.9915   0.03 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##                        Df     AIC      F Pr(>F)  
    ## + bmi                   1 -4.7009 2.1813  0.020 *
    ## + sex                   1 -3.7142 1.2910  0.200  
    ## + diet_BetaCaro_mcg     1 -3.6572 1.2416  0.235  
    ## + diet_Iron_mg          1 -3.5786 1.1738  0.250  
    ## + diet_TotSolFib_g      1 -3.5145 1.1188  0.295  
    ## + race                  4 -2.3632 1.0754  0.360  
    ## + diet_TotInsolFib_g    1 -3.4078 1.0280  0.420  
    ## + diet_condensed        4 -2.0980 1.0183  0.450  
    ## + diet_Mang_mg          1 -3.3173 0.9514  0.480  
    ## + diet_Copp_mg          1 -3.3694 0.9954  0.485  
    ## + diet_TotFib_g         1 -3.2877 0.9266  0.515  
    ## + diet_Chol_mg          1 -3.3414 0.9718  0.530  
    ## + diet_Disacc_g         1 -3.2191 0.8690  0.540  
    ## + diet_Fol.DFE_mcg_DFE  1 -3.2963 0.9338  0.555  
    ## + diet_OCarb_g          1 -3.2320 0.8798  0.650  
    ## + diet_MonSac_g         1 -3.0663 0.7421  0.770  
    ## + age                   1 -3.0237 0.7069  0.790  
    ## + diet_SatFat_g         1 -2.9919 0.6807  0.840  
    ## + diet_Vit_B1_mg        1 -2.7786 0.5070  0.900  
    ## + cohort                1 -2.8333 0.5512  0.925  
    ## + diet_Vit_B6_mg        1 -2.7271 0.4654  0.960  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Step: MB_asvtab_broccA ~ diet_Alc_g + bmi 
    ## 
    ##              Df     AIC      F Pr(>F)  
    ## - diet_Alc_g  1 -4.1330 2.1839  0.025 *
    ## - bmi         1 -4.1358 2.1813  0.015 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##                        Df     AIC      F Pr(>F)  
    ## + diet_TotInsolFib_g    1 -4.7703 1.5530  0.065 .
    ## + diet_BetaCaro_mcg     1 -4.4779 1.3183  0.165  
    ## + diet_TotSolFib_g      1 -4.5297 1.3594  0.190  
    ## + race                  4 -4.0951 1.1492  0.230  
    ## + diet_Iron_mg          1 -4.2668 1.1520  0.325  
    ## + sex                   1 -4.1864 1.0895  0.400  
    ## + diet_Mang_mg          1 -4.1525 1.0632  0.410  
    ## + diet_Chol_mg          1 -4.0456 0.9808  0.445  
    ## + diet_condensed        4 -3.4380 1.0186  0.455  
    ## + diet_Disacc_g         1 -4.0455 0.9807  0.480  
    ## + diet_Fol.DFE_mcg_DFE  1 -4.0341 0.9720  0.505  
    ## + bmi:diet_Alc_g        1 -3.9312 0.8933  0.520  
    ## + diet_TotFib_g         1 -4.0681 0.9981  0.525  
    ## + diet_Copp_mg          1 -4.0568 0.9894  0.530  
    ## + diet_MonSac_g         1 -3.9583 0.9140  0.560  
    ## + diet_Vit_B1_mg        1 -3.8456 0.8284  0.655  
    ## + diet_OCarb_g          1 -3.6397 0.6740  0.775  
    ## + diet_Vit_B6_mg        1 -3.7278 0.7398  0.805  
    ## + diet_SatFat_g         1 -3.6388 0.6733  0.840  
    ## + age                   1 -3.5465 0.6049  0.880  
    ## + cohort                1 -3.1615 0.3246  0.970  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ordibroccA
    ## Call: capscale(formula = MB_asvtab_broccA ~ diet_Alc_g + bmi, data =
    ## meta_groupA_brocc, distance = "bray")
    ## 
    ## -- Model Summary --
    ## 
    ##               Inertia Proportion Rank
    ## Total          0.6843     1.0000     
    ## Constrained    0.2086     0.3049    2
    ## Unconstrained  0.4756     0.6951   10
    ## 
    ## Inertia is squared Bray distance
    ## 
    ## -- Note --
    ## 
    ## Species scores projected from 'MB_asvtab_broccA'
    ## 
    ## -- Eigenvalues --
    ## 
    ## Eigenvalues for constrained axes:
    ##    CAP1    CAP2 
    ## 0.12166 0.08699 
    ## 
    ## Eigenvalues for unconstrained axes:
    ##    MDS1    MDS2    MDS3    MDS4    MDS5    MDS6    MDS7    MDS8    MDS9   MDS10 
    ## 0.11755 0.08670 0.06758 0.06391 0.04911 0.03375 0.02661 0.01683 0.01066 0.00293
    # best model is MB_asvtab_broccA ~ bmi + diet_Alc_g
    # model inertia was 0.6843 (low, not great); Inertia is squared Bray distance

    # calculate RDA (distance based)
    rda_broccA <- capscale(formula = MB_asvtab_broccA ~ bmi + diet_Alc_g, 
                          data = meta_groupA_brocc, 
                          distance = "bray")
    RsquareAdj(rda_broccA)
    ## $r.squared
    ## [1] 0.3049153
    ## 
    ## $adj.r.squared
    ## [1] 0.1658983
    # 0.17 of variation explained by this model
    summary(rda_broccA)
    ## 
    ## Call:
    ## capscale(formula = MB_asvtab_broccA ~ bmi + diet_Alc_g, data = meta_groupA_brocc,      distance = "bray") 
    ## 
    ## Partitioning of squared Bray distance:
    ##               Inertia Proportion
    ## Total          0.6843     1.0000
    ## Constrained    0.2086     0.3049
    ## Unconstrained  0.4756     0.6951
    ## 
    ## Eigenvalues, and their contribution to the squared Bray distance 
    ## 
    ## Importance of components:
    ##                         CAP1    CAP2   MDS1   MDS2    MDS3    MDS4    MDS5
    ## Eigenvalue            0.1217 0.08699 0.1176 0.0867 0.06758 0.06391 0.04911
    ## Proportion Explained  0.1778 0.12712 0.1718 0.1267 0.09877 0.09340 0.07177
    ## Cumulative Proportion 0.1778 0.30492 0.4767 0.6034 0.70217 0.79557 0.86734
    ##                          MDS6    MDS7    MDS8    MDS9    MDS10
    ## Eigenvalue            0.03375 0.02661 0.01683 0.01066 0.002934
    ## Proportion Explained  0.04932 0.03888 0.02460 0.01557 0.004288
    ## Cumulative Proportion 0.91666 0.95554 0.98014 0.99571 1.000000
    ## 
    ## Accumulated constrained eigenvalues
    ## Importance of components:
    ##                         CAP1    CAP2
    ## Eigenvalue            0.1217 0.08699
    ## Proportion Explained  0.5831 0.41691
    ## Cumulative Proportion 0.5831 1.00000

    plot(rda_broccA)

<img src="broccolipreMBpermanova_files/figure-markdown_strict/broccoli dbRDA run with best model from ordistep group A-1.png" width="98%" height="98%" />

    Perm_broccA <- adonis2(MB_asvtab_broccA ~ bmi + diet_Alc_g, 
                                data = meta_groupA_brocc,
                           method = "bray",
                           by = "terms")
    Perm_broccA
    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = MB_asvtab_broccA ~ bmi + diet_Alc_g, data = meta_groupA_brocc, method = "bray", by = "terms")
    ##            Df SumOfSqs      R2      F Pr(>F)   
    ## bmi         1  0.10477 0.15312 2.2028  0.008 **
    ## diet_Alc_g  1  0.10387 0.15180 2.1839  0.005 **
    ## Residual   10  0.47564 0.69508                 
    ## Total      12  0.68428 1.00000                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    p.adjust(Perm_broccA$`Pr(>F)`, method = "fdr")
    ## [1] 0.008 0.008    NA    NA

    # write.table(Perm_broccA,
    #             "../outputs/Perm_broccA.txt",
    #             sep = "\t")

# group B

    meta_groupB_brocc <- meta_groups_brocc %>%
      dplyr::filter(group == "B") %>%
      dplyr::select(-c(group, group_numb, ethnicity, relation, height_in, weight_lb)) %>%
      column_to_rownames("subject_id")
    MB_asvtab_broccB <- MB_asvtab_brocc %>%
      rownames_to_column("subject_id") %>%
      dplyr::filter(subject_id %in% rownames(meta_groupB_brocc)) %>%
      column_to_rownames("subject_id")
    Bbroccmod0 <- capscale(MB_asvtab_broccB ~ 1, meta_groupB_brocc, distance = "bray")  # Model with intercept only
    Bbroccmod1 <- capscale(MB_asvtab_broccB ~ . + .*., meta_groupB_brocc, distance = "bray")  # Model with all explanatory variables
    ordibroccB <- ordistep(Bbroccmod0, scope = formula(Bbroccmod1)) # this determines what the best model is to run RDA on
    ## 
    ## Start: MB_asvtab_broccB ~ 1 
    ## 
    ##                        Df    AIC      F Pr(>F)   
    ## + race                  3 13.197 1.9541  0.010 **
    ## + diet_Copp_mg          1 13.324 2.0277  0.020 * 
    ## + diet_SatFat_g         1 13.401 1.9509  0.040 * 
    ## + diet_condensed        4 14.685 1.5006  0.060 . 
    ## + diet_TotInsolFib_g    1 13.668 1.6856  0.110   
    ## + diet_TotSolFib_g      1 13.903 1.4568  0.120   
    ## + diet_Vit_B6_mg        1 13.899 1.4605  0.130   
    ## + diet_OCarb_g          1 14.033 1.3307  0.145   
    ## + sex                   1 14.132 1.2355  0.170   
    ## + diet_Chol_mg          1 14.091 1.2742  0.210   
    ## + diet_MonSac_g         1 14.260 1.1133  0.275   
    ## + diet_Iron_mg          1 14.318 1.0578  0.380   
    ## + diet_TotFib_g         1 14.445 0.9376  0.495   
    ## + diet_Alc_g            1 14.533 0.8546  0.530   
    ## + age                   1 14.507 0.8792  0.580   
    ## + bmi                   1 14.510 0.8761  0.580   
    ## + diet_BetaCaro_mcg     1 14.576 0.8141  0.595   
    ## + diet_Vit_B1_mg        1 14.589 0.8023  0.695   
    ## + diet_Disacc_g         1 14.648 0.7472  0.710   
    ## + diet_Mang_mg          1 14.627 0.7665  0.720   
    ## + cohort                1 14.652 0.7430  0.740   
    ## + diet_Fol.DFE_mcg_DFE  1 14.784 0.6207  0.885   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Step: MB_asvtab_broccB ~ race 
    ## 
    ##        Df    AIC      F Pr(>F)  
    ## - race  3 13.465 1.9541  0.015 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##                        Df    AIC      F Pr(>F)  
    ## + diet_Vit_B6_mg        1 12.595 2.0547  0.025 *
    ## + diet_SatFat_g         1 12.712 1.9562  0.035 *
    ## + diet_TotInsolFib_g    1 12.775 1.9036  0.060 .
    ## + diet_condensed        3 13.615 1.3660  0.125  
    ## + diet_TotSolFib_g      1 13.313 1.4596  0.140  
    ## + diet_OCarb_g          1 13.418 1.3738  0.150  
    ## + bmi                   1 13.531 1.2830  0.190  
    ## + diet_Vit_B1_mg        1 13.469 1.3328  0.205  
    ## + diet_MonSac_g         1 13.502 1.3064  0.210  
    ## + diet_Chol_mg          1 13.579 1.2441  0.230  
    ## + sex                   1 13.537 1.2779  0.280  
    ## + diet_Mang_mg          1 13.781 1.0832  0.310  
    ## + diet_Copp_mg          1 13.870 1.0128  0.355  
    ## + diet_Iron_mg          1 13.758 1.1018  0.425  
    ## + diet_TotFib_g         1 13.933 0.9633  0.430  
    ## + diet_Disacc_g         1 13.762 1.0982  0.435  
    ## + diet_Fol.DFE_mcg_DFE  1 13.989 0.9191  0.480  
    ## + diet_BetaCaro_mcg     1 14.120 0.8166  0.710  
    ## + age                   1 14.222 0.7373  0.730  
    ## + cohort                1 14.200 0.7546  0.735  
    ## + diet_Alc_g            1 14.660 0.4012  0.975  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Step: MB_asvtab_broccB ~ race + diet_Vit_B6_mg 
    ## 
    ##                  Df    AIC      F Pr(>F)   
    ## - diet_Vit_B6_mg  1 13.197 2.0547  0.010 **
    ## - race            3 13.899 2.1875  0.005 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##                        Df    AIC      F Pr(>F)  
    ## + diet_Vit_B1_mg        1 11.809 2.0532  0.030 *
    ## + diet_SatFat_g         1 12.278 1.6860  0.050 *
    ## + diet_TotInsolFib_g    1 12.448 1.5549  0.085 .
    ## + diet_condensed        3 12.180 1.4727  0.115  
    ## + bmi                   1 12.543 1.4824  0.120  
    ## + race:diet_Vit_B6_mg   1 12.286 1.6799  0.125  
    ## + diet_OCarb_g          1 12.605 1.4357  0.125  
    ## + diet_MonSac_g         1 12.734 1.3378  0.190  
    ## + diet_Chol_mg          1 12.878 1.2293  0.255  
    ## + diet_Disacc_g         1 12.910 1.2054  0.255  
    ## + diet_TotSolFib_g      1 12.771 1.3098  0.275  
    ## + diet_Iron_mg          1 12.955 1.1717  0.285  
    ## + sex                   1 12.965 1.1645  0.295  
    ## + diet_Mang_mg          1 12.957 1.1704  0.305  
    ## + diet_Copp_mg          1 13.069 1.0870  0.390  
    ## + diet_BetaCaro_mcg     1 13.287 0.9264  0.530  
    ## + diet_TotFib_g         1 13.340 0.8878  0.555  
    ## + age                   1 13.447 0.8100  0.660  
    ## + diet_Fol.DFE_mcg_DFE  1 13.550 0.7351  0.690  
    ## + cohort                1 13.901 0.4836  0.960  
    ## + diet_Alc_g            1 14.082 0.3561  0.990  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Step: MB_asvtab_broccB ~ race + diet_Vit_B6_mg + diet_Vit_B1_mg 
    ## 
    ##                  Df    AIC      F Pr(>F)   
    ## - diet_Vit_B1_mg  1 12.595 2.0532  0.020 * 
    ## - diet_Vit_B6_mg  1 13.469 2.7619  0.005 **
    ## - race            3 14.149 2.3879  0.005 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##                                 Df    AIC      F Pr(>F)  
    ## + race:diet_Vit_B6_mg            1 11.027 1.8923  0.035 *
    ## + race:diet_Vit_B1_mg            1 11.027 1.8923  0.045 *
    ## + diet_TotInsolFib_g             1 11.249 1.7308  0.090 .
    ## + diet_Copp_mg                   1 11.607 1.4746  0.140  
    ## + diet_Chol_mg                   1 11.747 1.3755  0.175  
    ## + diet_OCarb_g                   1 11.762 1.3652  0.180  
    ## + diet_condensed                 3 11.793 1.2416  0.205  
    ## + bmi                            1 11.840 1.3100  0.215  
    ## + sex                            1 11.867 1.2911  0.225  
    ## + diet_MonSac_g                  1 11.886 1.2783  0.235  
    ## + diet_TotSolFib_g               1 11.858 1.2976  0.250  
    ## + diet_Mang_mg                   1 12.055 1.1605  0.305  
    ## + diet_SatFat_g                  1 12.020 1.1850  0.320  
    ## + diet_Vit_B1_mg:diet_Vit_B6_mg  1 12.182 1.0729  0.405  
    ## + diet_BetaCaro_mcg              1 12.263 1.0171  0.480  
    ## + diet_TotFib_g                  1 12.335 0.9679  0.525  
    ## + diet_Disacc_g                  1 12.358 0.9525  0.540  
    ## + age                            1 12.636 0.7639  0.640  
    ## + diet_Iron_mg                   1 12.530 0.8357  0.670  
    ## + diet_Fol.DFE_mcg_DFE           1 12.631 0.7677  0.725  
    ## + cohort                         1 13.007 0.5173  0.950  
    ## + diet_Alc_g                     1 13.141 0.4290  0.950  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Step: MB_asvtab_broccB ~ race + diet_Vit_B6_mg + diet_Vit_B1_mg + race:diet_Vit_B6_mg 
    ## 
    ##                       Df    AIC      F Pr(>F)   
    ## - race:diet_Vit_B6_mg  1 11.809 1.8923   0.06 . 
    ## - diet_Vit_B1_mg       1 12.286 2.2455   0.01 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##                                 Df    AIC      F Pr(>F)
    ## + diet_TotInsolFib_g             1 10.604 1.4961  0.130
    ## + diet_Chol_mg                   1 10.549 1.5322  0.145
    ## + bmi                            1 10.791 1.3738  0.190
    ## + diet_SatFat_g                  1 10.963 1.2618  0.210
    ## + diet_Copp_mg                   1 11.187 1.1184  0.270
    ## + diet_OCarb_g                   1 11.124 1.1585  0.280
    ## + diet_Mang_mg                   1 11.183 1.1207  0.300
    ## + diet_condensed                 3 10.570 1.2141  0.310
    ## + diet_TotSolFib_g               1 10.805 1.3646  0.315
    ## + diet_BetaCaro_mcg              1 11.176 1.1252  0.335
    ## + diet_Vit_B1_mg:diet_Vit_B6_mg  1 11.224 1.0950  0.430
    ## + diet_TotFib_g                  1 11.377 0.9982  0.445
    ## + sex                            1 11.488 0.9281  0.525
    ## + diet_Fol.DFE_mcg_DFE           1 11.622 0.8442  0.565
    ## + diet_Iron_mg                   1 11.613 0.8495  0.585
    ## + age                            1 11.652 0.8252  0.660
    ## + diet_Disacc_g                  1 11.771 0.7515  0.735
    ## + diet_MonSac_g                  1 12.040 0.5861  0.805
    ## + cohort                         1 12.081 0.5616  0.910
    ## + diet_Alc_g                     1 12.293 0.4332  0.985
    ## + race:diet_Vit_B1_mg            0 11.027
    ordibroccB
    ## Call: capscale(formula = MB_asvtab_broccB ~ race + diet_Vit_B6_mg +
    ## diet_Vit_B1_mg + race:diet_Vit_B6_mg, data = meta_groupB_brocc, distance =
    ## "bray")
    ## 
    ## -- Model Summary --
    ## 
    ##                Inertia Proportion Rank
    ## Total          1.91939                
    ## RealTotal      1.92992    1.00000     
    ## Constrained    1.02727    0.53229    6
    ## Unconstrained  0.90265    0.46771   12
    ## Imaginary     -0.01053                
    ## 
    ## Inertia is squared Bray distance
    ## 
    ## -- Note --
    ## 
    ## Species scores projected from 'MB_asvtab_broccB'
    ## 
    ## Some constraints or conditions were aliased because they were redundant.
    ## This can happen if terms are linearly dependent (collinear):
    ## 'raceOther:diet_Vit_B6_mg', 'raceWhite:diet_Vit_B6_mg'
    ## 
    ## -- Eigenvalues --
    ## 
    ## Eigenvalues for constrained axes:
    ##   CAP1   CAP2   CAP3   CAP4   CAP5   CAP6 
    ## 0.4002 0.2540 0.1431 0.1273 0.0611 0.0415 
    ## 
    ## Eigenvalues for unconstrained axes:
    ##    MDS1    MDS2    MDS3    MDS4    MDS5    MDS6    MDS7    MDS8    MDS9   MDS10 
    ## 0.19994 0.17107 0.13066 0.10114 0.08478 0.06906 0.04064 0.03603 0.02845 0.02354 
    ##   MDS11   MDS12 
    ## 0.01254 0.00480
    # best model is MB_asvtab_broccB ~ race + diet_Vit_B6_mg + diet_Vit_B1_mg +
    # race:diet_Vit_B1_mg
    # model inertia was 1.91939; Inertia is squared Bray distance

    # calculate RDA (distance based)
    rda_broccB <- capscale(formula = MB_asvtab_broccB ~ race + diet_Vit_B6_mg + diet_Vit_B1_mg +
    race:diet_Vit_B1_mg, 
                          data = meta_groupB_brocc, 
                          distance = "bray")
    RsquareAdj(rda_broccB)
    ## $r.squared
    ## [1] 0.5322856
    ## 
    ## $adj.r.squared
    ## [1] 0.2984284
    # 0.3 of variation explained by this model
    # race other and race white are collinear with vitamin B1
    summary(rda_broccB)
    ## 
    ## Call:
    ## capscale(formula = MB_asvtab_broccB ~ race + diet_Vit_B6_mg +      diet_Vit_B1_mg + race:diet_Vit_B1_mg, data = meta_groupB_brocc,      distance = "bray") 
    ## 
    ## Partitioning of squared Bray distance:
    ##               Inertia Proportion
    ## Total          1.9299     1.0000
    ## Constrained    1.0273     0.5323
    ## Unconstrained  0.9027     0.4677
    ## 
    ## Eigenvalues, and their contribution to the squared Bray distance 
    ## 
    ## Importance of components:
    ##                         CAP1   CAP2    CAP3    CAP4    CAP5    CAP6   MDS1
    ## Eigenvalue            0.4002 0.2540 0.14314 0.12729 0.06111 0.04155 0.1999
    ## Proportion Explained  0.2074 0.1316 0.07417 0.06596 0.03166 0.02153 0.1036
    ## Cumulative Proportion 0.2074 0.3390 0.41314 0.47909 0.51076 0.53229 0.6359
    ##                          MDS2   MDS3    MDS4    MDS5    MDS6    MDS7    MDS8
    ## Eigenvalue            0.17107 0.1307 0.10114 0.08478 0.06906 0.04064 0.03603
    ## Proportion Explained  0.08864 0.0677 0.05241 0.04393 0.03578 0.02106 0.01867
    ## Cumulative Proportion 0.72453 0.7922 0.84464 0.88857 0.92435 0.94541 0.96408
    ##                          MDS9   MDS10    MDS11    MDS12
    ## Eigenvalue            0.02845 0.02354 0.012540 0.004795
    ## Proportion Explained  0.01474 0.01220 0.006497 0.002485
    ## Cumulative Proportion 0.97882 0.99102 0.997515 1.000000
    ## 
    ## Accumulated constrained eigenvalues
    ## Importance of components:
    ##                         CAP1   CAP2   CAP3   CAP4    CAP5    CAP6
    ## Eigenvalue            0.4002 0.2540 0.1431 0.1273 0.06111 0.04155
    ## Proportion Explained  0.3896 0.2473 0.1393 0.1239 0.05948 0.04045
    ## Cumulative Proportion 0.3896 0.6368 0.7762 0.9001 0.95955 1.00000

    plot(rda_broccB)

<img src="broccolipreMBpermanova_files/figure-markdown_strict/broccoli dbRDA run with best model from ordistep group B-1.png" width="98%" height="98%" />

    Perm_broccB <- adonis2(MB_asvtab_broccB ~ race + diet_Vit_B6_mg + diet_Vit_B1_mg +
    race:diet_Vit_B1_mg, 
                                data = meta_groupB_brocc,
                           method = "bray",
                           by = "terms")
    Perm_broccB
    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = MB_asvtab_broccB ~ race + diet_Vit_B6_mg + diet_Vit_B1_mg + race:diet_Vit_B1_mg, data = meta_groupB_brocc, method = "bray", by = "terms")
    ##                     Df SumOfSqs      R2      F Pr(>F)   
    ## race                 3  0.54190 0.28233 2.4271  0.003 **
    ## diet_Vit_B6_mg       1  0.17758 0.09252 2.3860  0.008 **
    ## diet_Vit_B1_mg       1  0.16450 0.08570 2.2103  0.014 * 
    ## race:diet_Vit_B1_mg  1  0.14232 0.07415 1.9123  0.044 * 
    ## Residual            12  0.89309 0.46530                 
    ## Total               18  1.91939 1.00000                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    p.adjust(Perm_broccB$`Pr(>F)`, method = "fdr")
    ## [1] 0.01200000 0.01600000 0.01866667 0.04400000         NA         NA

    # write.table(Perm_broccB,
    #             "../outputs/Perm_broccB.txt",
    #             sep = "\t")

# group C

    meta_groupC_brocc <- meta_groups_brocc %>%
      dplyr::filter(group == "C") %>%
      dplyr::select(-c(group, group_numb, ethnicity, relation, height_in, weight_lb)) %>%
      column_to_rownames("subject_id")
    MB_asvtab_groupC <- MB_asvtab_brocc %>%
      rownames_to_column("subject_id") %>%
      dplyr::filter(subject_id %in% rownames(meta_groupC_brocc)) %>%
      column_to_rownames("subject_id")
    Cbroccmod0 <- capscale(MB_asvtab_groupC ~ 1, meta_groupC_brocc, distance = "bray")  # Model with intercept only
    Cbroccmod1 <- capscale(MB_asvtab_groupC ~ . + .*., meta_groupC_brocc, distance = "bray")  # Model with all explanatory variables
    ordigroupC <- ordistep(Cbroccmod0, scope = formula(Cbroccmod1)) # this determines what the best model is to run RDA on
    ## 
    ## Start: MB_asvtab_groupC ~ 1
    ## 
    ##                        Df     AIC      F Pr(>F)   
    ## + diet_Fol.DFE_mcg_DFE  1 -5.1774 2.2007  0.010 **
    ## + diet_MonSac_g         1 -4.7534 1.7776  0.045 * 
    ## + diet_BetaCaro_mcg     1 -4.7166 1.7422  0.095 . 
    ## + cohort                1 -4.4980 1.5368  0.160   
    ## + race                  2 -4.6710 1.4825  0.185   
    ## + bmi                   1 -4.2413 1.3049  0.305   
    ## + diet_Copp_mg          1 -4.2004 1.2689  0.340   
    ## + sex                   1 -4.1034 1.1844  0.390   
    ## + age                   1 -3.9559 1.0585  0.435   
    ## + diet_SatFat_g         1 -3.8529 0.9724  0.460   
    ## + diet_TotSolFib_g      1 -3.7570 0.8935  0.490   
    ## + diet_TotFib_g         1 -3.8042 0.9322  0.545   
    ## + diet_TotInsolFib_g    1 -3.7463 0.8849  0.565   
    ## + diet_OCarb_g          1 -3.7780 0.9107  0.570   
    ## + diet_condensed        1 -3.4857 0.6772  0.730   
    ## + diet_Alc_g            1 -3.5330 0.7142  0.755   
    ## + diet_Vit_B6_mg        1 -3.4949 0.6844  0.760   
    ## + diet_Chol_mg          1 -3.4188 0.6254  0.775   
    ## + diet_Vit_B1_mg        1 -3.3879 0.6016  0.830   
    ## + diet_Mang_mg          1 -3.2263 0.4793  0.890   
    ## + diet_Disacc_g         1 -3.0820 0.3729  0.970   
    ## + diet_Iron_mg          1 -3.1404 0.4156  0.980   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Step: MB_asvtab_groupC ~ diet_Fol.DFE_mcg_DFE
    ## 
    ##                        Df     AIC      F Pr(>F)  
    ## - diet_Fol.DFE_mcg_DFE  1 -4.5473 2.2007  0.035 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##                      Df     AIC      F Pr(>F)
    ## + bmi                 1 -5.8353 1.6720  0.140
    ## + age                 1 -6.0147 1.8138  0.180
    ## + diet_MonSac_g       1 -5.4931 1.4130  0.225
    ## + cohort              1 -5.3708 1.3240  0.245
    ## + diet_TotSolFib_g    1 -5.4742 1.3991  0.265
    ## + sex                 1 -5.4352 1.3706  0.305
    ## + diet_TotInsolFib_g  1 -5.4707 1.3966  0.315
    ## + diet_BetaCaro_mcg   1 -5.0009 1.0654  0.360
    ## + diet_SatFat_g       1 -5.2489 1.2370  0.395
    ## + diet_Copp_mg        1 -4.9765 1.0490  0.505
    ## + diet_Vit_B6_mg      1 -4.5762 0.7876  0.540
    ## + race                2 -5.3694 1.0110  0.545
    ## + diet_condensed      1 -4.6160 0.8128  0.570
    ## + diet_Mang_mg        1 -4.4988 0.7391  0.630
    ## + diet_Alc_g          1 -4.4269 0.6946  0.640
    ## + diet_OCarb_g        1 -4.5522 0.7725  0.655
    ## + diet_Iron_mg        1 -4.3008 0.6177  0.745
    ## + diet_Chol_mg        1 -4.2775 0.6037  0.785
    ## + diet_Vit_B1_mg      1 -4.1620 0.5349  0.800
    ## + diet_Disacc_g       1 -4.1665 0.5376  0.840
    ## + diet_TotFib_g       1 -3.9338 0.4030  0.890
    ordigroupC
    ## Call: capscale(formula = MB_asvtab_groupC ~ diet_Fol.DFE_mcg_DFE, data =
    ## meta_groupC_brocc, distance = "bray")
    ## 
    ## -- Model Summary --
    ## 
    ##               Inertia Proportion Rank
    ## Total          0.4030     1.0000     
    ## Constrained    0.1430     0.3549    1
    ## Unconstrained  0.2600     0.6451    4
    ## 
    ## Inertia is squared Bray distance
    ## 
    ## -- Note --
    ## 
    ## Species scores projected from 'MB_asvtab_groupC'
    ## 
    ## -- Eigenvalues --
    ## 
    ## Eigenvalues for constrained axes:
    ##    CAP1 
    ## 0.14302 
    ## 
    ## Eigenvalues for unconstrained axes:
    ##    MDS1    MDS2    MDS3    MDS4 
    ## 0.12274 0.07196 0.04059 0.02467
    # best model is MB_asvtab_groupC ~ diet_Fol.DFE_mcg_DFE
    # model inertia was 0.4030 (low, not great); Inertia is squared Bray distance

    # calculate RDA (distance based)
    rda_groupC <- capscale(formula = MB_asvtab_groupC ~ diet_Fol.DFE_mcg_DFE, 
                          data = meta_groupC_brocc, 
                          distance = "bray")
    RsquareAdj(rda_groupC)
    ## $r.squared
    ## [1] 0.3549098
    ## 
    ## $adj.r.squared
    ## [1] 0.1936373
    # 0.19 of variation explained by this model
    summary(rda_groupC)
    ## 
    ## Call:
    ## capscale(formula = MB_asvtab_groupC ~ diet_Fol.DFE_mcg_DFE, data = meta_groupC_brocc,      distance = "bray") 
    ## 
    ## Partitioning of squared Bray distance:
    ##               Inertia Proportion
    ## Total           0.403     1.0000
    ## Constrained     0.143     0.3549
    ## Unconstrained   0.260     0.6451
    ## 
    ## Eigenvalues, and their contribution to the squared Bray distance 
    ## 
    ## Importance of components:
    ##                         CAP1   MDS1    MDS2    MDS3    MDS4
    ## Eigenvalue            0.1430 0.1227 0.07196 0.04059 0.02467
    ## Proportion Explained  0.3549 0.3046 0.17858 0.10072 0.06122
    ## Cumulative Proportion 0.3549 0.6595 0.83806 0.93878 1.00000
    ## 
    ## Accumulated constrained eigenvalues
    ## Importance of components:
    ##                        CAP1
    ## Eigenvalue            0.143
    ## Proportion Explained  1.000
    ## Cumulative Proportion 1.000

    plot(rda_groupC)

<img src="broccolipreMBpermanova_files/figure-markdown_strict/broccoli dbRDA run with best model from ordistep group C-1.png" width="98%" height="98%" />

    Perm_groupC <- adonis2(MB_asvtab_groupC ~ diet_Fol.DFE_mcg_DFE, 
                                data = meta_groupC_brocc,
                           method = "bray",
                           by = "terms")
    Perm_groupC
    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 719
    ## 
    ## adonis2(formula = MB_asvtab_groupC ~ diet_Fol.DFE_mcg_DFE, data = meta_groupC_brocc, method = "bray", by = "terms")
    ##                      Df SumOfSqs      R2      F  Pr(>F)  
    ## diet_Fol.DFE_mcg_DFE  1  0.14302 0.35491 2.2007 0.02083 *
    ## Residual              4  0.25995 0.64509                 
    ## Total                 5  0.40297 1.00000                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    p.adjust(Perm_groupC$`Pr(>F)`, method = "fdr")
    ## [1] 0.02083333         NA         NA

    # write.table(Perm_groupC,
    #             "../outputs/Perm_groupC.txt",
    #             sep = "\t")

## save image chunk

    save.image("/Users/alexieva/Documents/Projects/Analysis/broccoli project/02_pubAnalysis/Enterotyping/broccolipreMBpermanova.RData")
