############# r loop through all dietary components with each SFN ############# 

list_sfn <- c("SFN", "SFN_NIT")
# and use diet_compnts as the list of dietary components
# make an output file
model.0_lms <- c()
pvals_model.0_lms <- c()

for (a in 1:length(list_sfn)) {
  for (e in 1:length(diet_compnts)) {
    model.0 <- lm(paste0("log(", list_sfn[a], "+1) ~ ", 
                         diet_compnts[e], "+ I(", diet_compnts[e], "^2) + I(", 
                         diet_compnts[e], "^3) + I(", diet_compnts[e], "^4)"), 
                  data = med_input)
    modsum <- summary(model.0)
    coeff_tab <- as.data.frame(modsum$coefficients)
    model.0_lms[[paste(list_sfn[a], "~", diet_compnts[e])]] <- coeff_tab
    pvals_model.0_lms[[paste(list_sfn[a], "~", diet_compnts[e])]][["diet"]] <- coeff_tab[2,4]
    pvals_model.0_lms[[paste(list_sfn[a], "~", diet_compnts[e])]][["diet2"]] <- coeff_tab[3,4]
    pvals_model.0_lms[[paste(list_sfn[a], "~", diet_compnts[e])]][["diet3"]] <- coeff_tab[4,4]
    pvals_model.0_lms[[paste(list_sfn[a], "~", diet_compnts[e])]][["diet4"]] <- coeff_tab[5,4]
    e <- e + 1
  }
  a <- a + 1
}


############# r loop through all SFNs, dietary components, and ASVs ############# 
# Use list_sfn and diet_compnts as the list of dietary components
# make SFN taxa a list
taxa_list <- as.character(SFN_taxa$x)
# make an output file
med_res <- c()
pvals_med <- c()
exp_res_med <- c()
list_combos <- data.frame(sfns = c("SFN", "SFN", "SFN_NIT"),
                          diet = c("diet_Vit_B1_mg", "diet_TotFib_g", "diet_Iron_mg"))

for (a in 1:nrow(list_combos)) {
  for (f in 1:length(taxa_list)) {
    # add quartiles for each dietary component that will be used
    # transform ASV abundance and transform SFN abundance
    med_input <- med_input %>%
      dplyr::mutate(quarts = ntile(med_input[[list_combos[a,2]]], 4),
                    sfn_trans = log(med_input[[list_combos[a,1]]] + 1),
                    asv_trans = log(med_input[[taxa_list[f]]] + 1))
    
    # set models
    model.M <- lm(asv_trans ~ quarts,
                  med_input)
    model.Y <- lm(sfn_trans ~ quarts + asv_trans,
                  med_input)
    
    # mediate
    results <- mediate(model.M, model.Y, 
                       treat = "quarts", 
                       mediator = "asv_trans",
                       control.value = 1, treat.value = 4, # lowest and highest quartiles
                       boot = TRUE, sims = 500)
    
    # save results
    exp_res_med[[paste(list_combos[a,1], taxa_list[f], list_combos[a,2])]] <- results
    modsum <- summary(results)
    pvals_med[[paste(list_combos[a,1], taxa_list[f], list_combos[a,2])]] <- list(ACME_pval = modsum$d.avg.p,
                                                                                 ADE_pval = modsum$z.avg.p,
                                                                                 prop_mediated_pval = modsum$n.avg.p)
    
    if (modsum$d.avg.p <= 0.2) {
      coeff_tab <- list(ACME_est = modsum$d.avg,
                        ACME_pval = modsum$d.avg.p,
                        ADE_est = modsum$z.avg,
                        ADE_pval = modsum$z.avg.p,
                        prop_mediated_est = modsum$n.avg,
                        prop_mediated_pval = modsum$n.avg.p)
      
      med_res[[paste(list_combos[a,1], taxa_list[f], list_combos[a,2])]] <- coeff_tab
      
    } else { print(paste(list_combos[a,1], ",", taxa_list[f], ",",
                         list_combos[a,2], "are NOT significant"))}
    
    # counters
    f <- f + 1
  }
  a <- a + 1
}


