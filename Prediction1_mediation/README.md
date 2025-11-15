## libraries and directories

    library(dplyr)
    library(tibble)
    library(tidyr)
    library(ggplot2)
    library(mediation)
    library(ggfortify)
    library(ggpubr)
    library(tweedie)
    library(statmod)
    library(AICcmodavg)
    library(stringr)
    set.seed(3)

    setwd("/Users/alexieva/Documents/Projects/Analysis/broccoli project/02_pubAnalysis/Prediction1_mediation/")

    load("/Users/alexieva/Documents/Projects/Analysis/broccoli project/02_pubAnalysis/Prediction1_mediation/brocc_mediation.RData")

## Input data

    ## microbiome asv table; this is at time 0 as per previous exploratory data analysis
    # only taxa that were previously identified as correlated with SFN and SFN-NIT, whose taxa names are in this SFN_taxa vector
    SFN_taxa <- read.csv("/Users/alexieva/Documents/Projects/Analysis/broccoli project/01_explAnalysis/input_data/SFN_taxa.txt", header = T)

    MB_asvtab <- read.table("/Users/alexieva/Documents/Projects/Analysis/broccoli project/01_explAnalysis/input_data/MBasvtab_brocc_cleanT0.txt", header = TRUE) %>%
      dplyr::select(SFN_taxa$x) %>%
      rownames_to_column("subject_id")

    ## SFN data
    # using cumulative SFN, SFN-NIT, and SFN tot across all time points
    # note that there are a few alfalfa eaters here that get filtered out in the process of making combined data frame of med_input
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
      # %>% spread(key = "SFN_metabolites", value = "quantity")

    ## diet data
    # make a list of the "significant" diet factors we are interested in 
    diet_compnts <- c("diet_Chol_mg", "diet_SatFat_g", "diet_Alc_g", "diet_TotSolFib_g", 
                      "diet_MonSac_g", "diet_Disacc_g", "diet_BetaCaro_mcg", "diet_Vit_B1_mg",
                      "diet_Vit_B6_mg", "diet_Fol.DFE_mcg_DFE", "diet_Mang_mg", "diet_Copp_mg", 
                      "diet_TotFib_g", "diet_OCarb_g", "diet_Iron_mg", "diet_TotInsolFib_g")
    diet <- read.csv("/Users/alexieva/Documents/Projects/Analysis/broccoli project/00_data/Diet/New_Diet_Small.csv", header = TRUE) %>%
      dplyr::select(c("subject_id", all_of(diet_compnts)))

    # write.table(diet, "../input_data/dietData_clean.txt", row.names = F, quote = F)

    # match data based on subject id
    med_input <- MB_asvtab %>%
      inner_join(Urine_SFN, by = "subject_id") %>%
      inner_join(diet, by = "subject_id")

    # write.table(med_input, "../input_data/combinedData_clean.txt", row.names = F, quote = F)

## LEFT OFF: rerunning this below with insol fiber data frames

## check distributions of variables

    plot(density(MB_asvtab$Bifidobacterium))

<img src="brocc_mediation_files/figure-markdown_strict/check distrubutions of a couple ASVs we are interested in-1.png" width="98%" height="98%" />


    plot(density(MB_asvtab$Lachnospiraceae.NK4A136.group))

<img src="brocc_mediation_files/figure-markdown_strict/check distrubutions of a couple ASVs we are interested in-2.png" width="98%" height="98%" />


    plot(density(MB_asvtab$Ruminococcus))

<img src="brocc_mediation_files/figure-markdown_strict/check distrubutions of a couple ASVs we are interested in-3.png" width="98%" height="98%" />


    plot(density(MB_asvtab$Roseburia))

<img src="brocc_mediation_files/figure-markdown_strict/check distrubutions of a couple ASVs we are interested in-4.png" width="98%" height="98%" />


    plot(density(MB_asvtab$Blautia))

<img src="brocc_mediation_files/figure-markdown_strict/check distrubutions of a couple ASVs we are interested in-5.png" width="98%" height="98%" />

    d <- density(MB_asvtab$Blautia)
    quantile(MB_asvtab$Blautia)
    ##      0%     25%     50%     75%    100% 
    ##  763.00 1680.50 2099.50 3017.25 5206.00
    d$x[which.max(d$y)] # use the max peak of the curve for treat.value in mediate() later
    ## [1] 1929.751
    # this value is 1929.751

    # most have right skewed density distributions with long right side tails, as is usual for ASV abundances, Blautia is at least approximately normal

    plot(density(Urine_SFN$sum_SFN)) # kinda right skewed

<img src="brocc_mediation_files/figure-markdown_strict/check distributions of SFN and SFN-NIT, just to get an idea-1.png" width="98%" height="98%" />

    plot(density(Urine_SFN$sum_SFN_NIT)) # this one is more normal but has a right tail

<img src="brocc_mediation_files/figure-markdown_strict/check distributions of SFN and SFN-NIT, just to get an idea-2.png" width="98%" height="98%" />

    plot(density(Urine_SFN$sum_SFN_Tot)) # kinda right skewed

<img src="brocc_mediation_files/figure-markdown_strict/check distributions of SFN and SFN-NIT, just to get an idea-3.png" width="98%" height="98%" />


    hist(Urine_SFN$sum_SFN)

<img src="brocc_mediation_files/figure-markdown_strict/check distributions of SFN and SFN-NIT, just to get an idea-4.png" width="98%" height="98%" />

    hist(Urine_SFN$sum_SFN_NIT)

<img src="brocc_mediation_files/figure-markdown_strict/check distributions of SFN and SFN-NIT, just to get an idea-5.png" width="98%" height="98%" />

    hist(Urine_SFN$sum_SFN_Tot)

<img src="brocc_mediation_files/figure-markdown_strict/check distributions of SFN and SFN-NIT, just to get an idea-6.png" width="98%" height="98%" />

    plot(density(diet$diet_TotSolFib_g)) # very right skewed

<img src="brocc_mediation_files/figure-markdown_strict/check distributions of total and soluble fiber-1.png" width="98%" height="98%" />

    plot(density(diet$diet_TotFib_g)) # more normal

<img src="brocc_mediation_files/figure-markdown_strict/check distributions of total and soluble fiber-2.png" width="98%" height="98%" />

    plot(density(diet$diet_TotInsolFib_g)) # right skewed

<img src="brocc_mediation_files/figure-markdown_strict/check distributions of total and soluble fiber-3.png" width="98%" height="98%" />

    plot(density(diet$diet_OCarb_g)) # more normal but still slightly right skewed

<img src="brocc_mediation_files/figure-markdown_strict/check distributions of total and soluble fiber-4.png" width="98%" height="98%" />

    plot(density(diet$diet_Chol_mg)) # more normal but still slightly right skewed

<img src="brocc_mediation_files/figure-markdown_strict/check distributions of total and soluble fiber-5.png" width="98%" height="98%" />

## Run mediation analysis

-   the hypothetical model we tested was whether certain ASVs
    (identified via partial correlation analysis; ASVs correlated to SFN
    and SFN-NIT) mediate the SFN response (SFN and SFN-NIT) to diet
    (based on several dietary components, particularly with interest in
    fiber)

try on one set of each X, Y, and M, then make into for loop based on
visual model here:
<https://library.virginia.edu/data/articles/introduction-to-mediation-analysis>

X &lt;- dietary components

Y &lt;- SFN (3 types)

M &lt;- ASVs

### First, check model 0 is true

-   will use this output to also limit how many mediations I run later
    (i.e., only model 0 combinations of SFN and diet components will
    move forward to mediation analysis)

<!-- -->

    # default lm
    model.0.1 <- lm(sum_SFN ~ diet_TotFib_g, med_input)
    # lm with log of SFN with pseudocount added
    mod_log <- lm(log(sum_SFN+1) ~ diet_TotFib_g, med_input)
    summary(mod_log)
    ## 
    ## Call:
    ## lm(formula = log(sum_SFN + 1) ~ diet_TotFib_g, data = med_input)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -1.34576 -0.40833  0.03989  0.47716  1.19552 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)     2.4579     0.3432   7.162 2.02e-08 ***
    ## diet_TotFib_g -18.5176    31.3780  -0.590    0.559    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.7185 on 36 degrees of freedom
    ## Multiple R-squared:  0.009582,   Adjusted R-squared:  -0.01793 
    ## F-statistic: 0.3483 on 1 and 36 DF,  p-value: 0.5588
    plot(x = med_input$diet_TotFib_g, y = log(med_input$sum_SFN + 1))
    # lm with polynomial
    mod_logpoly <- lm(log(sum_SFN+1) ~ diet_TotFib_g + I(diet_TotFib_g^2) + I(diet_TotFib_g^3) + I(diet_TotFib_g^4), med_input)
    summary(mod_logpoly)
    ## 
    ## Call:
    ## lm(formula = log(sum_SFN + 1) ~ diet_TotFib_g + I(diet_TotFib_g^2) + 
    ##     I(diet_TotFib_g^3) + I(diet_TotFib_g^4), data = med_input)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -1.54348 -0.34172  0.01183  0.55878  1.20569 
    ## 
    ## Coefficients:
    ##                      Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)        -7.518e+00  5.983e+00  -1.257   0.2177  
    ## diet_TotFib_g       4.664e+03  2.553e+03   1.827   0.0767 .
    ## I(diet_TotFib_g^2) -7.415e+05  3.844e+05  -1.929   0.0624 .
    ## I(diet_TotFib_g^3)  4.767e+07  2.425e+07   1.966   0.0578 .
    ## I(diet_TotFib_g^4) -1.066e+09  5.432e+08  -1.963   0.0582 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.6969 on 33 degrees of freedom
    ## Multiple R-squared:  0.1457, Adjusted R-squared:  0.04215 
    ## F-statistic: 1.407 on 4 and 33 DF,  p-value: 0.2532
    plot(x = med_input$diet_TotFib_g, y = log(med_input$sum_SFN + 1))
    pred <- predict(mod_logpoly)
    ix <- sort(med_input$diet_TotFib_g, index.return = T)$ix
    lines(med_input$diet_TotFib_g[ix], pred[ix], col='red', lwd = 2)

<img src="brocc_mediation_files/figure-markdown_strict/choose a model for checking link between SFN/SFN-NIT and the different diet components-1.png" width="98%" height="98%" />

    # gamma glm with pseudocount added to SFN
    mod_gam <- glm((sum_SFN+1) ~ diet_TotFib_g, med_input, family = Gamma)
    summary(mod_gam)
    ## 
    ## Call:
    ## glm(formula = (sum_SFN + 1) ~ diet_TotFib_g, family = Gamma, 
    ##     data = med_input)
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)    0.05965    0.02604   2.291   0.0279 *
    ## diet_TotFib_g  2.29669    2.55339   0.899   0.3744  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Gamma family taken to be 0.4449801)
    ## 
    ##     Null deviance: 17.346  on 37  degrees of freedom
    ## Residual deviance: 16.971  on 36  degrees of freedom
    ## AIC: 258.34
    ## 
    ## Number of Fisher Scoring iterations: 5
    # inverse gaussian flm with pseudocount added to SFN
    mod_invgaus <- glm((sum_SFN+1) ~ diet_TotFib_g, med_input, 
                       family = inverse.gaussian(link = "log"))
    summary(mod_invgaus)
    ## 
    ## Call:
    ## glm(formula = (sum_SFN + 1) ~ diet_TotFib_g, family = inverse.gaussian(link = "log"), 
    ##     data = med_input)
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)     2.7325     0.3197   8.548 3.44e-10 ***
    ## diet_TotFib_g -23.4959    28.3264  -0.829    0.412    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for inverse.gaussian family taken to be 0.03690452)
    ## 
    ##     Null deviance: 1.9252  on 37  degrees of freedom
    ## Residual deviance: 1.8977  on 36  degrees of freedom
    ## AIC: 258.44
    ## 
    ## Number of Fisher Scoring iterations: 7
    # tweedie
    mod_tw <- glm((sum_SFN+1) ~ diet_TotFib_g, med_input,
                  family = tweedie(var.power = 1, link.power = 0))
    summary(mod_tw) # deviance super high
    ## 
    ## Call:
    ## glm(formula = (sum_SFN + 1) ~ diet_TotFib_g, family = tweedie(var.power = 1, 
    ##     link.power = 0), data = med_input)
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)     2.7691     0.3183   8.700 2.23e-10 ***
    ## diet_TotFib_g -27.0621    30.1536  -0.897    0.375    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Tweedie family taken to be 5.396297)
    ## 
    ##     Null deviance: 189.43  on 37  degrees of freedom
    ## Residual deviance: 184.98  on 36  degrees of freedom
    ## AIC: NA
    ## 
    ## Number of Fisher Scoring iterations: 5
    # negative binomial
    mod_negbin <- glm.nb((sum_SFN+1) ~ diet_TotFib_g, med_input)
    summary(mod_negbin)
    ## 
    ## Call:
    ## glm.nb(formula = (sum_SFN + 1) ~ diet_TotFib_g, data = med_input, 
    ##     init.theta = 3.046952046, link = log)
    ## 
    ## Coefficients:
    ##               Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)     2.7537     0.3063   8.991   <2e-16 ***
    ## diet_TotFib_g -25.5397    28.2008  -0.906    0.365    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(3.047) family taken to be 1)
    ## 
    ##     Null deviance: 40.028  on 37  degrees of freedom
    ## Residual deviance: 39.169  on 36  degrees of freedom
    ## AIC: 259.52
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  3.047 
    ##           Std. Err.:  0.860 
    ## 
    ##  2 x log-likelihood:  -253.519

    # compare with aic
    AIC(model.0.1, mod_log, mod_logpoly, mod_gam, mod_invgaus, mod_negbin)
    ##             df       AIC
    ## model.0.1    3 271.03013
    ## mod_log      3  86.65626
    ## mod_logpoly  6  87.03790
    ## mod_gam      3 258.34107
    ## mod_invgaus  3 258.43624
    ## mod_negbin   3 259.51943
    # lm with log(SFN+1) is the best model, but the polynomial is close behind
    # however, adjusted R^2 is higher in the logpoly model, so I should probably use that moving forward

    # model.0 <- lm(Y ~ X)
    model.0.1 <- lm(log(sum_SFN+1) ~ diet_TotFib_g + I(diet_TotFib_g^2) + I(diet_TotFib_g^3) + I(diet_TotFib_g^4), med_input)
    summary(model.0.1)
    ## 
    ## Call:
    ## lm(formula = log(sum_SFN + 1) ~ diet_TotFib_g + I(diet_TotFib_g^2) + 
    ##     I(diet_TotFib_g^3) + I(diet_TotFib_g^4), data = med_input)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -1.54348 -0.34172  0.01183  0.55878  1.20569 
    ## 
    ## Coefficients:
    ##                      Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)        -7.518e+00  5.983e+00  -1.257   0.2177  
    ## diet_TotFib_g       4.664e+03  2.553e+03   1.827   0.0767 .
    ## I(diet_TotFib_g^2) -7.415e+05  3.844e+05  -1.929   0.0624 .
    ## I(diet_TotFib_g^3)  4.767e+07  2.425e+07   1.966   0.0578 .
    ## I(diet_TotFib_g^4) -1.066e+09  5.432e+08  -1.963   0.0582 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.6969 on 33 degrees of freedom
    ## Multiple R-squared:  0.1457, Adjusted R-squared:  0.04215 
    ## F-statistic: 1.407 on 4 and 33 DF,  p-value: 0.2532
    # every term is almost significant

    # diagnostic plots
    autoplot(model.0.1)

<img src="brocc_mediation_files/figure-markdown_strict/check if total fiber affects SFN production-1.png" width="98%" height="98%" />

    # Homoskedastic. Residuals are more or less normally distributed. Residuals are mostly linear.

    # plot
    plot(x = med_input$diet_TotFib_g, y = log(med_input$sum_SFN + 1))
    pred <- predict(model.0.1)
    ix <- sort(med_input$diet_TotFib_g, index.return = T)$ix
    lines(med_input$diet_TotFib_g[ix], pred[ix], col='red', lwd = 2)

<img src="brocc_mediation_files/figure-markdown_strict/check if total fiber affects SFN production-2.png" width="98%" height="98%" />



    model.0.2 <- lm(log(sum_SFN_NIT+1) ~ diet_TotFib_g + I(diet_TotFib_g^2) + I(diet_TotFib_g^3) + I(diet_TotFib_g^4), med_input)
    summary(model.0.2)
    ## 
    ## Call:
    ## lm(formula = log(sum_SFN_NIT + 1) ~ diet_TotFib_g + I(diet_TotFib_g^2) + 
    ##     I(diet_TotFib_g^3) + I(diet_TotFib_g^4), data = med_input)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -1.72027 -0.42765 -0.01992  0.45185  1.19609 
    ## 
    ## Coefficients:
    ##                      Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)        -2.170e+00  5.911e+00  -0.367    0.716
    ## diet_TotFib_g       2.530e+03  2.522e+03   1.003    0.323
    ## I(diet_TotFib_g^2) -3.770e+05  3.798e+05  -0.993    0.328
    ## I(diet_TotFib_g^3)  2.249e+07  2.396e+07   0.939    0.355
    ## I(diet_TotFib_g^4) -4.636e+08  5.367e+08  -0.864    0.394
    ## 
    ## Residual standard error: 0.6886 on 33 degrees of freedom
    ## Multiple R-squared:  0.08722,    Adjusted R-squared:  -0.02342 
    ## F-statistic: 0.7883 on 4 and 33 DF,  p-value: 0.5411
    # no significance

    list_sfn <- c("sum_SFN", "sum_SFN_Cys", "sum_SFN_NAC", "sum_SFN_CG", 
                  "sum_SFN_GSH", "sum_SFN_NIT", "sum_SFN_Tot")
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

    pvals_model.0_lms <- t(as.data.frame(pvals_model.0_lms))
    # a few are significant with an alpha of 0.05

    # want to save insoluble fiber to show coauthors and maybe reviewers
    justinsols <- as.data.frame(pvals_model.0_lms) %>%
      filter(grepl("diet_TotInsolFib_g", rownames(pvals_model.0_lms)))
    write.table(justinsols, "../outputs/insol_pvals.txt", quote = F)

    hist(pvals_model.0_lms, breaks = 20)
    abline(v = quantile(pvals_model.0_lms, 0.05))

<img src="brocc_mediation_files/figure-markdown_strict/loop through all dietary components with each SFN-1.png" width="98%" height="98%" />


    # which ones are those?
    sign_dietSFNcombos <- pvals_model.0_lms[which(pvals_model.0_lms <= 0.05), ]
    sign_dietSFNcombos
    ##                sum_SFN...diet_Alc_g.diet 
    ##                              0.040577890 
    ##               sum_SFN...diet_Alc_g.diet2 
    ##                              0.041217022 
    ##            sum_SFN...diet_Vit_B1_mg.diet 
    ##                              0.011506227 
    ##           sum_SFN...diet_Vit_B1_mg.diet2 
    ##                              0.008644618 
    ##           sum_SFN...diet_Vit_B1_mg.diet3 
    ##                              0.006931764 
    ##           sum_SFN...diet_Vit_B1_mg.diet4 
    ##                              0.006240444 
    ##              sum_SFN...diet_Iron_mg.diet 
    ##                              0.013854151 
    ##             sum_SFN...diet_Iron_mg.diet2 
    ##                              0.015672449 
    ##             sum_SFN...diet_Iron_mg.diet3 
    ##                              0.018110412 
    ##             sum_SFN...diet_Iron_mg.diet4 
    ##                              0.021032034 
    ##            sum_SFN_Cys...diet_Alc_g.diet 
    ##                              0.002241525 
    ##           sum_SFN_Cys...diet_Alc_g.diet2 
    ##                              0.010586294 
    ##           sum_SFN_Cys...diet_Alc_g.diet3 
    ##                              0.049564738 
    ##      sum_SFN_Cys...diet_TotSolFib_g.diet 
    ##                              0.038840661 
    ##  sum_SFN_Cys...diet_Fol.DFE_mcg_DFE.diet 
    ##                              0.023293845 
    ## sum_SFN_Cys...diet_Fol.DFE_mcg_DFE.diet2 
    ##                              0.021059696 
    ## sum_SFN_Cys...diet_Fol.DFE_mcg_DFE.diet3 
    ##                              0.021803070 
    ## sum_SFN_Cys...diet_Fol.DFE_mcg_DFE.diet4 
    ##                              0.023830455 
    ##            sum_SFN_NAC...diet_Alc_g.diet 
    ##                              0.003842307 
    ##           sum_SFN_NAC...diet_Alc_g.diet2 
    ##                              0.010741173 
    ##           sum_SFN_NAC...diet_Alc_g.diet3 
    ##                              0.043078672 
    ##         sum_SFN_NAC...diet_TotFib_g.diet 
    ##                              0.046047626 
    ##        sum_SFN_NAC...diet_TotFib_g.diet2 
    ##                              0.046620371 
    ##           sum_SFN_CG...diet_Chol_mg.diet 
    ##                              0.038950624 
    ##          sum_SFN_CG...diet_Chol_mg.diet2 
    ##                              0.041030579 
    ##             sum_SFN_CG...diet_Alc_g.diet 
    ##                              0.018003057 
    ##            sum_SFN_CG...diet_Alc_g.diet2 
    ##                              0.049165400 
    ##        sum_SFN_GSH...diet_Vit_B1_mg.diet 
    ##                              0.010773993 
    ##       sum_SFN_GSH...diet_Vit_B1_mg.diet2 
    ##                              0.017851640 
    ##       sum_SFN_GSH...diet_Vit_B1_mg.diet3 
    ##                              0.025393989 
    ##       sum_SFN_GSH...diet_Vit_B1_mg.diet4 
    ##                              0.031546921 
    ##            sum_SFN_Tot...diet_Alc_g.diet 
    ##                              0.003884531 
    ##           sum_SFN_Tot...diet_Alc_g.diet2 
    ##                              0.010479773 
    ##           sum_SFN_Tot...diet_Alc_g.diet3 
    ##                              0.044113371 
    ##  sum_SFN_Tot...diet_Fol.DFE_mcg_DFE.diet 
    ##                              0.038037784 
    ## sum_SFN_Tot...diet_Fol.DFE_mcg_DFE.diet2 
    ##                              0.035250927 
    ## sum_SFN_Tot...diet_Fol.DFE_mcg_DFE.diet3 
    ##                              0.037257742 
    ## sum_SFN_Tot...diet_Fol.DFE_mcg_DFE.diet4 
    ##                              0.041233040

    # there are a ton of these that are significant

Pick model 0 as SFN and total dietary fiber to try code on

### Troubleshoot model Y and M that will go into mediation

    ## try to use a tweedie model
    model.Y.1 <- glm(sum_SFN ~ diet_TotFib_g + Blautia, med_input, 
                   family = tweedie(var.power = 1, link.power = 0))
    # no warnings
    summary(model.Y.1)
    # its ok, no warnings now, but the deviance has gone up, no AIC, and more overdispersed

    ## try scaling sfn
    # med_input$SFN_scaled <- round(((med_input$sum_SFN - min(med_input$sum_SFN)) / (max(med_input$sum_SFN) - min(med_input$sum_SFN)) * 100), digits = 0)
    # plot(density(med_input$SFN_scaled))
    # plot(density(med_input$sum_SFN))
    # 
    # model.Y.2 <- glm.nb(SFN_scaled ~ diet_TotFib_g + Blautia, med_input)
    # summary(model.Y.2)
    # plot(x = med_input$Blautia, y = med_input$SFN_scaled)
    # # the same as unscaled so no point in scaling

    ## try log transforming
    model.Y.3 <- glm.nb(log(sum_SFN+1) ~ diet_TotFib_g + Blautia, med_input)
    summary(model.Y.3) # nothing
    model.Y.3b <- glm.nb(sum_SFN ~ diet_TotFib_g + Blautia, med_input)

    ## try lm
    model.Y.4 <- lm(log(sum_SFN+1) ~ diet_TotFib_g + log(Blautia+1), med_input)
    summary(model.Y.4) # nothing

    ## inverse gaussian
    model.Y.5 <- glm(formula = log(sum_SFN + 1) ~ diet_TotFib_g + Blautia, 
                    family = inverse.gaussian(link = "log"), data = med_input)
    summary(model.Y.5)


    ## inv gaussian wo log sfn
    model.Y.5b <- glm(formula = sum_SFN ~ diet_TotFib_g + Blautia, 
                    family = inverse.gaussian(link = "log"), data = med_input)

    ## use AIC to evaluate best option
    AIC(model.Y.0, model.Y.2, model.Y.3, model.Y.3b, model.Y.4, model.Y.5, model.Y.5b)
    AICtweedie( model.Y.1, dispersion = 1, k = 2, verbose = TRUE) # so really not a good fit at all for the data
    # compare

    # model Y4 has a negative adj R2 despite having the lowest AIC, so Y5 is the best choice here

    plot(med_input$quarts_totfib, med_input$f_Erysipelotrichaceae_ASV487)
    plot(med_input$quarts_totfib, log(med_input$f_Erysipelotrichaceae_ASV487 + 1))

    ## basic lm
    model.M.0 <- lm(f_Erysipelotrichaceae_ASV487 ~ diet_TotFib_g, med_input)
    summary(model.M.0)

    ## glm negative binomial
    model.M.1 <- glm.nb(f_Erysipelotrichaceae_ASV487 ~ diet_TotFib_g, med_input)
    summary(model.M.1)

    ## log lm
    model.M.2 <- lm(log(f_Erysipelotrichaceae_ASV487+1) ~ diet_TotFib_g, med_input)
    summary(model.M.2)

    ## zero-infl neg bin
    library(pscl)
    model.M.3 <- zeroinfl(f_Erysipelotrichaceae_ASV487 ~ diet_TotFib_g, med_input, 
                          dist = "negbin")
    summary(model.M.3)

    ## evaluate
    AIC(model.M.0, model.M.1, model.M.2, model.M.3)
    # model M2 is best (tried with a couple ASVs)

    # add quartiles for X (diet_TotFib_g)
    med_input <- med_input %>%
      mutate(quarts_totfib = ntile(diet_TotFib_g, 4),
             SFN_trans = log(sum_SFN+1),
             SFN_Tot_trans = log(sum_SFN_Tot+1),
             Blautia_trans = log(Blautia))
    model.M <- lm(Blautia_trans ~ quarts_totfib, med_input)
    model.Y <- lm(SFN_trans ~ quarts_totfib + Blautia, med_input)

    results <- mediate(model.M, model.Y, treat = "quarts_totfib", mediator = "Blautia",
                       control.value = 1, treat.value = 4, # lowest and highest quartiles
                       boot = TRUE, sims = 500)
    summary(results) # not great but at least now I have a way to implement loop code
    ## 
    ## Causal Mediation Analysis 
    ## 
    ## Nonparametric Bootstrap Confidence Intervals with the Percentile Method
    ## 
    ##                 Estimate 95% CI Lower 95% CI Upper p-value
    ## ACME           -1.81e-05    -1.17e-04         0.00    0.62
    ## ADE            -2.05e-01    -7.69e-01         0.35    0.48
    ## Total Effect   -2.05e-01    -7.69e-01         0.35    0.48
    ## Prop. Mediated  8.85e-05    -1.70e-03         0.00    0.91
    ## 
    ## Sample Size Used: 38 
    ## 
    ## 
    ## Simulations: 500

## Implement loop to cycle through each relevant SFN-dietary combo and each ASV

-   only use significant SFN-diet combos from above

<!-- -->

    # Use list_sfn and diet_compnts as the list of dietary components
    # make SFN taxa a list
    taxa_list <- as.character(SFN_taxa$x)
    # make an output file
    med_res <- c()
    pvals_med <- c()
    exp_res_med <- c()
    list_combos <- as.data.frame(sign_dietSFNcombos) %>%
      rownames_to_column("comparison") %>%
      dplyr::filter(!grepl("diet2|diet3|diet4", comparison)) %>%
      mutate(comparison = str_remove(comparison, ".diet$"),
             sfn = stringr::str_split(comparison, "\\...", simplify = T)[, 1],
             diet = stringr::str_split(comparison, "\\...", simplify = T)[, 2],
             diet = replace(diet, diet == "diet_Fol", "diet_Fol.DFE_mcg_DFE")) %>%
      dplyr::select(-c("comparison", "sign_dietSFNcombos"))

    for (a in 1:nrow(list_combos)) {
        for (f in 1:length(taxa_list)) {
          # add quartiles for each dietary component that will be used
          # transform ASV abundance and transform SFN abundance
          med_input <- med_input %>%
            dplyr::mutate(quarts = ntile(med_input[[list_combos[a,2]]], 4),
                          sfn_trans = log(med_input[[list_combos[a,1]]] + 1),
                          asv_trans = log(med_input[[taxa_list[f]]] + 1),
                          asv = med_input[[taxa_list[f]]])

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
          modsum$p.adj <- p.adjust(modsum$d.avg.p, method = "fdr")
          
          if (modsum$p.adj <= 0.2) {
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
    ## [1] "sum_SFN , Agathobacter , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN , Blautia , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN , CAG.56 , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN , Colidextribacter , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN , Flavonifractor , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN , Incertae.Sedis , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN , Lachnospiraceae.NK4A136.group , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN , Turicibacter , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN , UCG.003 , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN , X.Eubacterium..fissicatena.group , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN , X.Eubacterium..xylanophilum.group , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN , X.Ruminococcus..gauvreauii.group , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN , X.Ruminococcus..gnavus.group , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN , f_Erysipelotrichaceae_ASV487 , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN , f_NA_ASV516 , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN , f_UCG.010_ASV667 , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN , Holdemania , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN , Ruminococcus , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN , Subdoligranulum , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN , Alistipes , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN , Collinsella , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN , Faecalibacterium , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN , NK4A214.group , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN , Roseburia , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN , Agathobacter , diet_Vit_B1_mg are NOT significant"
    ## [1] "sum_SFN , Blautia , diet_Vit_B1_mg are NOT significant"
    ## [1] "sum_SFN , CAG.56 , diet_Vit_B1_mg are NOT significant"
    ## [1] "sum_SFN , Colidextribacter , diet_Vit_B1_mg are NOT significant"
    ## [1] "sum_SFN , Flavonifractor , diet_Vit_B1_mg are NOT significant"
    ## [1] "sum_SFN , Incertae.Sedis , diet_Vit_B1_mg are NOT significant"
    ## [1] "sum_SFN , Lachnospiraceae.NK4A136.group , diet_Vit_B1_mg are NOT significant"
    ## [1] "sum_SFN , Turicibacter , diet_Vit_B1_mg are NOT significant"
    ## [1] "sum_SFN , UCG.003 , diet_Vit_B1_mg are NOT significant"
    ## [1] "sum_SFN , X.Eubacterium..fissicatena.group , diet_Vit_B1_mg are NOT significant"
    ## [1] "sum_SFN , X.Eubacterium..xylanophilum.group , diet_Vit_B1_mg are NOT significant"
    ## [1] "sum_SFN , X.Ruminococcus..gauvreauii.group , diet_Vit_B1_mg are NOT significant"
    ## [1] "sum_SFN , X.Ruminococcus..gnavus.group , diet_Vit_B1_mg are NOT significant"
    ## [1] "sum_SFN , f_Erysipelotrichaceae_ASV487 , diet_Vit_B1_mg are NOT significant"
    ## [1] "sum_SFN , f_NA_ASV516 , diet_Vit_B1_mg are NOT significant"
    ## [1] "sum_SFN , f_UCG.010_ASV667 , diet_Vit_B1_mg are NOT significant"
    ## [1] "sum_SFN , Holdemania , diet_Vit_B1_mg are NOT significant"
    ## [1] "sum_SFN , Ruminococcus , diet_Vit_B1_mg are NOT significant"
    ## [1] "sum_SFN , Subdoligranulum , diet_Vit_B1_mg are NOT significant"
    ## [1] "sum_SFN , Alistipes , diet_Vit_B1_mg are NOT significant"
    ## [1] "sum_SFN , Bifidobacterium , diet_Vit_B1_mg are NOT significant"
    ## [1] "sum_SFN , Faecalibacterium , diet_Vit_B1_mg are NOT significant"
    ## [1] "sum_SFN , NK4A214.group , diet_Vit_B1_mg are NOT significant"
    ## [1] "sum_SFN , Roseburia , diet_Vit_B1_mg are NOT significant"
    ## [1] "sum_SFN , Agathobacter , diet_Iron_mg are NOT significant"
    ## [1] "sum_SFN , Blautia , diet_Iron_mg are NOT significant"
    ## [1] "sum_SFN , CAG.56 , diet_Iron_mg are NOT significant"
    ## [1] "sum_SFN , Colidextribacter , diet_Iron_mg are NOT significant"
    ## [1] "sum_SFN , Flavonifractor , diet_Iron_mg are NOT significant"
    ## [1] "sum_SFN , Incertae.Sedis , diet_Iron_mg are NOT significant"
    ## [1] "sum_SFN , Lachnospiraceae.NK4A136.group , diet_Iron_mg are NOT significant"
    ## [1] "sum_SFN , Turicibacter , diet_Iron_mg are NOT significant"
    ## [1] "sum_SFN , UCG.003 , diet_Iron_mg are NOT significant"
    ## [1] "sum_SFN , X.Eubacterium..fissicatena.group , diet_Iron_mg are NOT significant"
    ## [1] "sum_SFN , X.Eubacterium..xylanophilum.group , diet_Iron_mg are NOT significant"
    ## [1] "sum_SFN , X.Ruminococcus..gauvreauii.group , diet_Iron_mg are NOT significant"
    ## [1] "sum_SFN , X.Ruminococcus..gnavus.group , diet_Iron_mg are NOT significant"
    ## [1] "sum_SFN , f_Erysipelotrichaceae_ASV487 , diet_Iron_mg are NOT significant"
    ## [1] "sum_SFN , f_NA_ASV516 , diet_Iron_mg are NOT significant"
    ## [1] "sum_SFN , f_UCG.010_ASV667 , diet_Iron_mg are NOT significant"
    ## [1] "sum_SFN , Holdemania , diet_Iron_mg are NOT significant"
    ## [1] "sum_SFN , Ruminococcus , diet_Iron_mg are NOT significant"
    ## [1] "sum_SFN , Subdoligranulum , diet_Iron_mg are NOT significant"
    ## [1] "sum_SFN , Alistipes , diet_Iron_mg are NOT significant"
    ## [1] "sum_SFN , Bifidobacterium , diet_Iron_mg are NOT significant"
    ## [1] "sum_SFN , Collinsella , diet_Iron_mg are NOT significant"
    ## [1] "sum_SFN , Faecalibacterium , diet_Iron_mg are NOT significant"
    ## [1] "sum_SFN , NK4A214.group , diet_Iron_mg are NOT significant"
    ## [1] "sum_SFN , Roseburia , diet_Iron_mg are NOT significant"
    ## [1] "sum_SFN_Cys , Agathobacter , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_Cys , Blautia , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_Cys , CAG.56 , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_Cys , Colidextribacter , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_Cys , Flavonifractor , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_Cys , Incertae.Sedis , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_Cys , Lachnospiraceae.NK4A136.group , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_Cys , Turicibacter , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_Cys , UCG.003 , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_Cys , X.Eubacterium..fissicatena.group , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_Cys , X.Eubacterium..xylanophilum.group , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_Cys , X.Ruminococcus..gauvreauii.group , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_Cys , X.Ruminococcus..gnavus.group , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_Cys , f_Erysipelotrichaceae_ASV487 , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_Cys , f_NA_ASV516 , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_Cys , f_UCG.010_ASV667 , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_Cys , Holdemania , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_Cys , Ruminococcus , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_Cys , Subdoligranulum , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_Cys , Alistipes , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_Cys , Bifidobacterium , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_Cys , Collinsella , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_Cys , NK4A214.group , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_Cys , Roseburia , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_Cys , Agathobacter , diet_TotSolFib_g are NOT significant"
    ## [1] "sum_SFN_Cys , Blautia , diet_TotSolFib_g are NOT significant"
    ## [1] "sum_SFN_Cys , CAG.56 , diet_TotSolFib_g are NOT significant"
    ## [1] "sum_SFN_Cys , Colidextribacter , diet_TotSolFib_g are NOT significant"
    ## [1] "sum_SFN_Cys , Flavonifractor , diet_TotSolFib_g are NOT significant"
    ## [1] "sum_SFN_Cys , Incertae.Sedis , diet_TotSolFib_g are NOT significant"
    ## [1] "sum_SFN_Cys , Turicibacter , diet_TotSolFib_g are NOT significant"
    ## [1] "sum_SFN_Cys , UCG.003 , diet_TotSolFib_g are NOT significant"
    ## [1] "sum_SFN_Cys , X.Eubacterium..fissicatena.group , diet_TotSolFib_g are NOT significant"
    ## [1] "sum_SFN_Cys , X.Eubacterium..xylanophilum.group , diet_TotSolFib_g are NOT significant"
    ## [1] "sum_SFN_Cys , X.Ruminococcus..gauvreauii.group , diet_TotSolFib_g are NOT significant"
    ## [1] "sum_SFN_Cys , X.Ruminococcus..gnavus.group , diet_TotSolFib_g are NOT significant"
    ## [1] "sum_SFN_Cys , f_Erysipelotrichaceae_ASV487 , diet_TotSolFib_g are NOT significant"
    ## [1] "sum_SFN_Cys , f_NA_ASV516 , diet_TotSolFib_g are NOT significant"
    ## [1] "sum_SFN_Cys , f_UCG.010_ASV667 , diet_TotSolFib_g are NOT significant"
    ## [1] "sum_SFN_Cys , Holdemania , diet_TotSolFib_g are NOT significant"
    ## [1] "sum_SFN_Cys , Ruminococcus , diet_TotSolFib_g are NOT significant"
    ## [1] "sum_SFN_Cys , Subdoligranulum , diet_TotSolFib_g are NOT significant"
    ## [1] "sum_SFN_Cys , Alistipes , diet_TotSolFib_g are NOT significant"
    ## [1] "sum_SFN_Cys , Bifidobacterium , diet_TotSolFib_g are NOT significant"
    ## [1] "sum_SFN_Cys , Collinsella , diet_TotSolFib_g are NOT significant"
    ## [1] "sum_SFN_Cys , NK4A214.group , diet_TotSolFib_g are NOT significant"
    ## [1] "sum_SFN_Cys , Roseburia , diet_TotSolFib_g are NOT significant"
    ## [1] "sum_SFN_Cys , Agathobacter , diet_Fol.DFE_mcg_DFE are NOT significant"
    ## [1] "sum_SFN_Cys , Blautia , diet_Fol.DFE_mcg_DFE are NOT significant"
    ## [1] "sum_SFN_Cys , CAG.56 , diet_Fol.DFE_mcg_DFE are NOT significant"
    ## [1] "sum_SFN_Cys , Colidextribacter , diet_Fol.DFE_mcg_DFE are NOT significant"
    ## [1] "sum_SFN_Cys , Flavonifractor , diet_Fol.DFE_mcg_DFE are NOT significant"
    ## [1] "sum_SFN_Cys , Incertae.Sedis , diet_Fol.DFE_mcg_DFE are NOT significant"
    ## [1] "sum_SFN_Cys , Lachnospiraceae.NK4A136.group , diet_Fol.DFE_mcg_DFE are NOT significant"
    ## [1] "sum_SFN_Cys , Turicibacter , diet_Fol.DFE_mcg_DFE are NOT significant"
    ## [1] "sum_SFN_Cys , UCG.003 , diet_Fol.DFE_mcg_DFE are NOT significant"
    ## [1] "sum_SFN_Cys , X.Eubacterium..fissicatena.group , diet_Fol.DFE_mcg_DFE are NOT significant"
    ## [1] "sum_SFN_Cys , X.Eubacterium..xylanophilum.group , diet_Fol.DFE_mcg_DFE are NOT significant"
    ## [1] "sum_SFN_Cys , X.Ruminococcus..gnavus.group , diet_Fol.DFE_mcg_DFE are NOT significant"
    ## [1] "sum_SFN_Cys , f_Erysipelotrichaceae_ASV487 , diet_Fol.DFE_mcg_DFE are NOT significant"
    ## [1] "sum_SFN_Cys , f_UCG.010_ASV667 , diet_Fol.DFE_mcg_DFE are NOT significant"
    ## [1] "sum_SFN_Cys , Holdemania , diet_Fol.DFE_mcg_DFE are NOT significant"
    ## [1] "sum_SFN_Cys , Ruminococcus , diet_Fol.DFE_mcg_DFE are NOT significant"
    ## [1] "sum_SFN_Cys , Subdoligranulum , diet_Fol.DFE_mcg_DFE are NOT significant"
    ## [1] "sum_SFN_Cys , Alistipes , diet_Fol.DFE_mcg_DFE are NOT significant"
    ## [1] "sum_SFN_Cys , Bifidobacterium , diet_Fol.DFE_mcg_DFE are NOT significant"
    ## [1] "sum_SFN_Cys , Collinsella , diet_Fol.DFE_mcg_DFE are NOT significant"
    ## [1] "sum_SFN_Cys , Faecalibacterium , diet_Fol.DFE_mcg_DFE are NOT significant"
    ## [1] "sum_SFN_Cys , NK4A214.group , diet_Fol.DFE_mcg_DFE are NOT significant"
    ## [1] "sum_SFN_Cys , Roseburia , diet_Fol.DFE_mcg_DFE are NOT significant"
    ## [1] "sum_SFN_NAC , Agathobacter , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_NAC , Blautia , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_NAC , CAG.56 , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_NAC , Colidextribacter , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_NAC , Flavonifractor , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_NAC , Incertae.Sedis , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_NAC , Lachnospiraceae.NK4A136.group , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_NAC , Turicibacter , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_NAC , UCG.003 , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_NAC , X.Eubacterium..fissicatena.group , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_NAC , X.Eubacterium..xylanophilum.group , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_NAC , X.Ruminococcus..gauvreauii.group , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_NAC , X.Ruminococcus..gnavus.group , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_NAC , f_Erysipelotrichaceae_ASV487 , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_NAC , f_NA_ASV516 , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_NAC , f_UCG.010_ASV667 , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_NAC , Holdemania , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_NAC , Ruminococcus , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_NAC , Subdoligranulum , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_NAC , Alistipes , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_NAC , Bifidobacterium , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_NAC , Collinsella , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_NAC , NK4A214.group , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_NAC , Roseburia , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_NAC , Agathobacter , diet_TotFib_g are NOT significant"
    ## [1] "sum_SFN_NAC , Blautia , diet_TotFib_g are NOT significant"
    ## [1] "sum_SFN_NAC , CAG.56 , diet_TotFib_g are NOT significant"
    ## [1] "sum_SFN_NAC , Colidextribacter , diet_TotFib_g are NOT significant"
    ## [1] "sum_SFN_NAC , Flavonifractor , diet_TotFib_g are NOT significant"
    ## [1] "sum_SFN_NAC , Incertae.Sedis , diet_TotFib_g are NOT significant"
    ## [1] "sum_SFN_NAC , Lachnospiraceae.NK4A136.group , diet_TotFib_g are NOT significant"
    ## [1] "sum_SFN_NAC , Turicibacter , diet_TotFib_g are NOT significant"
    ## [1] "sum_SFN_NAC , UCG.003 , diet_TotFib_g are NOT significant"
    ## [1] "sum_SFN_NAC , X.Eubacterium..fissicatena.group , diet_TotFib_g are NOT significant"
    ## [1] "sum_SFN_NAC , X.Eubacterium..xylanophilum.group , diet_TotFib_g are NOT significant"
    ## [1] "sum_SFN_NAC , X.Ruminococcus..gauvreauii.group , diet_TotFib_g are NOT significant"
    ## [1] "sum_SFN_NAC , X.Ruminococcus..gnavus.group , diet_TotFib_g are NOT significant"
    ## [1] "sum_SFN_NAC , f_Erysipelotrichaceae_ASV487 , diet_TotFib_g are NOT significant"
    ## [1] "sum_SFN_NAC , f_NA_ASV516 , diet_TotFib_g are NOT significant"
    ## [1] "sum_SFN_NAC , f_UCG.010_ASV667 , diet_TotFib_g are NOT significant"
    ## [1] "sum_SFN_NAC , Holdemania , diet_TotFib_g are NOT significant"
    ## [1] "sum_SFN_NAC , Ruminococcus , diet_TotFib_g are NOT significant"
    ## [1] "sum_SFN_NAC , Subdoligranulum , diet_TotFib_g are NOT significant"
    ## [1] "sum_SFN_NAC , Alistipes , diet_TotFib_g are NOT significant"
    ## [1] "sum_SFN_NAC , Bifidobacterium , diet_TotFib_g are NOT significant"
    ## [1] "sum_SFN_NAC , Faecalibacterium , diet_TotFib_g are NOT significant"
    ## [1] "sum_SFN_NAC , NK4A214.group , diet_TotFib_g are NOT significant"
    ## [1] "sum_SFN_NAC , Roseburia , diet_TotFib_g are NOT significant"
    ## [1] "sum_SFN_CG , Agathobacter , diet_Chol_mg are NOT significant"
    ## [1] "sum_SFN_CG , Blautia , diet_Chol_mg are NOT significant"
    ## [1] "sum_SFN_CG , CAG.56 , diet_Chol_mg are NOT significant"
    ## [1] "sum_SFN_CG , Colidextribacter , diet_Chol_mg are NOT significant"
    ## [1] "sum_SFN_CG , Flavonifractor , diet_Chol_mg are NOT significant"
    ## [1] "sum_SFN_CG , Incertae.Sedis , diet_Chol_mg are NOT significant"
    ## [1] "sum_SFN_CG , Lachnospiraceae.NK4A136.group , diet_Chol_mg are NOT significant"
    ## [1] "sum_SFN_CG , Turicibacter , diet_Chol_mg are NOT significant"
    ## [1] "sum_SFN_CG , UCG.003 , diet_Chol_mg are NOT significant"
    ## [1] "sum_SFN_CG , X.Eubacterium..fissicatena.group , diet_Chol_mg are NOT significant"
    ## [1] "sum_SFN_CG , X.Eubacterium..xylanophilum.group , diet_Chol_mg are NOT significant"
    ## [1] "sum_SFN_CG , X.Ruminococcus..gauvreauii.group , diet_Chol_mg are NOT significant"
    ## [1] "sum_SFN_CG , X.Ruminococcus..gnavus.group , diet_Chol_mg are NOT significant"
    ## [1] "sum_SFN_CG , f_Erysipelotrichaceae_ASV487 , diet_Chol_mg are NOT significant"
    ## [1] "sum_SFN_CG , f_NA_ASV516 , diet_Chol_mg are NOT significant"
    ## [1] "sum_SFN_CG , f_UCG.010_ASV667 , diet_Chol_mg are NOT significant"
    ## [1] "sum_SFN_CG , Holdemania , diet_Chol_mg are NOT significant"
    ## [1] "sum_SFN_CG , Ruminococcus , diet_Chol_mg are NOT significant"
    ## [1] "sum_SFN_CG , Alistipes , diet_Chol_mg are NOT significant"
    ## [1] "sum_SFN_CG , Bifidobacterium , diet_Chol_mg are NOT significant"
    ## [1] "sum_SFN_CG , Collinsella , diet_Chol_mg are NOT significant"
    ## [1] "sum_SFN_CG , Faecalibacterium , diet_Chol_mg are NOT significant"
    ## [1] "sum_SFN_CG , NK4A214.group , diet_Chol_mg are NOT significant"
    ## [1] "sum_SFN_CG , Roseburia , diet_Chol_mg are NOT significant"
    ## [1] "sum_SFN_CG , Agathobacter , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_CG , Blautia , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_CG , CAG.56 , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_CG , Colidextribacter , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_CG , Flavonifractor , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_CG , Incertae.Sedis , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_CG , Lachnospiraceae.NK4A136.group , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_CG , Turicibacter , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_CG , UCG.003 , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_CG , X.Eubacterium..fissicatena.group , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_CG , X.Eubacterium..xylanophilum.group , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_CG , X.Ruminococcus..gauvreauii.group , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_CG , X.Ruminococcus..gnavus.group , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_CG , f_Erysipelotrichaceae_ASV487 , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_CG , f_NA_ASV516 , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_CG , f_UCG.010_ASV667 , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_CG , Holdemania , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_CG , Ruminococcus , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_CG , Subdoligranulum , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_CG , Alistipes , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_CG , Bifidobacterium , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_CG , Collinsella , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_CG , Faecalibacterium , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_CG , NK4A214.group , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_CG , Roseburia , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_GSH , Agathobacter , diet_Vit_B1_mg are NOT significant"
    ## [1] "sum_SFN_GSH , Blautia , diet_Vit_B1_mg are NOT significant"
    ## [1] "sum_SFN_GSH , CAG.56 , diet_Vit_B1_mg are NOT significant"
    ## [1] "sum_SFN_GSH , Colidextribacter , diet_Vit_B1_mg are NOT significant"
    ## [1] "sum_SFN_GSH , Flavonifractor , diet_Vit_B1_mg are NOT significant"
    ## [1] "sum_SFN_GSH , Incertae.Sedis , diet_Vit_B1_mg are NOT significant"
    ## [1] "sum_SFN_GSH , Lachnospiraceae.NK4A136.group , diet_Vit_B1_mg are NOT significant"
    ## [1] "sum_SFN_GSH , Turicibacter , diet_Vit_B1_mg are NOT significant"
    ## [1] "sum_SFN_GSH , UCG.003 , diet_Vit_B1_mg are NOT significant"
    ## [1] "sum_SFN_GSH , X.Eubacterium..fissicatena.group , diet_Vit_B1_mg are NOT significant"
    ## [1] "sum_SFN_GSH , X.Eubacterium..xylanophilum.group , diet_Vit_B1_mg are NOT significant"
    ## [1] "sum_SFN_GSH , X.Ruminococcus..gauvreauii.group , diet_Vit_B1_mg are NOT significant"
    ## [1] "sum_SFN_GSH , f_Erysipelotrichaceae_ASV487 , diet_Vit_B1_mg are NOT significant"
    ## [1] "sum_SFN_GSH , f_NA_ASV516 , diet_Vit_B1_mg are NOT significant"
    ## [1] "sum_SFN_GSH , f_UCG.010_ASV667 , diet_Vit_B1_mg are NOT significant"
    ## [1] "sum_SFN_GSH , Holdemania , diet_Vit_B1_mg are NOT significant"
    ## [1] "sum_SFN_GSH , Ruminococcus , diet_Vit_B1_mg are NOT significant"
    ## [1] "sum_SFN_GSH , Subdoligranulum , diet_Vit_B1_mg are NOT significant"
    ## [1] "sum_SFN_GSH , Alistipes , diet_Vit_B1_mg are NOT significant"
    ## [1] "sum_SFN_GSH , Bifidobacterium , diet_Vit_B1_mg are NOT significant"
    ## [1] "sum_SFN_GSH , Collinsella , diet_Vit_B1_mg are NOT significant"
    ## [1] "sum_SFN_GSH , Faecalibacterium , diet_Vit_B1_mg are NOT significant"
    ## [1] "sum_SFN_GSH , NK4A214.group , diet_Vit_B1_mg are NOT significant"
    ## [1] "sum_SFN_GSH , Roseburia , diet_Vit_B1_mg are NOT significant"
    ## [1] "sum_SFN_Tot , Agathobacter , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_Tot , Blautia , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_Tot , CAG.56 , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_Tot , Colidextribacter , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_Tot , Flavonifractor , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_Tot , Incertae.Sedis , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_Tot , Lachnospiraceae.NK4A136.group , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_Tot , Turicibacter , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_Tot , UCG.003 , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_Tot , X.Eubacterium..fissicatena.group , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_Tot , X.Eubacterium..xylanophilum.group , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_Tot , X.Ruminococcus..gauvreauii.group , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_Tot , X.Ruminococcus..gnavus.group , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_Tot , f_Erysipelotrichaceae_ASV487 , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_Tot , f_NA_ASV516 , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_Tot , f_UCG.010_ASV667 , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_Tot , Holdemania , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_Tot , Ruminococcus , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_Tot , Subdoligranulum , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_Tot , Alistipes , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_Tot , Bifidobacterium , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_Tot , Collinsella , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_Tot , NK4A214.group , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_Tot , Roseburia , diet_Alc_g are NOT significant"
    ## [1] "sum_SFN_Tot , Agathobacter , diet_Fol.DFE_mcg_DFE are NOT significant"
    ## [1] "sum_SFN_Tot , Blautia , diet_Fol.DFE_mcg_DFE are NOT significant"
    ## [1] "sum_SFN_Tot , CAG.56 , diet_Fol.DFE_mcg_DFE are NOT significant"
    ## [1] "sum_SFN_Tot , Colidextribacter , diet_Fol.DFE_mcg_DFE are NOT significant"
    ## [1] "sum_SFN_Tot , Flavonifractor , diet_Fol.DFE_mcg_DFE are NOT significant"
    ## [1] "sum_SFN_Tot , Incertae.Sedis , diet_Fol.DFE_mcg_DFE are NOT significant"
    ## [1] "sum_SFN_Tot , Lachnospiraceae.NK4A136.group , diet_Fol.DFE_mcg_DFE are NOT significant"
    ## [1] "sum_SFN_Tot , Turicibacter , diet_Fol.DFE_mcg_DFE are NOT significant"
    ## [1] "sum_SFN_Tot , UCG.003 , diet_Fol.DFE_mcg_DFE are NOT significant"
    ## [1] "sum_SFN_Tot , X.Eubacterium..fissicatena.group , diet_Fol.DFE_mcg_DFE are NOT significant"
    ## [1] "sum_SFN_Tot , X.Eubacterium..xylanophilum.group , diet_Fol.DFE_mcg_DFE are NOT significant"
    ## [1] "sum_SFN_Tot , X.Ruminococcus..gauvreauii.group , diet_Fol.DFE_mcg_DFE are NOT significant"
    ## [1] "sum_SFN_Tot , X.Ruminococcus..gnavus.group , diet_Fol.DFE_mcg_DFE are NOT significant"
    ## [1] "sum_SFN_Tot , f_Erysipelotrichaceae_ASV487 , diet_Fol.DFE_mcg_DFE are NOT significant"
    ## [1] "sum_SFN_Tot , Holdemania , diet_Fol.DFE_mcg_DFE are NOT significant"
    ## [1] "sum_SFN_Tot , Ruminococcus , diet_Fol.DFE_mcg_DFE are NOT significant"
    ## [1] "sum_SFN_Tot , Subdoligranulum , diet_Fol.DFE_mcg_DFE are NOT significant"
    ## [1] "sum_SFN_Tot , Alistipes , diet_Fol.DFE_mcg_DFE are NOT significant"
    ## [1] "sum_SFN_Tot , Bifidobacterium , diet_Fol.DFE_mcg_DFE are NOT significant"
    ## [1] "sum_SFN_Tot , Faecalibacterium , diet_Fol.DFE_mcg_DFE are NOT significant"
    ## [1] "sum_SFN_Tot , NK4A214.group , diet_Fol.DFE_mcg_DFE are NOT significant"
    ## [1] "sum_SFN_Tot , Roseburia , diet_Fol.DFE_mcg_DFE are NOT significant"

    # saveRDS(med_res, "mediation_result_notconservative.Rdata")
    # saveRDS(pvals_med, "pvals_mediation_notconservative.Rdata")
    # saveRDS(exp_res_med, "expanded_mediation_results_notconservative.Rdata")

    # saveRDS(med_res, "mediation_result_signifsfndiets.Rdata")
    # saveRDS(pvals_med, "pvals_mediation_signifsfndiets.Rdata")
    # saveRDS(exp_res_med, "expanded_mediation_results_signifsfndiets.Rdata")

    saveRDS(med_res, "mediation_result_signifcumSFNs.Rdata")
    saveRDS(pvals_med, "pvals_mediation_signifcumSFNs.Rdata")
    saveRDS(exp_res_med, "expanded_mediation_results_signifcumSFNs.Rdata")

## Examine results

    signif <- do.call(rbind.data.frame, pvals_med[names(med_res)]) %>%
      rownames_to_column("names(med_res)")
    signif
    ##                                                       names(med_res) ACME_pval
    ## 1                                 sum_SFN Bifidobacterium diet_Alc_g     0.192
    ## 2                                 sum_SFN Collinsella diet_Vit_B1_mg     0.092
    ## 3                            sum_SFN_Cys Faecalibacterium diet_Alc_g     0.196
    ## 4         sum_SFN_Cys Lachnospiraceae.NK4A136.group diet_TotSolFib_g     0.200
    ## 5                      sum_SFN_Cys Faecalibacterium diet_TotSolFib_g     0.056
    ## 6  sum_SFN_Cys X.Ruminococcus..gauvreauii.group diet_Fol.DFE_mcg_DFE     0.144
    ## 7                       sum_SFN_Cys f_NA_ASV516 diet_Fol.DFE_mcg_DFE     0.128
    ## 8                            sum_SFN_NAC Faecalibacterium diet_Alc_g     0.172
    ## 9                              sum_SFN_NAC Collinsella diet_TotFib_g     0.072
    ## 10                           sum_SFN_CG Subdoligranulum diet_Chol_mg     0.056
    ## 11           sum_SFN_GSH X.Ruminococcus..gnavus.group diet_Vit_B1_mg     0.136
    ## 12                           sum_SFN_Tot Faecalibacterium diet_Alc_g     0.164
    ## 13                      sum_SFN_Tot f_NA_ASV516 diet_Fol.DFE_mcg_DFE     0.116
    ## 14                 sum_SFN_Tot f_UCG.010_ASV667 diet_Fol.DFE_mcg_DFE     0.152
    ## 15                      sum_SFN_Tot Collinsella diet_Fol.DFE_mcg_DFE     0.112
    ##    ADE_pval prop_mediated_pval
    ## 1     0.116              0.528
    ## 2     0.240              0.684
    ## 3     0.540              0.436
    ## 4     0.284              0.604
    ## 5     0.908              0.544
    ## 6     0.772              0.828
    ## 7     0.672              0.900
    ## 8     0.960              0.692
    ## 9     0.828              0.684
    ## 10    0.032              0.216
    ## 11    0.676              0.384
    ## 12    0.804              0.852
    ## 13    0.696              0.796
    ## 14    0.980              0.736
    ## 15    0.336              0.828

    # write.table(signif, "../outputs/signif_meds_cumsfn.txt")

    med_res
    ## $`sum_SFN Bifidobacterium diet_Alc_g`
    ## $`sum_SFN Bifidobacterium diet_Alc_g`$ACME_est
    ## [1] 0.1367869
    ## 
    ## $`sum_SFN Bifidobacterium diet_Alc_g`$ACME_pval
    ## [1] 0.192
    ## 
    ## $`sum_SFN Bifidobacterium diet_Alc_g`$ADE_est
    ## [1] -0.3789608
    ## 
    ## $`sum_SFN Bifidobacterium diet_Alc_g`$ADE_pval
    ## [1] 0.116
    ## 
    ## $`sum_SFN Bifidobacterium diet_Alc_g`$prop_mediated_est
    ## [1] -0.5648293
    ## 
    ## $`sum_SFN Bifidobacterium diet_Alc_g`$prop_mediated_pval
    ## [1] 0.528
    ## 
    ## 
    ## $`sum_SFN Collinsella diet_Vit_B1_mg`
    ## $`sum_SFN Collinsella diet_Vit_B1_mg`$ACME_est
    ## [1] -0.2192413
    ## 
    ## $`sum_SFN Collinsella diet_Vit_B1_mg`$ACME_pval
    ## [1] 0.092
    ## 
    ## $`sum_SFN Collinsella diet_Vit_B1_mg`$ADE_est
    ## [1] 0.3734043
    ## 
    ## $`sum_SFN Collinsella diet_Vit_B1_mg`$ADE_pval
    ## [1] 0.24
    ## 
    ## $`sum_SFN Collinsella diet_Vit_B1_mg`$prop_mediated_est
    ## [1] -1.422139
    ## 
    ## $`sum_SFN Collinsella diet_Vit_B1_mg`$prop_mediated_pval
    ## [1] 0.684
    ## 
    ## 
    ## $`sum_SFN_Cys Faecalibacterium diet_Alc_g`
    ## $`sum_SFN_Cys Faecalibacterium diet_Alc_g`$ACME_est
    ## [1] 0.1239119
    ## 
    ## $`sum_SFN_Cys Faecalibacterium diet_Alc_g`$ACME_pval
    ## [1] 0.196
    ## 
    ## $`sum_SFN_Cys Faecalibacterium diet_Alc_g`$ADE_est
    ## [1] 0.2144161
    ## 
    ## $`sum_SFN_Cys Faecalibacterium diet_Alc_g`$ADE_pval
    ## [1] 0.54
    ## 
    ## $`sum_SFN_Cys Faecalibacterium diet_Alc_g`$prop_mediated_est
    ## [1] 0.3662479
    ## 
    ## $`sum_SFN_Cys Faecalibacterium diet_Alc_g`$prop_mediated_pval
    ## [1] 0.436
    ## 
    ## 
    ## $`sum_SFN_Cys Lachnospiraceae.NK4A136.group diet_TotSolFib_g`
    ## $`sum_SFN_Cys Lachnospiraceae.NK4A136.group diet_TotSolFib_g`$ACME_est
    ## [1] -0.1265578
    ## 
    ## $`sum_SFN_Cys Lachnospiraceae.NK4A136.group diet_TotSolFib_g`$ACME_pval
    ## [1] 0.2
    ## 
    ## $`sum_SFN_Cys Lachnospiraceae.NK4A136.group diet_TotSolFib_g`$ADE_est
    ## [1] 0.3242522
    ## 
    ## $`sum_SFN_Cys Lachnospiraceae.NK4A136.group diet_TotSolFib_g`$ADE_pval
    ## [1] 0.284
    ## 
    ## $`sum_SFN_Cys Lachnospiraceae.NK4A136.group diet_TotSolFib_g`$prop_mediated_est
    ## [1] -0.6401692
    ## 
    ## $`sum_SFN_Cys Lachnospiraceae.NK4A136.group diet_TotSolFib_g`$prop_mediated_pval
    ## [1] 0.604
    ## 
    ## 
    ## $`sum_SFN_Cys Faecalibacterium diet_TotSolFib_g`
    ## $`sum_SFN_Cys Faecalibacterium diet_TotSolFib_g`$ACME_est
    ## [1] 0.2260267
    ## 
    ## $`sum_SFN_Cys Faecalibacterium diet_TotSolFib_g`$ACME_pval
    ## [1] 0.056
    ## 
    ## $`sum_SFN_Cys Faecalibacterium diet_TotSolFib_g`$ADE_est
    ## [1] -0.02833233
    ## 
    ## $`sum_SFN_Cys Faecalibacterium diet_TotSolFib_g`$ADE_pval
    ## [1] 0.908
    ## 
    ## $`sum_SFN_Cys Faecalibacterium diet_TotSolFib_g`$prop_mediated_est
    ## [1] 1.143314
    ## 
    ## $`sum_SFN_Cys Faecalibacterium diet_TotSolFib_g`$prop_mediated_pval
    ## [1] 0.544
    ## 
    ## 
    ## $`sum_SFN_Cys X.Ruminococcus..gauvreauii.group diet_Fol.DFE_mcg_DFE`
    ## $`sum_SFN_Cys X.Ruminococcus..gauvreauii.group diet_Fol.DFE_mcg_DFE`$ACME_est
    ## [1] 0.1523865
    ## 
    ## $`sum_SFN_Cys X.Ruminococcus..gauvreauii.group diet_Fol.DFE_mcg_DFE`$ACME_pval
    ## [1] 0.144
    ## 
    ## $`sum_SFN_Cys X.Ruminococcus..gauvreauii.group diet_Fol.DFE_mcg_DFE`$ADE_est
    ## [1] -0.09947036
    ## 
    ## $`sum_SFN_Cys X.Ruminococcus..gauvreauii.group diet_Fol.DFE_mcg_DFE`$ADE_pval
    ## [1] 0.772
    ## 
    ## $`sum_SFN_Cys X.Ruminococcus..gauvreauii.group diet_Fol.DFE_mcg_DFE`$prop_mediated_est
    ## [1] 2.879773
    ## 
    ## $`sum_SFN_Cys X.Ruminococcus..gauvreauii.group diet_Fol.DFE_mcg_DFE`$prop_mediated_pval
    ## [1] 0.828
    ## 
    ## 
    ## $`sum_SFN_Cys f_NA_ASV516 diet_Fol.DFE_mcg_DFE`
    ## $`sum_SFN_Cys f_NA_ASV516 diet_Fol.DFE_mcg_DFE`$ACME_est
    ## [1] 0.1850021
    ## 
    ## $`sum_SFN_Cys f_NA_ASV516 diet_Fol.DFE_mcg_DFE`$ACME_pval
    ## [1] 0.128
    ## 
    ## $`sum_SFN_Cys f_NA_ASV516 diet_Fol.DFE_mcg_DFE`$ADE_est
    ## [1] -0.1320859
    ## 
    ## $`sum_SFN_Cys f_NA_ASV516 diet_Fol.DFE_mcg_DFE`$ADE_pval
    ## [1] 0.672
    ## 
    ## $`sum_SFN_Cys f_NA_ASV516 diet_Fol.DFE_mcg_DFE`$prop_mediated_est
    ## [1] 3.496136
    ## 
    ## $`sum_SFN_Cys f_NA_ASV516 diet_Fol.DFE_mcg_DFE`$prop_mediated_pval
    ## [1] 0.9
    ## 
    ## 
    ## $`sum_SFN_NAC Faecalibacterium diet_Alc_g`
    ## $`sum_SFN_NAC Faecalibacterium diet_Alc_g`$ACME_est
    ## [1] 0.111018
    ## 
    ## $`sum_SFN_NAC Faecalibacterium diet_Alc_g`$ACME_pval
    ## [1] 0.172
    ## 
    ## $`sum_SFN_NAC Faecalibacterium diet_Alc_g`$ADE_est
    ## [1] -0.03740561
    ## 
    ## $`sum_SFN_NAC Faecalibacterium diet_Alc_g`$ADE_pval
    ## [1] 0.96
    ## 
    ## $`sum_SFN_NAC Faecalibacterium diet_Alc_g`$prop_mediated_est
    ## [1] 1.508143
    ## 
    ## $`sum_SFN_NAC Faecalibacterium diet_Alc_g`$prop_mediated_pval
    ## [1] 0.692
    ## 
    ## 
    ## $`sum_SFN_NAC Collinsella diet_TotFib_g`
    ## $`sum_SFN_NAC Collinsella diet_TotFib_g`$ACME_est
    ## [1] -0.1419928
    ## 
    ## $`sum_SFN_NAC Collinsella diet_TotFib_g`$ACME_pval
    ## [1] 0.072
    ## 
    ## $`sum_SFN_NAC Collinsella diet_TotFib_g`$ADE_est
    ## [1] 0.05842425
    ## 
    ## $`sum_SFN_NAC Collinsella diet_TotFib_g`$ADE_pval
    ## [1] 0.828
    ## 
    ## $`sum_SFN_NAC Collinsella diet_TotFib_g`$prop_mediated_est
    ## [1] 1.699118
    ## 
    ## $`sum_SFN_NAC Collinsella diet_TotFib_g`$prop_mediated_pval
    ## [1] 0.684
    ## 
    ## 
    ## $`sum_SFN_CG Subdoligranulum diet_Chol_mg`
    ## $`sum_SFN_CG Subdoligranulum diet_Chol_mg`$ACME_est
    ## [1] 0.02209599
    ## 
    ## $`sum_SFN_CG Subdoligranulum diet_Chol_mg`$ACME_pval
    ## [1] 0.056
    ## 
    ## $`sum_SFN_CG Subdoligranulum diet_Chol_mg`$ADE_est
    ## [1] -0.1135696
    ## 
    ## $`sum_SFN_CG Subdoligranulum diet_Chol_mg`$ADE_pval
    ## [1] 0.032
    ## 
    ## $`sum_SFN_CG Subdoligranulum diet_Chol_mg`$prop_mediated_est
    ## [1] -0.2415559
    ## 
    ## $`sum_SFN_CG Subdoligranulum diet_Chol_mg`$prop_mediated_pval
    ## [1] 0.216
    ## 
    ## 
    ## $`sum_SFN_GSH X.Ruminococcus..gnavus.group diet_Vit_B1_mg`
    ## $`sum_SFN_GSH X.Ruminococcus..gnavus.group diet_Vit_B1_mg`$ACME_est
    ## [1] -0.00955483
    ## 
    ## $`sum_SFN_GSH X.Ruminococcus..gnavus.group diet_Vit_B1_mg`$ACME_pval
    ## [1] 0.136
    ## 
    ## $`sum_SFN_GSH X.Ruminococcus..gnavus.group diet_Vit_B1_mg`$ADE_est
    ## [1] -0.01034253
    ## 
    ## $`sum_SFN_GSH X.Ruminococcus..gnavus.group diet_Vit_B1_mg`$ADE_pval
    ## [1] 0.676
    ## 
    ## $`sum_SFN_GSH X.Ruminococcus..gnavus.group diet_Vit_B1_mg`$prop_mediated_est
    ## [1] 0.480206
    ## 
    ## $`sum_SFN_GSH X.Ruminococcus..gnavus.group diet_Vit_B1_mg`$prop_mediated_pval
    ## [1] 0.384
    ## 
    ## 
    ## $`sum_SFN_Tot Faecalibacterium diet_Alc_g`
    ## $`sum_SFN_Tot Faecalibacterium diet_Alc_g`$ACME_est
    ## [1] 0.1082566
    ## 
    ## $`sum_SFN_Tot Faecalibacterium diet_Alc_g`$ACME_pval
    ## [1] 0.164
    ## 
    ## $`sum_SFN_Tot Faecalibacterium diet_Alc_g`$ADE_est
    ## [1] -0.06525081
    ## 
    ## $`sum_SFN_Tot Faecalibacterium diet_Alc_g`$ADE_pval
    ## [1] 0.804
    ## 
    ## $`sum_SFN_Tot Faecalibacterium diet_Alc_g`$prop_mediated_est
    ## [1] 2.517256
    ## 
    ## $`sum_SFN_Tot Faecalibacterium diet_Alc_g`$prop_mediated_pval
    ## [1] 0.852
    ## 
    ## 
    ## $`sum_SFN_Tot f_NA_ASV516 diet_Fol.DFE_mcg_DFE`
    ## $`sum_SFN_Tot f_NA_ASV516 diet_Fol.DFE_mcg_DFE`$ACME_est
    ## [1] 0.1387401
    ## 
    ## $`sum_SFN_Tot f_NA_ASV516 diet_Fol.DFE_mcg_DFE`$ACME_pval
    ## [1] 0.116
    ## 
    ## $`sum_SFN_Tot f_NA_ASV516 diet_Fol.DFE_mcg_DFE`$ADE_est
    ## [1] -0.06569426
    ## 
    ## $`sum_SFN_Tot f_NA_ASV516 diet_Fol.DFE_mcg_DFE`$ADE_pval
    ## [1] 0.696
    ## 
    ## $`sum_SFN_Tot f_NA_ASV516 diet_Fol.DFE_mcg_DFE`$prop_mediated_est
    ## [1] 1.899357
    ## 
    ## $`sum_SFN_Tot f_NA_ASV516 diet_Fol.DFE_mcg_DFE`$prop_mediated_pval
    ## [1] 0.796
    ## 
    ## 
    ## $`sum_SFN_Tot f_UCG.010_ASV667 diet_Fol.DFE_mcg_DFE`
    ## $`sum_SFN_Tot f_UCG.010_ASV667 diet_Fol.DFE_mcg_DFE`$ACME_est
    ## [1] 0.08795385
    ## 
    ## $`sum_SFN_Tot f_UCG.010_ASV667 diet_Fol.DFE_mcg_DFE`$ACME_pval
    ## [1] 0.152
    ## 
    ## $`sum_SFN_Tot f_UCG.010_ASV667 diet_Fol.DFE_mcg_DFE`$ADE_est
    ## [1] -0.01490804
    ## 
    ## $`sum_SFN_Tot f_UCG.010_ASV667 diet_Fol.DFE_mcg_DFE`$ADE_pval
    ## [1] 0.98
    ## 
    ## $`sum_SFN_Tot f_UCG.010_ASV667 diet_Fol.DFE_mcg_DFE`$prop_mediated_est
    ## [1] 1.204092
    ## 
    ## $`sum_SFN_Tot f_UCG.010_ASV667 diet_Fol.DFE_mcg_DFE`$prop_mediated_pval
    ## [1] 0.736
    ## 
    ## 
    ## $`sum_SFN_Tot Collinsella diet_Fol.DFE_mcg_DFE`
    ## $`sum_SFN_Tot Collinsella diet_Fol.DFE_mcg_DFE`$ACME_est
    ## [1] -0.1332669
    ## 
    ## $`sum_SFN_Tot Collinsella diet_Fol.DFE_mcg_DFE`$ACME_pval
    ## [1] 0.112
    ## 
    ## $`sum_SFN_Tot Collinsella diet_Fol.DFE_mcg_DFE`$ADE_est
    ## [1] 0.2063127
    ## 
    ## $`sum_SFN_Tot Collinsella diet_Fol.DFE_mcg_DFE`$ADE_pval
    ## [1] 0.336
    ## 
    ## $`sum_SFN_Tot Collinsella diet_Fol.DFE_mcg_DFE`$prop_mediated_est
    ## [1] -1.82443
    ## 
    ## $`sum_SFN_Tot Collinsella diet_Fol.DFE_mcg_DFE`$prop_mediated_pval
    ## [1] 0.828

    # write a for loop with the above code to make graphs of all of them
    medsignif_combos <- names(med_res) # 17 of these

    # make output plot list
    plot_list <- list()

    # for loop to make graphs
    for (i in 1:length(medsignif_combos)) {
      # first separate name to be three names
      name <- medsignif_combos[i]
      X <- strsplit(name, " ")[[1]][1]
      M <- strsplit(name, " ")[[1]][2]
      Y <- strsplit(name, " ")[[1]][3]
      
      # make a filtered data frame for ggplot to use
      regrplot <- med_input %>% 
        dplyr::select(all_of(M), all_of(X), all_of(Y))
      regrplot$trans_M <- log(regrplot[[M]]+1)
      regrplot$trans_X <- log(regrplot[[X]]+1)
      names(regrplot)[names(regrplot) == M] <- "M_vals"
      names(regrplot)[names(regrplot) == X] <- "X_vals"
      names(regrplot)[names(regrplot) == Y] <- "Y_vals"
      
      # make graph
      if(i == 1) {
        regrplot <- regrplot %>%
          dplyr::filter(regrplot$Y_vals != max(regrplot$Y_vals))
        }
      p <- ggplot(data = regrplot, aes(y = trans_X, x = trans_M)) +
        geom_point(aes(color = Y_vals)) +
        geom_smooth(se = TRUE, method = "lm") +
        labs(y = paste0(X, " quantity"),
           x = "ASV abundance",
           color = paste0("log transformed \n", Y, " quantity"), title = M) +
        theme(text = element_text(size = 14))
      
      # add to plot list
      plot_list[[i]] <- p
      # save
      ggsave(paste0("/Users/alexieva/Documents/Projects/Analysis/broccoli project/02_pubAnalysis/outputs/P1_medres/P1MedRes_plot_", i, ".png"),
             plot_list[[i]],
             height = 5, width = 7, dpi = 300)
      
      
    }

    g <- ggpubr::ggarrange(plotlist = plot_list,
              ncol = 5, nrow = 4,
              common.legend = T,
              legend = "bottom")

    ggsave("/Users/alexieva/Documents/Projects/Analysis/broccoli project/02_pubAnalysis/outputs/P1MedRes_plotALL.png",
            g,
            height = 14, width = 18, dpi = 300)

    # write a for loop with the above code to make graphs of all of them
    medsignif_combos <- names(med_res) # 17 of these

    # make output plot list
    plot_list <- list()

    # for loop to make graphs
    for (i in 1:length(medsignif_combos)) {
      # first separate name to be three names
      name <- medsignif_combos[i]
      X <- strsplit(name, " ")[[1]][1]
      M <- strsplit(name, " ")[[1]][2]
      Y <- strsplit(name, " ")[[1]][3]
      
      # make a filtered data frame for ggplot to use
      regrplot <- med_input %>% 
        dplyr::select(all_of(M), all_of(X), all_of(Y))
      regrplot$trans_M <- log(regrplot[[M]]+1)
      regrplot$trans_X <- log(regrplot[[X]]+1)
      names(regrplot)[names(regrplot) == M] <- "M_vals"
      names(regrplot)[names(regrplot) == X] <- "X_vals"
      names(regrplot)[names(regrplot) == Y] <- "Y_vals"
      
      # make graph
      if(i == 1) {
        regrplot <- regrplot %>%
          dplyr::filter(regrplot$Y_vals != max(regrplot$Y_vals))
        }
      p <- ggplot(data = regrplot, aes(y = trans_M, x = Y_vals)) +
        geom_point(aes(color = trans_X)) +
        geom_smooth(se = TRUE, method = "lm") +
        labs(y = "ASV abundance",
           x = paste0("log transformed \n", Y, " quantity"),
           color = paste0(X, " quantity"), title = M) +
        theme(text = element_text(size = 14))
      
      # add to plot list
      plot_list[[i]] <- p
      # save
      ggsave(paste0("/Users/alexieva/Documents/Projects/Analysis/broccoli project/02_pubAnalysis/outputs/P1_medres/P1MedRes_plot_", i, ".png"),
             plot_list[[i]],
             height = 5, width = 7, dpi = 300)
      
      
    }

    g <- ggpubr::ggarrange(plotlist = plot_list,
              ncol = 5, nrow = 4,
              common.legend = T,
              legend = "bottom")

    ggsave("/Users/alexieva/Documents/Projects/Analysis/broccoli project/02_pubAnalysis/outputs/P1MedRes_plotALL_reversed.png",
            g,
            height = 14, width = 18, dpi = 300)

## Re-run Collinsella result without Vit B1 outlier person

-   the significant combo from above, Collinsella - vitamin B1 - SFN,
    has an outlier in the vit B1 data (one inidividual just has
    extremely high B1 for some reason)
-   plan is to remove this outlier and rerun just that mediation and
    graph

<!-- -->

    med_input_B1filt <- med_input %>%
      dplyr::filter(diet_Vit_B1_mg != max(diet_Vit_B1_mg))

    med_input_B1filt2 <- med_input_B1filt %>%
            dplyr::mutate(quarts = ntile(diet_Vit_B1_mg, 4),
                          asv_trans = log(Collinsella+1),
                          sfn_trans = log(sum_SFN+1))

    model.M.B1 <- lm(asv_trans ~ quarts,
                        med_input_B1filt2)
    summary(model.M)
    ## 
    ## Call:
    ## lm(formula = asv_trans ~ quarts, data = med_input)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -1.0661 -0.4605  0.0160  0.3166  1.5162 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   4.9276     0.2604  18.922   <2e-16 ***
    ## quarts        0.1815     0.0968   1.875   0.0689 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.6664 on 36 degrees of freedom
    ## Multiple R-squared:  0.089,  Adjusted R-squared:  0.06369 
    ## F-statistic: 3.517 on 1 and 36 DF,  p-value: 0.06887

    model.Y.B1 <- lm(sfn_trans ~ quarts + asv_trans,
                        med_input_B1filt2)
    summary(model.Y)
    ## 
    ## Call:
    ## lm(formula = sfn_trans ~ quarts + asv_trans, data = med_input)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.80901 -0.48218  0.01936  0.38929  1.30601 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  5.50241    0.71547   7.691 5.02e-09 ***
    ## quarts       0.03052    0.08422   0.362    0.719    
    ## asv_trans   -0.03401    0.13840  -0.246    0.807    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.5534 on 35 degrees of freedom
    ## Multiple R-squared:  0.004328,   Adjusted R-squared:  -0.05257 
    ## F-statistic: 0.07606 on 2 and 35 DF,  p-value: 0.9269

    med_B1results <- mediate(model.M.B1, model.Y.B1, 
                             treat = "quarts", 
                             mediator = "asv_trans",
                             control.value = 1, treat.value = 4, # lowest and highest quartiles
                             boot = TRUE, sims = 500)


    modsumB1 <- summary(med_B1results)
    modsumB1
    ## 
    ## Causal Mediation Analysis 
    ## 
    ## Nonparametric Bootstrap Confidence Intervals with the Percentile Method
    ## 
    ##                Estimate 95% CI Lower 95% CI Upper p-value  
    ## ACME             -0.232       -0.626         0.01   0.064 .
    ## ADE               0.451       -0.181         1.02   0.160  
    ## Total Effect      0.219       -0.499         0.84   0.548  
    ## Prop. Mediated   -1.057       -6.299        17.55   0.596  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Sample Size Used: 37 
    ## 
    ## 
    ## Simulations: 500

    coeff_tab_medB1 <- list(ACME_est = modsumB1$d.avg,
               ACME_pval = p.adjust(modsumB1$d.avg.p, method = "fdr"),
               ADE_est = modsumB1$z.avg,
               ADE_pval = p.adjust(modsumB1$z.avg.p, method = "fdr"),
               prop_mediated_est = modsumB1$n.avg,
               prop_mediated_pval = p.adjust(modsumB1$n.avg.p, method = "fdr"))

## save image chunk

    save.image("/Users/alexieva/Documents/Projects/Analysis/broccoli project/02_pubAnalysis/Prediction1_mediation/brocc_mediation.RData")
