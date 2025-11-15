## libraries and directories

    library(dplyr)
    library(tibble)
    library(ggplot2)
    library(randomForest)
    set.seed(3)

    setwd("/Users/alexieva/Documents/Projects/Analysis/broccoli project/02_pubAnalysis/RandomForest/")

    load("/Users/alexieva/Documents/Projects/Analysis/broccoli project/02_pubAnalysis/RandomForest/brocc_randomforestSFNNIT.RData")

## Input files prep

We want an input file that has the “features” we want to test the
importance of vs the samples. The features we want there are the
microbial taxa followed by the SFN\_NIT data.

    ## input microbiome taxa table
    MB_asvtab <- read.delim("/Users/alexieva/Documents/Projects/Analysis/broccoli project/01_explAnalysis/input_data/MBasvtab_brocc_cleanT0.txt") %>%
      rownames_to_column("subject_id")
    # str(MB_asvtab)
    # MB_asvtab[,2:112] <- log(MB_asvtab[,2:112]+1)

    ## input SFN_NIT data but combine with asv tab
    SFN_NIT_tab <- read.delim("/Users/alexieva/Documents/Projects/Analysis/broccoli project/01_explAnalysis/input_data/Urine_SFN_cumul.txt")

    # no NA's?
    sum(!complete.cases(MB_asvtab))
    ## [1] 0
    sum(!complete.cases(SFN_NIT_tab))
    ## [1] 0
    # yes

    # for now test with SFN_NIT_cumulative data
    feat_tab <- MB_asvtab %>%
      inner_join(SFN_NIT_tab[c("subject_id", "sum_SFN_NIT")], by = "subject_id") %>%
      column_to_rownames("subject_id")

## graph the distribution of SFN-NIT and a test taxon

These were logged in the optimal model we used in the mediation analysis
so we intend to do that here for the random forest regressions as well.

    hist(feat_tab$sum_SFN_NIT) # non-normal

<img src="brocc_randomforest_SFNNIT_files/figure-markdown_strict/SFN-NIT and test taxon graphs-1.png" width="98%" height="98%" />

    hist(feat_tab$Dialister) # also non-normal

<img src="brocc_randomforest_SFNNIT_files/figure-markdown_strict/SFN-NIT and test taxon graphs-2.png" width="98%" height="98%" />


    hist(log(feat_tab$sum_SFN_NIT)) # more normal

<img src="brocc_randomforest_SFNNIT_files/figure-markdown_strict/SFN-NIT and test taxon graphs-3.png" width="98%" height="98%" />

    hist(log(feat_tab$Dialister)) # also more normal

<img src="brocc_randomforest_SFNNIT_files/figure-markdown_strict/SFN-NIT and test taxon graphs-4.png" width="98%" height="98%" />


    # let's just log the whole matrix (with pseudocount), as we did in the mediation analysis
    feat_tabl <- log(feat_tab + 1)

## Split into train and validation set of data

    # Training Set : Validation Set = 70 : 30 (random)
    train <- sample(nrow(feat_tabl), 0.7*nrow(feat_tabl), replace = FALSE)
    TrainSet <- feat_tabl[train,]
    ValidSet <- feat_tabl[-train,]
    summary(TrainSet)
    ##  Adlercreutzia     Agathobacter    Akkermansia      Alistipes    
    ##  Min.   :0.0000   Min.   :4.727   Min.   :0.000   Min.   :3.332  
    ##  1st Qu.:0.0000   1st Qu.:6.208   1st Qu.:0.000   1st Qu.:5.460  
    ##  Median :0.6931   Median :6.957   Median :0.000   Median :6.299  
    ##  Mean   :0.9238   Mean   :6.759   Mean   :1.479   Mean   :6.081  
    ##  3rd Qu.:1.5537   3rd Qu.:7.472   3rd Qu.:3.370   3rd Qu.:6.882  
    ##  Max.   :3.1781   Max.   :8.064   Max.   :6.136   Max.   :8.650  
    ##   Anaerostipes   Angelakisella     Bacteroides    Bifidobacterium
    ##  Min.   :2.303   Min.   :0.0000   Min.   :4.511   Min.   :0.000  
    ##  1st Qu.:5.522   1st Qu.:0.0000   1st Qu.:7.855   1st Qu.:5.087  
    ##  Median :5.948   Median :0.3466   Median :8.421   Median :5.768  
    ##  Mean   :5.893   Mean   :0.9346   Mean   :8.299   Mean   :5.622  
    ##  3rd Qu.:6.718   3rd Qu.:1.7462   3rd Qu.:9.057   3rd Qu.:7.281  
    ##  Max.   :7.743   Max.   :3.4012   Max.   :9.726   Max.   :8.241  
    ##    Bilophila        Blautia      Butyricicoccus   Butyricimonas   
    ##  Min.   :0.000   Min.   :6.639   Min.   :0.6931   Min.   :0.0000  
    ##  1st Qu.:1.727   1st Qu.:7.433   1st Qu.:4.2377   1st Qu.:0.0000  
    ##  Median :4.430   Median :7.583   Median :4.7791   Median :0.0000  
    ##  Mean   :3.637   Mean   :7.636   Mean   :4.4669   Mean   :0.6479  
    ##  3rd Qu.:5.100   3rd Qu.:7.875   3rd Qu.:4.9853   3rd Qu.:0.0000  
    ##  Max.   :6.526   Max.   :8.558   Max.   :5.7900   Max.   :3.4340  
    ##      CAG.56      Christensenellaceae.R.7.group Clostridium.sensu.stricto.1
    ##  Min.   :0.000   Min.   :0.000                 Min.   :0.0000             
    ##  1st Qu.:0.000   1st Qu.:2.136                 1st Qu.:0.1733             
    ##  Median :3.351   Median :3.658                 Median :2.8332             
    ##  Mean   :2.727   Mean   :3.810                 Mean   :2.5185             
    ##  3rd Qu.:5.074   3rd Qu.:5.804                 3rd Qu.:3.9167             
    ##  Max.   :5.808   Max.   :8.117                 Max.   :5.5984             
    ##  Colidextribacter  Collinsella     Coprococcus        DTU089      
    ##  Min.   :0.000    Min.   :0.000   Min.   :3.135   Min.   :0.0000  
    ##  1st Qu.:2.969    1st Qu.:3.330   1st Qu.:4.592   1st Qu.:0.0000  
    ##  Median :3.384    Median :3.951   Median :4.969   Median :0.6931  
    ##  Mean   :3.167    Mean   :3.708   Mean   :5.049   Mean   :0.8251  
    ##  3rd Qu.:3.749    3rd Qu.:4.840   3rd Qu.:5.721   3rd Qu.:1.3144  
    ##  Max.   :4.248    Max.   :5.778   Max.   :6.856   Max.   :2.6391  
    ##    Dialister         Dorea       Erysipelatoclostridium
    ##  Min.   :0.000   Min.   :3.045   Min.   :0.0000        
    ##  1st Qu.:0.000   1st Qu.:5.169   1st Qu.:0.0000        
    ##  Median :0.000   Median :5.622   Median :0.8959        
    ##  Mean   :1.607   Mean   :5.384   Mean   :1.3833        
    ##  3rd Qu.:3.984   3rd Qu.:5.967   3rd Qu.:2.7007        
    ##  Max.   :6.265   Max.   :6.211   Max.   :4.3944        
    ##  Erysipelotrichaceae.UCG.003 Escherichia.Shigella Faecalibacterium
    ##  Min.   :0.000               Min.   :0.0000       Min.   :5.971   
    ##  1st Qu.:3.377               1st Qu.:0.1733       1st Qu.:6.735   
    ##  Median :4.379               Median :1.4979       Median :7.460   
    ##  Mean   :4.265               Mean   :1.6104       Mean   :7.285   
    ##  3rd Qu.:5.712               3rd Qu.:2.4393       3rd Qu.:7.785   
    ##  Max.   :6.404               Max.   :4.8978       Max.   :8.103   
    ##  Family.XIII.AD3011.group Family.XIII.UCG.001 Flavonifractor  Fusicatenibacter
    ##  Min.   :0.0000           Min.   :0.0000      Min.   :0.000   Min.   :5.024   
    ##  1st Qu.:0.6931           1st Qu.:0.0000      1st Qu.:2.197   1st Qu.:5.579   
    ##  Median :1.6094           Median :0.0000      Median :3.061   Median :6.174   
    ##  Mean   :1.4480           Mean   :0.7152      Mean   :2.811   Mean   :6.176   
    ##  3rd Qu.:1.9459           3rd Qu.:1.6185      3rd Qu.:3.555   3rd Qu.:6.787   
    ##  Max.   :3.8712           Max.   :2.9444      Max.   :5.493   Max.   :7.435   
    ##  GCA.900066575    Haemophilus       Holdemania     Incertae.Sedis 
    ##  Min.   :0.000   Min.   :0.0000   Min.   :0.0000   Min.   :0.000  
    ##  1st Qu.:0.000   1st Qu.:0.1733   1st Qu.:0.0000   1st Qu.:2.404  
    ##  Median :1.386   Median :1.0397   Median :0.8959   Median :3.132  
    ##  Mean   :1.084   Mean   :1.5599   Mean   :0.9888   Mean   :3.125  
    ##  3rd Qu.:1.746   3rd Qu.:2.1678   3rd Qu.:1.6094   3rd Qu.:3.951  
    ##  Max.   :3.178   Max.   :5.5607   Max.   :2.7081   Max.   :5.442  
    ##  Intestinibacter Intestinimonas  Lachnoclostridium  Lachnospira   
    ##  Min.   :0.000   Min.   :0.000   Min.   :3.219     Min.   :0.000  
    ##  1st Qu.:1.488   1st Qu.:1.099   1st Qu.:4.435     1st Qu.:3.238  
    ##  Median :2.833   Median :1.869   Median :5.074     Median :4.848  
    ##  Mean   :2.886   Mean   :1.792   Mean   :5.130     Mean   :4.312  
    ##  3rd Qu.:4.540   3rd Qu.:2.555   3rd Qu.:5.772     3rd Qu.:5.324  
    ##  Max.   :5.394   Max.   :4.663   Max.   :6.494     Max.   :7.545  
    ##  Lachnospiraceae.AC2044.group Lachnospiraceae.FCS020.group
    ##  Min.   :0.0000               Min.   :0.000               
    ##  1st Qu.:0.0000               1st Qu.:1.226               
    ##  Median :0.0000               Median :2.565               
    ##  Mean   :0.7793               Mean   :2.161               
    ##  3rd Qu.:0.9972               3rd Qu.:3.175               
    ##  Max.   :5.1120               Max.   :4.277               
    ##  Lachnospiraceae.ND3007.group Lachnospiraceae.NK4A136.group
    ##  Min.   :0.000                Min.   :2.079                
    ##  1st Qu.:2.708                1st Qu.:3.314                
    ##  Median :4.625                Median :3.942                
    ##  Mean   :3.822                Mean   :4.051                
    ##  3rd Qu.:5.240                3rd Qu.:4.918                
    ##  Max.   :6.583                Max.   :5.861                
    ##  Lachnospiraceae.UCG.001 Lachnospiraceae.UCG.004 Lachnospiraceae.UCG.009
    ##  Min.   :0.0000          Min.   :0.000           Min.   :0.0000         
    ##  1st Qu.:0.0000          1st Qu.:3.628           1st Qu.:0.6931         
    ##  Median :0.8959          Median :4.412           Median :1.0986         
    ##  Mean   :1.6896          Mean   :4.215           Mean   :0.9783         
    ##  3rd Qu.:3.1945          3rd Qu.:5.019           3rd Qu.:1.3863         
    ##  Max.   :5.7746          Max.   :5.717           Max.   :2.1972         
    ##  Lachnospiraceae.UCG.010  Lachnotalea     Marvinbryantia  Methanobrevibacter
    ##  Min.   :0.000           Min.   :0.0000   Min.   :0.000   Min.   :0.0000    
    ##  1st Qu.:1.272           1st Qu.:0.0000   1st Qu.:1.655   1st Qu.:0.0000    
    ##  Median :2.817           Median :0.0000   Median :2.967   Median :0.0000    
    ##  Mean   :2.478           Mean   :0.5185   Mean   :2.643   Mean   :0.8807    
    ##  3rd Qu.:3.604           3rd Qu.:0.0000   3rd Qu.:3.575   3rd Qu.:0.0000    
    ##  Max.   :4.836           Max.   :4.2627   Max.   :5.501   Max.   :6.8222    
    ##    Monoglobus    NK4A214.group     Odoribacter     Oscillibacter  
    ##  Min.   :1.386   Min.   :0.0000   Min.   :0.0000   Min.   :0.000  
    ##  1st Qu.:3.495   1st Qu.:0.4479   1st Qu.:0.1733   1st Qu.:2.135  
    ##  Median :4.174   Median :3.2514   Median :2.9907   Median :2.917  
    ##  Mean   :3.915   Mean   :2.7603   Mean   :2.2979   Mean   :2.798  
    ##  3rd Qu.:4.655   3rd Qu.:4.3625   3rd Qu.:3.4411   3rd Qu.:3.526  
    ##  Max.   :5.561   Max.   :5.9584   Max.   :4.4886   Max.   :5.298  
    ##   Oscillospira     Oxalobacter     Parabacteroides Parasutterella  
    ##  Min.   :0.0000   Min.   :0.0000   Min.   :0.000   Min.   :0.0000  
    ##  1st Qu.:0.1733   1st Qu.:0.0000   1st Qu.:4.896   1st Qu.:0.4479  
    ##  Median :1.2425   Median :0.0000   Median :5.515   Median :2.7917  
    ##  Mean   :1.1665   Mean   :0.4746   Mean   :5.499   Mean   :2.9764  
    ##  3rd Qu.:1.7918   3rd Qu.:0.0000   3rd Qu.:6.240   3rd Qu.:5.0998  
    ##  Max.   :3.2189   Max.   :3.4340   Max.   :7.811   Max.   :7.7579  
    ##  Phascolarctobacterium  Prevotella_9      Romboutsia      Roseburia    
    ##  Min.   :0.0000        Min.   :0.0000   Min.   :0.000   Min.   :4.143  
    ##  1st Qu.:0.0000        1st Qu.:0.6931   1st Qu.:2.404   1st Qu.:4.892  
    ##  Median :0.3466        Median :1.3863   Median :4.254   Median :5.598  
    ##  Mean   :1.8145        Mean   :2.9777   Mean   :3.866   Mean   :5.406  
    ##  3rd Qu.:4.1267        3rd Qu.:5.5658   3rd Qu.:5.263   3rd Qu.:5.776  
    ##  Max.   :5.1761        Max.   :9.3712   Max.   :6.812   Max.   :6.988  
    ##   Ruminococcus   Shuttleworthia   Streptococcus   Subdoligranulum
    ##  Min.   :0.000   Min.   :0.0000   Min.   :0.000   Min.   :3.807  
    ##  1st Qu.:1.884   1st Qu.:0.0000   1st Qu.:2.109   1st Qu.:5.630  
    ##  Median :5.058   Median :0.0000   Median :2.826   Median :6.251  
    ##  Mean   :4.372   Mean   :0.5845   Mean   :2.879   Mean   :6.251  
    ##  3rd Qu.:6.417   3rd Qu.:0.0000   3rd Qu.:3.611   3rd Qu.:6.655  
    ##  Max.   :7.772   Max.   :4.1271   Max.   :4.852   Max.   :8.129  
    ##    Sutterella     Terrisporobacter  Turicibacter       UBA1819      
    ##  Min.   :0.0000   Min.   :0.0000   Min.   :0.0000   Min.   :0.0000  
    ##  1st Qu.:0.0000   1st Qu.:0.0000   1st Qu.:0.1733   1st Qu.:0.7945  
    ##  Median :0.6931   Median :0.0000   Median :1.8444   Median :1.6094  
    ##  Mean   :2.4621   Mean   :0.9604   Mean   :1.8441   Mean   :1.6146  
    ##  3rd Qu.:5.5328   3rd Qu.:0.6931   3rd Qu.:2.8332   3rd Qu.:2.1678  
    ##  Max.   :6.8701   Max.   :5.3279   Max.   :4.3041   Max.   :4.0943  
    ##     UCG.002         UCG.003          UCG.005         UCG.009      
    ##  Min.   :0.000   Min.   :0.0000   Min.   :0.000   Min.   :0.0000  
    ##  1st Qu.:1.609   1st Qu.:0.1733   1st Qu.:1.488   1st Qu.:0.0000  
    ##  Median :5.053   Median :3.2157   Median :3.927   Median :0.3466  
    ##  Mean   :3.978   Mean   :2.6366   Mean   :3.225   Mean   :0.8839  
    ##  3rd Qu.:5.821   3rd Qu.:4.2902   3rd Qu.:4.993   3rd Qu.:1.7462  
    ##  Max.   :7.003   Max.   :5.6204   Max.   :6.475   Max.   :3.2958  
    ##  X.Eubacterium..eligens.group X.Eubacterium..fissicatena.group
    ##  Min.   :0.000                Min.   :0.0000                  
    ##  1st Qu.:1.310                1st Qu.:0.0000                  
    ##  Median :3.377                Median :0.3466                  
    ##  Mean   :2.916                Mean   :0.6935                  
    ##  3rd Qu.:4.134                3rd Qu.:1.3863                  
    ##  Max.   :5.631                Max.   :2.0794                  
    ##  X.Eubacterium..hallii.group X.Eubacterium..ruminantium.group
    ##  Min.   :0.000               Min.   :0.0000                  
    ##  1st Qu.:5.439               1st Qu.:0.0000                  
    ##  Median :5.644               Median :0.3466                  
    ##  Mean   :5.228               Mean   :1.4806                  
    ##  3rd Qu.:5.852               3rd Qu.:2.9533                  
    ##  Max.   :6.628               Max.   :6.4489                  
    ##  X.Eubacterium..siraeum.group X.Eubacterium..ventriosum.group
    ##  Min.   :0.0000               Min.   :0.000                  
    ##  1st Qu.:0.0000               1st Qu.:2.119                  
    ##  Median :0.6931               Median :3.151                  
    ##  Mean   :1.9233               Mean   :3.025                  
    ##  3rd Qu.:3.6604               3rd Qu.:4.247                  
    ##  Max.   :6.0730               Max.   :5.598                  
    ##  X.Eubacterium..xylanophilum.group X.Ruminococcus..gauvreauii.group
    ##  Min.   :0.000                     Min.   :0.000                   
    ##  1st Qu.:0.000                     1st Qu.:3.023                   
    ##  Median :1.386                     Median :4.418                   
    ##  Mean   :1.778                     Mean   :3.654                   
    ##  3rd Qu.:3.624                     3rd Qu.:4.874                   
    ##  Max.   :4.745                     Max.   :5.930                   
    ##  X.Ruminococcus..gnavus.group X.Ruminococcus..torques.group
    ##  Min.   :0.0000               Min.   :2.079                
    ##  1st Qu.:0.1733               1st Qu.:4.556                
    ##  Median :1.5890               Median :5.447                
    ##  Mean   :1.8431               Mean   :5.100                
    ##  3rd Qu.:2.9953               3rd Qu.:5.850                
    ##  Max.   :5.4638               Max.   :6.756                
    ##  f_Christensenellaceae_ASV181 f_Coriobacteriales.Incertae.Sedis_ASV433
    ##  Min.   :0.0000               Min.   :0.000                           
    ##  1st Qu.:0.0000               1st Qu.:0.000                           
    ##  Median :0.0000               Median :0.000                           
    ##  Mean   :0.4426               Mean   :0.896                           
    ##  3rd Qu.:0.0000               3rd Qu.:1.482                           
    ##  Max.   :3.8712               Max.   :3.689                           
    ##  f_Eggerthellaceae_ASV319 f_Erysipelotrichaceae_ASV487 f_Lachnospiraceae_ASV117
    ##  Min.   :0.000            Min.   :0.0000               Min.   :0.0000          
    ##  1st Qu.:0.000            1st Qu.:0.0000               1st Qu.:0.0000          
    ##  Median :0.000            Median :0.0000               Median :0.6931          
    ##  Mean   :1.109            Mean   :0.5469               Mean   :0.9825          
    ##  3rd Qu.:2.545            3rd Qu.:0.0000               3rd Qu.:1.0986          
    ##  Max.   :5.088            Max.   :3.9703               Max.   :4.2195          
    ##  f_Lachnospiraceae_ASV239 f_Lachnospiraceae_ASV259 f_Lachnospiraceae_ASV275
    ##  Min.   :0.000            Min.   :0.0000           Min.   :0.000           
    ##  1st Qu.:1.655            1st Qu.:0.7945           1st Qu.:0.000           
    ##  Median :2.197            Median :2.0127           Median :1.792           
    ##  Mean   :2.162            Mean   :1.8018           Mean   :1.383           
    ##  3rd Qu.:2.818            3rd Qu.:2.6723           3rd Qu.:2.079           
    ##  Max.   :3.892            Max.   :3.3322           Max.   :3.807           
    ##  f_Lachnospiraceae_ASV283 f_Lachnospiraceae_ASV297 f_Lachnospiraceae_ASV326
    ##  Min.   :0.000            Min.   :0.0000           Min.   :0.000           
    ##  1st Qu.:0.000            1st Qu.:0.0000           1st Qu.:0.000           
    ##  Median :1.701            Median :0.0000           Median :1.869           
    ##  Mean   :1.400            Mean   :0.8688           Mean   :1.445           
    ##  3rd Qu.:2.303            3rd Qu.:1.8060           3rd Qu.:2.303           
    ##  Max.   :3.045            Max.   :3.8286           Max.   :3.296           
    ##  f_Lachnospiraceae_ASV384 f_Lachnospiraceae_ASV398 f_Lachnospiraceae_ASV410
    ##  Min.   :0.0000           Min.   :0.0000           Min.   :0.0000          
    ##  1st Qu.:0.6931           1st Qu.:0.0000           1st Qu.:0.0000          
    ##  Median :1.3863           Median :0.0000           Median :0.8959          
    ##  Mean   :1.3149           Mean   :0.7965           Mean   :1.1734          
    ##  3rd Qu.:1.9074           3rd Qu.:1.3144           3rd Qu.:2.0794          
    ##  Max.   :2.8904           Max.   :3.2581           Max.   :3.4340          
    ##  f_Lachnospiraceae_ASV514 f_Lachnospiraceae_ASV532 f_Lachnospiraceae_ASV555
    ##  Min.   :0.0000           Min.   :0.0000           Min.   :0.0000          
    ##  1st Qu.:0.0000           1st Qu.:0.0000           1st Qu.:0.0000          
    ##  Median :0.0000           Median :0.0000           Median :0.0000          
    ##  Mean   :0.5708           Mean   :0.7746           Mean   :0.5638          
    ##  3rd Qu.:0.5199           3rd Qu.:1.7462           3rd Qu.:0.6931          
    ##  Max.   :3.9512           Max.   :4.1109           Max.   :2.6391          
    ##  f_Lachnospiraceae_ASV558  f_NA_ASV128      f_NA_ASV156      f_NA_ASV184   
    ##  Min.   :0.0000           Min.   :0.0000   Min.   :0.0000   Min.   :0.000  
    ##  1st Qu.:0.0000           1st Qu.:0.0000   1st Qu.:0.0000   1st Qu.:0.000  
    ##  Median :0.3466           Median :0.0000   Median :0.0000   Median :0.000  
    ##  Mean   :0.8006           Mean   :0.7338   Mean   :0.6992   Mean   :1.182  
    ##  3rd Qu.:1.5537           3rd Qu.:0.0000   3rd Qu.:0.0000   3rd Qu.:3.045  
    ##  Max.   :2.8904           Max.   :6.5971   Max.   :4.7622   Max.   :5.198  
    ##   f_NA_ASV516      f_NA_ASV691    f_Oscillospiraceae_ASV141
    ##  Min.   :0.0000   Min.   :0.000   Min.   :0.000            
    ##  1st Qu.:0.0000   1st Qu.:0.000   1st Qu.:0.000            
    ##  Median :0.0000   Median :1.040   Median :3.387            
    ##  Mean   :0.5692   Mean   :1.059   Mean   :2.565            
    ##  3rd Qu.:0.9972   3rd Qu.:2.046   3rd Qu.:3.871            
    ##  Max.   :2.8332   Max.   :3.045   Max.   :5.293            
    ##  f_Ruminococcaceae_ASV327 f_Ruminococcaceae_ASV610 f_UCG.010_ASV245
    ##  Min.   :0.000            Min.   :0.0000           Min.   :0.0000  
    ##  1st Qu.:0.000            1st Qu.:0.0000           1st Qu.:0.0000  
    ##  Median :2.047            Median :0.0000           Median :0.0000  
    ##  Mean   :1.721            Mean   :0.8468           Mean   :0.8423  
    ##  3rd Qu.:2.983            3rd Qu.:1.6904           3rd Qu.:0.6931  
    ##  Max.   :3.584            Max.   :3.2958           Max.   :6.3936  
    ##  f_UCG.010_ASV667 f_.Eubacterium..coprostanoligenes.group_ASV203
    ##  Min.   :0.0000   Min.   :0.000                                 
    ##  1st Qu.:0.0000   1st Qu.:0.000                                 
    ##  Median :0.0000   Median :0.000                                 
    ##  Mean   :0.5682   Mean   :1.206                                 
    ##  3rd Qu.:0.8240   3rd Qu.:2.992                                 
    ##  Max.   :2.9444   Max.   :5.533                                 
    ##   sum_SFN_NIT   
    ##  Min.   :1.805  
    ##  1st Qu.:2.865  
    ##  Median :3.363  
    ##  Mean   :3.325  
    ##  3rd Qu.:3.824  
    ##  Max.   :4.785
    summary(ValidSet)
    ##  Adlercreutzia     Agathobacter    Akkermansia      Alistipes    
    ##  Min.   :0.0000   Min.   :3.497   Min.   :0.000   Min.   :0.000  
    ##  1st Qu.:0.5199   1st Qu.:6.111   1st Qu.:0.000   1st Qu.:5.878  
    ##  Median :1.0986   Median :6.821   Median :0.000   Median :6.223  
    ##  Mean   :1.2698   Mean   :6.466   Mean   :1.713   Mean   :5.926  
    ##  3rd Qu.:2.1352   3rd Qu.:7.058   3rd Qu.:1.895   3rd Qu.:6.957  
    ##  Max.   :2.7726   Max.   :7.923   Max.   :7.003   Max.   :7.773  
    ##   Anaerostipes   Angelakisella     Bacteroides    Bifidobacterium
    ##  Min.   :5.056   Min.   :0.0000   Min.   :6.395   Min.   :0.000  
    ##  1st Qu.:5.377   1st Qu.:0.0000   1st Qu.:7.696   1st Qu.:2.472  
    ##  Median :5.557   Median :0.0000   Median :8.362   Median :4.856  
    ##  Mean   :5.580   Mean   :0.2834   Mean   :8.163   Mean   :4.254  
    ##  3rd Qu.:5.636   3rd Qu.:0.0000   3rd Qu.:8.711   3rd Qu.:6.113  
    ##  Max.   :6.592   Max.   :1.7918   Max.   :9.079   Max.   :8.548  
    ##    Bilophila        Blautia      Butyricicoccus  Butyricimonas   
    ##  Min.   :0.000   Min.   :6.751   Min.   :0.000   Min.   :0.0000  
    ##  1st Qu.:3.767   1st Qu.:7.454   1st Qu.:4.164   1st Qu.:0.0000  
    ##  Median :4.430   Median :7.819   Median :4.663   Median :0.0000  
    ##  Mean   :3.981   Mean   :7.772   Mean   :4.320   Mean   :0.2071  
    ##  3rd Qu.:5.129   3rd Qu.:8.068   3rd Qu.:5.089   3rd Qu.:0.0000  
    ##  Max.   :5.252   Max.   :8.504   Max.   :5.796   Max.   :1.7918  
    ##      CAG.56      Christensenellaceae.R.7.group Clostridium.sensu.stricto.1
    ##  Min.   :0.000   Min.   :0.000                 Min.   :0.000              
    ##  1st Qu.:1.040   1st Qu.:2.649                 1st Qu.:0.000              
    ##  Median :4.429   Median :4.712                 Median :1.498              
    ##  Mean   :3.408   Mean   :4.153                 Mean   :2.276              
    ##  3rd Qu.:5.152   3rd Qu.:5.642                 3rd Qu.:4.329              
    ##  Max.   :6.122   Max.   :7.426                 Max.   :5.580              
    ##  Colidextribacter  Collinsella     Coprococcus        DTU089     
    ##  Min.   :0.000    Min.   :0.000   Min.   :3.761   Min.   :0.000  
    ##  1st Qu.:3.005    1st Qu.:2.599   1st Qu.:4.822   1st Qu.:0.000  
    ##  Median :3.584    Median :4.297   Median :5.531   Median :1.099  
    ##  Mean   :3.305    Mean   :3.369   Mean   :5.437   Mean   :1.149  
    ##  3rd Qu.:4.047    3rd Qu.:4.686   3rd Qu.:6.171   3rd Qu.:2.109  
    ##  Max.   :4.754    Max.   :5.283   Max.   :6.829   Max.   :2.708  
    ##    Dialister         Dorea       Erysipelatoclostridium
    ##  Min.   :0.000   Min.   :3.258   Min.   :0.0000        
    ##  1st Qu.:0.000   1st Qu.:4.700   1st Qu.:0.0000        
    ##  Median :3.332   Median :5.145   Median :0.0000        
    ##  Mean   :2.500   Mean   :4.984   Mean   :0.6885        
    ##  3rd Qu.:4.615   3rd Qu.:5.399   3rd Qu.:1.5596        
    ##  Max.   :4.905   Max.   :5.759   Max.   :2.3979        
    ##  Erysipelotrichaceae.UCG.003 Escherichia.Shigella Faecalibacterium
    ##  Min.   :0.000               Min.   :0.0000       Min.   :6.122   
    ##  1st Qu.:3.643               1st Qu.:0.0000       1st Qu.:7.167   
    ##  Median :4.376               Median :0.6931       Median :7.415   
    ##  Mean   :4.043               Mean   :1.0782       Mean   :7.381   
    ##  3rd Qu.:5.144               3rd Qu.:1.8303       3rd Qu.:7.821   
    ##  Max.   :6.260               Max.   :3.8067       Max.   :8.396   
    ##  Family.XIII.AD3011.group Family.XIII.UCG.001 Flavonifractor  Fusicatenibacter
    ##  Min.   :0.000            Min.   :0.0000      Min.   :0.000   Min.   :5.075   
    ##  1st Qu.:0.000            1st Qu.:0.0000      1st Qu.:1.099   1st Qu.:6.130   
    ##  Median :2.398            Median :0.0000      Median :1.869   Median :6.241   
    ##  Mean   :1.866            Mean   :0.6246      Mean   :2.174   Mean   :6.340   
    ##  3rd Qu.:2.846            3rd Qu.:0.4479      3rd Qu.:3.324   3rd Qu.:6.523   
    ##  Max.   :4.431            Max.   :2.9957      Max.   :4.043   Max.   :7.340   
    ##  GCA.900066575    Haemophilus      Holdemania     Incertae.Sedis 
    ##  Min.   :0.000   Min.   :0.000   Min.   :0.0000   Min.   :0.000  
    ##  1st Qu.:0.000   1st Qu.:0.000   1st Qu.:0.0000   1st Qu.:2.969  
    ##  Median :1.994   Median :1.151   Median :0.6931   Median :3.596  
    ##  Mean   :1.531   Mean   :1.760   Mean   :0.9286   Mean   :3.798  
    ##  3rd Qu.:2.348   3rd Qu.:3.156   3rd Qu.:1.6550   3rd Qu.:4.967  
    ##  Max.   :3.367   Max.   :5.094   Max.   :2.4849   Max.   :6.413  
    ##  Intestinibacter Intestinimonas  Lachnoclostridium  Lachnospira   
    ##  Min.   :0.000   Min.   :0.000   Min.   :3.664     Min.   :0.000  
    ##  1st Qu.:1.380   1st Qu.:1.554   1st Qu.:4.608     1st Qu.:4.225  
    ##  Median :3.003   Median :2.250   Median :5.225     Median :4.627  
    ##  Mean   :2.743   Mean   :2.111   Mean   :5.083     Mean   :4.353  
    ##  3rd Qu.:3.942   3rd Qu.:2.821   3rd Qu.:5.524     3rd Qu.:5.024  
    ##  Max.   :5.198   Max.   :3.829   Max.   :6.223     Max.   :6.230  
    ##  Lachnospiraceae.AC2044.group Lachnospiraceae.FCS020.group
    ##  Min.   :0.0000               Min.   :0.000               
    ##  1st Qu.:0.0000               1st Qu.:1.862               
    ##  Median :0.0000               Median :2.453               
    ##  Mean   :0.8633               Mean   :2.403               
    ##  3rd Qu.:1.5596               3rd Qu.:3.229               
    ##  Max.   :4.0604               Max.   :4.277               
    ##  Lachnospiraceae.ND3007.group Lachnospiraceae.NK4A136.group
    ##  Min.   :0.000                Min.   :3.045                
    ##  1st Qu.:4.266                1st Qu.:3.604                
    ##  Median :4.958                Median :4.454                
    ##  Mean   :4.495                Mean   :4.654                
    ##  3rd Qu.:5.129                3rd Qu.:5.510                
    ##  Max.   :5.796                Max.   :6.547                
    ##  Lachnospiraceae.UCG.001 Lachnospiraceae.UCG.004 Lachnospiraceae.UCG.009
    ##  Min.   :0.000           Min.   :1.386           Min.   :0.0000         
    ##  1st Qu.:0.000           1st Qu.:3.734           1st Qu.:0.0000         
    ##  Median :2.674           Median :4.388           Median :0.8959         
    ##  Mean   :2.430           Mean   :4.025           Mean   :1.0170         
    ##  3rd Qu.:3.829           3rd Qu.:4.582           3rd Qu.:1.4421         
    ##  Max.   :5.663           Max.   :5.460           Max.   :3.2581         
    ##  Lachnospiraceae.UCG.010  Lachnotalea     Marvinbryantia  Methanobrevibacter
    ##  Min.   :0.0000          Min.   :0.0000   Min.   :1.386   Min.   :0.000     
    ##  1st Qu.:0.9972          1st Qu.:0.0000   1st Qu.:2.652   1st Qu.:0.000     
    ##  Median :2.4452          Median :0.0000   Median :3.219   Median :0.000     
    ##  Mean   :2.2719          Mean   :0.8438   Mean   :3.243   Mean   :1.560     
    ##  3rd Qu.:3.2957          3rd Qu.:1.8637   3rd Qu.:3.707   3rd Qu.:3.814     
    ##  Max.   :4.8828          Max.   :2.5649   Max.   :5.100   Max.   :5.442     
    ##    Monoglobus    NK4A214.group    Odoribacter    Oscillibacter  
    ##  Min.   :0.000   Min.   :0.000   Min.   :0.000   Min.   :0.000  
    ##  1st Qu.:3.353   1st Qu.:0.000   1st Qu.:3.054   1st Qu.:2.348  
    ##  Median :4.201   Median :3.290   Median :3.350   Median :3.043  
    ##  Mean   :3.932   Mean   :2.309   Mean   :2.839   Mean   :2.881  
    ##  3rd Qu.:4.840   3rd Qu.:3.763   3rd Qu.:3.495   3rd Qu.:3.922  
    ##  Max.   :6.752   Max.   :5.481   Max.   :3.871   Max.   :4.331  
    ##   Oscillospira     Oxalobacter     Parabacteroides Parasutterella 
    ##  Min.   :0.0000   Min.   :0.0000   Min.   :3.091   Min.   :0.000  
    ##  1st Qu.:0.0000   1st Qu.:0.0000   1st Qu.:4.961   1st Qu.:3.404  
    ##  Median :0.0000   Median :0.0000   Median :5.954   Median :5.446  
    ##  Mean   :0.8187   Mean   :0.7723   Mean   :5.576   Mean   :4.493  
    ##  3rd Qu.:1.6550   3rd Qu.:1.9459   3rd Qu.:6.331   3rd Qu.:5.952  
    ##  Max.   :3.0910   Max.   :2.8904   Max.   :6.590   Max.   :6.304  
    ##  Phascolarctobacterium  Prevotella_9      Romboutsia      Roseburia    
    ##  Min.   :0.000         Min.   :0.0000   Min.   :0.000   Min.   :4.043  
    ##  1st Qu.:0.000         1st Qu.:0.0000   1st Qu.:2.973   1st Qu.:4.921  
    ##  Median :0.973         Median :0.6931   Median :3.839   Median :5.258  
    ##  Mean   :2.051         Mean   :1.9464   Mean   :3.410   Mean   :5.297  
    ##  3rd Qu.:4.137         3rd Qu.:0.8664   3rd Qu.:4.278   3rd Qu.:5.743  
    ##  Max.   :5.011         Max.   :9.7906   Max.   :6.207   Max.   :6.430  
    ##   Ruminococcus    Shuttleworthia   Streptococcus   Subdoligranulum
    ##  Min.   :0.6931   Min.   :0.0000   Min.   :1.099   Min.   :4.970  
    ##  1st Qu.:2.4830   1st Qu.:0.0000   1st Qu.:2.655   1st Qu.:5.882  
    ##  Median :5.3523   Median :0.0000   Median :4.047   Median :6.305  
    ##  Mean   :4.5392   Mean   :0.4479   Mean   :3.522   Mean   :6.157  
    ##  3rd Qu.:6.3773   3rd Qu.:0.0000   3rd Qu.:4.374   3rd Qu.:6.515  
    ##  Max.   :7.7129   Max.   :3.5835   Max.   :4.997   Max.   :7.031  
    ##    Sutterella    Terrisporobacter  Turicibacter       UBA1819      
    ##  Min.   :0.000   Min.   :0.0000   Min.   :0.0000   Min.   :0.0000  
    ##  1st Qu.:0.000   1st Qu.:0.0000   1st Qu.:0.5199   1st Qu.:0.9972  
    ##  Median :0.000   Median :0.3466   Median :3.2551   Median :1.6094  
    ##  Mean   :1.547   Mean   :0.8316   Mean   :2.5590   Mean   :1.4551  
    ##  3rd Qu.:2.971   3rd Qu.:1.0063   3rd Qu.:3.8967   3rd Qu.:2.0087  
    ##  Max.   :5.919   Max.   :4.0073   Max.   :4.9273   Max.   :2.5649  
    ##     UCG.002         UCG.003         UCG.005         UCG.009      
    ##  Min.   :0.000   Min.   :0.000   Min.   :0.000   Min.   :0.0000  
    ##  1st Qu.:4.596   1st Qu.:3.244   1st Qu.:2.808   1st Qu.:0.0000  
    ##  Median :5.178   Median :4.282   Median :3.543   Median :0.0000  
    ##  Mean   :4.511   Mean   :3.482   Mean   :3.582   Mean   :0.2409  
    ##  3rd Qu.:5.926   3rd Qu.:4.527   3rd Qu.:5.347   3rd Qu.:0.1733  
    ##  Max.   :6.737   Max.   :5.050   Max.   :6.292   Max.   :1.0986  
    ##  X.Eubacterium..eligens.group X.Eubacterium..fissicatena.group
    ##  Min.   :0.6931               Min.   :0.0000                  
    ##  1st Qu.:2.4632               1st Qu.:0.0000                  
    ##  Median :3.5488               Median :0.8959                  
    ##  Mean   :3.5124               Mean   :1.1857                  
    ##  3rd Qu.:4.6956               3rd Qu.:2.3026                  
    ##  Max.   :5.7991               Max.   :2.9957                  
    ##  X.Eubacterium..hallii.group X.Eubacterium..ruminantium.group
    ##  Min.   :0.000               Min.   :0.0000                  
    ##  1st Qu.:4.928               1st Qu.:0.0000                  
    ##  Median :5.527               Median :0.3466                  
    ##  Mean   :5.012               Mean   :2.0242                  
    ##  3rd Qu.:5.687               3rd Qu.:4.1123                  
    ##  Max.   :6.240               Max.   :6.8522                  
    ##  X.Eubacterium..siraeum.group X.Eubacterium..ventriosum.group
    ##  Min.   :0.000                Min.   :0.000                  
    ##  1st Qu.:0.000                1st Qu.:3.137                  
    ##  Median :2.359                Median :3.661                  
    ##  Mean   :2.377                Mean   :3.496                  
    ##  3rd Qu.:3.914                3rd Qu.:4.390                  
    ##  Max.   :5.749                Max.   :5.342                  
    ##  X.Eubacterium..xylanophilum.group X.Ruminococcus..gauvreauii.group
    ##  Min.   :0.000                     Min.   :0.000                   
    ##  1st Qu.:0.000                     1st Qu.:2.551                   
    ##  Median :2.975                     Median :4.602                   
    ##  Mean   :2.707                     Mean   :3.543                   
    ##  3rd Qu.:4.411                     3rd Qu.:4.818                   
    ##  Max.   :5.932                     Max.   :5.568                   
    ##  X.Ruminococcus..gnavus.group X.Ruminococcus..torques.group
    ##  Min.   :0.000                Min.   :1.099                
    ##  1st Qu.:0.000                1st Qu.:4.605                
    ##  Median :1.869                Median :5.584                
    ##  Mean   :1.783                Mean   :5.130                
    ##  3rd Qu.:3.288                3rd Qu.:6.011                
    ##  Max.   :4.205                Max.   :6.483                
    ##  f_Christensenellaceae_ASV181 f_Coriobacteriales.Incertae.Sedis_ASV433
    ##  Min.   :0.000                Min.   :0.0000                          
    ##  1st Qu.:0.000                1st Qu.:0.0000                          
    ##  Median :0.000                Median :0.3466                          
    ##  Mean   :1.336                Mean   :0.6114                          
    ##  3rd Qu.:3.400                3rd Qu.:0.7945                          
    ##  Max.   :5.273                Max.   :2.7726                          
    ##  f_Eggerthellaceae_ASV319 f_Erysipelotrichaceae_ASV487 f_Lachnospiraceae_ASV117
    ##  Min.   :0.000            Min.   :0.000                Min.   :0.000           
    ##  1st Qu.:0.000            1st Qu.:0.000                1st Qu.:0.000           
    ##  Median :0.973            Median :0.000                Median :0.000           
    ##  Mean   :1.516            Mean   :1.114                Mean   :1.334           
    ##  3rd Qu.:2.969            3rd Qu.:2.368                3rd Qu.:2.194           
    ##  Max.   :4.234            Max.   :3.296                Max.   :5.743           
    ##  f_Lachnospiraceae_ASV239 f_Lachnospiraceae_ASV259 f_Lachnospiraceae_ASV275
    ##  Min.   :0.000            Min.   :0.0000           Min.   :0.000           
    ##  1st Qu.:0.824            1st Qu.:0.5199           1st Qu.:0.000           
    ##  Median :2.191            Median :2.8315           Median :1.903           
    ##  Mean   :1.952            Mean   :2.1436           Mean   :1.407           
    ##  3rd Qu.:3.097            3rd Qu.:3.3758           3rd Qu.:2.326           
    ##  Max.   :3.850            Max.   :3.5553           Max.   :3.689           
    ##  f_Lachnospiraceae_ASV283 f_Lachnospiraceae_ASV297 f_Lachnospiraceae_ASV326
    ##  Min.   :0.000            Min.   :0.0000           Min.   :1.099           
    ##  1st Qu.:1.517            1st Qu.:0.5199           1st Qu.:1.746           
    ##  Median :2.636            Median :1.9560           Median :1.936           
    ##  Mean   :2.091            Mean   :1.6904           Mean   :2.154           
    ##  3rd Qu.:2.981            3rd Qu.:2.5234           3rd Qu.:2.476           
    ##  Max.   :3.135            Max.   :4.1589           Max.   :3.738           
    ##  f_Lachnospiraceae_ASV384 f_Lachnospiraceae_ASV398 f_Lachnospiraceae_ASV410
    ##  Min.   :0.0000           Min.   :0.0000           Min.   :0.0000          
    ##  1st Qu.:0.5199           1st Qu.:0.0000           1st Qu.:0.0000          
    ##  Median :1.3863           Median :0.0000           Median :0.6931          
    ##  Mean   :1.3874           Mean   :0.9402           Mean   :0.7584          
    ##  3rd Qu.:2.1089           3rd Qu.:1.8066           3rd Qu.:1.4421          
    ##  Max.   :3.7136           Max.   :3.6376           Max.   :2.0794          
    ##  f_Lachnospiraceae_ASV514 f_Lachnospiraceae_ASV532 f_Lachnospiraceae_ASV555
    ##  Min.   :0.0000           Min.   :0.0000           Min.   :0.0000          
    ##  1st Qu.:0.0000           1st Qu.:0.0000           1st Qu.:0.5199          
    ##  Median :0.0000           Median :0.0000           Median :1.3863          
    ##  Mean   :0.7134           Mean   :0.2137           Mean   :1.3242          
    ##  3rd Qu.:1.2263           3rd Qu.:0.0000           3rd Qu.:1.9793          
    ##  Max.   :3.3673           Max.   :2.5649           Max.   :2.9444          
    ##  f_Lachnospiraceae_ASV558  f_NA_ASV128      f_NA_ASV156      f_NA_ASV184    
    ##  Min.   :0.0000           Min.   :0.0000   Min.   :0.0000   Min.   :0.0000  
    ##  1st Qu.:0.5199           1st Qu.:0.0000   1st Qu.:0.0000   1st Qu.:0.0000  
    ##  Median :1.3540           Median :0.3466   Median :0.0000   Median :0.0000  
    ##  Mean   :1.1994           Mean   :1.5694   Mean   :0.2648   Mean   :0.9745  
    ##  3rd Qu.:1.6550           3rd Qu.:3.5318   3rd Qu.:0.1733   3rd Qu.:1.6479  
    ##  Max.   :2.4849           Max.   :5.8201   Max.   :1.7918   Max.   :3.9318  
    ##   f_NA_ASV516      f_NA_ASV691    f_Oscillospiraceae_ASV141
    ##  Min.   :0.0000   Min.   :0.000   Min.   :0.000            
    ##  1st Qu.:0.0000   1st Qu.:0.000   1st Qu.:0.000            
    ##  Median :0.0000   Median :1.386   Median :0.000            
    ##  Mean   :0.5821   Mean   :1.057   Mean   :1.712            
    ##  3rd Qu.:0.5199   3rd Qu.:1.792   3rd Qu.:4.123            
    ##  Max.   :2.7081   Max.   :2.079   Max.   :4.700            
    ##  f_Ruminococcaceae_ASV327 f_Ruminococcaceae_ASV610 f_UCG.010_ASV245
    ##  Min.   :0.0000           Min.   :0.0000           Min.   :0.000   
    ##  1st Qu.:0.0000           1st Qu.:0.0000           1st Qu.:0.000   
    ##  Median :0.0000           Median :0.8959           Median :0.000   
    ##  Mean   :0.9023           Mean   :1.0354           Mean   :1.238   
    ##  3rd Qu.:1.2810           3rd Qu.:1.7918           3rd Qu.:2.767   
    ##  Max.   :3.2189           Max.   :2.4849           Max.   :3.829   
    ##  f_UCG.010_ASV667 f_.Eubacterium..coprostanoligenes.group_ASV203
    ##  Min.   :0.0000   Min.   :0.000                                 
    ##  1st Qu.:0.0000   1st Qu.:0.000                                 
    ##  Median :0.0000   Median :0.000                                 
    ##  Mean   :0.6772   Mean   :1.230                                 
    ##  3rd Qu.:1.3104   3rd Qu.:2.673                                 
    ##  Max.   :3.1355   Max.   :4.443                                 
    ##   sum_SFN_NIT   
    ##  Min.   :2.688  
    ##  1st Qu.:3.337  
    ##  Median :3.714  
    ##  Mean   :3.699  
    ##  3rd Qu.:4.133  
    ##  Max.   :4.700

## create random forest model to test which taxa are associated with cumulative SFN\_NIT

-   default uses 500 trees
-   regression model

<!-- -->

    # Create a Random Forest model with default parameters
    model1 <- randomForest(sum_SFN_NIT ~ ., data = TrainSet, importance = TRUE)
    model1
    ## 
    ## Call:
    ##  randomForest(formula = sum_SFN_NIT ~ ., data = TrainSet, importance = TRUE) 
    ##                Type of random forest: regression
    ##                      Number of trees: 500
    ## No. of variables tried at each split: 37
    ## 
    ##           Mean of squared residuals: 0.5492153
    ##                     % Var explained: -25.01
    round(importance(model1), 2)
    ##                                                %IncMSE IncNodePurity
    ## Adlercreutzia                                     0.04          0.04
    ## Agathobacter                                      0.98          0.14
    ## Akkermansia                                      -1.48          0.01
    ## Alistipes                                        -0.28          0.24
    ## Anaerostipes                                      0.98          0.21
    ## Angelakisella                                    -0.03          0.01
    ## Bacteroides                                       0.52          0.21
    ## Bifidobacterium                                  -1.08          0.11
    ## Bilophila                                        -3.04          0.23
    ## Blautia                                           0.32          0.14
    ## Butyricicoccus                                   -1.35          0.09
    ## Butyricimonas                                    -0.78          0.09
    ## CAG.56                                           -2.03          0.08
    ## Christensenellaceae.R.7.group                    -0.62          0.05
    ## Clostridium.sensu.stricto.1                      -1.76          0.20
    ## Colidextribacter                                  0.75          0.11
    ## Collinsella                                      -1.54          0.04
    ## Coprococcus                                      -1.80          0.16
    ## DTU089                                            0.36          0.02
    ## Dialister                                         0.93          0.33
    ## Dorea                                             0.71          0.43
    ## Erysipelatoclostridium                           -0.53          0.02
    ## Erysipelotrichaceae.UCG.003                      -1.06          0.06
    ## Escherichia.Shigella                             -0.56          0.07
    ## Faecalibacterium                                 -1.67          0.12
    ## Family.XIII.AD3011.group                         -1.05          0.02
    ## Family.XIII.UCG.001                               0.54          0.01
    ## Flavonifractor                                   -2.18          0.04
    ## Fusicatenibacter                                 -0.55          0.18
    ## GCA.900066575                                    -0.71          0.01
    ## Haemophilus                                      -0.78          0.31
    ## Holdemania                                       -0.35          0.04
    ## Incertae.Sedis                                   -0.50          0.03
    ## Intestinibacter                                  -1.07          0.03
    ## Intestinimonas                                   -1.01          0.06
    ## Lachnoclostridium                                -0.23          0.07
    ## Lachnospira                                       1.24          0.14
    ## Lachnospiraceae.AC2044.group                     -1.13          0.00
    ## Lachnospiraceae.FCS020.group                     -0.91          0.08
    ## Lachnospiraceae.ND3007.group                     -1.25          0.05
    ## Lachnospiraceae.NK4A136.group                    -0.17          0.07
    ## Lachnospiraceae.UCG.001                          -2.10          0.07
    ## Lachnospiraceae.UCG.004                          -0.69          0.09
    ## Lachnospiraceae.UCG.009                           0.37          0.06
    ## Lachnospiraceae.UCG.010                          -1.07          0.21
    ## Lachnotalea                                      -1.78          0.12
    ## Marvinbryantia                                    0.22          0.05
    ## Methanobrevibacter                                1.00          0.00
    ## Monoglobus                                       -0.96          0.04
    ## NK4A214.group                                    -1.86          0.06
    ## Odoribacter                                       0.01          0.10
    ## Oscillibacter                                    -0.50          0.03
    ## Oscillospira                                      2.62          0.50
    ## Oxalobacter                                       0.00          0.00
    ## Parabacteroides                                  -0.22          0.13
    ## Parasutterella                                   -2.29          0.26
    ## Phascolarctobacterium                             2.60          0.54
    ## Prevotella_9                                      2.70          0.07
    ## Romboutsia                                       -1.68          0.10
    ## Roseburia                                         0.06          0.24
    ## Ruminococcus                                      0.80          0.13
    ## Shuttleworthia                                    1.39          0.06
    ## Streptococcus                                     0.96          0.04
    ## Subdoligranulum                                   0.03          0.16
    ## Sutterella                                        0.29          0.04
    ## Terrisporobacter                                 -1.70          0.09
    ## Turicibacter                                     -0.12          0.17
    ## UBA1819                                          -0.02          0.06
    ## UCG.002                                           0.69          0.02
    ## UCG.003                                          -0.26          0.04
    ## UCG.005                                          -0.43          0.04
    ## UCG.009                                          -1.01          0.01
    ## X.Eubacterium..eligens.group                     -0.58          0.05
    ## X.Eubacterium..fissicatena.group                 -0.86          0.03
    ## X.Eubacterium..hallii.group                       3.19          0.27
    ## X.Eubacterium..ruminantium.group                 -1.03          0.03
    ## X.Eubacterium..siraeum.group                      1.03          0.02
    ## X.Eubacterium..ventriosum.group                   0.44          0.10
    ## X.Eubacterium..xylanophilum.group                -0.43          0.04
    ## X.Ruminococcus..gauvreauii.group                  0.52          0.24
    ## X.Ruminococcus..gnavus.group                     -1.06          0.08
    ## X.Ruminococcus..torques.group                     0.23          0.07
    ## f_Christensenellaceae_ASV181                     -1.00          0.00
    ## f_Coriobacteriales.Incertae.Sedis_ASV433         -0.83          0.01
    ## f_Eggerthellaceae_ASV319                         -2.10          0.12
    ## f_Erysipelotrichaceae_ASV487                      0.26          0.02
    ## f_Lachnospiraceae_ASV117                         -0.44          0.03
    ## f_Lachnospiraceae_ASV239                         -0.46          0.29
    ## f_Lachnospiraceae_ASV259                         -0.36          0.04
    ## f_Lachnospiraceae_ASV275                         -0.08          0.03
    ## f_Lachnospiraceae_ASV283                         -0.78          0.03
    ## f_Lachnospiraceae_ASV297                         -1.09          0.00
    ## f_Lachnospiraceae_ASV326                         -0.31          0.38
    ## f_Lachnospiraceae_ASV384                         -0.14          0.07
    ## f_Lachnospiraceae_ASV398                         -0.29          0.02
    ## f_Lachnospiraceae_ASV410                         -0.84          0.08
    ## f_Lachnospiraceae_ASV514                         -0.73          0.03
    ## f_Lachnospiraceae_ASV532                         -0.01          0.07
    ## f_Lachnospiraceae_ASV555                         -0.75          0.01
    ## f_Lachnospiraceae_ASV558                          0.85          0.10
    ## f_NA_ASV128                                       0.92          0.00
    ## f_NA_ASV156                                       0.04          0.11
    ## f_NA_ASV184                                      -1.14          0.04
    ## f_NA_ASV516                                       1.42          0.00
    ## f_NA_ASV691                                       1.20          0.10
    ## f_Oscillospiraceae_ASV141                        -0.14          0.01
    ## f_Ruminococcaceae_ASV327                          0.82          0.04
    ## f_Ruminococcaceae_ASV610                         -1.15          0.01
    ## f_UCG.010_ASV245                                 -0.88          0.00
    ## f_UCG.010_ASV667                                  1.91          0.03
    ## f_.Eubacterium..coprostanoligenes.group_ASV203   -1.00          0.00

    # Get variable importance from the model fit
    ImpData <- as.data.frame(importance(model1))
    ImpData$Var.Names <- row.names(ImpData)

    ggplot(ImpData, aes(x=Var.Names, y=`%IncMSE`)) +
      geom_segment( aes(x=Var.Names, xend=Var.Names, y=0, yend=`%IncMSE`), color="skyblue") +
      geom_point(aes(size = IncNodePurity), color="blue", alpha=0.6) +
      theme_light() +
      coord_flip() +
      theme(
        legend.position="bottom",
        panel.grid.major.y = element_blank(),
        panel.border = element_blank(),
        axis.ticks.y = element_blank()
      )

<img src="brocc_randomforest_SFNNIT_files/figure-markdown_strict/create default random forest model on training set-1.png" width="98%" height="98%" />

    # since some of these are negative make sure later we use absolute valyue MSE greater than 1

    # ggsave("import_graph_trainSFNNIT.png",
    #        dpi = 300, height = 13,
    #        width = 8)

    # calculate RMSE of train data
    sqrt(model1$mse[length(model1$mse)]) 
    ## [1] 0.7410906
    # 0.7410906
    # want to be low; SFN-NIT logged ranges up to 5

    # calculate RMSE of test data
    predValues <- predict(model1, ValidSet)
    sqrt(mean((ValidSet$sum_SFN_NIT - predValues)^2))
    ## [1] 0.7835532
    # 0.7835532

    R2 <- 1 - (sum((ValidSet$sum_SFN_NIT-predValues)^2)/sum((ValidSet$sum_SFN_NIT-mean(ValidSet$sum_SFN_NIT))^2))
    # -0.612301
    # not a great fit for data? but also found an article that describes that it *could* 
    # be negative because we used different sets of data to train and test so the sums of 
    # squares are not exactly comparable; this is probably not the best metric for fit then
    # https://towardsdatascience.com/explaining-negative-r-squared-17894ca26321/ 

    # find number of trees that produced lowest test MSE
    which.min(model1$mse) # 25
    ## [1] 25

    # RMSE of best model
    sqrt(model1$mse[which.min(model1$mse)])
    ## [1] 0.7045018
    # 0.7045018; pretty good

    #plot the test MSE by number of trees
    plot(model1)

<img src="brocc_randomforest_SFNNIT_files/figure-markdown_strict/create default random forest model on training set-2.png" width="98%" height="98%" />

    model2 <- randomForest(formula = sum_SFN_NIT ~ . , data = feat_tabl, importance = TRUE)
    model2 # not much different than train model
    ## 
    ## Call:
    ##  randomForest(formula = sum_SFN_NIT ~ ., data = feat_tabl, importance = TRUE) 
    ##                Type of random forest: regression
    ##                      Number of trees: 500
    ## No. of variables tried at each split: 37
    ## 
    ##           Mean of squared residuals: 0.547447
    ##                     % Var explained: -21.35

    # find number of trees that produced lowest test MSE
    which.min(model2$mse)
    ## [1] 188
    # 126 trees best, but this changes every time I run it

    # RMSE of best model
    sqrt(model2$mse[which.min(model2$mse)])
    ## [1] 0.7141419
    # 0.7368217 is best mse

    #plot the test MSE by number of trees
    plot(model2)

<img src="brocc_randomforest_SFNNIT_files/figure-markdown_strict/run full model-1.png" width="98%" height="98%" />


    # Get variable importance from the model fit
    ImpData <- as.data.frame(importance(model2))
    ImpData$Var.Names <- row.names(ImpData)

    ggplot(ImpData, aes(x=Var.Names, y=`%IncMSE`)) +
      geom_segment( aes(x=Var.Names, xend=Var.Names, y=0, yend=`%IncMSE`), color="skyblue") +
      geom_point(aes(size = IncNodePurity), color="blue", alpha=0.6) +
      theme_light() +
      coord_flip() +
      theme(
        legend.position="bottom",
        panel.grid.major.y = element_blank(),
        panel.border = element_blank(),
        axis.ticks.y = element_blank()
      )

<img src="brocc_randomforest_SFNNIT_files/figure-markdown_strict/run full model-2.png" width="98%" height="98%" />


    # ggsave("import_graph.png",
    #        dpi = 300, height = 13,
    #        width = 8)

optimize model: - ntreeTry: The number of trees to build. - mtryStart:
The starting number of predictor variables to consider at each split. -
stepFactor: The factor to increase by until the out-of-bag estimated
error stops improving by a certain amount. - improve: The amount that
the out-of-bag error needs to improve by to keep increasing the step
factor.

    model_tuned <- tuneRF(
                   x=feat_tabl[,-112], #define predictor variables
                   y=feat_tabl$sum_SFN_NIT, #define response variable
                   ntreeTry=500,
                   mtryStart=3, 
                   stepFactor=1.5,
                   improve=0.01,
                   trace=FALSE #don't show real-time progress
                   )
    ## 0.02166608 0.01 
    ## 0.003352413 0.01

<img src="brocc_randomforest_SFNNIT_files/figure-markdown_strict/optimize/tune the rf model-1.png" width="98%" height="98%" />

    # 3 predictors (mtry) at each split will give best OOB

    # run model with best trees
    model3 <- randomForest(formula = sum_SFN_NIT ~ . , data = feat_tabl, 
                           importance = TRUE, ntree = 500, mtry = 4)
    model3 # not much different than train model
    ## 
    ## Call:
    ##  randomForest(formula = sum_SFN_NIT ~ ., data = feat_tabl, importance = TRUE,      ntree = 500, mtry = 4) 
    ##                Type of random forest: regression
    ##                      Number of trees: 500
    ## No. of variables tried at each split: 4
    ## 
    ##           Mean of squared residuals: 0.5110997
    ##                     % Var explained: -13.29

    # RMSE of best model
    sqrt(model3$mse[which.min(model3$mse)])
    ## [1] 0.7051657
    # 0.7108695
    # compare to range of Dialister abundance for e.g.
    range(feat_tabl$Dialister)
    ## [1] 0.000000 6.265301
    # so the error is not too bad

    #plot the test MSE by number of trees
    plot(model3)

<img src="brocc_randomforest_SFNNIT_files/figure-markdown_strict/rerun rf model with new information about best number of predictors-1.png" width="98%" height="98%" />


    # Get variable importance from the model fit
    ImpData <- as.data.frame(importance(model3)) %>%
      rename_at('%IncMSE', ~'percIncMSE') %>%
      dplyr::filter(percIncMSE > 1)
    ImpData$Var.Names <- row.names(ImpData)

    mod3Impgraph <- ggplot(ImpData, aes(x=Var.Names, y=`percIncMSE`)) +
      geom_segment( aes(x=Var.Names, xend=Var.Names, y=0, yend=`percIncMSE`), color="skyblue") +
      geom_point(aes(size = IncNodePurity), color="blue", alpha=0.6) +
      theme_light() +
      coord_flip() +
      theme(legend.position="bottom",
        panel.grid.major.y = element_blank(),
        panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        text = element_text(size = 14)) +
      labs(y = "Percent increase in MSE", x = "", size = "Increase in Node Purity")
    mod3Impgraph

<img src="brocc_randomforest_SFNNIT_files/figure-markdown_strict/rerun rf model with new information about best number of predictors-2.png" width="98%" height="98%" />


    # ggsave("import_graph_optimtestmod_SFNNIT.png",
    #        dpi = 300, height = 9,
    #        width = 8)

    # psuedo r squared
    predValues <- predict(model3, ValidSet)
    R2 <- 1 - (sum((ValidSet$sum_SFN_NIT-predValues)^2)/sum((ValidSet$sum_SFN_NIT-mean(ValidSet$sum_SFN_NIT))^2))
    R2
    ## [1] 0.7546221
    # pretty ok
    # 0.7452101

    # list the taxa with an perincmse higher than 2 and node purity inc higher than 0.1
    top_tax <- ImpData %>%
      dplyr::filter(percIncMSE > 2 & IncNodePurity > 0.1)
    # note this uses logged tax abundances and SFN
    toptax_graph <- feat_tabl %>%
      dplyr::select(c(rownames(top_tax), "sum_SFN_NIT"))

    # stats
    summary(lm(sum_SFN_NIT ~ Dialister, data = toptax_graph)) # nope
    ## 
    ## Call:
    ## lm(formula = sum_SFN_NIT ~ Dialister, data = toptax_graph)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -1.78403 -0.40465  0.09618  0.47904  1.12636 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  3.32906    0.14352  23.196   <2e-16 ***
    ## Dialister    0.06017    0.04900   1.228    0.227    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.6761 on 36 degrees of freedom
    ## Multiple R-squared:  0.0402, Adjusted R-squared:  0.01354 
    ## F-statistic: 1.508 on 1 and 36 DF,  p-value: 0.2274
    summary(lm(sum_SFN_NIT ~ Lachnospiraceae.UCG.009, data = toptax_graph)) # almost
    ## 
    ## Call:
    ## lm(formula = sum_SFN_NIT ~ Lachnospiraceae.UCG.009, data = toptax_graph)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -1.40023 -0.43979  0.02987  0.50879  1.17229 
    ## 
    ## Coefficients:
    ##                         Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)               3.2051     0.1669  19.199   <2e-16 ***
    ## Lachnospiraceae.UCG.009   0.2399     0.1294   1.854   0.0719 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.6593 on 36 degrees of freedom
    ## Multiple R-squared:  0.08718,    Adjusted R-squared:  0.06182 
    ## F-statistic: 3.438 on 1 and 36 DF,  p-value: 0.07192
    summary(lm(sum_SFN_NIT ~ Oscillospira, data = toptax_graph)) # yes
    ## 
    ## Call:
    ## lm(formula = sum_SFN_NIT ~ Oscillospira, data = toptax_graph)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -1.01104 -0.53079 -0.02775  0.33548  1.39339 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   3.74911    0.14533  25.797  < 2e-16 ***
    ## Oscillospira -0.28992    0.09922  -2.922  0.00597 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.6204 on 36 degrees of freedom
    ## Multiple R-squared:  0.1917, Adjusted R-squared:  0.1693 
    ## F-statistic: 8.539 on 1 and 36 DF,  p-value: 0.005973
    summary(lm(sum_SFN_NIT ~ f_NA_ASV156, data = toptax_graph)) # yes
    ## 
    ## Call:
    ## lm(formula = sum_SFN_NIT ~ f_NA_ASV156, data = toptax_graph)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -1.54177 -0.47578 -0.08091  0.38393  1.43861 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  3.34661    0.11599   28.85   <2e-16 ***
    ## f_NA_ASV156  0.17105    0.08386    2.04   0.0488 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.6533 on 36 degrees of freedom
    ## Multiple R-squared:  0.1036, Adjusted R-squared:  0.07869 
    ## F-statistic:  4.16 on 1 and 36 DF,  p-value: 0.04877

    # make a little table of these to add as annotations to graph later
    lachno <- summary(lm(sum_SFN_NIT ~ Lachnospiraceae.UCG.009, data = toptax_graph)) # almost
    osci <- summary(lm(sum_SFN_NIT ~ Oscillospira, data = toptax_graph)) # yes
    NAasv <- summary(lm(sum_SFN_NIT ~ f_NA_ASV156, data = toptax_graph)) # yes
    stat_tab <- data.frame(Taxon = c(rownames(lachno$coefficients)[2],
                                     rownames(osci$coefficients)[2],
                                     rownames(NAasv$coefficients)[2]),
                           Est = round(c(lachno$coefficients[2,1],
                                     osci$coefficients[2,1],
                                     NAasv$coefficients[2,1]), digits = 1),
                           Error = round(c(lachno$coefficients[2,2],
                                     osci$coefficients[2,2],
                                     NAasv$coefficients[2,2]), digits = 2),
                           pval = round(c(lachno$coefficients[2,4],
                                     osci$coefficients[2,4],
                                     NAasv$coefficients[2,4]), digits = 2))
    # turn into text object for.graph
    text_annot <- data.frame(label = c(paste0("Estimate = ", stat_tab$Est[1],
                                                    " +/- ", stat_tab$Error[1], "\n",
                                                    "p = ", stat_tab$pval[1]),
                                       paste0("Estimate = ", stat_tab$Est[2],
                                                    " +/- ", stat_tab$Error[2], "\n",
                                                    "p = ", stat_tab$pval[2]),
                                       paste0("Estimate = ", stat_tab$Est[3],
                                                    " +/- ", stat_tab$Error[3], "\n",
                                                    "p = ", stat_tab$pval[3])),
                             Taxon = stat_tab$Taxon)

    # change format for graphing
    toptax_graph_long <- toptax_graph %>%
      tidyr::gather(key = "Taxon", value = "abund", -sum_SFN_NIT)

    # graph
    taxnit_plot <- ggplot(data = toptax_graph_long,
                           aes(x = abund, y = sum_SFN_NIT)) +
      geom_point() +
      geom_smooth(method = 'lm', color = "black") +
      labs(x = "Log taxon abundance", y = "log cumulative SFN-nitrile in urine") +
      theme(text = element_text(size = 14)) +
      facet_wrap("Taxon")

    taxnit_plot

<img src="brocc_randomforest_SFNNIT_files/figure-markdown_strict/graph and linear regression of dialister and SFN-NIT-1.png" width="98%" height="98%" />


    # ggsave("TOP5NIT_lmplot.jpg", taxnit_plot,
    #        dpi = 300, width = 7, height = 5)

    # different version of graph
    # change format for graphing
    toptax_graph_long2 <- toptax_graph %>%
      dplyr::select(c("Oscillospira", "f_NA_ASV156", "Lachnospiraceae.UCG.009", "sum_SFN_NIT")) %>%
      tidyr::gather(key = "Taxon", value = "abund", -sum_SFN_NIT)
    # graph
    taxnit_plot2 <- ggplot(data = toptax_graph_long2,
                           aes(x = abund, y = sum_SFN_NIT)) +
      geom_point() +
      geom_smooth(method = 'lm', color = "black") +
      labs(x = "Log taxon abundance", y = "log cumulative SFN-nitrile in urine") +
      theme(text = element_text(size = 14)) +
      facet_wrap("Taxon") +
      geom_text(data = text_annot, mapping = aes(x = 2.5, y = 2, label = label))

    taxnit_plot2

<img src="brocc_randomforest_SFNNIT_files/figure-markdown_strict/graph and linear regression of dialister and SFN-NIT-2.png" width="98%" height="98%" />


    # ggsave("TOP3NIT_lmplot.jpg", taxnit_plot2,
    #        dpi = 300, width = 9, height = 4)

    library(ggpubr)
    implm <- ggarrange(mod3Impgraph, taxnit_plot2,
                       nrow = 2, ncol = 1,
                       heights = c(8, 4),
                       labels = c("A", "B"))
    implm

<img src="brocc_randomforest_SFNNIT_files/figure-markdown_strict/make a paneled figure of importance graph and regressions-1.png" width="98%" height="98%" />


    # ggsave("imp_and_lm_rfplot.jpg", implm,
    #        dpi = 300, width = 9, height = 12)

## random forest model to test whether the presence or absence of certain taxa is associated with cumulative SFN\_NIT

    metadata <- read.delim("../input_data/metadata_wgroups_brocc.txt",
                           header = T, sep = "\t") %>%
      dplyr::select(c("subject_id", "group"))
    MB_asvtab_pa <- column_to_rownames(MB_asvtab, "subject_id")
    MB_asvtab_pa <- as.data.frame(ifelse(MB_asvtab_pa == 0, "A", "P")) %>%
      rownames_to_column("subject_id") %>%
      inner_join(SFN_NIT_tab[c("subject_id", "sum_SFN_NIT")], by = "subject_id") %>%
      column_to_rownames("subject_id")

    ## separate into each group
    MBasvtab_pa_A <- dplyr::filter(MB_asvtab_pa %>% rownames_to_column("subject_id"), 
                                   subject_id %in% metadata$subject_id[metadata$group == "A"])
    MBasvtab_pa_B <- dplyr::filter(MB_asvtab_pa %>% rownames_to_column("subject_id"), 
                                   subject_id %in% metadata$subject_id[metadata$group == "B"])
    MBasvtab_pa_C <- dplyr::filter(MB_asvtab_pa %>% rownames_to_column("subject_id"), 
                                   subject_id %in% metadata$subject_id[metadata$group == "C"])

    # Training Set : Validation Set = 70 : 30 (random)
    list <- list(MBasvtab_pa_A, MBasvtab_pa_B, MBasvtab_pa_C)

    for (i in list) {
      train <- sample(nrow(i), 0.7*nrow(i), replace = FALSE)
      TrainSet <- i[train,]
      ValidSet <- i[-train,]
      
    }

    # for whole data
    train <- sample(nrow(MB_asvtab_pa), 0.7*nrow(MB_asvtab_pa), replace = FALSE)
    TrainSet <- MB_asvtab_pa[train,]
    ValidSet <- MB_asvtab_pa[-train,]
    summary(TrainSet)
    ##  Adlercreutzia      Agathobacter       Akkermansia         Alistipes        
    ##  Length:26          Length:26          Length:26          Length:26         
    ##  Class :character   Class :character   Class :character   Class :character  
    ##  Mode  :character   Mode  :character   Mode  :character   Mode  :character  
    ##                                                                             
    ##                                                                             
    ##                                                                             
    ##  Anaerostipes       Angelakisella      Bacteroides        Bifidobacterium   
    ##  Length:26          Length:26          Length:26          Length:26         
    ##  Class :character   Class :character   Class :character   Class :character  
    ##  Mode  :character   Mode  :character   Mode  :character   Mode  :character  
    ##                                                                             
    ##                                                                             
    ##                                                                             
    ##   Bilophila           Blautia          Butyricicoccus     Butyricimonas     
    ##  Length:26          Length:26          Length:26          Length:26         
    ##  Class :character   Class :character   Class :character   Class :character  
    ##  Mode  :character   Mode  :character   Mode  :character   Mode  :character  
    ##                                                                             
    ##                                                                             
    ##                                                                             
    ##     CAG.56          Christensenellaceae.R.7.group Clostridium.sensu.stricto.1
    ##  Length:26          Length:26                     Length:26                  
    ##  Class :character   Class :character              Class :character           
    ##  Mode  :character   Mode  :character              Mode  :character           
    ##                                                                              
    ##                                                                              
    ##                                                                              
    ##  Colidextribacter   Collinsella        Coprococcus           DTU089         
    ##  Length:26          Length:26          Length:26          Length:26         
    ##  Class :character   Class :character   Class :character   Class :character  
    ##  Mode  :character   Mode  :character   Mode  :character   Mode  :character  
    ##                                                                             
    ##                                                                             
    ##                                                                             
    ##   Dialister            Dorea           Erysipelatoclostridium
    ##  Length:26          Length:26          Length:26             
    ##  Class :character   Class :character   Class :character      
    ##  Mode  :character   Mode  :character   Mode  :character      
    ##                                                              
    ##                                                              
    ##                                                              
    ##  Erysipelotrichaceae.UCG.003 Escherichia.Shigella Faecalibacterium  
    ##  Length:26                   Length:26            Length:26         
    ##  Class :character            Class :character     Class :character  
    ##  Mode  :character            Mode  :character     Mode  :character  
    ##                                                                     
    ##                                                                     
    ##                                                                     
    ##  Family.XIII.AD3011.group Family.XIII.UCG.001 Flavonifractor    
    ##  Length:26                Length:26           Length:26         
    ##  Class :character         Class :character    Class :character  
    ##  Mode  :character         Mode  :character    Mode  :character  
    ##                                                                 
    ##                                                                 
    ##                                                                 
    ##  Fusicatenibacter   GCA.900066575      Haemophilus         Holdemania       
    ##  Length:26          Length:26          Length:26          Length:26         
    ##  Class :character   Class :character   Class :character   Class :character  
    ##  Mode  :character   Mode  :character   Mode  :character   Mode  :character  
    ##                                                                             
    ##                                                                             
    ##                                                                             
    ##  Incertae.Sedis     Intestinibacter    Intestinimonas     Lachnoclostridium 
    ##  Length:26          Length:26          Length:26          Length:26         
    ##  Class :character   Class :character   Class :character   Class :character  
    ##  Mode  :character   Mode  :character   Mode  :character   Mode  :character  
    ##                                                                             
    ##                                                                             
    ##                                                                             
    ##  Lachnospira        Lachnospiraceae.AC2044.group Lachnospiraceae.FCS020.group
    ##  Length:26          Length:26                    Length:26                   
    ##  Class :character   Class :character             Class :character            
    ##  Mode  :character   Mode  :character             Mode  :character            
    ##                                                                              
    ##                                                                              
    ##                                                                              
    ##  Lachnospiraceae.ND3007.group Lachnospiraceae.NK4A136.group
    ##  Length:26                    Length:26                    
    ##  Class :character             Class :character             
    ##  Mode  :character             Mode  :character             
    ##                                                            
    ##                                                            
    ##                                                            
    ##  Lachnospiraceae.UCG.001 Lachnospiraceae.UCG.004 Lachnospiraceae.UCG.009
    ##  Length:26               Length:26               Length:26              
    ##  Class :character        Class :character        Class :character       
    ##  Mode  :character        Mode  :character        Mode  :character       
    ##                                                                         
    ##                                                                         
    ##                                                                         
    ##  Lachnospiraceae.UCG.010 Lachnotalea        Marvinbryantia    
    ##  Length:26               Length:26          Length:26         
    ##  Class :character        Class :character   Class :character  
    ##  Mode  :character        Mode  :character   Mode  :character  
    ##                                                               
    ##                                                               
    ##                                                               
    ##  Methanobrevibacter  Monoglobus        NK4A214.group      Odoribacter       
    ##  Length:26          Length:26          Length:26          Length:26         
    ##  Class :character   Class :character   Class :character   Class :character  
    ##  Mode  :character   Mode  :character   Mode  :character   Mode  :character  
    ##                                                                             
    ##                                                                             
    ##                                                                             
    ##  Oscillibacter      Oscillospira       Oxalobacter        Parabacteroides   
    ##  Length:26          Length:26          Length:26          Length:26         
    ##  Class :character   Class :character   Class :character   Class :character  
    ##  Mode  :character   Mode  :character   Mode  :character   Mode  :character  
    ##                                                                             
    ##                                                                             
    ##                                                                             
    ##  Parasutterella     Phascolarctobacterium Prevotella_9        Romboutsia       
    ##  Length:26          Length:26             Length:26          Length:26         
    ##  Class :character   Class :character      Class :character   Class :character  
    ##  Mode  :character   Mode  :character      Mode  :character   Mode  :character  
    ##                                                                                
    ##                                                                                
    ##                                                                                
    ##   Roseburia         Ruminococcus       Shuttleworthia     Streptococcus     
    ##  Length:26          Length:26          Length:26          Length:26         
    ##  Class :character   Class :character   Class :character   Class :character  
    ##  Mode  :character   Mode  :character   Mode  :character   Mode  :character  
    ##                                                                             
    ##                                                                             
    ##                                                                             
    ##  Subdoligranulum     Sutterella        Terrisporobacter   Turicibacter      
    ##  Length:26          Length:26          Length:26          Length:26         
    ##  Class :character   Class :character   Class :character   Class :character  
    ##  Mode  :character   Mode  :character   Mode  :character   Mode  :character  
    ##                                                                             
    ##                                                                             
    ##                                                                             
    ##    UBA1819            UCG.002            UCG.003            UCG.005         
    ##  Length:26          Length:26          Length:26          Length:26         
    ##  Class :character   Class :character   Class :character   Class :character  
    ##  Mode  :character   Mode  :character   Mode  :character   Mode  :character  
    ##                                                                             
    ##                                                                             
    ##                                                                             
    ##    UCG.009          X.Eubacterium..eligens.group
    ##  Length:26          Length:26                   
    ##  Class :character   Class :character            
    ##  Mode  :character   Mode  :character            
    ##                                                 
    ##                                                 
    ##                                                 
    ##  X.Eubacterium..fissicatena.group X.Eubacterium..hallii.group
    ##  Length:26                        Length:26                  
    ##  Class :character                 Class :character           
    ##  Mode  :character                 Mode  :character           
    ##                                                              
    ##                                                              
    ##                                                              
    ##  X.Eubacterium..ruminantium.group X.Eubacterium..siraeum.group
    ##  Length:26                        Length:26                   
    ##  Class :character                 Class :character            
    ##  Mode  :character                 Mode  :character            
    ##                                                               
    ##                                                               
    ##                                                               
    ##  X.Eubacterium..ventriosum.group X.Eubacterium..xylanophilum.group
    ##  Length:26                       Length:26                        
    ##  Class :character                Class :character                 
    ##  Mode  :character                Mode  :character                 
    ##                                                                   
    ##                                                                   
    ##                                                                   
    ##  X.Ruminococcus..gauvreauii.group X.Ruminococcus..gnavus.group
    ##  Length:26                        Length:26                   
    ##  Class :character                 Class :character            
    ##  Mode  :character                 Mode  :character            
    ##                                                               
    ##                                                               
    ##                                                               
    ##  X.Ruminococcus..torques.group f_Christensenellaceae_ASV181
    ##  Length:26                     Length:26                   
    ##  Class :character              Class :character            
    ##  Mode  :character              Mode  :character            
    ##                                                            
    ##                                                            
    ##                                                            
    ##  f_Coriobacteriales.Incertae.Sedis_ASV433 f_Eggerthellaceae_ASV319
    ##  Length:26                                Length:26               
    ##  Class :character                         Class :character        
    ##  Mode  :character                         Mode  :character        
    ##                                                                   
    ##                                                                   
    ##                                                                   
    ##  f_Erysipelotrichaceae_ASV487 f_Lachnospiraceae_ASV117 f_Lachnospiraceae_ASV239
    ##  Length:26                    Length:26                Length:26               
    ##  Class :character             Class :character         Class :character        
    ##  Mode  :character             Mode  :character         Mode  :character        
    ##                                                                                
    ##                                                                                
    ##                                                                                
    ##  f_Lachnospiraceae_ASV259 f_Lachnospiraceae_ASV275 f_Lachnospiraceae_ASV283
    ##  Length:26                Length:26                Length:26               
    ##  Class :character         Class :character         Class :character        
    ##  Mode  :character         Mode  :character         Mode  :character        
    ##                                                                            
    ##                                                                            
    ##                                                                            
    ##  f_Lachnospiraceae_ASV297 f_Lachnospiraceae_ASV326 f_Lachnospiraceae_ASV384
    ##  Length:26                Length:26                Length:26               
    ##  Class :character         Class :character         Class :character        
    ##  Mode  :character         Mode  :character         Mode  :character        
    ##                                                                            
    ##                                                                            
    ##                                                                            
    ##  f_Lachnospiraceae_ASV398 f_Lachnospiraceae_ASV410 f_Lachnospiraceae_ASV514
    ##  Length:26                Length:26                Length:26               
    ##  Class :character         Class :character         Class :character        
    ##  Mode  :character         Mode  :character         Mode  :character        
    ##                                                                            
    ##                                                                            
    ##                                                                            
    ##  f_Lachnospiraceae_ASV532 f_Lachnospiraceae_ASV555 f_Lachnospiraceae_ASV558
    ##  Length:26                Length:26                Length:26               
    ##  Class :character         Class :character         Class :character        
    ##  Mode  :character         Mode  :character         Mode  :character        
    ##                                                                            
    ##                                                                            
    ##                                                                            
    ##  f_NA_ASV128        f_NA_ASV156        f_NA_ASV184        f_NA_ASV516       
    ##  Length:26          Length:26          Length:26          Length:26         
    ##  Class :character   Class :character   Class :character   Class :character  
    ##  Mode  :character   Mode  :character   Mode  :character   Mode  :character  
    ##                                                                             
    ##                                                                             
    ##                                                                             
    ##  f_NA_ASV691        f_Oscillospiraceae_ASV141 f_Ruminococcaceae_ASV327
    ##  Length:26          Length:26                 Length:26               
    ##  Class :character   Class :character          Class :character        
    ##  Mode  :character   Mode  :character          Mode  :character        
    ##                                                                       
    ##                                                                       
    ##                                                                       
    ##  f_Ruminococcaceae_ASV610 f_UCG.010_ASV245   f_UCG.010_ASV667  
    ##  Length:26                Length:26          Length:26         
    ##  Class :character         Class :character   Class :character  
    ##  Mode  :character         Mode  :character   Mode  :character  
    ##                                                                
    ##                                                                
    ##                                                                
    ##  f_.Eubacterium..coprostanoligenes.group_ASV203  sum_SFN_NIT     
    ##  Length:26                                      Min.   :  9.089  
    ##  Class :character                               1st Qu.: 16.575  
    ##  Mode  :character                               Median : 26.694  
    ##                                                 Mean   : 32.950  
    ##                                                 3rd Qu.: 41.113  
    ##                                                 Max.   :118.728
    summary(ValidSet)
    ##  Adlercreutzia      Agathobacter       Akkermansia         Alistipes        
    ##  Length:12          Length:12          Length:12          Length:12         
    ##  Class :character   Class :character   Class :character   Class :character  
    ##  Mode  :character   Mode  :character   Mode  :character   Mode  :character  
    ##                                                                             
    ##                                                                             
    ##                                                                             
    ##  Anaerostipes       Angelakisella      Bacteroides        Bifidobacterium   
    ##  Length:12          Length:12          Length:12          Length:12         
    ##  Class :character   Class :character   Class :character   Class :character  
    ##  Mode  :character   Mode  :character   Mode  :character   Mode  :character  
    ##                                                                             
    ##                                                                             
    ##                                                                             
    ##   Bilophila           Blautia          Butyricicoccus     Butyricimonas     
    ##  Length:12          Length:12          Length:12          Length:12         
    ##  Class :character   Class :character   Class :character   Class :character  
    ##  Mode  :character   Mode  :character   Mode  :character   Mode  :character  
    ##                                                                             
    ##                                                                             
    ##                                                                             
    ##     CAG.56          Christensenellaceae.R.7.group Clostridium.sensu.stricto.1
    ##  Length:12          Length:12                     Length:12                  
    ##  Class :character   Class :character              Class :character           
    ##  Mode  :character   Mode  :character              Mode  :character           
    ##                                                                              
    ##                                                                              
    ##                                                                              
    ##  Colidextribacter   Collinsella        Coprococcus           DTU089         
    ##  Length:12          Length:12          Length:12          Length:12         
    ##  Class :character   Class :character   Class :character   Class :character  
    ##  Mode  :character   Mode  :character   Mode  :character   Mode  :character  
    ##                                                                             
    ##                                                                             
    ##                                                                             
    ##   Dialister            Dorea           Erysipelatoclostridium
    ##  Length:12          Length:12          Length:12             
    ##  Class :character   Class :character   Class :character      
    ##  Mode  :character   Mode  :character   Mode  :character      
    ##                                                              
    ##                                                              
    ##                                                              
    ##  Erysipelotrichaceae.UCG.003 Escherichia.Shigella Faecalibacterium  
    ##  Length:12                   Length:12            Length:12         
    ##  Class :character            Class :character     Class :character  
    ##  Mode  :character            Mode  :character     Mode  :character  
    ##                                                                     
    ##                                                                     
    ##                                                                     
    ##  Family.XIII.AD3011.group Family.XIII.UCG.001 Flavonifractor    
    ##  Length:12                Length:12           Length:12         
    ##  Class :character         Class :character    Class :character  
    ##  Mode  :character         Mode  :character    Mode  :character  
    ##                                                                 
    ##                                                                 
    ##                                                                 
    ##  Fusicatenibacter   GCA.900066575      Haemophilus         Holdemania       
    ##  Length:12          Length:12          Length:12          Length:12         
    ##  Class :character   Class :character   Class :character   Class :character  
    ##  Mode  :character   Mode  :character   Mode  :character   Mode  :character  
    ##                                                                             
    ##                                                                             
    ##                                                                             
    ##  Incertae.Sedis     Intestinibacter    Intestinimonas     Lachnoclostridium 
    ##  Length:12          Length:12          Length:12          Length:12         
    ##  Class :character   Class :character   Class :character   Class :character  
    ##  Mode  :character   Mode  :character   Mode  :character   Mode  :character  
    ##                                                                             
    ##                                                                             
    ##                                                                             
    ##  Lachnospira        Lachnospiraceae.AC2044.group Lachnospiraceae.FCS020.group
    ##  Length:12          Length:12                    Length:12                   
    ##  Class :character   Class :character             Class :character            
    ##  Mode  :character   Mode  :character             Mode  :character            
    ##                                                                              
    ##                                                                              
    ##                                                                              
    ##  Lachnospiraceae.ND3007.group Lachnospiraceae.NK4A136.group
    ##  Length:12                    Length:12                    
    ##  Class :character             Class :character             
    ##  Mode  :character             Mode  :character             
    ##                                                            
    ##                                                            
    ##                                                            
    ##  Lachnospiraceae.UCG.001 Lachnospiraceae.UCG.004 Lachnospiraceae.UCG.009
    ##  Length:12               Length:12               Length:12              
    ##  Class :character        Class :character        Class :character       
    ##  Mode  :character        Mode  :character        Mode  :character       
    ##                                                                         
    ##                                                                         
    ##                                                                         
    ##  Lachnospiraceae.UCG.010 Lachnotalea        Marvinbryantia    
    ##  Length:12               Length:12          Length:12         
    ##  Class :character        Class :character   Class :character  
    ##  Mode  :character        Mode  :character   Mode  :character  
    ##                                                               
    ##                                                               
    ##                                                               
    ##  Methanobrevibacter  Monoglobus        NK4A214.group      Odoribacter       
    ##  Length:12          Length:12          Length:12          Length:12         
    ##  Class :character   Class :character   Class :character   Class :character  
    ##  Mode  :character   Mode  :character   Mode  :character   Mode  :character  
    ##                                                                             
    ##                                                                             
    ##                                                                             
    ##  Oscillibacter      Oscillospira       Oxalobacter        Parabacteroides   
    ##  Length:12          Length:12          Length:12          Length:12         
    ##  Class :character   Class :character   Class :character   Class :character  
    ##  Mode  :character   Mode  :character   Mode  :character   Mode  :character  
    ##                                                                             
    ##                                                                             
    ##                                                                             
    ##  Parasutterella     Phascolarctobacterium Prevotella_9        Romboutsia       
    ##  Length:12          Length:12             Length:12          Length:12         
    ##  Class :character   Class :character      Class :character   Class :character  
    ##  Mode  :character   Mode  :character      Mode  :character   Mode  :character  
    ##                                                                                
    ##                                                                                
    ##                                                                                
    ##   Roseburia         Ruminococcus       Shuttleworthia     Streptococcus     
    ##  Length:12          Length:12          Length:12          Length:12         
    ##  Class :character   Class :character   Class :character   Class :character  
    ##  Mode  :character   Mode  :character   Mode  :character   Mode  :character  
    ##                                                                             
    ##                                                                             
    ##                                                                             
    ##  Subdoligranulum     Sutterella        Terrisporobacter   Turicibacter      
    ##  Length:12          Length:12          Length:12          Length:12         
    ##  Class :character   Class :character   Class :character   Class :character  
    ##  Mode  :character   Mode  :character   Mode  :character   Mode  :character  
    ##                                                                             
    ##                                                                             
    ##                                                                             
    ##    UBA1819            UCG.002            UCG.003            UCG.005         
    ##  Length:12          Length:12          Length:12          Length:12         
    ##  Class :character   Class :character   Class :character   Class :character  
    ##  Mode  :character   Mode  :character   Mode  :character   Mode  :character  
    ##                                                                             
    ##                                                                             
    ##                                                                             
    ##    UCG.009          X.Eubacterium..eligens.group
    ##  Length:12          Length:12                   
    ##  Class :character   Class :character            
    ##  Mode  :character   Mode  :character            
    ##                                                 
    ##                                                 
    ##                                                 
    ##  X.Eubacterium..fissicatena.group X.Eubacterium..hallii.group
    ##  Length:12                        Length:12                  
    ##  Class :character                 Class :character           
    ##  Mode  :character                 Mode  :character           
    ##                                                              
    ##                                                              
    ##                                                              
    ##  X.Eubacterium..ruminantium.group X.Eubacterium..siraeum.group
    ##  Length:12                        Length:12                   
    ##  Class :character                 Class :character            
    ##  Mode  :character                 Mode  :character            
    ##                                                               
    ##                                                               
    ##                                                               
    ##  X.Eubacterium..ventriosum.group X.Eubacterium..xylanophilum.group
    ##  Length:12                       Length:12                        
    ##  Class :character                Class :character                 
    ##  Mode  :character                Mode  :character                 
    ##                                                                   
    ##                                                                   
    ##                                                                   
    ##  X.Ruminococcus..gauvreauii.group X.Ruminococcus..gnavus.group
    ##  Length:12                        Length:12                   
    ##  Class :character                 Class :character            
    ##  Mode  :character                 Mode  :character            
    ##                                                               
    ##                                                               
    ##                                                               
    ##  X.Ruminococcus..torques.group f_Christensenellaceae_ASV181
    ##  Length:12                     Length:12                   
    ##  Class :character              Class :character            
    ##  Mode  :character              Mode  :character            
    ##                                                            
    ##                                                            
    ##                                                            
    ##  f_Coriobacteriales.Incertae.Sedis_ASV433 f_Eggerthellaceae_ASV319
    ##  Length:12                                Length:12               
    ##  Class :character                         Class :character        
    ##  Mode  :character                         Mode  :character        
    ##                                                                   
    ##                                                                   
    ##                                                                   
    ##  f_Erysipelotrichaceae_ASV487 f_Lachnospiraceae_ASV117 f_Lachnospiraceae_ASV239
    ##  Length:12                    Length:12                Length:12               
    ##  Class :character             Class :character         Class :character        
    ##  Mode  :character             Mode  :character         Mode  :character        
    ##                                                                                
    ##                                                                                
    ##                                                                                
    ##  f_Lachnospiraceae_ASV259 f_Lachnospiraceae_ASV275 f_Lachnospiraceae_ASV283
    ##  Length:12                Length:12                Length:12               
    ##  Class :character         Class :character         Class :character        
    ##  Mode  :character         Mode  :character         Mode  :character        
    ##                                                                            
    ##                                                                            
    ##                                                                            
    ##  f_Lachnospiraceae_ASV297 f_Lachnospiraceae_ASV326 f_Lachnospiraceae_ASV384
    ##  Length:12                Length:12                Length:12               
    ##  Class :character         Class :character         Class :character        
    ##  Mode  :character         Mode  :character         Mode  :character        
    ##                                                                            
    ##                                                                            
    ##                                                                            
    ##  f_Lachnospiraceae_ASV398 f_Lachnospiraceae_ASV410 f_Lachnospiraceae_ASV514
    ##  Length:12                Length:12                Length:12               
    ##  Class :character         Class :character         Class :character        
    ##  Mode  :character         Mode  :character         Mode  :character        
    ##                                                                            
    ##                                                                            
    ##                                                                            
    ##  f_Lachnospiraceae_ASV532 f_Lachnospiraceae_ASV555 f_Lachnospiraceae_ASV558
    ##  Length:12                Length:12                Length:12               
    ##  Class :character         Class :character         Class :character        
    ##  Mode  :character         Mode  :character         Mode  :character        
    ##                                                                            
    ##                                                                            
    ##                                                                            
    ##  f_NA_ASV128        f_NA_ASV156        f_NA_ASV184        f_NA_ASV516       
    ##  Length:12          Length:12          Length:12          Length:12         
    ##  Class :character   Class :character   Class :character   Class :character  
    ##  Mode  :character   Mode  :character   Mode  :character   Mode  :character  
    ##                                                                             
    ##                                                                             
    ##                                                                             
    ##  f_NA_ASV691        f_Oscillospiraceae_ASV141 f_Ruminococcaceae_ASV327
    ##  Length:12          Length:12                 Length:12               
    ##  Class :character   Class :character          Class :character        
    ##  Mode  :character   Mode  :character          Mode  :character        
    ##                                                                       
    ##                                                                       
    ##                                                                       
    ##  f_Ruminococcaceae_ASV610 f_UCG.010_ASV245   f_UCG.010_ASV667  
    ##  Length:12                Length:12          Length:12         
    ##  Class :character         Class :character   Class :character  
    ##  Mode  :character         Mode  :character   Mode  :character  
    ##                                                                
    ##                                                                
    ##                                                                
    ##  f_.Eubacterium..coprostanoligenes.group_ASV203  sum_SFN_NIT     
    ##  Length:12                                      Min.   :  5.079  
    ##  Class :character                               1st Qu.: 30.363  
    ##  Mode  :character                               Median : 45.622  
    ##                                                 Mean   : 48.387  
    ##                                                 3rd Qu.: 60.816  
    ##                                                 Max.   :108.909

    # Create a Random Forest model with default parameters
    model1 <- randomForest(sum_SFN_NIT ~ ., data = TrainSet, importance = TRUE, type = "class")
    model1
    ## 
    ## Call:
    ##  randomForest(formula = sum_SFN_NIT ~ ., data = TrainSet, importance = TRUE,      type = "class") 
    ##                Type of random forest: regression
    ##                      Number of trees: 500
    ## No. of variables tried at each split: 37
    ## 
    ##           Mean of squared residuals: 770.7486
    ##                     % Var explained: -38.47

    R2 <- 1 - (sum((ValidSet$sum_SFN_NIT-predValues)^2)/sum((ValidSet$sum_SFN_NIT-mean(ValidSet$sum_SFN_NIT))^2))
    # negative, not a great fit for data

    # find number of trees that produced lowest test MSE
    which.min(model1$mse)
    ## [1] 496
    # 500

    # RMSE of best model
    sqrt(model1$mse[which.min(model1$mse)])
    ## [1] 27.74876
    # 27.4

    model2 <- randomForest(formula = sum_SFN_NIT ~ . , data = MB_asvtab_pa, importance = TRUE)
    model2 # not much different than train model
    ## 
    ## Call:
    ##  randomForest(formula = sum_SFN_NIT ~ ., data = MB_asvtab_pa,      importance = TRUE) 
    ##                Type of random forest: regression
    ##                      Number of trees: 500
    ## No. of variables tried at each split: 37
    ## 
    ##           Mean of squared residuals: 716.0811
    ##                     % Var explained: -3.54

    # find number of trees that produced lowest test MSE
    which.min(model2$mse)
    ## [1] 1
    # 106 trees best

    # RMSE of best model
    sqrt(model1$mse[which.min(model2$mse)])
    ## [1] 51.35442
    # 28 is best mse

    #plot the test MSE by number of trees
    plot(model2)

<img src="brocc_randomforest_SFNNIT_files/figure-markdown_strict/run full model 2-1.png" width="98%" height="98%" />

optimize model: - ntreeTry: The number of trees to build. - mtryStart:
The starting number of predictor variables to consider at each split. -
stepFactor: The factor to increase by until the out-of-bag estimated
error stops improving by a certain amount. - improve: The amount that
the out-of-bag error needs to improve by to keep increasing the step
factor.

    model_tuned <- tuneRF(
                   x=MB_asvtab_pa[,-112], #define predictor variables
                   y=MB_asvtab_pa$sum_SFN_NIT, #define response variable
                   ntreeTry=500,
                   mtryStart=2, 
                   stepFactor=1.5,
                   improve=0.01,
                   trace=FALSE #don't show real-time progress
                   )
    ## 0.01778237 0.01 
    ## 0.02529919 0.01 
    ## 0.0302794 0.01 
    ## -0.02843199 0.01

<img src="brocc_randomforest_SFNNIT_files/figure-markdown_strict/optimize/tune the rf model 2-1.png" width="98%" height="98%" />

    # 2 predictors at each split will give best OOB

    model3 <- randomForest(formula = sum_SFN_NIT ~ . , data = MB_asvtab_pa, importance = TRUE, mtry = 2)
    model3 # not much different than train model
    ## 
    ## Call:
    ##  randomForest(formula = sum_SFN_NIT ~ ., data = MB_asvtab_pa,      importance = TRUE, mtry = 2) 
    ##                Type of random forest: regression
    ##                      Number of trees: 500
    ## No. of variables tried at each split: 2
    ## 
    ##           Mean of squared residuals: 732.7673
    ##                     % Var explained: -5.95

    # find number of trees that produced lowest test MSE
    which.min(model3$mse)
    ## [1] 298
    # 29 trees best, but this changes every time I run it

    # RMSE of best model
    sqrt(model3$mse[which.min(model3$mse)])
    ## [1] 26.73269
    # 7.8 is best mse

    # r squared
    predValues <- predict(model3, ValidSet)
    R2 <- 1 - (sum((ValidSet$sum_SFN_NIT-predValues)^2)/sum((ValidSet$sum_SFN_NIT-mean(ValidSet$sum_SFN_NIT))^2))
    # pretty mid

    #plot the test MSE by number of trees
    plot(model3)

<img src="brocc_randomforest_SFNNIT_files/figure-markdown_strict/rerun rf model with new information about best number of predictors 2-1.png" width="98%" height="98%" />


    # Get variable importance from the model fit
    ImpData <- as.data.frame(model3[["importance"]]) %>%
      rename_at('%IncMSE', ~'percIncMSE') %>%
      dplyr::filter(percIncMSE > 0)
    ImpData$Var.Names <- row.names(ImpData)

    ggplot(ImpData, aes(x=Var.Names, y=`percIncMSE`)) +
      geom_segment( aes(x=Var.Names, xend=Var.Names, y=0, yend=`percIncMSE`), color="skyblue") +
      geom_point(aes(size = IncNodePurity), color="blue", alpha=0.6) +
      theme_light() +
      coord_flip() +
      theme(
        legend.position="bottom",
        panel.grid.major.y = element_blank(),
        panel.border = element_blank(),
        axis.ticks.y = element_blank()
      )

<img src="brocc_randomforest_SFNNIT_files/figure-markdown_strict/rerun rf model with new information about best number of predictors 2-2.png" width="98%" height="98%" />


    # ggsave("import_graph_optimtestmod_presabstaxa.png",
    #        dpi = 300, height = 8,
    #        width = 8)

    # another verison of importance plot
    varImpPlot(model3, 
               sort=FALSE, 
               main="Variable Importance Plot")

<img src="brocc_randomforest_SFNNIT_files/figure-markdown_strict/rerun rf model with new information about best number of predictors 2-3.png" width="98%" height="98%" />

## random forest model to test whether microbiome grouping is more important than dietary params in determining SFN\_NIT metrics in broccoli consumers

    MB_wgroups <- read.delim("meta_groups_brocc.txt", header = T, sep = " ") %>%
      dplyr::select(-c("group_numb", "sample", "treatment", "relation")) %>%
      inner_join(SFN_NIT_tab[c("subject_id", "sum_SFN_NIT")], by = "subject_id") %>%
      column_to_rownames("subject_id")

    # Training Set : Validation Set = 70 : 30 (random)
    train <- sample(nrow(MB_wgroups), 0.7*nrow(MB_wgroups), replace = FALSE)
    TrainSet <- MB_wgroups[train,]
    ValidSet <- MB_wgroups[-train,]
    summary(TrainSet)
    ##      cohort          sex                race            ethnicity        
    ##  Min.   :1.000   Length:26          Length:26          Length:26         
    ##  1st Qu.:3.000   Class :character   Class :character   Class :character  
    ##  Median :4.000   Mode  :character   Mode  :character   Mode  :character  
    ##  Mean   :4.423                                                           
    ##  3rd Qu.:6.750                                                           
    ##  Max.   :8.000                                                           
    ##       age          height_in       weight_lb          bmi       
    ##  Min.   :18.00   Min.   :60.00   Min.   :115.0   Min.   :20.10  
    ##  1st Qu.:25.25   1st Qu.:65.00   1st Qu.:134.2   1st Qu.:21.60  
    ##  Median :28.00   Median :67.00   Median :146.7   Median :23.60  
    ##  Mean   :32.27   Mean   :66.85   Mean   :151.7   Mean   :23.85  
    ##  3rd Qu.:39.75   3rd Qu.:68.00   3rd Qu.:168.5   3rd Qu.:26.18  
    ##  Max.   :59.00   Max.   :75.00   Max.   :196.0   Max.   :28.10  
    ##      diet           sample_issue    diet_condensed        group          
    ##  Length:26          Mode :logical   Length:26          Length:26         
    ##  Class :character   FALSE:25        Class :character   Class :character  
    ##  Mode  :character   TRUE :1         Mode  :character   Mode  :character  
    ##                                                                          
    ##                                                                          
    ##                                                                          
    ##   diet_Chol_mg     diet_SatFat_g        diet_Alc_g       diet_TotSolFib_g   
    ##  Min.   :0.02535   Min.   :0.006994   Min.   :0.000000   Min.   :0.000e+00  
    ##  1st Qu.:0.09972   1st Qu.:0.010507   1st Qu.:0.000000   1st Qu.:8.245e-05  
    ##  Median :0.11716   Median :0.012615   Median :0.001890   Median :1.632e-04  
    ##  Mean   :0.12131   Mean   :0.012637   Mean   :0.002962   Mean   :2.902e-04  
    ##  3rd Qu.:0.15248   3rd Qu.:0.013946   3rd Qu.:0.004616   3rd Qu.:3.082e-04  
    ##  Max.   :0.24234   Max.   :0.019871   Max.   :0.013543   Max.   :1.182e-03  
    ##  diet_MonSac_g       diet_Disacc_g       diet_BetaCaro_mcg diet_Vit_B1_mg     
    ##  Min.   :0.0004014   Min.   :0.0004052   Min.   :0.07519   Min.   :0.0002232  
    ##  1st Qu.:0.0022233   1st Qu.:0.0011052   1st Qu.:0.40927   1st Qu.:0.0003599  
    ##  Median :0.0038450   Median :0.0016050   Median :0.74136   Median :0.0004230  
    ##  Mean   :0.0053200   Mean   :0.0025290   Mean   :1.44676   Mean   :0.0004509  
    ##  3rd Qu.:0.0058547   3rd Qu.:0.0028405   3rd Qu.:1.64145   3rd Qu.:0.0005459  
    ##  Max.   :0.0217423   Max.   :0.0147332   Max.   :7.22214   Max.   :0.0009221  
    ##  diet_Vit_B6_mg      diet_Fol.DFE_mcg_DFE  diet_Mang_mg      
    ##  Min.   :0.0002290   Min.   :0.03001      Min.   :9.916e-05  
    ##  1st Qu.:0.0003911   1st Qu.:0.10209      1st Qu.:2.793e-04  
    ##  Median :0.0004746   Median :0.13759      Median :4.157e-04  
    ##  Mean   :0.0005255   Mean   :0.15030      Mean   :6.895e-04  
    ##  3rd Qu.:0.0006274   3rd Qu.:0.17630      3rd Qu.:7.397e-04  
    ##  Max.   :0.0010802   Max.   :0.31817      Max.   :2.606e-03  
    ##   diet_Copp_mg       diet_TotFib_g       diet_OCarb_g      diet_Iron_mg     
    ##  Min.   :0.0001737   Min.   :0.005930   Min.   :0.04542   Min.   :0.004600  
    ##  1st Qu.:0.0002458   1st Qu.:0.007669   1st Qu.:0.05896   1st Qu.:0.005559  
    ##  Median :0.0003404   Median :0.009891   Median :0.06495   Median :0.006086  
    ##  Mean   :0.0003823   Mean   :0.010361   Mean   :0.06787   Mean   :0.006322  
    ##  3rd Qu.:0.0004850   3rd Qu.:0.012335   3rd Qu.:0.07276   3rd Qu.:0.006826  
    ##  Max.   :0.0007740   Max.   :0.018292   Max.   :0.12612   Max.   :0.011428  
    ##   sum_SFN_NIT     
    ##  Min.   :  9.089  
    ##  1st Qu.: 16.575  
    ##  Median : 36.169  
    ##  Mean   : 38.854  
    ##  3rd Qu.: 48.526  
    ##  Max.   :118.728
    summary(ValidSet)
    ##      cohort          sex                race            ethnicity        
    ##  Min.   :2.000   Length:12          Length:12          Length:12         
    ##  1st Qu.:3.750   Class :character   Class :character   Class :character  
    ##  Median :6.000   Mode  :character   Mode  :character   Mode  :character  
    ##  Mean   :5.333                                                           
    ##  3rd Qu.:7.000                                                           
    ##  Max.   :8.000                                                           
    ##       age          height_in       weight_lb          bmi       
    ##  Min.   :22.00   Min.   :61.00   Min.   :113.0   Min.   :20.00  
    ##  1st Qu.:25.50   1st Qu.:63.75   1st Qu.:134.2   1st Qu.:22.77  
    ##  Median :29.50   Median :67.50   Median :159.3   Median :24.60  
    ##  Mean   :35.00   Mean   :66.75   Mean   :153.5   Mean   :24.29  
    ##  3rd Qu.:44.75   3rd Qu.:69.00   3rd Qu.:166.6   3rd Qu.:25.68  
    ##  Max.   :58.00   Max.   :72.00   Max.   :200.0   Max.   :27.70  
    ##      diet           sample_issue    diet_condensed        group          
    ##  Length:12          Mode :logical   Length:12          Length:12         
    ##  Class :character   FALSE:12        Class :character   Class :character  
    ##  Mode  :character                   Mode  :character   Mode  :character  
    ##                                                                          
    ##                                                                          
    ##                                                                          
    ##   diet_Chol_mg     diet_SatFat_g        diet_Alc_g       diet_TotSolFib_g   
    ##  Min.   :0.03253   Min.   :0.008948   Min.   :0.000000   Min.   :0.000e+00  
    ##  1st Qu.:0.10672   1st Qu.:0.012093   1st Qu.:0.000000   1st Qu.:1.624e-05  
    ##  Median :0.13846   Median :0.014206   Median :0.000000   Median :1.483e-04  
    ##  Mean   :0.15660   Mean   :0.014804   Mean   :0.002400   Mean   :2.879e-04  
    ##  3rd Qu.:0.18629   3rd Qu.:0.016396   3rd Qu.:0.003165   3rd Qu.:4.531e-04  
    ##  Max.   :0.39274   Max.   :0.026479   Max.   :0.011471   Max.   :1.022e-03  
    ##  diet_MonSac_g       diet_Disacc_g       diet_BetaCaro_mcg diet_Vit_B1_mg     
    ##  Min.   :0.0003911   Min.   :8.932e-05   Min.   :0.09368   Min.   :0.0001685  
    ##  1st Qu.:0.0012335   1st Qu.:9.640e-04   1st Qu.:0.19316   1st Qu.:0.0002630  
    ##  Median :0.0022568   Median :1.295e-03   Median :0.98601   Median :0.0003388  
    ##  Mean   :0.0029231   Mean   :2.334e-03   Mean   :1.53058   Mean   :0.0004979  
    ##  3rd Qu.:0.0038292   3rd Qu.:2.153e-03   3rd Qu.:1.85476   3rd Qu.:0.0004671  
    ##  Max.   :0.0072088   Max.   :9.533e-03   Max.   :6.41777   Max.   :0.0018019  
    ##  diet_Vit_B6_mg      diet_Fol.DFE_mcg_DFE  diet_Mang_mg      
    ##  Min.   :0.0003271   Min.   :0.04229      Min.   :6.507e-05  
    ##  1st Qu.:0.0004196   1st Qu.:0.09003      1st Qu.:1.823e-04  
    ##  Median :0.0004916   Median :0.10561      Median :4.398e-04  
    ##  Mean   :0.0005555   Mean   :0.12523      Mean   :5.538e-04  
    ##  3rd Qu.:0.0006153   3rd Qu.:0.13227      3rd Qu.:6.990e-04  
    ##  Max.   :0.0011872   Max.   :0.29975      Max.   :2.284e-03  
    ##   diet_Copp_mg       diet_TotFib_g       diet_OCarb_g      diet_Iron_mg     
    ##  Min.   :0.0001474   Min.   :0.004311   Min.   :0.03178   Min.   :0.004710  
    ##  1st Qu.:0.0002127   1st Qu.:0.006980   1st Qu.:0.04939   1st Qu.:0.005109  
    ##  Median :0.0002644   Median :0.008980   Median :0.05203   Median :0.005507  
    ##  Mean   :0.0003244   Mean   :0.010128   Mean   :0.05569   Mean   :0.005972  
    ##  3rd Qu.:0.0003471   3rd Qu.:0.012238   3rd Qu.:0.06159   3rd Qu.:0.005794  
    ##  Max.   :0.0009170   Max.   :0.018270   Max.   :0.08983   Max.   :0.010683  
    ##   sum_SFN_NIT    
    ##  Min.   : 5.079  
    ##  1st Qu.:23.158  
    ##  Median :29.613  
    ##  Mean   :35.594  
    ##  3rd Qu.:38.529  
    ##  Max.   :87.613

    # Create a Random Forest model with default parameters
    model1 <- randomForest(sum_SFN_NIT ~ ., data = TrainSet, importance = TRUE)
    model1
    ## 
    ## Call:
    ##  randomForest(formula = sum_SFN_NIT ~ ., data = TrainSet, importance = TRUE) 
    ##                Type of random forest: regression
    ##                      Number of trees: 500
    ## No. of variables tried at each split: 9
    ## 
    ##           Mean of squared residuals: 771.5091
    ##                     % Var explained: 0.49

    R2 <- 1 - (sum((ValidSet$sum_SFN_NIT-predValues)^2)/sum((ValidSet$sum_SFN_NIT-mean(ValidSet$sum_SFN_NIT))^2))
    # negative, not a great fit for data

    # find number of trees that produced lowest test MSE
    which.min(model1$mse)
    ## [1] 34
    # 1

    # RMSE of best model
    sqrt(model1$mse[which.min(model1$mse)])
    ## [1] 27.27936
    # 7.1

    model2 <- randomForest(formula = sum_SFN_NIT ~ . , data = MB_wgroups, importance = TRUE)
    model2 # not much different than train model
    ## 
    ## Call:
    ##  randomForest(formula = sum_SFN_NIT ~ ., data = MB_wgroups, importance = TRUE) 
    ##                Type of random forest: regression
    ##                      Number of trees: 500
    ## No. of variables tried at each split: 9
    ## 
    ##           Mean of squared residuals: 732.7818
    ##                     % Var explained: -5.96

    # find number of trees that produced lowest test MSE
    which.min(model2$mse)
    ## [1] 67
    # 136 trees best

    # RMSE of best model
    sqrt(model1$mse[which.min(model2$mse)])
    ## [1] 28.39737
    # 9.7 is best mse

    #plot the test MSE by number of trees
    plot(model2)

<img src="brocc_randomforest_SFNNIT_files/figure-markdown_strict/run full model 3-1.png" width="98%" height="98%" />

optimize model: - ntreeTry: The number of trees to build. - mtryStart:
The starting number of predictor variables to consider at each split. -
stepFactor: The factor to increase by until the out-of-bag estimated
error stops improving by a certain amount. - improve: The amount that
the out-of-bag error needs to improve by to keep increasing the step
factor.

    model_tuned <- tuneRF(
                   x=MB_wgroups[,-1], #define predictor variables
                   y=MB_wgroups$sum_SFN_NIT, #define response variable
                   ntreeTry=500,
                   mtryStart=2, 
                   stepFactor=1.5,
                   improve=0.01,
                   trace=FALSE #don't show real-time progress
                   )
    ## 0.08759485 0.01 
    ## 0.08921163 0.01 
    ## 0.2295429 0.01 
    ## 0.2200969 0.01 
    ## 0.309115 0.01 
    ## 0.4296999 0.01 
    ## 0.473212 0.01

<img src="brocc_randomforest_SFNNIT_files/figure-markdown_strict/optimize/tune the rf model 3-1.png" width="98%" height="98%" />

    # 27 predictors at each split will give best OOB

    model3 <- randomForest(formula = sum_SFN_NIT ~ . , data = MB_wgroups, importance = TRUE, mtry = 27)
    model3 # not much different than train model
    ## 
    ## Call:
    ##  randomForest(formula = sum_SFN_NIT ~ ., data = MB_wgroups, importance = TRUE,      mtry = 27) 
    ##                Type of random forest: regression
    ##                      Number of trees: 500
    ## No. of variables tried at each split: 27
    ## 
    ##           Mean of squared residuals: 755.901
    ##                     % Var explained: -9.3

    # find number of trees that produced lowest test MSE
    which.min(model3$mse)
    ## [1] 346
    # 1 tree best, but this changes every time I run it

    # RMSE of best model
    sqrt(model3$mse[which.min(model3$mse)])
    ## [1] 27.41026
    # 5.5 is best mse

    # r squared
    predValues <- predict(model3, ValidSet)
    R2 <- 1 - (sum((ValidSet$sum_SFN_NIT-predValues)^2)/sum((ValidSet$sum_SFN_NIT-mean(ValidSet$sum_SFN_NIT))^2))
    # pretty high

    #plot the test MSE by number of trees
    plot(model3)

<img src="brocc_randomforest_SFNNIT_files/figure-markdown_strict/rerun rf model with new information about best number of predictors 3-1.png" width="98%" height="98%" />


    # Get variable importance from the model fit
    ImpData <- as.data.frame(model3[["importance"]]) %>%
      rename_at('%IncMSE', ~'percIncMSE') %>%
      dplyr::filter(percIncMSE > 0)
    ImpData$Var.Names <- row.names(ImpData)

    ggplot(ImpData, aes(x=Var.Names, y=`percIncMSE`)) +
      geom_segment( aes(x=Var.Names, xend=Var.Names, y=0, yend=`percIncMSE`), color="skyblue") +
      geom_point(aes(size = IncNodePurity), color="blue", alpha=0.6) +
      theme_light() +
      coord_flip() +
      theme(
        legend.position="bottom",
        panel.grid.major.y = element_blank(),
        panel.border = element_blank(),
        axis.ticks.y = element_blank()
      )

<img src="brocc_randomforest_SFNNIT_files/figure-markdown_strict/rerun rf model with new information about best number of predictors 3-2.png" width="98%" height="98%" />


    # ggsave("importgraph_optimtestmod_AllvarsMBgroupSFN_NIT.png",
    #        dpi = 300, height = 4,
    #        width = 5)

    # another verison of importance plot
    varImpPlot(model3, 
               sort=FALSE, 
               main="Variable Importance Plot")

<img src="brocc_randomforest_SFNNIT_files/figure-markdown_strict/rerun rf model with new information about best number of predictors 3-3.png" width="98%" height="98%" />

## save environment

    save.image("brocc_randomforestSFNNIT.RData")
