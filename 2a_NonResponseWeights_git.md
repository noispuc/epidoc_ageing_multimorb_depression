1 Correction for non-response using inverted probability weighting
================
Tomoe Gusberti, Dr. Eng
2025-03-20

``` r
load(file='./Data/Data_cleanedTreated_wide.RData')#data_wide,demog,

data_wide$NonResp=is.na(data_wide$Idade_T3)
data_wide=left_join(data_wide,demog)
```

    ## Joining with `by = join_by(ID, Escolaridade)`

``` r
data_wide$NUTII <- factor(data_wide$NUTII, 
                            levels = c(1:7),
                            labels=c('Norte',  'Centro', 'Lisboa', 'Alentejo', 'Algarve', 'Açores', 'Madeira'))


data_wide$SEXO <- factor(data_wide$SEXO, levels = c(1,2),labels=c('Male','Female'))
data_wide$Estado_Civil <- factor(data_wide$Estado_Civil, 
                            levels = c(1:5),
                            labels=c('Single', 'Married', 'Divorced',
                                     'Widowed', 'in a stable union'
))

data_wide$Escolaridade <- factor(data_wide$Escolaridade, 
                            levels = c(1:5),
                            labels=c('>12 years (university level)', 
                                     '10 ~ 12 years study)',
                                     '9 years study (elementary level)',
                                     '4 years study (primary level)',
                                     '< 4 years study'

))
data_wide$EmploymentStatus_T0<-factor(data_wide$EmploymentStatus_T0,labels=(c('Active Worker','unemployed','retired','other')))


data_wide$MMorbCat_T3<-factor(data_wide$MMorbCat_T3,levels=c('0','1','2','3','+3'))
data_wide$MMorbCat0<-factor(data_wide$MMorbCat0,levels=c('0','1','2','3','+3'))
data_wide$DepressaoCat_T3 <- factor(data_wide$DepressaoCat_T3, labels = c('normal','mild','moderate',"severe"))
data_wide$Depressao0Cat <- factor(data_wide$Depressao0Cat, labels = c('normal','mild','moderate',"severe"))

data_wide$familySize0 <- factor(data_wide$familySize0, levels = c('LivingAlone','SmallFamily','MediumFamily',"LargeFamily"))
data_wide$FamilyMembers0<-as.numeric(data_wide$FamilyMembers0)
data_wide$FamilyMembers_T3<-as.numeric(data_wide$FamilyMembers_T3)

data_wide$ageCat0 <- factor(data_wide$ageCat0, levels = c('Youth_<40','Adult_40_49','Adult_50_59',"Adult_60_69",'Adult_70_79','Adult_80+'))
save(data_wide,demog,file='./Data/Data_cleanedTreated_wide2.RData')#data_wide,demog,
```

\#checking balance

``` r
library('cobalt')
```

    ##  cobalt (Version 4.6.0, Build Date: 2025-04-15)

``` r
bal.plot(NonResp ~ Depressao0Cat+ageCat0+SEXO+familySize0+EmploymentStatus_T0,
        data = data_wide, estimand = "ATE", thresholds = c(m = .05))
```

    ## Warning: Missing values exist in the covariates. Displayed values omit these
    ## observations.

    ## No `var.name` was provided. Displaying balance for Depressao0Cat.

![](2a_NonResponseWeights_git_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
ggsave('./Inspect_balance/balPlot.png')
```

    ## Saving 7 x 5 in image

``` r
love.plot(NonResp ~ Depressao0Cat+ageCat0+SEXO+familySize0+EmploymentStatus_T0,
        data = data_wide, estimand = "ATE", thresholds = c(m = .05))
```

    ## Warning: Missing values exist in the covariates. Displayed values omit these
    ## observations.

![](2a_NonResponseWeights_git_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
ggsave('./Inspect_balance/lovePlot.png')
```

    ## Saving 7 x 5 in image

# weights usign WieghtIt

``` r
library("WeightIt")
W.out <- weightit(NonResp ~ Depressao0Cat+ageCat0+SEXO+MMorbCat0+familySize0+Smoke0+Alcohol0+Exercise0+EmploymentStatus_T0+Escolaridade,
        data = data_wide, estimand = "ATE", method = "ipt")
```

    ## Warning: Missing values are present in the covariates. See
    ## `?WeightIt::method_ipt` for information on how these are handled.

    ## Warning in stode(y, times, func, parms = parms, ...): error during
    ## factorisation of matrix (dgefa); singular matrix

    ## Warning in stode(y, times, func, parms = parms, ...): steady-state not reached

    ## Warning: The optimization failed to converge; consider using fewer covariates or
    ## a different link function.

``` r
W.out #print the output
```

    ## A weightit object
    ##  - method: "ipt" (inverse probability tilting)
    ##  - number of obs.: 10661
    ##  - sampling weights: none
    ##  - treatment: 2-category
    ##  - estimand: ATE
    ##  - covariates: Depressao0Cat, ageCat0, SEXO, MMorbCat0, familySize0, Smoke0, Alcohol0, Exercise0, EmploymentStatus_T0, Escolaridade
    ##  - missingness method: missingness indicators

``` r
summary(W.out)
```

    ##                   Summary of weights
    ## 
    ## - Weight ranges:
    ## 
    ##            Min                                   Max
    ## treated 1.0000 |--|                           2.8405
    ## control 1.5391  |--------------------------| 15.8213
    ## 
    ## - Units with the 5 most extreme weights by group:
    ##                                                 
    ##            10538    5320    5771    1615      31
    ##  treated  2.6888  2.6888  2.6996   2.766  2.8405
    ##             1135    2608    9267    4799    8798
    ##  control 11.6032 11.8012 13.1892 14.1603 15.8213
    ## 
    ## - Weight statistics:
    ## 
    ##         Coef of Var   MAD Entropy # Zeros
    ## treated       0.173 0.134   0.014       0
    ## control       0.399 0.244   0.061       0
    ## 
    ## - Effective Sample Sizes:
    ## 
    ##            Control Treated
    ## Unweighted 3757.   6904.  
    ## Weighted   3241.74 6704.02

``` r
data_wide$W=W.out[["weights"]]
data_wide$ps=W.out[["ps"]]
data_wide%>%group_by(NonResp,Depressao0Cat)%>%summarize(mean(W),sd(W))
```

    ## `summarise()` has grouped output by 'NonResp'. You can override using the
    ## `.groups` argument.

    ## # A tibble: 8 × 4
    ## # Groups:   NonResp [2]
    ##   NonResp Depressao0Cat `mean(W)` `sd(W)`
    ##   <lgl>   <fct>             <dbl>   <dbl>
    ## 1 FALSE   normal             2.74   0.965
    ## 2 FALSE   mild               3.36   1.67 
    ## 3 FALSE   moderate           3.26   1.70 
    ## 4 FALSE   severe             3.89   1.64 
    ## 5 TRUE    normal             1.58   0.266
    ## 6 TRUE    mild               1.43   0.226
    ## 7 TRUE    moderate           1.44   0.246
    ## 8 TRUE    severe             1.31   0.191

# by logit glm

``` r
library(broom)
# Logit model to predict net use
model <- glm(NonResp ~ Depressao0Cat+ageCat0+SEXO+MMorbCat0+familySize0+Smoke0+Alcohol0+Exercise0+EmploymentStatus_T0+Escolaridade,
                         family = binomial(link = "logit"),
                         data = data_wide)

# Generate propensity scores and IPWs
data_ipw <- augment_columns(model, data_wide,
                                type.predict = "response") %>% 
  rename(propensity = .fitted) %>% 
  mutate(ipw = (NonResp / propensity) + ((1 - NonResp) / (1 - propensity)))
```

``` r
summary(model)
```

    ## 
    ## Call:
    ## glm(formula = NonResp ~ Depressao0Cat + ageCat0 + SEXO + MMorbCat0 + 
    ##     familySize0 + Smoke0 + Alcohol0 + Exercise0 + EmploymentStatus_T0 + 
    ##     Escolaridade, family = binomial(link = "logit"), data = data_wide)
    ## 
    ## Coefficients:
    ##                                              Estimate Std. Error z value
    ## (Intercept)                                   0.39515    0.25772   1.533
    ## Depressao0Catmild                             0.21527    0.20229   1.064
    ## Depressao0Catmoderate                        -0.19126    0.29348  -0.652
    ## Depressao0Catsevere                           0.25752    0.45823   0.562
    ## ageCat0Adult_40_49                           -0.34407    0.12606  -2.729
    ## ageCat0Adult_50_59                           -0.29753    0.15378  -1.935
    ## ageCat0Adult_60_69                           -0.46071    0.23568  -1.955
    ## ageCat0Adult_70_79                           -0.09253    0.37609  -0.246
    ## ageCat0Adult_80+                              1.34537    1.06160   1.267
    ## SEXOFemale                                   -0.10263    0.10964  -0.936
    ## MMorbCat01                                   -0.21493    0.11941  -1.800
    ## MMorbCat02                                   -0.34152    0.16800  -2.033
    ## MMorbCat03                                   -0.30796    0.25699  -1.198
    ## MMorbCat0+3                                  -0.15652    0.36913  -0.424
    ## familySize0SmallFamily                       -0.13967    0.14435  -0.968
    ## familySize0MediumFamily                      -0.05985    0.16268  -0.368
    ## familySize0LargeFamily                        0.39689    0.58172   0.682
    ## Smoke0Yes                                     0.24510    0.14741   1.663
    ## Alcohol0Ocasionally                          -0.03779    0.12316  -0.307
    ## Alcohol0Daily                                -0.15809    0.15021  -1.052
    ## Exercise0TRUE                                 0.07883    0.10733   0.734
    ## EmploymentStatus_T0unemployed                 0.18468    0.15838   1.166
    ## EmploymentStatus_T0retired                    0.47491    0.21873   2.171
    ## EmploymentStatus_T0other                      0.44051    0.13740   3.206
    ## Escolaridade10 ~ 12 years study)              0.29538    0.14021   2.107
    ## Escolaridade9 years study (elementary level)  0.69431    0.14695   4.725
    ## Escolaridade4 years study (primary level)     0.59199    0.16699   3.545
    ## Escolaridade< 4 years study                   0.99917    0.35787   2.792
    ##                                              Pr(>|z|)    
    ## (Intercept)                                  0.125219    
    ## Depressao0Catmild                            0.287259    
    ## Depressao0Catmoderate                        0.514603    
    ## Depressao0Catsevere                          0.574118    
    ## ageCat0Adult_40_49                           0.006344 ** 
    ## ageCat0Adult_50_59                           0.053014 .  
    ## ageCat0Adult_60_69                           0.050602 .  
    ## ageCat0Adult_70_79                           0.805660    
    ## ageCat0Adult_80+                             0.205048    
    ## SEXOFemale                                   0.349237    
    ## MMorbCat01                                   0.071877 .  
    ## MMorbCat02                                   0.042066 *  
    ## MMorbCat03                                   0.230790    
    ## MMorbCat0+3                                  0.671543    
    ## familySize0SmallFamily                       0.333250    
    ## familySize0MediumFamily                      0.712948    
    ## familySize0LargeFamily                       0.495067    
    ## Smoke0Yes                                    0.096361 .  
    ## Alcohol0Ocasionally                          0.758988    
    ## Alcohol0Daily                                0.292603    
    ## Exercise0TRUE                                0.462662    
    ## EmploymentStatus_T0unemployed                0.243597    
    ## EmploymentStatus_T0retired                   0.029917 *  
    ## EmploymentStatus_T0other                     0.001346 ** 
    ## Escolaridade10 ~ 12 years study)             0.035147 *  
    ## Escolaridade9 years study (elementary level)  2.3e-06 ***
    ## Escolaridade4 years study (primary level)    0.000392 ***
    ## Escolaridade< 4 years study                  0.005238 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 2623.6  on 2087  degrees of freedom
    ## Residual deviance: 2536.0  on 2060  degrees of freedom
    ##   (8573 observations deleted due to missingness)
    ## AIC: 2592
    ## 
    ## Number of Fisher Scoring iterations: 5

``` r
data_ipw%>%group_by(NonResp,Depressao0Cat)%>%summarize(mean(ipw),sd(ipw),mean(Idade_T3))
```

    ## `summarise()` has grouped output by 'NonResp'. You can override using the
    ## `.groups` argument.

    ## # A tibble: 8 × 5
    ## # Groups:   NonResp [2]
    ##   NonResp Depressao0Cat `mean(ipw)` `sd(ipw)` `mean(Idade_T3)`
    ##   <lgl>   <fct>               <dbl>     <dbl>            <dbl>
    ## 1 FALSE   normal               3.06     1.03              51.6
    ## 2 FALSE   mild                 3.72     1.11              59.4
    ## 3 FALSE   moderate             2.93     0.845             54.5
    ## 4 FALSE   severe               3.24     0.324             55.4
    ## 5 TRUE    normal               1.49     0.223             NA  
    ## 6 TRUE    mild                 1.37     0.167             NA  
    ## 7 TRUE    moderate             1.50     0.259             NA  
    ## 8 TRUE    severe               1.33     0.161             NA

``` r
data_bal<-data_ipw%>%filter(NonResp==FALSE)
save(data_bal,file='./Data/Data_ipwBalancedT3_wide.RData')
```

``` r
devtools::session_info()
```

    ## ─ Session info ───────────────────────────────────────────────────────────────
    ##  setting  value
    ##  version  R version 4.5.2 (2025-10-31)
    ##  os       Ubuntu 24.04.4 LTS
    ##  system   x86_64, linux-gnu
    ##  ui       X11
    ##  language (EN)
    ##  collate  en_US.UTF-8
    ##  ctype    en_US.UTF-8
    ##  tz       America/Sao_Paulo
    ##  date     2026-04-08
    ##  pandoc   3.6.3 @ /usr/lib/rstudio/resources/app/bin/quarto/bin/tools/x86_64/ (via rmarkdown)
    ##  quarto   1.8.25 @ /usr/lib/rstudio/resources/app/bin/quarto/bin/quarto
    ## 
    ## ─ Packages ───────────────────────────────────────────────────────────────────
    ##  ! package      * version date (UTC) lib source
    ##  P backports      1.5.0   2024-05-23 [?] RSPM
    ##  P broom        * 1.0.9   2025-07-28 [?] RSPM
    ##  P cachem         1.1.0   2024-05-16 [?] CRAN (R 4.5.2)
    ##  P chk            0.10.0  2025-01-24 [?] CRAN (R 4.5.2)
    ##  P cli            3.6.5   2025-04-23 [?] CRAN (R 4.5.2)
    ##  P cobalt       * 4.6.0   2025-04-15 [?] CRAN (R 4.5.2)
    ##  P crayon         1.5.3   2024-06-20 [?] CRAN (R 4.5.2)
    ##  P devtools       2.5.0   2026-03-14 [?] CRAN (R 4.5.2)
    ##  P digest         0.6.37  2024-08-19 [?] RSPM
    ##  P dplyr        * 1.1.4   2023-11-17 [?] CRAN (R 4.5.2)
    ##  P ellipsis       0.3.2   2021-04-29 [?] RSPM
    ##  P evaluate       1.0.4   2025-06-18 [?] RSPM
    ##  P farver         2.1.2   2024-05-13 [?] CRAN (R 4.5.2)
    ##  P fastmap        1.2.0   2024-05-15 [?] CRAN (R 4.5.2)
    ##  P fs             2.0.1   2026-03-24 [?] RSPM
    ##  P generics       0.1.4   2025-05-09 [?] CRAN (R 4.5.2)
    ##  P ggplot2      * 4.0.0   2025-09-11 [?] RSPM
    ##  P glue           1.8.0   2024-09-30 [?] CRAN (R 4.5.2)
    ##  P gtable         0.3.6   2024-10-25 [?] CRAN (R 4.5.2)
    ##  P htmltools      0.5.8.1 2024-04-04 [?] RSPM
    ##  P knitr          1.50    2025-03-16 [?] RSPM
    ##  P labeling       0.4.3   2023-08-29 [?] CRAN (R 4.5.2)
    ##  P lattice        0.22-6  2024-03-20 [?] RSPM
    ##  P lifecycle      1.0.5   2026-01-08 [?] CRAN (R 4.5.2)
    ##  P magrittr       2.0.3   2022-03-30 [?] RSPM
    ##  P memoise        2.0.1   2021-11-26 [?] CRAN (R 4.5.2)
    ##  P pillar         1.11.0  2025-07-04 [?] RSPM
    ##  P pkgbuild       1.4.8   2025-05-26 [?] CRAN (R 4.5.2)
    ##  P pkgconfig      2.0.3   2019-09-22 [?] CRAN (R 4.5.2)
    ##  P pkgload        1.5.1   2026-04-01 [?] CRAN (R 4.5.2)
    ##  P purrr          1.2.1   2026-01-09 [?] CRAN (R 4.5.2)
    ##  P R6             2.6.1   2025-02-15 [?] CRAN (R 4.5.2)
    ##  P ragg           1.5.2   2026-03-23 [?] CRAN (R 4.5.2)
    ##  P RColorBrewer   1.1-3   2022-04-03 [?] CRAN (R 4.5.2)
    ##    renv           1.1.5   2025-07-24 [1] RSPM
    ##  P rlang          1.2.0   2026-04-06 [?] RSPM
    ##  P rmarkdown      2.29    2024-11-04 [?] RSPM
    ##  P rootSolve      1.8.2.4 2023-09-21 [?] CRAN (R 4.5.2)
    ##  P rstudioapi     0.17.1  2024-10-22 [?] RSPM
    ##  P S7             0.2.0   2024-11-07 [?] RSPM
    ##  P sandwich       3.1-1   2024-09-15 [?] CRAN (R 4.5.2)
    ##  P scales         1.4.0   2025-04-24 [?] CRAN (R 4.5.2)
    ##  P sessioninfo    1.2.3   2025-02-05 [?] CRAN (R 4.5.2)
    ##  P systemfonts    1.3.1   2025-10-01 [?] RSPM
    ##  P textshaping    1.0.4   2025-10-10 [?] RSPM
    ##  P tibble         3.3.0   2025-06-08 [?] RSPM
    ##  P tidyr          1.3.1   2024-01-24 [?] RSPM
    ##  P tidyselect     1.2.1   2024-03-11 [?] CRAN (R 4.5.2)
    ##  P usethis        3.2.1   2025-09-06 [?] CRAN (R 4.5.2)
    ##  P utf8           1.2.6   2025-06-08 [?] CRAN (R 4.5.2)
    ##  P vctrs          0.7.2   2026-03-21 [?] CRAN (R 4.5.2)
    ##  P WeightIt     * 1.4.0   2025-02-24 [?] RSPM
    ##  P withr          3.0.2   2024-10-28 [?] CRAN (R 4.5.2)
    ##  P xfun           0.53    2025-08-19 [?] CRAN (R 4.5.2)
    ##  P yaml           2.3.10  2024-07-26 [?] RSPM
    ##  P zoo            1.8-14  2025-04-10 [?] RSPM
    ## 
    ##  [1] /mnt/Data/Projects/RutePortugal/renv/library/linux-ubuntu-noble/R-4.5/x86_64-pc-linux-gnu
    ##  [2] /home/tomoe/.cache/R/renv/sandbox/linux-ubuntu-noble/R-4.5/x86_64-pc-linux-gnu/9a444a72
    ## 
    ##  * ── Packages attached to the search path.
    ##  P ── Loaded and on-disk path mismatch.
    ## 
    ## ──────────────────────────────────────────────────────────────────────────────
