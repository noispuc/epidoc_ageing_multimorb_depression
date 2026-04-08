2 Cross-Lagged Structural Equation analysis
================
Tomoe Gusberti, Dr. Eng
2025-03-20

# PREPARE DATA

- depression levels- from imputed (multiple) data
- multimorbidity from original scores?
- inverted propensity weight - from ipwBalanced (imputation does not
  change it)

``` r
load(file='./Data/Data_Imp20.RData')
```

``` r
load(file='./Data/Data_ipwBalancedT3_wide.RData')
data_bal<-data_bal%>%mutate(depression=ifelse(Depressao0Cat=='normal','non-Depressive',
                                      ifelse(Depressao0Cat%in%c('mild',
                                      'moderate','severe'),'Depressive',NA)))%>%
    rename('Dep_EpiDoc1' := Depressao_T0, 'Dep_EpiDoc2':= Depressao_T1,
           'Dep_EpiDoc4':= Depressao_T3,
           'Age':=Idade_T0)
```

## merging and creating as mids

``` r
library(miceadds)
```

    ## Loading required package: mice

    ## 
    ## Attaching package: 'mice'

    ## The following object is masked from 'package:stats':
    ## 
    ##     filter

    ## The following objects are masked from 'package:base':
    ## 
    ##     cbind, rbind

    ## * miceadds 3.17-44 (2024-01-08 19:08:24)

``` r
#datlist2mids

data_mi=list()
for (i in 1:imp20$m){
  dt<-complete(imp20,i)
  dt<-left_join(data_bal%>%
                  select(ID,ipw,depression,MMorb_T0,MMorb_T1,MMorb_T3,MMorb0)
                ,dt) %>%
    mutate(
      MMorbidity=ifelse(MMorb0>=2,'>=2 diseases','<2 diseases'),
      EducationalLevel=ifelse(Escolaridade>2,
                              '<=9 years','>9years')
      )%>%
    rename('Dep_EpiDoc1' := Pontuacao_Depressao_T0, 'Dep_EpiDoc2':= Pontuacao_Depressao_T1,
           'Dep_EpiDoc4':= Pontuacao_Depressao_T3,
           'MM_EpiDoc1':=MMorb_T0,
           'MM_EpiDoc2':=MMorb_T1,
           'MM_EpiDoc4':=MMorb_T3,
           'Age':=Idade_T0,
           'Education':=Escolaridade)%>%mutate(
             Sex=ifelse(SEXO==1,'Male','Female')
           )
  # prep interaction terms
  dt<-dt%>%mutate(
    MM2_ed=scale(Education*MM_EpiDoc2)%>%as.numeric(),
    MM1_ed=scale(Education*MM_EpiDoc1)%>%as.numeric(),
    d2_ed=scale(Education*Dep_EpiDoc2)%>%as.numeric(),
    d1_ed=scale(Education*Dep_EpiDoc1)%>%as.numeric(),
    age_ed=scale(Education*Age)%>%as.numeric()
  )
  data_mi[[i]]<-dt
}
```

    ## Joining with `by = join_by(ID)`

    ## Joining with `by = join_by(ID)`
    ## Joining with `by = join_by(ID)`
    ## Joining with `by = join_by(ID)`
    ## Joining with `by = join_by(ID)`
    ## Joining with `by = join_by(ID)`
    ## Joining with `by = join_by(ID)`
    ## Joining with `by = join_by(ID)`
    ## Joining with `by = join_by(ID)`
    ## Joining with `by = join_by(ID)`
    ## Joining with `by = join_by(ID)`
    ## Joining with `by = join_by(ID)`
    ## Joining with `by = join_by(ID)`
    ## Joining with `by = join_by(ID)`
    ## Joining with `by = join_by(ID)`
    ## Joining with `by = join_by(ID)`
    ## Joining with `by = join_by(ID)`
    ## Joining with `by = join_by(ID)`
    ## Joining with `by = join_by(ID)`
    ## Joining with `by = join_by(ID)`

``` r
data.mi=datlist2mids(data_mi)
```

    ## Warning: Number of logged events: 80

``` r
save(data.mi,file='./Data/Data_Imp20_forSEM.RData')
```

# Lagged SEM model

## ondas 0 e 1

``` r
model01<-'
Dep_EpiDoc2~MM_EpiDoc1+Dep_EpiDoc1
MM_EpiDoc2~MM_EpiDoc1+Dep_EpiDoc1
MM_EpiDoc1~~Dep_EpiDoc1
MM_EpiDoc2~~Dep_EpiDoc2
'
res<-lavaan.mi::sem.mi(data=data.mi,std.ov=TRUE,
         model=model01,
         estimator='MLR',likelihood = "wishart",
           sampling.weights='ipw',
         orthogonal=TRUE
         )
summary(res,fit.meas=TRUE)
```

    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.

    ## lavaan.mi object fit to 20 imputed data sets using:
    ##  - lavaan    (0.6-19)
    ##  - lavaan.mi (0.1-0)
    ## See class?lavaan.mi help page for available methods. 
    ## 
    ## Convergence information:
    ## The model converged on 20 imputed data sets.
    ## Standard errors were available for all imputations.
    ## 
    ##   Estimator                                         ML
    ##   Optimization method                           NLMINB
    ##   Number of model parameters                        10
    ## 
    ##   Number of observations                           672
    ##   Sampling weights variable                        ipw
    ## 
    ## Model Test User Model:
    ## 
    ##                                                     Standard      Scaled
    ##   Test statistic                                       0.000          NA
    ##   Degrees of freedom                                       0           0
    ##   P-value                                              1.000          NA
    ##   Average scaling correction factor                                   NA
    ##   Pooling method                                          D4            
    ##     Pooled statistic                              "standard"            
    ##     "yuan.bentler.mplus" correction applied            AFTER     pooling
    ## 
    ## Model Test Baseline Model:
    ## 
    ##   Test statistic                               308.174     192.692
    ##   Degrees of freedom                                 6           6
    ##   P-value                                        0.000       0.000
    ##   Scaling correction factor                                  1.599
    ## 
    ## User Model versus Baseline Model:
    ## 
    ##   Comparative Fit Index (CFI)                    1.000          NA
    ##   Tucker-Lewis Index (TLI)                       1.000          NA
    ##                                                                   
    ##   Robust Comparative Fit Index (CFI)                            NA
    ##   Robust Tucker-Lewis Index (TLI)                               NA
    ## 
    ## Loglikelihood and Information Criteria:
    ## 
    ##   Loglikelihood user model (H0)              -3617.151   -3617.151
    ##   Loglikelihood unrestricted model (H1)      -3617.165   -3617.165
    ##                                                                   
    ##   Akaike (AIC)                                7254.301    7254.301
    ##   Bayesian (BIC)                              7299.404    7299.404
    ##   Sample-size adjusted Bayesian (SABIC)       7267.653    7267.653
    ## 
    ## Root Mean Square Error of Approximation:
    ## 
    ##   RMSEA                                          0.000          NA
    ##   90 Percent confidence interval - lower         0.000          NA
    ##   90 Percent confidence interval - upper         0.000          NA
    ##   P-value H_0: RMSEA <= 0.050                       NA          NA
    ##   P-value H_0: RMSEA >= 0.080                       NA          NA
    ##                                                                   
    ##   Robust RMSEA                                               0.000
    ##   90 Percent confidence interval - lower                        NA
    ##   90 Percent confidence interval - upper                        NA
    ##   P-value H_0: Robust RMSEA <= 0.050                            NA
    ##   P-value H_0: Robust RMSEA >= 0.080                            NA
    ## 
    ## Standardized Root Mean Square Residual:
    ## 
    ##   SRMR                                           0.013       0.013
    ## 
    ## Parameter Estimates:
    ## 
    ##   Standard errors                                    Sandwich
    ##   Information bread                                  Observed
    ##   Observed information based on                       Hessian
    ##                                                              
    ##   Pooled across imputations              Rubin's (1987) rules
    ##   Augment within-imputation variance     Scale by average RIV
    ##   Wald test for pooled parameters          t(df) distribution
    ## 
    ##   Pooled t statistics with df >= 1000 are displayed with
    ##   df = Inf(inity) to save space. Although the t distribution
    ##   with large df closely approximates a standard normal
    ##   distribution, exact df for reporting these t tests can be
    ##   obtained from parameterEstimates.mi() 
    ## 
    ## 
    ## Regressions:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|)
    ##   Dep_EpiDoc2 ~                                                
    ##     MM_EpiDoc1        0.088    0.047    1.872  746.271    0.062
    ##     Dep_EpiDoc1       0.437    0.056    7.854  598.426    0.000
    ##   MM_EpiDoc2 ~                                                 
    ##     MM_EpiDoc1        0.334    0.077    4.368      Inf    0.000
    ##     Dep_EpiDoc1       0.093    0.049    1.914      Inf    0.056
    ## 
    ## Covariances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|)
    ##   MM_EpiDoc1 ~~                                                
    ##     Dep_EpiDoc1       0.356    0.070    5.097      Inf    0.000
    ##  .Dep_EpiDoc2 ~~                                               
    ##    .MM_EpiDoc2        0.018    0.039    0.459  335.630    0.647
    ## 
    ## Variances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|)
    ##    .Dep_EpiDoc2       0.792    0.076   10.378  429.862    0.000
    ##    .MM_EpiDoc2        0.741    0.103    7.214      Inf    0.000
    ##     MM_EpiDoc1        1.018    0.101   10.061      Inf    0.000
    ##     Dep_EpiDoc1       1.034    0.095   10.848      Inf    0.000

## ondas 1, 2 e 3

``` r
model013<-'
# evolução
Dep_EpiDoc4~Dep_EpiDoc1+Dep_EpiDoc2
MM_EpiDoc4~MM_EpiDoc2+MM_EpiDoc1
Dep_EpiDoc2~Dep_EpiDoc1
MM_EpiDoc2~MM_EpiDoc1

# causalidade - c/ temporalidade/lag
Dep_EpiDoc4~MM_EpiDoc2+MM_EpiDoc1
MM_EpiDoc4~Dep_EpiDoc2+Dep_EpiDoc1
Dep_EpiDoc2~MM_EpiDoc1
MM_EpiDoc2~Dep_EpiDoc1

# corrrelação mesmo t
MM_EpiDoc1~~Dep_EpiDoc1
MM_EpiDoc2~~Dep_EpiDoc2
MM_EpiDoc4~~Dep_EpiDoc4
'
res<-lavaan.mi::sem.mi(data=data.mi,std.ov=TRUE,
         model=model013,
         estimator='MLR',likelihood = "wishart", 
           sampling.weights='ipw',
         orthogonal=TRUE
         )
summary(res,fit.meas=TRUE)
```

    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.

    ## lavaan.mi object fit to 20 imputed data sets using:
    ##  - lavaan    (0.6-19)
    ##  - lavaan.mi (0.1-0)
    ## See class?lavaan.mi help page for available methods. 
    ## 
    ## Convergence information:
    ## The model converged on 20 imputed data sets.
    ## Standard errors were available for all imputations.
    ## 
    ##   Estimator                                         ML
    ##   Optimization method                           NLMINB
    ##   Number of model parameters                        21
    ## 
    ##   Number of observations                           672
    ##   Sampling weights variable                        ipw
    ## 
    ## Model Test User Model:
    ## 
    ##                                                     Standard      Scaled
    ##   Test statistic                                       0.000          NA
    ##   Degrees of freedom                                       0           0
    ##   P-value                                              1.000          NA
    ##   Average scaling correction factor                                   NA
    ##   Pooling method                                          D4            
    ##     Pooled statistic                              "standard"            
    ##     "yuan.bentler.mplus" correction applied            AFTER     pooling
    ## 
    ## Model Test Baseline Model:
    ## 
    ##   Test statistic                               722.112     500.728
    ##   Degrees of freedom                                15          15
    ##   P-value                                        0.000       0.000
    ##   Scaling correction factor                                  1.442
    ## 
    ## User Model versus Baseline Model:
    ## 
    ##   Comparative Fit Index (CFI)                    1.000          NA
    ##   Tucker-Lewis Index (TLI)                       1.000          NA
    ##                                                                   
    ##   Robust Comparative Fit Index (CFI)                            NA
    ##   Robust Tucker-Lewis Index (TLI)                               NA
    ## 
    ## Loglikelihood and Information Criteria:
    ## 
    ##   Loglikelihood user model (H0)              -5292.824   -5292.824
    ##   Loglikelihood unrestricted model (H1)      -5292.838   -5292.838
    ##                                                                   
    ##   Akaike (AIC)                               10627.648   10627.648
    ##   Bayesian (BIC)                             10722.363   10722.363
    ##   Sample-size adjusted Bayesian (SABIC)      10655.686   10655.686
    ## 
    ## Root Mean Square Error of Approximation:
    ## 
    ##   RMSEA                                          0.000          NA
    ##   90 Percent confidence interval - lower         0.000          NA
    ##   90 Percent confidence interval - upper         0.000          NA
    ##   P-value H_0: RMSEA <= 0.050                       NA          NA
    ##   P-value H_0: RMSEA >= 0.080                       NA          NA
    ##                                                                   
    ##   Robust RMSEA                                               0.000
    ##   90 Percent confidence interval - lower                        NA
    ##   90 Percent confidence interval - upper                        NA
    ##   P-value H_0: Robust RMSEA <= 0.050                            NA
    ##   P-value H_0: Robust RMSEA >= 0.080                            NA
    ## 
    ## Standardized Root Mean Square Residual:
    ## 
    ##   SRMR                                           0.013       0.013
    ## 
    ## Parameter Estimates:
    ## 
    ##   Standard errors                                    Sandwich
    ##   Information bread                                  Observed
    ##   Observed information based on                       Hessian
    ##                                                              
    ##   Pooled across imputations              Rubin's (1987) rules
    ##   Augment within-imputation variance     Scale by average RIV
    ##   Wald test for pooled parameters          t(df) distribution
    ## 
    ##   Pooled t statistics with df >= 1000 are displayed with
    ##   df = Inf(inity) to save space. Although the t distribution
    ##   with large df closely approximates a standard normal
    ##   distribution, exact df for reporting these t tests can be
    ##   obtained from parameterEstimates.mi() 
    ## 
    ## 
    ## Regressions:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|)
    ##   Dep_EpiDoc4 ~                                                
    ##     Dep_EpiDoc1       0.280    0.057    4.879      Inf    0.000
    ##     Dep_EpiDoc2       0.242    0.057    4.243  659.070    0.000
    ##   MM_EpiDoc4 ~                                                 
    ##     MM_EpiDoc2        0.133    0.049    2.690      Inf    0.007
    ##     MM_EpiDoc1        0.470    0.045   10.558      Inf    0.000
    ##   Dep_EpiDoc2 ~                                                
    ##     Dep_EpiDoc1       0.437    0.056    7.839  598.426    0.000
    ##   MM_EpiDoc2 ~                                                 
    ##     MM_EpiDoc1        0.334    0.077    4.360      Inf    0.000
    ##   Dep_EpiDoc4 ~                                                
    ##     MM_EpiDoc2        0.113    0.044    2.570      Inf    0.010
    ##     MM_EpiDoc1       -0.029    0.046   -0.629  885.577    0.529
    ##   MM_EpiDoc4 ~                                                 
    ##     Dep_EpiDoc2       0.070    0.042    1.656      Inf    0.098
    ##     Dep_EpiDoc1       0.055    0.040    1.378      Inf    0.168
    ##   Dep_EpiDoc2 ~                                                
    ##     MM_EpiDoc1        0.088    0.047    1.869  746.271    0.062
    ##   MM_EpiDoc2 ~                                                 
    ##     Dep_EpiDoc1       0.093    0.049    1.910      Inf    0.056
    ## 
    ## Covariances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|)
    ##   Dep_EpiDoc1 ~~                                               
    ##     MM_EpiDoc1        0.356    0.070    5.088      Inf    0.000
    ##  .Dep_EpiDoc2 ~~                                               
    ##    .MM_EpiDoc2        0.018    0.039    0.458  335.630    0.647
    ##  .Dep_EpiDoc4 ~~                                               
    ##    .MM_EpiDoc4        0.039    0.032    1.224  537.259    0.222
    ## 
    ## Variances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|)
    ##    .Dep_EpiDoc4       0.786    0.070   11.168      Inf    0.000
    ##    .MM_EpiDoc4        0.635    0.046   13.924      Inf    0.000
    ##    .Dep_EpiDoc2       0.792    0.076   10.359  429.862    0.000
    ##    .MM_EpiDoc2        0.741    0.103    7.201      Inf    0.000
    ##     Dep_EpiDoc1       1.034    0.096   10.828      Inf    0.000
    ##     MM_EpiDoc1        1.018    0.101   10.043      Inf    0.000

# adding age

``` r
modelAll0<-'
# evolução
Dep_EpiDoc4~dd03*Dep_EpiDoc1+dd13*Dep_EpiDoc2
MM_EpiDoc4~MM_EpiDoc2+MM_EpiDoc1
Dep_EpiDoc2~dd01*Dep_EpiDoc1
MM_EpiDoc2~MM_EpiDoc1

# causalidade - c/ temporalidade/lag
Dep_EpiDoc4~MM_EpiDoc2+MM_EpiDoc1
#Dep_EpiDoc4~MM_EpiDoc2
MM_EpiDoc4~Dep_EpiDoc2+Dep_EpiDoc1
Dep_EpiDoc2~dm1*MM_EpiDoc1
MM_EpiDoc2~md1*Dep_EpiDoc1

# corrrelação mesmo t
MM_EpiDoc1~~Dep_EpiDoc1
MM_EpiDoc2~~Dep_EpiDoc2
MM_EpiDoc4~~Dep_EpiDoc4

# causalidade - idade
Dep_EpiDoc4~d3*Age
MM_EpiDoc4~m3*Age
Dep_EpiDoc2~d1*Age
MM_EpiDoc2~m1*Age
Dep_EpiDoc1~d0*Age
MM_EpiDoc1~m0*Age

age_dep_indirect:=d0*dd01*dd13+d0*dd03
age_dep_direct:=d3

dep0_direct:=dd03
dep0_indirect:=dd01*dd13
dep_total:=dd03+dd13*dd01

# hypotheses 3
age_dep1_direct:=d1
age_dep1_m_indirect:=m0*dm1

'
resAll0<-lavaan.mi::sem.mi(data=data.mi,
         model=modelAll0,
         estimator='MLR',likelihood = "wishart", 
           sampling.weights='ipw',
         orthogonal=TRUE
         )
summary(resAll0,fit.meas=TRUE,standardized=TRUE,rsquare=TRUE,ci=TRUE)
```

    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.

    ## lavaan.mi object fit to 20 imputed data sets using:
    ##  - lavaan    (0.6-19)
    ##  - lavaan.mi (0.1-0)
    ## See class?lavaan.mi help page for available methods. 
    ## 
    ## Convergence information:
    ## The model converged on 20 imputed data sets.
    ## Standard errors were available for all imputations.
    ## 
    ##   Estimator                                         ML
    ##   Optimization method                           NLMINB
    ##   Number of model parameters                        27
    ## 
    ##   Number of observations                           672
    ##   Sampling weights variable                        ipw
    ## 
    ## Model Test User Model:
    ## 
    ##                                                     Standard      Scaled
    ##   Test statistic                                       0.000          NA
    ##   Degrees of freedom                                       0           0
    ##   P-value                                              1.000          NA
    ##   Average scaling correction factor                                   NA
    ##   Pooling method                                          D4            
    ##     Pooled statistic                              "standard"            
    ##     "yuan.bentler.mplus" correction applied            AFTER     pooling
    ## 
    ## Model Test Baseline Model:
    ## 
    ##   Test statistic                               965.949     712.936
    ##   Degrees of freedom                                21          21
    ##   P-value                                        0.000       0.000
    ##   Scaling correction factor                                  1.355
    ## 
    ## User Model versus Baseline Model:
    ## 
    ##   Comparative Fit Index (CFI)                    1.000          NA
    ##   Tucker-Lewis Index (TLI)                       1.000          NA
    ##                                                                   
    ##   Robust Comparative Fit Index (CFI)                            NA
    ##   Robust Tucker-Lewis Index (TLI)                               NA
    ## 
    ## Loglikelihood and Information Criteria:
    ## 
    ##   Loglikelihood user model (H0)              -7426.648   -7426.648
    ##   Loglikelihood unrestricted model (H1)      -7426.697   -7426.697
    ##                                                                   
    ##   Akaike (AIC)                               14907.296   14907.296
    ##   Bayesian (BIC)                             15029.073   15029.073
    ##   Sample-size adjusted Bayesian (SABIC)      14943.346   14943.346
    ## 
    ## Root Mean Square Error of Approximation:
    ## 
    ##   RMSEA                                          0.000          NA
    ##   90 Percent confidence interval - lower         0.000          NA
    ##   90 Percent confidence interval - upper         0.000          NA
    ##   P-value H_0: RMSEA <= 0.050                       NA          NA
    ##   P-value H_0: RMSEA >= 0.080                       NA          NA
    ##                                                                   
    ##   Robust RMSEA                                               0.000
    ##   90 Percent confidence interval - lower                        NA
    ##   90 Percent confidence interval - upper                        NA
    ##   P-value H_0: Robust RMSEA <= 0.050                            NA
    ##   P-value H_0: Robust RMSEA >= 0.080                            NA
    ## 
    ## Standardized Root Mean Square Residual:
    ## 
    ##   SRMR                                           0.017       0.017
    ## 
    ## Parameter Estimates:
    ## 
    ##   Standard errors                                    Sandwich
    ##   Information bread                                  Observed
    ##   Observed information based on                       Hessian
    ##                                                              
    ##   Pooled across imputations              Rubin's (1987) rules
    ##   Augment within-imputation variance     Scale by average RIV
    ##   Wald test for pooled parameters          t(df) distribution
    ## 
    ##   Pooled t statistics with df >= 1000 are displayed with
    ##   df = Inf(inity) to save space. Although the t distribution
    ##   with large df closely approximates a standard normal
    ##   distribution, exact df for reporting these t tests can be
    ##   obtained from parameterEstimates.mi() 
    ## 
    ## 
    ## Regressions:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##   Dep_EpiDoc4 ~                                                         
    ##     Dp_EpD1 (dd03)    0.265    0.056    4.768      Inf    0.000    0.156
    ##     Dp_EpD2 (dd13)    0.230    0.056    4.077  817.933    0.000    0.119
    ##   MM_EpiDoc4 ~                                                          
    ##     MM_EpD2           0.295    0.096    3.071      Inf    0.002    0.107
    ##     MM_EpD1           0.475    0.063    7.598      Inf    0.000    0.353
    ##   Dep_EpiDoc2 ~                                                         
    ##     Dp_EpD1 (dd01)    0.419    0.055    7.582  382.591    0.000    0.311
    ##   MM_EpiDoc2 ~                                                          
    ##     MM_EpD1           0.229    0.055    4.181      Inf    0.000    0.122
    ##   Dep_EpiDoc4 ~                                                         
    ##     MM_EpD2           0.588    0.224    2.626      Inf    0.009    0.149
    ##     MM_EpD1          -0.184    0.170   -1.084  943.034    0.279   -0.518
    ##   MM_EpiDoc4 ~                                                          
    ##     Dp_EpD2           0.019    0.016    1.145  392.238    0.253   -0.013
    ##     Dp_EpD1           0.019    0.016    1.137      Inf    0.256   -0.013
    ##   Dep_EpiDoc2 ~                                                         
    ##     MM_EpD1  (dm1)    0.123    0.180    0.686  389.573    0.493   -0.230
    ##   MM_EpiDoc2 ~                                                          
    ##     Dp_EpD1  (md1)    0.018    0.009    1.938      Inf    0.053   -0.000
    ##   Dep_EpiDoc4 ~                                                         
    ##     Age       (d3)    0.014    0.012    1.194  199.715    0.234   -0.009
    ##   MM_EpiDoc4 ~                                                          
    ##     Age       (m3)    0.022    0.006    3.784      Inf    0.000    0.010
    ##   Dep_EpiDoc2 ~                                                         
    ##     Age       (d1)    0.027    0.013    1.977  122.712    0.050   -0.000
    ##   MM_EpiDoc2 ~                                                          
    ##     Age       (m1)   -0.002    0.002   -1.165      Inf    0.244   -0.007
    ##   Dep_EpiDoc1 ~                                                         
    ##     Age       (d0)    0.061    0.014    4.250      Inf    0.000    0.033
    ##   MM_EpiDoc1 ~                                                          
    ##     Age       (m0)    0.038    0.004    9.557      Inf    0.000    0.030
    ##  ci.upper   Std.lv  Std.all
    ##                            
    ##     0.374    0.265    0.280
    ##     0.341    0.230    0.237
    ##                            
    ##     0.483    0.295    0.138
    ##     0.598    0.475    0.372
    ##                            
    ##     0.528    0.419    0.431
    ##                            
    ##     0.337    0.229    0.384
    ##                            
    ##     1.027    0.588    0.109
    ##     0.150   -0.184   -0.057
    ##                            
    ##     0.051    0.019    0.049
    ##     0.050    0.019    0.050
    ##                            
    ##     0.477    0.123    0.037
    ##                            
    ##     0.037    0.018    0.105
    ##                            
    ##     0.038    0.014    0.059
    ##                            
    ##     0.033    0.022    0.225
    ##                            
    ##     0.053    0.027    0.106
    ##                            
    ##     0.002   -0.002   -0.054
    ##                            
    ##     0.089    0.061    0.237
    ##                            
    ##     0.045    0.038    0.497
    ## 
    ## Covariances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##  .Dep_EpiDoc1 ~~                                                        
    ##    .MM_EpiDoc1        0.782    0.199    3.922      Inf    0.000    0.391
    ##  .Dep_EpiDoc2 ~~                                                        
    ##    .MM_EpiDoc2        0.046    0.082    0.565  305.320    0.572   -0.115
    ##  .Dep_EpiDoc4 ~~                                                        
    ##    .MM_EpiDoc4        0.120    0.128    0.935  430.056    0.350   -0.132
    ##  ci.upper   Std.lv  Std.all
    ##                            
    ##     1.173    0.782    0.272
    ##                            
    ##     0.208    0.046    0.029
    ##                            
    ##     0.372    0.120    0.042
    ## 
    ## Variances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##    .Dep_EpiDoc4       8.044    0.723   11.122      Inf    0.000    6.625
    ##    .MM_EpiDoc4        1.001    0.072   13.811      Inf    0.000    0.859
    ##    .Dep_EpiDoc2       8.429    0.800   10.536  192.641    0.000    6.851
    ##    .MM_EpiDoc2        0.298    0.042    7.123      Inf    0.000    0.216
    ##    .Dep_EpiDoc1      11.012    1.074   10.258      Inf    0.000    8.908
    ##    .MM_EpiDoc1        0.753    0.076    9.963      Inf    0.000    0.605
    ##  ci.upper   Std.lv  Std.all
    ##     9.462    8.044    0.772
    ##     1.143    1.001    0.614
    ##    10.007    8.429    0.765
    ##     0.380    0.298    0.834
    ##    13.117   11.012    0.944
    ##     0.901    0.753    0.753
    ## 
    ## R-Square:
    ##                    Estimate
    ##     Dep_EpiDoc4       0.228
    ##     MM_EpiDoc4        0.386
    ##     Dep_EpiDoc2       0.235
    ##     MM_EpiDoc2        0.166
    ##     Dep_EpiDoc1       0.056
    ##     MM_EpiDoc1        0.247
    ## 
    ## Defined Parameters:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##     age_dep_indrct    0.022    0.006    3.680      Inf    0.000    0.010
    ##     age_dep_direct    0.014    0.012    1.194  199.715    0.234   -0.009
    ##     dep0_direct       0.265    0.056    4.768      Inf    0.000    0.156
    ##     dep0_indirect     0.097    0.028    3.513  522.799    0.000    0.043
    ##     dep_total         0.362    0.053    6.881      Inf    0.000    0.259
    ##     age_dep1_dirct    0.027    0.013    1.977  122.712    0.050   -0.000
    ##     ag_dp1_m_ndrct    0.005    0.007    0.689  382.968    0.491   -0.009
    ##  ci.upper   Std.lv  Std.all
    ##     0.034    0.022    0.091
    ##     0.038    0.014    0.059
    ##     0.374    0.265    0.280
    ##     0.151    0.097    0.102
    ##     0.465    0.362    0.383
    ##     0.053    0.027    0.106
    ##     0.018    0.005    0.018

``` r
source('./code/parTableToCSV.R')
parTableToCSV(resAll0,path='./Results/LaggedSEM_mi20/g_SEMED')
```

    ## $partable
    ##                  label     est    se      t           df pvalue ci.lower
    ## 1                 dd03   0.265 0.056  4.768 1.681122e+03  0.000    0.156
    ## 2                 dd13   0.230 0.056  4.077 8.179330e+02  0.000    0.119
    ## 3                        0.295 0.096  3.071 1.415974e+07  0.002    0.107
    ## 4                        0.475 0.063  7.598 9.071777e+06  0.000    0.353
    ## 5                 dd01   0.419 0.055  7.582 3.825910e+02  0.000    0.311
    ## 6                        0.229 0.055  4.181 9.422294e+36  0.000    0.122
    ## 7                        0.588 0.224  2.626 1.365779e+03  0.009    0.149
    ## 8                       -0.184 0.170 -1.084 9.430340e+02  0.279   -0.518
    ## 9                        0.019 0.016  1.145 3.922380e+02  0.253   -0.013
    ## 10                       0.019 0.016  1.137 7.331101e+03  0.256   -0.013
    ## 11                 dm1   0.123 0.180  0.686 3.895730e+02  0.493   -0.230
    ## 12                 md1   0.018 0.009  1.938 3.272508e+29  0.053    0.000
    ## 13                       0.782 0.199  3.922 3.954019e+29  0.000    0.391
    ## 14                       0.046 0.082  0.565 3.053200e+02  0.572   -0.115
    ## 15                       0.120 0.128  0.935 4.300560e+02  0.350   -0.132
    ## 16                  d3   0.014 0.012  1.194 1.997150e+02  0.234   -0.009
    ## 17                  m3   0.022 0.006  3.784 1.900788e+07  0.000    0.010
    ## 18                  d1   0.027 0.013  1.977 1.227120e+02  0.050    0.000
    ## 19                  m1  -0.002 0.002 -1.165 1.151641e+29  0.244   -0.007
    ## 20                  d0   0.061 0.014  4.250 9.794115e+28  0.000    0.033
    ## 21                  m0   0.038 0.004  9.557 2.535496e+28  0.000    0.030
    ## 22                       8.044 0.723 11.122 1.531384e+03  0.000    6.625
    ## 23                       1.001 0.072 13.811 1.456266e+07  0.000    0.859
    ## 24                       8.429 0.800 10.536 1.926410e+02  0.000    6.851
    ## 25                       0.298 0.042  7.123 5.806690e+30  0.000    0.216
    ## 26                      11.012 1.074 10.258 1.414457e+30  0.000    8.908
    ## 27                       0.753 0.076  9.963 4.403520e+29  0.000    0.605
    ## 28                     175.024 0.000     NA           NA     NA  175.024
    ## 29    age_dep_indirect   0.022 0.006  3.680 7.452287e+04  0.000    0.010
    ## 30      age_dep_direct   0.014 0.012  1.194 1.997150e+02  0.234   -0.009
    ## 31         dep0_direct   0.265 0.056  4.768 1.681122e+03  0.000    0.156
    ## 32       dep0_indirect   0.097 0.028  3.513 5.227990e+02  0.000    0.043
    ## 33           dep_total   0.362 0.053  6.881 7.142625e+03  0.000    0.259
    ## 34     age_dep1_direct   0.027 0.013  1.977 1.227120e+02  0.050    0.000
    ## 35 age_dep1_m_indirect   0.005 0.007  0.689 3.829680e+02  0.491   -0.009
    ##    ci.upper  std.lv std.all std.nox                                 label2
    ## 1     0.374   0.265   0.280   0.280                Dep_EpiDoc4~Dep_EpiDoc1
    ## 2     0.341   0.230   0.237   0.237                Dep_EpiDoc4~Dep_EpiDoc2
    ## 3     0.483   0.295   0.138   0.138                  MM_EpiDoc4~MM_EpiDoc2
    ## 4     0.598   0.475   0.372   0.372                  MM_EpiDoc4~MM_EpiDoc1
    ## 5     0.528   0.419   0.431   0.431                Dep_EpiDoc2~Dep_EpiDoc1
    ## 6     0.337   0.229   0.384   0.384                  MM_EpiDoc2~MM_EpiDoc1
    ## 7     1.027   0.588   0.109   0.109                 Dep_EpiDoc4~MM_EpiDoc2
    ## 8     0.150  -0.184  -0.057  -0.057                 Dep_EpiDoc4~MM_EpiDoc1
    ## 9     0.051   0.019   0.049   0.049                 MM_EpiDoc4~Dep_EpiDoc2
    ## 10    0.050   0.019   0.050   0.050                 MM_EpiDoc4~Dep_EpiDoc1
    ## 11    0.477   0.123   0.037   0.037                 Dep_EpiDoc2~MM_EpiDoc1
    ## 12    0.037   0.018   0.105   0.105                 MM_EpiDoc2~Dep_EpiDoc1
    ## 13    1.173   0.782   0.272   0.272                Dep_EpiDoc1~~MM_EpiDoc1
    ## 14    0.208   0.046   0.029   0.029                Dep_EpiDoc2~~MM_EpiDoc2
    ## 15    0.372   0.120   0.042   0.042                Dep_EpiDoc4~~MM_EpiDoc4
    ## 16    0.038   0.014   0.059   0.004                        Dep_EpiDoc4~Age
    ## 17    0.033   0.022   0.225   0.017                         MM_EpiDoc4~Age
    ## 18    0.053   0.027   0.106   0.008                        Dep_EpiDoc2~Age
    ## 19    0.002  -0.002  -0.054  -0.004                         MM_EpiDoc2~Age
    ## 20    0.089   0.061   0.237   0.018                        Dep_EpiDoc1~Age
    ## 21    0.045   0.038   0.497   0.038                         MM_EpiDoc1~Age
    ## 22    9.462   8.044   0.772   0.772               Dep_EpiDoc4~~Dep_EpiDoc4
    ## 23    1.143   1.001   0.614   0.614                 MM_EpiDoc4~~MM_EpiDoc4
    ## 24   10.007   8.429   0.765   0.765               Dep_EpiDoc2~~Dep_EpiDoc2
    ## 25    0.380   0.298   0.834   0.834                 MM_EpiDoc2~~MM_EpiDoc2
    ## 26   13.117  11.012   0.944   0.944               Dep_EpiDoc1~~Dep_EpiDoc1
    ## 27    0.901   0.753   0.753   0.753                 MM_EpiDoc1~~MM_EpiDoc1
    ## 28  175.024 175.024   1.000 175.024                               Age~~Age
    ## 29    0.034   0.022   0.091   0.007 age_dep_indirect:=d0*dd01*dd13+d0*dd03
    ## 30    0.038   0.014   0.059   0.004                     age_dep_direct:=d3
    ## 31    0.374   0.265   0.280   0.280                      dep0_direct:=dd03
    ## 32    0.151   0.097   0.102   0.102               dep0_indirect:=dd01*dd13
    ## 33    0.465   0.362   0.383   0.383              dep_total:=dd03+dd13*dd01
    ## 34    0.053   0.027   0.106   0.008                    age_dep1_direct:=d1
    ## 35    0.018   0.005   0.018   0.001            age_dep1_m_indirect:=m0*dm1
    ##        estimate_p                CI
    ## 1       0.265 (0)     [0.156:0.374]
    ## 2        0.23 (0)     [0.119:0.341]
    ## 3   0.295 (0.002)     [0.107:0.483]
    ## 4       0.475 (0)     [0.353:0.598]
    ## 5       0.419 (0)     [0.311:0.528]
    ## 6       0.229 (0)     [0.122:0.337]
    ## 7   0.588 (0.009)     [0.149:1.027]
    ## 8  -0.184 (0.279)     [-0.518:0.15]
    ## 9   0.019 (0.253)    [-0.013:0.051]
    ## 10  0.019 (0.256)     [-0.013:0.05]
    ## 11  0.123 (0.493)     [-0.23:0.477]
    ## 12  0.018 (0.053)         [0:0.037]
    ## 13      0.782 (0)     [0.391:1.173]
    ## 14  0.046 (0.572)    [-0.115:0.208]
    ## 15    0.12 (0.35)    [-0.132:0.372]
    ## 16  0.014 (0.234)    [-0.009:0.038]
    ## 17      0.022 (0)      [0.01:0.033]
    ## 18   0.027 (0.05)         [0:0.053]
    ## 19 -0.002 (0.244)    [-0.007:0.002]
    ## 20      0.061 (0)     [0.033:0.089]
    ## 21      0.038 (0)      [0.03:0.045]
    ## 22      8.044 (0)     [6.625:9.462]
    ## 23      1.001 (0)     [0.859:1.143]
    ## 24      8.429 (0)    [6.851:10.007]
    ## 25      0.298 (0)      [0.216:0.38]
    ## 26     11.012 (0)    [8.908:13.117]
    ## 27      0.753 (0)     [0.605:0.901]
    ## 28   175.024 (NA) [175.024:175.024]
    ## 29      0.022 (0)      [0.01:0.034]
    ## 30  0.014 (0.234)    [-0.009:0.038]
    ## 31      0.265 (0)     [0.156:0.374]
    ## 32      0.097 (0)     [0.043:0.151]
    ## 33      0.362 (0)     [0.259:0.465]
    ## 34   0.027 (0.05)         [0:0.053]
    ## 35  0.005 (0.491)    [-0.009:0.018]

# adjusted and education

``` r
modelAll0_adj<-'
# evolução
Dep_EpiDoc4~dd03*Dep_EpiDoc1+dd13*Dep_EpiDoc2
MM_EpiDoc4~mm13*MM_EpiDoc2+mm03*MM_EpiDoc1
Dep_EpiDoc2~dd01*Dep_EpiDoc1
MM_EpiDoc2~mm01*MM_EpiDoc1

# causalidade - c/ temporalidade/lag
Dep_EpiDoc4~md13*MM_EpiDoc2+md03*MM_EpiDoc1
#Dep_EpiDoc4~md13*MM_EpiDoc2
MM_EpiDoc4~dm13*Dep_EpiDoc2+dm03*Dep_EpiDoc1
Dep_EpiDoc2~dm01*MM_EpiDoc1
MM_EpiDoc2~md01*Dep_EpiDoc1

# corrrelação mesmo t
MM_EpiDoc1~~Dep_EpiDoc1
MM_EpiDoc2~~Dep_EpiDoc2
MM_EpiDoc4~~Dep_EpiDoc4

# causalidade - idade
Dep_EpiDoc4~d3*Age
MM_EpiDoc4~m3*Age
Dep_EpiDoc2~d1*Age
MM_EpiDoc2~m1*Age
Dep_EpiDoc1~d0*Age
MM_EpiDoc1~m0*Age

# educação
Dep_EpiDoc4~de3*Education
MM_EpiDoc4~me3*Education
Dep_EpiDoc2~de1*Education
MM_EpiDoc2~me1*Education
Dep_EpiDoc1~de0*Education
MM_EpiDoc1~me0*Education

#age effects
age_dep_indirect:=d0*dd01*dd13+d0*dd03+d1*dd13
age_depM_indirect:=m0*md03+m1*md13+d0*dm01*md13+d0*dm01*md13
age_direct:=d3
age_total:=age_dep_indirect+age_depM_indirect


## indirect effect of educaiton on dep4
ed_dep_indirect:=de0*dd01*dd13+de1*dd13
ed_depM_indirect:=de1*dm01*md13
ed_direct:=de3

ed_total:=ed_dep_indirect+ed_depM_indirect

## cumulative effect d1 on dep4
dep0_direct:=dd03
dep0_indirect:=dd01*dd13+dd01*dd13
dep0_m_indirect:=dm01*md13+dm01*md13
dep_total:=dep0_direct+dep0_indirect+dep0_m_indirect

##indirect effect MM1 on dep4
m0_direct:=md03
m0_indirect:=md01*dd13+mm01*md13
m0_total:=m0_direct+m0_indirect


'
res_adj<-lavaan.mi::sem.mi(data=data.mi,
         model=modelAll0_adj,
         estimator = "ML", test = "yuan.bentler.mplus",
         likelihood = "wishart", missing = "FIML", 
           sampling.weights='ipw',
         orthogonal=TRUE
         )
summary(res_adj,fit.meas=TRUE)
```

    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.

    ## Warning: lavaan->lav_lavaan_step11_estoptim():  
    ##    Model estimation FAILED! Returning starting values.
    ## Warning: lavaan->lav_lavaan_step11_estoptim():  
    ##    Model estimation FAILED! Returning starting values.

    ## lavaan.mi object fit to 20 imputed data sets using:
    ##  - lavaan    (0.6-19)
    ##  - lavaan.mi (0.1-0)
    ## See class?lavaan.mi help page for available methods. 
    ## 
    ## Convergence information:
    ## The model converged on 20 imputed data sets.
    ## Standard errors were available for all imputations.
    ## 
    ##   Estimator                                         ML
    ##   Optimization method                           NLMINB
    ##   Number of model parameters                        39
    ## 
    ##   Number of observations                           672
    ##   Number of missing patterns                         1
    ##   Sampling weights variable                        ipw
    ## 
    ## Model Test User Model:
    ## 
    ##                                                     Standard      Scaled
    ##   Test statistic                                       0.000          NA
    ##   Degrees of freedom                                       0           0
    ##   P-value                                              1.000          NA
    ##   Average scaling correction factor                                   NA
    ##   Pooling method                                          D4            
    ##     Pooled statistic                              "standard"            
    ##     "yuan.bentler.mplus" correction applied            AFTER     pooling
    ## 
    ## Model Test Baseline Model:
    ## 
    ##   Test statistic                              1026.436     799.284
    ##   Degrees of freedom                                27          27
    ##   P-value                                        0.000       0.000
    ##   Scaling correction factor                                  1.284
    ## 
    ## User Model versus Baseline Model:
    ## 
    ##   Comparative Fit Index (CFI)                    1.000          NA
    ##   Tucker-Lewis Index (TLI)                       1.000          NA
    ##                                                                   
    ##   Robust Comparative Fit Index (CFI)                         0.992
    ##   Robust Tucker-Lewis Index (TLI)                            1.000
    ## 
    ## Loglikelihood and Information Criteria:
    ## 
    ##   Loglikelihood user model (H0)              -7421.263   -7421.263
    ##   Loglikelihood unrestricted model (H1)      -7421.323   -7421.323
    ##                                                                   
    ##   Akaike (AIC)                               14920.526   14920.526
    ##   Bayesian (BIC)                             15096.426   15096.426
    ##   Sample-size adjusted Bayesian (SABIC)      14972.598   14972.598
    ## 
    ## Root Mean Square Error of Approximation:
    ## 
    ##   RMSEA                                          0.000          NA
    ##   90 Percent confidence interval - lower         0.000          NA
    ##   90 Percent confidence interval - upper         0.000          NA
    ##   P-value H_0: RMSEA <= 0.050                       NA          NA
    ##   P-value H_0: RMSEA >= 0.080                       NA          NA
    ##                                                                   
    ##   Robust RMSEA                                               0.000
    ##   90 Percent confidence interval - lower                     0.000
    ##   90 Percent confidence interval - upper                     0.000
    ##   P-value H_0: Robust RMSEA <= 0.050                            NA
    ##   P-value H_0: Robust RMSEA >= 0.080                            NA
    ## 
    ## Standardized Root Mean Square Residual:
    ## 
    ##   SRMR                                           0.014       0.014
    ## 
    ## Parameter Estimates:
    ## 
    ##   Standard errors                                    Sandwich
    ##   Information bread                                  Observed
    ##   Observed information based on                       Hessian
    ##                                                              
    ##   Pooled across imputations              Rubin's (1987) rules
    ##   Augment within-imputation variance     Scale by average RIV
    ##   Wald test for pooled parameters          t(df) distribution
    ## 
    ##   Pooled t statistics with df >= 1000 are displayed with
    ##   df = Inf(inity) to save space. Although the t distribution
    ##   with large df closely approximates a standard normal
    ##   distribution, exact df for reporting these t tests can be
    ##   obtained from parameterEstimates.mi() 
    ## 
    ## 
    ## Regressions:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|)
    ##   Dep_EpiDoc4 ~                                                
    ##     Dp_EpD1 (dd03)    0.261    0.055    4.769      Inf    0.000
    ##     Dp_EpD2 (dd13)    0.225    0.056    3.981  742.959    0.000
    ##   MM_EpiDoc4 ~                                                 
    ##     MM_EpD2 (mm13)    0.295    0.095    3.108      Inf    0.002
    ##     MM_EpD1 (mm03)    0.475    0.062    7.667      Inf    0.000
    ##   Dep_EpiDoc2 ~                                                
    ##     Dp_EpD1 (dd01)    0.398    0.055    7.185  349.103    0.000
    ##   MM_EpiDoc2 ~                                                 
    ##     MM_EpD1 (mm01)    0.228    0.054    4.207      Inf    0.000
    ##   Dep_EpiDoc4 ~                                                
    ##     MM_EpD2 (md13)    0.598    0.222    2.691      Inf    0.007
    ##     MM_EpD1 (md03)   -0.182    0.168   -1.078  947.239    0.281
    ##   MM_EpiDoc4 ~                                                 
    ##     Dp_EpD2 (dm13)    0.019    0.016    1.150  363.038    0.251
    ##     Dp_EpD1 (dm03)    0.018    0.016    1.130      Inf    0.258
    ##   Dep_EpiDoc2 ~                                                
    ##     MM_EpD1 (dm01)    0.139    0.174    0.801  362.394    0.424
    ##   MM_EpiDoc2 ~                                                 
    ##     Dp_EpD1 (md01)    0.020    0.009    2.093      Inf    0.036
    ##   Dep_EpiDoc4 ~                                                
    ##     Age       (d3)    0.011    0.012    0.919  244.178    0.359
    ##   MM_EpiDoc4 ~                                                 
    ##     Age       (m3)    0.022    0.006    3.828      Inf    0.000
    ##   Dep_EpiDoc2 ~                                                
    ##     Age       (d1)    0.014    0.013    1.095  121.526    0.276
    ##   MM_EpiDoc2 ~                                                 
    ##     Age       (m1)   -0.002    0.002   -0.743      Inf    0.458
    ##   Dep_EpiDoc1 ~                                                
    ##     Age       (d0)    0.045    0.015    2.974      Inf    0.003
    ##   MM_EpiDoc1 ~                                                 
    ##     Age       (m0)    0.037    0.004    9.553      Inf    0.000
    ##   Dep_EpiDoc4 ~                                                
    ##     Educatn  (de3)    0.112    0.113    0.989      Inf    0.323
    ##   MM_EpiDoc4 ~                                                 
    ##     Educatn  (me3)    0.003    0.041    0.081      Inf    0.935
    ##   Dep_EpiDoc2 ~                                                
    ##     Educatn  (de1)    0.407    0.123    3.316  306.762    0.001
    ##   MM_EpiDoc2 ~                                                 
    ##     Educatn  (me1)   -0.025    0.020   -1.231      Inf    0.218
    ##   Dep_EpiDoc1 ~                                                
    ##     Educatn  (de0)    0.513    0.139    3.693      Inf    0.000
    ##   MM_EpiDoc1 ~                                                 
    ##     Educatn  (me0)    0.011    0.036    0.306      Inf    0.760
    ## 
    ## Covariances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|)
    ##  .Dep_EpiDoc1 ~~                                               
    ##    .MM_EpiDoc1        0.776    0.193    4.016      Inf    0.000
    ##  .Dep_EpiDoc2 ~~                                               
    ##    .MM_EpiDoc2        0.057    0.080    0.711  276.239    0.478
    ##  .Dep_EpiDoc4 ~~                                               
    ##    .MM_EpiDoc4        0.120    0.127    0.944  426.291    0.345
    ## 
    ## Intercepts:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|)
    ##    .Dep_EpiDoc4       0.709    0.468    1.516  288.267    0.131
    ##    .MM_EpiDoc4       -0.388    0.210   -1.842      Inf    0.065
    ##    .Dep_EpiDoc2       0.018    0.516    0.035  175.169    0.972
    ##    .MM_EpiDoc2        0.096    0.071    1.348      Inf    0.178
    ##    .Dep_EpiDoc1      -0.220    0.598   -0.368      Inf    0.713
    ##    .MM_EpiDoc1       -0.965    0.158   -6.091      Inf    0.000
    ## 
    ## Variances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|)
    ##    .Dep_EpiDoc4       8.031    0.712   11.277      Inf    0.000
    ##    .MM_EpiDoc4        1.001    0.072   13.946      Inf    0.000
    ##    .Dep_EpiDoc2       8.255    0.782   10.561  208.549    0.000
    ##    .MM_EpiDoc2        0.297    0.041    7.177      Inf    0.000
    ##    .Dep_EpiDoc1      10.734    1.046   10.258      Inf    0.000
    ##    .MM_EpiDoc1        0.753    0.075   10.090      Inf    0.000
    ## 
    ## Defined Parameters:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|)
    ##     age_dep_indrct    0.019    0.006    3.038      Inf    0.002
    ##     age_depM_ndrct   -0.000    0.014   -0.018      Inf    0.986
    ##     age_direct        0.011    0.012    0.919  244.178    0.359
    ##     age_total         0.019    0.015    1.221      Inf    0.222
    ##     ed_dep_indirct    0.137    0.048    2.875  790.357    0.004
    ##     ed_depM_indrct    0.034    0.054    0.631  341.346    0.529
    ##     ed_direct         0.112    0.113    0.989      Inf    0.323
    ##     ed_total          0.171    0.075    2.282  603.479    0.023
    ##     dep0_direct       0.261    0.055    4.769      Inf    0.000
    ##     dep0_indirect     0.179    0.052    3.429  457.400    0.001
    ##     dep0_m_indirct    0.168    0.251    0.671  393.876    0.503
    ##     dep_total         0.609    0.252    2.417  401.613    0.016
    ##     m0_direct        -0.182    0.168   -1.078  947.239    0.281
    ##     m0_indirect       0.141    0.059    2.408      Inf    0.016
    ##     m0_total         -0.041    0.168   -0.243      Inf    0.808

``` r
parameterEstimates.mi(res_adj)
```

    ## 
    ## Rubin's (1987) rules were used to pool point and SE estimates across 20 imputed data sets, and to calculate degrees of freedom for each parameter's t test and CI.
    ## 
    ## 
    ##                  lhs op                                       rhs
    ## 1        Dep_EpiDoc4  ~                               Dep_EpiDoc1
    ## 2        Dep_EpiDoc4  ~                               Dep_EpiDoc2
    ## 3         MM_EpiDoc4  ~                                MM_EpiDoc2
    ## 4         MM_EpiDoc4  ~                                MM_EpiDoc1
    ## 5        Dep_EpiDoc2  ~                               Dep_EpiDoc1
    ## 6         MM_EpiDoc2  ~                                MM_EpiDoc1
    ## 7        Dep_EpiDoc4  ~                                MM_EpiDoc2
    ## 8        Dep_EpiDoc4  ~                                MM_EpiDoc1
    ## 9         MM_EpiDoc4  ~                               Dep_EpiDoc2
    ## 10        MM_EpiDoc4  ~                               Dep_EpiDoc1
    ## 11       Dep_EpiDoc2  ~                                MM_EpiDoc1
    ## 12        MM_EpiDoc2  ~                               Dep_EpiDoc1
    ## 13       Dep_EpiDoc1 ~~                                MM_EpiDoc1
    ## 14       Dep_EpiDoc2 ~~                                MM_EpiDoc2
    ## 15       Dep_EpiDoc4 ~~                                MM_EpiDoc4
    ## 16       Dep_EpiDoc4  ~                                       Age
    ## 17        MM_EpiDoc4  ~                                       Age
    ## 18       Dep_EpiDoc2  ~                                       Age
    ## 19        MM_EpiDoc2  ~                                       Age
    ## 20       Dep_EpiDoc1  ~                                       Age
    ## 21        MM_EpiDoc1  ~                                       Age
    ## 22       Dep_EpiDoc4  ~                                 Education
    ## 23        MM_EpiDoc4  ~                                 Education
    ## 24       Dep_EpiDoc2  ~                                 Education
    ## 25        MM_EpiDoc2  ~                                 Education
    ## 26       Dep_EpiDoc1  ~                                 Education
    ## 27        MM_EpiDoc1  ~                                 Education
    ## 28       Dep_EpiDoc4 ~~                               Dep_EpiDoc4
    ## 29        MM_EpiDoc4 ~~                                MM_EpiDoc4
    ## 30       Dep_EpiDoc2 ~~                               Dep_EpiDoc2
    ## 31        MM_EpiDoc2 ~~                                MM_EpiDoc2
    ## 32       Dep_EpiDoc1 ~~                               Dep_EpiDoc1
    ## 33        MM_EpiDoc1 ~~                                MM_EpiDoc1
    ## 34               Age ~~                                       Age
    ## 35               Age ~~                                 Education
    ## 36         Education ~~                                 Education
    ## 37       Dep_EpiDoc4 ~1                                          
    ## 38        MM_EpiDoc4 ~1                                          
    ## 39       Dep_EpiDoc2 ~1                                          
    ## 40        MM_EpiDoc2 ~1                                          
    ## 41       Dep_EpiDoc1 ~1                                          
    ## 42        MM_EpiDoc1 ~1                                          
    ## 43               Age ~1                                          
    ## 44         Education ~1                                          
    ## 45  age_dep_indirect :=              d0*dd01*dd13+d0*dd03+d1*dd13
    ## 46 age_depM_indirect := m0*md03+m1*md13+d0*dm01*md13+d0*dm01*md13
    ## 47        age_direct :=                                        d3
    ## 48         age_total :=        age_dep_indirect+age_depM_indirect
    ## 49   ed_dep_indirect :=                    de0*dd01*dd13+de1*dd13
    ## 50  ed_depM_indirect :=                             de1*dm01*md13
    ## 51         ed_direct :=                                       de3
    ## 52          ed_total :=          ed_dep_indirect+ed_depM_indirect
    ## 53       dep0_direct :=                                      dd03
    ## 54     dep0_indirect :=                       dd01*dd13+dd01*dd13
    ## 55   dep0_m_indirect :=                       dm01*md13+dm01*md13
    ## 56         dep_total := dep0_direct+dep0_indirect+dep0_m_indirect
    ## 57         m0_direct :=                                      md03
    ## 58       m0_indirect :=                       md01*dd13+mm01*md13
    ## 59          m0_total :=                     m0_direct+m0_indirect
    ##                label     est    se      t           df pvalue ci.lower ci.upper
    ## 1               dd03   0.261 0.055  4.769 1.704662e+03  0.000    0.154    0.369
    ## 2               dd13   0.225 0.056  3.981 7.429590e+02  0.000    0.114    0.335
    ## 3               mm13   0.295 0.095  3.108 1.166900e+07  0.002    0.109    0.481
    ## 4               mm03   0.475 0.062  7.667 8.258328e+06  0.000    0.354    0.597
    ## 5               dd01   0.398 0.055  7.185 3.491030e+02  0.000    0.289    0.507
    ## 6               mm01   0.228 0.054  4.207 4.368787e+31  0.000    0.122    0.335
    ## 7               md13   0.598 0.222  2.691 1.451251e+03  0.007    0.162    1.035
    ## 8               md03  -0.182 0.168 -1.078 9.472390e+02  0.281   -0.512    0.149
    ## 9               dm13   0.019 0.016  1.150 3.630380e+02  0.251   -0.013    0.050
    ## 10              dm03   0.018 0.016  1.130 8.557981e+03  0.258   -0.014    0.050
    ## 11              dm01   0.139 0.174  0.801 3.623940e+02  0.424   -0.203    0.481
    ## 12              md01   0.020 0.009  2.093 3.724402e+29  0.036    0.001    0.038
    ## 13                     0.776 0.193  4.016 5.012095e+29  0.000    0.398    1.155
    ## 14                     0.057 0.080  0.711 2.762390e+02  0.478   -0.100    0.214
    ## 15                     0.120 0.127  0.944 4.262910e+02  0.345   -0.129    0.369
    ## 16                d3   0.011 0.012  0.919 2.441780e+02  0.359   -0.013    0.035
    ## 17                m3   0.022 0.006  3.828 5.996417e+07  0.000    0.011    0.033
    ## 18                d1   0.014 0.013  1.095 1.215260e+02  0.276   -0.012    0.040
    ## 19                m1  -0.002 0.002 -0.743 1.257277e+29  0.458   -0.006    0.003
    ## 20                d0   0.045 0.015  2.974 1.015398e+32  0.003    0.015    0.074
    ## 21                m0   0.037 0.004  9.553 4.678267e+31  0.000    0.030    0.045
    ## 22               de3   0.112 0.113  0.989 5.384078e+03  0.323   -0.110    0.333
    ## 23               me3   0.003 0.041  0.081 2.104496e+05  0.935   -0.076    0.083
    ## 24               de1   0.407 0.123  3.316 3.067620e+02  0.001    0.165    0.648
    ## 25               me1  -0.025 0.020 -1.231 1.015776e+28  0.218   -0.065    0.015
    ## 26               de0   0.513 0.139  3.693 6.354181e+40  0.000    0.241    0.785
    ## 27               me0   0.011 0.036  0.306 2.853486e+40  0.760   -0.059    0.081
    ## 28                     8.031 0.712 11.277 1.593760e+03  0.000    6.634    9.427
    ## 29                     1.001 0.072 13.946 1.514246e+07  0.000    0.860    1.141
    ## 30                     8.255 0.782 10.561 2.085490e+02  0.000    6.714    9.796
    ## 31                     0.297 0.041  7.177 1.097986e+31  0.000    0.216    0.378
    ## 32                    10.734 1.046 10.258 5.299726e+29  0.000    8.683   12.785
    ## 33                     0.753 0.075 10.090 7.369901e+28  0.000    0.607    0.899
    ## 34                   175.024 0.000     NA           NA     NA  175.024  175.024
    ## 35                     5.613 0.000     NA           NA     NA    5.613    5.613
    ## 36                     1.237 0.000     NA           NA     NA    1.237    1.237
    ## 37                     0.709 0.468  1.516 2.882670e+02  0.131   -0.212    1.629
    ## 38                    -0.388 0.210 -1.842 3.428393e+07  0.065   -0.800    0.025
    ## 39                     0.018 0.516  0.035 1.751690e+02  0.972   -1.000    1.035
    ## 40                     0.096 0.071  1.348 8.544992e+33  0.178   -0.043    0.235
    ## 41                    -0.220 0.598 -0.368 1.293655e+45  0.713   -1.393    0.953
    ## 42                    -0.965 0.158 -6.091 6.496677e+44  0.000   -1.276   -0.655
    ## 43                    43.446 0.000     NA           NA     NA   43.446   43.446
    ## 44                     2.643 0.000     NA           NA     NA    2.643    2.643
    ## 45  age_dep_indirect   0.019 0.006  3.038 1.007801e+03  0.002    0.007    0.031
    ## 46 age_depM_indirect   0.000 0.014 -0.018 1.061837e+03  0.986   -0.028    0.027
    ## 47        age_direct   0.011 0.012  0.919 2.441780e+02  0.359   -0.013    0.035
    ## 48         age_total   0.019 0.015  1.221 2.932985e+03  0.222   -0.011    0.049
    ## 49   ed_dep_indirect   0.137 0.048  2.875 7.903570e+02  0.004    0.044    0.231
    ## 50  ed_depM_indirect   0.034 0.054  0.631 3.413460e+02  0.529   -0.072    0.139
    ## 51         ed_direct   0.112 0.113  0.989 5.384078e+03  0.323   -0.110    0.333
    ## 52          ed_total   0.171 0.075  2.282 6.034790e+02  0.023    0.024    0.318
    ## 53       dep0_direct   0.261 0.055  4.769 1.704662e+03  0.000    0.154    0.369
    ## 54     dep0_indirect   0.179 0.052  3.429 4.574000e+02  0.001    0.076    0.282
    ## 55   dep0_m_indirect   0.168 0.251  0.671 3.938760e+02  0.503   -0.325    0.661
    ## 56         dep_total   0.609 0.252  2.417 4.016130e+02  0.016    0.114    1.103
    ## 57         m0_direct  -0.182 0.168 -1.078 9.472390e+02  0.281   -0.512    0.149
    ## 58       m0_indirect   0.141 0.059  2.408 2.189804e+03  0.016    0.026    0.256
    ## 59          m0_total  -0.041 0.168 -0.243 1.507570e+03  0.808   -0.370    0.288

``` r
anova(res_adj)
```

    ## Robust correction can only be applied to pooled chi-squared statistic, not F statistic. "asymptotic" was switched to TRUE.

    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.

    ## Test statistic(s) pooled using the D4 pooling method.
    ##   Pooled statistic: "standard"
    ##   Method to robustify pooled statistic:  "yuan.bentler.mplus"
    ## 
    ##           Df Chisq Chisq diff Df diff Pr(>Chisq) RIV FMI
    ## Saturated  0     0                                      
    ## object     0     0                  0              0   0

``` r
source('./code/parTableToCSV.R')
res<-parTableToCSV(res_adj,path=('./Results/LaggedSEM_mi20_final'))
res$partable
```

    ##                label     est    se      t           df pvalue ci.lower ci.upper
    ## 1               dd03   0.261 0.055  4.769 1.704662e+03  0.000    0.154    0.369
    ## 2               dd13   0.225 0.056  3.981 7.429590e+02  0.000    0.114    0.335
    ## 3               mm13   0.295 0.095  3.108 1.166900e+07  0.002    0.109    0.481
    ## 4               mm03   0.475 0.062  7.667 8.258328e+06  0.000    0.354    0.597
    ## 5               dd01   0.398 0.055  7.185 3.491030e+02  0.000    0.289    0.507
    ## 6               mm01   0.228 0.054  4.207 4.368787e+31  0.000    0.122    0.335
    ## 7               md13   0.598 0.222  2.691 1.451251e+03  0.007    0.162    1.035
    ## 8               md03  -0.182 0.168 -1.078 9.472390e+02  0.281   -0.512    0.149
    ## 9               dm13   0.019 0.016  1.150 3.630380e+02  0.251   -0.013    0.050
    ## 10              dm03   0.018 0.016  1.130 8.557981e+03  0.258   -0.014    0.050
    ## 11              dm01   0.139 0.174  0.801 3.623940e+02  0.424   -0.203    0.481
    ## 12              md01   0.020 0.009  2.093 3.724402e+29  0.036    0.001    0.038
    ## 13                     0.776 0.193  4.016 5.012095e+29  0.000    0.398    1.155
    ## 14                     0.057 0.080  0.711 2.762390e+02  0.478   -0.100    0.214
    ## 15                     0.120 0.127  0.944 4.262910e+02  0.345   -0.129    0.369
    ## 16                d3   0.011 0.012  0.919 2.441780e+02  0.359   -0.013    0.035
    ## 17                m3   0.022 0.006  3.828 5.996417e+07  0.000    0.011    0.033
    ## 18                d1   0.014 0.013  1.095 1.215260e+02  0.276   -0.012    0.040
    ## 19                m1  -0.002 0.002 -0.743 1.257277e+29  0.458   -0.006    0.003
    ## 20                d0   0.045 0.015  2.974 1.015398e+32  0.003    0.015    0.074
    ## 21                m0   0.037 0.004  9.553 4.678267e+31  0.000    0.030    0.045
    ## 22               de3   0.112 0.113  0.989 5.384078e+03  0.323   -0.110    0.333
    ## 23               me3   0.003 0.041  0.081 2.104496e+05  0.935   -0.076    0.083
    ## 24               de1   0.407 0.123  3.316 3.067620e+02  0.001    0.165    0.648
    ## 25               me1  -0.025 0.020 -1.231 1.015776e+28  0.218   -0.065    0.015
    ## 26               de0   0.513 0.139  3.693 6.354181e+40  0.000    0.241    0.785
    ## 27               me0   0.011 0.036  0.306 2.853486e+40  0.760   -0.059    0.081
    ## 28                     8.031 0.712 11.277 1.593760e+03  0.000    6.634    9.427
    ## 29                     1.001 0.072 13.946 1.514246e+07  0.000    0.860    1.141
    ## 30                     8.255 0.782 10.561 2.085490e+02  0.000    6.714    9.796
    ## 31                     0.297 0.041  7.177 1.097986e+31  0.000    0.216    0.378
    ## 32                    10.734 1.046 10.258 5.299726e+29  0.000    8.683   12.785
    ## 33                     0.753 0.075 10.090 7.369901e+28  0.000    0.607    0.899
    ## 34                   175.024 0.000     NA           NA     NA  175.024  175.024
    ## 35                     5.613 0.000     NA           NA     NA    5.613    5.613
    ## 36                     1.237 0.000     NA           NA     NA    1.237    1.237
    ## 37                     0.709 0.468  1.516 2.882670e+02  0.131   -0.212    1.629
    ## 38                    -0.388 0.210 -1.842 3.428393e+07  0.065   -0.800    0.025
    ## 39                     0.018 0.516  0.035 1.751690e+02  0.972   -1.000    1.035
    ## 40                     0.096 0.071  1.348 8.544992e+33  0.178   -0.043    0.235
    ## 41                    -0.220 0.598 -0.368 1.293655e+45  0.713   -1.393    0.953
    ## 42                    -0.965 0.158 -6.091 6.496677e+44  0.000   -1.276   -0.655
    ## 43                    43.446 0.000     NA           NA     NA   43.446   43.446
    ## 44                     2.643 0.000     NA           NA     NA    2.643    2.643
    ## 45  age_dep_indirect   0.019 0.006  3.038 1.007801e+03  0.002    0.007    0.031
    ## 46 age_depM_indirect   0.000 0.014 -0.018 1.061837e+03  0.986   -0.028    0.027
    ## 47        age_direct   0.011 0.012  0.919 2.441780e+02  0.359   -0.013    0.035
    ## 48         age_total   0.019 0.015  1.221 2.932985e+03  0.222   -0.011    0.049
    ## 49   ed_dep_indirect   0.137 0.048  2.875 7.903570e+02  0.004    0.044    0.231
    ## 50  ed_depM_indirect   0.034 0.054  0.631 3.413460e+02  0.529   -0.072    0.139
    ## 51         ed_direct   0.112 0.113  0.989 5.384078e+03  0.323   -0.110    0.333
    ## 52          ed_total   0.171 0.075  2.282 6.034790e+02  0.023    0.024    0.318
    ## 53       dep0_direct   0.261 0.055  4.769 1.704662e+03  0.000    0.154    0.369
    ## 54     dep0_indirect   0.179 0.052  3.429 4.574000e+02  0.001    0.076    0.282
    ## 55   dep0_m_indirect   0.168 0.251  0.671 3.938760e+02  0.503   -0.325    0.661
    ## 56         dep_total   0.609 0.252  2.417 4.016130e+02  0.016    0.114    1.103
    ## 57         m0_direct  -0.182 0.168 -1.078 9.472390e+02  0.281   -0.512    0.149
    ## 58       m0_indirect   0.141 0.059  2.408 2.189804e+03  0.016    0.026    0.256
    ## 59          m0_total  -0.041 0.168 -0.243 1.507570e+03  0.808   -0.370    0.288
    ##     std.lv std.all std.nox
    ## 1    0.261   0.277   0.277
    ## 2    0.225   0.231   0.231
    ## 3    0.295   0.138   0.138
    ## 4    0.475   0.372   0.372
    ## 5    0.398   0.410   0.410
    ## 6    0.228   0.382   0.382
    ## 7    0.598   0.111   0.111
    ## 8   -0.182  -0.056  -0.056
    ## 9    0.019   0.048   0.048
    ## 10   0.018   0.049   0.049
    ## 11   0.139   0.042   0.042
    ## 12   0.020   0.112   0.112
    ## 13   0.776   0.273   0.273
    ## 14   0.057   0.036   0.036
    ## 15   0.120   0.042   0.042
    ## 16   0.011   0.046   0.003
    ## 17   0.022   0.224   0.017
    ## 18   0.014   0.057   0.004
    ## 19  -0.002  -0.037  -0.003
    ## 20   0.045   0.173   0.013
    ## 21   0.037   0.492   0.037
    ## 22   0.112   0.039   0.035
    ## 23   0.003   0.003   0.003
    ## 24   0.407   0.136   0.122
    ## 25  -0.025  -0.046  -0.042
    ## 26   0.513   0.167   0.150
    ## 27   0.011   0.012   0.011
    ## 28   8.031   0.771   0.771
    ## 29   1.001   0.614   0.614
    ## 30   8.255   0.749   0.749
    ## 31   0.297   0.832   0.832
    ## 32  10.734   0.920   0.920
    ## 33   0.753   0.753   0.753
    ## 34 175.024   1.000 175.024
    ## 35   5.613   0.382   5.613
    ## 36   1.237   1.000   1.237
    ## 37   0.709   0.220   0.220
    ## 38  -0.388  -0.304  -0.304
    ## 39   0.018   0.005   0.005
    ## 40   0.096   0.160   0.160
    ## 41  -0.220  -0.064  -0.064
    ## 42  -0.965  -0.965  -0.965
    ## 43  43.446   3.284  43.446
    ## 44   2.643   2.376   2.643
    ## 45   0.019   0.078   0.006
    ## 46   0.000  -0.030  -0.002
    ## 47   0.011   0.046   0.003
    ## 48   0.019   0.047   0.004
    ## 49   0.137   0.047   0.043
    ## 50   0.034   0.001   0.001
    ## 51   0.112   0.039   0.035
    ## 52   0.171   0.048   0.043
    ## 53   0.261   0.277   0.277
    ## 54   0.179   0.189   0.189
    ## 55   0.167   0.009   0.009
    ## 56   0.607   0.475   0.475
    ## 57  -0.182  -0.056  -0.056
    ## 58   0.141   0.068   0.068
    ## 59  -0.041   0.012   0.012
    ##                                                          label2     estimate_p
    ## 1                                       Dep_EpiDoc4~Dep_EpiDoc1      0.261 (0)
    ## 2                                       Dep_EpiDoc4~Dep_EpiDoc2      0.225 (0)
    ## 3                                         MM_EpiDoc4~MM_EpiDoc2  0.295 (0.002)
    ## 4                                         MM_EpiDoc4~MM_EpiDoc1      0.475 (0)
    ## 5                                       Dep_EpiDoc2~Dep_EpiDoc1      0.398 (0)
    ## 6                                         MM_EpiDoc2~MM_EpiDoc1      0.228 (0)
    ## 7                                        Dep_EpiDoc4~MM_EpiDoc2  0.598 (0.007)
    ## 8                                        Dep_EpiDoc4~MM_EpiDoc1 -0.182 (0.281)
    ## 9                                        MM_EpiDoc4~Dep_EpiDoc2  0.019 (0.251)
    ## 10                                       MM_EpiDoc4~Dep_EpiDoc1  0.018 (0.258)
    ## 11                                       Dep_EpiDoc2~MM_EpiDoc1  0.139 (0.424)
    ## 12                                       MM_EpiDoc2~Dep_EpiDoc1   0.02 (0.036)
    ## 13                                      Dep_EpiDoc1~~MM_EpiDoc1      0.776 (0)
    ## 14                                      Dep_EpiDoc2~~MM_EpiDoc2  0.057 (0.478)
    ## 15                                      Dep_EpiDoc4~~MM_EpiDoc4   0.12 (0.345)
    ## 16                                              Dep_EpiDoc4~Age  0.011 (0.359)
    ## 17                                               MM_EpiDoc4~Age      0.022 (0)
    ## 18                                              Dep_EpiDoc2~Age  0.014 (0.276)
    ## 19                                               MM_EpiDoc2~Age -0.002 (0.458)
    ## 20                                              Dep_EpiDoc1~Age  0.045 (0.003)
    ## 21                                               MM_EpiDoc1~Age      0.037 (0)
    ## 22                                        Dep_EpiDoc4~Education  0.112 (0.323)
    ## 23                                         MM_EpiDoc4~Education  0.003 (0.935)
    ## 24                                        Dep_EpiDoc2~Education  0.407 (0.001)
    ## 25                                         MM_EpiDoc2~Education -0.025 (0.218)
    ## 26                                        Dep_EpiDoc1~Education      0.513 (0)
    ## 27                                         MM_EpiDoc1~Education   0.011 (0.76)
    ## 28                                     Dep_EpiDoc4~~Dep_EpiDoc4      8.031 (0)
    ## 29                                       MM_EpiDoc4~~MM_EpiDoc4      1.001 (0)
    ## 30                                     Dep_EpiDoc2~~Dep_EpiDoc2      8.255 (0)
    ## 31                                       MM_EpiDoc2~~MM_EpiDoc2      0.297 (0)
    ## 32                                     Dep_EpiDoc1~~Dep_EpiDoc1     10.734 (0)
    ## 33                                       MM_EpiDoc1~~MM_EpiDoc1      0.753 (0)
    ## 34                                                     Age~~Age   175.024 (NA)
    ## 35                                               Age~~Education     5.613 (NA)
    ## 36                                         Education~~Education     1.237 (NA)
    ## 37                                                Dep_EpiDoc4~1  0.709 (0.131)
    ## 38                                                 MM_EpiDoc4~1 -0.388 (0.065)
    ## 39                                                Dep_EpiDoc2~1  0.018 (0.972)
    ## 40                                                 MM_EpiDoc2~1  0.096 (0.178)
    ## 41                                                Dep_EpiDoc1~1  -0.22 (0.713)
    ## 42                                                 MM_EpiDoc1~1     -0.965 (0)
    ## 43                                                        Age~1    43.446 (NA)
    ## 44                                                  Education~1     2.643 (NA)
    ## 45               age_dep_indirect:=d0*dd01*dd13+d0*dd03+d1*dd13  0.019 (0.002)
    ## 46 age_depM_indirect:=m0*md03+m1*md13+d0*dm01*md13+d0*dm01*md13      0 (0.986)
    ## 47                                               age_direct:=d3  0.011 (0.359)
    ## 48                age_total:=age_dep_indirect+age_depM_indirect  0.019 (0.222)
    ## 49                      ed_dep_indirect:=de0*dd01*dd13+de1*dd13  0.137 (0.004)
    ## 50                              ed_depM_indirect:=de1*dm01*md13  0.034 (0.529)
    ## 51                                               ed_direct:=de3  0.112 (0.323)
    ## 52                   ed_total:=ed_dep_indirect+ed_depM_indirect  0.171 (0.023)
    ## 53                                            dep0_direct:=dd03      0.261 (0)
    ## 54                           dep0_indirect:=dd01*dd13+dd01*dd13  0.179 (0.001)
    ## 55                         dep0_m_indirect:=dm01*md13+dm01*md13  0.168 (0.503)
    ## 56         dep_total:=dep0_direct+dep0_indirect+dep0_m_indirect  0.609 (0.016)
    ## 57                                              m0_direct:=md03 -0.182 (0.281)
    ## 58                             m0_indirect:=md01*dd13+mm01*md13  0.141 (0.016)
    ## 59                              m0_total:=m0_direct+m0_indirect -0.041 (0.808)
    ##                   CI
    ## 1      [0.154:0.369]
    ## 2      [0.114:0.335]
    ## 3      [0.109:0.481]
    ## 4      [0.354:0.597]
    ## 5      [0.289:0.507]
    ## 6      [0.122:0.335]
    ## 7      [0.162:1.035]
    ## 8     [-0.512:0.149]
    ## 9      [-0.013:0.05]
    ## 10     [-0.014:0.05]
    ## 11    [-0.203:0.481]
    ## 12     [0.001:0.038]
    ## 13     [0.398:1.155]
    ## 14      [-0.1:0.214]
    ## 15    [-0.129:0.369]
    ## 16    [-0.013:0.035]
    ## 17     [0.011:0.033]
    ## 18     [-0.012:0.04]
    ## 19    [-0.006:0.003]
    ## 20     [0.015:0.074]
    ## 21      [0.03:0.045]
    ## 22     [-0.11:0.333]
    ## 23    [-0.076:0.083]
    ## 24     [0.165:0.648]
    ## 25    [-0.065:0.015]
    ## 26     [0.241:0.785]
    ## 27    [-0.059:0.081]
    ## 28     [6.634:9.427]
    ## 29      [0.86:1.141]
    ## 30     [6.714:9.796]
    ## 31     [0.216:0.378]
    ## 32    [8.683:12.785]
    ## 33     [0.607:0.899]
    ## 34 [175.024:175.024]
    ## 35     [5.613:5.613]
    ## 36     [1.237:1.237]
    ## 37    [-0.212:1.629]
    ## 38      [-0.8:0.025]
    ## 39        [-1:1.035]
    ## 40    [-0.043:0.235]
    ## 41    [-1.393:0.953]
    ## 42   [-1.276:-0.655]
    ## 43   [43.446:43.446]
    ## 44     [2.643:2.643]
    ## 45     [0.007:0.031]
    ## 46    [-0.028:0.027]
    ## 47    [-0.013:0.035]
    ## 48    [-0.011:0.049]
    ## 49     [0.044:0.231]
    ## 50    [-0.072:0.139]
    ## 51     [-0.11:0.333]
    ## 52     [0.024:0.318]
    ## 53     [0.154:0.369]
    ## 54     [0.076:0.282]
    ## 55    [-0.325:0.661]
    ## 56     [0.114:1.103]
    ## 57    [-0.512:0.149]
    ## 58     [0.026:0.256]
    ## 59     [-0.37:0.288]

# Multi-group assessment

``` r
modelfree<-'
# evolução
Dep_EpiDoc4~Dep_EpiDoc1+Dep_EpiDoc2
MM_EpiDoc4~MM_EpiDoc2+MM_EpiDoc1
Dep_EpiDoc2~Dep_EpiDoc1
MM_EpiDoc2~MM_EpiDoc1

# causalidade - c/ temporalidade/lag
Dep_EpiDoc4~MM_EpiDoc2+MM_EpiDoc1
Dep_EpiDoc4~MM_EpiDoc2
MM_EpiDoc4~Dep_EpiDoc2+Dep_EpiDoc1
Dep_EpiDoc2~MM_EpiDoc1
MM_EpiDoc2~Dep_EpiDoc1

# corrrelação mesmo t
MM_EpiDoc1~~Dep_EpiDoc1
MM_EpiDoc2~~Dep_EpiDoc2
MM_EpiDoc4~~Dep_EpiDoc4

# causalidade - idade
Dep_EpiDoc4~0*Age
MM_EpiDoc4~Age
Dep_EpiDoc2~Age
MM_EpiDoc2~Age
Dep_EpiDoc1~Age
MM_EpiDoc1~Age

# educação
Dep_EpiDoc4~Education
Dep_EpiDoc2~Education
Dep_EpiDoc1~Education
MM_EpiDoc4~Education
MM_EpiDoc2~Education
MM_EpiDoc1~Education
'
```

``` r
resFre_g<-lavaan.mi::sem.mi(data=data.mi,
         model=modelfree,
         estimator = "ML", test = "yuan.bentler.mplus",
         likelihood = "wishart",
           sampling.weights='ipw',group='Sex',
         orthogonal=TRUE
         )
summary(resFre_g,fit.meas=TRUE,standardized=TRUE,rsquare=TRUE,ci=TRUE)
```

    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.

    ## lavaan.mi object fit to 20 imputed data sets using:
    ##  - lavaan    (0.6-19)
    ##  - lavaan.mi (0.1-0)
    ## See class?lavaan.mi help page for available methods. 
    ## 
    ## Convergence information:
    ## The model converged on 20 imputed data sets.
    ## Standard errors were available for all imputations.
    ## 
    ##   Estimator                                         ML
    ##   Optimization method                           NLMINB
    ##   Number of model parameters                        76
    ## 
    ##   Number of observations per group:                   
    ##     Male                                           355
    ##     Female                                         317
    ##   Sampling weights variable                        ipw
    ## 
    ## Model Test User Model:
    ## 
    ##                                                     Standard      Scaled
    ##   Test statistic                                       2.313       2.132
    ##   Degrees of freedom                                       2           2
    ##   P-value                                              0.315       0.344
    ##   Average scaling correction factor                                1.085
    ##   Pooling method                                          D4            
    ##     Pooled statistic                              "standard"            
    ##     "yuan.bentler.mplus" correction applied            AFTER     pooling
    ## 
    ## Model Test Baseline Model:
    ## 
    ##   Test statistic                              1048.982     846.230
    ##   Degrees of freedom                                54          54
    ##   P-value                                        0.000       0.000
    ##   Scaling correction factor                                  1.240
    ## 
    ## User Model versus Baseline Model:
    ## 
    ##   Comparative Fit Index (CFI)                    1.000       1.000
    ##   Tucker-Lewis Index (TLI)                       0.992       0.996
    ##                                                                   
    ##   Robust Comparative Fit Index (CFI)                         1.000
    ##   Robust Tucker-Lewis Index (TLI)                            0.996
    ## 
    ## Loglikelihood and Information Criteria:
    ## 
    ##   Loglikelihood user model (H0)              -7326.069   -7326.069
    ##   Scaling correction factor                                  1.736
    ##       for the MLR correction                                      
    ##   Loglikelihood unrestricted model (H1)      -7325.594   -7325.594
    ##   Scaling correction factor                                  1.717
    ##       for the MLR correction                                      
    ##                                                                   
    ##   Akaike (AIC)                               14804.139   14804.139
    ##   Bayesian (BIC)                             15146.919   15146.919
    ##   Sample-size adjusted Bayesian (SABIC)      14905.612   14905.612
    ## 
    ## Root Mean Square Error of Approximation:
    ## 
    ##   RMSEA                                          0.022       0.021
    ##   90 Percent confidence interval - lower         0.000       0.000
    ##   90 Percent confidence interval - upper         0.113       0.112
    ##   P-value H_0: RMSEA <= 0.050                    0.569       0.573
    ##   P-value H_0: RMSEA >= 0.080                    0.198       0.194
    ##                                                                   
    ##   Robust RMSEA                                               0.015
    ##   90 Percent confidence interval - lower                     0.000
    ##   90 Percent confidence interval - upper                     0.115
    ##   P-value H_0: Robust RMSEA <= 0.050                         0.580
    ##   P-value H_0: Robust RMSEA >= 0.080                         0.200
    ## 
    ## Standardized Root Mean Square Residual:
    ## 
    ##   SRMR                                           0.020       0.020
    ## 
    ## Parameter Estimates:
    ## 
    ##   Standard errors                                    Sandwich
    ##   Information bread                                  Observed
    ##   Observed information based on                       Hessian
    ##                                                              
    ##   Pooled across imputations              Rubin's (1987) rules
    ##   Augment within-imputation variance     Scale by average RIV
    ##   Wald test for pooled parameters          t(df) distribution
    ## 
    ##   Pooled t statistics with df >= 1000 are displayed with
    ##   df = Inf(inity) to save space. Although the t distribution
    ##   with large df closely approximates a standard normal
    ##   distribution, exact df for reporting these t tests can be
    ##   obtained from parameterEstimates.mi() 
    ## 
    ## 
    ## 
    ## Group 1 [Male]:
    ## 
    ## Regressions:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##   Dep_EpiDoc4 ~                                                         
    ##     Dep_EpiDoc1       0.164    0.079    2.079  527.391    0.038    0.009
    ##     Dep_EpiDoc2       0.234    0.085    2.761  327.279    0.006    0.067
    ##   MM_EpiDoc4 ~                                                          
    ##     MM_EpiDoc2        0.368    0.155    2.381      Inf    0.017    0.065
    ##     MM_EpiDoc1        0.564    0.089    6.314      Inf    0.000    0.389
    ##   Dep_EpiDoc2 ~                                                         
    ##     Dep_EpiDoc1       0.407    0.070    5.812  356.369    0.000    0.269
    ##   MM_EpiDoc2 ~                                                          
    ##     MM_EpiDoc1        0.102    0.047    2.195      Inf    0.028    0.011
    ##   Dep_EpiDoc4 ~                                                         
    ##     MM_EpiDoc2        0.530    0.292    1.817      Inf    0.069   -0.042
    ##     MM_EpiDoc1        0.046    0.208    0.223      Inf    0.824   -0.362
    ##   MM_EpiDoc4 ~                                                          
    ##     Dep_EpiDoc2       0.036    0.030    1.199  181.697    0.232   -0.023
    ##     Dep_EpiDoc1       0.005    0.029    0.162      Inf    0.871   -0.052
    ##   Dep_EpiDoc2 ~                                                         
    ##     MM_EpiDoc1        0.281    0.224    1.256  272.102    0.210   -0.160
    ##   MM_EpiDoc2 ~                                                          
    ##     Dep_EpiDoc1       0.021    0.015    1.363      Inf    0.173   -0.009
    ##   Dep_EpiDoc4 ~                                                         
    ##     Age               0.000                                        0.000
    ##   MM_EpiDoc4 ~                                                          
    ##     Age               0.018    0.007    2.456      Inf    0.014    0.004
    ##   Dep_EpiDoc2 ~                                                         
    ##     Age               0.001    0.016    0.033   77.311    0.974   -0.031
    ##   MM_EpiDoc2 ~                                                          
    ##     Age               0.002    0.002    0.671      Inf    0.502   -0.003
    ##   Dep_EpiDoc1 ~                                                         
    ##     Age               0.041    0.019    2.188      Inf    0.029    0.004
    ##   MM_EpiDoc1 ~                                                          
    ##     Age               0.031    0.005    6.080      Inf    0.000    0.021
    ##   Dep_EpiDoc4 ~                                                         
    ##     Education         0.162    0.166    0.977  599.381    0.329   -0.163
    ##   Dep_EpiDoc2 ~                                                         
    ##     Education         0.530    0.146    3.619  451.180    0.000    0.242
    ##   Dep_EpiDoc1 ~                                                         
    ##     Education         0.495    0.165    3.004      Inf    0.003    0.172
    ##   MM_EpiDoc4 ~                                                          
    ##     Education        -0.006    0.057   -0.097      Inf    0.923   -0.117
    ##   MM_EpiDoc2 ~                                                          
    ##     Education        -0.045    0.029   -1.524      Inf    0.128   -0.103
    ##   MM_EpiDoc1 ~                                                          
    ##     Education         0.031    0.054    0.574      Inf    0.566   -0.075
    ##  ci.upper   Std.lv  Std.all
    ##                            
    ##     0.319    0.164    0.174
    ##     0.400    0.234    0.241
    ##                            
    ##     0.671    0.368    0.131
    ##     0.739    0.564    0.424
    ##                            
    ##     0.544    0.407    0.419
    ##                            
    ##     0.194    0.102    0.216
    ##                            
    ##     1.102    0.530    0.083
    ##     0.455    0.046    0.015
    ##                            
    ##     0.095    0.036    0.083
    ##     0.061    0.005    0.011
    ##                            
    ##     0.722    0.281    0.091
    ##                            
    ##     0.050    0.021    0.139
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.032    0.018    0.203
    ##                            
    ##     0.032    0.001    0.003
    ##                            
    ##     0.006    0.002    0.049
    ##                            
    ##     0.078    0.041    0.191
    ##                            
    ##     0.041    0.031    0.465
    ##                            
    ##     0.487    0.162    0.060
    ##                            
    ##     0.818    0.530    0.192
    ##                            
    ##     0.819    0.495    0.175
    ##                            
    ##     0.106   -0.006   -0.005
    ##                            
    ##     0.013   -0.045   -0.107
    ##                            
    ##     0.138    0.031    0.035
    ## 
    ## Covariances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##  .Dep_EpiDoc1 ~~                                                        
    ##    .MM_EpiDoc1        0.548    0.216    2.540      Inf    0.011    0.125
    ##  .Dep_EpiDoc2 ~~                                                        
    ##    .MM_EpiDoc2        0.047    0.051    0.921  982.443    0.357   -0.053
    ##  .Dep_EpiDoc4 ~~                                                        
    ##    .MM_EpiDoc4        0.215    0.169    1.267  195.734    0.207   -0.119
    ##  ci.upper   Std.lv  Std.all
    ##                            
    ##     0.971    0.548    0.226
    ##                            
    ##     0.148    0.047    0.044
    ##                            
    ##     0.549    0.215    0.083
    ## 
    ## Intercepts:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##    .Dep_EpiDoc4       1.073    0.372    2.882  828.814    0.004    0.342
    ##    .MM_EpiDoc4       -0.336    0.270   -1.242      Inf    0.214   -0.866
    ##    .Dep_EpiDoc2      -0.013    0.633   -0.021  100.906    0.983   -1.269
    ##    .MM_EpiDoc2        0.041    0.055    0.745      Inf    0.456   -0.066
    ##    .Dep_EpiDoc1      -0.575    0.722   -0.796      Inf    0.426   -1.990
    ##    .MM_EpiDoc1       -0.818    0.190   -4.296      Inf    0.000   -1.191
    ##  ci.upper   Std.lv  Std.all
    ##     1.804    1.073    0.372
    ##     0.194   -0.336   -0.264
    ##     1.243   -0.013   -0.004
    ##     0.148    0.041    0.090
    ##     0.840   -0.575   -0.188
    ##    -0.445   -0.818   -0.856
    ## 
    ## Variances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##    .Dep_EpiDoc4       6.879    0.992    6.937      Inf    0.000    4.934
    ##    .MM_EpiDoc4        0.971    0.106    9.148      Inf    0.000    0.763
    ##    .Dep_EpiDoc2       6.187    0.865    7.154  101.737    0.000    4.472
    ##    .MM_EpiDoc2        0.187    0.051    3.650      Inf    0.000    0.086
    ##    .Dep_EpiDoc1       8.427    1.162    7.250      Inf    0.000    6.149
    ##    .MM_EpiDoc1        0.700    0.086    8.105      Inf    0.000    0.531
    ##  ci.upper   Std.lv  Std.all
    ##     8.824    6.879    0.828
    ##     1.179    0.971    0.602
    ##     7.903    6.187    0.704
    ##     0.287    0.187    0.911
    ##    10.705    8.427    0.903
    ##     0.870    0.700    0.768
    ## 
    ## R-Square:
    ##                    Estimate
    ##     Dep_EpiDoc4       0.172
    ##     MM_EpiDoc4        0.398
    ##     Dep_EpiDoc2       0.296
    ##     MM_EpiDoc2        0.089
    ##     Dep_EpiDoc1       0.097
    ##     MM_EpiDoc1        0.232
    ## 
    ## 
    ## Group 2 [Female]:
    ## 
    ## Regressions:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##   Dep_EpiDoc4 ~                                                         
    ##     Dep_EpiDoc1       0.342    0.075    4.555      Inf    0.000    0.195
    ##     Dep_EpiDoc2       0.220    0.075    2.941  943.731    0.003    0.073
    ##   MM_EpiDoc4 ~                                                          
    ##     MM_EpiDoc2        0.282    0.114    2.481      Inf    0.013    0.059
    ##     MM_EpiDoc1        0.350    0.088    4.000      Inf    0.000    0.179
    ##   Dep_EpiDoc2 ~                                                         
    ##     Dep_EpiDoc1       0.374    0.084    4.477  814.636    0.000    0.210
    ##   MM_EpiDoc2 ~                                                          
    ##     MM_EpiDoc1        0.366    0.097    3.776      Inf    0.000    0.176
    ##   Dep_EpiDoc4 ~                                                         
    ##     MM_EpiDoc2        0.669    0.307    2.180      Inf    0.029    0.067
    ##     MM_EpiDoc1       -0.320    0.218   -1.471      Inf    0.141   -0.747
    ##   MM_EpiDoc4 ~                                                          
    ##     Dep_EpiDoc2      -0.001    0.017   -0.042  495.077    0.967   -0.034
    ##     Dep_EpiDoc1       0.025    0.018    1.401      Inf    0.161   -0.010
    ##   Dep_EpiDoc2 ~                                                         
    ##     MM_EpiDoc1       -0.088    0.273   -0.322  585.338    0.748   -0.623
    ##   MM_EpiDoc2 ~                                                          
    ##     Dep_EpiDoc1       0.007    0.013    0.561      Inf    0.575   -0.018
    ##   Dep_EpiDoc4 ~                                                         
    ##     Age               0.000                                        0.000
    ##   MM_EpiDoc4 ~                                                          
    ##     Age               0.030    0.007    4.322      Inf    0.000    0.016
    ##   Dep_EpiDoc2 ~                                                         
    ##     Age               0.040    0.022    1.866  689.305    0.062   -0.002
    ##   MM_EpiDoc2 ~                                                          
    ##     Age              -0.005    0.005   -1.185      Inf    0.236   -0.014
    ##   Dep_EpiDoc1 ~                                                         
    ##     Age               0.055    0.020    2.761      Inf    0.006    0.016
    ##   MM_EpiDoc1 ~                                                          
    ##     Age               0.049    0.005    9.572      Inf    0.000    0.039
    ##   Dep_EpiDoc4 ~                                                         
    ##     Education         0.163    0.166    0.978      Inf    0.328   -0.163
    ##   Dep_EpiDoc2 ~                                                         
    ##     Education         0.373    0.210    1.774  309.152    0.077   -0.041
    ##   Dep_EpiDoc1 ~                                                         
    ##     Education         0.777    0.213    3.644      Inf    0.000    0.359
    ##   MM_EpiDoc4 ~                                                          
    ##     Education         0.045    0.060    0.746      Inf    0.456   -0.073
    ##   MM_EpiDoc2 ~                                                          
    ##     Education         0.028    0.033    0.853      Inf    0.394   -0.036
    ##   MM_EpiDoc1 ~                                                          
    ##     Education         0.025    0.050    0.504      Inf    0.614   -0.073
    ##  ci.upper   Std.lv  Std.all
    ##                            
    ##     0.489    0.342    0.357
    ##     0.366    0.220    0.225
    ##                            
    ##     0.504    0.282    0.162
    ##     0.522    0.350    0.290
    ##                            
    ##     0.538    0.374    0.380
    ##                            
    ##     0.556    0.366    0.526
    ##                            
    ##     1.270    0.669    0.136
    ##     0.107   -0.320   -0.094
    ##                            
    ##     0.033   -0.001   -0.002
    ##     0.060    0.025    0.074
    ##                            
    ##     0.448   -0.088   -0.025
    ##                            
    ##     0.033    0.007    0.038
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.043    0.030    0.270
    ##                            
    ##     0.083    0.040    0.127
    ##                            
    ##     0.004   -0.005   -0.086
    ##                            
    ##     0.095    0.055    0.172
    ##                            
    ##     0.059    0.049    0.537
    ##                            
    ##     0.489    0.163    0.050
    ##                            
    ##     0.786    0.373    0.112
    ##                            
    ##     1.195    0.777    0.230
    ##                            
    ##     0.163    0.045    0.039
    ##                            
    ##     0.092    0.028    0.042
    ##                            
    ##     0.124    0.025    0.027
    ## 
    ## Covariances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##  .Dep_EpiDoc1 ~~                                                        
    ##    .MM_EpiDoc1        0.921    0.309    2.976      Inf    0.003    0.314
    ##  .Dep_EpiDoc2 ~~                                                        
    ##    .MM_EpiDoc2        0.081    0.154    0.525  273.474    0.600   -0.223
    ##  .Dep_EpiDoc4 ~~                                                        
    ##    .MM_EpiDoc4       -0.026    0.189   -0.137      Inf    0.891   -0.397
    ##  ci.upper   Std.lv  Std.all
    ##                            
    ##     1.527    0.921    0.293
    ##                            
    ##     0.385    0.081    0.039
    ##                            
    ##     0.345   -0.026   -0.009
    ## 
    ## Intercepts:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##    .Dep_EpiDoc4       0.968    0.377    2.566      Inf    0.010    0.228
    ##    .MM_EpiDoc4       -0.611    0.243   -2.514      Inf    0.012   -1.087
    ##    .Dep_EpiDoc2      -0.517    0.825   -0.627  381.594    0.531   -2.139
    ##    .MM_EpiDoc2        0.133    0.159    0.834      Inf    0.404   -0.179
    ##    .Dep_EpiDoc1      -0.568    0.799   -0.710      Inf    0.477   -2.134
    ##    .MM_EpiDoc1       -1.396    0.226   -6.169      Inf    0.000   -1.840
    ##  ci.upper   Std.lv  Std.all
    ##     1.708    0.968    0.269
    ##    -0.135   -0.611   -0.480
    ##     1.105   -0.517   -0.140
    ##     0.445    0.133    0.181
    ##     0.999   -0.568   -0.151
    ##    -0.953   -1.396   -1.325
    ## 
    ## Variances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##    .Dep_EpiDoc4       9.339    1.049    8.903      Inf    0.000    7.281
    ##    .MM_EpiDoc4        0.997    0.084   11.825      Inf    0.000    0.832
    ##    .Dep_EpiDoc2      10.657    1.325    8.044  331.704    0.000    8.051
    ##    .MM_EpiDoc2        0.400    0.061    6.522      Inf    0.000    0.280
    ##    .Dep_EpiDoc1      12.629    1.736    7.276      Inf    0.000    9.227
    ##    .MM_EpiDoc1        0.781    0.121    6.445      Inf    0.000    0.543
    ##  ci.upper   Std.lv  Std.all
    ##    11.396    9.339    0.722
    ##     1.162    0.997    0.615
    ##    13.263   10.657    0.782
    ##     0.520    0.400    0.742
    ##    16.031   12.629    0.895
    ##     1.018    0.781    0.703
    ## 
    ## R-Square:
    ##                    Estimate
    ##     Dep_EpiDoc4       0.278
    ##     MM_EpiDoc4        0.385
    ##     Dep_EpiDoc2       0.218
    ##     MM_EpiDoc2        0.258
    ##     Dep_EpiDoc1       0.105
    ##     MM_EpiDoc1        0.297

## theoretical restrictions on age, auto-regressive or, cross-lagged components

### restriction on residual cov

``` r
model_covR<-'
# evolução
Dep_EpiDoc4~Dep_EpiDoc1+Dep_EpiDoc2
MM_EpiDoc4~MM_EpiDoc2+MM_EpiDoc1
Dep_EpiDoc2~Dep_EpiDoc1
MM_EpiDoc2~MM_EpiDoc1 

# causalidade - c/ temporalidade/lag
Dep_EpiDoc4~MM_EpiDoc2+MM_EpiDoc1
MM_EpiDoc4~Dep_EpiDoc2+Dep_EpiDoc1
Dep_EpiDoc2~MM_EpiDoc1
MM_EpiDoc2~Dep_EpiDoc1

# corrrelação mesmo t
MM_EpiDoc1~~Dep_EpiDoc1
MM_EpiDoc2~~0*Dep_EpiDoc2
MM_EpiDoc4~~0*Dep_EpiDoc4

# causalidade - idade
Dep_EpiDoc4~0*Age
Dep_EpiDoc2~Age
Dep_EpiDoc1~Age

MM_EpiDoc4~Age
MM_EpiDoc2~Age
MM_EpiDoc1~Age

# educação
Dep_EpiDoc4~Education
Dep_EpiDoc2~Education
Dep_EpiDoc1~Education
MM_EpiDoc4~Education
MM_EpiDoc2~Education
MM_EpiDoc1~Education
'


res_covR_g<-lavaan.mi::sem.mi(data=data.mi,
         model=model_covR,
         estimator = "ML", test = "yuan.bentler.mplus",
         likelihood = "wishart",
           sampling.weights='ipw',group='Sex',
         orthogonal=TRUE
         )
summary(res_covR_g,fit.meas=TRUE,standardized=TRUE,rsquare=TRUE,ci=TRUE)
```

    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.

    ## lavaan.mi object fit to 20 imputed data sets using:
    ##  - lavaan    (0.6-19)
    ##  - lavaan.mi (0.1-0)
    ## See class?lavaan.mi help page for available methods. 
    ## 
    ## Convergence information:
    ## The model converged on 20 imputed data sets.
    ## Standard errors were available for all imputations.
    ## 
    ##   Estimator                                         ML
    ##   Optimization method                           NLMINB
    ##   Number of model parameters                        72
    ## 
    ##   Number of observations per group:                   
    ##     Male                                           355
    ##     Female                                         317
    ##   Sampling weights variable                        ipw
    ## 
    ## Model Test User Model:
    ## 
    ##                                                     Standard      Scaled
    ##   Test statistic                                       5.172       4.610
    ##   Degrees of freedom                                       6           6
    ##   P-value                                              0.522       0.595
    ##   Average scaling correction factor                                1.122
    ##   Pooling method                                          D4            
    ##     Pooled statistic                              "standard"            
    ##     "yuan.bentler.mplus" correction applied            AFTER     pooling
    ## 
    ## Model Test Baseline Model:
    ## 
    ##   Test statistic                              1048.982     846.230
    ##   Degrees of freedom                                54          54
    ##   P-value                                        0.000       0.000
    ##   Scaling correction factor                                  1.240
    ## 
    ## User Model versus Baseline Model:
    ## 
    ##   Comparative Fit Index (CFI)                    1.000       1.000
    ##   Tucker-Lewis Index (TLI)                       1.007       1.016
    ##                                                                   
    ##   Robust Comparative Fit Index (CFI)                         1.000
    ##   Robust Tucker-Lewis Index (TLI)                            1.014
    ## 
    ## Loglikelihood and Information Criteria:
    ## 
    ##   Loglikelihood user model (H0)              -7327.716   -7327.716
    ##   Scaling correction factor                                  1.765
    ##       for the MLR correction                                      
    ##   Loglikelihood unrestricted model (H1)      -7325.594   -7325.594
    ##   Scaling correction factor                                  1.717
    ##       for the MLR correction                                      
    ##                                                                   
    ##   Akaike (AIC)                               14799.432   14799.432
    ##   Bayesian (BIC)                             15124.171   15124.171
    ##   Sample-size adjusted Bayesian (SABIC)      14895.565   14895.565
    ## 
    ## Root Mean Square Error of Approximation:
    ## 
    ##   RMSEA                                          0.000       0.000
    ##   90 Percent confidence interval - lower         0.000       0.000
    ##   90 Percent confidence interval - upper         0.065       0.057
    ##   P-value H_0: RMSEA <= 0.050                    0.864       0.919
    ##   P-value H_0: RMSEA >= 0.080                    0.014       0.005
    ##                                                                   
    ##   Robust RMSEA                                               0.000
    ##   90 Percent confidence interval - lower                     0.000
    ##   90 Percent confidence interval - upper                     0.065
    ##   P-value H_0: Robust RMSEA <= 0.050                         0.879
    ##   P-value H_0: Robust RMSEA >= 0.080                         0.015
    ## 
    ## Standardized Root Mean Square Residual:
    ## 
    ##   SRMR                                           0.021       0.021
    ## 
    ## Parameter Estimates:
    ## 
    ##   Standard errors                                    Sandwich
    ##   Information bread                                  Observed
    ##   Observed information based on                       Hessian
    ##                                                              
    ##   Pooled across imputations              Rubin's (1987) rules
    ##   Augment within-imputation variance     Scale by average RIV
    ##   Wald test for pooled parameters          t(df) distribution
    ## 
    ##   Pooled t statistics with df >= 1000 are displayed with
    ##   df = Inf(inity) to save space. Although the t distribution
    ##   with large df closely approximates a standard normal
    ##   distribution, exact df for reporting these t tests can be
    ##   obtained from parameterEstimates.mi() 
    ## 
    ## 
    ## 
    ## Group 1 [Male]:
    ## 
    ## Regressions:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##   Dep_EpiDoc4 ~                                                         
    ##     Dep_EpiDoc1       0.164    0.078    2.091  527.387    0.037    0.010
    ##     Dep_EpiDoc2       0.234    0.084    2.777  327.279    0.006    0.068
    ##   MM_EpiDoc4 ~                                                          
    ##     MM_EpiDoc2        0.367    0.153    2.394      Inf    0.017    0.067
    ##     MM_EpiDoc1        0.560    0.089    6.278      Inf    0.000    0.385
    ##   Dep_EpiDoc2 ~                                                         
    ##     Dep_EpiDoc1       0.407    0.070    5.845  356.369    0.000    0.270
    ##   MM_EpiDoc2 ~                                                          
    ##     MM_EpiDoc1        0.102    0.046    2.207      Inf    0.027    0.011
    ##   Dep_EpiDoc4 ~                                                         
    ##     MM_EpiDoc2        0.530    0.290    1.827      Inf    0.068   -0.039
    ##     MM_EpiDoc1        0.046    0.207    0.224      Inf    0.823   -0.360
    ##   MM_EpiDoc4 ~                                                          
    ##     Dep_EpiDoc2       0.036    0.030    1.206  179.079    0.230   -0.023
    ##     Dep_EpiDoc1       0.005    0.029    0.157      Inf    0.875   -0.052
    ##   Dep_EpiDoc2 ~                                                         
    ##     MM_EpiDoc1        0.281    0.223    1.263  272.102    0.208   -0.157
    ##   MM_EpiDoc2 ~                                                          
    ##     Dep_EpiDoc1       0.021    0.015    1.371      Inf    0.170   -0.009
    ##   Dep_EpiDoc4 ~                                                         
    ##     Age               0.000                                        0.000
    ##   Dep_EpiDoc2 ~                                                         
    ##     Age               0.001    0.016    0.033   77.311    0.973   -0.031
    ##   Dep_EpiDoc1 ~                                                         
    ##     Age               0.041    0.019    2.200      Inf    0.028    0.004
    ##   MM_EpiDoc4 ~                                                          
    ##     Age               0.019    0.007    2.537      Inf    0.011    0.004
    ##   MM_EpiDoc2 ~                                                          
    ##     Age               0.002    0.002    0.675      Inf    0.500   -0.003
    ##   MM_EpiDoc1 ~                                                          
    ##     Age               0.031    0.005    6.115      Inf    0.000    0.021
    ##   Dep_EpiDoc4 ~                                                         
    ##     Education         0.162    0.165    0.983  599.378    0.326   -0.161
    ##   Dep_EpiDoc2 ~                                                         
    ##     Education         0.530    0.146    3.640  451.180    0.000    0.244
    ##   Dep_EpiDoc1 ~                                                         
    ##     Education         0.495    0.164    3.021      Inf    0.003    0.174
    ##   MM_EpiDoc4 ~                                                          
    ##     Education        -0.008    0.056   -0.146      Inf    0.884   -0.118
    ##   MM_EpiDoc2 ~                                                          
    ##     Education        -0.045    0.029   -1.532      Inf    0.125   -0.102
    ##   MM_EpiDoc1 ~                                                          
    ##     Education         0.031    0.054    0.578      Inf    0.564   -0.075
    ##  ci.upper   Std.lv  Std.all
    ##                            
    ##     0.318    0.164    0.174
    ##     0.399    0.234    0.241
    ##                            
    ##     0.668    0.367    0.131
    ##     0.735    0.560    0.421
    ##                            
    ##     0.543    0.407    0.419
    ##                            
    ##     0.193    0.102    0.216
    ##                            
    ##     1.099    0.530    0.083
    ##     0.453    0.046    0.015
    ##                            
    ##     0.094    0.036    0.083
    ##     0.061    0.005    0.011
    ##                            
    ##     0.720    0.281    0.091
    ##                            
    ##     0.050    0.021    0.139
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.032    0.001    0.003
    ##                            
    ##     0.077    0.041    0.191
    ##                            
    ##     0.033    0.019    0.209
    ##                            
    ##     0.006    0.002    0.049
    ##                            
    ##     0.041    0.031    0.465
    ##                            
    ##     0.485    0.162    0.061
    ##                            
    ##     0.816    0.530    0.192
    ##                            
    ##     0.817    0.495    0.175
    ##                            
    ##     0.102   -0.008   -0.007
    ##                            
    ##     0.013   -0.045   -0.107
    ##                            
    ##     0.137    0.031    0.035
    ## 
    ## Covariances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##  .Dep_EpiDoc1 ~~                                                        
    ##    .MM_EpiDoc1        0.548    0.214    2.555      Inf    0.011    0.128
    ##  .Dep_EpiDoc2 ~~                                                        
    ##    .MM_EpiDoc2        0.000                                        0.000
    ##  .Dep_EpiDoc4 ~~                                                        
    ##    .MM_EpiDoc4        0.000                                        0.000
    ##  ci.upper   Std.lv  Std.all
    ##                            
    ##     0.968    0.548    0.226
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.000    0.000    0.000
    ## 
    ## Intercepts:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##    .Dep_EpiDoc4       1.073    0.370    2.898  828.811    0.004    0.346
    ##    .MM_EpiDoc4       -0.352    0.271   -1.299      Inf    0.194   -0.883
    ##    .Dep_EpiDoc2      -0.013    0.630   -0.021  100.906    0.983   -1.262
    ##    .MM_EpiDoc2        0.041    0.054    0.749      Inf    0.454   -0.066
    ##    .Dep_EpiDoc1      -0.575    0.718   -0.801      Inf    0.423   -1.982
    ##    .MM_EpiDoc1       -0.818    0.189   -4.320      Inf    0.000   -1.189
    ##  ci.upper   Std.lv  Std.all
    ##     1.800    1.073    0.373
    ##     0.179   -0.352   -0.277
    ##     1.236   -0.013   -0.004
    ##     0.147    0.041    0.090
    ##     0.832   -0.575   -0.188
    ##    -0.447   -0.818   -0.856
    ## 
    ## Variances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##    .Dep_EpiDoc4       6.879    0.986    6.976      Inf    0.000    4.945
    ##    .MM_EpiDoc4        0.971    0.106    9.200      Inf    0.000    0.764
    ##    .Dep_EpiDoc2       6.187    0.860    7.194  101.737    0.000    4.481
    ##    .MM_EpiDoc2        0.187    0.051    3.671      Inf    0.000    0.087
    ##    .Dep_EpiDoc1       8.427    1.156    7.291      Inf    0.000    6.162
    ##    .MM_EpiDoc1        0.700    0.086    8.151      Inf    0.000    0.532
    ##  ci.upper   Std.lv  Std.all
    ##     8.813    6.879    0.830
    ##     1.178    0.971    0.601
    ##     7.894    6.187    0.704
    ##     0.287    0.187    0.911
    ##    10.692    8.427    0.903
    ##     0.869    0.700    0.768
    ## 
    ## R-Square:
    ##                    Estimate
    ##     Dep_EpiDoc4       0.170
    ##     MM_EpiDoc4        0.399
    ##     Dep_EpiDoc2       0.296
    ##     MM_EpiDoc2        0.089
    ##     Dep_EpiDoc1       0.097
    ##     MM_EpiDoc1        0.232
    ## 
    ## 
    ## Group 2 [Female]:
    ## 
    ## Regressions:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##   Dep_EpiDoc4 ~                                                         
    ##     Dep_EpiDoc1       0.342    0.075    4.581      Inf    0.000    0.195
    ##     Dep_EpiDoc2       0.220    0.074    2.958  943.724    0.003    0.074
    ##   MM_EpiDoc4 ~                                                          
    ##     MM_EpiDoc2        0.282    0.113    2.494      Inf    0.013    0.060
    ##     MM_EpiDoc1        0.350    0.087    4.024      Inf    0.000    0.180
    ##   Dep_EpiDoc2 ~                                                         
    ##     Dep_EpiDoc1       0.374    0.083    4.502  814.637    0.000    0.211
    ##   MM_EpiDoc2 ~                                                          
    ##     MM_EpiDoc1        0.366    0.096    3.798      Inf    0.000    0.177
    ##   Dep_EpiDoc4 ~                                                         
    ##     MM_EpiDoc2        0.669    0.305    2.192      Inf    0.028    0.070
    ##     MM_EpiDoc1       -0.320    0.216   -1.479      Inf    0.139   -0.745
    ##   MM_EpiDoc4 ~                                                          
    ##     Dep_EpiDoc2      -0.001    0.017   -0.042  494.987    0.967   -0.034
    ##     Dep_EpiDoc1       0.025    0.018    1.409      Inf    0.159   -0.010
    ##   Dep_EpiDoc2 ~                                                         
    ##     MM_EpiDoc1       -0.088    0.271   -0.323  585.338    0.746   -0.620
    ##   MM_EpiDoc2 ~                                                          
    ##     Dep_EpiDoc1       0.007    0.013    0.564      Inf    0.573   -0.018
    ##   Dep_EpiDoc4 ~                                                         
    ##     Age               0.000                                        0.000
    ##   Dep_EpiDoc2 ~                                                         
    ##     Age               0.040    0.021    1.877  689.305    0.061   -0.002
    ##   Dep_EpiDoc1 ~                                                         
    ##     Age               0.055    0.020    2.777      Inf    0.005    0.016
    ##   MM_EpiDoc4 ~                                                          
    ##     Age               0.030    0.007    4.346      Inf    0.000    0.016
    ##   MM_EpiDoc2 ~                                                          
    ##     Age              -0.005    0.005   -1.192      Inf    0.233   -0.014
    ##   MM_EpiDoc1 ~                                                          
    ##     Age               0.049    0.005    9.626      Inf    0.000    0.039
    ##   Dep_EpiDoc4 ~                                                         
    ##     Education         0.163    0.165    0.984      Inf    0.325   -0.162
    ##   Dep_EpiDoc2 ~                                                         
    ##     Education         0.373    0.209    1.784  309.152    0.075   -0.038
    ##   Dep_EpiDoc1 ~                                                         
    ##     Education         0.777    0.212    3.665      Inf    0.000    0.361
    ##   MM_EpiDoc4 ~                                                          
    ##     Education         0.045    0.060    0.750      Inf    0.453   -0.072
    ##   MM_EpiDoc2 ~                                                          
    ##     Education         0.028    0.033    0.858      Inf    0.391   -0.036
    ##   MM_EpiDoc1 ~                                                          
    ##     Education         0.025    0.050    0.507      Inf    0.612   -0.073
    ##  ci.upper   Std.lv  Std.all
    ##                            
    ##     0.488    0.342    0.357
    ##     0.365    0.220    0.226
    ##                            
    ##     0.503    0.282    0.162
    ##     0.521    0.350    0.290
    ##                            
    ##     0.537    0.374    0.380
    ##                            
    ##     0.555    0.366    0.526
    ##                            
    ##     1.267    0.669    0.137
    ##     0.104   -0.320   -0.094
    ##                            
    ##     0.033   -0.001   -0.002
    ##     0.060    0.025    0.074
    ##                            
    ##     0.445   -0.088   -0.025
    ##                            
    ##     0.033    0.007    0.038
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.082    0.040    0.127
    ##                            
    ##     0.094    0.055    0.172
    ##                            
    ##     0.043    0.030    0.270
    ##                            
    ##     0.003   -0.005   -0.086
    ##                            
    ##     0.059    0.049    0.537
    ##                            
    ##     0.487    0.163    0.050
    ##                            
    ##     0.784    0.373    0.112
    ##                            
    ##     1.193    0.777    0.230
    ##                            
    ##     0.162    0.045    0.039
    ##                            
    ##     0.092    0.028    0.042
    ##                            
    ##     0.123    0.025    0.027
    ## 
    ## Covariances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##  .Dep_EpiDoc1 ~~                                                        
    ##    .MM_EpiDoc1        0.921    0.308    2.993      Inf    0.003    0.318
    ##  .Dep_EpiDoc2 ~~                                                        
    ##    .MM_EpiDoc2        0.000                                        0.000
    ##  .Dep_EpiDoc4 ~~                                                        
    ##    .MM_EpiDoc4        0.000                                        0.000
    ##  ci.upper   Std.lv  Std.all
    ##                            
    ##     1.524    0.921    0.293
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.000    0.000    0.000
    ## 
    ## Intercepts:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##    .Dep_EpiDoc4       0.968    0.375    2.580      Inf    0.010    0.233
    ##    .MM_EpiDoc4       -0.611    0.242   -2.527      Inf    0.012   -1.084
    ##    .Dep_EpiDoc2      -0.517    0.820   -0.630  381.594    0.529   -2.130
    ##    .MM_EpiDoc2        0.133    0.158    0.839      Inf    0.401   -0.177
    ##    .Dep_EpiDoc1      -0.568    0.795   -0.714      Inf    0.475   -2.125
    ##    .MM_EpiDoc1       -1.396    0.225   -6.204      Inf    0.000   -1.838
    ##  ci.upper   Std.lv  Std.all
    ##     1.704    0.968    0.270
    ##    -0.137   -0.611   -0.479
    ##     1.096   -0.517   -0.140
    ##     0.443    0.133    0.181
    ##     0.990   -0.568   -0.151
    ##    -0.955   -1.396   -1.325
    ## 
    ## Variances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##    .Dep_EpiDoc4       9.339    1.043    8.954      Inf    0.000    7.293
    ##    .MM_EpiDoc4        0.997    0.084   11.892      Inf    0.000    0.833
    ##    .Dep_EpiDoc2      10.657    1.317    8.090  331.703    0.000    8.066
    ##    .MM_EpiDoc2        0.400    0.061    6.558      Inf    0.000    0.281
    ##    .Dep_EpiDoc1      12.629    1.726    7.317      Inf    0.000    9.247
    ##    .MM_EpiDoc1        0.781    0.120    6.482      Inf    0.000    0.545
    ##  ci.upper   Std.lv  Std.all
    ##    11.385    9.339    0.723
    ##     1.161    0.997    0.615
    ##    13.249   10.657    0.782
    ##     0.520    0.400    0.742
    ##    16.012   12.629    0.895
    ##     1.017    0.781    0.703
    ## 
    ## R-Square:
    ##                    Estimate
    ##     Dep_EpiDoc4       0.277
    ##     MM_EpiDoc4        0.385
    ##     Dep_EpiDoc2       0.218
    ##     MM_EpiDoc2        0.258
    ##     Dep_EpiDoc1       0.105
    ##     MM_EpiDoc1        0.297

``` r
anova(resFre_g,  res_covR_g)
```

    ## Robust correction can only be applied to pooled chi-squared statistic, not F statistic. "asymptotic" was switched to TRUE.

    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.
    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.
    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.
    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.

    ## Test statistic(s) pooled using the D4 pooling method.
    ##   Pooled statistic: "standard"
    ##   Method to robustify pooled statistic:  "yuan.bentler.mplus"
    ## 
    ##            Df  Chisq Chisq diff Df diff Pr(>Chisq)    RIV     FMI
    ## resFre_g    2 2.3128                                             
    ## res_covR_g  6 5.1723     2.4926       4    0.64597 0.3435 0.25568

### restriction on age direct effect

- we free equality on cross-lagged effect… but if gender (and contextual
  issues matter, time point 4 would differ also… )
- equal among group (doubt/hypotheses under test here)

#### fully restrictive - on equality

``` r
model_covAge_R<-'
# evolução
Dep_EpiDoc4~Dep_EpiDoc1+Dep_EpiDoc2
MM_EpiDoc4~MM_EpiDoc2+MM_EpiDoc1
Dep_EpiDoc2~Dep_EpiDoc1
MM_EpiDoc2~MM_EpiDoc1 

# causalidade - c/ temporalidade/lag
Dep_EpiDoc4~MM_EpiDoc2+MM_EpiDoc1
MM_EpiDoc4~Dep_EpiDoc2+Dep_EpiDoc1
Dep_EpiDoc2~MM_EpiDoc1
MM_EpiDoc2~Dep_EpiDoc1

# corrrelação mesmo t
MM_EpiDoc1~~Dep_EpiDoc1
MM_EpiDoc2~~0*Dep_EpiDoc2
MM_EpiDoc4~~0*Dep_EpiDoc4

# causalidade - idade
Dep_EpiDoc4~0*Age
Dep_EpiDoc2~d2*Age
Dep_EpiDoc1~d1*Age

MM_EpiDoc4~m4*Age
MM_EpiDoc2~m2*Age
MM_EpiDoc1~m1*Age

# educação
Dep_EpiDoc4~Education
Dep_EpiDoc2~Education
Dep_EpiDoc1~Education
MM_EpiDoc4~Education
MM_EpiDoc2~Education
MM_EpiDoc1~Education
'


res_covAge_g<-lavaan.mi::sem.mi(data=data.mi,
         model=model_covAge_R,
         estimator = "ML", test = "yuan.bentler.mplus",
         likelihood = "wishart",
           sampling.weights='ipw',group='Sex',
         orthogonal=TRUE
         )
```

    ## Warning: lavaan->lavParTable():  
    ##    using a single label per parameter in a multiple group setting implies 
    ##    imposing equality constraints across all the groups; If this is not 
    ##    intended, either remove the label(s), or use a vector of labels (one for 
    ##    each group); See the Multiple groups section in the man page of 
    ##    model.syntax.

``` r
summary(res_covAge_g,fit.meas=TRUE,standardized=TRUE,rsquare=TRUE,ci=TRUE)
```

    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.

    ## lavaan.mi object fit to 20 imputed data sets using:
    ##  - lavaan    (0.6-19)
    ##  - lavaan.mi (0.1-0)
    ## See class?lavaan.mi help page for available methods. 
    ## 
    ## Convergence information:
    ## The model converged on 20 imputed data sets.
    ## Standard errors were available for all imputations.
    ## 
    ##   Estimator                                         ML
    ##   Optimization method                           NLMINB
    ##   Number of model parameters                        72
    ##   Number of equality constraints                     5
    ## 
    ##   Number of observations per group:                   
    ##     Male                                           355
    ##     Female                                         317
    ##   Sampling weights variable                        ipw
    ## 
    ## Model Test User Model:
    ## 
    ##                                                     Standard      Scaled
    ##   Test statistic                                      19.269      16.587
    ##   Degrees of freedom                                      11          11
    ##   P-value                                              0.056       0.121
    ##   Average scaling correction factor                                1.162
    ##   Pooling method                                          D4            
    ##     Pooled statistic                              "standard"            
    ##     "yuan.bentler.mplus" correction applied            AFTER     pooling
    ## 
    ## Model Test Baseline Model:
    ## 
    ##   Test statistic                              1048.982     846.230
    ##   Degrees of freedom                                54          54
    ##   P-value                                        0.000       0.000
    ##   Scaling correction factor                                  1.240
    ## 
    ## User Model versus Baseline Model:
    ## 
    ##   Comparative Fit Index (CFI)                    0.992       0.993
    ##   Tucker-Lewis Index (TLI)                       0.959       0.965
    ##                                                                   
    ##   Robust Comparative Fit Index (CFI)                         0.993
    ##   Robust Tucker-Lewis Index (TLI)                            0.968
    ## 
    ## Loglikelihood and Information Criteria:
    ## 
    ##   Loglikelihood user model (H0)              -7333.985   -7333.985
    ##   Scaling correction factor                                  1.685
    ##       for the MLR correction                                      
    ##   Loglikelihood unrestricted model (H1)      -7325.594   -7325.594
    ##   Scaling correction factor                                  1.717
    ##       for the MLR correction                                      
    ##                                                                   
    ##   Akaike (AIC)                               14801.970   14801.970
    ##   Bayesian (BIC)                             15104.157   15104.157
    ##   Sample-size adjusted Bayesian (SABIC)      14891.427   14891.427
    ## 
    ## Root Mean Square Error of Approximation:
    ## 
    ##   RMSEA                                          0.047       0.040
    ##   90 Percent confidence interval - lower         0.000       0.000
    ##   90 Percent confidence interval - upper         0.082       0.073
    ##   P-value H_0: RMSEA <= 0.050                    0.503       0.654
    ##   P-value H_0: RMSEA >= 0.080                    0.060       0.021
    ##                                                                   
    ##   Robust RMSEA                                               0.042
    ##   90 Percent confidence interval - lower                     0.000
    ##   90 Percent confidence interval - upper                     0.081
    ##   P-value H_0: Robust RMSEA <= 0.050                         0.584
    ##   P-value H_0: Robust RMSEA >= 0.080                         0.054
    ## 
    ## Standardized Root Mean Square Residual:
    ## 
    ##   SRMR                                           0.036       0.036
    ## 
    ## Parameter Estimates:
    ## 
    ##   Standard errors                                    Sandwich
    ##   Information bread                                  Observed
    ##   Observed information based on                       Hessian
    ##                                                              
    ##   Pooled across imputations              Rubin's (1987) rules
    ##   Augment within-imputation variance     Scale by average RIV
    ##   Wald test for pooled parameters          t(df) distribution
    ## 
    ##   Pooled t statistics with df >= 1000 are displayed with
    ##   df = Inf(inity) to save space. Although the t distribution
    ##   with large df closely approximates a standard normal
    ##   distribution, exact df for reporting these t tests can be
    ##   obtained from parameterEstimates.mi() 
    ## 
    ## 
    ## 
    ## Group 1 [Male]:
    ## 
    ## Regressions:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##   Dep_EpiDoc4 ~                                                         
    ##     Dep_EpDc1         0.164    0.078    2.091  527.384    0.037    0.010
    ##     Dep_EpDc2         0.234    0.084    2.776  327.277    0.006    0.068
    ##   MM_EpiDoc4 ~                                                          
    ##     MM_EpiDc2         0.363    0.152    2.389      Inf    0.017    0.065
    ##     MM_EpiDc1         0.540    0.087    6.218      Inf    0.000    0.370
    ##   Dep_EpiDoc2 ~                                                         
    ##     Dep_EpDc1         0.404    0.070    5.781  382.954    0.000    0.267
    ##   MM_EpiDoc2 ~                                                          
    ##     MM_EpiDc1         0.110    0.048    2.279      Inf    0.023    0.015
    ##   Dep_EpiDoc4 ~                                                         
    ##     MM_EpiDc2         0.530    0.290    1.826      Inf    0.068   -0.039
    ##     MM_EpiDc1         0.046    0.207    0.224      Inf    0.823   -0.360
    ##   MM_EpiDoc4 ~                                                          
    ##     Dep_EpDc2         0.036    0.030    1.191  165.239    0.236   -0.024
    ##     Dep_EpDc1         0.004    0.029    0.124      Inf    0.902   -0.053
    ##   Dep_EpiDoc2 ~                                                         
    ##     MM_EpiDc1         0.230    0.221    1.041  358.236    0.299   -0.204
    ##   MM_EpiDoc2 ~                                                          
    ##     Dep_EpDc1         0.021    0.015    1.408      Inf    0.159   -0.008
    ##   Dep_EpiDoc4 ~                                                         
    ##     Age               0.000                                        0.000
    ##   Dep_EpiDoc2 ~                                                         
    ##     Age       (d2)    0.010    0.013    0.739   88.941    0.462   -0.016
    ##   Dep_EpiDoc1 ~                                                         
    ##     Age       (d1)    0.045    0.015    3.039      Inf    0.002    0.016
    ##   MM_EpiDoc4 ~                                                          
    ##     Age       (m4)    0.022    0.006    3.991      Inf    0.000    0.011
    ##   MM_EpiDoc2 ~                                                          
    ##     Age       (m2)    0.000    0.002    0.097      Inf    0.923   -0.004
    ##   MM_EpiDoc1 ~                                                          
    ##     Age       (m1)    0.037    0.004    9.156      Inf    0.000    0.029
    ##   Dep_EpiDoc4 ~                                                         
    ##     Education         0.162    0.165    0.983  599.403    0.326   -0.162
    ##   Dep_EpiDoc2 ~                                                         
    ##     Education         0.490    0.146    3.350  335.373    0.001    0.202
    ##   Dep_EpiDoc1 ~                                                         
    ##     Education         0.474    0.159    2.988      Inf    0.003    0.163
    ##   MM_EpiDoc4 ~                                                          
    ##     Education        -0.024    0.057   -0.431      Inf    0.666   -0.135
    ##   MM_EpiDoc2 ~                                                          
    ##     Education        -0.039    0.028   -1.390      Inf    0.164   -0.094
    ##   MM_EpiDoc1 ~                                                          
    ##     Education        -0.004    0.051   -0.078      Inf    0.938   -0.105
    ##  ci.upper   Std.lv  Std.all
    ##                            
    ##     0.318    0.164    0.174
    ##     0.399    0.234    0.241
    ##                            
    ##     0.662    0.363    0.127
    ##     0.710    0.540    0.412
    ##                            
    ##     0.542    0.404    0.416
    ##                            
    ##     0.205    0.110    0.240
    ##                            
    ##     1.099    0.530    0.084
    ##     0.453    0.046    0.016
    ##                            
    ##     0.095    0.036    0.082
    ##     0.060    0.004    0.008
    ##                            
    ##     0.663    0.230    0.077
    ##                            
    ##     0.050    0.021    0.142
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.035    0.010    0.046
    ##                            
    ##     0.073    0.045    0.208
    ##                            
    ##     0.033    0.022    0.245
    ##                            
    ##     0.004    0.000    0.006
    ##                            
    ##     0.045    0.037    0.534
    ##                            
    ##     0.485    0.162    0.060
    ##                            
    ##     0.777    0.490    0.177
    ##                            
    ##     0.785    0.474    0.167
    ##                            
    ##     0.086   -0.024   -0.020
    ##                            
    ##     0.016   -0.039   -0.092
    ##                            
    ##     0.097   -0.004   -0.004
    ## 
    ## Covariances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##  .Dep_EpiDoc1 ~~                                                        
    ##    .MM_EpiDoc1        0.551    0.216    2.550      Inf    0.011    0.128
    ##  .Dep_EpiDoc2 ~~                                                        
    ##    .MM_EpiDoc2        0.000                                        0.000
    ##  .Dep_EpiDoc4 ~~                                                        
    ##    .MM_EpiDoc4        0.000                                        0.000
    ##  ci.upper   Std.lv  Std.all
    ##                            
    ##     0.975    0.551    0.226
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.000    0.000    0.000
    ## 
    ## Intercepts:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##    .Dep_EpiDoc4       1.073    0.370    2.897  828.777    0.004    0.346
    ##    .MM_EpiDoc4       -0.451    0.232   -1.942      Inf    0.052   -0.905
    ##    .Dep_EpiDoc2      -0.261    0.566   -0.461  122.540    0.646   -1.382
    ##    .MM_EpiDoc2        0.078    0.055    1.421      Inf    0.155   -0.030
    ##    .Dep_EpiDoc1      -0.675    0.628   -1.075      Inf    0.282   -1.905
    ##    .MM_EpiDoc1       -0.984    0.175   -5.608      Inf    0.000   -1.328
    ##  ci.upper   Std.lv  Std.all
    ##     1.800    1.073    0.372
    ##     0.004   -0.451   -0.347
    ##     0.860   -0.261   -0.088
    ##     0.185    0.078    0.171
    ##     0.555   -0.675   -0.220
    ##    -0.640   -0.984   -0.991
    ## 
    ## Variances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##    .Dep_EpiDoc4       6.879    0.986    6.974      Inf    0.000    4.945
    ##    .MM_EpiDoc4        0.973    0.106    9.179      Inf    0.000    0.765
    ##    .Dep_EpiDoc2       6.201    0.862    7.190  103.268    0.000    4.490
    ##    .MM_EpiDoc2        0.187    0.051    3.664      Inf    0.000    0.087
    ##    .Dep_EpiDoc1       8.429    1.158    7.280      Inf    0.000    6.160
    ##    .MM_EpiDoc1        0.706    0.086    8.178      Inf    0.000    0.537
    ##  ci.upper   Std.lv  Std.all
    ##     8.814    6.879    0.828
    ##     1.180    0.973    0.576
    ##     7.911    6.201    0.700
    ##     0.287    0.187    0.906
    ##    10.698    8.429    0.898
    ##     0.875    0.706    0.717
    ## 
    ## R-Square:
    ##                    Estimate
    ##     Dep_EpiDoc4       0.172
    ##     MM_EpiDoc4        0.424
    ##     Dep_EpiDoc2       0.300
    ##     MM_EpiDoc2        0.094
    ##     Dep_EpiDoc1       0.102
    ##     MM_EpiDoc1        0.283
    ## 
    ## 
    ## Group 2 [Female]:
    ## 
    ## Regressions:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##   Dep_EpiDoc4 ~                                                         
    ##     Dep_EpDc1         0.342    0.075    4.579      Inf    0.000    0.195
    ##     Dep_EpDc2         0.220    0.074    2.957  943.733    0.003    0.074
    ##   MM_EpiDoc4 ~                                                          
    ##     MM_EpiDc2         0.272    0.113    2.410      Inf    0.016    0.051
    ##     MM_EpiDc1         0.394    0.078    5.053      Inf    0.000    0.241
    ##   Dep_EpiDoc2 ~                                                         
    ##     Dep_EpDc1         0.373    0.084    4.451  844.400    0.000    0.209
    ##   MM_EpiDoc2 ~                                                          
    ##     MM_EpiDc1         0.334    0.085    3.929      Inf    0.000    0.168
    ##   Dep_EpiDoc4 ~                                                         
    ##     MM_EpiDc2         0.669    0.305    2.192      Inf    0.029    0.070
    ##     MM_EpiDc1        -0.320    0.217   -1.479      Inf    0.139   -0.745
    ##   MM_EpiDoc4 ~                                                          
    ##     Dep_EpDc2         0.002    0.017    0.108  590.349    0.914   -0.031
    ##     Dep_EpDc1         0.024    0.018    1.355      Inf    0.176   -0.011
    ##   Dep_EpiDoc2 ~                                                         
    ##     MM_EpiDc1         0.086    0.246    0.350  301.262    0.727   -0.399
    ##   MM_EpiDoc2 ~                                                          
    ##     Dep_EpDc1         0.007    0.013    0.567      Inf    0.570   -0.018
    ##   Dep_EpiDoc4 ~                                                         
    ##     Age               0.000                                        0.000
    ##   Dep_EpiDoc2 ~                                                         
    ##     Age       (d2)    0.010    0.013    0.739   88.941    0.462   -0.016
    ##   Dep_EpiDoc1 ~                                                         
    ##     Age       (d1)    0.045    0.015    3.039      Inf    0.002    0.016
    ##   MM_EpiDoc4 ~                                                          
    ##     Age       (m4)    0.022    0.006    3.991      Inf    0.000    0.011
    ##   MM_EpiDoc2 ~                                                          
    ##     Age       (m2)    0.000    0.002    0.097      Inf    0.923   -0.004
    ##   MM_EpiDoc1 ~                                                          
    ##     Age       (m1)    0.037    0.004    9.156      Inf    0.000    0.029
    ##   Dep_EpiDoc4 ~                                                         
    ##     Education         0.163    0.166    0.983      Inf    0.325   -0.162
    ##   Dep_EpiDoc2 ~                                                         
    ##     Education         0.435    0.211    2.064  296.394    0.040    0.020
    ##   Dep_EpiDoc1 ~                                                         
    ##     Education         0.809    0.203    3.978      Inf    0.000    0.411
    ##   MM_EpiDoc4 ~                                                          
    ##     Education         0.059    0.060    0.971      Inf    0.331   -0.060
    ##   MM_EpiDoc2 ~                                                          
    ##     Education         0.017    0.032    0.513      Inf    0.608   -0.047
    ##   MM_EpiDoc1 ~                                                          
    ##     Education         0.059    0.056    1.064      Inf    0.287   -0.050
    ##  ci.upper   Std.lv  Std.all
    ##                            
    ##     0.488    0.342    0.356
    ##     0.365    0.220    0.225
    ##                            
    ##     0.494    0.272    0.159
    ##     0.547    0.394    0.319
    ##                            
    ##     0.538    0.373    0.379
    ##                            
    ##     0.501    0.334    0.462
    ##                            
    ##     1.267    0.669    0.135
    ##     0.104   -0.320   -0.089
    ##                            
    ##     0.035    0.002    0.005
    ##     0.059    0.024    0.073
    ##                            
    ##     0.571    0.086    0.023
    ##                            
    ##     0.033    0.007    0.038
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.035    0.010    0.030
    ##                            
    ##     0.073    0.045    0.139
    ##                            
    ##     0.033    0.022    0.209
    ##                            
    ##     0.004    0.000    0.003
    ##                            
    ##     0.045    0.037    0.431
    ##                            
    ##     0.487    0.163    0.050
    ##                            
    ##     0.850    0.435    0.132
    ##                            
    ##     1.208    0.809    0.241
    ##                            
    ##     0.177    0.059    0.053
    ##                            
    ##     0.080    0.017    0.025
    ##                            
    ##     0.169    0.059    0.066
    ## 
    ## Covariances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##  .Dep_EpiDoc1 ~~                                                        
    ##    .MM_EpiDoc1        0.936    0.317    2.956      Inf    0.003    0.315
    ##  .Dep_EpiDoc2 ~~                                                        
    ##    .MM_EpiDoc2        0.000                                        0.000
    ##  .Dep_EpiDoc4 ~~                                                        
    ##    .MM_EpiDoc4        0.000                                        0.000
    ##  ci.upper   Std.lv  Std.all
    ##                            
    ##     1.557    0.936    0.295
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.000    0.000    0.000
    ## 
    ## Intercepts:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##    .Dep_EpiDoc4       0.968    0.375    2.580      Inf    0.010    0.232
    ##    .MM_EpiDoc4       -0.370    0.218   -1.697      Inf    0.090   -0.797
    ##    .Dep_EpiDoc2       0.513    0.609    0.843  147.075    0.401   -0.690
    ##    .MM_EpiDoc2       -0.055    0.105   -0.526      Inf    0.599   -0.260
    ##    .Dep_EpiDoc1      -0.185    0.689   -0.268      Inf    0.788   -1.535
    ##    .MM_EpiDoc1       -0.990    0.185   -5.358      Inf    0.000   -1.352
    ##  ci.upper   Std.lv  Std.all
    ##     1.704    0.968    0.270
    ##     0.057   -0.370   -0.298
    ##     1.717    0.513    0.139
    ##     0.150   -0.055   -0.076
    ##     1.165   -0.185   -0.049
    ##    -0.628   -0.990   -0.987
    ## 
    ## Variances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##    .Dep_EpiDoc4       9.339    1.043    8.951      Inf    0.000    7.292
    ##    .MM_EpiDoc4        1.002    0.085   11.833      Inf    0.000    0.836
    ##    .Dep_EpiDoc2      10.759    1.348    7.980  379.008    0.000    8.108
    ##    .MM_EpiDoc2        0.403    0.062    6.509      Inf    0.000    0.282
    ##    .Dep_EpiDoc1      12.644    1.733    7.294      Inf    0.000    9.247
    ##    .MM_EpiDoc1        0.797    0.126    6.336      Inf    0.000    0.551
    ##  ci.upper   Std.lv  Std.all
    ##    11.385    9.339    0.724
    ##     1.167    1.002    0.651
    ##    13.410   10.759    0.794
    ##     0.524    0.403    0.765
    ##    16.042   12.644    0.904
    ##     1.044    0.797    0.794
    ## 
    ## R-Square:
    ##                    Estimate
    ##     Dep_EpiDoc4       0.276
    ##     MM_EpiDoc4        0.349
    ##     Dep_EpiDoc2       0.206
    ##     MM_EpiDoc2        0.235
    ##     Dep_EpiDoc1       0.096
    ##     MM_EpiDoc1        0.206

``` r
anova(resFre_g,  res_covR_g,res_covAge_g)
```

    ## Robust correction can only be applied to pooled chi-squared statistic, not F statistic. "asymptotic" was switched to TRUE.

    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.

    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.

    ## Robust correction can only be applied to pooled chi-squared statistic, not F statistic. "asymptotic" was switched to TRUE.

    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.
    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.
    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.
    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.
    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.

    ## Test statistic(s) pooled using the D4 pooling method.
    ##   Pooled statistic: "standard"
    ##   Method to robustify pooled statistic:  "yuan.bentler.mplus"
    ## 
    ##              Df   Chisq Chisq diff Df diff Pr(>Chisq)     RIV     FMI
    ## resFre_g      2  2.3128                                              
    ## res_covR_g    6  5.1723     2.4926       4    0.64597 0.34350 0.25568
    ## res_covAge_g 11 19.2691    12.5958       5    0.02748 0.11802 0.10556

difference is significant, constraining to make the effect of age being
the same among groups leads to worsening of the fit.

``` r
source('./code/parTableToCSV.R')
res<-parTableToCSV(res_covAge_g,path=('./Results/LaggedSEM_mi20/gSEX_covAge'))
```

    ## Joining with `by = join_by(group)`

    ## [1] Group2     label      estimate_p CI         Group     
    ## <0 rows> (or 0-length row.names)

``` r
res$gPartable
```

    ##                      label2 estimate_p_Male estimate_p_Female           CI_Male
    ## 1   Dep_EpiDoc4~Dep_EpiDoc1   0.164 (0.037)         0.342 (0)      [0.01:0.318]
    ## 2   Dep_EpiDoc4~Dep_EpiDoc2   0.234 (0.006)      0.22 (0.003)     [0.068:0.399]
    ## 3     MM_EpiDoc4~MM_EpiDoc2   0.363 (0.017)     0.272 (0.016)     [0.065:0.662]
    ## 4     MM_EpiDoc4~MM_EpiDoc1        0.54 (0)         0.394 (0)       [0.37:0.71]
    ## 5   Dep_EpiDoc2~Dep_EpiDoc1       0.404 (0)         0.373 (0)     [0.267:0.542]
    ## 6     MM_EpiDoc2~MM_EpiDoc1    0.11 (0.023)         0.334 (0)     [0.015:0.205]
    ## 7    Dep_EpiDoc4~MM_EpiDoc2    0.53 (0.068)     0.669 (0.029)    [-0.039:1.099]
    ## 8    Dep_EpiDoc4~MM_EpiDoc1   0.046 (0.823)     -0.32 (0.139)     [-0.36:0.453]
    ## 9    MM_EpiDoc4~Dep_EpiDoc2   0.036 (0.236)     0.002 (0.914)    [-0.024:0.095]
    ## 10   MM_EpiDoc4~Dep_EpiDoc1   0.004 (0.902)     0.024 (0.176)     [-0.053:0.06]
    ## 11   Dep_EpiDoc2~MM_EpiDoc1    0.23 (0.299)     0.086 (0.727)    [-0.204:0.663]
    ## 12   MM_EpiDoc2~Dep_EpiDoc1   0.021 (0.159)      0.007 (0.57)     [-0.008:0.05]
    ## 13  Dep_EpiDoc1~~MM_EpiDoc1   0.551 (0.011)     0.936 (0.003)     [0.128:0.975]
    ## 14  Dep_EpiDoc2~~MM_EpiDoc2          0 (NA)            0 (NA)             [0:0]
    ## 15  Dep_EpiDoc4~~MM_EpiDoc4          0 (NA)            0 (NA)             [0:0]
    ## 16          Dep_EpiDoc4~Age          0 (NA)            0 (NA)             [0:0]
    ## 17          Dep_EpiDoc2~Age    0.01 (0.462)      0.01 (0.462)    [-0.016:0.035]
    ## 18          Dep_EpiDoc1~Age   0.045 (0.002)     0.045 (0.002)     [0.016:0.073]
    ## 19           MM_EpiDoc4~Age       0.022 (0)         0.022 (0)     [0.011:0.033]
    ## 20           MM_EpiDoc2~Age       0 (0.923)         0 (0.923)    [-0.004:0.004]
    ## 21           MM_EpiDoc1~Age       0.037 (0)         0.037 (0)     [0.029:0.045]
    ## 22    Dep_EpiDoc4~Education   0.162 (0.326)     0.163 (0.325)    [-0.162:0.485]
    ## 23    Dep_EpiDoc2~Education    0.49 (0.001)      0.435 (0.04)     [0.202:0.777]
    ## 24    Dep_EpiDoc1~Education   0.474 (0.003)         0.809 (0)     [0.163:0.785]
    ## 25     MM_EpiDoc4~Education  -0.024 (0.666)     0.059 (0.331)    [-0.135:0.086]
    ## 26     MM_EpiDoc2~Education  -0.039 (0.164)     0.017 (0.608)    [-0.094:0.016]
    ## 27     MM_EpiDoc1~Education  -0.004 (0.938)     0.059 (0.287)    [-0.105:0.097]
    ## 28 Dep_EpiDoc4~~Dep_EpiDoc4       6.879 (0)         9.339 (0)     [4.945:8.814]
    ## 29   MM_EpiDoc4~~MM_EpiDoc4       0.973 (0)         1.002 (0)      [0.765:1.18]
    ## 30 Dep_EpiDoc2~~Dep_EpiDoc2       6.201 (0)        10.759 (0)      [4.49:7.911]
    ## 31   MM_EpiDoc2~~MM_EpiDoc2       0.187 (0)         0.403 (0)     [0.087:0.287]
    ## 32 Dep_EpiDoc1~~Dep_EpiDoc1       8.429 (0)        12.644 (0)     [6.16:10.698]
    ## 33   MM_EpiDoc1~~MM_EpiDoc1       0.706 (0)         0.797 (0)     [0.537:0.875]
    ## 34                 Age~~Age    204.021 (NA)      135.726 (NA) [204.021:204.021]
    ## 35           Age~~Education      6.794 (NA)         3.68 (NA)     [6.794:6.794]
    ## 36     Education~~Education      1.159 (NA)         1.24 (NA)     [1.159:1.159]
    ## 37            Dep_EpiDoc4~1   1.073 (0.004)      0.968 (0.01)       [0.346:1.8]
    ## 38             MM_EpiDoc4~1  -0.451 (0.052)      -0.37 (0.09)    [-0.905:0.004]
    ## 39            Dep_EpiDoc2~1  -0.261 (0.646)     0.513 (0.401)     [-1.382:0.86]
    ## 40             MM_EpiDoc2~1   0.078 (0.155)    -0.055 (0.599)     [-0.03:0.185]
    ## 41            Dep_EpiDoc1~1  -0.675 (0.282)    -0.185 (0.788)    [-1.905:0.555]
    ## 42             MM_EpiDoc1~1      -0.984 (0)         -0.99 (0)    [-1.328:-0.64]
    ## 43                    Age~1     44.188 (NA)       42.481 (NA)   [44.188:44.188]
    ## 44              Education~1      2.823 (NA)        2.409 (NA)     [2.823:2.823]
    ##            CI_Female
    ## 1      [0.195:0.488]
    ## 2      [0.074:0.365]
    ## 3      [0.051:0.494]
    ## 4      [0.241:0.547]
    ## 5      [0.209:0.538]
    ## 6      [0.168:0.501]
    ## 7       [0.07:1.267]
    ## 8     [-0.745:0.104]
    ## 9     [-0.031:0.035]
    ## 10    [-0.011:0.059]
    ## 11    [-0.399:0.571]
    ## 12    [-0.018:0.033]
    ## 13     [0.315:1.557]
    ## 14             [0:0]
    ## 15             [0:0]
    ## 16             [0:0]
    ## 17    [-0.016:0.035]
    ## 18     [0.016:0.073]
    ## 19     [0.011:0.033]
    ## 20    [-0.004:0.004]
    ## 21     [0.029:0.045]
    ## 22    [-0.162:0.487]
    ## 23       [0.02:0.85]
    ## 24     [0.411:1.208]
    ## 25     [-0.06:0.177]
    ## 26     [-0.047:0.08]
    ## 27     [-0.05:0.169]
    ## 28    [7.292:11.385]
    ## 29     [0.836:1.167]
    ## 30     [8.108:13.41]
    ## 31     [0.282:0.524]
    ## 32    [9.247:16.042]
    ## 33     [0.551:1.044]
    ## 34 [135.726:135.726]
    ## 35       [3.68:3.68]
    ## 36       [1.24:1.24]
    ## 37     [0.232:1.704]
    ## 38    [-0.797:0.057]
    ## 39     [-0.69:1.717]
    ## 40      [-0.26:0.15]
    ## 41    [-1.535:1.165]
    ## 42   [-1.352:-0.628]
    ## 43   [42.481:42.481]
    ## 44     [2.409:2.409]

#### freeing up MM

- we free equality on cross-lagged effect… but if gender (and contextual
  issues matter, time point 4 would differ also… )
- equal among group for dep (doubt/hypotheses under test here)

``` r
model_covAge_R<-'
# evolução
Dep_EpiDoc4~Dep_EpiDoc1+Dep_EpiDoc2
MM_EpiDoc4~MM_EpiDoc2+MM_EpiDoc1
Dep_EpiDoc2~Dep_EpiDoc1
MM_EpiDoc2~MM_EpiDoc1 

# causalidade - c/ temporalidade/lag
Dep_EpiDoc4~MM_EpiDoc2+MM_EpiDoc1
MM_EpiDoc4~Dep_EpiDoc2+Dep_EpiDoc1
Dep_EpiDoc2~MM_EpiDoc1
MM_EpiDoc2~Dep_EpiDoc1


# corrrelação mesmo t
MM_EpiDoc1~~Dep_EpiDoc1
MM_EpiDoc2~~0*Dep_EpiDoc2
MM_EpiDoc4~~0*Dep_EpiDoc4

# causalidade - idade
Dep_EpiDoc4~0*Age
Dep_EpiDoc2~ad2*Age
Dep_EpiDoc1~ad1*Age

MM_EpiDoc4~c(am4a,am4b)*Age
MM_EpiDoc2~c(am2a,am2b)*Age
MM_EpiDoc1~c(am1a,am1b)*Age



# educação
Dep_EpiDoc4~Education
Dep_EpiDoc2~Education
Dep_EpiDoc1~Education
MM_EpiDoc4~Education
MM_EpiDoc2~Education
MM_EpiDoc1~Education
'


res_covAge_g<-lavaan.mi::sem.mi(data=data.mi,
         model=model_covAge_R,
         estimator = "ML", test = "yuan.bentler.mplus",
         likelihood = "wishart",
           sampling.weights='ipw',group='Sex',
         orthogonal=TRUE
         )
```

    ## Warning: lavaan->lavParTable():  
    ##    using a single label per parameter in a multiple group setting implies 
    ##    imposing equality constraints across all the groups; If this is not 
    ##    intended, either remove the label(s), or use a vector of labels (one for 
    ##    each group); See the Multiple groups section in the man page of 
    ##    model.syntax.

``` r
summary(res_covAge_g,fit.meas=TRUE,standardized=TRUE,rsquare=TRUE,ci=TRUE)
```

    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.

    ## lavaan.mi object fit to 20 imputed data sets using:
    ##  - lavaan    (0.6-19)
    ##  - lavaan.mi (0.1-0)
    ## See class?lavaan.mi help page for available methods. 
    ## 
    ## Convergence information:
    ## The model converged on 20 imputed data sets.
    ## Standard errors were available for all imputations.
    ## 
    ##   Estimator                                         ML
    ##   Optimization method                           NLMINB
    ##   Number of model parameters                        72
    ##   Number of equality constraints                     2
    ## 
    ##   Number of observations per group:                   
    ##     Male                                           355
    ##     Female                                         317
    ##   Sampling weights variable                        ipw
    ## 
    ## Model Test User Model:
    ## 
    ##                                                     Standard      Scaled
    ##   Test statistic                                       7.701       6.748
    ##   Degrees of freedom                                       8           8
    ##   P-value                                              0.463       0.564
    ##   Average scaling correction factor                                1.141
    ##   Pooling method                                          D4            
    ##     Pooled statistic                              "standard"            
    ##     "yuan.bentler.mplus" correction applied            AFTER     pooling
    ## 
    ## Model Test Baseline Model:
    ## 
    ##   Test statistic                              1048.982     846.230
    ##   Degrees of freedom                                54          54
    ##   P-value                                        0.000       0.000
    ##   Scaling correction factor                                  1.240
    ## 
    ## User Model versus Baseline Model:
    ## 
    ##   Comparative Fit Index (CFI)                    1.000       1.000
    ##   Tucker-Lewis Index (TLI)                       1.002       1.011
    ##                                                                   
    ##   Robust Comparative Fit Index (CFI)                         1.000
    ##   Robust Tucker-Lewis Index (TLI)                            1.010
    ## 
    ## Loglikelihood and Information Criteria:
    ## 
    ##   Loglikelihood user model (H0)              -7329.078   -7329.078
    ##   Scaling correction factor                                  1.736
    ##       for the MLR correction                                      
    ##   Loglikelihood unrestricted model (H1)      -7325.594   -7325.594
    ##   Scaling correction factor                                  1.717
    ##       for the MLR correction                                      
    ##                                                                   
    ##   Akaike (AIC)                               14798.155   14798.155
    ##   Bayesian (BIC)                             15113.873   15113.873
    ##   Sample-size adjusted Bayesian (SABIC)      14891.617   14891.617
    ## 
    ## Root Mean Square Error of Approximation:
    ## 
    ##   RMSEA                                          0.000       0.000
    ##   90 Percent confidence interval - lower         0.000       0.000
    ##   90 Percent confidence interval - upper         0.063       0.055
    ##   P-value H_0: RMSEA <= 0.050                    0.872       0.922
    ##   P-value H_0: RMSEA >= 0.080                    0.008       0.003
    ##                                                                   
    ##   Robust RMSEA                                               0.000
    ##   90 Percent confidence interval - lower                     0.000
    ##   90 Percent confidence interval - upper                     0.061
    ##   P-value H_0: Robust RMSEA <= 0.050                         0.892
    ##   P-value H_0: Robust RMSEA >= 0.080                         0.008
    ## 
    ## Standardized Root Mean Square Residual:
    ## 
    ##   SRMR                                           0.025       0.025
    ## 
    ## Parameter Estimates:
    ## 
    ##   Standard errors                                    Sandwich
    ##   Information bread                                  Observed
    ##   Observed information based on                       Hessian
    ##                                                              
    ##   Pooled across imputations              Rubin's (1987) rules
    ##   Augment within-imputation variance     Scale by average RIV
    ##   Wald test for pooled parameters          t(df) distribution
    ## 
    ##   Pooled t statistics with df >= 1000 are displayed with
    ##   df = Inf(inity) to save space. Although the t distribution
    ##   with large df closely approximates a standard normal
    ##   distribution, exact df for reporting these t tests can be
    ##   obtained from parameterEstimates.mi() 
    ## 
    ## 
    ## 
    ## Group 1 [Male]:
    ## 
    ## Regressions:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##   Dep_EpiDoc4 ~                                                         
    ##     Dp_EpD1           0.164    0.078    2.094  527.385    0.037    0.010
    ##     Dp_EpD2           0.234    0.084    2.780  327.278    0.006    0.068
    ##   MM_EpiDoc4 ~                                                          
    ##     MM_EpD2           0.367    0.153    2.397      Inf    0.017    0.067
    ##     MM_EpD1           0.560    0.089    6.286      Inf    0.000    0.386
    ##   Dep_EpiDoc2 ~                                                         
    ##     Dp_EpD1           0.404    0.070    5.791  382.956    0.000    0.267
    ##   MM_EpiDoc2 ~                                                          
    ##     MM_EpD1           0.102    0.046    2.210      Inf    0.027    0.012
    ##   Dep_EpiDoc4 ~                                                         
    ##     MM_EpD2           0.530    0.290    1.829      Inf    0.067   -0.038
    ##     MM_EpD1           0.046    0.207    0.225      Inf    0.822   -0.359
    ##   MM_EpiDoc4 ~                                                          
    ##     Dp_EpD2           0.036    0.030    1.207  179.079    0.229   -0.023
    ##     Dp_EpD1           0.005    0.029    0.157      Inf    0.875   -0.052
    ##   Dep_EpiDoc2 ~                                                         
    ##     MM_EpD1           0.230    0.220    1.043  358.233    0.298   -0.203
    ##   MM_EpiDoc2 ~                                                          
    ##     Dp_EpD1           0.021    0.015    1.373      Inf    0.170   -0.009
    ##   Dep_EpiDoc4 ~                                                         
    ##     Age               0.000                                        0.000
    ##   Dep_EpiDoc2 ~                                                         
    ##     Age      (ad2)    0.010    0.013    0.740   88.941    0.461   -0.016
    ##   Dep_EpiDoc1 ~                                                         
    ##     Age      (ad1)    0.045    0.015    3.096      Inf    0.002    0.017
    ##   MM_EpiDoc4 ~                                                          
    ##     Age     (am4a)    0.019    0.007    2.541      Inf    0.011    0.004
    ##   MM_EpiDoc2 ~                                                          
    ##     Age     (am2a)    0.002    0.002    0.676      Inf    0.499   -0.003
    ##   MM_EpiDoc1 ~                                                          
    ##     Age     (am1a)    0.031    0.005    6.380      Inf    0.000    0.022
    ##   Dep_EpiDoc4 ~                                                         
    ##     Educatn           0.162    0.164    0.984  599.365    0.325   -0.161
    ##   Dep_EpiDoc2 ~                                                         
    ##     Educatn           0.490    0.146    3.356  335.395    0.001    0.203
    ##   Dep_EpiDoc1 ~                                                         
    ##     Educatn           0.472    0.159    2.975      Inf    0.003    0.161
    ##   MM_EpiDoc4 ~                                                          
    ##     Educatn          -0.008    0.056   -0.146      Inf    0.884   -0.118
    ##   MM_EpiDoc2 ~                                                          
    ##     Educatn          -0.045    0.029   -1.534      Inf    0.125   -0.102
    ##   MM_EpiDoc1 ~                                                          
    ##     Educatn           0.030    0.054    0.553      Inf    0.580   -0.075
    ##  ci.upper   Std.lv  Std.all
    ##                            
    ##     0.318    0.164    0.174
    ##     0.399    0.234    0.241
    ##                            
    ##     0.668    0.367    0.131
    ##     0.735    0.560    0.421
    ##                            
    ##     0.541    0.404    0.417
    ##                            
    ##     0.193    0.102    0.216
    ##                            
    ##     1.098    0.530    0.083
    ##     0.452    0.046    0.015
    ##                            
    ##     0.094    0.036    0.083
    ##     0.061    0.005    0.011
    ##                            
    ##     0.663    0.230    0.074
    ##                            
    ##     0.050    0.021    0.140
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.035    0.010    0.046
    ##                            
    ##     0.073    0.045    0.210
    ##                            
    ##     0.033    0.019    0.209
    ##                            
    ##     0.006    0.002    0.049
    ##                            
    ##     0.041    0.031    0.468
    ##                            
    ##     0.485    0.162    0.060
    ##                            
    ##     0.777    0.490    0.177
    ##                            
    ##     0.782    0.472    0.166
    ##                            
    ##     0.102   -0.008   -0.007
    ##                            
    ##     0.012   -0.045   -0.107
    ##                            
    ##     0.135    0.030    0.033
    ## 
    ## Covariances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##  .Dep_EpiDoc1 ~~                                                        
    ##    .MM_EpiDoc1        0.548    0.214    2.558      Inf    0.011    0.128
    ##  .Dep_EpiDoc2 ~~                                                        
    ##    .MM_EpiDoc2        0.000                                        0.000
    ##  .Dep_EpiDoc4 ~~                                                        
    ##    .MM_EpiDoc4        0.000                                        0.000
    ##  ci.upper   Std.lv  Std.all
    ##                            
    ##     0.968    0.548    0.226
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.000    0.000    0.000
    ## 
    ## Intercepts:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##    .Dep_EpiDoc4       1.073    0.370    2.902  828.763    0.004    0.347
    ##    .MM_EpiDoc4       -0.352    0.271   -1.301      Inf    0.193   -0.882
    ##    .Dep_EpiDoc2      -0.261    0.565   -0.461  122.543    0.645   -1.380
    ##    .MM_EpiDoc2        0.041    0.054    0.750      Inf    0.453   -0.066
    ##    .Dep_EpiDoc1      -0.687    0.623   -1.102      Inf    0.270   -1.909
    ##    .MM_EpiDoc1       -0.825    0.186   -4.444      Inf    0.000   -1.189
    ##  ci.upper   Std.lv  Std.all
    ##     1.799    1.073    0.372
    ##     0.178   -0.352   -0.276
    ##     0.858   -0.261   -0.088
    ##     0.147    0.041    0.090
    ##     0.535   -0.687   -0.224
    ##    -0.461   -0.825   -0.863
    ## 
    ## Variances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##    .Dep_EpiDoc4       6.879    0.985    6.986      Inf    0.000    4.948
    ##    .MM_EpiDoc4        0.971    0.105    9.213      Inf    0.000    0.764
    ##    .Dep_EpiDoc2       6.201    0.861    7.202  103.268    0.000    4.493
    ##    .MM_EpiDoc2        0.187    0.051    3.676      Inf    0.000    0.087
    ##    .Dep_EpiDoc1       8.430    1.156    7.292      Inf    0.000    6.164
    ##    .MM_EpiDoc1        0.700    0.086    8.162      Inf    0.000    0.532
    ##  ci.upper   Std.lv  Std.all
    ##     8.811    6.879    0.829
    ##     1.177    0.971    0.599
    ##     7.908    6.201    0.701
    ##     0.286    0.187    0.910
    ##    10.695    8.430    0.898
    ##     0.868    0.700    0.766
    ## 
    ## R-Square:
    ##                    Estimate
    ##     Dep_EpiDoc4       0.171
    ##     MM_EpiDoc4        0.401
    ##     Dep_EpiDoc2       0.299
    ##     MM_EpiDoc2        0.090
    ##     Dep_EpiDoc1       0.102
    ##     MM_EpiDoc1        0.234
    ## 
    ## 
    ## Group 2 [Female]:
    ## 
    ## Regressions:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##   Dep_EpiDoc4 ~                                                         
    ##     Dp_EpD1           0.342    0.074    4.587      Inf    0.000    0.196
    ##     Dp_EpD2           0.220    0.074    2.962  943.735    0.003    0.074
    ##   MM_EpiDoc4 ~                                                          
    ##     MM_EpD2           0.282    0.113    2.498      Inf    0.012    0.061
    ##     MM_EpD1           0.350    0.087    4.030      Inf    0.000    0.180
    ##   Dep_EpiDoc2 ~                                                         
    ##     Dp_EpD1           0.373    0.084    4.459  844.400    0.000    0.209
    ##   MM_EpiDoc2 ~                                                          
    ##     MM_EpD1           0.366    0.096    3.803      Inf    0.000    0.177
    ##   Dep_EpiDoc4 ~                                                         
    ##     MM_EpD2           0.669    0.305    2.195      Inf    0.028    0.071
    ##     MM_EpD1          -0.320    0.216   -1.481      Inf    0.139   -0.744
    ##   MM_EpiDoc4 ~                                                          
    ##     Dp_EpD2          -0.001    0.017   -0.042  494.986    0.967   -0.034
    ##     Dp_EpD1           0.025    0.018    1.411      Inf    0.158   -0.010
    ##   Dep_EpiDoc2 ~                                                         
    ##     MM_EpD1           0.086    0.246    0.350  301.262    0.726   -0.398
    ##   MM_EpiDoc2 ~                                                          
    ##     Dp_EpD1           0.007    0.013    0.565      Inf    0.572   -0.018
    ##   Dep_EpiDoc4 ~                                                         
    ##     Age               0.000                                        0.000
    ##   Dep_EpiDoc2 ~                                                         
    ##     Age      (ad2)    0.010    0.013    0.740   88.941    0.461   -0.016
    ##   Dep_EpiDoc1 ~                                                         
    ##     Age      (ad1)    0.045    0.015    3.096      Inf    0.002    0.017
    ##   MM_EpiDoc4 ~                                                          
    ##     Age     (am4b)    0.030    0.007    4.352      Inf    0.000    0.016
    ##   MM_EpiDoc2 ~                                                          
    ##     Age     (am2b)   -0.005    0.005   -1.193      Inf    0.233   -0.014
    ##   MM_EpiDoc1 ~                                                          
    ##     Age     (am1b)    0.048    0.005    9.497      Inf    0.000    0.038
    ##   Dep_EpiDoc4 ~                                                         
    ##     Educatn           0.163    0.165    0.985      Inf    0.325   -0.161
    ##   Dep_EpiDoc2 ~                                                         
    ##     Educatn           0.435    0.210    2.067  296.399    0.040    0.021
    ##   Dep_EpiDoc1 ~                                                         
    ##     Educatn           0.808    0.203    3.987      Inf    0.000    0.411
    ##   MM_EpiDoc4 ~                                                          
    ##     Educatn           0.045    0.060    0.752      Inf    0.452   -0.072
    ##   MM_EpiDoc2 ~                                                          
    ##     Educatn           0.028    0.033    0.859      Inf    0.390   -0.036
    ##   MM_EpiDoc1 ~                                                          
    ##     Educatn           0.028    0.050    0.549      Inf    0.583   -0.071
    ##  ci.upper   Std.lv  Std.all
    ##                            
    ##     0.488    0.342    0.356
    ##     0.365    0.220    0.225
    ##                            
    ##     0.503    0.282    0.163
    ##     0.521    0.350    0.290
    ##                            
    ##     0.538    0.373    0.379
    ##                            
    ##     0.555    0.366    0.524
    ##                            
    ##     1.266    0.669    0.137
    ##     0.104   -0.320   -0.094
    ##                            
    ##     0.033   -0.001   -0.002
    ##     0.060    0.025    0.074
    ##                            
    ##     0.570    0.086    0.025
    ##                            
    ##     0.033    0.007    0.038
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.035    0.010    0.030
    ##                            
    ##     0.073    0.045    0.140
    ##                            
    ##     0.043    0.030    0.271
    ##                            
    ##     0.003   -0.005   -0.086
    ##                            
    ##     0.058    0.048    0.531
    ##                            
    ##     0.487    0.163    0.050
    ##                            
    ##     0.849    0.435    0.132
    ##                            
    ##     1.205    0.808    0.240
    ##                            
    ##     0.162    0.045    0.039
    ##                            
    ##     0.092    0.028    0.042
    ##                            
    ##     0.126    0.028    0.029
    ## 
    ## Covariances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##  .Dep_EpiDoc1 ~~                                                        
    ##    .MM_EpiDoc1        0.922    0.308    2.991      Inf    0.003    0.318
    ##  .Dep_EpiDoc2 ~~                                                        
    ##    .MM_EpiDoc2        0.000                                        0.000
    ##  .Dep_EpiDoc4 ~~                                                        
    ##    .MM_EpiDoc4        0.000                                        0.000
    ##  ci.upper   Std.lv  Std.all
    ##                            
    ##     1.525    0.922    0.293
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.000    0.000    0.000
    ## 
    ## Intercepts:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##    .Dep_EpiDoc4       0.968    0.375    2.584      Inf    0.010    0.234
    ##    .MM_EpiDoc4       -0.611    0.241   -2.530      Inf    0.011   -1.084
    ##    .Dep_EpiDoc2       0.513    0.608    0.844  147.075    0.400   -0.688
    ##    .MM_EpiDoc2        0.133    0.158    0.840      Inf    0.401   -0.177
    ##    .Dep_EpiDoc1      -0.201    0.687   -0.292      Inf    0.770   -1.547
    ##    .MM_EpiDoc1       -1.370    0.225   -6.082      Inf    0.000   -1.811
    ##  ci.upper   Std.lv  Std.all
    ##     1.703    0.968    0.270
    ##    -0.138   -0.611   -0.481
    ##     1.715    0.513    0.139
    ##     0.443    0.133    0.181
    ##     1.146   -0.201   -0.054
    ##    -0.928   -1.370   -1.305
    ## 
    ## Variances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##    .Dep_EpiDoc4       9.339    1.042    8.966      Inf    0.000    7.295
    ##    .MM_EpiDoc4        0.997    0.084   11.908      Inf    0.000    0.833
    ##    .Dep_EpiDoc2      10.759    1.346    7.994  379.009    0.000    8.112
    ##    .MM_EpiDoc2        0.400    0.061    6.568      Inf    0.000    0.281
    ##    .Dep_EpiDoc1      12.643    1.730    7.310      Inf    0.000    9.253
    ##    .MM_EpiDoc1        0.781    0.120    6.489      Inf    0.000    0.545
    ##  ci.upper   Std.lv  Std.all
    ##    11.382    9.339    0.724
    ##     1.161    0.997    0.618
    ##    13.405   10.759    0.794
    ##     0.520    0.400    0.744
    ##    16.033   12.643    0.903
    ##     1.017    0.781    0.708
    ## 
    ## R-Square:
    ##                    Estimate
    ##     Dep_EpiDoc4       0.276
    ##     MM_EpiDoc4        0.382
    ##     Dep_EpiDoc2       0.206
    ##     MM_EpiDoc2        0.256
    ##     Dep_EpiDoc1       0.097
    ##     MM_EpiDoc1        0.292

``` r
anova(resFre_g,  res_covR_g,res_covAge_g)
```

    ## Robust correction can only be applied to pooled chi-squared statistic, not F statistic. "asymptotic" was switched to TRUE.

    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.

    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.

    ## Robust correction can only be applied to pooled chi-squared statistic, not F statistic. "asymptotic" was switched to TRUE.

    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.
    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.
    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.
    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.
    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.

    ## Test statistic(s) pooled using the D4 pooling method.
    ##   Pooled statistic: "standard"
    ##   Method to robustify pooled statistic:  "yuan.bentler.mplus"
    ## 
    ##              Df  Chisq Chisq diff Df diff Pr(>Chisq)     RIV     FMI
    ## resFre_g      2 2.3128                                              
    ## res_covR_g    6 5.1723     2.4926       4    0.64597 0.34350 0.25568
    ## res_covAge_g  8 7.7007     2.1111       2    0.34799 0.35321 0.26102

difference is still significant, constraining to make the effect of age
being the same among groups leads to worsening of the fit.

``` r
source('./code/parTableToCSV.R')
res<-parTableToCSV(res_covAge_g,path=('./Results/LaggedSEM_mi20/gSEX_covAge'))
```

    ## Joining with `by = join_by(group)`

    ## [1] Group2     label      estimate_p CI         Group     
    ## <0 rows> (or 0-length row.names)

``` r
res$gPartable
```

    ##                      label2 estimate_p_Male estimate_p_Female           CI_Male
    ## 1   Dep_EpiDoc4~Dep_EpiDoc1   0.164 (0.037)         0.342 (0)      [0.01:0.318]
    ## 2   Dep_EpiDoc4~Dep_EpiDoc2   0.234 (0.006)      0.22 (0.003)     [0.068:0.399]
    ## 3     MM_EpiDoc4~MM_EpiDoc2   0.367 (0.017)     0.282 (0.012)     [0.067:0.668]
    ## 4     MM_EpiDoc4~MM_EpiDoc1        0.56 (0)          0.35 (0)     [0.386:0.735]
    ## 5   Dep_EpiDoc2~Dep_EpiDoc1       0.404 (0)         0.373 (0)     [0.267:0.541]
    ## 6     MM_EpiDoc2~MM_EpiDoc1   0.102 (0.027)         0.366 (0)     [0.012:0.193]
    ## 7    Dep_EpiDoc4~MM_EpiDoc2    0.53 (0.067)     0.669 (0.028)    [-0.038:1.098]
    ## 8    Dep_EpiDoc4~MM_EpiDoc1   0.046 (0.822)     -0.32 (0.139)    [-0.359:0.452]
    ## 9    MM_EpiDoc4~Dep_EpiDoc2   0.036 (0.229)    -0.001 (0.967)    [-0.023:0.094]
    ## 10   MM_EpiDoc4~Dep_EpiDoc1   0.005 (0.875)     0.025 (0.158)    [-0.052:0.061]
    ## 11   Dep_EpiDoc2~MM_EpiDoc1    0.23 (0.298)     0.086 (0.726)    [-0.203:0.663]
    ## 12   MM_EpiDoc2~Dep_EpiDoc1    0.021 (0.17)     0.007 (0.572)     [-0.009:0.05]
    ## 13  Dep_EpiDoc1~~MM_EpiDoc1   0.548 (0.011)     0.922 (0.003)     [0.128:0.968]
    ## 14  Dep_EpiDoc2~~MM_EpiDoc2          0 (NA)            0 (NA)             [0:0]
    ## 15  Dep_EpiDoc4~~MM_EpiDoc4          0 (NA)            0 (NA)             [0:0]
    ## 16          Dep_EpiDoc4~Age          0 (NA)            0 (NA)             [0:0]
    ## 17          Dep_EpiDoc2~Age    0.01 (0.461)      0.01 (0.461)    [-0.016:0.035]
    ## 18          Dep_EpiDoc1~Age   0.045 (0.002)     0.045 (0.002)     [0.017:0.073]
    ## 19           MM_EpiDoc4~Age   0.019 (0.011)          0.03 (0)     [0.004:0.033]
    ## 20           MM_EpiDoc2~Age   0.002 (0.499)    -0.005 (0.233)    [-0.003:0.006]
    ## 21           MM_EpiDoc1~Age       0.031 (0)         0.048 (0)     [0.022:0.041]
    ## 22    Dep_EpiDoc4~Education   0.162 (0.325)     0.163 (0.325)    [-0.161:0.485]
    ## 23    Dep_EpiDoc2~Education    0.49 (0.001)      0.435 (0.04)     [0.203:0.777]
    ## 24    Dep_EpiDoc1~Education   0.472 (0.003)         0.808 (0)     [0.161:0.782]
    ## 25     MM_EpiDoc4~Education  -0.008 (0.884)     0.045 (0.452)    [-0.118:0.102]
    ## 26     MM_EpiDoc2~Education  -0.045 (0.125)      0.028 (0.39)    [-0.102:0.012]
    ## 27     MM_EpiDoc1~Education     0.03 (0.58)     0.028 (0.583)    [-0.075:0.135]
    ## 28 Dep_EpiDoc4~~Dep_EpiDoc4       6.879 (0)         9.339 (0)     [4.948:8.811]
    ## 29   MM_EpiDoc4~~MM_EpiDoc4       0.971 (0)         0.997 (0)     [0.764:1.177]
    ## 30 Dep_EpiDoc2~~Dep_EpiDoc2       6.201 (0)        10.759 (0)     [4.493:7.908]
    ## 31   MM_EpiDoc2~~MM_EpiDoc2       0.187 (0)           0.4 (0)     [0.087:0.286]
    ## 32 Dep_EpiDoc1~~Dep_EpiDoc1        8.43 (0)        12.643 (0)    [6.164:10.695]
    ## 33   MM_EpiDoc1~~MM_EpiDoc1         0.7 (0)         0.781 (0)     [0.532:0.868]
    ## 34                 Age~~Age    204.021 (NA)      135.726 (NA) [204.021:204.021]
    ## 35           Age~~Education      6.794 (NA)         3.68 (NA)     [6.794:6.794]
    ## 36     Education~~Education      1.159 (NA)         1.24 (NA)     [1.159:1.159]
    ## 37            Dep_EpiDoc4~1   1.073 (0.004)      0.968 (0.01)     [0.347:1.799]
    ## 38             MM_EpiDoc4~1  -0.352 (0.193)    -0.611 (0.011)    [-0.882:0.178]
    ## 39            Dep_EpiDoc2~1  -0.261 (0.645)       0.513 (0.4)     [-1.38:0.858]
    ## 40             MM_EpiDoc2~1   0.041 (0.453)     0.133 (0.401)    [-0.066:0.147]
    ## 41            Dep_EpiDoc1~1   -0.687 (0.27)     -0.201 (0.77)    [-1.909:0.535]
    ## 42             MM_EpiDoc1~1      -0.825 (0)         -1.37 (0)   [-1.189:-0.461]
    ## 43                    Age~1     44.188 (NA)       42.481 (NA)   [44.188:44.188]
    ## 44              Education~1      2.823 (NA)        2.409 (NA)     [2.823:2.823]
    ##            CI_Female
    ## 1      [0.196:0.488]
    ## 2      [0.074:0.365]
    ## 3      [0.061:0.503]
    ## 4       [0.18:0.521]
    ## 5      [0.209:0.538]
    ## 6      [0.177:0.555]
    ## 7      [0.071:1.266]
    ## 8     [-0.744:0.104]
    ## 9     [-0.034:0.033]
    ## 10      [-0.01:0.06]
    ## 11     [-0.398:0.57]
    ## 12    [-0.018:0.033]
    ## 13     [0.318:1.525]
    ## 14             [0:0]
    ## 15             [0:0]
    ## 16             [0:0]
    ## 17    [-0.016:0.035]
    ## 18     [0.017:0.073]
    ## 19     [0.016:0.043]
    ## 20    [-0.014:0.003]
    ## 21     [0.038:0.058]
    ## 22    [-0.161:0.487]
    ## 23     [0.021:0.849]
    ## 24     [0.411:1.205]
    ## 25    [-0.072:0.162]
    ## 26    [-0.036:0.092]
    ## 27    [-0.071:0.126]
    ## 28    [7.295:11.382]
    ## 29     [0.833:1.161]
    ## 30    [8.112:13.405]
    ## 31      [0.281:0.52]
    ## 32    [9.253:16.033]
    ## 33     [0.545:1.017]
    ## 34 [135.726:135.726]
    ## 35       [3.68:3.68]
    ## 36       [1.24:1.24]
    ## 37     [0.234:1.703]
    ## 38   [-1.084:-0.138]
    ## 39    [-0.688:1.715]
    ## 40    [-0.177:0.443]
    ## 41    [-1.547:1.146]
    ## 42   [-1.811:-0.928]
    ## 43   [42.481:42.481]
    ## 44     [2.409:2.409]

``` r
model_covAge_R<-'
# evolução
Dep_EpiDoc4~Dep_EpiDoc1+Dep_EpiDoc2
MM_EpiDoc4~MM_EpiDoc2+MM_EpiDoc1
Dep_EpiDoc2~Dep_EpiDoc1
MM_EpiDoc2~MM_EpiDoc1 

# causalidade - c/ temporalidade/lag
Dep_EpiDoc4~MM_EpiDoc2+MM_EpiDoc1
MM_EpiDoc4~Dep_EpiDoc2+Dep_EpiDoc1
Dep_EpiDoc2~MM_EpiDoc1
MM_EpiDoc2~Dep_EpiDoc1


# corrrelação mesmo t
MM_EpiDoc1~~Dep_EpiDoc1
MM_EpiDoc2~~0*Dep_EpiDoc2
MM_EpiDoc4~~0*Dep_EpiDoc4

# causalidade - idade
Dep_EpiDoc4~0*Age
Dep_EpiDoc2~ad2*Age
Dep_EpiDoc1~ad1*Age

MM_EpiDoc4~c(am4a,am4b)*Age
MM_EpiDoc2~c(am2a,am2b)*Age
MM_EpiDoc1~c(am1a,am1b)*Age



# educação
Dep_EpiDoc4~Education
Dep_EpiDoc2~Education
Dep_EpiDoc1~Education
MM_EpiDoc4~Education
MM_EpiDoc2~Education
MM_EpiDoc1~Education
'


res_covAge_g<-lavaan.mi::sem.mi(data=data.mi,
         model=model_covAge_R,
         estimator = "ML", test = "yuan.bentler.mplus",
         likelihood = "wishart",
           sampling.weights='ipw',group='Sex',
         orthogonal=TRUE
         )
```

    ## Warning: lavaan->lavParTable():  
    ##    using a single label per parameter in a multiple group setting implies 
    ##    imposing equality constraints across all the groups; If this is not 
    ##    intended, either remove the label(s), or use a vector of labels (one for 
    ##    each group); See the Multiple groups section in the man page of 
    ##    model.syntax.

``` r
summary(res_covAge_g,fit.meas=TRUE,standardized=TRUE,rsquare=TRUE,ci=TRUE)
```

    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.

    ## lavaan.mi object fit to 20 imputed data sets using:
    ##  - lavaan    (0.6-19)
    ##  - lavaan.mi (0.1-0)
    ## See class?lavaan.mi help page for available methods. 
    ## 
    ## Convergence information:
    ## The model converged on 20 imputed data sets.
    ## Standard errors were available for all imputations.
    ## 
    ##   Estimator                                         ML
    ##   Optimization method                           NLMINB
    ##   Number of model parameters                        72
    ##   Number of equality constraints                     2
    ## 
    ##   Number of observations per group:                   
    ##     Male                                           355
    ##     Female                                         317
    ##   Sampling weights variable                        ipw
    ## 
    ## Model Test User Model:
    ## 
    ##                                                     Standard      Scaled
    ##   Test statistic                                       7.701       6.748
    ##   Degrees of freedom                                       8           8
    ##   P-value                                              0.463       0.564
    ##   Average scaling correction factor                                1.141
    ##   Pooling method                                          D4            
    ##     Pooled statistic                              "standard"            
    ##     "yuan.bentler.mplus" correction applied            AFTER     pooling
    ## 
    ## Model Test Baseline Model:
    ## 
    ##   Test statistic                              1048.982     846.230
    ##   Degrees of freedom                                54          54
    ##   P-value                                        0.000       0.000
    ##   Scaling correction factor                                  1.240
    ## 
    ## User Model versus Baseline Model:
    ## 
    ##   Comparative Fit Index (CFI)                    1.000       1.000
    ##   Tucker-Lewis Index (TLI)                       1.002       1.011
    ##                                                                   
    ##   Robust Comparative Fit Index (CFI)                         1.000
    ##   Robust Tucker-Lewis Index (TLI)                            1.010
    ## 
    ## Loglikelihood and Information Criteria:
    ## 
    ##   Loglikelihood user model (H0)              -7329.078   -7329.078
    ##   Scaling correction factor                                  1.736
    ##       for the MLR correction                                      
    ##   Loglikelihood unrestricted model (H1)      -7325.594   -7325.594
    ##   Scaling correction factor                                  1.717
    ##       for the MLR correction                                      
    ##                                                                   
    ##   Akaike (AIC)                               14798.155   14798.155
    ##   Bayesian (BIC)                             15113.873   15113.873
    ##   Sample-size adjusted Bayesian (SABIC)      14891.617   14891.617
    ## 
    ## Root Mean Square Error of Approximation:
    ## 
    ##   RMSEA                                          0.000       0.000
    ##   90 Percent confidence interval - lower         0.000       0.000
    ##   90 Percent confidence interval - upper         0.063       0.055
    ##   P-value H_0: RMSEA <= 0.050                    0.872       0.922
    ##   P-value H_0: RMSEA >= 0.080                    0.008       0.003
    ##                                                                   
    ##   Robust RMSEA                                               0.000
    ##   90 Percent confidence interval - lower                     0.000
    ##   90 Percent confidence interval - upper                     0.061
    ##   P-value H_0: Robust RMSEA <= 0.050                         0.892
    ##   P-value H_0: Robust RMSEA >= 0.080                         0.008
    ## 
    ## Standardized Root Mean Square Residual:
    ## 
    ##   SRMR                                           0.025       0.025
    ## 
    ## Parameter Estimates:
    ## 
    ##   Standard errors                                    Sandwich
    ##   Information bread                                  Observed
    ##   Observed information based on                       Hessian
    ##                                                              
    ##   Pooled across imputations              Rubin's (1987) rules
    ##   Augment within-imputation variance     Scale by average RIV
    ##   Wald test for pooled parameters          t(df) distribution
    ## 
    ##   Pooled t statistics with df >= 1000 are displayed with
    ##   df = Inf(inity) to save space. Although the t distribution
    ##   with large df closely approximates a standard normal
    ##   distribution, exact df for reporting these t tests can be
    ##   obtained from parameterEstimates.mi() 
    ## 
    ## 
    ## 
    ## Group 1 [Male]:
    ## 
    ## Regressions:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##   Dep_EpiDoc4 ~                                                         
    ##     Dp_EpD1           0.164    0.078    2.094  527.385    0.037    0.010
    ##     Dp_EpD2           0.234    0.084    2.780  327.278    0.006    0.068
    ##   MM_EpiDoc4 ~                                                          
    ##     MM_EpD2           0.367    0.153    2.397      Inf    0.017    0.067
    ##     MM_EpD1           0.560    0.089    6.286      Inf    0.000    0.386
    ##   Dep_EpiDoc2 ~                                                         
    ##     Dp_EpD1           0.404    0.070    5.791  382.956    0.000    0.267
    ##   MM_EpiDoc2 ~                                                          
    ##     MM_EpD1           0.102    0.046    2.210      Inf    0.027    0.012
    ##   Dep_EpiDoc4 ~                                                         
    ##     MM_EpD2           0.530    0.290    1.829      Inf    0.067   -0.038
    ##     MM_EpD1           0.046    0.207    0.225      Inf    0.822   -0.359
    ##   MM_EpiDoc4 ~                                                          
    ##     Dp_EpD2           0.036    0.030    1.207  179.079    0.229   -0.023
    ##     Dp_EpD1           0.005    0.029    0.157      Inf    0.875   -0.052
    ##   Dep_EpiDoc2 ~                                                         
    ##     MM_EpD1           0.230    0.220    1.043  358.233    0.298   -0.203
    ##   MM_EpiDoc2 ~                                                          
    ##     Dp_EpD1           0.021    0.015    1.373      Inf    0.170   -0.009
    ##   Dep_EpiDoc4 ~                                                         
    ##     Age               0.000                                        0.000
    ##   Dep_EpiDoc2 ~                                                         
    ##     Age      (ad2)    0.010    0.013    0.740   88.941    0.461   -0.016
    ##   Dep_EpiDoc1 ~                                                         
    ##     Age      (ad1)    0.045    0.015    3.096      Inf    0.002    0.017
    ##   MM_EpiDoc4 ~                                                          
    ##     Age     (am4a)    0.019    0.007    2.541      Inf    0.011    0.004
    ##   MM_EpiDoc2 ~                                                          
    ##     Age     (am2a)    0.002    0.002    0.676      Inf    0.499   -0.003
    ##   MM_EpiDoc1 ~                                                          
    ##     Age     (am1a)    0.031    0.005    6.380      Inf    0.000    0.022
    ##   Dep_EpiDoc4 ~                                                         
    ##     Educatn           0.162    0.164    0.984  599.365    0.325   -0.161
    ##   Dep_EpiDoc2 ~                                                         
    ##     Educatn           0.490    0.146    3.356  335.395    0.001    0.203
    ##   Dep_EpiDoc1 ~                                                         
    ##     Educatn           0.472    0.159    2.975      Inf    0.003    0.161
    ##   MM_EpiDoc4 ~                                                          
    ##     Educatn          -0.008    0.056   -0.146      Inf    0.884   -0.118
    ##   MM_EpiDoc2 ~                                                          
    ##     Educatn          -0.045    0.029   -1.534      Inf    0.125   -0.102
    ##   MM_EpiDoc1 ~                                                          
    ##     Educatn           0.030    0.054    0.553      Inf    0.580   -0.075
    ##  ci.upper   Std.lv  Std.all
    ##                            
    ##     0.318    0.164    0.174
    ##     0.399    0.234    0.241
    ##                            
    ##     0.668    0.367    0.131
    ##     0.735    0.560    0.421
    ##                            
    ##     0.541    0.404    0.417
    ##                            
    ##     0.193    0.102    0.216
    ##                            
    ##     1.098    0.530    0.083
    ##     0.452    0.046    0.015
    ##                            
    ##     0.094    0.036    0.083
    ##     0.061    0.005    0.011
    ##                            
    ##     0.663    0.230    0.074
    ##                            
    ##     0.050    0.021    0.140
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.035    0.010    0.046
    ##                            
    ##     0.073    0.045    0.210
    ##                            
    ##     0.033    0.019    0.209
    ##                            
    ##     0.006    0.002    0.049
    ##                            
    ##     0.041    0.031    0.468
    ##                            
    ##     0.485    0.162    0.060
    ##                            
    ##     0.777    0.490    0.177
    ##                            
    ##     0.782    0.472    0.166
    ##                            
    ##     0.102   -0.008   -0.007
    ##                            
    ##     0.012   -0.045   -0.107
    ##                            
    ##     0.135    0.030    0.033
    ## 
    ## Covariances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##  .Dep_EpiDoc1 ~~                                                        
    ##    .MM_EpiDoc1        0.548    0.214    2.558      Inf    0.011    0.128
    ##  .Dep_EpiDoc2 ~~                                                        
    ##    .MM_EpiDoc2        0.000                                        0.000
    ##  .Dep_EpiDoc4 ~~                                                        
    ##    .MM_EpiDoc4        0.000                                        0.000
    ##  ci.upper   Std.lv  Std.all
    ##                            
    ##     0.968    0.548    0.226
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.000    0.000    0.000
    ## 
    ## Intercepts:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##    .Dep_EpiDoc4       1.073    0.370    2.902  828.763    0.004    0.347
    ##    .MM_EpiDoc4       -0.352    0.271   -1.301      Inf    0.193   -0.882
    ##    .Dep_EpiDoc2      -0.261    0.565   -0.461  122.543    0.645   -1.380
    ##    .MM_EpiDoc2        0.041    0.054    0.750      Inf    0.453   -0.066
    ##    .Dep_EpiDoc1      -0.687    0.623   -1.102      Inf    0.270   -1.909
    ##    .MM_EpiDoc1       -0.825    0.186   -4.444      Inf    0.000   -1.189
    ##  ci.upper   Std.lv  Std.all
    ##     1.799    1.073    0.372
    ##     0.178   -0.352   -0.276
    ##     0.858   -0.261   -0.088
    ##     0.147    0.041    0.090
    ##     0.535   -0.687   -0.224
    ##    -0.461   -0.825   -0.863
    ## 
    ## Variances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##    .Dep_EpiDoc4       6.879    0.985    6.986      Inf    0.000    4.948
    ##    .MM_EpiDoc4        0.971    0.105    9.213      Inf    0.000    0.764
    ##    .Dep_EpiDoc2       6.201    0.861    7.202  103.268    0.000    4.493
    ##    .MM_EpiDoc2        0.187    0.051    3.676      Inf    0.000    0.087
    ##    .Dep_EpiDoc1       8.430    1.156    7.292      Inf    0.000    6.164
    ##    .MM_EpiDoc1        0.700    0.086    8.162      Inf    0.000    0.532
    ##  ci.upper   Std.lv  Std.all
    ##     8.811    6.879    0.829
    ##     1.177    0.971    0.599
    ##     7.908    6.201    0.701
    ##     0.286    0.187    0.910
    ##    10.695    8.430    0.898
    ##     0.868    0.700    0.766
    ## 
    ## R-Square:
    ##                    Estimate
    ##     Dep_EpiDoc4       0.171
    ##     MM_EpiDoc4        0.401
    ##     Dep_EpiDoc2       0.299
    ##     MM_EpiDoc2        0.090
    ##     Dep_EpiDoc1       0.102
    ##     MM_EpiDoc1        0.234
    ## 
    ## 
    ## Group 2 [Female]:
    ## 
    ## Regressions:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##   Dep_EpiDoc4 ~                                                         
    ##     Dp_EpD1           0.342    0.074    4.587      Inf    0.000    0.196
    ##     Dp_EpD2           0.220    0.074    2.962  943.735    0.003    0.074
    ##   MM_EpiDoc4 ~                                                          
    ##     MM_EpD2           0.282    0.113    2.498      Inf    0.012    0.061
    ##     MM_EpD1           0.350    0.087    4.030      Inf    0.000    0.180
    ##   Dep_EpiDoc2 ~                                                         
    ##     Dp_EpD1           0.373    0.084    4.459  844.400    0.000    0.209
    ##   MM_EpiDoc2 ~                                                          
    ##     MM_EpD1           0.366    0.096    3.803      Inf    0.000    0.177
    ##   Dep_EpiDoc4 ~                                                         
    ##     MM_EpD2           0.669    0.305    2.195      Inf    0.028    0.071
    ##     MM_EpD1          -0.320    0.216   -1.481      Inf    0.139   -0.744
    ##   MM_EpiDoc4 ~                                                          
    ##     Dp_EpD2          -0.001    0.017   -0.042  494.986    0.967   -0.034
    ##     Dp_EpD1           0.025    0.018    1.411      Inf    0.158   -0.010
    ##   Dep_EpiDoc2 ~                                                         
    ##     MM_EpD1           0.086    0.246    0.350  301.262    0.726   -0.398
    ##   MM_EpiDoc2 ~                                                          
    ##     Dp_EpD1           0.007    0.013    0.565      Inf    0.572   -0.018
    ##   Dep_EpiDoc4 ~                                                         
    ##     Age               0.000                                        0.000
    ##   Dep_EpiDoc2 ~                                                         
    ##     Age      (ad2)    0.010    0.013    0.740   88.941    0.461   -0.016
    ##   Dep_EpiDoc1 ~                                                         
    ##     Age      (ad1)    0.045    0.015    3.096      Inf    0.002    0.017
    ##   MM_EpiDoc4 ~                                                          
    ##     Age     (am4b)    0.030    0.007    4.352      Inf    0.000    0.016
    ##   MM_EpiDoc2 ~                                                          
    ##     Age     (am2b)   -0.005    0.005   -1.193      Inf    0.233   -0.014
    ##   MM_EpiDoc1 ~                                                          
    ##     Age     (am1b)    0.048    0.005    9.497      Inf    0.000    0.038
    ##   Dep_EpiDoc4 ~                                                         
    ##     Educatn           0.163    0.165    0.985      Inf    0.325   -0.161
    ##   Dep_EpiDoc2 ~                                                         
    ##     Educatn           0.435    0.210    2.067  296.399    0.040    0.021
    ##   Dep_EpiDoc1 ~                                                         
    ##     Educatn           0.808    0.203    3.987      Inf    0.000    0.411
    ##   MM_EpiDoc4 ~                                                          
    ##     Educatn           0.045    0.060    0.752      Inf    0.452   -0.072
    ##   MM_EpiDoc2 ~                                                          
    ##     Educatn           0.028    0.033    0.859      Inf    0.390   -0.036
    ##   MM_EpiDoc1 ~                                                          
    ##     Educatn           0.028    0.050    0.549      Inf    0.583   -0.071
    ##  ci.upper   Std.lv  Std.all
    ##                            
    ##     0.488    0.342    0.356
    ##     0.365    0.220    0.225
    ##                            
    ##     0.503    0.282    0.163
    ##     0.521    0.350    0.290
    ##                            
    ##     0.538    0.373    0.379
    ##                            
    ##     0.555    0.366    0.524
    ##                            
    ##     1.266    0.669    0.137
    ##     0.104   -0.320   -0.094
    ##                            
    ##     0.033   -0.001   -0.002
    ##     0.060    0.025    0.074
    ##                            
    ##     0.570    0.086    0.025
    ##                            
    ##     0.033    0.007    0.038
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.035    0.010    0.030
    ##                            
    ##     0.073    0.045    0.140
    ##                            
    ##     0.043    0.030    0.271
    ##                            
    ##     0.003   -0.005   -0.086
    ##                            
    ##     0.058    0.048    0.531
    ##                            
    ##     0.487    0.163    0.050
    ##                            
    ##     0.849    0.435    0.132
    ##                            
    ##     1.205    0.808    0.240
    ##                            
    ##     0.162    0.045    0.039
    ##                            
    ##     0.092    0.028    0.042
    ##                            
    ##     0.126    0.028    0.029
    ## 
    ## Covariances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##  .Dep_EpiDoc1 ~~                                                        
    ##    .MM_EpiDoc1        0.922    0.308    2.991      Inf    0.003    0.318
    ##  .Dep_EpiDoc2 ~~                                                        
    ##    .MM_EpiDoc2        0.000                                        0.000
    ##  .Dep_EpiDoc4 ~~                                                        
    ##    .MM_EpiDoc4        0.000                                        0.000
    ##  ci.upper   Std.lv  Std.all
    ##                            
    ##     1.525    0.922    0.293
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.000    0.000    0.000
    ## 
    ## Intercepts:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##    .Dep_EpiDoc4       0.968    0.375    2.584      Inf    0.010    0.234
    ##    .MM_EpiDoc4       -0.611    0.241   -2.530      Inf    0.011   -1.084
    ##    .Dep_EpiDoc2       0.513    0.608    0.844  147.075    0.400   -0.688
    ##    .MM_EpiDoc2        0.133    0.158    0.840      Inf    0.401   -0.177
    ##    .Dep_EpiDoc1      -0.201    0.687   -0.292      Inf    0.770   -1.547
    ##    .MM_EpiDoc1       -1.370    0.225   -6.082      Inf    0.000   -1.811
    ##  ci.upper   Std.lv  Std.all
    ##     1.703    0.968    0.270
    ##    -0.138   -0.611   -0.481
    ##     1.715    0.513    0.139
    ##     0.443    0.133    0.181
    ##     1.146   -0.201   -0.054
    ##    -0.928   -1.370   -1.305
    ## 
    ## Variances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##    .Dep_EpiDoc4       9.339    1.042    8.966      Inf    0.000    7.295
    ##    .MM_EpiDoc4        0.997    0.084   11.908      Inf    0.000    0.833
    ##    .Dep_EpiDoc2      10.759    1.346    7.994  379.009    0.000    8.112
    ##    .MM_EpiDoc2        0.400    0.061    6.568      Inf    0.000    0.281
    ##    .Dep_EpiDoc1      12.643    1.730    7.310      Inf    0.000    9.253
    ##    .MM_EpiDoc1        0.781    0.120    6.489      Inf    0.000    0.545
    ##  ci.upper   Std.lv  Std.all
    ##    11.382    9.339    0.724
    ##     1.161    0.997    0.618
    ##    13.405   10.759    0.794
    ##     0.520    0.400    0.744
    ##    16.033   12.643    0.903
    ##     1.017    0.781    0.708
    ## 
    ## R-Square:
    ##                    Estimate
    ##     Dep_EpiDoc4       0.276
    ##     MM_EpiDoc4        0.382
    ##     Dep_EpiDoc2       0.206
    ##     MM_EpiDoc2        0.256
    ##     Dep_EpiDoc1       0.097
    ##     MM_EpiDoc1        0.292

``` r
anova(resFre_g,  res_covR_g,res_covAge_g)
```

    ## Robust correction can only be applied to pooled chi-squared statistic, not F statistic. "asymptotic" was switched to TRUE.

    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.

    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.

    ## Robust correction can only be applied to pooled chi-squared statistic, not F statistic. "asymptotic" was switched to TRUE.

    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.
    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.
    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.
    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.
    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.

    ## Test statistic(s) pooled using the D4 pooling method.
    ##   Pooled statistic: "standard"
    ##   Method to robustify pooled statistic:  "yuan.bentler.mplus"
    ## 
    ##              Df  Chisq Chisq diff Df diff Pr(>Chisq)     RIV     FMI
    ## resFre_g      2 2.3128                                              
    ## res_covR_g    6 5.1723     2.4926       4    0.64597 0.34350 0.25568
    ## res_covAge_g  8 7.7007     2.1111       2    0.34799 0.35321 0.26102

difference is still significant, constraining to make the effect of age
being the same among groups leads to worsening of the fit.

``` r
source('./code/parTableToCSV.R')
res<-parTableToCSV(res_covAge_g,path=('./Results/LaggedSEM_mi20/gSEX_covAge'))
```

    ## Joining with `by = join_by(group)`

    ## [1] Group2     label      estimate_p CI         Group     
    ## <0 rows> (or 0-length row.names)

``` r
res$gPartable
```

    ##                      label2 estimate_p_Male estimate_p_Female           CI_Male
    ## 1   Dep_EpiDoc4~Dep_EpiDoc1   0.164 (0.037)         0.342 (0)      [0.01:0.318]
    ## 2   Dep_EpiDoc4~Dep_EpiDoc2   0.234 (0.006)      0.22 (0.003)     [0.068:0.399]
    ## 3     MM_EpiDoc4~MM_EpiDoc2   0.367 (0.017)     0.282 (0.012)     [0.067:0.668]
    ## 4     MM_EpiDoc4~MM_EpiDoc1        0.56 (0)          0.35 (0)     [0.386:0.735]
    ## 5   Dep_EpiDoc2~Dep_EpiDoc1       0.404 (0)         0.373 (0)     [0.267:0.541]
    ## 6     MM_EpiDoc2~MM_EpiDoc1   0.102 (0.027)         0.366 (0)     [0.012:0.193]
    ## 7    Dep_EpiDoc4~MM_EpiDoc2    0.53 (0.067)     0.669 (0.028)    [-0.038:1.098]
    ## 8    Dep_EpiDoc4~MM_EpiDoc1   0.046 (0.822)     -0.32 (0.139)    [-0.359:0.452]
    ## 9    MM_EpiDoc4~Dep_EpiDoc2   0.036 (0.229)    -0.001 (0.967)    [-0.023:0.094]
    ## 10   MM_EpiDoc4~Dep_EpiDoc1   0.005 (0.875)     0.025 (0.158)    [-0.052:0.061]
    ## 11   Dep_EpiDoc2~MM_EpiDoc1    0.23 (0.298)     0.086 (0.726)    [-0.203:0.663]
    ## 12   MM_EpiDoc2~Dep_EpiDoc1    0.021 (0.17)     0.007 (0.572)     [-0.009:0.05]
    ## 13  Dep_EpiDoc1~~MM_EpiDoc1   0.548 (0.011)     0.922 (0.003)     [0.128:0.968]
    ## 14  Dep_EpiDoc2~~MM_EpiDoc2          0 (NA)            0 (NA)             [0:0]
    ## 15  Dep_EpiDoc4~~MM_EpiDoc4          0 (NA)            0 (NA)             [0:0]
    ## 16          Dep_EpiDoc4~Age          0 (NA)            0 (NA)             [0:0]
    ## 17          Dep_EpiDoc2~Age    0.01 (0.461)      0.01 (0.461)    [-0.016:0.035]
    ## 18          Dep_EpiDoc1~Age   0.045 (0.002)     0.045 (0.002)     [0.017:0.073]
    ## 19           MM_EpiDoc4~Age   0.019 (0.011)          0.03 (0)     [0.004:0.033]
    ## 20           MM_EpiDoc2~Age   0.002 (0.499)    -0.005 (0.233)    [-0.003:0.006]
    ## 21           MM_EpiDoc1~Age       0.031 (0)         0.048 (0)     [0.022:0.041]
    ## 22    Dep_EpiDoc4~Education   0.162 (0.325)     0.163 (0.325)    [-0.161:0.485]
    ## 23    Dep_EpiDoc2~Education    0.49 (0.001)      0.435 (0.04)     [0.203:0.777]
    ## 24    Dep_EpiDoc1~Education   0.472 (0.003)         0.808 (0)     [0.161:0.782]
    ## 25     MM_EpiDoc4~Education  -0.008 (0.884)     0.045 (0.452)    [-0.118:0.102]
    ## 26     MM_EpiDoc2~Education  -0.045 (0.125)      0.028 (0.39)    [-0.102:0.012]
    ## 27     MM_EpiDoc1~Education     0.03 (0.58)     0.028 (0.583)    [-0.075:0.135]
    ## 28 Dep_EpiDoc4~~Dep_EpiDoc4       6.879 (0)         9.339 (0)     [4.948:8.811]
    ## 29   MM_EpiDoc4~~MM_EpiDoc4       0.971 (0)         0.997 (0)     [0.764:1.177]
    ## 30 Dep_EpiDoc2~~Dep_EpiDoc2       6.201 (0)        10.759 (0)     [4.493:7.908]
    ## 31   MM_EpiDoc2~~MM_EpiDoc2       0.187 (0)           0.4 (0)     [0.087:0.286]
    ## 32 Dep_EpiDoc1~~Dep_EpiDoc1        8.43 (0)        12.643 (0)    [6.164:10.695]
    ## 33   MM_EpiDoc1~~MM_EpiDoc1         0.7 (0)         0.781 (0)     [0.532:0.868]
    ## 34                 Age~~Age    204.021 (NA)      135.726 (NA) [204.021:204.021]
    ## 35           Age~~Education      6.794 (NA)         3.68 (NA)     [6.794:6.794]
    ## 36     Education~~Education      1.159 (NA)         1.24 (NA)     [1.159:1.159]
    ## 37            Dep_EpiDoc4~1   1.073 (0.004)      0.968 (0.01)     [0.347:1.799]
    ## 38             MM_EpiDoc4~1  -0.352 (0.193)    -0.611 (0.011)    [-0.882:0.178]
    ## 39            Dep_EpiDoc2~1  -0.261 (0.645)       0.513 (0.4)     [-1.38:0.858]
    ## 40             MM_EpiDoc2~1   0.041 (0.453)     0.133 (0.401)    [-0.066:0.147]
    ## 41            Dep_EpiDoc1~1   -0.687 (0.27)     -0.201 (0.77)    [-1.909:0.535]
    ## 42             MM_EpiDoc1~1      -0.825 (0)         -1.37 (0)   [-1.189:-0.461]
    ## 43                    Age~1     44.188 (NA)       42.481 (NA)   [44.188:44.188]
    ## 44              Education~1      2.823 (NA)        2.409 (NA)     [2.823:2.823]
    ##            CI_Female
    ## 1      [0.196:0.488]
    ## 2      [0.074:0.365]
    ## 3      [0.061:0.503]
    ## 4       [0.18:0.521]
    ## 5      [0.209:0.538]
    ## 6      [0.177:0.555]
    ## 7      [0.071:1.266]
    ## 8     [-0.744:0.104]
    ## 9     [-0.034:0.033]
    ## 10      [-0.01:0.06]
    ## 11     [-0.398:0.57]
    ## 12    [-0.018:0.033]
    ## 13     [0.318:1.525]
    ## 14             [0:0]
    ## 15             [0:0]
    ## 16             [0:0]
    ## 17    [-0.016:0.035]
    ## 18     [0.017:0.073]
    ## 19     [0.016:0.043]
    ## 20    [-0.014:0.003]
    ## 21     [0.038:0.058]
    ## 22    [-0.161:0.487]
    ## 23     [0.021:0.849]
    ## 24     [0.411:1.205]
    ## 25    [-0.072:0.162]
    ## 26    [-0.036:0.092]
    ## 27    [-0.071:0.126]
    ## 28    [7.295:11.382]
    ## 29     [0.833:1.161]
    ## 30    [8.112:13.405]
    ## 31      [0.281:0.52]
    ## 32    [9.253:16.033]
    ## 33     [0.545:1.017]
    ## 34 [135.726:135.726]
    ## 35       [3.68:3.68]
    ## 36       [1.24:1.24]
    ## 37     [0.234:1.703]
    ## 38   [-1.084:-0.138]
    ## 39    [-0.688:1.715]
    ## 40    [-0.177:0.443]
    ## 41    [-1.547:1.146]
    ## 42   [-1.811:-0.928]
    ## 43   [42.481:42.481]
    ## 44     [2.409:2.409]

### restriction on education

#### restrictions on NO direct effect

``` r
model_covAgeEd_R<-'
# evolução
Dep_EpiDoc4~Dep_EpiDoc1+Dep_EpiDoc2
MM_EpiDoc4~MM_EpiDoc2+MM_EpiDoc1
Dep_EpiDoc2~Dep_EpiDoc1
MM_EpiDoc2~MM_EpiDoc1 

# causalidade - c/ temporalidade/lag
Dep_EpiDoc4~MM_EpiDoc2+MM_EpiDoc1
MM_EpiDoc4~Dep_EpiDoc2+Dep_EpiDoc1
Dep_EpiDoc2~MM_EpiDoc1
MM_EpiDoc2~Dep_EpiDoc1


# corrrelação mesmo t
MM_EpiDoc1~~Dep_EpiDoc1
MM_EpiDoc2~~0*Dep_EpiDoc2
MM_EpiDoc4~~0*Dep_EpiDoc4

# causalidade - idade
Dep_EpiDoc4~0*Age
Dep_EpiDoc2~ad2*Age
Dep_EpiDoc1~ad1*Age

MM_EpiDoc4~c(am4a,am4b)*Age
MM_EpiDoc2~c(am2a,am2b)*Age
MM_EpiDoc1~c(am1a,am1b)*Age



# educação
Dep_EpiDoc4~0*Education
Dep_EpiDoc2~Education
Dep_EpiDoc1~Education
MM_EpiDoc4~0*Education
MM_EpiDoc2~Education
MM_EpiDoc1~Education
'


res_covAgeEd_g<-lavaan.mi::sem.mi(data=data.mi,
         model=model_covAgeEd_R,
         estimator = "ML", test = "yuan.bentler.mplus",
         likelihood = "wishart",
           sampling.weights='ipw',group='Sex',
         orthogonal=TRUE
         )
```

    ## Warning: lavaan->lavParTable():  
    ##    using a single label per parameter in a multiple group setting implies 
    ##    imposing equality constraints across all the groups; If this is not 
    ##    intended, either remove the label(s), or use a vector of labels (one for 
    ##    each group); See the Multiple groups section in the man page of 
    ##    model.syntax.

``` r
summary(res_covAgeEd_g,fit.meas=TRUE,standardized=TRUE,rsquare=TRUE,ci=TRUE)
```

    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.

    ## lavaan.mi object fit to 20 imputed data sets using:
    ##  - lavaan    (0.6-19)
    ##  - lavaan.mi (0.1-0)
    ## See class?lavaan.mi help page for available methods. 
    ## 
    ## Convergence information:
    ## The model converged on 20 imputed data sets.
    ## Standard errors were available for all imputations.
    ## 
    ##   Estimator                                         ML
    ##   Optimization method                           NLMINB
    ##   Number of model parameters                        68
    ##   Number of equality constraints                     2
    ## 
    ##   Number of observations per group:                   
    ##     Male                                           355
    ##     Female                                         317
    ##   Sampling weights variable                        ipw
    ## 
    ## Model Test User Model:
    ## 
    ##                                                     Standard      Scaled
    ##   Test statistic                                      10.619       9.674
    ##   Degrees of freedom                                      12          12
    ##   P-value                                              0.562       0.645
    ##   Average scaling correction factor                                1.098
    ##   Pooling method                                          D4            
    ##     Pooled statistic                              "standard"            
    ##     "yuan.bentler.mplus" correction applied            AFTER     pooling
    ## 
    ## Model Test Baseline Model:
    ## 
    ##   Test statistic                              1048.982     846.230
    ##   Degrees of freedom                                54          54
    ##   P-value                                        0.000       0.000
    ##   Scaling correction factor                                  1.240
    ## 
    ## User Model versus Baseline Model:
    ## 
    ##   Comparative Fit Index (CFI)                    1.000       1.000
    ##   Tucker-Lewis Index (TLI)                       1.006       1.013
    ##                                                                   
    ##   Robust Comparative Fit Index (CFI)                         1.000
    ##   Robust Tucker-Lewis Index (TLI)                            1.012
    ## 
    ## Loglikelihood and Information Criteria:
    ## 
    ##   Loglikelihood user model (H0)              -7330.378   -7330.378
    ##   Scaling correction factor                                  1.779
    ##       for the MLR correction                                      
    ##   Loglikelihood unrestricted model (H1)      -7325.594   -7325.594
    ##   Scaling correction factor                                  1.717
    ##       for the MLR correction                                      
    ##                                                                   
    ##   Akaike (AIC)                               14792.756   14792.756
    ##   Bayesian (BIC)                             15090.433   15090.433
    ##   Sample-size adjusted Bayesian (SABIC)      14880.878   14880.878
    ## 
    ## Root Mean Square Error of Approximation:
    ## 
    ##   RMSEA                                          0.000       0.000
    ##   90 Percent confidence interval - lower         0.000       0.000
    ##   90 Percent confidence interval - upper         0.050       0.045
    ##   P-value H_0: RMSEA <= 0.050                    0.948       0.969
    ##   P-value H_0: RMSEA >= 0.080                    0.001       0.000
    ##                                                                   
    ##   Robust RMSEA                                               0.000
    ##   90 Percent confidence interval - lower                     0.000
    ##   90 Percent confidence interval - upper                     0.048
    ##   P-value H_0: Robust RMSEA <= 0.050                         0.957
    ##   P-value H_0: Robust RMSEA >= 0.080                         0.001
    ## 
    ## Standardized Root Mean Square Residual:
    ## 
    ##   SRMR                                           0.028       0.028
    ## 
    ## Parameter Estimates:
    ## 
    ##   Standard errors                                    Sandwich
    ##   Information bread                                  Observed
    ##   Observed information based on                       Hessian
    ##                                                              
    ##   Pooled across imputations              Rubin's (1987) rules
    ##   Augment within-imputation variance     Scale by average RIV
    ##   Wald test for pooled parameters          t(df) distribution
    ## 
    ##   Pooled t statistics with df >= 1000 are displayed with
    ##   df = Inf(inity) to save space. Although the t distribution
    ##   with large df closely approximates a standard normal
    ##   distribution, exact df for reporting these t tests can be
    ##   obtained from parameterEstimates.mi() 
    ## 
    ## 
    ## 
    ## Group 1 [Male]:
    ## 
    ## Regressions:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##   Dep_EpiDoc4 ~                                                         
    ##     Dp_EpD1           0.170    0.079    2.163  563.968    0.031    0.016
    ##     Dp_EpD2           0.247    0.080    3.086  334.567    0.002    0.090
    ##   MM_EpiDoc4 ~                                                          
    ##     MM_EpD2           0.369    0.154    2.404      Inf    0.016    0.068
    ##     MM_EpD1           0.560    0.089    6.280      Inf    0.000    0.385
    ##   Dep_EpiDoc2 ~                                                         
    ##     Dp_EpD1           0.404    0.070    5.786  382.957    0.000    0.267
    ##   MM_EpiDoc2 ~                                                          
    ##     MM_EpD1           0.102    0.046    2.208      Inf    0.027    0.012
    ##   Dep_EpiDoc4 ~                                                         
    ##     MM_EpD2           0.494    0.285    1.736      Inf    0.083   -0.064
    ##     MM_EpD1           0.077    0.198    0.390      Inf    0.697   -0.311
    ##   MM_EpiDoc4 ~                                                          
    ##     Dp_EpD2           0.035    0.030    1.172  210.440    0.243   -0.024
    ##     Dp_EpD1           0.004    0.029    0.153      Inf    0.879   -0.052
    ##   Dep_EpiDoc2 ~                                                         
    ##     MM_EpD1           0.230    0.220    1.042  358.230    0.298   -0.204
    ##   MM_EpiDoc2 ~                                                          
    ##     Dp_EpD1           0.021    0.015    1.372      Inf    0.170   -0.009
    ##   Dep_EpiDoc4 ~                                                         
    ##     Age               0.000                                        0.000
    ##   Dep_EpiDoc2 ~                                                         
    ##     Age      (ad2)    0.010    0.013    0.739   88.940    0.462   -0.016
    ##   Dep_EpiDoc1 ~                                                         
    ##     Age      (ad1)    0.045    0.015    3.093      Inf    0.002    0.016
    ##   MM_EpiDoc4 ~                                                          
    ##     Age     (am4a)    0.018    0.007    2.540      Inf    0.011    0.004
    ##   MM_EpiDoc2 ~                                                          
    ##     Age     (am2a)    0.002    0.002    0.675      Inf    0.500   -0.003
    ##   MM_EpiDoc1 ~                                                          
    ##     Age     (am1a)    0.031    0.005    6.374      Inf    0.000    0.022
    ##   Dep_EpiDoc4 ~                                                         
    ##     Educatn           0.000                                        0.000
    ##   Dep_EpiDoc2 ~                                                         
    ##     Educatn           0.490    0.146    3.353  335.391    0.001    0.202
    ##   Dep_EpiDoc1 ~                                                         
    ##     Educatn           0.472    0.159    2.973      Inf    0.003    0.161
    ##   MM_EpiDoc4 ~                                                          
    ##     Educatn           0.000                                        0.000
    ##   MM_EpiDoc2 ~                                                          
    ##     Educatn          -0.045    0.029   -1.533      Inf    0.125   -0.102
    ##   MM_EpiDoc1 ~                                                          
    ##     Educatn           0.030    0.054    0.553      Inf    0.580   -0.075
    ##  ci.upper   Std.lv  Std.all
    ##                            
    ##     0.325    0.170    0.181
    ##     0.405    0.247    0.255
    ##                            
    ##     0.671    0.369    0.132
    ##     0.735    0.560    0.421
    ##                            
    ##     0.541    0.404    0.417
    ##                            
    ##     0.193    0.102    0.216
    ##                            
    ##     1.052    0.494    0.078
    ##     0.465    0.077    0.026
    ##                            
    ##     0.094    0.035    0.082
    ##     0.061    0.004    0.011
    ##                            
    ##     0.663    0.230    0.074
    ##                            
    ##     0.050    0.021    0.140
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.035    0.010    0.046
    ##                            
    ##     0.073    0.045    0.210
    ##                            
    ##     0.033    0.018    0.206
    ##                            
    ##     0.006    0.002    0.049
    ##                            
    ##     0.041    0.031    0.468
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.777    0.490    0.177
    ##                            
    ##     0.783    0.472    0.166
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.013   -0.045   -0.107
    ##                            
    ##     0.135    0.030    0.033
    ## 
    ## Covariances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##  .Dep_EpiDoc1 ~~                                                        
    ##    .MM_EpiDoc1        0.548    0.214    2.556      Inf    0.011    0.128
    ##  .Dep_EpiDoc2 ~~                                                        
    ##    .MM_EpiDoc2        0.000                                        0.000
    ##  .Dep_EpiDoc4 ~~                                                        
    ##    .MM_EpiDoc4        0.000                                        0.000
    ##  ci.upper   Std.lv  Std.all
    ##                            
    ##     0.968    0.548    0.226
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.000    0.000    0.000
    ## 
    ## Intercepts:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##    .Dep_EpiDoc4       1.461    0.212    6.875  401.448    0.000    1.043
    ##    .MM_EpiDoc4       -0.362    0.258   -1.403      Inf    0.161   -0.868
    ##    .Dep_EpiDoc2      -0.261    0.566   -0.461  122.542    0.646   -1.381
    ##    .MM_EpiDoc2        0.041    0.054    0.750      Inf    0.453   -0.066
    ##    .Dep_EpiDoc1      -0.687    0.624   -1.101      Inf    0.271   -1.910
    ##    .MM_EpiDoc1       -0.825    0.186   -4.440      Inf    0.000   -1.189
    ##  ci.upper   Std.lv  Std.all
    ##     1.878    1.461    0.507
    ##     0.144   -0.362   -0.284
    ##     0.859   -0.261   -0.088
    ##     0.147    0.041    0.090
    ##     0.536   -0.687   -0.224
    ##    -0.461   -0.825   -0.863
    ## 
    ## Variances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##    .Dep_EpiDoc4       6.910    1.001    6.903      Inf    0.000    4.947
    ##    .MM_EpiDoc4        0.971    0.106    9.201      Inf    0.000    0.764
    ##    .Dep_EpiDoc2       6.201    0.862    7.196  103.268    0.000    4.492
    ##    .MM_EpiDoc2        0.187    0.051    3.673      Inf    0.000    0.087
    ##    .Dep_EpiDoc1       8.430    1.157    7.286      Inf    0.000    6.162
    ##    .MM_EpiDoc1        0.700    0.086    8.155      Inf    0.000    0.532
    ##  ci.upper   Std.lv  Std.all
    ##     8.874    6.910    0.832
    ##     1.178    0.971    0.599
    ##     7.910    6.201    0.701
    ##     0.286    0.187    0.910
    ##    10.697    8.430    0.898
    ##     0.869    0.700    0.766
    ## 
    ## R-Square:
    ##                    Estimate
    ##     Dep_EpiDoc4       0.168
    ##     MM_EpiDoc4        0.401
    ##     Dep_EpiDoc2       0.299
    ##     MM_EpiDoc2        0.090
    ##     Dep_EpiDoc1       0.102
    ##     MM_EpiDoc1        0.234
    ## 
    ## 
    ## Group 2 [Female]:
    ## 
    ## Regressions:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##   Dep_EpiDoc4 ~                                                         
    ##     Dp_EpD1           0.351    0.074    4.721      Inf    0.000    0.205
    ##     Dp_EpD2           0.227    0.075    3.031      Inf    0.002    0.080
    ##   MM_EpiDoc4 ~                                                          
    ##     MM_EpD2           0.285    0.112    2.530      Inf    0.011    0.064
    ##     MM_EpD1           0.347    0.087    3.977      Inf    0.000    0.176
    ##   Dep_EpiDoc2 ~                                                         
    ##     Dp_EpD1           0.373    0.084    4.455  844.406    0.000    0.209
    ##   MM_EpiDoc2 ~                                                          
    ##     MM_EpD1           0.366    0.096    3.800      Inf    0.000    0.177
    ##   Dep_EpiDoc4 ~                                                         
    ##     MM_EpD2           0.675    0.306    2.205      Inf    0.028    0.075
    ##     MM_EpD1          -0.309    0.215   -1.437      Inf    0.151   -0.730
    ##   MM_EpiDoc4 ~                                                          
    ##     Dp_EpD2           0.001    0.017    0.057  492.345    0.954   -0.032
    ##     Dp_EpD1           0.028    0.018    1.565      Inf    0.118   -0.007
    ##   Dep_EpiDoc2 ~                                                         
    ##     MM_EpD1           0.086    0.246    0.350  301.271    0.727   -0.398
    ##   MM_EpiDoc2 ~                                                          
    ##     Dp_EpD1           0.007    0.013    0.565      Inf    0.572   -0.018
    ##   Dep_EpiDoc4 ~                                                         
    ##     Age               0.000                                        0.000
    ##   Dep_EpiDoc2 ~                                                         
    ##     Age      (ad2)    0.010    0.013    0.739   88.940    0.462   -0.016
    ##   Dep_EpiDoc1 ~                                                         
    ##     Age      (ad1)    0.045    0.015    3.093      Inf    0.002    0.016
    ##   MM_EpiDoc4 ~                                                          
    ##     Age     (am4b)    0.031    0.007    4.513      Inf    0.000    0.017
    ##   MM_EpiDoc2 ~                                                          
    ##     Age     (am2b)   -0.005    0.005   -1.192      Inf    0.233   -0.014
    ##   MM_EpiDoc1 ~                                                          
    ##     Age     (am1b)    0.048    0.005    9.489      Inf    0.000    0.038
    ##   Dep_EpiDoc4 ~                                                         
    ##     Educatn           0.000                                        0.000
    ##   Dep_EpiDoc2 ~                                                         
    ##     Educatn           0.435    0.211    2.066  296.395    0.040    0.021
    ##   Dep_EpiDoc1 ~                                                         
    ##     Educatn           0.808    0.203    3.984      Inf    0.000    0.410
    ##   MM_EpiDoc4 ~                                                          
    ##     Educatn           0.000                                        0.000
    ##   MM_EpiDoc2 ~                                                          
    ##     Educatn           0.028    0.033    0.858      Inf    0.391   -0.036
    ##   MM_EpiDoc1 ~                                                          
    ##     Educatn           0.028    0.050    0.549      Inf    0.583   -0.071
    ##  ci.upper   Std.lv  Std.all
    ##                            
    ##     0.496    0.351    0.365
    ##     0.374    0.227    0.233
    ##                            
    ##     0.505    0.285    0.164
    ##     0.519    0.347    0.287
    ##                            
    ##     0.538    0.373    0.379
    ##                            
    ##     0.555    0.366    0.524
    ##                            
    ##     1.275    0.675    0.138
    ##     0.113   -0.309   -0.090
    ##                            
    ##     0.034    0.001    0.003
    ##     0.062    0.028    0.081
    ##                            
    ##     0.570    0.086    0.025
    ##                            
    ##     0.033    0.007    0.038
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.035    0.010    0.030
    ##                            
    ##     0.073    0.045    0.140
    ##                            
    ##     0.044    0.031    0.280
    ##                            
    ##     0.003   -0.005   -0.086
    ##                            
    ##     0.058    0.048    0.531
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.849    0.435    0.132
    ##                            
    ##     1.205    0.808    0.240
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.092    0.028    0.042
    ##                            
    ##     0.126    0.028    0.029
    ## 
    ## Covariances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##  .Dep_EpiDoc1 ~~                                                        
    ##    .MM_EpiDoc1        0.922    0.308    2.989      Inf    0.003    0.317
    ##  .Dep_EpiDoc2 ~~                                                        
    ##    .MM_EpiDoc2        0.000                                        0.000
    ##  .Dep_EpiDoc4 ~~                                                        
    ##    .MM_EpiDoc4        0.000                                        0.000
    ##  ci.upper   Std.lv  Std.all
    ##                            
    ##     1.526    0.922    0.293
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.000    0.000    0.000
    ## 
    ## Intercepts:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##    .Dep_EpiDoc4       1.293    0.246    5.257      Inf    0.000    0.811
    ##    .MM_EpiDoc4       -0.559    0.230   -2.431      Inf    0.015   -1.009
    ##    .Dep_EpiDoc2       0.513    0.608    0.843  147.075    0.400   -0.689
    ##    .MM_EpiDoc2        0.133    0.158    0.839      Inf    0.401   -0.177
    ##    .Dep_EpiDoc1      -0.201    0.688   -0.292      Inf    0.770   -1.549
    ##    .MM_EpiDoc1       -1.370    0.225   -6.077      Inf    0.000   -1.811
    ##  ci.upper   Std.lv  Std.all
    ##     1.775    1.293    0.360
    ##    -0.108   -0.559   -0.440
    ##     1.716    0.513    0.139
    ##     0.443    0.133    0.181
    ##     1.147   -0.201   -0.054
    ##    -0.928   -1.370   -1.305
    ## 
    ## Variances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##    .Dep_EpiDoc4       9.370    1.053    8.902      Inf    0.000    7.305
    ##    .MM_EpiDoc4        0.999    0.085   11.791      Inf    0.000    0.833
    ##    .Dep_EpiDoc2      10.759    1.347    7.987  379.009    0.000    8.110
    ##    .MM_EpiDoc2        0.400    0.061    6.562      Inf    0.000    0.281
    ##    .Dep_EpiDoc1      12.643    1.731    7.304      Inf    0.000    9.250
    ##    .MM_EpiDoc1        0.781    0.120    6.483      Inf    0.000    0.545
    ##  ci.upper   Std.lv  Std.all
    ##    11.435    9.370    0.727
    ##     1.165    0.999    0.619
    ##    13.407   10.759    0.794
    ##     0.520    0.400    0.744
    ##    16.036   12.643    0.903
    ##     1.017    0.781    0.708
    ## 
    ## R-Square:
    ##                    Estimate
    ##     Dep_EpiDoc4       0.273
    ##     MM_EpiDoc4        0.381
    ##     Dep_EpiDoc2       0.206
    ##     MM_EpiDoc2        0.256
    ##     Dep_EpiDoc1       0.097
    ##     MM_EpiDoc1        0.292

``` r
anova(resFre_g,  res_covR_g,res_covAge_g,res_covAgeEd_g)
```

    ## Robust correction can only be applied to pooled chi-squared statistic, not F statistic. "asymptotic" was switched to TRUE.

    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.

    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.

    ## Robust correction can only be applied to pooled chi-squared statistic, not F statistic. "asymptotic" was switched to TRUE.

    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.
    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.

    ## Robust correction can only be applied to pooled chi-squared statistic, not F statistic. "asymptotic" was switched to TRUE.

    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.
    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.
    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.
    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.
    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.
    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.

    ## Test statistic(s) pooled using the D4 pooling method.
    ##   Pooled statistic: "standard"
    ##   Method to robustify pooled statistic:  "yuan.bentler.mplus"
    ## 
    ##                Df   Chisq Chisq diff Df diff Pr(>Chisq)     RIV      FMI
    ## resFre_g        2  2.3128                                               
    ## res_covR_g      6  5.1723     2.4926       4    0.64597 0.34350 0.255676
    ## res_covAge_g    8  7.7007     2.1111       2    0.34799 0.35321 0.261018
    ## res_covAgeEd_g 12 10.6192     2.7331       4    0.60343 0.08429 0.077737

difference is still significant, constraining to make the effect of age
being the same among groups leads to worsening of the fit.

``` r
source('./code/parTableToCSV.R')
res<-parTableToCSV(res_covAgeEd_g,path=('./Results/LaggedSEM_mi20/gSEX_covAgeEd'))
```

    ## Joining with `by = join_by(group)`

    ## [1] Group2     label      estimate_p CI         Group     
    ## <0 rows> (or 0-length row.names)

``` r
res$gPartable
```

    ##                      label2 estimate_p_Male estimate_p_Female           CI_Male
    ## 1   Dep_EpiDoc4~Dep_EpiDoc1    0.17 (0.031)         0.351 (0)     [0.016:0.325]
    ## 2   Dep_EpiDoc4~Dep_EpiDoc2   0.247 (0.002)     0.227 (0.002)      [0.09:0.405]
    ## 3     MM_EpiDoc4~MM_EpiDoc2   0.369 (0.016)     0.285 (0.011)     [0.068:0.671]
    ## 4     MM_EpiDoc4~MM_EpiDoc1        0.56 (0)         0.347 (0)     [0.385:0.735]
    ## 5   Dep_EpiDoc2~Dep_EpiDoc1       0.404 (0)         0.373 (0)     [0.267:0.541]
    ## 6     MM_EpiDoc2~MM_EpiDoc1   0.102 (0.027)         0.366 (0)     [0.012:0.193]
    ## 7    Dep_EpiDoc4~MM_EpiDoc2   0.494 (0.083)     0.675 (0.028)    [-0.064:1.052]
    ## 8    Dep_EpiDoc4~MM_EpiDoc1   0.077 (0.697)    -0.309 (0.151)    [-0.311:0.465]
    ## 9    MM_EpiDoc4~Dep_EpiDoc2   0.035 (0.243)     0.001 (0.954)    [-0.024:0.094]
    ## 10   MM_EpiDoc4~Dep_EpiDoc1   0.004 (0.879)     0.028 (0.118)    [-0.052:0.061]
    ## 11   Dep_EpiDoc2~MM_EpiDoc1    0.23 (0.298)     0.086 (0.727)    [-0.204:0.663]
    ## 12   MM_EpiDoc2~Dep_EpiDoc1    0.021 (0.17)     0.007 (0.572)     [-0.009:0.05]
    ## 13  Dep_EpiDoc1~~MM_EpiDoc1   0.548 (0.011)     0.922 (0.003)     [0.128:0.968]
    ## 14  Dep_EpiDoc2~~MM_EpiDoc2          0 (NA)            0 (NA)             [0:0]
    ## 15  Dep_EpiDoc4~~MM_EpiDoc4          0 (NA)            0 (NA)             [0:0]
    ## 16          Dep_EpiDoc4~Age          0 (NA)            0 (NA)             [0:0]
    ## 17          Dep_EpiDoc2~Age    0.01 (0.462)      0.01 (0.462)    [-0.016:0.035]
    ## 18          Dep_EpiDoc1~Age   0.045 (0.002)     0.045 (0.002)     [0.016:0.073]
    ## 19           MM_EpiDoc4~Age   0.018 (0.011)         0.031 (0)     [0.004:0.033]
    ## 20           MM_EpiDoc2~Age     0.002 (0.5)    -0.005 (0.233)    [-0.003:0.006]
    ## 21           MM_EpiDoc1~Age       0.031 (0)         0.048 (0)     [0.022:0.041]
    ## 22    Dep_EpiDoc4~Education          0 (NA)            0 (NA)             [0:0]
    ## 23    Dep_EpiDoc2~Education    0.49 (0.001)      0.435 (0.04)     [0.202:0.777]
    ## 24    Dep_EpiDoc1~Education   0.472 (0.003)         0.808 (0)     [0.161:0.783]
    ## 25     MM_EpiDoc4~Education          0 (NA)            0 (NA)             [0:0]
    ## 26     MM_EpiDoc2~Education  -0.045 (0.125)     0.028 (0.391)    [-0.102:0.013]
    ## 27     MM_EpiDoc1~Education     0.03 (0.58)     0.028 (0.583)    [-0.075:0.135]
    ## 28 Dep_EpiDoc4~~Dep_EpiDoc4        6.91 (0)          9.37 (0)     [4.947:8.874]
    ## 29   MM_EpiDoc4~~MM_EpiDoc4       0.971 (0)         0.999 (0)     [0.764:1.178]
    ## 30 Dep_EpiDoc2~~Dep_EpiDoc2       6.201 (0)        10.759 (0)      [4.492:7.91]
    ## 31   MM_EpiDoc2~~MM_EpiDoc2       0.187 (0)           0.4 (0)     [0.087:0.286]
    ## 32 Dep_EpiDoc1~~Dep_EpiDoc1        8.43 (0)        12.643 (0)    [6.162:10.697]
    ## 33   MM_EpiDoc1~~MM_EpiDoc1         0.7 (0)         0.781 (0)     [0.532:0.869]
    ## 34                 Age~~Age    204.021 (NA)      135.726 (NA) [204.021:204.021]
    ## 35           Age~~Education      6.794 (NA)         3.68 (NA)     [6.794:6.794]
    ## 36     Education~~Education      1.159 (NA)         1.24 (NA)     [1.159:1.159]
    ## 37            Dep_EpiDoc4~1       1.461 (0)         1.293 (0)     [1.043:1.878]
    ## 38             MM_EpiDoc4~1  -0.362 (0.161)    -0.559 (0.015)    [-0.868:0.144]
    ## 39            Dep_EpiDoc2~1  -0.261 (0.646)       0.513 (0.4)    [-1.381:0.859]
    ## 40             MM_EpiDoc2~1   0.041 (0.453)     0.133 (0.401)    [-0.066:0.147]
    ## 41            Dep_EpiDoc1~1  -0.687 (0.271)     -0.201 (0.77)     [-1.91:0.536]
    ## 42             MM_EpiDoc1~1      -0.825 (0)         -1.37 (0)   [-1.189:-0.461]
    ## 43                    Age~1     44.188 (NA)       42.481 (NA)   [44.188:44.188]
    ## 44              Education~1      2.823 (NA)        2.409 (NA)     [2.823:2.823]
    ##            CI_Female
    ## 1      [0.205:0.496]
    ## 2       [0.08:0.374]
    ## 3      [0.064:0.505]
    ## 4      [0.176:0.519]
    ## 5      [0.209:0.538]
    ## 6      [0.177:0.555]
    ## 7      [0.075:1.275]
    ## 8      [-0.73:0.113]
    ## 9     [-0.032:0.034]
    ## 10    [-0.007:0.062]
    ## 11     [-0.398:0.57]
    ## 12    [-0.018:0.033]
    ## 13     [0.317:1.526]
    ## 14             [0:0]
    ## 15             [0:0]
    ## 16             [0:0]
    ## 17    [-0.016:0.035]
    ## 18     [0.016:0.073]
    ## 19     [0.017:0.044]
    ## 20    [-0.014:0.003]
    ## 21     [0.038:0.058]
    ## 22             [0:0]
    ## 23     [0.021:0.849]
    ## 24      [0.41:1.205]
    ## 25             [0:0]
    ## 26    [-0.036:0.092]
    ## 27    [-0.071:0.126]
    ## 28    [7.305:11.435]
    ## 29     [0.833:1.165]
    ## 30     [8.11:13.407]
    ## 31      [0.281:0.52]
    ## 32     [9.25:16.036]
    ## 33     [0.545:1.017]
    ## 34 [135.726:135.726]
    ## 35       [3.68:3.68]
    ## 36       [1.24:1.24]
    ## 37     [0.811:1.775]
    ## 38   [-1.009:-0.108]
    ## 39    [-0.689:1.716]
    ## 40    [-0.177:0.443]
    ## 41    [-1.549:1.147]
    ## 42   [-1.811:-0.928]
    ## 43   [42.481:42.481]
    ## 44     [2.409:2.409]

#### restrictions on NO effect MM

``` r
model_covAgeEd2_R<-'
# evolução
Dep_EpiDoc4~c(dd14a,dd14b)*Dep_EpiDoc1+c(dd24a,dd24b)*Dep_EpiDoc2
Dep_EpiDoc2~c(dd12a,dd12b)*Dep_EpiDoc1

MM_EpiDoc4~c(mm24a,mm24b)*MM_EpiDoc2+c(mm14a,mm14b)*MM_EpiDoc1
MM_EpiDoc2~c(mm12a,mm12b)*MM_EpiDoc1 

# causalidade - c/ temporalidade/lag
Dep_EpiDoc4~c(md24a,md24b)*MM_EpiDoc2+c(md14a,md14b)*MM_EpiDoc1
Dep_EpiDoc2~c(md12a,md12b)*MM_EpiDoc1

MM_EpiDoc4~c(dm24a,dm24b)*Dep_EpiDoc2+c(dm14a,dm14b)*Dep_EpiDoc1
MM_EpiDoc2~c(dm12a,dm12b)*Dep_EpiDoc1


# corrrelação mesmo t
MM_EpiDoc1~~Dep_EpiDoc1
MM_EpiDoc2~~0*Dep_EpiDoc2
MM_EpiDoc4~~0*Dep_EpiDoc4

# causalidade - idade
Dep_EpiDoc4~0*Age
Dep_EpiDoc2~ad2*Age
Dep_EpiDoc1~ad1*Age

MM_EpiDoc4~c(am4a,am4b)*Age
MM_EpiDoc2~c(am2a,am2b)*Age
MM_EpiDoc1~c(am1a,am1b)*Age



# educação
Dep_EpiDoc4~0*Education
Dep_EpiDoc2~c(ed2a,ed2b)*Education
Dep_EpiDoc1~c(ed1a,ed1b)*Education
MM_EpiDoc4~0*Education
MM_EpiDoc2~0*Education
MM_EpiDoc1~0*Education


# effect on multimorbidity


# effects on depression =============
##indirect effect age on dep4
ada:=ad1*dd14a+ad2*dd24a+ad1*dd12a*dd24a
adb:=ad1*dd14b+ad2*dd24b+ad1*dd12b*dd24b
amda:=am1a*md14a+am2a*md24a+am1a*mm12a*md24a+ad1*dm12a*md24a
amdb:=am1b*md14b+am2b*md24b+am1b*mm12b*md24b+ad1*dm12b*md24b


adta:=ada+amda
adtb:=adb+amdb

##indirect effect MM1 on dep4
mda:=md12a*dd24a+mm12a*md24a+md14a
mdb:=md12b*dd24b+mm12b*md24b+md14b

## cumulative effect d1 on dep4
dda:=dm12a*md24a+dd12a*dd24a+dd14a
ddb:=dm12b*md24b+dd12b*dd24b+dd14b

## indirect effect of educaiton on dep4
eda:=ed1a*dd12a*dd24a+ed2a*dd24a
edb:=ed1b*dd12b*dd24b+ed2b*dd24b
emda:=ed1a*dm12a*md24a
emdb:=ed1b*dm12b*md24b

edta:=eda+amda
edtb:=edb+emdb

'


res_covAgeEd2_g<-lavaan.mi::sem.mi(data=data.mi,
         model=model_covAgeEd2_R,
         estimator = "ML", test = "yuan.bentler.mplus",
         likelihood = "wishart",
           sampling.weights='ipw',group='Sex',
         orthogonal=TRUE
         )
```

    ## Warning: lavaan->lavParTable():  
    ##    using a single label per parameter in a multiple group setting implies 
    ##    imposing equality constraints across all the groups; If this is not 
    ##    intended, either remove the label(s), or use a vector of labels (one for 
    ##    each group); See the Multiple groups section in the man page of 
    ##    model.syntax.

``` r
summary(res_covAgeEd2_g,fit.meas=TRUE,standardized=TRUE,rsquare=TRUE,ci=TRUE)
```

    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.

    ## lavaan.mi object fit to 20 imputed data sets using:
    ##  - lavaan    (0.6-19)
    ##  - lavaan.mi (0.1-0)
    ## See class?lavaan.mi help page for available methods. 
    ## 
    ## Convergence information:
    ## The model converged on 20 imputed data sets.
    ## Standard errors were available for all imputations.
    ## 
    ##   Estimator                                         ML
    ##   Optimization method                           NLMINB
    ##   Number of model parameters                        64
    ##   Number of equality constraints                     2
    ## 
    ##   Number of observations per group:                   
    ##     Male                                           355
    ##     Female                                         317
    ##   Sampling weights variable                        ipw
    ## 
    ## Model Test User Model:
    ## 
    ##                                                     Standard      Scaled
    ##   Test statistic                                      15.467      14.210
    ##   Degrees of freedom                                      16          16
    ##   P-value                                              0.491       0.583
    ##   Average scaling correction factor                                1.088
    ##   Pooling method                                          D4            
    ##     Pooled statistic                              "standard"            
    ##     "yuan.bentler.mplus" correction applied            AFTER     pooling
    ## 
    ## Model Test Baseline Model:
    ## 
    ##   Test statistic                              1048.982     846.230
    ##   Degrees of freedom                                54          54
    ##   P-value                                        0.000       0.000
    ##   Scaling correction factor                                  1.240
    ## 
    ## User Model versus Baseline Model:
    ## 
    ##   Comparative Fit Index (CFI)                    1.000       1.000
    ##   Tucker-Lewis Index (TLI)                       1.002       1.008
    ##                                                                   
    ##   Robust Comparative Fit Index (CFI)                         1.000
    ##   Robust Tucker-Lewis Index (TLI)                            1.007
    ## 
    ## Loglikelihood and Information Criteria:
    ## 
    ##   Loglikelihood user model (H0)              -7333.788   -7333.788
    ##   Scaling correction factor                                  1.824
    ##       for the MLR correction                                      
    ##   Loglikelihood unrestricted model (H1)      -7325.594   -7325.594
    ##   Scaling correction factor                                  1.717
    ##       for the MLR correction                                      
    ##                                                                   
    ##   Akaike (AIC)                               14791.577   14791.577
    ##   Bayesian (BIC)                             15071.213   15071.213
    ##   Sample-size adjusted Bayesian (SABIC)      14874.358   14874.358
    ## 
    ## Root Mean Square Error of Approximation:
    ## 
    ##   RMSEA                                          0.000       0.000
    ##   90 Percent confidence interval - lower         0.000       0.000
    ##   90 Percent confidence interval - upper         0.049       0.044
    ##   P-value H_0: RMSEA <= 0.050                    0.956       0.975
    ##   P-value H_0: RMSEA >= 0.080                    0.000       0.000
    ##                                                                   
    ##   Robust RMSEA                                               0.000
    ##   90 Percent confidence interval - lower                     0.000
    ##   90 Percent confidence interval - upper                     0.047
    ##   P-value H_0: Robust RMSEA <= 0.050                         0.964
    ##   P-value H_0: Robust RMSEA >= 0.080                         0.000
    ## 
    ## Standardized Root Mean Square Residual:
    ## 
    ##   SRMR                                           0.030       0.030
    ## 
    ## Parameter Estimates:
    ## 
    ##   Standard errors                                    Sandwich
    ##   Information bread                                  Observed
    ##   Observed information based on                       Hessian
    ##                                                              
    ##   Pooled across imputations              Rubin's (1987) rules
    ##   Augment within-imputation variance     Scale by average RIV
    ##   Wald test for pooled parameters          t(df) distribution
    ## 
    ##   Pooled t statistics with df >= 1000 are displayed with
    ##   df = Inf(inity) to save space. Although the t distribution
    ##   with large df closely approximates a standard normal
    ##   distribution, exact df for reporting these t tests can be
    ##   obtained from parameterEstimates.mi() 
    ## 
    ## 
    ## 
    ## Group 1 [Male]:
    ## 
    ## Regressions:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##   Dep_EpiDoc4 ~                                                         
    ##     Dp_EpD1 (dd14)    0.170    0.079    2.158  563.980    0.031    0.015
    ##     Dp_EpD2 (dd24)    0.247    0.080    3.079  334.572    0.002    0.089
    ##   Dep_EpiDoc2 ~                                                         
    ##     Dp_EpD1 (dd12)    0.404    0.070    5.772  382.956    0.000    0.266
    ##   MM_EpiDoc4 ~                                                          
    ##     MM_EpD2 (mm24)    0.369    0.154    2.398      Inf    0.016    0.068
    ##     MM_EpD1 (mm14)    0.560    0.089    6.265      Inf    0.000    0.385
    ##   MM_EpiDoc2 ~                                                          
    ##     MM_EpD1 (mm12)    0.102    0.047    2.183      Inf    0.029    0.010
    ##   Dep_EpiDoc4 ~                                                         
    ##     MM_EpD2 (md24)    0.494    0.285    1.732      Inf    0.083   -0.065
    ##     MM_EpD1 (md14)    0.077    0.198    0.389      Inf    0.697   -0.312
    ##   Dep_EpiDoc2 ~                                                         
    ##     MM_EpD1 (md12)    0.230    0.221    1.039  358.228    0.299   -0.205
    ##   MM_EpiDoc4 ~                                                          
    ##     Dp_EpD2 (dm24)    0.035    0.030    1.169  210.440    0.244   -0.024
    ##     Dp_EpD1 (dm14)    0.004    0.029    0.152      Inf    0.879   -0.052
    ##   MM_EpiDoc2 ~                                                          
    ##     Dp_EpD1 (dm12)    0.018    0.015    1.206      Inf    0.228   -0.011
    ##   Dep_EpiDoc4 ~                                                         
    ##     Age               0.000                                        0.000
    ##   Dep_EpiDoc2 ~                                                         
    ##     Age      (ad2)    0.010    0.013    0.737   88.940    0.463   -0.016
    ##   Dep_EpiDoc1 ~                                                         
    ##     Age      (ad1)    0.046    0.014    3.165      Inf    0.002    0.017
    ##   MM_EpiDoc4 ~                                                          
    ##     Age     (am4a)    0.018    0.007    2.534      Inf    0.011    0.004
    ##   MM_EpiDoc2 ~                                                          
    ##     Age     (am2a)    0.000    0.002    0.108      Inf    0.914   -0.003
    ##   MM_EpiDoc1 ~                                                          
    ##     Age     (am1a)    0.032    0.005    7.009      Inf    0.000    0.023
    ##   Dep_EpiDoc4 ~                                                         
    ##     Educatn           0.000                                        0.000
    ##   Dep_EpiDoc2 ~                                                         
    ##     Educatn (ed2a)    0.490    0.146    3.345  335.383    0.001    0.202
    ##   Dep_EpiDoc1 ~                                                         
    ##     Educatn (ed1a)    0.448    0.146    3.061      Inf    0.002    0.161
    ##   MM_EpiDoc4 ~                                                          
    ##     Educatn           0.000                                        0.000
    ##   MM_EpiDoc2 ~                                                          
    ##     Educatn           0.000                                        0.000
    ##   MM_EpiDoc1 ~                                                          
    ##     Educatn           0.000                                        0.000
    ##  ci.upper   Std.lv  Std.all
    ##                            
    ##     0.325    0.170    0.181
    ##     0.405    0.247    0.255
    ##                            
    ##     0.542    0.404    0.417
    ##                            
    ##     0.671    0.369    0.131
    ##     0.736    0.560    0.421
    ##                            
    ##     0.194    0.102    0.216
    ##                            
    ##     1.053    0.494    0.078
    ##     0.466    0.077    0.026
    ##                            
    ##     0.664    0.230    0.074
    ##                            
    ##     0.094    0.035    0.082
    ##     0.061    0.004    0.011
    ##                            
    ##     0.048    0.018    0.123
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.035    0.010    0.046
    ##                            
    ##     0.074    0.046    0.214
    ##                            
    ##     0.033    0.018    0.206
    ##                            
    ##     0.004    0.000    0.006
    ##                            
    ##     0.041    0.032    0.483
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.777    0.490    0.177
    ##                            
    ##     0.735    0.448    0.158
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.000    0.000    0.000
    ## 
    ## Covariances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##  .Dep_EpiDoc1 ~~                                                        
    ##    .MM_EpiDoc1        0.549    0.216    2.543      Inf    0.011    0.126
    ##  .Dep_EpiDoc2 ~~                                                        
    ##    .MM_EpiDoc2        0.000                                        0.000
    ##  .Dep_EpiDoc4 ~~                                                        
    ##    .MM_EpiDoc4        0.000                                        0.000
    ##  ci.upper   Std.lv  Std.all
    ##                            
    ##     0.973    0.549    0.226
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.000    0.000    0.000
    ## 
    ## Intercepts:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##    .Dep_EpiDoc4       1.461    0.213    6.859  401.449    0.000    1.042
    ##    .Dep_EpiDoc2      -0.261    0.567   -0.460  122.542    0.646   -1.383
    ##    .MM_EpiDoc4       -0.362    0.259   -1.400      Inf    0.162   -0.869
    ##    .MM_EpiDoc2       -0.020    0.063   -0.311      Inf    0.756   -0.144
    ##    .Dep_EpiDoc1      -0.657    0.620   -1.059      Inf    0.290   -1.872
    ##    .MM_EpiDoc1       -0.785    0.176   -4.469      Inf    0.000   -1.129
    ##  ci.upper   Std.lv  Std.all
    ##     1.879    1.461    0.507
    ##     0.862   -0.261   -0.088
    ##     0.145   -0.362   -0.285
    ##     0.105   -0.020   -0.044
    ##     0.559   -0.657   -0.214
    ##    -0.441   -0.785   -0.821
    ## 
    ## Variances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##    .Dep_EpiDoc4       6.910    1.003    6.887      Inf    0.000    4.942
    ##    .Dep_EpiDoc2       6.201    0.864    7.178  103.268    0.000    4.487
    ##    .MM_EpiDoc4        0.971    0.106    9.179      Inf    0.000    0.764
    ##    .MM_EpiDoc2        0.189    0.052    3.638      Inf    0.000    0.087
    ##    .Dep_EpiDoc1       8.431    1.161    7.264      Inf    0.000    6.156
    ##    .MM_EpiDoc1        0.701    0.087    8.072      Inf    0.000    0.531
    ##  ci.upper   Std.lv  Std.all
    ##     8.878    6.910    0.832
    ##     7.914    6.201    0.703
    ##     1.178    0.971    0.599
    ##     0.290    0.189    0.919
    ##    10.706    8.431    0.900
    ##     0.871    0.701    0.767
    ## 
    ## R-Square:
    ##                    Estimate
    ##     Dep_EpiDoc4       0.168
    ##     Dep_EpiDoc2       0.297
    ##     MM_EpiDoc4        0.401
    ##     MM_EpiDoc2        0.081
    ##     Dep_EpiDoc1       0.100
    ##     MM_EpiDoc1        0.233
    ## 
    ## 
    ## Group 2 [Female]:
    ## 
    ## Regressions:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##   Dep_EpiDoc4 ~                                                         
    ##     Dp_EpD1 (dd14)    0.351    0.074    4.710      Inf    0.000    0.205
    ##     Dp_EpD2 (dd24)    0.227    0.075    3.024      Inf    0.003    0.080
    ##   Dep_EpiDoc2 ~                                                         
    ##     Dp_EpD1 (dd12)    0.373    0.084    4.444  844.392    0.000    0.208
    ##   MM_EpiDoc4 ~                                                          
    ##     MM_EpD2 (mm24)    0.285    0.113    2.524      Inf    0.012    0.064
    ##     MM_EpD1 (mm14)    0.347    0.088    3.968      Inf    0.000    0.176
    ##   MM_EpiDoc2 ~                                                          
    ##     MM_EpD1 (mm12)    0.365    0.096    3.803      Inf    0.000    0.177
    ##   Dep_EpiDoc4 ~                                                         
    ##     MM_EpD2 (md24)    0.675    0.307    2.199      Inf    0.028    0.073
    ##     MM_EpD1 (md14)   -0.309    0.215   -1.433      Inf    0.152   -0.731
    ##   Dep_EpiDoc2 ~                                                         
    ##     MM_EpD1 (md12)    0.086    0.247    0.349  301.269    0.727   -0.399
    ##   MM_EpiDoc4 ~                                                          
    ##     Dp_EpD2 (dm24)    0.001    0.017    0.057  492.345    0.955   -0.032
    ##     Dp_EpD1 (dm14)    0.028    0.018    1.562      Inf    0.118   -0.007
    ##   MM_EpiDoc2 ~                                                          
    ##     Dp_EpD1 (dm12)    0.009    0.013    0.739      Inf    0.460   -0.015
    ##   Dep_EpiDoc4 ~                                                         
    ##     Age               0.000                                        0.000
    ##   Dep_EpiDoc2 ~                                                         
    ##     Age      (ad2)    0.010    0.013    0.737   88.940    0.463   -0.016
    ##   Dep_EpiDoc1 ~                                                         
    ##     Age      (ad1)    0.046    0.014    3.165      Inf    0.002    0.017
    ##   MM_EpiDoc4 ~                                                          
    ##     Age     (am4b)    0.031    0.007    4.502      Inf    0.000    0.017
    ##   MM_EpiDoc2 ~                                                          
    ##     Age     (am2b)   -0.005    0.004   -1.059      Inf    0.290   -0.013
    ##   MM_EpiDoc1 ~                                                          
    ##     Age     (am1b)    0.049    0.005    9.116      Inf    0.000    0.038
    ##   Dep_EpiDoc4 ~                                                         
    ##     Educatn           0.000                                        0.000
    ##   Dep_EpiDoc2 ~                                                         
    ##     Educatn (ed2b)    0.435    0.211    2.061  296.395    0.040    0.020
    ##   Dep_EpiDoc1 ~                                                         
    ##     Educatn (ed1b)    0.776    0.197    3.935      Inf    0.000    0.389
    ##   MM_EpiDoc4 ~                                                          
    ##     Educatn           0.000                                        0.000
    ##   MM_EpiDoc2 ~                                                          
    ##     Educatn           0.000                                        0.000
    ##   MM_EpiDoc1 ~                                                          
    ##     Educatn           0.000                                        0.000
    ##  ci.upper   Std.lv  Std.all
    ##                            
    ##     0.496    0.351    0.365
    ##     0.374    0.227    0.233
    ##                            
    ##     0.538    0.373    0.379
    ##                            
    ##     0.506    0.285    0.164
    ##     0.519    0.347    0.287
    ##                            
    ##     0.553    0.365    0.522
    ##                            
    ##     1.277    0.675    0.138
    ##     0.114   -0.309   -0.090
    ##                            
    ##     0.572    0.086    0.025
    ##                            
    ##     0.034    0.001    0.003
    ##     0.062    0.028    0.081
    ##                            
    ##     0.034    0.009    0.047
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.035    0.010    0.030
    ##                            
    ##     0.074    0.046    0.143
    ##                            
    ##     0.044    0.031    0.280
    ##                            
    ##     0.004   -0.005   -0.075
    ##                            
    ##     0.059    0.049    0.539
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.850    0.435    0.132
    ##                            
    ##     1.162    0.776    0.231
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.000    0.000    0.000
    ## 
    ## Covariances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##  .Dep_EpiDoc1 ~~                                                        
    ##    .MM_EpiDoc1        0.922    0.309    2.979      Inf    0.003    0.315
    ##  .Dep_EpiDoc2 ~~                                                        
    ##    .MM_EpiDoc2        0.000                                        0.000
    ##  .Dep_EpiDoc4 ~~                                                        
    ##    .MM_EpiDoc4        0.000                                        0.000
    ##  ci.upper   Std.lv  Std.all
    ##                            
    ##     1.528    0.922    0.293
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.000    0.000    0.000
    ## 
    ## Intercepts:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##    .Dep_EpiDoc4       1.293    0.247    5.244      Inf    0.000    0.809
    ##    .Dep_EpiDoc2       0.513    0.610    0.841  147.074    0.401   -0.692
    ##    .MM_EpiDoc4       -0.559    0.230   -2.425      Inf    0.015   -1.010
    ##    .MM_EpiDoc2        0.165    0.152    1.084      Inf    0.278   -0.133
    ##    .Dep_EpiDoc1      -0.157    0.685   -0.229      Inf    0.819   -1.500
    ##    .MM_EpiDoc1       -1.335    0.201   -6.648      Inf    0.000   -1.728
    ##  ci.upper   Std.lv  Std.all
    ##     1.776    1.293    0.360
    ##     1.719    0.513    0.140
    ##    -0.107   -0.559   -0.440
    ##     0.464    0.165    0.225
    ##     1.186   -0.157   -0.042
    ##    -0.941   -1.335   -1.272
    ## 
    ## Variances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##    .Dep_EpiDoc4       9.370    1.055    8.880      Inf    0.000    7.300
    ##    .Dep_EpiDoc2      10.759    1.350    7.968  379.009    0.000    8.104
    ##    .MM_EpiDoc4        0.999    0.085   11.762      Inf    0.000    0.833
    ##    .MM_EpiDoc2        0.401    0.061    6.524      Inf    0.000    0.280
    ##    .Dep_EpiDoc1      12.642    1.735    7.286      Inf    0.000    9.241
    ##    .MM_EpiDoc1        0.782    0.122    6.411      Inf    0.000    0.543
    ##  ci.upper   Std.lv  Std.all
    ##    11.440    9.370    0.728
    ##    13.414   10.759    0.795
    ##     1.166    0.999    0.620
    ##     0.521    0.401    0.746
    ##    16.043   12.642    0.907
    ##     1.020    0.782    0.709
    ## 
    ## R-Square:
    ##                    Estimate
    ##     Dep_EpiDoc4       0.272
    ##     Dep_EpiDoc2       0.205
    ##     MM_EpiDoc4        0.380
    ##     MM_EpiDoc2        0.254
    ##     Dep_EpiDoc1       0.093
    ##     MM_EpiDoc1        0.291
    ## 
    ## Defined Parameters:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##     ada               0.015    0.006    2.518  417.171    0.012    0.003
    ##     adb               0.022    0.007    3.027      Inf    0.003    0.008
    ##     amda              0.005    0.006    0.725      Inf    0.468   -0.008
    ##     amdb             -0.006    0.011   -0.541      Inf    0.589   -0.027
    ##     adta              0.019    0.007    2.784  439.633    0.006    0.006
    ##     adtb              0.016    0.012    1.311      Inf    0.190   -0.008
    ##     mda               0.183    0.211    0.867      Inf    0.386   -0.230
    ##     mdb              -0.044    0.226   -0.196      Inf    0.844   -0.487
    ##     dda               0.279    0.079    3.551      Inf    0.000    0.125
    ##     ddb               0.442    0.067    6.625      Inf    0.000    0.311
    ##     eda               0.166    0.066    2.523  287.587    0.012    0.036
    ##     edb               0.164    0.079    2.066      Inf    0.039    0.008
    ##     emda              0.004    0.005    0.846      Inf    0.397   -0.005
    ##     emdb              0.005    0.008    0.637      Inf    0.524   -0.010
    ##     edta              0.171    0.065    2.614  290.232    0.009    0.042
    ##     edtb              0.168    0.081    2.086      Inf    0.037    0.010
    ##  ci.upper   Std.lv  Std.all
    ##     0.026    0.015    0.073
    ##     0.037    0.022    0.108
    ##     0.017    0.005    0.023
    ##     0.016   -0.006   -0.019
    ##     0.033    0.019    0.096
    ##     0.041    0.016    0.089
    ##     0.596    0.184    0.061
    ##     0.399   -0.043   -0.013
    ##     0.433    0.279    0.296
    ##     0.572    0.442    0.459
    ##     0.296    0.166    0.062
    ##     0.319    0.165    0.051
    ##     0.013    0.004    0.002
    ##     0.020    0.005    0.002
    ##     0.299    0.171    0.085
    ##     0.327    0.169    0.053

``` r
anova(resFre_g,  res_covR_g,res_covAge_g,res_covAgeEd_g,res_covAgeEd2_g)
```

    ## Robust correction can only be applied to pooled chi-squared statistic, not F statistic. "asymptotic" was switched to TRUE.

    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.

    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.

    ## Robust correction can only be applied to pooled chi-squared statistic, not F statistic. "asymptotic" was switched to TRUE.

    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.
    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.

    ## Robust correction can only be applied to pooled chi-squared statistic, not F statistic. "asymptotic" was switched to TRUE.

    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.
    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.

    ## Robust correction can only be applied to pooled chi-squared statistic, not F statistic. "asymptotic" was switched to TRUE.

    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.
    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.
    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.
    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.
    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.
    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.
    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.

    ## Test statistic(s) pooled using the D4 pooling method.
    ##   Pooled statistic: "standard"
    ##   Method to robustify pooled statistic:  "yuan.bentler.mplus"
    ## 
    ##                 Df   Chisq Chisq diff Df diff Pr(>Chisq)     RIV      FMI
    ## resFre_g         2  2.3128                                               
    ## res_covR_g       6  5.1723     2.4926       4    0.64597 0.34350 0.255676
    ## res_covAge_g     8  7.7007     2.1111       2    0.34799 0.35321 0.261018
    ## res_covAgeEd_g  12 10.6192     2.7331       4    0.60343 0.08429 0.077737
    ## res_covAgeEd2_g 16 15.4667     4.8023       4    0.30819 0.00000 0.000000

difference is still significant, constraining to make the effect of age
being the same among groups leads to worsening of the fit.

``` r
source('./code/parTableToCSV.R')
res<-parTableToCSV(res_covAgeEd2_g,path=('./Results/LaggedSEM_mi20/gSEX_covAgeEd'))
```

    ## Joining with `by = join_by(group)`

    ## Warning: There was 1 warning in `mutate()`.
    ## ℹ In argument: `Group2 = as.numeric(substr(label, 2, 2))`.
    ## Caused by warning:
    ## ! NAs introduced by coercion

    ##    Group2 label     estimate_p             CI Group
    ## 1      NA     a  0.015 (0.012)  [0.003:0.026]  <NA>
    ## 2      NA     b  0.022 (0.003)  [0.008:0.037]  <NA>
    ## 3      NA    da  0.005 (0.468) [-0.008:0.017]  <NA>
    ## 4      NA    db -0.006 (0.589) [-0.027:0.016]  <NA>
    ## 5      NA    ta  0.019 (0.006)  [0.006:0.033]  <NA>
    ## 6      NA    tb   0.016 (0.19) [-0.008:0.041]  <NA>
    ## 7      NA     a  0.183 (0.386)  [-0.23:0.596]  <NA>
    ## 8      NA     b -0.044 (0.844) [-0.487:0.399]  <NA>
    ## 9      NA     a      0.279 (0)  [0.125:0.433]  <NA>
    ## 10     NA     b      0.442 (0)  [0.311:0.572]  <NA>
    ## 11     NA     a  0.166 (0.012)  [0.036:0.296]  <NA>
    ## 12     NA     b  0.164 (0.039)  [0.008:0.319]  <NA>
    ## 13     NA    da  0.004 (0.397) [-0.005:0.013]  <NA>
    ## 14     NA    db  0.005 (0.524)   [-0.01:0.02]  <NA>
    ## 15     NA    ta  0.171 (0.009)  [0.042:0.299]  <NA>
    ## 16     NA    tb  0.168 (0.037)   [0.01:0.327]  <NA>

    ## Warning in reshapeWide(data, idvar = idvar, timevar = timevar, varying =
    ## varying, : there are records with missing times, which will be dropped.

    ## Warning in reshapeWide(data, idvar = idvar, timevar = timevar, varying =
    ## varying, : multiple rows match for Group=NA: first taken

``` r
res$gPartable
```

    ##                      label2 estimate_p_Male estimate_p_Female           CI_Male
    ## 1   Dep_EpiDoc4~Dep_EpiDoc1    0.17 (0.031)         0.351 (0)     [0.015:0.325]
    ## 2   Dep_EpiDoc4~Dep_EpiDoc2   0.247 (0.002)     0.227 (0.003)     [0.089:0.405]
    ## 3   Dep_EpiDoc2~Dep_EpiDoc1       0.404 (0)         0.373 (0)     [0.266:0.542]
    ## 4     MM_EpiDoc4~MM_EpiDoc2   0.369 (0.016)     0.285 (0.012)     [0.068:0.671]
    ## 5     MM_EpiDoc4~MM_EpiDoc1        0.56 (0)         0.347 (0)     [0.385:0.736]
    ## 6     MM_EpiDoc2~MM_EpiDoc1   0.102 (0.029)         0.365 (0)      [0.01:0.194]
    ## 7    Dep_EpiDoc4~MM_EpiDoc2   0.494 (0.083)     0.675 (0.028)    [-0.065:1.053]
    ## 8    Dep_EpiDoc4~MM_EpiDoc1   0.077 (0.697)    -0.309 (0.152)    [-0.312:0.466]
    ## 9    Dep_EpiDoc2~MM_EpiDoc1    0.23 (0.299)     0.086 (0.727)    [-0.205:0.664]
    ## 10   MM_EpiDoc4~Dep_EpiDoc2   0.035 (0.244)     0.001 (0.955)    [-0.024:0.094]
    ## 11   MM_EpiDoc4~Dep_EpiDoc1   0.004 (0.879)     0.028 (0.118)    [-0.052:0.061]
    ## 12   MM_EpiDoc2~Dep_EpiDoc1   0.018 (0.228)      0.009 (0.46)    [-0.011:0.048]
    ## 13  Dep_EpiDoc1~~MM_EpiDoc1   0.549 (0.011)     0.922 (0.003)     [0.126:0.973]
    ## 14  Dep_EpiDoc2~~MM_EpiDoc2          0 (NA)            0 (NA)             [0:0]
    ## 15  Dep_EpiDoc4~~MM_EpiDoc4          0 (NA)            0 (NA)             [0:0]
    ## 16          Dep_EpiDoc4~Age          0 (NA)            0 (NA)             [0:0]
    ## 17          Dep_EpiDoc2~Age    0.01 (0.463)      0.01 (0.463)    [-0.016:0.035]
    ## 18          Dep_EpiDoc1~Age   0.046 (0.002)     0.046 (0.002)     [0.017:0.074]
    ## 19           MM_EpiDoc4~Age   0.018 (0.011)         0.031 (0)     [0.004:0.033]
    ## 20           MM_EpiDoc2~Age       0 (0.914)     -0.005 (0.29)    [-0.003:0.004]
    ## 21           MM_EpiDoc1~Age       0.032 (0)         0.049 (0)     [0.023:0.041]
    ## 22    Dep_EpiDoc4~Education          0 (NA)            0 (NA)             [0:0]
    ## 23    Dep_EpiDoc2~Education    0.49 (0.001)      0.435 (0.04)     [0.202:0.777]
    ## 24    Dep_EpiDoc1~Education   0.448 (0.002)         0.776 (0)     [0.161:0.735]
    ## 25     MM_EpiDoc4~Education          0 (NA)            0 (NA)             [0:0]
    ## 26     MM_EpiDoc2~Education          0 (NA)            0 (NA)             [0:0]
    ## 27     MM_EpiDoc1~Education          0 (NA)            0 (NA)             [0:0]
    ## 28 Dep_EpiDoc4~~Dep_EpiDoc4        6.91 (0)          9.37 (0)     [4.942:8.878]
    ## 29 Dep_EpiDoc2~~Dep_EpiDoc2       6.201 (0)        10.759 (0)     [4.487:7.914]
    ## 30   MM_EpiDoc4~~MM_EpiDoc4       0.971 (0)         0.999 (0)     [0.764:1.178]
    ## 31   MM_EpiDoc2~~MM_EpiDoc2       0.189 (0)         0.401 (0)      [0.087:0.29]
    ## 32 Dep_EpiDoc1~~Dep_EpiDoc1       8.431 (0)        12.642 (0)    [6.156:10.706]
    ## 33   MM_EpiDoc1~~MM_EpiDoc1       0.701 (0)         0.782 (0)     [0.531:0.871]
    ## 34                 Age~~Age    204.021 (NA)      135.726 (NA) [204.021:204.021]
    ## 35           Age~~Education      6.794 (NA)         3.68 (NA)     [6.794:6.794]
    ## 36     Education~~Education      1.159 (NA)         1.24 (NA)     [1.159:1.159]
    ## 37            Dep_EpiDoc4~1       1.461 (0)         1.293 (0)     [1.042:1.879]
    ## 38            Dep_EpiDoc2~1  -0.261 (0.646)     0.513 (0.401)    [-1.383:0.862]
    ## 39             MM_EpiDoc4~1  -0.362 (0.162)    -0.559 (0.015)    [-0.869:0.145]
    ## 40             MM_EpiDoc2~1   -0.02 (0.756)     0.165 (0.278)    [-0.144:0.105]
    ## 41            Dep_EpiDoc1~1   -0.657 (0.29)    -0.157 (0.819)    [-1.872:0.559]
    ## 42             MM_EpiDoc1~1      -0.785 (0)        -1.335 (0)   [-1.129:-0.441]
    ## 43                    Age~1     44.188 (NA)       42.481 (NA)   [44.188:44.188]
    ## 44              Education~1      2.823 (NA)        2.409 (NA)     [2.823:2.823]
    ##            CI_Female
    ## 1      [0.205:0.496]
    ## 2       [0.08:0.374]
    ## 3      [0.208:0.538]
    ## 4      [0.064:0.506]
    ## 5      [0.176:0.519]
    ## 6      [0.177:0.553]
    ## 7      [0.073:1.277]
    ## 8     [-0.731:0.114]
    ## 9     [-0.399:0.572]
    ## 10    [-0.032:0.034]
    ## 11    [-0.007:0.062]
    ## 12    [-0.015:0.034]
    ## 13     [0.315:1.528]
    ## 14             [0:0]
    ## 15             [0:0]
    ## 16             [0:0]
    ## 17    [-0.016:0.035]
    ## 18     [0.017:0.074]
    ## 19     [0.017:0.044]
    ## 20    [-0.013:0.004]
    ## 21     [0.038:0.059]
    ## 22             [0:0]
    ## 23       [0.02:0.85]
    ## 24     [0.389:1.162]
    ## 25             [0:0]
    ## 26             [0:0]
    ## 27             [0:0]
    ## 28       [7.3:11.44]
    ## 29    [8.104:13.414]
    ## 30     [0.833:1.166]
    ## 31      [0.28:0.521]
    ## 32    [9.241:16.043]
    ## 33      [0.543:1.02]
    ## 34 [135.726:135.726]
    ## 35       [3.68:3.68]
    ## 36       [1.24:1.24]
    ## 37     [0.809:1.776]
    ## 38    [-0.692:1.719]
    ## 39    [-1.01:-0.107]
    ## 40    [-0.133:0.464]
    ## 41      [-1.5:1.186]
    ## 42   [-1.728:-0.941]
    ## 43   [42.481:42.481]
    ## 44     [2.409:2.409]

# testing moderation

``` r
model_mod<-'
# evolução
Dep_EpiDoc4~Dep_EpiDoc1+Dep_EpiDoc2
Dep_EpiDoc2~Dep_EpiDoc1

MM_EpiDoc4~MM_EpiDoc2+MM_EpiDoc1
MM_EpiDoc2~MM_EpiDoc1 

# causalidade - c/ temporalidade/lag
Dep_EpiDoc4~MM_EpiDoc2+MM_EpiDoc1
Dep_EpiDoc2~MM_EpiDoc1

MM_EpiDoc4~Dep_EpiDoc2+Dep_EpiDoc1
MM_EpiDoc2~Dep_EpiDoc1



# corrrelação mesmo t
MM_EpiDoc1~~Dep_EpiDoc1
MM_EpiDoc2~~0*Dep_EpiDoc2
MM_EpiDoc4~~0*Dep_EpiDoc4

# causalidade - idade
Dep_EpiDoc4~0*Age
Dep_EpiDoc2~Age
Dep_EpiDoc1~Age

MM_EpiDoc4~0*Age
MM_EpiDoc2~Age
MM_EpiDoc1~Age


# educação
Dep_EpiDoc4~0*Education
Dep_EpiDoc2~Education
Dep_EpiDoc1~Education
MM_EpiDoc4~0*Education
MM_EpiDoc2~Education
MM_EpiDoc1~Education


# interaction terms
#Dep_EpiDoc4~MM2_ed+MM1_ed
#Dep_EpiDoc2~MM1_ed



MM_EpiDoc4~d2_ed+d2_ed
MM_EpiDoc2~d1_ed



MM_EpiDoc4~age_ed
MM_EpiDoc2~age_ed
MM_EpiDoc1~age_ed

# correlation among interaction terms and its original var
#MM2_ed~~Education+MM_EpiDoc2+age_ed+Age+MM1_ed+MM_EpiDoc1+d2_ed+d1_ed
#MM1_ed~~Education+MM_EpiDoc1+age_ed+Age+MM_EpiDoc2+d2_ed+d1_ed
d2_ed~~Education+Dep_EpiDoc2+age_ed+d1_ed+Age+Dep_EpiDoc1
d1_ed~~Education+Dep_EpiDoc1+age_ed+Age+Dep_EpiDoc2
age_ed~~Education+Age
Education~~Age

'


res_mod<-lavaan.mi::sem.mi(data=data.mi,
         model=model_mod,
         estimator = "ML", test = "yuan.bentler.mplus",
         likelihood = "wishart",
         sampling.weights='ipw',group='Sex',
         orthogonal=FALSE
         )
```

    ## Warning: lavaan->lav_model_vcov():  
    ##    The variance-covariance matrix of the estimated parameters (vcov) does not 
    ##    appear to be positive definite! The smallest eigenvalue (= -1.048120e-14) 
    ##    is smaller than zero. This may be a symptom that the model is not 
    ##    identified.
    ## Warning: lavaan->lav_model_vcov():  
    ##    The variance-covariance matrix of the estimated parameters (vcov) does not 
    ##    appear to be positive definite! The smallest eigenvalue (= -1.048120e-14) 
    ##    is smaller than zero. This may be a symptom that the model is not 
    ##    identified.

    ## Warning: lavaan->lav_model_vcov():  
    ##    The variance-covariance matrix of the estimated parameters (vcov) does not 
    ##    appear to be positive definite! The smallest eigenvalue (= -1.594779e-15) 
    ##    is smaller than zero. This may be a symptom that the model is not 
    ##    identified.

    ## Warning: lavaan->lav_model_vcov():  
    ##    The variance-covariance matrix of the estimated parameters (vcov) does not 
    ##    appear to be positive definite! The smallest eigenvalue (= -2.678873e-14) 
    ##    is smaller than zero. This may be a symptom that the model is not 
    ##    identified.

    ## Warning: lavaan->lav_model_vcov():  
    ##    The variance-covariance matrix of the estimated parameters (vcov) does not 
    ##    appear to be positive definite! The smallest eigenvalue (= -9.995965e-15) 
    ##    is smaller than zero. This may be a symptom that the model is not 
    ##    identified.

    ## Warning: lavaan->lav_model_vcov():  
    ##    The variance-covariance matrix of the estimated parameters (vcov) does not 
    ##    appear to be positive definite! The smallest eigenvalue (= -7.502046e-15) 
    ##    is smaller than zero. This may be a symptom that the model is not 
    ##    identified.

    ## Warning: lavaan->lav_model_vcov():  
    ##    The variance-covariance matrix of the estimated parameters (vcov) does not 
    ##    appear to be positive definite! The smallest eigenvalue (= -1.758755e-15) 
    ##    is smaller than zero. This may be a symptom that the model is not 
    ##    identified.

    ## Warning: lavaan->lav_model_vcov():  
    ##    The variance-covariance matrix of the estimated parameters (vcov) does not 
    ##    appear to be positive definite! The smallest eigenvalue (= -4.613177e-15) 
    ##    is smaller than zero. This may be a symptom that the model is not 
    ##    identified.

    ## Warning: lavaan->lav_model_vcov():  
    ##    The variance-covariance matrix of the estimated parameters (vcov) does not 
    ##    appear to be positive definite! The smallest eigenvalue (= -6.559778e-15) 
    ##    is smaller than zero. This may be a symptom that the model is not 
    ##    identified.

    ## Warning: lavaan->lav_model_vcov():  
    ##    The variance-covariance matrix of the estimated parameters (vcov) does not 
    ##    appear to be positive definite! The smallest eigenvalue (= -7.202657e-15) 
    ##    is smaller than zero. This may be a symptom that the model is not 
    ##    identified.

    ## Warning: lavaan->lav_model_vcov():  
    ##    The variance-covariance matrix of the estimated parameters (vcov) does not 
    ##    appear to be positive definite! The smallest eigenvalue (= -6.917265e-15) 
    ##    is smaller than zero. This may be a symptom that the model is not 
    ##    identified.

    ## Warning: lavaan->lav_model_vcov():  
    ##    The variance-covariance matrix of the estimated parameters (vcov) does not 
    ##    appear to be positive definite! The smallest eigenvalue (= -3.426988e-15) 
    ##    is smaller than zero. This may be a symptom that the model is not 
    ##    identified.

    ## Warning: lavaan->lav_model_vcov():  
    ##    The variance-covariance matrix of the estimated parameters (vcov) does not 
    ##    appear to be positive definite! The smallest eigenvalue (= -5.528988e-16) 
    ##    is smaller than zero. This may be a symptom that the model is not 
    ##    identified.

    ## Warning: lavaan->lav_model_vcov():  
    ##    The variance-covariance matrix of the estimated parameters (vcov) does not 
    ##    appear to be positive definite! The smallest eigenvalue (= -6.201479e-16) 
    ##    is smaller than zero. This may be a symptom that the model is not 
    ##    identified.

    ## Warning: lavaan->lav_model_vcov():  
    ##    The variance-covariance matrix of the estimated parameters (vcov) does not 
    ##    appear to be positive definite! The smallest eigenvalue (= -1.061326e-14) 
    ##    is smaller than zero. This may be a symptom that the model is not 
    ##    identified.

    ## Warning: lavaan->lav_model_vcov():  
    ##    The variance-covariance matrix of the estimated parameters (vcov) does not 
    ##    appear to be positive definite! The smallest eigenvalue (= -2.850728e-17) 
    ##    is smaller than zero. This may be a symptom that the model is not 
    ##    identified.

    ## Warning: lavaan->lav_model_vcov():  
    ##    The variance-covariance matrix of the estimated parameters (vcov) does not 
    ##    appear to be positive definite! The smallest eigenvalue (= -4.132677e-15) 
    ##    is smaller than zero. This may be a symptom that the model is not 
    ##    identified.

    ## Warning: lavaan->lav_model_vcov():  
    ##    The variance-covariance matrix of the estimated parameters (vcov) does not 
    ##    appear to be positive definite! The smallest eigenvalue (= -9.506669e-16) 
    ##    is smaller than zero. This may be a symptom that the model is not 
    ##    identified.

    ## Warning: lavaan->lav_model_vcov():  
    ##    The variance-covariance matrix of the estimated parameters (vcov) does not 
    ##    appear to be positive definite! The smallest eigenvalue (= -5.663202e-17) 
    ##    is smaller than zero. This may be a symptom that the model is not 
    ##    identified.

    ## Warning: lavaan->lav_model_vcov():  
    ##    The variance-covariance matrix of the estimated parameters (vcov) does not 
    ##    appear to be positive definite! The smallest eigenvalue (= -4.468342e-15) 
    ##    is smaller than zero. This may be a symptom that the model is not 
    ##    identified.

    ## Warning: lavaan->lav_model_vcov():  
    ##    The variance-covariance matrix of the estimated parameters (vcov) does not 
    ##    appear to be positive definite! The smallest eigenvalue (= -2.423886e-14) 
    ##    is smaller than zero. This may be a symptom that the model is not 
    ##    identified.

``` r
summary(res_mod,fit.meas=TRUE,standardized=TRUE,rsquare=TRUE,ci=TRUE)
```

    ## Warning: lavaan->lav_model_vcov():  
    ##    The variance-covariance matrix of the estimated parameters (vcov) does not 
    ##    appear to be positive definite! The smallest eigenvalue (= -1.048120e-14) 
    ##    is smaller than zero. This may be a symptom that the model is not 
    ##    identified.

    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.

    ## Warning in vcov_lavaan_mi(object, scale.W = scale.W, omit.imps = omit.imps): Could not invert within-imputation covariance matrix. Generalized inverse used instead.
    ## It may be safer to set `scale.W = FALSE' (and `asymptotic = TRUE').

    ## Warning: lavaan->lav_model_vcov():  
    ##    The variance-covariance matrix of the estimated parameters (vcov) does not 
    ##    appear to be positive definite! The smallest eigenvalue (= -1.048120e-14) 
    ##    is smaller than zero. This may be a symptom that the model is not 
    ##    identified.

    ## Warning in vcov_lavaan_mi(object, omit.imps = omit.imps): Could not invert within-imputation covariance matrix. Generalized inverse used instead.
    ## It may be safer to set `scale.W = FALSE' (and `asymptotic = TRUE').

    ## Warning: lavaan->lav_model_vcov():  
    ##    The variance-covariance matrix of the estimated parameters (vcov) does not 
    ##    appear to be positive definite! The smallest eigenvalue (= -1.048120e-14) 
    ##    is smaller than zero. This may be a symptom that the model is not 
    ##    identified.

    ## Warning in vcov_lavaan_mi(object, omit.imps = omit.imps): Could not invert within-imputation covariance matrix. Generalized inverse used instead.
    ## It may be safer to set `scale.W = FALSE' (and `asymptotic = TRUE').

    ## Warning: lavaan->lav_model_vcov():  
    ##    The variance-covariance matrix of the estimated parameters (vcov) does not 
    ##    appear to be positive definite! The smallest eigenvalue (= -1.048120e-14) 
    ##    is smaller than zero. This may be a symptom that the model is not 
    ##    identified.

    ## Warning in vcov_lavaan_mi(object, omit.imps = omit.imps): Could not invert within-imputation covariance matrix. Generalized inverse used instead.
    ## It may be safer to set `scale.W = FALSE' (and `asymptotic = TRUE').

    ## lavaan.mi object fit to 20 imputed data sets using:
    ##  - lavaan    (0.6-19)
    ##  - lavaan.mi (0.1-0)
    ## See class?lavaan.mi help page for available methods. 
    ## 
    ## Convergence information:
    ## The model converged on 20 imputed data sets.
    ## Standard errors were available for all imputations.
    ## 
    ##   Estimator                                         ML
    ##   Optimization method                           NLMINB
    ##   Number of model parameters                       124
    ## 
    ##   Number of observations per group:                   
    ##     Male                                           355
    ##     Female                                         317
    ##   Sampling weights variable                        ipw
    ## 
    ## Model Test User Model:
    ## 
    ##                                                     Standard      Scaled
    ##   Test statistic                                      87.890      70.930
    ##   Degrees of freedom                                      30          30
    ##   P-value                                              0.000       0.000
    ##   Average scaling correction factor                                1.239
    ##   Pooling method                                          D4            
    ##     Pooled statistic                              "standard"            
    ##     "yuan.bentler.mplus" correction applied            AFTER     pooling
    ## 
    ## Model Test Baseline Model:
    ## 
    ##   Test statistic                              5413.458    3535.216
    ##   Degrees of freedom                               110         110
    ##   P-value                                        0.000       0.000
    ##   Scaling correction factor                                  1.531
    ## 
    ## User Model versus Baseline Model:
    ## 
    ##   Comparative Fit Index (CFI)                    0.989       0.988
    ##   Tucker-Lewis Index (TLI)                       0.960       0.956
    ##                                                                   
    ##   Robust Comparative Fit Index (CFI)                         0.990
    ##   Robust Tucker-Lewis Index (TLI)                            0.965
    ## 
    ## Loglikelihood and Information Criteria:
    ## 
    ##   Loglikelihood user model (H0)             -11039.319  -11039.319
    ##   Scaling correction factor                                  1.975
    ##       for the MLR correction                                      
    ##   Loglikelihood unrestricted model (H1)     -10983.802  -10983.802
    ##   Scaling correction factor                                  1.819
    ##       for the MLR correction                                      
    ##                                                                   
    ##   Akaike (AIC)                               22326.638   22326.638
    ##   Bayesian (BIC)                             22885.910   22885.910
    ##   Sample-size adjusted Bayesian (SABIC)      22492.200   22492.200
    ## 
    ## Root Mean Square Error of Approximation:
    ## 
    ##   RMSEA                                          0.076       0.067
    ##   90 Percent confidence interval - lower         0.058       0.049
    ##   90 Percent confidence interval - upper         0.095       0.084
    ##   P-value H_0: RMSEA <= 0.050                    0.011       0.056
    ##   P-value H_0: RMSEA >= 0.080                    0.377       0.111
    ##                                                                   
    ##   Robust RMSEA                                               0.071
    ##   90 Percent confidence interval - lower                     0.050
    ##   90 Percent confidence interval - upper                     0.093
    ##   P-value H_0: Robust RMSEA <= 0.050                         0.052
    ##   P-value H_0: Robust RMSEA >= 0.080                         0.261
    ## 
    ## Standardized Root Mean Square Residual:
    ## 
    ##   SRMR                                           0.052       0.052
    ## 
    ## Parameter Estimates:
    ## 
    ##   Standard errors                                    Sandwich
    ##   Information bread                                  Observed
    ##   Observed information based on                       Hessian
    ##                                                              
    ##   Pooled across imputations              Rubin's (1987) rules
    ##   Augment within-imputation variance     Scale by average RIV
    ##   Wald test for pooled parameters          t(df) distribution
    ## 
    ##   Pooled t statistics with df >= 1000 are displayed with
    ##   df = Inf(inity) to save space. Although the t distribution
    ##   with large df closely approximates a standard normal
    ##   distribution, exact df for reporting these t tests can be
    ##   obtained from parameterEstimates.mi() 
    ## 
    ## 
    ## 
    ## Group 1 [Male]:
    ## 
    ## Regressions:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##   Dep_EpiDoc4 ~                                                         
    ##     Dep_EpiDoc1       0.170    0.081    2.098  563.969    0.036    0.011
    ##     Dep_EpiDoc2       0.247    0.083    2.993  334.574    0.003    0.085
    ##   Dep_EpiDoc2 ~                                                         
    ##     Dep_EpiDoc1       0.425    0.071    5.948  279.870    0.000    0.284
    ##   MM_EpiDoc4 ~                                                          
    ##     MM_EpiDoc2        0.409    0.173    2.368      Inf    0.018    0.071
    ##     MM_EpiDoc1        0.623    0.090    6.913      Inf    0.000    0.447
    ##   MM_EpiDoc2 ~                                                          
    ##     MM_EpiDoc1        0.101    0.046    2.170      Inf    0.030    0.010
    ##   Dep_EpiDoc4 ~                                                         
    ##     MM_EpiDoc2        0.494    0.293    1.684      Inf    0.092   -0.081
    ##     MM_EpiDoc1        0.077    0.204    0.378      Inf    0.705   -0.323
    ##   Dep_EpiDoc2 ~                                                         
    ##     MM_EpiDoc1        0.003    0.067    0.052      Inf    0.959   -0.127
    ##   MM_EpiDoc4 ~                                                          
    ##     Dep_EpiDoc2      -0.046    0.071   -0.653      Inf    0.514   -0.186
    ##     Dep_EpiDoc1       0.006    0.029    0.200      Inf    0.841   -0.050
    ##   MM_EpiDoc2 ~                                                          
    ##     Dep_EpiDoc1       0.010    0.038    0.258      Inf    0.796   -0.065
    ##   Dep_EpiDoc4 ~                                                         
    ##     Age               0.000                                        0.000
    ##   Dep_EpiDoc2 ~                                                         
    ##     Age               0.008    0.015    0.574   89.440    0.567   -0.021
    ##   Dep_EpiDoc1 ~                                                         
    ##     Age               0.041    0.019    2.135      Inf    0.033    0.003
    ##   MM_EpiDoc4 ~                                                          
    ##     Age               0.000                                        0.000
    ##   MM_EpiDoc2 ~                                                          
    ##     Age               0.016    0.008    1.958      Inf    0.050   -0.000
    ##   MM_EpiDoc1 ~                                                          
    ##     Age               0.032    0.010    3.070      Inf    0.002    0.012
    ##   Dep_EpiDoc4 ~                                                         
    ##     Education         0.000                                        0.000
    ##   Dep_EpiDoc2 ~                                                         
    ##     Education         0.529    0.151    3.497  463.569    0.001    0.232
    ##   Dep_EpiDoc1 ~                                                         
    ##     Education         0.495    0.169    2.931      Inf    0.003    0.164
    ##   MM_EpiDoc4 ~                                                          
    ##     Education         0.000                                        0.000
    ##   MM_EpiDoc2 ~                                                          
    ##     Education         0.162    0.089    1.820      Inf    0.069   -0.013
    ##   MM_EpiDoc1 ~                                                          
    ##     Education         0.047    0.160    0.292      Inf    0.770   -0.267
    ##   MM_EpiDoc4 ~                                                          
    ##     d2_ed             0.260    0.269    0.968      Inf    0.333   -0.267
    ##   MM_EpiDoc2 ~                                                          
    ##     d1_ed             0.037    0.157    0.236      Inf    0.813   -0.271
    ##   MM_EpiDoc4 ~                                                          
    ##     age_ed            0.105    0.098    1.070      Inf    0.285   -0.087
    ##   MM_EpiDoc2 ~                                                          
    ##     age_ed           -0.343    0.168   -2.040      Inf    0.041   -0.672
    ##   MM_EpiDoc1 ~                                                          
    ##     age_ed           -0.025    0.267   -0.093      Inf    0.926   -0.548
    ##  ci.upper   Std.lv  Std.all
    ##                            
    ##     0.329    0.170    0.181
    ##     0.410    0.247    0.256
    ##                            
    ##     0.566    0.425    0.438
    ##                            
    ##     0.748    0.409    0.146
    ##     0.800    0.623    0.473
    ##                            
    ##     0.192    0.101    0.214
    ##                            
    ##     1.069    0.494    0.077
    ##     0.477    0.077    0.026
    ##                            
    ##     0.134    0.003    0.001
    ##                            
    ##     0.093   -0.046   -0.110
    ##     0.062    0.006    0.014
    ##                            
    ##     0.084    0.010    0.067
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.037    0.008    0.040
    ##                            
    ##     0.079    0.041    0.191
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.032    0.016    0.512
    ##                            
    ##     0.053    0.032    0.481
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.827    0.529    0.192
    ##                            
    ##     0.827    0.495    0.174
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.337    0.162    0.387
    ##                            
    ##     0.360    0.047    0.053
    ##                            
    ##     0.788    0.260    0.210
    ##                            
    ##     0.345    0.037    0.083
    ##                            
    ##     0.297    0.105    0.092
    ##                            
    ##    -0.013   -0.343   -0.842
    ##                            
    ##     0.498   -0.025   -0.029
    ## 
    ## Covariances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##  .Dep_EpiDoc1 ~~                                                        
    ##    .MM_EpiDoc1       -0.013    0.058   -0.216      Inf    0.829   -0.126
    ##  .Dep_EpiDoc2 ~~                                                        
    ##    .MM_EpiDoc2        0.000                                        0.000
    ##  .Dep_EpiDoc4 ~~                                                        
    ##    .MM_EpiDoc4        0.000                                        0.000
    ##   Education ~~                                                          
    ##     d2_ed             0.551    0.090    6.137      Inf    0.000    0.375
    ##  .Dep_EpiDoc2 ~~                                                        
    ##     d2_ed             1.841    0.277    6.653  158.833    0.000    1.294
    ##   d2_ed ~~                                                              
    ##     age_ed            0.546    0.107    5.116  414.492    0.000    0.336
    ##     d1_ed             0.594    0.138    4.309      Inf    0.000    0.324
    ##   Age ~~                                                                
    ##     d2_ed             4.365    1.221    3.576  170.410    0.000    1.956
    ##  .Dep_EpiDoc1 ~~                                                        
    ##     d2_ed             1.062    0.300    3.546  888.649    0.000    0.474
    ##   Education ~~                                                          
    ##     d1_ed             0.480    0.074    6.504      Inf    0.000    0.336
    ##  .Dep_EpiDoc1 ~~                                                        
    ##     d1_ed             2.487    0.404    6.163      Inf    0.000    1.696
    ##   d1_ed ~~                                                              
    ##     age_ed            0.515    0.082    6.249      Inf    0.000    0.353
    ##   Age ~~                                                                
    ##     d1_ed             4.754    1.063    4.473      Inf    0.000    2.671
    ##  .Dep_EpiDoc2 ~~                                                        
    ##     d1_ed             0.015    0.044    0.348  455.806    0.728   -0.071
    ##   Education ~~                                                          
    ##     age_ed            1.013    0.079   12.879      Inf    0.000    0.859
    ##   Age ~~                                                                
    ##     age_ed           12.895    1.757    7.337      Inf    0.000    9.450
    ##     Education         6.794    1.016    6.686      Inf    0.000    4.802
    ##  ci.upper   Std.lv  Std.all
    ##                            
    ##     0.101   -0.013   -0.005
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.728    0.551    0.506
    ##                            
    ##     2.387    1.841    0.727
    ##                            
    ##     0.756    0.546    0.488
    ##     0.864    0.594    0.580
    ##                            
    ##     6.774    4.365    0.302
    ##                            
    ##     1.650    1.062    0.361
    ##                            
    ##     0.625    0.480    0.442
    ##                            
    ##     3.278    2.487    0.847
    ##                            
    ##     0.676    0.515    0.461
    ##                            
    ##     6.838    4.754    0.329
    ##                            
    ##     0.102    0.015    0.006
    ##                            
    ##     1.167    1.013    0.851
    ##                            
    ##    16.339   12.895    0.816
    ##     8.786    6.794    0.442
    ## 
    ## Intercepts:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##    .Dep_EpiDoc4       1.461    0.219    6.668  401.475    0.000    1.030
    ##    .Dep_EpiDoc2      -0.230    0.611   -0.376  113.038    0.707   -1.440
    ##    .MM_EpiDoc4        0.584    0.207    2.827      Inf    0.005    0.179
    ##    .MM_EpiDoc2       -1.066    0.591   -1.805      Inf    0.071   -2.225
    ##    .Dep_EpiDoc1      -0.575    0.740   -0.777      Inf    0.437   -2.025
    ##    .MM_EpiDoc1       -0.901    0.790   -1.142      Inf    0.254   -2.449
    ##     Age              44.188    1.006   43.940      Inf    0.000   42.217
    ##     Education         2.823    0.064   43.810      Inf    0.000    2.697
    ##     d2_ed             0.059    0.067    0.884      Inf    0.377   -0.072
    ##     d1_ed             0.001    0.064    0.014      Inf    0.989   -0.125
    ##     age_ed            0.267    0.076    3.509      Inf    0.000    0.118
    ##  ci.upper   Std.lv  Std.all
    ##     1.891    1.461    0.509
    ##     0.980   -0.230   -0.077
    ##     0.989    0.584    0.464
    ##     0.092   -1.066   -2.369
    ##     0.875   -0.575   -0.188
    ##     0.646   -0.901   -0.944
    ##    46.159   44.188    3.094
    ##     2.949    2.823    2.622
    ##     0.190    0.059    0.058
    ##     0.127    0.001    0.001
    ##     0.416    0.267    0.242
    ## 
    ## Variances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##    .Dep_EpiDoc4       6.910    1.032    6.695      Inf    0.000    4.886
    ##    .Dep_EpiDoc2       6.248    0.886    7.056  108.076    0.000    4.493
    ##    .MM_EpiDoc4        0.994    0.108    9.172      Inf    0.000    0.782
    ##    .MM_EpiDoc2        0.182    0.051    3.563      Inf    0.000    0.082
    ##    .Dep_EpiDoc1       8.447    1.208    6.992      Inf    0.000    6.079
    ##    .MM_EpiDoc1        0.700    0.088    7.913      Inf    0.000    0.527
    ##     Age             204.022   23.394    8.721      Inf    0.000  158.170
    ##     Education         1.159    0.073   15.954      Inf    0.000    1.017
    ##     d2_ed             1.025    0.182    5.617  930.617    0.000    0.667
    ##     d1_ed             1.021    0.177    5.761      Inf    0.000    0.673
    ##     age_ed            1.223    0.136    8.997      Inf    0.000    0.957
    ##  ci.upper   Std.lv  Std.all
    ##     8.935    6.910    0.838
    ##     8.004    6.248    0.709
    ##     1.207    0.994    0.628
    ##     0.282    0.182    0.898
    ##    10.815    8.447    0.904
    ##     0.874    0.700    0.768
    ##   249.873  204.022    1.000
    ##     1.302    1.159    1.000
    ##     1.383    1.025    1.000
    ##     1.368    1.021    1.000
    ##     1.489    1.223    1.000
    ## 
    ## R-Square:
    ##                    Estimate
    ##     Dep_EpiDoc4       0.162
    ##     Dep_EpiDoc2       0.291
    ##     MM_EpiDoc4        0.372
    ##     MM_EpiDoc2        0.102
    ##     Dep_EpiDoc1       0.096
    ##     MM_EpiDoc1        0.232
    ##     Age               0.000
    ##     Education         0.000
    ##     d2_ed             0.000
    ##     d1_ed             0.000
    ##     age_ed            0.000
    ## 
    ## 
    ## Group 2 [Female]:
    ## 
    ## Regressions:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##   Dep_EpiDoc4 ~                                                         
    ##     Dep_EpiDoc1       0.351    0.077    4.579      Inf    0.000    0.200
    ##     Dep_EpiDoc2       0.227    0.077    2.940      Inf    0.003    0.076
    ##   Dep_EpiDoc2 ~                                                         
    ##     Dep_EpiDoc1       0.360    0.085    4.261  697.240    0.000    0.194
    ##   MM_EpiDoc4 ~                                                          
    ##     MM_EpiDoc2        0.253    0.120    2.110      Inf    0.035    0.018
    ##     MM_EpiDoc1        0.450    0.081    5.557      Inf    0.000    0.292
    ##   MM_EpiDoc2 ~                                                          
    ##     MM_EpiDoc1        0.359    0.098    3.685      Inf    0.000    0.168
    ##   Dep_EpiDoc4 ~                                                         
    ##     MM_EpiDoc2        0.675    0.316    2.138      Inf    0.033    0.056
    ##     MM_EpiDoc1       -0.309    0.222   -1.393      Inf    0.164   -0.743
    ##   Dep_EpiDoc2 ~                                                         
    ##     MM_EpiDoc1        0.114    0.090    1.265  344.420    0.207   -0.063
    ##   MM_EpiDoc4 ~                                                          
    ##     Dep_EpiDoc2       0.054    0.036    1.501  378.444    0.134   -0.017
    ##     Dep_EpiDoc1       0.020    0.019    1.068      Inf    0.286   -0.017
    ##   MM_EpiDoc2 ~                                                          
    ##     Dep_EpiDoc1       0.007    0.031    0.234      Inf    0.815   -0.054
    ##   Dep_EpiDoc4 ~                                                         
    ##     Age               0.000                                        0.000
    ##   Dep_EpiDoc2 ~                                                         
    ##     Age               0.031    0.019    1.606  429.226    0.109   -0.007
    ##   Dep_EpiDoc1 ~                                                         
    ##     Age               0.055    0.021    2.695      Inf    0.007    0.015
    ##   MM_EpiDoc4 ~                                                          
    ##     Age               0.000                                        0.000
    ##   MM_EpiDoc2 ~                                                          
    ##     Age              -0.012    0.010   -1.218      Inf    0.223   -0.031
    ##   MM_EpiDoc1 ~                                                          
    ##     Age               0.026    0.011    2.275      Inf    0.023    0.004
    ##   Dep_EpiDoc4 ~                                                         
    ##     Education         0.000                                        0.000
    ##   Dep_EpiDoc2 ~                                                         
    ##     Education         0.378    0.218    1.736  331.944    0.083   -0.050
    ##   Dep_EpiDoc1 ~                                                         
    ##     Education         0.777    0.218    3.556      Inf    0.000    0.349
    ##   MM_EpiDoc4 ~                                                          
    ##     Education         0.000                                        0.000
    ##   MM_EpiDoc2 ~                                                          
    ##     Education        -0.099    0.172   -0.578      Inf    0.563   -0.437
    ##   MM_EpiDoc1 ~                                                          
    ##     Education        -0.385    0.227   -1.699      Inf    0.089   -0.830
    ##   MM_EpiDoc4 ~                                                          
    ##     d2_ed            -0.212    0.143   -1.480  336.280    0.140   -0.495
    ##   MM_EpiDoc2 ~                                                          
    ##     d1_ed             0.006    0.140    0.043      Inf    0.966   -0.268
    ##   MM_EpiDoc4 ~                                                          
    ##     age_ed            0.320    0.114    2.796      Inf    0.005    0.095
    ##   MM_EpiDoc2 ~                                                          
    ##     age_ed            0.198    0.294    0.674      Inf    0.500   -0.378
    ##   MM_EpiDoc1 ~                                                          
    ##     age_ed            0.655    0.391    1.675      Inf    0.094   -0.111
    ##  ci.upper   Std.lv  Std.all
    ##                            
    ##     0.501    0.351    0.364
    ##     0.379    0.227    0.233
    ##                            
    ##     0.527    0.360    0.365
    ##                            
    ##     0.488    0.253    0.146
    ##     0.609    0.450    0.375
    ##                            
    ##     0.551    0.359    0.519
    ##                            
    ##     1.294    0.675    0.137
    ##     0.126   -0.309   -0.091
    ##                            
    ##     0.290    0.114    0.032
    ##                            
    ##     0.125    0.054    0.158
    ##     0.057    0.020    0.059
    ##                            
    ##     0.068    0.007    0.037
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.070    0.031    0.099
    ##                            
    ##     0.096    0.055    0.173
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.007   -0.012   -0.190
    ##                            
    ##     0.049    0.026    0.289
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.807    0.378    0.114
    ##                            
    ##     1.205    0.777    0.232
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.238   -0.099   -0.152
    ##                            
    ##     0.059   -0.385   -0.407
    ##                            
    ##     0.070   -0.212   -0.186
    ##                            
    ##     0.280    0.006    0.009
    ##                            
    ##     0.544    0.320    0.231
    ##                            
    ##     0.775    0.198    0.248
    ##                            
    ##     1.420    0.655    0.568
    ## 
    ## Covariances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##  .Dep_EpiDoc1 ~~                                                        
    ##    .MM_EpiDoc1        0.086    0.097    0.893      Inf    0.372   -0.103
    ##  .Dep_EpiDoc2 ~~                                                        
    ##    .MM_EpiDoc2        0.000                                        0.000
    ##  .Dep_EpiDoc4 ~~                                                        
    ##    .MM_EpiDoc4        0.000                                        0.000
    ##   Education ~~                                                          
    ##     d2_ed             0.673    0.110    6.141      Inf    0.000    0.458
    ##  .Dep_EpiDoc2 ~~                                                        
    ##     d2_ed             2.554    0.346    7.373  606.048    0.000    1.874
    ##   d2_ed ~~                                                              
    ##     age_ed            0.551    0.084    6.570  826.019    0.000    0.386
    ##     d1_ed             0.701    0.150    4.680      Inf    0.000    0.407
    ##   Age ~~                                                                
    ##     d2_ed             3.444    0.777    4.432  434.134    0.000    1.916
    ##  .Dep_EpiDoc1 ~~                                                        
    ##     d2_ed             0.990    0.266    3.726      Inf    0.000    0.469
    ##   Education ~~                                                          
    ##     d1_ed             0.719    0.109    6.592      Inf    0.000    0.505
    ##  .Dep_EpiDoc1 ~~                                                        
    ##     d1_ed             2.862    0.371    7.706      Inf    0.000    2.134
    ##   d1_ed ~~                                                              
    ##     age_ed            0.578    0.088    6.606      Inf    0.000    0.407
    ##   Age ~~                                                                
    ##     d1_ed             3.324    0.840    3.959      Inf    0.000    1.678
    ##  .Dep_EpiDoc2 ~~                                                        
    ##     d1_ed            -0.053    0.105   -0.502      Inf    0.616   -0.258
    ##   Education ~~                                                          
    ##     age_ed            0.904    0.078   11.556      Inf    0.000    0.751
    ##   Age ~~                                                                
    ##     age_ed            6.967    0.809    8.614      Inf    0.000    5.382
    ##     Education         3.680    0.796    4.623      Inf    0.000    2.120
    ##  ci.upper   Std.lv  Std.all
    ##                            
    ##     0.275    0.086    0.028
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.887    0.673    0.545
    ##                            
    ##     3.235    2.554    0.704
    ##                            
    ##     0.715    0.551    0.543
    ##     0.994    0.701    0.571
    ##                            
    ##     4.971    3.444    0.267
    ##                            
    ##     1.512    0.990    0.253
    ##                            
    ##     0.933    0.719    0.584
    ##                            
    ##     3.590    2.862    0.734
    ##                            
    ##     0.750    0.578    0.572
    ##                            
    ##     4.969    3.324    0.258
    ##                            
    ##     0.153   -0.053   -0.015
    ##                            
    ##     1.058    0.904    0.888
    ##                            
    ##     8.552    6.967    0.654
    ##     5.241    3.680    0.284
    ## 
    ## Intercepts:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##    .Dep_EpiDoc4       1.293    0.254    5.098      Inf    0.000    0.796
    ##    .Dep_EpiDoc2      -0.243    0.777   -0.313  208.463    0.754   -1.775
    ##    .MM_EpiDoc4        0.569    0.123    4.637      Inf    0.000    0.328
    ##    .MM_EpiDoc2        0.740    0.838    0.883      Inf    0.377   -0.902
    ##    .Dep_EpiDoc1      -0.568    0.819   -0.693      Inf    0.488   -2.173
    ##    .MM_EpiDoc1        0.610    1.055    0.578      Inf    0.563   -1.459
    ##     Age              42.481    0.759   55.998      Inf    0.000   40.994
    ##     Education         2.409    0.072   33.389      Inf    0.000    2.267
    ##     d2_ed             0.094    0.073    1.286      Inf    0.199   -0.049
    ##     d1_ed             0.154    0.074    2.079      Inf    0.038    0.009
    ##     age_ed           -0.096    0.059   -1.637      Inf    0.102   -0.211
    ##  ci.upper   Std.lv  Std.all
    ##     1.790    1.293    0.360
    ##     1.288   -0.243   -0.066
    ##     0.810    0.569    0.450
    ##     2.383    0.740    1.013
    ##     1.037   -0.568   -0.152
    ##     2.678    0.610    0.578
    ##    43.968   42.481    3.646
    ##     2.550    2.409    2.163
    ##     0.238    0.094    0.085
    ##     0.298    0.154    0.139
    ##     0.019   -0.096   -0.105
    ## 
    ## Variances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##    .Dep_EpiDoc4       9.370    1.085    8.633      Inf    0.000    7.241
    ##    .Dep_EpiDoc2      10.697    1.377    7.768  336.106    0.000    7.988
    ##    .MM_EpiDoc4        1.042    0.088   11.829      Inf    0.000    0.869
    ##    .MM_EpiDoc2        0.399    0.063    6.325      Inf    0.000    0.275
    ##    .Dep_EpiDoc1      12.435    1.726    7.204      Inf    0.000    9.052
    ##    .MM_EpiDoc1        0.768    0.121    6.374      Inf    0.000    0.532
    ##     Age             135.726   10.223   13.277      Inf    0.000  115.689
    ##     Education         1.240    0.102   12.178      Inf    0.000    1.041
    ##     d2_ed             1.230    0.207    5.929      Inf    0.000    0.823
    ##     d1_ed             1.223    0.185    6.614      Inf    0.000    0.861
    ##     age_ed            0.836    0.080   10.440      Inf    0.000    0.679
    ##  ci.upper   Std.lv  Std.all
    ##    11.499    9.370    0.726
    ##    13.405   10.697    0.787
    ##     1.214    1.042    0.650
    ##     0.523    0.399    0.747
    ##    15.818   12.435    0.894
    ##     1.005    0.768    0.691
    ##   155.762  135.726    1.000
    ##     1.440    1.240    1.000
    ##     1.637    1.230    1.000
    ##     1.586    1.223    1.000
    ##     0.993    0.836    1.000
    ## 
    ## R-Square:
    ##                    Estimate
    ##     Dep_EpiDoc4       0.274
    ##     Dep_EpiDoc2       0.213
    ##     MM_EpiDoc4        0.350
    ##     MM_EpiDoc2        0.253
    ##     Dep_EpiDoc1       0.106
    ##     MM_EpiDoc1        0.309
    ##     Age               0.000
    ##     Education         0.000
    ##     d2_ed             0.000
    ##     d1_ed             0.000
    ##     age_ed            0.000

``` r
source('./code/parTableToCSV.R')
res<-parTableToCSV(res_mod,path=('./Results/LaggedSEM_mi20/Moderation'))
```

    ## Warning in vcov_lavaan_mi(object, scale.W = scale.W, omit.imps = omit.imps): Could not invert within-imputation covariance matrix. Generalized inverse used instead.
    ## It may be safer to set `scale.W = FALSE' (and `asymptotic = TRUE').

    ## Warning: lavaan->lav_model_vcov():  
    ##    The variance-covariance matrix of the estimated parameters (vcov) does not 
    ##    appear to be positive definite! The smallest eigenvalue (= -1.048120e-14) 
    ##    is smaller than zero. This may be a symptom that the model is not 
    ##    identified.

    ## Warning in vcov_lavaan_mi(object, omit.imps = omit.imps): Could not invert within-imputation covariance matrix. Generalized inverse used instead.
    ## It may be safer to set `scale.W = FALSE' (and `asymptotic = TRUE').

    ## Warning: lavaan->lav_model_vcov():  
    ##    The variance-covariance matrix of the estimated parameters (vcov) does not 
    ##    appear to be positive definite! The smallest eigenvalue (= -1.048120e-14) 
    ##    is smaller than zero. This may be a symptom that the model is not 
    ##    identified.

    ## Warning in vcov_lavaan_mi(object, omit.imps = omit.imps): Could not invert within-imputation covariance matrix. Generalized inverse used instead.
    ## It may be safer to set `scale.W = FALSE' (and `asymptotic = TRUE').

    ## Joining with `by = join_by(group)`

``` r
res$gPartable
```

    ##                      label2 estimate_p_Male estimate_p_Female          CI_Male
    ## 1   Dep_EpiDoc4~Dep_EpiDoc1    0.17 (0.036)         0.351 (0)    [0.011:0.329]
    ## 2   Dep_EpiDoc4~Dep_EpiDoc2   0.247 (0.003)     0.227 (0.003)     [0.085:0.41]
    ## 3   Dep_EpiDoc2~Dep_EpiDoc1       0.425 (0)          0.36 (0)    [0.284:0.566]
    ## 4     MM_EpiDoc4~MM_EpiDoc2   0.409 (0.018)     0.253 (0.035)    [0.071:0.748]
    ## 5     MM_EpiDoc4~MM_EpiDoc1       0.623 (0)          0.45 (0)      [0.447:0.8]
    ## 6     MM_EpiDoc2~MM_EpiDoc1    0.101 (0.03)         0.359 (0)     [0.01:0.192]
    ## 7    Dep_EpiDoc4~MM_EpiDoc2   0.494 (0.092)     0.675 (0.033)   [-0.081:1.069]
    ## 8    Dep_EpiDoc4~MM_EpiDoc1   0.077 (0.705)    -0.309 (0.164)   [-0.323:0.477]
    ## 9    Dep_EpiDoc2~MM_EpiDoc1   0.003 (0.959)     0.114 (0.207)   [-0.127:0.134]
    ## 10   MM_EpiDoc4~Dep_EpiDoc2  -0.046 (0.514)     0.054 (0.134)   [-0.186:0.093]
    ## 11   MM_EpiDoc4~Dep_EpiDoc1   0.006 (0.841)      0.02 (0.286)    [-0.05:0.062]
    ## 12   MM_EpiDoc2~Dep_EpiDoc1    0.01 (0.796)     0.007 (0.815)   [-0.065:0.084]
    ## 13  Dep_EpiDoc1~~MM_EpiDoc1  -0.013 (0.829)     0.086 (0.372)   [-0.126:0.101]
    ## 14  Dep_EpiDoc2~~MM_EpiDoc2          0 (NA)            0 (NA)            [0:0]
    ## 15  Dep_EpiDoc4~~MM_EpiDoc4          0 (NA)            0 (NA)            [0:0]
    ## 16          Dep_EpiDoc4~Age          0 (NA)            0 (NA)            [0:0]
    ## 17          Dep_EpiDoc2~Age   0.008 (0.567)     0.031 (0.109)   [-0.021:0.037]
    ## 18          Dep_EpiDoc1~Age   0.041 (0.033)     0.055 (0.007)    [0.003:0.079]
    ## 19           MM_EpiDoc4~Age          0 (NA)            0 (NA)            [0:0]
    ## 20           MM_EpiDoc2~Age    0.016 (0.05)    -0.012 (0.223)        [0:0.032]
    ## 21           MM_EpiDoc1~Age   0.032 (0.002)     0.026 (0.023)    [0.012:0.053]
    ## 22    Dep_EpiDoc4~Education          0 (NA)            0 (NA)            [0:0]
    ## 23    Dep_EpiDoc2~Education   0.529 (0.001)     0.378 (0.083)    [0.232:0.827]
    ## 24    Dep_EpiDoc1~Education   0.495 (0.003)         0.777 (0)    [0.164:0.827]
    ## 25     MM_EpiDoc4~Education          0 (NA)            0 (NA)            [0:0]
    ## 26     MM_EpiDoc2~Education   0.162 (0.069)    -0.099 (0.563)   [-0.013:0.337]
    ## 27     MM_EpiDoc1~Education    0.047 (0.77)    -0.385 (0.089)    [-0.267:0.36]
    ## 28         MM_EpiDoc4~d2_ed    0.26 (0.333)     -0.212 (0.14)   [-0.267:0.788]
    ## 29         MM_EpiDoc2~d1_ed   0.037 (0.813)     0.006 (0.966)   [-0.271:0.345]
    ## 30        MM_EpiDoc4~age_ed   0.105 (0.285)      0.32 (0.005)   [-0.087:0.297]
    ## 31        MM_EpiDoc2~age_ed  -0.343 (0.041)       0.198 (0.5)  [-0.672:-0.013]
    ## 32        MM_EpiDoc1~age_ed  -0.025 (0.926)     0.655 (0.094)   [-0.548:0.498]
    ## 33         Education~~d2_ed       0.551 (0)         0.673 (0)    [0.375:0.728]
    ## 34       Dep_EpiDoc2~~d2_ed       1.841 (0)         2.554 (0)    [1.294:2.387]
    ## 35            d2_ed~~age_ed       0.546 (0)         0.551 (0)    [0.336:0.756]
    ## 36             d2_ed~~d1_ed       0.594 (0)         0.701 (0)    [0.324:0.864]
    ## 37               Age~~d2_ed       4.365 (0)         3.444 (0)    [1.956:6.774]
    ## 38       Dep_EpiDoc1~~d2_ed       1.062 (0)          0.99 (0)     [0.474:1.65]
    ## 39         Education~~d1_ed        0.48 (0)         0.719 (0)    [0.336:0.625]
    ## 40       Dep_EpiDoc1~~d1_ed       2.487 (0)         2.862 (0)    [1.696:3.278]
    ## 41            d1_ed~~age_ed       0.515 (0)         0.578 (0)    [0.353:0.676]
    ## 42               Age~~d1_ed       4.754 (0)         3.324 (0)    [2.671:6.838]
    ## 43       Dep_EpiDoc2~~d1_ed   0.015 (0.728)    -0.053 (0.616)   [-0.071:0.102]
    ## 44        Education~~age_ed       1.013 (0)         0.904 (0)    [0.859:1.167]
    ## 45              Age~~age_ed      12.895 (0)         6.967 (0)    [9.45:16.339]
    ## 46           Age~~Education       6.794 (0)          3.68 (0)    [4.802:8.786]
    ## 47 Dep_EpiDoc4~~Dep_EpiDoc4        6.91 (0)          9.37 (0)    [4.886:8.935]
    ## 48 Dep_EpiDoc2~~Dep_EpiDoc2       6.248 (0)        10.697 (0)    [4.493:8.004]
    ## 49   MM_EpiDoc4~~MM_EpiDoc4       0.994 (0)         1.042 (0)    [0.782:1.207]
    ## 50   MM_EpiDoc2~~MM_EpiDoc2       0.182 (0)         0.399 (0)    [0.082:0.282]
    ## 51 Dep_EpiDoc1~~Dep_EpiDoc1       8.447 (0)        12.435 (0)   [6.079:10.815]
    ## 52   MM_EpiDoc1~~MM_EpiDoc1         0.7 (0)         0.768 (0)    [0.527:0.874]
    ## 53                 Age~~Age     204.022 (0)       135.726 (0) [158.17:249.873]
    ## 54     Education~~Education       1.159 (0)          1.24 (0)    [1.017:1.302]
    ## 55             d2_ed~~d2_ed       1.025 (0)          1.23 (0)    [0.667:1.383]
    ## 56             d1_ed~~d1_ed       1.021 (0)         1.223 (0)    [0.673:1.368]
    ## 57           age_ed~~age_ed       1.223 (0)         0.836 (0)    [0.957:1.489]
    ## 58            Dep_EpiDoc4~1       1.461 (0)         1.293 (0)     [1.03:1.891]
    ## 59            Dep_EpiDoc2~1   -0.23 (0.707)    -0.243 (0.754)     [-1.44:0.98]
    ## 60             MM_EpiDoc4~1   0.584 (0.005)         0.569 (0)    [0.179:0.989]
    ## 61             MM_EpiDoc2~1  -1.066 (0.071)      0.74 (0.377)   [-2.225:0.092]
    ## 62            Dep_EpiDoc1~1  -0.575 (0.437)    -0.568 (0.488)   [-2.025:0.875]
    ## 63             MM_EpiDoc1~1  -0.901 (0.254)      0.61 (0.563)   [-2.449:0.646]
    ## 64                    Age~1      44.188 (0)        42.481 (0)  [42.217:46.159]
    ## 65              Education~1       2.823 (0)         2.409 (0)    [2.697:2.949]
    ## 66                  d2_ed~1   0.059 (0.377)     0.094 (0.199)    [-0.072:0.19]
    ## 67                  d1_ed~1   0.001 (0.989)     0.154 (0.038)   [-0.125:0.127]
    ## 68                 age_ed~1       0.267 (0)    -0.096 (0.102)    [0.118:0.416]
    ##            CI_Female
    ## 1        [0.2:0.501]
    ## 2      [0.076:0.379]
    ## 3      [0.194:0.527]
    ## 4      [0.018:0.488]
    ## 5      [0.292:0.609]
    ## 6      [0.168:0.551]
    ## 7      [0.056:1.294]
    ## 8     [-0.743:0.126]
    ## 9      [-0.063:0.29]
    ## 10    [-0.017:0.125]
    ## 11    [-0.017:0.057]
    ## 12    [-0.054:0.068]
    ## 13    [-0.103:0.275]
    ## 14             [0:0]
    ## 15             [0:0]
    ## 16             [0:0]
    ## 17     [-0.007:0.07]
    ## 18     [0.015:0.096]
    ## 19             [0:0]
    ## 20    [-0.031:0.007]
    ## 21     [0.004:0.049]
    ## 22             [0:0]
    ## 23     [-0.05:0.807]
    ## 24     [0.349:1.205]
    ## 25             [0:0]
    ## 26    [-0.437:0.238]
    ## 27     [-0.83:0.059]
    ## 28     [-0.495:0.07]
    ## 29     [-0.268:0.28]
    ## 30     [0.095:0.544]
    ## 31    [-0.378:0.775]
    ## 32     [-0.111:1.42]
    ## 33     [0.458:0.887]
    ## 34     [1.874:3.235]
    ## 35     [0.386:0.715]
    ## 36     [0.407:0.994]
    ## 37     [1.916:4.971]
    ## 38     [0.469:1.512]
    ## 39     [0.505:0.933]
    ## 40      [2.134:3.59]
    ## 41      [0.407:0.75]
    ## 42     [1.678:4.969]
    ## 43    [-0.258:0.153]
    ## 44     [0.751:1.058]
    ## 45     [5.382:8.552]
    ## 46      [2.12:5.241]
    ## 47    [7.241:11.499]
    ## 48    [7.988:13.405]
    ## 49     [0.869:1.214]
    ## 50     [0.275:0.523]
    ## 51    [9.052:15.818]
    ## 52     [0.532:1.005]
    ## 53 [115.689:155.762]
    ## 54      [1.041:1.44]
    ## 55     [0.823:1.637]
    ## 56     [0.861:1.586]
    ## 57     [0.679:0.993]
    ## 58      [0.796:1.79]
    ## 59    [-1.775:1.288]
    ## 60      [0.328:0.81]
    ## 61    [-0.902:2.383]
    ## 62    [-2.173:1.037]
    ## 63    [-1.459:2.678]
    ## 64   [40.994:43.968]
    ## 65      [2.267:2.55]
    ## 66    [-0.049:0.238]
    ## 67     [0.009:0.298]
    ## 68    [-0.211:0.019]

# testing moderation with quant

``` r
model_mod<-'
# evolução
Dep_EpiDoc4~c(dd14a,dd14b)*Dep_EpiDoc1+c(dd24a,dd24b)*Dep_EpiDoc2
Dep_EpiDoc2~c(dd12a,dd12b)*Dep_EpiDoc1

MM_EpiDoc4~c(mm24a,mm24b)*MM_EpiDoc2+c(mm14a,mm14b)*MM_EpiDoc1
MM_EpiDoc2~c(mm12a,mm12b)*MM_EpiDoc1 

# causalidade - c/ temporalidade/lag
Dep_EpiDoc4~c(md24a,md24b)*MM_EpiDoc2+c(md14a,md14b)*MM_EpiDoc1
Dep_EpiDoc2~c(md12a,md12b)*MM_EpiDoc1

MM_EpiDoc4~c(dm24a,dm24b)*Dep_EpiDoc2+c(dm14a,dm14b)*Dep_EpiDoc1
MM_EpiDoc2~c(dm12a,dm12b)*Dep_EpiDoc1


# corrrelação mesmo t
MM_EpiDoc1~~Dep_EpiDoc1
MM_EpiDoc2~~0*Dep_EpiDoc2
MM_EpiDoc4~~0*Dep_EpiDoc4

# causalidade - idade
Dep_EpiDoc4~0*Age
Dep_EpiDoc2~ad2*Age
Dep_EpiDoc1~ad1*Age

MM_EpiDoc4~c(am4a,am4b)*Age
MM_EpiDoc2~c(am2a,am2b)*Age
MM_EpiDoc1~c(am1a,am1b)*Age


# educação
Dep_EpiDoc4~0*Education
Dep_EpiDoc2~c(ed2a,ed2b)*Education
Dep_EpiDoc1~c(ed1a,ed1b)*Education
MM_EpiDoc4~0*Education
MM_EpiDoc2~0*Education
MM_EpiDoc1~0*Education



# effects on depression =============
##indirect effect age on dep4
ada:=ad1*dd14a+ad2*dd24a+ad1*dd12a*dd24a
adb:=ad1*dd14b+ad2*dd24b+ad1*dd12b*dd24b
amda:=am1a*md14a+am2a*md24a+am1a*mm12a*md24a+ad1*dm12a*md24a
amdb:=am1b*md14b+am2b*md24b+am1b*mm12b*md24b+ad1*dm12b*md24b


adta:=ada+amda
adtb:=adb+amdb

##indirect effect MM1 on dep4
mda:=md12a*dd24a+mm12a*md24a+md14a
mdb:=md12b*dd24b+mm12b*md24b+md14b

## cumulative effect d1 on dep4
dda:=dm12a*md24a+dd12a*dd24a+dd14a
ddb:=dm12b*md24b+dd12b*dd24b+dd14b

## indirect effect of educaiton on dep4
eda:=ed1a*dd12a*dd24a+ed2a*dd24a
edb:=ed1b*dd12b*dd24b+ed2b*dd24b
emda:=ed1a*dm12a*md24a
emdb:=ed1b*dm12b*md24b

edta:=eda+amda
edtb:=edb+emdb


# interaction terms =========
Dep_EpiDoc4~MM2_ed+MM1_ed
Dep_EpiDoc2~MM1_ed



MM_EpiDoc4~d2_ed+d2_ed
MM_EpiDoc2~d1_ed



MM_EpiDoc4~age_ed
MM_EpiDoc2~age_ed
MM_EpiDoc1~age_ed

# correlation among interaction terms and its original var
MM2_ed~~Education+MM_EpiDoc2+age_ed+Age+MM1_ed+MM_EpiDoc1+d2_ed+d1_ed
MM1_ed~~Education+MM_EpiDoc1+age_ed+Age+MM_EpiDoc2+d2_ed+d1_ed
d2_ed~~Education+Dep_EpiDoc2+age_ed+d1_ed+Age+Dep_EpiDoc1
d1_ed~~Education+Dep_EpiDoc1+age_ed+Age+Dep_EpiDoc2
age_ed~~Education+Age
Education~~Age



'


res_mod<-lavaan.mi::sem.mi(data=data.mi,
         model=model_mod,
         estimator = "ML", test = "yuan.bentler.mplus",
         likelihood = "wishart",
           sampling.weights='ipw',group='Sex',
         orthogonal=FALSE
         )
```

    ## Warning: lavaan->lavParTable():  
    ##    using a single label per parameter in a multiple group setting implies 
    ##    imposing equality constraints across all the groups; If this is not 
    ##    intended, either remove the label(s), or use a vector of labels (one for 
    ##    each group); See the Multiple groups section in the man page of 
    ##    model.syntax.

    ## Warning: lavaan->lav_model_vcov():  
    ##    The variance-covariance matrix of the estimated parameters (vcov) does not 
    ##    appear to be positive definite! The smallest eigenvalue (= -6.239058e-14) 
    ##    is smaller than zero. This may be a symptom that the model is not 
    ##    identified.
    ## Warning: lavaan->lav_model_vcov():  
    ##    The variance-covariance matrix of the estimated parameters (vcov) does not 
    ##    appear to be positive definite! The smallest eigenvalue (= -6.239058e-14) 
    ##    is smaller than zero. This may be a symptom that the model is not 
    ##    identified.

    ## Warning: lavaan->lav_model_vcov():  
    ##    The variance-covariance matrix of the estimated parameters (vcov) does not 
    ##    appear to be positive definite! The smallest eigenvalue (= -1.087566e-14) 
    ##    is smaller than zero. This may be a symptom that the model is not 
    ##    identified.

    ## Warning: lavaan->lav_model_vcov():  
    ##    The variance-covariance matrix of the estimated parameters (vcov) does not 
    ##    appear to be positive definite! The smallest eigenvalue (= -4.646157e-14) 
    ##    is smaller than zero. This may be a symptom that the model is not 
    ##    identified.

    ## Warning: lavaan->lav_model_vcov():  
    ##    The variance-covariance matrix of the estimated parameters (vcov) does not 
    ##    appear to be positive definite! The smallest eigenvalue (= -2.288876e-15) 
    ##    is smaller than zero. This may be a symptom that the model is not 
    ##    identified.

    ## Warning: lavaan->lav_model_vcov():  
    ##    The variance-covariance matrix of the estimated parameters (vcov) does not 
    ##    appear to be positive definite! The smallest eigenvalue (= -7.805493e-15) 
    ##    is smaller than zero. This may be a symptom that the model is not 
    ##    identified.

    ## Warning: lavaan->lav_model_vcov():  
    ##    The variance-covariance matrix of the estimated parameters (vcov) does not 
    ##    appear to be positive definite! The smallest eigenvalue (= -9.456429e-14) 
    ##    is smaller than zero. This may be a symptom that the model is not 
    ##    identified.

    ## Warning: lavaan->lav_model_vcov():  
    ##    The variance-covariance matrix of the estimated parameters (vcov) does not 
    ##    appear to be positive definite! The smallest eigenvalue (= -3.281801e-14) 
    ##    is smaller than zero. This may be a symptom that the model is not 
    ##    identified.

    ## Warning: lavaan->lav_model_vcov():  
    ##    The variance-covariance matrix of the estimated parameters (vcov) does not 
    ##    appear to be positive definite! The smallest eigenvalue (= -7.016358e-14) 
    ##    is smaller than zero. This may be a symptom that the model is not 
    ##    identified.

    ## Warning: lavaan->lav_model_vcov():  
    ##    The variance-covariance matrix of the estimated parameters (vcov) does not 
    ##    appear to be positive definite! The smallest eigenvalue (= -7.172195e-15) 
    ##    is smaller than zero. This may be a symptom that the model is not 
    ##    identified.

    ## Warning: lavaan->lav_model_vcov():  
    ##    The variance-covariance matrix of the estimated parameters (vcov) does not 
    ##    appear to be positive definite! The smallest eigenvalue (= -1.290474e-14) 
    ##    is smaller than zero. This may be a symptom that the model is not 
    ##    identified.

    ## Warning: lavaan->lav_model_vcov():  
    ##    The variance-covariance matrix of the estimated parameters (vcov) does not 
    ##    appear to be positive definite! The smallest eigenvalue (= -4.367246e-14) 
    ##    is smaller than zero. This may be a symptom that the model is not 
    ##    identified.

    ## Warning: lavaan->lav_model_vcov():  
    ##    The variance-covariance matrix of the estimated parameters (vcov) does not 
    ##    appear to be positive definite! The smallest eigenvalue (= -1.042837e-14) 
    ##    is smaller than zero. This may be a symptom that the model is not 
    ##    identified.

    ## Warning: lavaan->lav_model_vcov():  
    ##    The variance-covariance matrix of the estimated parameters (vcov) does not 
    ##    appear to be positive definite! The smallest eigenvalue (= -8.136477e-14) 
    ##    is smaller than zero. This may be a symptom that the model is not 
    ##    identified.

    ## Warning: lavaan->lav_model_vcov():  
    ##    The variance-covariance matrix of the estimated parameters (vcov) does not 
    ##    appear to be positive definite! The smallest eigenvalue (= -2.370710e-14) 
    ##    is smaller than zero. This may be a symptom that the model is not 
    ##    identified.

    ## Warning: lavaan->lav_model_vcov():  
    ##    The variance-covariance matrix of the estimated parameters (vcov) does not 
    ##    appear to be positive definite! The smallest eigenvalue (= -4.551505e-15) 
    ##    is smaller than zero. This may be a symptom that the model is not 
    ##    identified.

    ## Warning: lavaan->lav_model_vcov():  
    ##    The variance-covariance matrix of the estimated parameters (vcov) does not 
    ##    appear to be positive definite! The smallest eigenvalue (= -1.445364e-14) 
    ##    is smaller than zero. This may be a symptom that the model is not 
    ##    identified.

    ## Warning: lavaan->lav_model_vcov():  
    ##    The variance-covariance matrix of the estimated parameters (vcov) does not 
    ##    appear to be positive definite! The smallest eigenvalue (= -1.764425e-14) 
    ##    is smaller than zero. This may be a symptom that the model is not 
    ##    identified.

    ## Warning: lavaan->lav_model_vcov():  
    ##    The variance-covariance matrix of the estimated parameters (vcov) does not 
    ##    appear to be positive definite! The smallest eigenvalue (= -1.448614e-14) 
    ##    is smaller than zero. This may be a symptom that the model is not 
    ##    identified.

    ## Warning: lavaan->lav_model_vcov():  
    ##    The variance-covariance matrix of the estimated parameters (vcov) does not 
    ##    appear to be positive definite! The smallest eigenvalue (= -2.846134e-14) 
    ##    is smaller than zero. This may be a symptom that the model is not 
    ##    identified.

    ## Warning: lavaan->lav_model_vcov():  
    ##    The variance-covariance matrix of the estimated parameters (vcov) does not 
    ##    appear to be positive definite! The smallest eigenvalue (= -4.510984e-14) 
    ##    is smaller than zero. This may be a symptom that the model is not 
    ##    identified.

``` r
summary(res_mod,fit.meas=TRUE,standardized=TRUE,rsquare=TRUE,ci=TRUE)
```

    ## Warning: lavaan->lav_model_vcov():  
    ##    The variance-covariance matrix of the estimated parameters (vcov) does not 
    ##    appear to be positive definite! The smallest eigenvalue (= -6.239058e-14) 
    ##    is smaller than zero. This may be a symptom that the model is not 
    ##    identified.

    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.

    ## Warning: lavaan->lav_model_vcov():  
    ##    The variance-covariance matrix of the estimated parameters (vcov) does not 
    ##    appear to be positive definite! The smallest eigenvalue (= -6.239058e-14) 
    ##    is smaller than zero. This may be a symptom that the model is not 
    ##    identified.
    ## Warning: lavaan->lav_model_vcov():  
    ##    The variance-covariance matrix of the estimated parameters (vcov) does not 
    ##    appear to be positive definite! The smallest eigenvalue (= -6.239058e-14) 
    ##    is smaller than zero. This may be a symptom that the model is not 
    ##    identified.
    ## Warning: lavaan->lav_model_vcov():  
    ##    The variance-covariance matrix of the estimated parameters (vcov) does not 
    ##    appear to be positive definite! The smallest eigenvalue (= -6.239058e-14) 
    ##    is smaller than zero. This may be a symptom that the model is not 
    ##    identified.

    ## lavaan.mi object fit to 20 imputed data sets using:
    ##  - lavaan    (0.6-19)
    ##  - lavaan.mi (0.1-0)
    ## See class?lavaan.mi help page for available methods. 
    ## 
    ## Convergence information:
    ## The model converged on 20 imputed data sets.
    ## Standard errors were available for all imputations.
    ## 
    ##   Estimator                                         ML
    ##   Optimization method                           NLMINB
    ##   Number of model parameters                       166
    ##   Number of equality constraints                     2
    ## 
    ##   Number of observations per group:                   
    ##     Male                                           355
    ##     Female                                         317
    ##   Sampling weights variable                        ipw
    ## 
    ## Model Test User Model:
    ## 
    ##                                                     Standard      Scaled
    ##   Test statistic                                      95.659      78.513
    ##   Degrees of freedom                                      44          44
    ##   P-value                                              0.000       0.001
    ##   Average scaling correction factor                                1.218
    ##   Pooling method                                          D4            
    ##     Pooled statistic                              "standard"            
    ##     "yuan.bentler.mplus" correction applied            AFTER     pooling
    ## 
    ## Model Test Baseline Model:
    ## 
    ##   Test statistic                              8387.828    4647.954
    ##   Degrees of freedom                               156         156
    ##   P-value                                        0.000       0.000
    ##   Scaling correction factor                                  1.805
    ## 
    ## User Model versus Baseline Model:
    ## 
    ##   Comparative Fit Index (CFI)                    0.994       0.992
    ##   Tucker-Lewis Index (TLI)                       0.978       0.973
    ##                                                                   
    ##   Robust Comparative Fit Index (CFI)                         0.995
    ##   Robust Tucker-Lewis Index (TLI)                            0.982
    ## 
    ## Loglikelihood and Information Criteria:
    ## 
    ##   Loglikelihood user model (H0)             -11239.356  -11239.356
    ##   Scaling correction factor                                  2.392
    ##       for the MLR correction                                      
    ##   Loglikelihood unrestricted model (H1)     -11187.670  -11187.670
    ##   Scaling correction factor                                  2.162
    ##       for the MLR correction                                      
    ##                                                                   
    ##   Akaike (AIC)                               22806.712   22806.712
    ##   Bayesian (BIC)                             23546.394   23546.394
    ##   Sample-size adjusted Bayesian (SABIC)      23025.681   23025.681
    ## 
    ## Root Mean Square Error of Approximation:
    ## 
    ##   RMSEA                                          0.059       0.049
    ##   90 Percent confidence interval - lower         0.043       0.033
    ##   90 Percent confidence interval - upper         0.075       0.065
    ##   P-value H_0: RMSEA <= 0.050                    0.165       0.510
    ##   P-value H_0: RMSEA >= 0.080                    0.016       0.000
    ##                                                                   
    ##   Robust RMSEA                                               0.053
    ##   90 Percent confidence interval - lower                     0.034
    ##   90 Percent confidence interval - upper                     0.072
    ##   P-value H_0: Robust RMSEA <= 0.050                         0.363
    ##   P-value H_0: Robust RMSEA >= 0.080                         0.009
    ## 
    ## Standardized Root Mean Square Residual:
    ## 
    ##   SRMR                                           0.067       0.067
    ## 
    ## Parameter Estimates:
    ## 
    ##   Standard errors                                    Sandwich
    ##   Information bread                                  Observed
    ##   Observed information based on                       Hessian
    ##                                                              
    ##   Pooled across imputations              Rubin's (1987) rules
    ##   Augment within-imputation variance     Scale by average RIV
    ##   Wald test for pooled parameters          t(df) distribution
    ## 
    ##   Pooled t statistics with df >= 1000 are displayed with
    ##   df = Inf(inity) to save space. Although the t distribution
    ##   with large df closely approximates a standard normal
    ##   distribution, exact df for reporting these t tests can be
    ##   obtained from parameterEstimates.mi() 
    ## 
    ## 
    ## 
    ## Group 1 [Male]:
    ## 
    ## Regressions:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##   Dep_EpiDoc4 ~                                                         
    ##     Dp_EpD1 (dd14)    0.175    0.082    2.128  614.521    0.034    0.013
    ##     Dp_EpD2 (dd24)    0.255    0.082    3.101  339.667    0.002    0.093
    ##   Dep_EpiDoc2 ~                                                         
    ##     Dp_EpD1 (dd12)    0.420    0.071    5.896  324.512    0.000    0.280
    ##   MM_EpiDoc4 ~                                                          
    ##     MM_EpD2 (mm24)    0.345    0.153    2.249      Inf    0.025    0.044
    ##     MM_EpD1 (mm14)    0.549    0.092    5.934      Inf    0.000    0.367
    ##   MM_EpiDoc2 ~                                                          
    ##     MM_EpD1 (mm12)    0.120    0.050    2.384      Inf    0.017    0.021
    ##   Dep_EpiDoc4 ~                                                         
    ##     MM_EpD2 (md24)   -0.170    0.674   -0.252      Inf    0.801   -1.493
    ##     MM_EpD1 (md14)    0.656    0.584    1.123      Inf    0.262   -0.489
    ##   Dep_EpiDoc2 ~                                                         
    ##     MM_EpD1 (md12)   -0.266    0.586   -0.454      Inf    0.650   -1.416
    ##   MM_EpiDoc4 ~                                                          
    ##     Dp_EpD2 (dm24)   -0.156    0.073   -2.135      Inf    0.033   -0.300
    ##     Dp_EpD1 (dm14)    0.006    0.029    0.212      Inf    0.832   -0.050
    ##   MM_EpiDoc2 ~                                                          
    ##     Dp_EpD1 (dm12)   -0.014    0.051   -0.283      Inf    0.777   -0.114
    ##   Dep_EpiDoc4 ~                                                         
    ##     Age               0.000                                        0.000
    ##   Dep_EpiDoc2 ~                                                         
    ##     Age      (ad2)    0.014    0.012    1.186  115.143    0.238   -0.010
    ##   Dep_EpiDoc1 ~                                                         
    ##     Age      (ad1)    0.045    0.015    3.042      Inf    0.002    0.016
    ##   MM_EpiDoc4 ~                                                          
    ##     Age     (am4a)    0.033    0.009    3.906      Inf    0.000    0.017
    ##   MM_EpiDoc2 ~                                                          
    ##     Age     (am2a)    0.007    0.004    1.645      Inf    0.100   -0.001
    ##   MM_EpiDoc1 ~                                                          
    ##     Age     (am1a)    0.029    0.007    4.243      Inf    0.000    0.016
    ##   Dep_EpiDoc4 ~                                                         
    ##     Educatn           0.000                                        0.000
    ##   Dep_EpiDoc2 ~                                                         
    ##     Educatn (ed2a)    0.455    0.155    2.927  384.885    0.004    0.149
    ##   Dep_EpiDoc1 ~                                                         
    ##     Educatn (ed1a)    0.472    0.161    2.923      Inf    0.003    0.155
    ##   MM_EpiDoc4 ~                                                          
    ##     Educatn           0.000                                        0.000
    ##   MM_EpiDoc2 ~                                                          
    ##     Educatn           0.000                                        0.000
    ##   MM_EpiDoc1 ~                                                          
    ##     Educatn           0.000                                        0.000
    ##   Dep_EpiDoc4 ~                                                         
    ##     MM2_ed            0.399    0.378    1.055      Inf    0.292   -0.343
    ##     MM1_ed           -0.573    0.546   -1.050  609.393    0.294   -1.645
    ##   Dep_EpiDoc2 ~                                                         
    ##     MM1_ed            0.269    0.545    0.493      Inf    0.622   -0.800
    ##   MM_EpiDoc4 ~                                                          
    ##     d2_ed             0.661    0.277    2.387  603.340    0.017    0.117
    ##   MM_EpiDoc2 ~                                                          
    ##     d1_ed             0.045    0.164    0.272      Inf    0.786   -0.277
    ##   MM_EpiDoc4 ~                                                          
    ##     age_ed           -0.300    0.134   -2.245  848.648    0.025   -0.563
    ##   MM_EpiDoc2 ~                                                          
    ##     age_ed           -0.098    0.058   -1.693      Inf    0.090   -0.212
    ##   MM_EpiDoc1 ~                                                          
    ##     age_ed            0.042    0.091    0.465      Inf    0.642   -0.136
    ##  ci.upper   Std.lv  Std.all
    ##                            
    ##     0.336    0.175    0.187
    ##     0.417    0.255    0.265
    ##                            
    ##     0.560    0.420    0.433
    ##                            
    ##     0.645    0.345    0.124
    ##     0.730    0.549    0.416
    ##                            
    ##     0.219    0.120    0.254
    ##                            
    ##     1.153   -0.170   -0.027
    ##     1.801    0.656    0.219
    ##                            
    ##     0.884   -0.266   -0.085
    ##                            
    ##    -0.013   -0.156   -0.369
    ##     0.062    0.006    0.015
    ##                            
    ##     0.085   -0.014   -0.097
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.038    0.014    0.069
    ##                            
    ##     0.074    0.045    0.209
    ##                            
    ##     0.050    0.033    0.378
    ##                            
    ##     0.015    0.007    0.210
    ##                            
    ##     0.043    0.029    0.439
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.761    0.455    0.165
    ##                            
    ##     0.788    0.472    0.166
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     1.141    0.399    0.101
    ##     0.499   -0.573   -0.205
    ##                            
    ##     1.338    0.269    0.093
    ##                            
    ##     1.205    0.661    0.531
    ##                            
    ##     0.366    0.045    0.099
    ##                            
    ##    -0.038   -0.300   -0.263
    ##                            
    ##     0.016   -0.098   -0.240
    ##                            
    ##     0.220    0.042    0.049
    ## 
    ## Covariances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##  .Dep_EpiDoc1 ~~                                                        
    ##    .MM_EpiDoc1       -0.047    0.057   -0.830      Inf    0.407   -0.159
    ##  .Dep_EpiDoc2 ~~                                                        
    ##    .MM_EpiDoc2        0.000                                        0.000
    ##  .Dep_EpiDoc4 ~~                                                        
    ##    .MM_EpiDoc4        0.000                                        0.000
    ##   Education ~~                                                          
    ##     MM2_ed            0.029    0.045    0.643      Inf    0.520   -0.059
    ##  .MM_EpiDoc2 ~~                                                         
    ##     MM2_ed            0.285    0.095    3.011      Inf    0.003    0.099
    ##   MM2_ed ~~                                                             
    ##     age_ed            0.077    0.051    1.517      Inf    0.129   -0.023
    ##   Age ~~                                                                
    ##     MM2_ed            1.259    0.641    1.964      Inf    0.050    0.002
    ##   MM2_ed ~~                                                             
    ##     MM1_ed            0.170    0.089    1.905      Inf    0.057   -0.005
    ##  .MM_EpiDoc1 ~~                                                         
    ##     MM2_ed            0.130    0.061    2.118      Inf    0.034    0.010
    ##   MM2_ed ~~                                                             
    ##     d2_ed             0.032    0.032    1.021      Inf    0.307   -0.030
    ##     d1_ed             0.045    0.037    1.199      Inf    0.230   -0.028
    ##   Education ~~                                                          
    ##     MM1_ed            0.429    0.088    4.888      Inf    0.000    0.257
    ##  .MM_EpiDoc1 ~~                                                         
    ##     MM1_ed            0.688    0.101    6.799      Inf    0.000    0.490
    ##   MM1_ed ~~                                                             
    ##     age_ed            0.601    0.100    6.020      Inf    0.000    0.405
    ##   Age ~~                                                                
    ##     MM1_ed            7.187    1.122    6.407      Inf    0.000    4.988
    ##  .MM_EpiDoc2 ~~                                                         
    ##     MM1_ed           -0.009    0.009   -1.011      Inf    0.312   -0.028
    ##   MM1_ed ~~                                                             
    ##     d2_ed             0.288    0.086    3.362      Inf    0.001    0.120
    ##     d1_ed             0.282    0.078    3.606      Inf    0.000    0.129
    ##   Education ~~                                                          
    ##     d2_ed             0.551    0.089    6.222      Inf    0.000    0.378
    ##  .Dep_EpiDoc2 ~~                                                        
    ##     d2_ed             1.832    0.273    6.698  159.237    0.000    1.292
    ##   d2_ed ~~                                                              
    ##     age_ed            0.562    0.104    5.391  619.505    0.000    0.357
    ##     d1_ed             0.588    0.133    4.425      Inf    0.000    0.327
    ##   Age ~~                                                                
    ##     d2_ed             4.687    1.108    4.229  252.885    0.000    2.504
    ##  .Dep_EpiDoc1 ~~                                                        
    ##     d2_ed             1.040    0.290    3.585  987.528    0.000    0.471
    ##   Education ~~                                                          
    ##     d1_ed             0.480    0.073    6.594      Inf    0.000    0.338
    ##  .Dep_EpiDoc1 ~~                                                        
    ##     d1_ed             2.465    0.382    6.453      Inf    0.000    1.717
    ##   d1_ed ~~                                                              
    ##     age_ed            0.523    0.077    6.777      Inf    0.000    0.372
    ##   Age ~~                                                                
    ##     d1_ed             4.951    0.902    5.491      Inf    0.000    3.184
    ##  .Dep_EpiDoc2 ~~                                                        
    ##     d1_ed             0.002    0.040    0.041  479.578    0.967   -0.077
    ##   Education ~~                                                          
    ##     age_ed            1.013    0.078   13.057      Inf    0.000    0.861
    ##   Age ~~                                                                
    ##     age_ed           12.895    1.733    7.439      Inf    0.000    9.497
    ##     Education         6.794    1.002    6.779      Inf    0.000    4.830
    ##  ci.upper   Std.lv  Std.all
    ##                            
    ##     0.064   -0.047   -0.019
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.116    0.029    0.037
    ##                            
    ##     0.471    0.285    0.908
    ##                            
    ##     0.177    0.077    0.097
    ##                            
    ##     2.516    1.259    0.122
    ##                            
    ##     0.346    0.170    0.230
    ##                            
    ##     0.250    0.130    0.214
    ##                            
    ##     0.094    0.032    0.044
    ##     0.118    0.045    0.061
    ##                            
    ##     0.601    0.429    0.389
    ##                            
    ##     0.887    0.688    0.800
    ##                            
    ##     0.797    0.601    0.530
    ##                            
    ##     9.386    7.187    0.491
    ##                            
    ##     0.009   -0.009   -0.021
    ##                            
    ##     0.457    0.288    0.278
    ##     0.435    0.282    0.273
    ##                            
    ##     0.725    0.551    0.506
    ##                            
    ##     2.372    1.832    0.724
    ##                            
    ##     0.767    0.562    0.503
    ##     0.848    0.588    0.578
    ##                            
    ##     6.869    4.687    0.324
    ##                            
    ##     1.609    1.040    0.354
    ##                            
    ##     0.623    0.480    0.444
    ##                            
    ##     3.214    2.465    0.843
    ##                            
    ##     0.674    0.523    0.470
    ##                            
    ##     6.718    4.951    0.345
    ##                            
    ##     0.081    0.002    0.001
    ##                            
    ##     1.165    1.013    0.851
    ##                            
    ##    16.292   12.895    0.816
    ##     8.758    6.794    0.442
    ## 
    ## Intercepts:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##    .Dep_EpiDoc4       1.178    0.449    2.624  854.096    0.009    0.297
    ##    .Dep_EpiDoc2      -0.100    0.659   -0.152  313.730    0.879   -1.396
    ##    .MM_EpiDoc4       -0.448    0.300   -1.495      Inf    0.135   -1.036
    ##    .MM_EpiDoc2       -0.206    0.212   -0.971      Inf    0.332   -0.621
    ##    .Dep_EpiDoc1      -0.688    0.633   -1.086      Inf    0.277   -1.930
    ##    .MM_EpiDoc1       -0.668    0.281   -2.380      Inf    0.017   -1.219
    ##     Age              44.188    0.992   44.547      Inf    0.000   42.244
    ##     Education         2.823    0.064   44.415      Inf    0.000    2.698
    ##     MM2_ed           -0.124    0.040   -3.143      Inf    0.002   -0.202
    ##     MM1_ed            0.011    0.065    0.162      Inf    0.872   -0.117
    ##     d2_ed             0.059    0.066    0.897      Inf    0.370   -0.070
    ##     d1_ed             0.001    0.063    0.014      Inf    0.989   -0.123
    ##     age_ed            0.267    0.075    3.558      Inf    0.000    0.120
    ##  ci.upper   Std.lv  Std.all
    ##     2.060    1.178    0.412
    ##     1.196   -0.100   -0.034
    ##     0.139   -0.448   -0.356
    ##     0.210   -0.206   -0.453
    ##     0.553   -0.688   -0.224
    ##    -0.118   -0.668   -0.699
    ##    46.132   44.188    3.094
    ##     2.947    2.823    2.622
    ##    -0.047   -0.124   -0.172
    ##     0.138    0.011    0.010
    ##     0.188    0.059    0.058
    ##     0.125    0.001    0.001
    ##     0.414    0.267    0.242
    ## 
    ## Variances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##    .Dep_EpiDoc4       6.871    1.014    6.777      Inf    0.000    4.882
    ##    .Dep_EpiDoc2       6.251    0.874    7.152  110.017    0.000    4.519
    ##    .MM_EpiDoc4        0.939    0.102    9.188      Inf    0.000    0.739
    ##    .MM_EpiDoc2        0.189    0.053    3.578      Inf    0.000    0.085
    ##    .Dep_EpiDoc1       8.462    1.192    7.097      Inf    0.000    6.125
    ##    .MM_EpiDoc1        0.704    0.089    7.895      Inf    0.000    0.530
    ##     Age             204.022   23.075    8.842      Inf    0.000  158.795
    ##     Education         1.159    0.072   16.175      Inf    0.000    1.019
    ##     MM2_ed            0.522    0.219    2.381      Inf    0.017    0.092
    ##     MM1_ed            1.051    0.190    5.543      Inf    0.000    0.679
    ##     d2_ed             1.023    0.182    5.631  882.750    0.000    0.666
    ##     d1_ed             1.011    0.167    6.070      Inf    0.000    0.685
    ##     age_ed            1.223    0.134    9.122      Inf    0.000    0.960
    ##  ci.upper   Std.lv  Std.all
    ##     8.860    6.871    0.839
    ##     7.983    6.251    0.706
    ##     1.140    0.939    0.591
    ##     0.292    0.189    0.916
    ##    10.799    8.462    0.898
    ##     0.879    0.704    0.769
    ##   249.248  204.022    1.000
    ##     1.300    1.159    1.000
    ##     0.952    0.522    1.000
    ##     1.423    1.051    1.000
    ##     1.379    1.023    1.000
    ##     1.338    1.011    1.000
    ##     1.486    1.223    1.000
    ## 
    ## R-Square:
    ##                    Estimate
    ##     Dep_EpiDoc4       0.161
    ##     Dep_EpiDoc2       0.294
    ##     MM_EpiDoc4        0.409
    ##     MM_EpiDoc2        0.084
    ##     Dep_EpiDoc1       0.102
    ##     MM_EpiDoc1        0.231
    ##     Age               0.000
    ##     Education         0.000
    ##     MM2_ed            0.000
    ##     MM1_ed            0.000
    ##     d2_ed             0.000
    ##     d1_ed             0.000
    ##     age_ed            0.000
    ## 
    ## 
    ## Group 2 [Female]:
    ## 
    ## Regressions:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##   Dep_EpiDoc4 ~                                                         
    ##     Dp_EpD1 (dd14)    0.347    0.074    4.690      Inf    0.000    0.202
    ##     Dp_EpD2 (dd24)    0.226    0.076    2.985  959.237    0.003    0.077
    ##   Dep_EpiDoc2 ~                                                         
    ##     Dp_EpD1 (dd12)    0.365    0.083    4.394  604.104    0.000    0.202
    ##   MM_EpiDoc4 ~                                                          
    ##     MM_EpD2 (mm24)    0.282    0.115    2.454      Inf    0.014    0.057
    ##     MM_EpD1 (mm14)    0.348    0.089    3.888      Inf    0.000    0.172
    ##   MM_EpiDoc2 ~                                                          
    ##     MM_EpD1 (mm12)    0.373    0.097    3.861      Inf    0.000    0.184
    ##   Dep_EpiDoc4 ~                                                         
    ##     MM_EpD2 (md24)    0.125    0.936    0.133  375.383    0.894   -1.715
    ##     MM_EpD1 (md14)   -0.677    0.501   -1.352  404.150    0.177   -1.661
    ##   Dep_EpiDoc2 ~                                                         
    ##     MM_EpD1 (md12)    0.489    0.676    0.722  162.507    0.471   -0.847
    ##   MM_EpiDoc4 ~                                                          
    ##     Dp_EpD2 (dm24)    0.004    0.036    0.103  243.620    0.918   -0.067
    ##     Dp_EpD1 (dm14)    0.026    0.018    1.461      Inf    0.144   -0.009
    ##   MM_EpiDoc2 ~                                                          
    ##     Dp_EpD1 (dm12)    0.001    0.035    0.022      Inf    0.983   -0.068
    ##   Dep_EpiDoc4 ~                                                         
    ##     Age               0.000                                        0.000
    ##   Dep_EpiDoc2 ~                                                         
    ##     Age      (ad2)    0.014    0.012    1.186  115.143    0.238   -0.010
    ##   Dep_EpiDoc1 ~                                                         
    ##     Age      (ad1)    0.045    0.015    3.042      Inf    0.002    0.016
    ##   MM_EpiDoc4 ~                                                          
    ##     Age     (am4b)    0.028    0.008    3.666      Inf    0.000    0.013
    ##   MM_EpiDoc2 ~                                                          
    ##     Age     (am2b)   -0.008    0.005   -1.504      Inf    0.132   -0.018
    ##   MM_EpiDoc1 ~                                                          
    ##     Age     (am1b)    0.045    0.005    8.647      Inf    0.000    0.035
    ##   Dep_EpiDoc4 ~                                                         
    ##     Educatn           0.000                                        0.000
    ##   Dep_EpiDoc2 ~                                                         
    ##     Educatn (ed2b)    0.507    0.235    2.159  523.092    0.031    0.046
    ##   Dep_EpiDoc1 ~                                                         
    ##     Educatn (ed1b)    0.809    0.206    3.923      Inf    0.000    0.405
    ##   MM_EpiDoc4 ~                                                          
    ##     Educatn           0.000                                        0.000
    ##   MM_EpiDoc2 ~                                                          
    ##     Educatn           0.000                                        0.000
    ##   MM_EpiDoc1 ~                                                          
    ##     Educatn           0.000                                        0.000
    ##   Dep_EpiDoc4 ~                                                         
    ##     MM2_ed            0.361    0.579    0.623  245.209    0.534   -0.780
    ##     MM1_ed            0.397    0.494    0.803  170.292    0.423   -0.579
    ##   Dep_EpiDoc2 ~                                                         
    ##     MM1_ed           -0.420    0.685   -0.613  150.532    0.541   -1.774
    ##   MM_EpiDoc4 ~                                                          
    ##     d2_ed            -0.015    0.145   -0.105  249.842    0.916   -0.300
    ##   MM_EpiDoc2 ~                                                          
    ##     d1_ed            -0.003    0.140   -0.022      Inf    0.982   -0.277
    ##   MM_EpiDoc4 ~                                                          
    ##     age_ed            0.047    0.130    0.363      Inf    0.717   -0.207
    ##   MM_EpiDoc2 ~                                                          
    ##     age_ed            0.063    0.083    0.756      Inf    0.450   -0.100
    ##   MM_EpiDoc1 ~                                                          
    ##     age_ed            0.078    0.090    0.872      Inf    0.383   -0.098
    ##  ci.upper   Std.lv  Std.all
    ##                            
    ##     0.491    0.347    0.361
    ##     0.374    0.226    0.231
    ##                            
    ##     0.528    0.365    0.371
    ##                            
    ##     0.508    0.282    0.164
    ##     0.523    0.348    0.290
    ##                            
    ##     0.562    0.373    0.535
    ##                            
    ##     1.965    0.125    0.026
    ##     0.307   -0.677   -0.199
    ##                            
    ##     1.825    0.489    0.140
    ##                            
    ##     0.074    0.004    0.011
    ##     0.062    0.026    0.078
    ##                            
    ##     0.069    0.001    0.004
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.038    0.014    0.046
    ##                            
    ##     0.074    0.045    0.140
    ##                            
    ##     0.044    0.028    0.262
    ##                            
    ##     0.002   -0.008   -0.121
    ##                            
    ##     0.055    0.045    0.500
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.969    0.507    0.154
    ##                            
    ##     1.213    0.809    0.241
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     1.501    0.361    0.119
    ##     1.373    0.397    0.115
    ##                            
    ##     0.934   -0.420   -0.119
    ##                            
    ##     0.270   -0.015   -0.013
    ##                            
    ##     0.270   -0.003   -0.005
    ##                            
    ##     0.301    0.047    0.034
    ##                            
    ##     0.226    0.063    0.078
    ##                            
    ##     0.255    0.078    0.068
    ## 
    ## Covariances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##  .Dep_EpiDoc1 ~~                                                        
    ##    .MM_EpiDoc1        0.036    0.091    0.393      Inf    0.695   -0.142
    ##  .Dep_EpiDoc2 ~~                                                        
    ##    .MM_EpiDoc2        0.000                                        0.000
    ##  .Dep_EpiDoc4 ~~                                                        
    ##    .MM_EpiDoc4        0.000                                        0.000
    ##   Education ~~                                                          
    ##     MM2_ed            0.308    0.102    3.006      Inf    0.003    0.107
    ##  .MM_EpiDoc2 ~~                                                         
    ##     MM2_ed            0.598    0.105    5.675      Inf    0.000    0.391
    ##   MM2_ed ~~                                                             
    ##     age_ed            0.329    0.103    3.194      Inf    0.001    0.127
    ##   Age ~~                                                                
    ##     MM2_ed            3.390    1.038    3.266      Inf    0.001    1.356
    ##   MM2_ed ~~                                                             
    ##     MM1_ed            0.638    0.203    3.137      Inf    0.002    0.239
    ##  .MM_EpiDoc1 ~~                                                         
    ##     MM2_ed            0.427    0.139    3.062      Inf    0.002    0.154
    ##   MM2_ed ~~                                                             
    ##     d2_ed             0.193    0.072    2.703      Inf    0.007    0.053
    ##     d1_ed             0.210    0.076    2.759      Inf    0.006    0.061
    ##   Education ~~                                                          
    ##     MM1_ed            0.458    0.110    4.179      Inf    0.000    0.243
    ##  .MM_EpiDoc1 ~~                                                         
    ##     MM1_ed            0.686    0.128    5.347      Inf    0.000    0.434
    ##   MM1_ed ~~                                                             
    ##     age_ed            0.538    0.114    4.714      Inf    0.000    0.314
    ##   Age ~~                                                                
    ##     MM1_ed            6.205    1.057    5.868      Inf    0.000    4.132
    ##  .MM_EpiDoc2 ~~                                                         
    ##     MM1_ed           -0.006    0.018   -0.342      Inf    0.733   -0.042
    ##   MM1_ed ~~                                                             
    ##     d2_ed             0.294    0.086    3.429      Inf    0.001    0.126
    ##     d1_ed             0.331    0.102    3.254      Inf    0.001    0.132
    ##   Education ~~                                                          
    ##     d2_ed             0.673    0.108    6.226      Inf    0.000    0.461
    ##  .Dep_EpiDoc2 ~~                                                        
    ##     d2_ed             2.564    0.344    7.458  705.285    0.000    1.889
    ##   d2_ed ~~                                                              
    ##     age_ed            0.528    0.081    6.488      Inf    0.000    0.368
    ##     d1_ed             0.697    0.148    4.716      Inf    0.000    0.407
    ##   Age ~~                                                                
    ##     d2_ed             2.887    0.653    4.421  498.365    0.000    1.604
    ##  .Dep_EpiDoc1 ~~                                                        
    ##     d2_ed             1.013    0.262    3.871  889.770    0.000    0.499
    ##   Education ~~                                                          
    ##     d1_ed             0.719    0.108    6.683      Inf    0.000    0.508
    ##  .Dep_EpiDoc1 ~~                                                        
    ##     d1_ed             2.879    0.371    7.754      Inf    0.000    2.151
    ##   d1_ed ~~                                                              
    ##     age_ed            0.569    0.085    6.668      Inf    0.000    0.402
    ##   Age ~~                                                                
    ##     d1_ed             3.036    0.735    4.132      Inf    0.000    1.596
    ##  .Dep_EpiDoc2 ~~                                                        
    ##     d1_ed            -0.037    0.102   -0.360      Inf    0.719   -0.236
    ##   Education ~~                                                          
    ##     age_ed            0.904    0.077   11.715      Inf    0.000    0.753
    ##   Age ~~                                                                
    ##     age_ed            6.967    0.798    8.733      Inf    0.000    5.403
    ##     Education         3.680    0.785    4.687      Inf    0.000    2.141
    ##  ci.upper   Std.lv  Std.all
    ##                            
    ##     0.213    0.036    0.011
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.509    0.308    0.234
    ##                            
    ##     0.805    0.598    0.801
    ##                            
    ##     0.531    0.329    0.305
    ##                            
    ##     5.424    3.390    0.246
    ##                            
    ##     1.037    0.638    0.518
    ##                            
    ##     0.700    0.427    0.410
    ##                            
    ##     0.334    0.193    0.148
    ##     0.359    0.210    0.161
    ##                            
    ##     0.672    0.458    0.394
    ##                            
    ##     0.937    0.686    0.746
    ##                            
    ##     0.762    0.538    0.564
    ##                            
    ##     8.277    6.205    0.511
    ##                            
    ##     0.030   -0.006   -0.010
    ##                            
    ##     0.461    0.294    0.254
    ##     0.531    0.331    0.288
    ##                            
    ##     0.884    0.673    0.546
    ##                            
    ##     3.240    2.564    0.708
    ##                            
    ##     0.688    0.528    0.522
    ##     0.987    0.697    0.571
    ##                            
    ##     4.171    2.887    0.224
    ##                            
    ##     1.526    1.013    0.258
    ##                            
    ##     0.930    0.719    0.585
    ##                            
    ##     3.607    2.879    0.736
    ##                            
    ##     0.736    0.569    0.564
    ##                            
    ##     4.476    3.036    0.236
    ##                            
    ##     0.163   -0.037   -0.010
    ##                            
    ##     1.056    0.904    0.888
    ##                            
    ##     8.530    6.967    0.654
    ##     5.220    3.680    0.284
    ## 
    ## Intercepts:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##    .Dep_EpiDoc4       1.695    0.421    4.027  667.165    0.000    0.869
    ##    .Dep_EpiDoc2      -0.135    0.889   -0.152  239.699    0.880   -1.886
    ##    .MM_EpiDoc4       -0.468    0.302   -1.551      Inf    0.121   -1.059
    ##    .MM_EpiDoc2        0.320    0.183    1.753      Inf    0.080   -0.038
    ##    .Dep_EpiDoc1      -0.203    0.699   -0.291      Inf    0.771   -1.574
    ##    .MM_EpiDoc1       -1.185    0.222   -5.349      Inf    0.000   -1.619
    ##     Age              42.481    0.748   56.772      Inf    0.000   41.015
    ##     Education         2.409    0.071   33.850      Inf    0.000    2.269
    ##     MM2_ed            0.105    0.071    1.477      Inf    0.140   -0.034
    ##     MM1_ed           -0.019    0.068   -0.282      Inf    0.778   -0.152
    ##     d2_ed             0.094    0.072    1.304      Inf    0.192   -0.047
    ##     d1_ed             0.154    0.073    2.107      Inf    0.035    0.011
    ##     age_ed           -0.096    0.058   -1.659      Inf    0.097   -0.210
    ##  ci.upper   Std.lv  Std.all
    ##     2.522    1.695    0.473
    ##     1.616   -0.135   -0.037
    ##     0.123   -0.468   -0.370
    ##     0.678    0.320    0.436
    ##     1.167   -0.203   -0.054
    ##    -0.751   -1.185   -1.125
    ##    43.948   42.481    3.646
    ##     2.548    2.409    2.163
    ##     0.245    0.105    0.089
    ##     0.113   -0.019   -0.018
    ##     0.236    0.094    0.085
    ##     0.296    0.154    0.139
    ##     0.017   -0.096   -0.105
    ## 
    ## Variances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##    .Dep_EpiDoc4       9.282    1.058    8.770      Inf    0.000    7.206
    ##    .Dep_EpiDoc2      10.705    1.363    7.853  365.201    0.000    8.024
    ##    .MM_EpiDoc4        0.997    0.086   11.663      Inf    0.000    0.830
    ##    .MM_EpiDoc2        0.400    0.062    6.465      Inf    0.000    0.279
    ##    .Dep_EpiDoc1      12.586    1.721    7.312      Inf    0.000    9.212
    ##    .MM_EpiDoc1        0.777    0.121    6.400      Inf    0.000    0.539
    ##     Age             135.726   10.084   13.460      Inf    0.000  115.962
    ##     Education         1.240    0.100   12.346      Inf    0.000    1.043
    ##     MM2_ed            1.394    0.360    3.877      Inf    0.000    0.689
    ##     MM1_ed            1.087    0.242    4.495      Inf    0.000    0.613
    ##     d2_ed             1.224    0.205    5.970      Inf    0.000    0.822
    ##     d1_ed             1.217    0.182    6.700      Inf    0.000    0.861
    ##     age_ed            0.836    0.079   10.584      Inf    0.000    0.681
    ##  ci.upper   Std.lv  Std.all
    ##    11.358    9.282    0.721
    ##    13.385   10.705    0.795
    ##     1.165    0.997    0.625
    ##     0.522    0.400    0.743
    ##    15.959   12.586    0.903
    ##     1.015    0.777    0.700
    ##   155.489  135.726    1.000
    ##     1.437    1.240    1.000
    ##     2.099    1.394    1.000
    ##     1.561    1.087    1.000
    ##     1.626    1.224    1.000
    ##     1.573    1.217    1.000
    ##     0.991    0.836    1.000
    ## 
    ## R-Square:
    ##                    Estimate
    ##     Dep_EpiDoc4       0.279
    ##     Dep_EpiDoc2       0.205
    ##     MM_EpiDoc4        0.375
    ##     MM_EpiDoc2        0.257
    ##     Dep_EpiDoc1       0.097
    ##     MM_EpiDoc1        0.300
    ##     Age               0.000
    ##     Education         0.000
    ##     MM2_ed            0.000
    ##     MM1_ed            0.000
    ##     d2_ed             0.000
    ##     d1_ed             0.000
    ##     age_ed            0.000
    ## 
    ## Defined Parameters:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##     ada               0.016    0.006    2.610  671.686    0.009    0.004
    ##     adb               0.023    0.007    3.163      Inf    0.002    0.009
    ##     amda              0.018    0.017    1.063      Inf    0.288   -0.015
    ##     amdb             -0.029    0.021   -1.376      Inf    0.169   -0.071
    ##     adta              0.034    0.019    1.785      Inf    0.074   -0.003
    ##     adtb             -0.007    0.022   -0.313      Inf    0.755   -0.050
    ##     mda               0.569    0.616    0.923      Inf    0.356   -0.640
    ##     mdb              -0.518    0.509   -1.018      Inf    0.309   -1.515
    ##     dda               0.284    0.089    3.188      Inf    0.001    0.109
    ##     ddb               0.428    0.068    6.285      Inf    0.000    0.295
    ##     eda               0.167    0.070    2.394  319.012    0.017    0.030
    ##     edb               0.180    0.084    2.149      Inf    0.032    0.016
    ##     emda              0.001    0.011    0.109      Inf    0.913   -0.020
    ##     emdb             -0.001    0.003   -0.194      Inf    0.846   -0.006
    ##     edta              0.185    0.074    2.494  328.657    0.013    0.039
    ##     edtb              0.180    0.084    2.141      Inf    0.032    0.015
    ##  ci.upper   Std.lv  Std.all
    ##     0.029    0.016    0.082
    ##     0.037    0.023    0.109
    ##     0.050    0.018    0.088
    ##     0.013   -0.029   -0.096
    ##     0.071    0.034    0.170
    ##     0.036   -0.007    0.014
    ##     1.777    0.567    0.190
    ##     0.479   -0.520   -0.153
    ##     0.459    0.284    0.305
    ##     0.562    0.429    0.446
    ##     0.304    0.167    0.063
    ##     0.345    0.181    0.056
    ##     0.022    0.001    0.000
    ##     0.005    0.000    0.000
    ##     0.330    0.184    0.151
    ##     0.344    0.181    0.056

``` r
source('./code/parTableToCSV.R')
res<-parTableToCSV(res_mod,path=('./Results/LaggedSEM_mi20/Moderation'))
```

    ## Warning: lavaan->lav_model_vcov():  
    ##    The variance-covariance matrix of the estimated parameters (vcov) does not 
    ##    appear to be positive definite! The smallest eigenvalue (= -6.239058e-14) 
    ##    is smaller than zero. This may be a symptom that the model is not 
    ##    identified.
    ## Warning: lavaan->lav_model_vcov():  
    ##    The variance-covariance matrix of the estimated parameters (vcov) does not 
    ##    appear to be positive definite! The smallest eigenvalue (= -6.239058e-14) 
    ##    is smaller than zero. This may be a symptom that the model is not 
    ##    identified.

    ## Joining with `by = join_by(group)`

    ## Warning: There was 1 warning in `mutate()`.
    ## ℹ In argument: `Group2 = as.numeric(substr(label, 2, 2))`.
    ## Caused by warning:
    ## ! NAs introduced by coercion

    ##    Group2 label     estimate_p             CI Group
    ## 1      NA     a  0.016 (0.009)  [0.004:0.029]  <NA>
    ## 2      NA     b  0.023 (0.002)  [0.009:0.037]  <NA>
    ## 3      NA    da  0.018 (0.288)  [-0.015:0.05]  <NA>
    ## 4      NA    db -0.029 (0.169) [-0.071:0.013]  <NA>
    ## 5      NA    ta  0.034 (0.074) [-0.003:0.071]  <NA>
    ## 6      NA    tb -0.007 (0.755)  [-0.05:0.036]  <NA>
    ## 7      NA     a  0.569 (0.356)  [-0.64:1.777]  <NA>
    ## 8      NA     b -0.518 (0.309) [-1.515:0.479]  <NA>
    ## 9      NA     a  0.284 (0.001)  [0.109:0.459]  <NA>
    ## 10     NA     b      0.428 (0)  [0.295:0.562]  <NA>
    ## 11     NA     a  0.167 (0.017)   [0.03:0.304]  <NA>
    ## 12     NA     b   0.18 (0.032)  [0.016:0.345]  <NA>
    ## 13     NA    da  0.001 (0.913)  [-0.02:0.022]  <NA>
    ## 14     NA    db -0.001 (0.846) [-0.006:0.005]  <NA>
    ## 15     NA    ta  0.185 (0.013)   [0.039:0.33]  <NA>
    ## 16     NA    tb   0.18 (0.032)  [0.015:0.344]  <NA>

    ## Warning in reshapeWide(data, idvar = idvar, timevar = timevar, varying =
    ## varying, : there are records with missing times, which will be dropped.

    ## Warning in reshapeWide(data, idvar = idvar, timevar = timevar, varying =
    ## varying, : multiple rows match for Group=NA: first taken

``` r
res$gPartable
```

    ##                      label2 estimate_p_Male estimate_p_Female           CI_Male
    ## 1   Dep_EpiDoc4~Dep_EpiDoc1   0.175 (0.034)         0.347 (0)     [0.013:0.336]
    ## 2   Dep_EpiDoc4~Dep_EpiDoc2   0.255 (0.002)     0.226 (0.003)     [0.093:0.417]
    ## 3   Dep_EpiDoc2~Dep_EpiDoc1        0.42 (0)         0.365 (0)       [0.28:0.56]
    ## 4     MM_EpiDoc4~MM_EpiDoc2   0.345 (0.025)     0.282 (0.014)     [0.044:0.645]
    ## 5     MM_EpiDoc4~MM_EpiDoc1       0.549 (0)         0.348 (0)      [0.367:0.73]
    ## 6     MM_EpiDoc2~MM_EpiDoc1    0.12 (0.017)         0.373 (0)     [0.021:0.219]
    ## 7    Dep_EpiDoc4~MM_EpiDoc2   -0.17 (0.801)     0.125 (0.894)    [-1.493:1.153]
    ## 8    Dep_EpiDoc4~MM_EpiDoc1   0.656 (0.262)    -0.677 (0.177)    [-0.489:1.801]
    ## 9    Dep_EpiDoc2~MM_EpiDoc1   -0.266 (0.65)     0.489 (0.471)    [-1.416:0.884]
    ## 10   MM_EpiDoc4~Dep_EpiDoc2  -0.156 (0.033)     0.004 (0.918)     [-0.3:-0.013]
    ## 11   MM_EpiDoc4~Dep_EpiDoc1   0.006 (0.832)     0.026 (0.144)     [-0.05:0.062]
    ## 12   MM_EpiDoc2~Dep_EpiDoc1  -0.014 (0.777)     0.001 (0.983)    [-0.114:0.085]
    ## 13  Dep_EpiDoc1~~MM_EpiDoc1  -0.047 (0.407)     0.036 (0.695)    [-0.159:0.064]
    ## 14  Dep_EpiDoc2~~MM_EpiDoc2          0 (NA)            0 (NA)             [0:0]
    ## 15  Dep_EpiDoc4~~MM_EpiDoc4          0 (NA)            0 (NA)             [0:0]
    ## 16          Dep_EpiDoc4~Age          0 (NA)            0 (NA)             [0:0]
    ## 17          Dep_EpiDoc2~Age   0.014 (0.238)     0.014 (0.238)     [-0.01:0.038]
    ## 18          Dep_EpiDoc1~Age   0.045 (0.002)     0.045 (0.002)     [0.016:0.074]
    ## 19           MM_EpiDoc4~Age       0.033 (0)         0.028 (0)      [0.017:0.05]
    ## 20           MM_EpiDoc2~Age     0.007 (0.1)    -0.008 (0.132)    [-0.001:0.015]
    ## 21           MM_EpiDoc1~Age       0.029 (0)         0.045 (0)     [0.016:0.043]
    ## 22    Dep_EpiDoc4~Education          0 (NA)            0 (NA)             [0:0]
    ## 23    Dep_EpiDoc2~Education   0.455 (0.004)     0.507 (0.031)     [0.149:0.761]
    ## 24    Dep_EpiDoc1~Education   0.472 (0.003)         0.809 (0)     [0.155:0.788]
    ## 25     MM_EpiDoc4~Education          0 (NA)            0 (NA)             [0:0]
    ## 26     MM_EpiDoc2~Education          0 (NA)            0 (NA)             [0:0]
    ## 27     MM_EpiDoc1~Education          0 (NA)            0 (NA)             [0:0]
    ## 28       Dep_EpiDoc4~MM2_ed   0.399 (0.292)     0.361 (0.534)    [-0.343:1.141]
    ## 29       Dep_EpiDoc4~MM1_ed  -0.573 (0.294)     0.397 (0.423)    [-1.645:0.499]
    ## 30       Dep_EpiDoc2~MM1_ed   0.269 (0.622)     -0.42 (0.541)      [-0.8:1.338]
    ## 31         MM_EpiDoc4~d2_ed   0.661 (0.017)    -0.015 (0.916)     [0.117:1.205]
    ## 32         MM_EpiDoc2~d1_ed   0.045 (0.786)    -0.003 (0.982)    [-0.277:0.366]
    ## 33        MM_EpiDoc4~age_ed    -0.3 (0.025)     0.047 (0.717)   [-0.563:-0.038]
    ## 34        MM_EpiDoc2~age_ed   -0.098 (0.09)      0.063 (0.45)    [-0.212:0.016]
    ## 35        MM_EpiDoc1~age_ed   0.042 (0.642)     0.078 (0.383)     [-0.136:0.22]
    ## 36        Education~~MM2_ed    0.029 (0.52)     0.308 (0.003)    [-0.059:0.116]
    ## 37       MM_EpiDoc2~~MM2_ed   0.285 (0.003)         0.598 (0)     [0.099:0.471]
    ## 38           MM2_ed~~age_ed   0.077 (0.129)     0.329 (0.001)    [-0.023:0.177]
    ## 39              Age~~MM2_ed    1.259 (0.05)      3.39 (0.001)     [0.002:2.516]
    ## 40           MM2_ed~~MM1_ed    0.17 (0.057)     0.638 (0.002)    [-0.005:0.346]
    ## 41       MM_EpiDoc1~~MM2_ed    0.13 (0.034)     0.427 (0.002)       [0.01:0.25]
    ## 42            MM2_ed~~d2_ed   0.032 (0.307)     0.193 (0.007)     [-0.03:0.094]
    ## 43            MM2_ed~~d1_ed    0.045 (0.23)      0.21 (0.006)    [-0.028:0.118]
    ## 44        Education~~MM1_ed       0.429 (0)         0.458 (0)     [0.257:0.601]
    ## 45       MM_EpiDoc1~~MM1_ed       0.688 (0)         0.686 (0)      [0.49:0.887]
    ## 46           MM1_ed~~age_ed       0.601 (0)         0.538 (0)     [0.405:0.797]
    ## 47              Age~~MM1_ed       7.187 (0)         6.205 (0)     [4.988:9.386]
    ## 48       MM_EpiDoc2~~MM1_ed  -0.009 (0.312)    -0.006 (0.733)    [-0.028:0.009]
    ## 49            MM1_ed~~d2_ed   0.288 (0.001)     0.294 (0.001)      [0.12:0.457]
    ## 50            MM1_ed~~d1_ed       0.282 (0)     0.331 (0.001)     [0.129:0.435]
    ## 51         Education~~d2_ed       0.551 (0)         0.673 (0)     [0.378:0.725]
    ## 52       Dep_EpiDoc2~~d2_ed       1.832 (0)         2.564 (0)     [1.292:2.372]
    ## 53            d2_ed~~age_ed       0.562 (0)         0.528 (0)     [0.357:0.767]
    ## 54             d2_ed~~d1_ed       0.588 (0)         0.697 (0)     [0.327:0.848]
    ## 55               Age~~d2_ed       4.687 (0)         2.887 (0)     [2.504:6.869]
    ## 56       Dep_EpiDoc1~~d2_ed        1.04 (0)         1.013 (0)     [0.471:1.609]
    ## 57         Education~~d1_ed        0.48 (0)         0.719 (0)     [0.338:0.623]
    ## 58       Dep_EpiDoc1~~d1_ed       2.465 (0)         2.879 (0)     [1.717:3.214]
    ## 59            d1_ed~~age_ed       0.523 (0)         0.569 (0)     [0.372:0.674]
    ## 60               Age~~d1_ed       4.951 (0)         3.036 (0)     [3.184:6.718]
    ## 61       Dep_EpiDoc2~~d1_ed   0.002 (0.967)    -0.037 (0.719)    [-0.077:0.081]
    ## 62        Education~~age_ed       1.013 (0)         0.904 (0)     [0.861:1.165]
    ## 63              Age~~age_ed      12.895 (0)         6.967 (0)    [9.497:16.292]
    ## 64           Age~~Education       6.794 (0)          3.68 (0)      [4.83:8.758]
    ## 65 Dep_EpiDoc4~~Dep_EpiDoc4       6.871 (0)         9.282 (0)      [4.882:8.86]
    ## 66 Dep_EpiDoc2~~Dep_EpiDoc2       6.251 (0)        10.705 (0)     [4.519:7.983]
    ## 67   MM_EpiDoc4~~MM_EpiDoc4       0.939 (0)         0.997 (0)      [0.739:1.14]
    ## 68   MM_EpiDoc2~~MM_EpiDoc2       0.189 (0)           0.4 (0)     [0.085:0.292]
    ## 69 Dep_EpiDoc1~~Dep_EpiDoc1       8.462 (0)        12.586 (0)    [6.125:10.799]
    ## 70   MM_EpiDoc1~~MM_EpiDoc1       0.704 (0)         0.777 (0)      [0.53:0.879]
    ## 71                 Age~~Age     204.022 (0)       135.726 (0) [158.795:249.248]
    ## 72     Education~~Education       1.159 (0)          1.24 (0)       [1.019:1.3]
    ## 73           MM2_ed~~MM2_ed   0.522 (0.017)         1.394 (0)     [0.092:0.952]
    ## 74           MM1_ed~~MM1_ed       1.051 (0)         1.087 (0)     [0.679:1.423]
    ## 75             d2_ed~~d2_ed       1.023 (0)         1.224 (0)     [0.666:1.379]
    ## 76             d1_ed~~d1_ed       1.011 (0)         1.217 (0)     [0.685:1.338]
    ## 77           age_ed~~age_ed       1.223 (0)         0.836 (0)      [0.96:1.486]
    ## 78            Dep_EpiDoc4~1   1.178 (0.009)         1.695 (0)      [0.297:2.06]
    ## 79            Dep_EpiDoc2~1    -0.1 (0.879)     -0.135 (0.88)    [-1.396:1.196]
    ## 80             MM_EpiDoc4~1  -0.448 (0.135)    -0.468 (0.121)    [-1.036:0.139]
    ## 81             MM_EpiDoc2~1  -0.206 (0.332)       0.32 (0.08)     [-0.621:0.21]
    ## 82            Dep_EpiDoc1~1  -0.688 (0.277)    -0.203 (0.771)     [-1.93:0.553]
    ## 83             MM_EpiDoc1~1  -0.668 (0.017)        -1.185 (0)   [-1.219:-0.118]
    ## 84                    Age~1      44.188 (0)        42.481 (0)   [42.244:46.132]
    ## 85              Education~1       2.823 (0)         2.409 (0)     [2.698:2.947]
    ## 86                 MM2_ed~1  -0.124 (0.002)      0.105 (0.14)   [-0.202:-0.047]
    ## 87                 MM1_ed~1   0.011 (0.872)    -0.019 (0.778)    [-0.117:0.138]
    ## 88                  d2_ed~1    0.059 (0.37)     0.094 (0.192)     [-0.07:0.188]
    ## 89                  d1_ed~1   0.001 (0.989)     0.154 (0.035)    [-0.123:0.125]
    ## 90                 age_ed~1       0.267 (0)    -0.096 (0.097)      [0.12:0.414]
    ##            CI_Female
    ## 1      [0.202:0.491]
    ## 2      [0.077:0.374]
    ## 3      [0.202:0.528]
    ## 4      [0.057:0.508]
    ## 5      [0.172:0.523]
    ## 6      [0.184:0.562]
    ## 7     [-1.715:1.965]
    ## 8     [-1.661:0.307]
    ## 9     [-0.847:1.825]
    ## 10    [-0.067:0.074]
    ## 11    [-0.009:0.062]
    ## 12    [-0.068:0.069]
    ## 13    [-0.142:0.213]
    ## 14             [0:0]
    ## 15             [0:0]
    ## 16             [0:0]
    ## 17     [-0.01:0.038]
    ## 18     [0.016:0.074]
    ## 19     [0.013:0.044]
    ## 20    [-0.018:0.002]
    ## 21     [0.035:0.055]
    ## 22             [0:0]
    ## 23     [0.046:0.969]
    ## 24     [0.405:1.213]
    ## 25             [0:0]
    ## 26             [0:0]
    ## 27             [0:0]
    ## 28     [-0.78:1.501]
    ## 29    [-0.579:1.373]
    ## 30    [-1.774:0.934]
    ## 31       [-0.3:0.27]
    ## 32     [-0.277:0.27]
    ## 33    [-0.207:0.301]
    ## 34      [-0.1:0.226]
    ## 35    [-0.098:0.255]
    ## 36     [0.107:0.509]
    ## 37     [0.391:0.805]
    ## 38     [0.127:0.531]
    ## 39     [1.356:5.424]
    ## 40     [0.239:1.037]
    ## 41       [0.154:0.7]
    ## 42     [0.053:0.334]
    ## 43     [0.061:0.359]
    ## 44     [0.243:0.672]
    ## 45     [0.434:0.937]
    ## 46     [0.314:0.762]
    ## 47     [4.132:8.277]
    ## 48     [-0.042:0.03]
    ## 49     [0.126:0.461]
    ## 50     [0.132:0.531]
    ## 51     [0.461:0.884]
    ## 52      [1.889:3.24]
    ## 53     [0.368:0.688]
    ## 54     [0.407:0.987]
    ## 55     [1.604:4.171]
    ## 56     [0.499:1.526]
    ## 57      [0.508:0.93]
    ## 58     [2.151:3.607]
    ## 59     [0.402:0.736]
    ## 60     [1.596:4.476]
    ## 61    [-0.236:0.163]
    ## 62     [0.753:1.056]
    ## 63      [5.403:8.53]
    ## 64      [2.141:5.22]
    ## 65    [7.206:11.358]
    ## 66    [8.024:13.385]
    ## 67      [0.83:1.165]
    ## 68     [0.279:0.522]
    ## 69    [9.212:15.959]
    ## 70     [0.539:1.015]
    ## 71 [115.962:155.489]
    ## 72     [1.043:1.437]
    ## 73     [0.689:2.099]
    ## 74     [0.613:1.561]
    ## 75     [0.822:1.626]
    ## 76     [0.861:1.573]
    ## 77     [0.681:0.991]
    ## 78     [0.869:2.522]
    ## 79    [-1.886:1.616]
    ## 80    [-1.059:0.123]
    ## 81    [-0.038:0.678]
    ## 82    [-1.574:1.167]
    ## 83   [-1.619:-0.751]
    ## 84   [41.015:43.948]
    ## 85     [2.269:2.548]
    ## 86    [-0.034:0.245]
    ## 87    [-0.152:0.113]
    ## 88    [-0.047:0.236]
    ## 89     [0.011:0.296]
    ## 90     [-0.21:0.017]

## full data

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
    ##  P boot           1.3-31  2024-08-28 [?] RSPM
    ##  P broom          1.0.9   2025-07-28 [?] RSPM
    ##  P cachem         1.1.0   2024-05-16 [?] CRAN (R 4.5.2)
    ##  P cli            3.6.5   2025-04-23 [?] CRAN (R 4.5.2)
    ##  P codetools      0.2-20  2024-03-31 [?] RSPM
    ##  P DBI            1.2.3   2024-06-02 [?] RSPM
    ##  P devtools       2.5.0   2026-03-14 [?] CRAN (R 4.5.2)
    ##  P digest         0.6.37  2024-08-19 [?] RSPM
    ##  P dplyr        * 1.1.4   2023-11-17 [?] CRAN (R 4.5.2)
    ##  P ellipsis       0.3.2   2021-04-29 [?] RSPM
    ##  P evaluate       1.0.4   2025-06-18 [?] RSPM
    ##  P farver         2.1.2   2024-05-13 [?] CRAN (R 4.5.2)
    ##  P fastmap        1.2.0   2024-05-15 [?] CRAN (R 4.5.2)
    ##  P foreach        1.5.2   2022-02-02 [?] RSPM
    ##  P fs             2.0.1   2026-03-24 [?] RSPM
    ##  P generics       0.1.4   2025-05-09 [?] CRAN (R 4.5.2)
    ##  P ggplot2      * 4.0.0   2025-09-11 [?] RSPM
    ##  P glmnet         4.1-10  2025-07-17 [?] RSPM
    ##  P glue           1.8.0   2024-09-30 [?] CRAN (R 4.5.2)
    ##  P gtable         0.3.6   2024-10-25 [?] CRAN (R 4.5.2)
    ##  P htmltools      0.5.8.1 2024-04-04 [?] RSPM
    ##  P iterators      1.0.14  2022-02-05 [?] RSPM
    ##  P jomo           2.7-6   2023-04-15 [?] RSPM
    ##  P knitr          1.50    2025-03-16 [?] RSPM
    ##  P lattice        0.22-6  2024-03-20 [?] RSPM
    ##  P lavaan       * 0.6-19  2024-09-26 [?] CRAN (R 4.5.2)
    ##  P lavaan.mi    * 0.1-0   2025-03-10 [?] RSPM
    ##  P lifecycle      1.0.5   2026-01-08 [?] CRAN (R 4.5.2)
    ##  P lme4           1.1-37  2025-03-26 [?] RSPM
    ##  P magrittr       2.0.3   2022-03-30 [?] RSPM
    ##  P MASS           7.3-61  2024-06-13 [?] RSPM
    ##  P Matrix         1.7-1   2024-10-18 [?] RSPM
    ##  P memoise        2.0.1   2021-11-26 [?] CRAN (R 4.5.2)
    ##  P mice         * 3.18.0  2025-05-27 [?] RSPM
    ##  P miceadds     * 3.17-44 2024-01-09 [?] RSPM
    ##  P minqa          1.2.8   2024-08-17 [?] RSPM
    ##  P mitml          0.4-5   2023-03-08 [?] RSPM
    ##  P mitools        2.4     2019-04-26 [?] RSPM
    ##  P mnormt         2.1.1   2022-09-26 [?] RSPM
    ##  P nlme           3.1-166 2024-08-14 [?] RSPM
    ##  P nloptr         2.2.1   2025-03-17 [?] RSPM
    ##  P nnet           7.3-19  2023-05-03 [?] RSPM
    ##  P pan            1.9     2023-12-07 [?] RSPM
    ##  P pbivnorm       0.6.0   2015-01-23 [?] RSPM
    ##  P pillar         1.11.0  2025-07-04 [?] RSPM
    ##  P pkgbuild       1.4.8   2025-05-26 [?] CRAN (R 4.5.2)
    ##  P pkgconfig      2.0.3   2019-09-22 [?] CRAN (R 4.5.2)
    ##  P pkgload        1.5.1   2026-04-01 [?] CRAN (R 4.5.2)
    ##  P purrr          1.2.1   2026-01-09 [?] CRAN (R 4.5.2)
    ##  P quadprog       1.5-8   2019-11-20 [?] RSPM
    ##  P R6             2.6.1   2025-02-15 [?] CRAN (R 4.5.2)
    ##  P rbibutils      2.4.1   2026-01-21 [?] RSPM
    ##  P RColorBrewer   1.1-3   2022-04-03 [?] CRAN (R 4.5.2)
    ##  P Rcpp           1.1.0   2025-07-02 [?] RSPM
    ##  P Rdpack         2.6.4   2025-04-09 [?] RSPM
    ##  P reformulas     0.4.4   2026-02-02 [?] RSPM
    ##    renv           1.1.5   2025-07-24 [1] RSPM
    ##  P rlang          1.2.0   2026-04-06 [?] RSPM
    ##  P rmarkdown      2.29    2024-11-04 [?] RSPM
    ##  P rpart          4.1.23  2023-12-05 [?] RSPM
    ##  P rstudioapi     0.17.1  2024-10-22 [?] RSPM
    ##  P S7             0.2.0   2024-11-07 [?] RSPM
    ##  P scales         1.4.0   2025-04-24 [?] CRAN (R 4.5.2)
    ##  P sessioninfo    1.2.3   2025-02-05 [?] CRAN (R 4.5.2)
    ##  P shape          1.4.6.1 2024-02-23 [?] RSPM
    ##  P survival       3.7-0   2024-06-05 [?] RSPM
    ##  P tibble         3.3.0   2025-06-08 [?] RSPM
    ##  P tidyr          1.3.1   2024-01-24 [?] RSPM
    ##  P tidyselect     1.2.1   2024-03-11 [?] CRAN (R 4.5.2)
    ##  P usethis        3.2.1   2025-09-06 [?] CRAN (R 4.5.2)
    ##  P vctrs          0.7.2   2026-03-21 [?] CRAN (R 4.5.2)
    ##  P withr          3.0.2   2024-10-28 [?] CRAN (R 4.5.2)
    ##  P xfun           0.53    2025-08-19 [?] CRAN (R 4.5.2)
    ##  P yaml           2.3.10  2024-07-26 [?] RSPM
    ## 
    ##  [1] /mnt/Data/Projects/RutePortugal/renv/library/linux-ubuntu-noble/R-4.5/x86_64-pc-linux-gnu
    ##  [2] /home/tomoe/.cache/R/renv/sandbox/linux-ubuntu-noble/R-4.5/x86_64-pc-linux-gnu/9a444a72
    ## 
    ##  * ── Packages attached to the search path.
    ##  P ── Loaded and on-disk path mismatch.
    ## 
    ## ──────────────────────────────────────────────────────────────────────────────
