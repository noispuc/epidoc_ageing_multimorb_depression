4 Multi-group assessment with Cross-Lagged Structural Equation Modelling
analysis
================
Tomoe Gusberti, Dr. Eng
2025-03-20

load data

``` r
load(file='./Data/Data_Imp20_forSEM.RData')#data.mi
GVAr='Sex'
folderName=paste0('MGroup_',GVAr)
```

# Defining models

## model with (almost) no restriction (almost saturated model)

- we just let DepEpidoc4~agE be 0, based on previous results

``` r
Model_sat<-'
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

# causalidade - idade
Dep_EpiDoc4~0*Age
MM_EpiDoc4~Age
Dep_EpiDoc2~Age
MM_EpiDoc2~Age
Dep_EpiDoc1~Age
MM_EpiDoc1~Age

# Escolaridade
Dep_EpiDoc4~Education
Dep_EpiDoc2~Education
Dep_EpiDoc1~Education
MM_EpiDoc4+MM_EpiDoc2+MM_EpiDoc1~Education
'
```

``` r
resSatModel<-lavaan.mi::sem.mi(data=data.mi,
         model=Model_sat,
         estimator = "ML", test = "yuan.bentler.mplus",
         likelihood = "wishart",
           sampling.weights='ipw',group=GVAr,
         orthogonal=TRUE
         )
summary(resSatModel,fit.meas=TRUE,standardized=TRUE,rsquare=TRUE,ci=TRUE)
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

``` r
partable<-parameterEstimates.mi(resSatModel,standardized =TRUE)
partable<-partable%>%mutate(
  label2=paste0(lhs,op,rhs),
  estimate_p=paste0(round(est,3),' (',round(pvalue,3),')'),
  CI=paste0('[',round(ci.lower,3),':',round(ci.upper,3),']')
)%>%select(-c(lhs,op,rhs))
gPartable<-reshape(partable%>%filter(!block==0), idvar = c('label2'),
                                                   timevar = "group", direction = "wide",sep='_')
gPartable%>%select(label2,estimate_p_1,estimate_p_2,CI_1,CI_2)
```

    ##                      label2   estimate_p_1   estimate_p_2              CI_1
    ## 1   Dep_EpiDoc4~Dep_EpiDoc1  0.164 (0.038)      0.342 (0)     [0.009:0.319]
    ## 2   Dep_EpiDoc4~Dep_EpiDoc2  0.234 (0.006)   0.22 (0.003)       [0.067:0.4]
    ## 3     MM_EpiDoc4~MM_EpiDoc2  0.368 (0.017)  0.282 (0.013)     [0.065:0.671]
    ## 4     MM_EpiDoc4~MM_EpiDoc1      0.564 (0)       0.35 (0)     [0.389:0.739]
    ## 5   Dep_EpiDoc2~Dep_EpiDoc1      0.407 (0)      0.374 (0)     [0.269:0.544]
    ## 6     MM_EpiDoc2~MM_EpiDoc1  0.102 (0.028)      0.366 (0)     [0.011:0.194]
    ## 7    Dep_EpiDoc4~MM_EpiDoc2   0.53 (0.069)  0.669 (0.029)    [-0.042:1.102]
    ## 8    Dep_EpiDoc4~MM_EpiDoc1  0.046 (0.824)  -0.32 (0.141)    [-0.362:0.455]
    ## 9    MM_EpiDoc4~Dep_EpiDoc2  0.036 (0.232) -0.001 (0.967)    [-0.023:0.095]
    ## 10   MM_EpiDoc4~Dep_EpiDoc1  0.005 (0.871)  0.025 (0.161)    [-0.052:0.061]
    ## 11   Dep_EpiDoc2~MM_EpiDoc1   0.281 (0.21) -0.088 (0.748)     [-0.16:0.722]
    ## 12   MM_EpiDoc2~Dep_EpiDoc1  0.021 (0.173)  0.007 (0.575)     [-0.009:0.05]
    ## 13  Dep_EpiDoc1~~MM_EpiDoc1  0.548 (0.011)  0.921 (0.003)     [0.125:0.971]
    ## 14  Dep_EpiDoc2~~MM_EpiDoc2  0.047 (0.357)    0.081 (0.6)    [-0.053:0.148]
    ## 15  Dep_EpiDoc4~~MM_EpiDoc4  0.215 (0.207) -0.026 (0.891)    [-0.119:0.549]
    ## 16          Dep_EpiDoc4~Age         0 (NA)         0 (NA)             [0:0]
    ## 17           MM_EpiDoc4~Age  0.018 (0.014)       0.03 (0)     [0.004:0.032]
    ## 18          Dep_EpiDoc2~Age  0.001 (0.974)   0.04 (0.062)    [-0.031:0.032]
    ## 19           MM_EpiDoc2~Age  0.002 (0.502) -0.005 (0.236)    [-0.003:0.006]
    ## 20          Dep_EpiDoc1~Age  0.041 (0.029)  0.055 (0.006)     [0.004:0.078]
    ## 21           MM_EpiDoc1~Age      0.031 (0)      0.049 (0)     [0.021:0.041]
    ## 22    Dep_EpiDoc4~Education  0.162 (0.329)  0.163 (0.328)    [-0.163:0.487]
    ## 23    Dep_EpiDoc2~Education       0.53 (0)  0.373 (0.077)     [0.242:0.818]
    ## 24    Dep_EpiDoc1~Education  0.495 (0.003)      0.777 (0)     [0.172:0.819]
    ## 25     MM_EpiDoc4~Education -0.006 (0.923)  0.045 (0.456)    [-0.117:0.106]
    ## 26     MM_EpiDoc2~Education -0.045 (0.128)  0.028 (0.394)    [-0.103:0.013]
    ## 27     MM_EpiDoc1~Education  0.031 (0.566)  0.025 (0.614)    [-0.075:0.138]
    ## 28 Dep_EpiDoc4~~Dep_EpiDoc4      6.879 (0)      9.339 (0)     [4.934:8.824]
    ## 29   MM_EpiDoc4~~MM_EpiDoc4      0.971 (0)      0.997 (0)     [0.763:1.179]
    ## 30 Dep_EpiDoc2~~Dep_EpiDoc2      6.187 (0)     10.657 (0)     [4.472:7.903]
    ## 31   MM_EpiDoc2~~MM_EpiDoc2      0.187 (0)        0.4 (0)     [0.086:0.287]
    ## 32 Dep_EpiDoc1~~Dep_EpiDoc1      8.427 (0)     12.629 (0)    [6.149:10.705]
    ## 33   MM_EpiDoc1~~MM_EpiDoc1        0.7 (0)      0.781 (0)      [0.531:0.87]
    ## 34                 Age~~Age   204.021 (NA)   135.726 (NA) [204.021:204.021]
    ## 35           Age~~Education     6.794 (NA)      3.68 (NA)     [6.794:6.794]
    ## 36     Education~~Education     1.159 (NA)      1.24 (NA)     [1.159:1.159]
    ## 37            Dep_EpiDoc4~1  1.073 (0.004)   0.968 (0.01)     [0.342:1.804]
    ## 38             MM_EpiDoc4~1 -0.336 (0.214) -0.611 (0.012)    [-0.866:0.194]
    ## 39            Dep_EpiDoc2~1 -0.013 (0.983) -0.517 (0.531)    [-1.269:1.243]
    ## 40             MM_EpiDoc2~1  0.041 (0.456)  0.133 (0.404)    [-0.066:0.148]
    ## 41            Dep_EpiDoc1~1 -0.575 (0.426) -0.568 (0.477)      [-1.99:0.84]
    ## 42             MM_EpiDoc1~1     -0.818 (0)     -1.396 (0)   [-1.191:-0.445]
    ## 43                    Age~1    44.188 (NA)    42.481 (NA)   [44.188:44.188]
    ## 44              Education~1     2.823 (NA)     2.409 (NA)     [2.823:2.823]
    ##                 CI_2
    ## 1      [0.195:0.489]
    ## 2      [0.073:0.366]
    ## 3      [0.059:0.504]
    ## 4      [0.179:0.522]
    ## 5       [0.21:0.538]
    ## 6      [0.176:0.556]
    ## 7       [0.067:1.27]
    ## 8     [-0.747:0.107]
    ## 9     [-0.034:0.033]
    ## 10      [-0.01:0.06]
    ## 11    [-0.623:0.448]
    ## 12    [-0.018:0.033]
    ## 13     [0.314:1.527]
    ## 14    [-0.223:0.385]
    ## 15    [-0.397:0.345]
    ## 16             [0:0]
    ## 17     [0.016:0.043]
    ## 18    [-0.002:0.083]
    ## 19    [-0.014:0.004]
    ## 20     [0.016:0.095]
    ## 21     [0.039:0.059]
    ## 22    [-0.163:0.489]
    ## 23    [-0.041:0.786]
    ## 24     [0.359:1.195]
    ## 25    [-0.073:0.163]
    ## 26    [-0.036:0.092]
    ## 27    [-0.073:0.124]
    ## 28    [7.281:11.396]
    ## 29     [0.832:1.162]
    ## 30    [8.051:13.263]
    ## 31       [0.28:0.52]
    ## 32    [9.227:16.031]
    ## 33     [0.543:1.018]
    ## 34 [135.726:135.726]
    ## 35       [3.68:3.68]
    ## 36       [1.24:1.24]
    ## 37     [0.228:1.708]
    ## 38   [-1.087:-0.135]
    ## 39    [-2.139:1.105]
    ## 40    [-0.179:0.445]
    ## 41    [-2.134:0.999]
    ## 42    [-1.84:-0.953]
    ## 43   [42.481:42.481]
    ## 44     [2.409:2.409]

# the extreme test: strong inv

``` r
resAllEd_g_strongInv<-lavaan.mi::sem.mi(data=data.mi,
         model=Model_sat,
         estimator = "ML", test = "yuan.bentler.mplus",
         likelihood = "wishart",
           sampling.weights='ipw',group=GVAr,
         orthogonal=TRUE,
         group.equal = c("intercepts", "loadings")
         )

# model comparison tests
lavaan.mi::lavTestLRT.mi(resSatModel,  resAllEd_g_strongInv)
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
    ##                      Df  Chisq Chisq diff Df diff Pr(>Chisq)     RIV     FMI
    ## resSatModel           2 2.3128                                              
    ## resAllEd_g_strongInv  8 8.5760     6.4922       6    0.37037 0.11798 0.10553

# testing partial constrains

## model with restriction on residualCov

``` r
Model_RCov<-'
# evolução
#Dep_EpiDoc4~Dep_EpiDoc1+Dep_EpiDoc2
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
MM_EpiDoc1~~md1*Dep_EpiDoc1
MM_EpiDoc2~~0* Dep_EpiDoc2
MM_EpiDoc4~~Dep_EpiDoc4

# causalidade - idade
Dep_EpiDoc4~0*Age
MM_EpiDoc4~Age
Dep_EpiDoc2~Age
MM_EpiDoc2~Age
Dep_EpiDoc1~Age
MM_EpiDoc1~Age

# Escolaridade
Dep_EpiDoc4~Education
Dep_EpiDoc2~Education
Dep_EpiDoc1~Education
MM_EpiDoc4+MM_EpiDoc2+MM_EpiDoc1~Education
'
```

``` r
resRCov<-lavaan.mi::sem.mi(data=data.mi,
         model=Model_RCov,
         estimator = "ML", test = "yuan.bentler.mplus",
         likelihood = "wishart",
           sampling.weights='ipw',group=GVAr,
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
summary(resRCov,fit.meas=TRUE,standardized=TRUE,rsquare=TRUE,ci=TRUE)
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
    ##   Number of model parameters                        74
    ##   Number of equality constraints                     1
    ## 
    ##   Number of observations per group:                   
    ##     Male                                           355
    ##     Female                                         317
    ##   Sampling weights variable                        ipw
    ## 
    ## Model Test User Model:
    ## 
    ##                                                     Standard      Scaled
    ##   Test statistic                                       5.514       3.966
    ##   Degrees of freedom                                       5           5
    ##   P-value                                              0.356       0.554
    ##   Average scaling correction factor                                1.390
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
    ##   Comparative Fit Index (CFI)                    0.999       1.000
    ##   Tucker-Lewis Index (TLI)                       0.994       1.014
    ##                                                                   
    ##   Robust Comparative Fit Index (CFI)                         1.000
    ##   Robust Tucker-Lewis Index (TLI)                            1.016
    ## 
    ## Loglikelihood and Information Criteria:
    ## 
    ##   Loglikelihood user model (H0)              -7328.602   -7328.602
    ##   Scaling correction factor                                  1.720
    ##       for the MLR correction                                      
    ##   Loglikelihood unrestricted model (H1)      -7325.594   -7325.594
    ##   Scaling correction factor                                  1.717
    ##       for the MLR correction                                      
    ##                                                                   
    ##   Akaike (AIC)                               14803.203   14803.203
    ##   Bayesian (BIC)                             15132.452   15132.452
    ##   Sample-size adjusted Bayesian (SABIC)      14900.671   14900.671
    ## 
    ## Root Mean Square Error of Approximation:
    ## 
    ##   RMSEA                                          0.018       0.000
    ##   90 Percent confidence interval - lower         0.000       0.000
    ##   90 Percent confidence interval - upper         0.080       0.061
    ##   P-value H_0: RMSEA <= 0.050                    0.737       0.894
    ##   P-value H_0: RMSEA >= 0.080                    0.048       0.008
    ##                                                                   
    ##   Robust RMSEA                                               0.000
    ##   90 Percent confidence interval - lower                     0.000
    ##   90 Percent confidence interval - upper                     0.080
    ##   P-value H_0: Robust RMSEA <= 0.050                         0.803
    ##   P-value H_0: Robust RMSEA >= 0.080                         0.049
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
    ##     Dep_EpiDoc1       0.164    0.079    2.086  527.390    0.037    0.010
    ##     Dep_EpiDoc2       0.234    0.084    2.770  327.278    0.006    0.068
    ##   MM_EpiDoc4 ~                                                          
    ##     MM_EpiDoc2        0.368    0.154    2.388      Inf    0.017    0.066
    ##     MM_EpiDoc1        0.564    0.089    6.335      Inf    0.000    0.389
    ##   Dep_EpiDoc2 ~                                                         
    ##     Dep_EpiDoc1       0.407    0.070    5.831  356.369    0.000    0.270
    ##   MM_EpiDoc2 ~                                                          
    ##     MM_EpiDoc1        0.102    0.046    2.202      Inf    0.028    0.011
    ##   Dep_EpiDoc4 ~                                                         
    ##     MM_EpiDoc2        0.530    0.291    1.823      Inf    0.068   -0.040
    ##     MM_EpiDoc1        0.046    0.208    0.224      Inf    0.823   -0.361
    ##   MM_EpiDoc4 ~                                                          
    ##     Dep_EpiDoc2       0.036    0.030    1.203  181.697    0.230   -0.023
    ##     Dep_EpiDoc1       0.005    0.029    0.163      Inf    0.871   -0.052
    ##   Dep_EpiDoc2 ~                                                         
    ##     MM_EpiDoc1        0.281    0.223    1.260  272.102    0.209   -0.158
    ##   MM_EpiDoc2 ~                                                          
    ##     Dep_EpiDoc1       0.021    0.015    1.368      Inf    0.171   -0.009
    ##   Dep_EpiDoc4 ~                                                         
    ##     Age               0.000                                        0.000
    ##   MM_EpiDoc4 ~                                                          
    ##     Age               0.018    0.007    2.464      Inf    0.014    0.004
    ##   Dep_EpiDoc2 ~                                                         
    ##     Age               0.001    0.016    0.033   77.311    0.974   -0.031
    ##   MM_EpiDoc2 ~                                                          
    ##     Age               0.002    0.002    0.673      Inf    0.501   -0.003
    ##   Dep_EpiDoc1 ~                                                         
    ##     Age               0.041    0.019    2.195      Inf    0.028    0.004
    ##   MM_EpiDoc1 ~                                                          
    ##     Age               0.031    0.005    6.100      Inf    0.000    0.021
    ##   Dep_EpiDoc4 ~                                                         
    ##     Education         0.162    0.165    0.981  599.382    0.327   -0.162
    ##   Dep_EpiDoc2 ~                                                         
    ##     Education         0.530    0.146    3.631  451.180    0.000    0.243
    ##   Dep_EpiDoc1 ~                                                         
    ##     Education         0.495    0.164    3.014      Inf    0.003    0.173
    ##   MM_EpiDoc4 ~                                                          
    ##     Education        -0.006    0.057   -0.097      Inf    0.923   -0.117
    ##   MM_EpiDoc2 ~                                                          
    ##     Education        -0.045    0.029   -1.529      Inf    0.126   -0.103
    ##   MM_EpiDoc1 ~                                                          
    ##     Education         0.031    0.054    0.576      Inf    0.564   -0.075
    ##  ci.upper   Std.lv  Std.all
    ##                            
    ##     0.318    0.164    0.176
    ##     0.400    0.234    0.241
    ##                            
    ##     0.670    0.368    0.131
    ##     0.738    0.564    0.427
    ##                            
    ##     0.544    0.407    0.422
    ##                            
    ##     0.193    0.102    0.217
    ##                            
    ##     1.100    0.530    0.083
    ##     0.454    0.046    0.016
    ##                            
    ##     0.094    0.036    0.083
    ##     0.061    0.005    0.011
    ##                            
    ##     0.721    0.281    0.091
    ##                            
    ##     0.050    0.021    0.141
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.032    0.018    0.202
    ##                            
    ##     0.032    0.001    0.002
    ##                            
    ##     0.006    0.002    0.049
    ##                            
    ##     0.077    0.041    0.189
    ##                            
    ##     0.041    0.031    0.460
    ##                            
    ##     0.486    0.162    0.060
    ##                            
    ##     0.817    0.530    0.192
    ##                            
    ##     0.818    0.495    0.173
    ##                            
    ##     0.106   -0.006   -0.005
    ##                            
    ##     0.013   -0.045   -0.107
    ##                            
    ##     0.137    0.031    0.035
    ## 
    ## Covariances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##  .Dep_EpiDoc1 ~~                                                        
    ##    .MM_EpDc1 (md1)    0.677    0.179    3.779      Inf    0.000    0.326
    ##  .Dep_EpiDoc2 ~~                                                        
    ##    .MM_EpDc2          0.000                                        0.000
    ##  .Dep_EpiDoc4 ~~                                                        
    ##    .MM_EpDc4          0.215    0.169    1.271  195.736    0.205   -0.118
    ##  ci.upper   Std.lv  Std.all
    ##                            
    ##     1.029    0.677    0.271
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.547    0.215    0.083
    ## 
    ## Intercepts:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##    .Dep_EpiDoc4       1.073    0.371    2.891  828.813    0.004    0.345
    ##    .MM_EpiDoc4       -0.336    0.269   -1.246      Inf    0.213   -0.864
    ##    .Dep_EpiDoc2      -0.013    0.631   -0.021  100.906    0.983   -1.265
    ##    .MM_EpiDoc2        0.041    0.054    0.748      Inf    0.455   -0.066
    ##    .Dep_EpiDoc1      -0.575    0.720   -0.799      Inf    0.424   -1.985
    ##    .MM_EpiDoc1       -0.818    0.190   -4.310      Inf    0.000   -1.190
    ##  ci.upper   Std.lv  Std.all
    ##     1.802    1.073    0.372
    ##     0.192   -0.336   -0.263
    ##     1.238   -0.013   -0.004
    ##     0.147    0.041    0.090
    ##     0.835   -0.575   -0.186
    ##    -0.446   -0.818   -0.848
    ## 
    ## Variances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##    .Dep_EpiDoc4       6.879    0.988    6.960      Inf    0.000    4.941
    ##    .MM_EpiDoc4        0.971    0.106    9.178      Inf    0.000    0.764
    ##    .Dep_EpiDoc2       6.187    0.862    7.177  101.737    0.000    4.477
    ##    .MM_EpiDoc2        0.187    0.051    3.662      Inf    0.000    0.087
    ##    .Dep_EpiDoc1       8.654    1.159    7.465      Inf    0.000    6.382
    ##    .MM_EpiDoc1        0.719    0.086    8.386      Inf    0.000    0.551
    ##  ci.upper   Std.lv  Std.all
    ##     8.818    6.879    0.827
    ##     1.178    0.971    0.598
    ##     7.898    6.187    0.698
    ##     0.287    0.187    0.907
    ##    10.926    8.654    0.906
    ##     0.887    0.719    0.773
    ## 
    ## R-Square:
    ##                    Estimate
    ##     Dep_EpiDoc4       0.173
    ##     MM_EpiDoc4        0.402
    ##     Dep_EpiDoc2       0.302
    ##     MM_EpiDoc2        0.093
    ##     Dep_EpiDoc1       0.094
    ##     MM_EpiDoc1        0.227
    ## 
    ## 
    ## Group 2 [Female]:
    ## 
    ## Regressions:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##   Dep_EpiDoc4 ~                                                         
    ##     Dep_EpiDoc1       0.342    0.075    4.570      Inf    0.000    0.195
    ##     Dep_EpiDoc2       0.220    0.074    2.951  943.731    0.003    0.074
    ##   MM_EpiDoc4 ~                                                          
    ##     MM_EpiDoc2        0.282    0.113    2.489      Inf    0.013    0.060
    ##     MM_EpiDoc1        0.350    0.087    4.013      Inf    0.000    0.179
    ##   Dep_EpiDoc2 ~                                                         
    ##     Dep_EpiDoc1       0.374    0.083    4.492  814.637    0.000    0.210
    ##   MM_EpiDoc2 ~                                                          
    ##     MM_EpiDoc1        0.366    0.097    3.789      Inf    0.000    0.177
    ##   Dep_EpiDoc4 ~                                                         
    ##     MM_EpiDoc2        0.669    0.306    2.187      Inf    0.029    0.069
    ##     MM_EpiDoc1       -0.320    0.217   -1.476      Inf    0.140   -0.746
    ##   MM_EpiDoc4 ~                                                          
    ##     Dep_EpiDoc2      -0.001    0.017   -0.042  495.075    0.967   -0.034
    ##     Dep_EpiDoc1       0.025    0.018    1.405      Inf    0.160   -0.010
    ##   Dep_EpiDoc2 ~                                                         
    ##     MM_EpiDoc1       -0.088    0.272   -0.323  585.338    0.747   -0.621
    ##   MM_EpiDoc2 ~                                                          
    ##     Dep_EpiDoc1       0.007    0.013    0.563      Inf    0.573   -0.018
    ##   Dep_EpiDoc4 ~                                                         
    ##     Age               0.000                                        0.000
    ##   MM_EpiDoc4 ~                                                          
    ##     Age               0.030    0.007    4.336      Inf    0.000    0.016
    ##   Dep_EpiDoc2 ~                                                         
    ##     Age               0.040    0.022    1.873  689.305    0.062   -0.002
    ##   MM_EpiDoc2 ~                                                          
    ##     Age              -0.005    0.005   -1.189      Inf    0.234   -0.014
    ##   Dep_EpiDoc1 ~                                                         
    ##     Age               0.055    0.020    2.771      Inf    0.006    0.016
    ##   MM_EpiDoc1 ~                                                          
    ##     Age               0.049    0.005    9.604      Inf    0.000    0.039
    ##   Dep_EpiDoc4 ~                                                         
    ##     Education         0.163    0.166    0.981      Inf    0.326   -0.162
    ##   Dep_EpiDoc2 ~                                                         
    ##     Education         0.373    0.209    1.780  309.152    0.076   -0.039
    ##   Dep_EpiDoc1 ~                                                         
    ##     Education         0.777    0.213    3.656      Inf    0.000    0.361
    ##   MM_EpiDoc4 ~                                                          
    ##     Education         0.045    0.060    0.748      Inf    0.454   -0.073
    ##   MM_EpiDoc2 ~                                                          
    ##     Education         0.028    0.033    0.856      Inf    0.392   -0.036
    ##   MM_EpiDoc1 ~                                                          
    ##     Education         0.025    0.050    0.506      Inf    0.613   -0.073
    ##  ci.upper   Std.lv  Std.all
    ##                            
    ##     0.488    0.342    0.353
    ##     0.366    0.220    0.226
    ##                            
    ##     0.503    0.282    0.162
    ##     0.521    0.350    0.288
    ##                            
    ##     0.537    0.374    0.375
    ##                            
    ##     0.555    0.366    0.522
    ##                            
    ##     1.268    0.669    0.136
    ##     0.105   -0.320   -0.093
    ##                            
    ##     0.033   -0.001   -0.002
    ##     0.060    0.025    0.073
    ##                            
    ##     0.446   -0.088   -0.025
    ##                            
    ##     0.033    0.007    0.037
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.043    0.030    0.271
    ##                            
    ##     0.083    0.040    0.127
    ##                            
    ##     0.004   -0.005   -0.086
    ##                            
    ##     0.095    0.055    0.174
    ##                            
    ##     0.059    0.049    0.543
    ##                            
    ##     0.488    0.163    0.051
    ##                            
    ##     0.785    0.373    0.113
    ##                            
    ##     1.194    0.777    0.234
    ##                            
    ##     0.163    0.045    0.039
    ##                            
    ##     0.092    0.028    0.043
    ##                            
    ##     0.123    0.025    0.027
    ## 
    ## Covariances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##  .Dep_EpiDoc1 ~~                                                        
    ##    .MM_EpDc1 (md1)    0.677    0.179    3.779      Inf    0.000    0.326
    ##  .Dep_EpiDoc2 ~~                                                        
    ##    .MM_EpDc2          0.000                                        0.000
    ##  .Dep_EpiDoc4 ~~                                                        
    ##    .MM_EpDc4         -0.026    0.189   -0.138      Inf    0.890   -0.396
    ##  ci.upper   Std.lv  Std.all
    ##                            
    ##     1.029    0.677    0.223
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.344   -0.026   -0.009
    ## 
    ## Intercepts:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##    .Dep_EpiDoc4       0.968    0.376    2.574      Inf    0.010    0.231
    ##    .MM_EpiDoc4       -0.611    0.242   -2.522      Inf    0.012   -1.086
    ##    .Dep_EpiDoc2      -0.517    0.822   -0.629  381.594    0.530   -2.133
    ##    .MM_EpiDoc2        0.133    0.159    0.837      Inf    0.403   -0.178
    ##    .Dep_EpiDoc1      -0.568    0.797   -0.713      Inf    0.476   -2.129
    ##    .MM_EpiDoc1       -1.396    0.226   -6.190      Inf    0.000   -1.839
    ##  ci.upper   Std.lv  Std.all
    ##     1.706    0.968    0.270
    ##    -0.136   -0.611   -0.481
    ##     1.099   -0.517   -0.140
    ##     0.444    0.133    0.182
    ##     0.994   -0.568   -0.153
    ##    -0.954   -1.396   -1.340
    ## 
    ## Variances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##    .Dep_EpiDoc4       9.339    1.045    8.933      Inf    0.000    7.288
    ##    .MM_EpiDoc4        0.997    0.084   11.864      Inf    0.000    0.832
    ##    .Dep_EpiDoc2      10.657    1.320    8.071  331.705    0.000    8.060
    ##    .MM_EpiDoc2        0.400    0.061    6.543      Inf    0.000    0.280
    ##    .Dep_EpiDoc1      12.214    1.524    8.014      Inf    0.000    9.227
    ##    .MM_EpiDoc1        0.755    0.107    7.050      Inf    0.000    0.545
    ##  ci.upper   Std.lv  Std.all
    ##    11.390    9.339    0.727
    ##     1.162    0.997    0.619
    ##    13.255   10.657    0.785
    ##     0.520    0.400    0.749
    ##    15.201   12.214    0.892
    ##     0.965    0.755    0.696
    ## 
    ## R-Square:
    ##                    Estimate
    ##     Dep_EpiDoc4       0.273
    ##     MM_EpiDoc4        0.381
    ##     Dep_EpiDoc2       0.215
    ##     MM_EpiDoc2        0.251
    ##     Dep_EpiDoc1       0.108
    ##     MM_EpiDoc1        0.304

### Test significance free model x resRCov

``` r
lavaan.mi::lavTestLRT.mi(resRCov,resSatModel )
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
    ##             Df  Chisq Chisq diff Df diff Pr(>Chisq)     RIV     FMI
    ## resSatModel  2 2.3128                                              
    ## resRCov      5 5.5145     1.9997       3    0.57247 0.21946 0.17996

## model with restriction on AR components

``` r
Model_AR<-'
# evolução
#Dep_EpiDoc4~dd03*Dep_EpiDoc1+dd13*Dep_EpiDoc2
Dep_EpiDoc4~Dep_EpiDoc1+dd13*Dep_EpiDoc2
MM_EpiDoc4~mm14*MM_EpiDoc2+mm04*MM_EpiDoc1
Dep_EpiDoc2~dd01*Dep_EpiDoc1
MM_EpiDoc2~mm01*MM_EpiDoc1 # freed up considering mi

# causalidade - c/ temporalidade/lag
Dep_EpiDoc4~MM_EpiDoc2+MM_EpiDoc1
MM_EpiDoc4~Dep_EpiDoc2+Dep_EpiDoc1
Dep_EpiDoc2~MM_EpiDoc1
MM_EpiDoc2~Dep_EpiDoc1

# corrrelação mesmo t
MM_EpiDoc1~~0*Dep_EpiDoc1
MM_EpiDoc2~~0* Dep_EpiDoc2
MM_EpiDoc4~~Dep_EpiDoc4

# causalidade - idade
Dep_EpiDoc4~0*Age
MM_EpiDoc4~Age
Dep_EpiDoc2~Age
MM_EpiDoc2~Age
Dep_EpiDoc1~Age
MM_EpiDoc1~Age

# Escolaridade
Dep_EpiDoc4~Education
Dep_EpiDoc2~Education
Dep_EpiDoc1~Education
MM_EpiDoc4+MM_EpiDoc2+MM_EpiDoc1~Education
'
```

``` r
resAR<-lavaan.mi::sem.mi(data=data.mi,
         model=Model_AR,
         estimator = "ML", test = "yuan.bentler.mplus",
         likelihood = "wishart",
           sampling.weights='ipw',group=GVAr,
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
summary(resAR,fit.meas=TRUE,standardized=TRUE,rsquare=TRUE,ci=TRUE)
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
    ##   Test statistic                                      70.295      44.097
    ##   Degrees of freedom                                      11          11
    ##   P-value                                              0.000       0.000
    ##   Average scaling correction factor                                1.594
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
    ##   Comparative Fit Index (CFI)                    0.940       0.958
    ##   Tucker-Lewis Index (TLI)                       0.707       0.795
    ##                                                                   
    ##   Robust Comparative Fit Index (CFI)                         0.946
    ##   Robust Tucker-Lewis Index (TLI)                            0.736
    ## 
    ## Loglikelihood and Information Criteria:
    ## 
    ##   Loglikelihood user model (H0)              -7361.000   -7361.000
    ##   Scaling correction factor                                  1.620
    ##       for the MLR correction                                      
    ##   Loglikelihood unrestricted model (H1)      -7325.594   -7325.594
    ##   Scaling correction factor                                  1.717
    ##       for the MLR correction                                      
    ##                                                                   
    ##   Akaike (AIC)                               14855.999   14855.999
    ##   Bayesian (BIC)                             15158.187   15158.187
    ##   Sample-size adjusted Bayesian (SABIC)      14945.456   14945.456
    ## 
    ## Root Mean Square Error of Approximation:
    ## 
    ##   RMSEA                                          0.127       0.096
    ##   90 Percent confidence interval - lower         0.099       0.073
    ##   90 Percent confidence interval - upper         0.156       0.119
    ##   P-value H_0: RMSEA <= 0.050                    0.000       0.001
    ##   P-value H_0: RMSEA >= 0.080                    0.997       0.877
    ##                                                                   
    ##   Robust RMSEA                                               0.120
    ##   90 Percent confidence interval - lower                     0.084
    ##   90 Percent confidence interval - upper                     0.158
    ##   P-value H_0: Robust RMSEA <= 0.050                         0.001
    ##   P-value H_0: Robust RMSEA >= 0.080                         0.965
    ## 
    ## Standardized Root Mean Square Residual:
    ## 
    ##   SRMR                                           0.055       0.055
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
    ##     Dp_EpD1           0.168    0.078    2.160  764.760    0.031    0.015
    ##     Dp_EpD2 (dd13)    0.228    0.056    4.047  557.171    0.000    0.117
    ##   MM_EpiDoc4 ~                                                          
    ##     MM_EpD2 (mm14)    0.279    0.094    2.955      Inf    0.003    0.094
    ##     MM_EpD1 (mm04)    0.474    0.063    7.482      Inf    0.000    0.349
    ##   Dep_EpiDoc2 ~                                                         
    ##     Dp_EpD1 (dd01)    0.394    0.054    7.321  335.639    0.000    0.288
    ##   MM_EpiDoc2 ~                                                          
    ##     MM_EpD1 (mm01)    0.173    0.054    3.231      Inf    0.001    0.068
    ##   Dep_EpiDoc4 ~                                                         
    ##     MM_EpD2           0.513    0.292    1.755      Inf    0.079   -0.060
    ##     MM_EpD1           0.032    0.209    0.153      Inf    0.879   -0.378
    ##   MM_EpiDoc4 ~                                                          
    ##     Dp_EpD2           0.039    0.030    1.324  159.587    0.188   -0.019
    ##     Dp_EpD1           0.011    0.030    0.385      Inf    0.701   -0.047
    ##   Dep_EpiDoc2 ~                                                         
    ##     MM_EpD1           0.291    0.224    1.301  254.835    0.194   -0.150
    ##   MM_EpiDoc2 ~                                                          
    ##     Dp_EpD1           0.016    0.015    1.084      Inf    0.278   -0.013
    ##   Dep_EpiDoc4 ~                                                         
    ##     Age               0.000                                        0.000
    ##   MM_EpiDoc4 ~                                                          
    ##     Age               0.021    0.007    2.897      Inf    0.004    0.007
    ##   Dep_EpiDoc2 ~                                                         
    ##     Age               0.001    0.016    0.047   78.786    0.963   -0.030
    ##   MM_EpiDoc2 ~                                                          
    ##     Age              -0.000    0.002   -0.198      Inf    0.843   -0.005
    ##   Dep_EpiDoc1 ~                                                         
    ##     Age               0.041    0.019    2.194      Inf    0.028    0.004
    ##   MM_EpiDoc1 ~                                                          
    ##     Age               0.031    0.005    6.098      Inf    0.000    0.021
    ##   Dep_EpiDoc4 ~                                                         
    ##     Educatn           0.167    0.160    1.044  610.416    0.297   -0.147
    ##   Dep_EpiDoc2 ~                                                         
    ##     Educatn           0.536    0.144    3.712  480.882    0.000    0.252
    ##   Dep_EpiDoc1 ~                                                         
    ##     Educatn           0.495    0.164    3.013      Inf    0.003    0.173
    ##   MM_EpiDoc4 ~                                                          
    ##     Educatn          -0.012    0.057   -0.202      Inf    0.840   -0.123
    ##   MM_EpiDoc2 ~                                                          
    ##     Educatn          -0.045    0.030   -1.501      Inf    0.133   -0.103
    ##   MM_EpiDoc1 ~                                                          
    ##     Educatn           0.031    0.054    0.576      Inf    0.565   -0.075
    ##  ci.upper   Std.lv  Std.all
    ##                            
    ##     0.321    0.168    0.179
    ##     0.339    0.228    0.233
    ##                            
    ##     0.464    0.279    0.105
    ##     0.598    0.474    0.365
    ##                            
    ##     0.500    0.394    0.410
    ##                            
    ##     0.279    0.173    0.354
    ##                            
    ##     1.087    0.513    0.084
    ##     0.441    0.032    0.011
    ##                            
    ##     0.098    0.039    0.093
    ##     0.070    0.011    0.028
    ##                            
    ##     0.732    0.291    0.095
    ##                            
    ##     0.045    0.016    0.105
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.035    0.021    0.242
    ##                            
    ##     0.032    0.001    0.004
    ##                            
    ##     0.004   -0.000   -0.014
    ##                            
    ##     0.078    0.041    0.191
    ##                            
    ##     0.041    0.031    0.465
    ##                            
    ##     0.482    0.167    0.063
    ##                            
    ##     0.820    0.536    0.197
    ##                            
    ##     0.818    0.495    0.175
    ##                            
    ##     0.100   -0.012   -0.010
    ##                            
    ##     0.014   -0.045   -0.103
    ##                            
    ##     0.137    0.031    0.035
    ## 
    ## Covariances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##  .Dep_EpiDoc1 ~~                                                        
    ##    .MM_EpiDoc1        0.000                                        0.000
    ##  .Dep_EpiDoc2 ~~                                                        
    ##    .MM_EpiDoc2        0.000                                        0.000
    ##  .Dep_EpiDoc4 ~~                                                        
    ##    .MM_EpiDoc4        0.208    0.171    1.218  188.859    0.225   -0.129
    ##  ci.upper   Std.lv  Std.all
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.545    0.208    0.080
    ## 
    ## Intercepts:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##    .Dep_EpiDoc4       1.075    0.371    2.900  867.249    0.004    0.347
    ##    .MM_EpiDoc4       -0.409    0.271   -1.510      Inf    0.131   -0.940
    ##    .Dep_EpiDoc2      -0.012    0.633   -0.020  101.285    0.984   -1.268
    ##    .MM_EpiDoc2        0.096    0.058    1.663      Inf    0.096   -0.017
    ##    .Dep_EpiDoc1      -0.575    0.720   -0.799      Inf    0.424   -1.986
    ##    .MM_EpiDoc1       -0.818    0.190   -4.308      Inf    0.000   -1.190
    ##  ci.upper   Std.lv  Std.all
    ##     1.803    1.075    0.375
    ##     0.122   -0.409   -0.330
    ##     1.243   -0.012   -0.004
    ##     0.209    0.096    0.206
    ##     0.836   -0.575   -0.188
    ##    -0.446   -0.818   -0.856
    ## 
    ## Variances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##    .Dep_EpiDoc4       6.884    0.991    6.949      Inf    0.000    4.942
    ##    .MM_EpiDoc4        0.978    0.107    9.180      Inf    0.000    0.770
    ##    .Dep_EpiDoc2       6.191    0.864    7.164  101.905    0.000    4.477
    ##    .MM_EpiDoc2        0.190    0.052    3.683      Inf    0.000    0.089
    ##    .Dep_EpiDoc1       8.427    1.159    7.271      Inf    0.000    6.155
    ##    .MM_EpiDoc1        0.700    0.086    8.128      Inf    0.000    0.531
    ##  ci.upper   Std.lv  Std.all
    ##     8.827    6.884    0.838
    ##     1.187    0.978    0.638
    ##     7.905    6.191    0.721
    ##     0.291    0.190    0.870
    ##    10.698    8.427    0.903
    ##     0.869    0.700    0.768
    ## 
    ## R-Square:
    ##                    Estimate
    ##     Dep_EpiDoc4       0.162
    ##     MM_EpiDoc4        0.362
    ##     Dep_EpiDoc2       0.279
    ##     MM_EpiDoc2        0.130
    ##     Dep_EpiDoc1       0.097
    ##     MM_EpiDoc1        0.232
    ## 
    ## 
    ## Group 2 [Female]:
    ## 
    ## Regressions:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##   Dep_EpiDoc4 ~                                                         
    ##     Dp_EpD1           0.339    0.071    4.773      Inf    0.000    0.200
    ##     Dp_EpD2 (dd13)    0.228    0.056    4.047  557.171    0.000    0.117
    ##   MM_EpiDoc4 ~                                                          
    ##     MM_EpD2 (mm14)    0.279    0.094    2.955      Inf    0.003    0.094
    ##     MM_EpD1 (mm04)    0.474    0.063    7.482      Inf    0.000    0.349
    ##   Dep_EpiDoc2 ~                                                         
    ##     Dp_EpD1 (dd01)    0.394    0.054    7.321  335.639    0.000    0.288
    ##   MM_EpiDoc2 ~                                                          
    ##     MM_EpD1 (mm01)    0.173    0.054    3.231      Inf    0.001    0.068
    ##   Dep_EpiDoc4 ~                                                         
    ##     MM_EpD2           0.665    0.305    2.179      Inf    0.029    0.066
    ##     MM_EpD1          -0.324    0.216   -1.500      Inf    0.134   -0.748
    ##   MM_EpiDoc4 ~                                                          
    ##     Dp_EpD2          -0.000    0.017   -0.000  533.366    1.000   -0.033
    ##     Dp_EpD1           0.016    0.018    0.902      Inf    0.367   -0.019
    ##   Dep_EpiDoc2 ~                                                         
    ##     MM_EpD1          -0.111    0.275   -0.405  611.454    0.685   -0.650
    ##   MM_EpiDoc2 ~                                                          
    ##     Dp_EpD1           0.021    0.015    1.461      Inf    0.144   -0.007
    ##   Dep_EpiDoc4 ~                                                         
    ##     Age               0.000                                        0.000
    ##   MM_EpiDoc4 ~                                                          
    ##     Age               0.024    0.006    3.967      Inf    0.000    0.012
    ##   Dep_EpiDoc2 ~                                                         
    ##     Age               0.040    0.022    1.870  692.829    0.062   -0.002
    ##   MM_EpiDoc2 ~                                                          
    ##     Age               0.003    0.004    0.761      Inf    0.447   -0.005
    ##   Dep_EpiDoc1 ~                                                         
    ##     Age               0.055    0.020    2.769      Inf    0.006    0.016
    ##   MM_EpiDoc1 ~                                                          
    ##     Age               0.049    0.005    9.600      Inf    0.000    0.039
    ##   Dep_EpiDoc4 ~                                                         
    ##     Educatn           0.158    0.168    0.943      Inf    0.346   -0.171
    ##   Dep_EpiDoc2 ~                                                         
    ##     Educatn           0.358    0.204    1.758  298.000    0.080   -0.043
    ##   Dep_EpiDoc1 ~                                                         
    ##     Educatn           0.777    0.213    3.655      Inf    0.000    0.360
    ##   MM_EpiDoc4 ~                                                          
    ##     Educatn           0.049    0.061    0.803      Inf    0.422   -0.070
    ##   MM_EpiDoc2 ~                                                          
    ##     Educatn           0.022    0.033    0.665      Inf    0.506   -0.043
    ##   MM_EpiDoc1 ~                                                          
    ##     Educatn           0.025    0.050    0.505      Inf    0.613   -0.073
    ##  ci.upper   Std.lv  Std.all
    ##                            
    ##     0.478    0.339    0.350
    ##     0.339    0.228    0.234
    ##                            
    ##     0.464    0.279    0.151
    ##     0.598    0.474    0.387
    ##                            
    ##     0.500    0.394    0.397
    ##                            
    ##     0.279    0.173    0.263
    ##                            
    ##     1.263    0.665    0.127
    ##     0.100   -0.324   -0.094
    ##                            
    ##     0.033   -0.000   -0.000
    ##     0.051    0.016    0.046
    ##                            
    ##     0.428   -0.111   -0.031
    ##                            
    ##     0.050    0.021    0.116
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.036    0.024    0.217
    ##                            
    ##     0.083    0.040    0.126
    ##                            
    ##     0.011    0.003    0.053
    ##                            
    ##     0.095    0.055    0.172
    ##                            
    ##     0.059    0.049    0.537
    ##                            
    ##     0.488    0.158    0.048
    ##                            
    ##     0.759    0.358    0.107
    ##                            
    ##     1.194    0.777    0.230
    ##                            
    ##     0.167    0.049    0.042
    ##                            
    ##     0.086    0.022    0.035
    ##                            
    ##     0.123    0.025    0.027
    ## 
    ## Covariances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##  .Dep_EpiDoc1 ~~                                                        
    ##    .MM_EpiDoc1        0.000                                        0.000
    ##  .Dep_EpiDoc2 ~~                                                        
    ##    .MM_EpiDoc2        0.000                                        0.000
    ##  .Dep_EpiDoc4 ~~                                                        
    ##    .MM_EpiDoc4       -0.027    0.191   -0.141      Inf    0.888   -0.400
    ##  ci.upper   Std.lv  Std.all
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.347   -0.027   -0.009
    ## 
    ## Intercepts:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##    .Dep_EpiDoc4       0.965    0.373    2.583      Inf    0.010    0.232
    ##    .MM_EpiDoc4       -0.445    0.221   -2.016      Inf    0.044   -0.877
    ##    .Dep_EpiDoc2      -0.539    0.833   -0.646  382.529    0.518   -2.177
    ##    .MM_EpiDoc2       -0.128    0.171   -0.750      Inf    0.453   -0.464
    ##    .Dep_EpiDoc1      -0.568    0.797   -0.712      Inf    0.476   -2.130
    ##    .MM_EpiDoc1       -1.396    0.226   -6.187      Inf    0.000   -1.839
    ##  ci.upper   Std.lv  Std.all
    ##     1.697    0.965    0.265
    ##    -0.012   -0.445   -0.345
    ##     1.100   -0.539   -0.144
    ##     0.207   -0.128   -0.184
    ##     0.994   -0.568   -0.151
    ##    -0.954   -1.396   -1.325
    ## 
    ## Variances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##    .Dep_EpiDoc4       9.348    1.048    8.920      Inf    0.000    7.292
    ##    .MM_EpiDoc4        1.008    0.086   11.649      Inf    0.000    0.838
    ##    .Dep_EpiDoc2      10.668    1.322    8.072  337.091    0.000    8.068
    ##    .MM_EpiDoc2        0.427    0.074    5.791      Inf    0.000    0.282
    ##    .Dep_EpiDoc1      12.629    1.731    7.297      Inf    0.000    9.237
    ##    .MM_EpiDoc1        0.781    0.121    6.464      Inf    0.000    0.544
    ##  ci.upper   Std.lv  Std.all
    ##    11.404    9.348    0.706
    ##     1.177    1.008    0.607
    ##    13.268   10.668    0.768
    ##     0.571    0.427    0.881
    ##    16.022   12.629    0.895
    ##     1.017    0.781    0.703
    ## 
    ## R-Square:
    ##                    Estimate
    ##     Dep_EpiDoc4       0.294
    ##     MM_EpiDoc4        0.393
    ##     Dep_EpiDoc2       0.232
    ##     MM_EpiDoc2        0.119
    ##     Dep_EpiDoc1       0.105
    ##     MM_EpiDoc1        0.297

### Test significance free model x AR stricted

``` r
lavaan.mi::lavTestLRT.mi(resAR,resSatModel )
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
    ##             Df   Chisq Chisq diff Df diff Pr(>Chisq)     RIV     FMI
    ## resSatModel  2  2.3128                                              
    ## resAR       11 70.2953     41.191       9 4.6174e-06 0.12821 0.11364

``` r
partable<-parameterEstimates.mi(resAR,standardized =TRUE)
partable<-partable%>%mutate(
  label2=paste0(lhs,op,rhs),
  estimate_p=paste0(round(est,3),' (',round(pvalue,3),')'),
  CI=paste0('[',round(ci.lower,3),':',round(ci.upper,3),']')
)%>%select(-c(lhs,op,rhs))
gPartable<-reshape(partable%>%filter(!block==0), idvar = c('label2'),
                                                   timevar = "group", direction = "wide",sep='_')
gPartable%>%select(label2,estimate_p_1,estimate_p_2,CI_1,CI_2)
```

    ##                      label2   estimate_p_1   estimate_p_2              CI_1
    ## 1   Dep_EpiDoc4~Dep_EpiDoc1  0.168 (0.031)      0.339 (0)     [0.015:0.321]
    ## 2   Dep_EpiDoc4~Dep_EpiDoc2      0.228 (0)      0.228 (0)     [0.117:0.339]
    ## 3     MM_EpiDoc4~MM_EpiDoc2  0.279 (0.003)  0.279 (0.003)     [0.094:0.464]
    ## 4     MM_EpiDoc4~MM_EpiDoc1      0.474 (0)      0.474 (0)     [0.349:0.598]
    ## 5   Dep_EpiDoc2~Dep_EpiDoc1      0.394 (0)      0.394 (0)       [0.288:0.5]
    ## 6     MM_EpiDoc2~MM_EpiDoc1  0.173 (0.001)  0.173 (0.001)     [0.068:0.279]
    ## 7    Dep_EpiDoc4~MM_EpiDoc2  0.513 (0.079)  0.665 (0.029)     [-0.06:1.087]
    ## 8    Dep_EpiDoc4~MM_EpiDoc1  0.032 (0.879) -0.324 (0.134)    [-0.378:0.441]
    ## 9    MM_EpiDoc4~Dep_EpiDoc2  0.039 (0.188)          0 (1)    [-0.019:0.098]
    ## 10   MM_EpiDoc4~Dep_EpiDoc1  0.011 (0.701)  0.016 (0.367)     [-0.047:0.07]
    ## 11   Dep_EpiDoc2~MM_EpiDoc1  0.291 (0.194) -0.111 (0.685)     [-0.15:0.732]
    ## 12   MM_EpiDoc2~Dep_EpiDoc1  0.016 (0.278)  0.021 (0.144)    [-0.013:0.045]
    ## 13  Dep_EpiDoc1~~MM_EpiDoc1         0 (NA)         0 (NA)             [0:0]
    ## 14  Dep_EpiDoc2~~MM_EpiDoc2         0 (NA)         0 (NA)             [0:0]
    ## 15  Dep_EpiDoc4~~MM_EpiDoc4  0.208 (0.225) -0.027 (0.888)    [-0.129:0.545]
    ## 16          Dep_EpiDoc4~Age         0 (NA)         0 (NA)             [0:0]
    ## 17           MM_EpiDoc4~Age  0.021 (0.004)      0.024 (0)     [0.007:0.035]
    ## 18          Dep_EpiDoc2~Age  0.001 (0.963)   0.04 (0.062)     [-0.03:0.032]
    ## 19           MM_EpiDoc2~Age      0 (0.843)  0.003 (0.447)    [-0.005:0.004]
    ## 20          Dep_EpiDoc1~Age  0.041 (0.028)  0.055 (0.006)     [0.004:0.078]
    ## 21           MM_EpiDoc1~Age      0.031 (0)      0.049 (0)     [0.021:0.041]
    ## 22    Dep_EpiDoc4~Education  0.167 (0.297)  0.158 (0.346)    [-0.147:0.482]
    ## 23    Dep_EpiDoc2~Education      0.536 (0)   0.358 (0.08)      [0.252:0.82]
    ## 24    Dep_EpiDoc1~Education  0.495 (0.003)      0.777 (0)     [0.173:0.818]
    ## 25     MM_EpiDoc4~Education  -0.012 (0.84)  0.049 (0.422)      [-0.123:0.1]
    ## 26     MM_EpiDoc2~Education -0.045 (0.133)  0.022 (0.506)    [-0.103:0.014]
    ## 27     MM_EpiDoc1~Education  0.031 (0.565)  0.025 (0.613)    [-0.075:0.137]
    ## 28 Dep_EpiDoc4~~Dep_EpiDoc4      6.884 (0)      9.348 (0)     [4.942:8.827]
    ## 29   MM_EpiDoc4~~MM_EpiDoc4      0.978 (0)      1.008 (0)      [0.77:1.187]
    ## 30 Dep_EpiDoc2~~Dep_EpiDoc2      6.191 (0)     10.668 (0)     [4.477:7.905]
    ## 31   MM_EpiDoc2~~MM_EpiDoc2       0.19 (0)      0.427 (0)     [0.089:0.291]
    ## 32 Dep_EpiDoc1~~Dep_EpiDoc1      8.427 (0)     12.629 (0)    [6.155:10.698]
    ## 33   MM_EpiDoc1~~MM_EpiDoc1        0.7 (0)      0.781 (0)     [0.531:0.869]
    ## 34                 Age~~Age   204.021 (NA)   135.726 (NA) [204.021:204.021]
    ## 35           Age~~Education     6.794 (NA)      3.68 (NA)     [6.794:6.794]
    ## 36     Education~~Education     1.159 (NA)      1.24 (NA)     [1.159:1.159]
    ## 37            Dep_EpiDoc4~1  1.075 (0.004)   0.965 (0.01)     [0.347:1.803]
    ## 38             MM_EpiDoc4~1 -0.409 (0.131) -0.445 (0.044)     [-0.94:0.122]
    ## 39            Dep_EpiDoc2~1 -0.012 (0.984) -0.539 (0.518)    [-1.268:1.243]
    ## 40             MM_EpiDoc2~1  0.096 (0.096) -0.128 (0.453)    [-0.017:0.209]
    ## 41            Dep_EpiDoc1~1 -0.575 (0.424) -0.568 (0.476)    [-1.986:0.836]
    ## 42             MM_EpiDoc1~1     -0.818 (0)     -1.396 (0)    [-1.19:-0.446]
    ## 43                    Age~1    44.188 (NA)    42.481 (NA)   [44.188:44.188]
    ## 44              Education~1     2.823 (NA)     2.409 (NA)     [2.823:2.823]
    ##                 CI_2
    ## 1        [0.2:0.478]
    ## 2      [0.117:0.339]
    ## 3      [0.094:0.464]
    ## 4      [0.349:0.598]
    ## 5        [0.288:0.5]
    ## 6      [0.068:0.279]
    ## 7      [0.066:1.263]
    ## 8       [-0.748:0.1]
    ## 9     [-0.033:0.033]
    ## 10    [-0.019:0.051]
    ## 11     [-0.65:0.428]
    ## 12     [-0.007:0.05]
    ## 13             [0:0]
    ## 14             [0:0]
    ## 15      [-0.4:0.347]
    ## 16             [0:0]
    ## 17     [0.012:0.036]
    ## 18    [-0.002:0.083]
    ## 19    [-0.005:0.011]
    ## 20     [0.016:0.095]
    ## 21     [0.039:0.059]
    ## 22    [-0.171:0.488]
    ## 23    [-0.043:0.759]
    ## 24      [0.36:1.194]
    ## 25     [-0.07:0.167]
    ## 26    [-0.043:0.086]
    ## 27    [-0.073:0.123]
    ## 28    [7.292:11.404]
    ## 29     [0.838:1.177]
    ## 30    [8.068:13.268]
    ## 31     [0.282:0.571]
    ## 32    [9.237:16.022]
    ## 33     [0.544:1.017]
    ## 34 [135.726:135.726]
    ## 35       [3.68:3.68]
    ## 36       [1.24:1.24]
    ## 37     [0.232:1.697]
    ## 38   [-0.877:-0.012]
    ## 39      [-2.177:1.1]
    ## 40    [-0.464:0.207]
    ## 41     [-2.13:0.994]
    ## 42   [-1.839:-0.954]
    ## 43   [42.481:42.481]
    ## 44     [2.409:2.409]

## model with restriction on AR and Bi-dir

``` r
Model_AR_BiD<-'
# evolução
# evolução
#Dep_EpiDoc4~dd03*Dep_EpiDoc1+dd13*Dep_EpiDoc2
Dep_EpiDoc4~Dep_EpiDoc1+dd13*Dep_EpiDoc2
MM_EpiDoc4~mm14*MM_EpiDoc2+mm04*MM_EpiDoc1
Dep_EpiDoc2~dd01*Dep_EpiDoc1
MM_EpiDoc2~mm01*MM_EpiDoc1 

# causalidade - c/ temporalidade/lag
Dep_EpiDoc4~MM_EpiDoc2+0*MM_EpiDoc1
Dep_EpiDoc2~dm21*MM_EpiDoc1
MM_EpiDoc4~0*Dep_EpiDoc2+Dep_EpiDoc1 
MM_EpiDoc2~Dep_EpiDoc1

# corrrelação mesmo t
MM_EpiDoc1~~0*Dep_EpiDoc1
MM_EpiDoc2~~0* Dep_EpiDoc2
MM_EpiDoc4~~Dep_EpiDoc4


# causalidade - idade
Dep_EpiDoc4~Age
MM_EpiDoc4~Age
Dep_EpiDoc2~Age
MM_EpiDoc2~Age
Dep_EpiDoc1~Age
MM_EpiDoc1~Age

# Escolaridade
Dep_EpiDoc4~Education
Dep_EpiDoc2~Education
Dep_EpiDoc1~Education
MM_EpiDoc4+MM_EpiDoc2+MM_EpiDoc1~Education

'
```

``` r
resAR_BiD<-lavaan.mi::sem.mi(data=data.mi,
         model=Model_AR_BiD,
         estimator = "ML", test = "yuan.bentler.mplus",
         likelihood = "wishart",
           sampling.weights='ipw',group=GVAr,
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
summary(resAR_BiD,fit.meas=TRUE,standardized=TRUE,rsquare=TRUE,ci=TRUE)
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
    ##   Number of model parameters                        70
    ##   Number of equality constraints                     6
    ## 
    ##   Number of observations per group:                   
    ##     Male                                           355
    ##     Female                                         317
    ##   Sampling weights variable                        ipw
    ## 
    ## Model Test User Model:
    ## 
    ##                                                     Standard      Scaled
    ##   Test statistic                                      70.543      45.345
    ##   Degrees of freedom                                      14          14
    ##   P-value                                              0.000       0.000
    ##   Average scaling correction factor                                1.556
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
    ##   Comparative Fit Index (CFI)                    0.943       0.960
    ##   Tucker-Lewis Index (TLI)                       0.781       0.847
    ##                                                                   
    ##   Robust Comparative Fit Index (CFI)                         0.950
    ##   Robust Tucker-Lewis Index (TLI)                            0.808
    ## 
    ## Loglikelihood and Information Criteria:
    ## 
    ##   Loglikelihood user model (H0)              -7363.705   -7363.705
    ##   Scaling correction factor                                  1.618
    ##       for the MLR correction                                      
    ##   Loglikelihood unrestricted model (H1)      -7325.594   -7325.594
    ##   Scaling correction factor                                  1.717
    ##       for the MLR correction                                      
    ##                                                                   
    ##   Akaike (AIC)                               14855.409   14855.409
    ##   Bayesian (BIC)                             15144.066   15144.066
    ##   Sample-size adjusted Bayesian (SABIC)      14940.861   14940.861
    ## 
    ## Root Mean Square Error of Approximation:
    ## 
    ##   RMSEA                                          0.110       0.085
    ##   90 Percent confidence interval - lower         0.085       0.064
    ##   90 Percent confidence interval - upper         0.136       0.107
    ##   P-value H_0: RMSEA <= 0.050                    0.000       0.004
    ##   P-value H_0: RMSEA >= 0.080                    0.976       0.668
    ##                                                                   
    ##   Robust RMSEA                                               0.102
    ##   90 Percent confidence interval - lower                     0.070
    ##   90 Percent confidence interval - upper                     0.136
    ##   P-value H_0: Robust RMSEA <= 0.050                         0.006
    ##   P-value H_0: Robust RMSEA >= 0.080                         0.876
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
    ##     Dp_EpD1           0.159    0.074    2.142  521.136    0.033    0.013
    ##     Dp_EpD2 (dd13)    0.225    0.056    3.991  560.655    0.000    0.114
    ##   MM_EpiDoc4 ~                                                          
    ##     MM_EpD2 (mm14)    0.281    0.094    2.983      Inf    0.003    0.096
    ##     MM_EpD1 (mm04)    0.479    0.063    7.610      Inf    0.000    0.356
    ##   Dep_EpiDoc2 ~                                                         
    ##     Dp_EpD1 (dd01)    0.391    0.054    7.228  345.812    0.000    0.285
    ##   MM_EpiDoc2 ~                                                          
    ##     MM_EpD1 (mm01)    0.173    0.054    3.240      Inf    0.001    0.068
    ##   Dep_EpiDoc4 ~                                                         
    ##     MM_EpD2           0.458    0.299    1.534      Inf    0.125   -0.127
    ##     MM_EpD1           0.000                                        0.000
    ##   Dep_EpiDoc2 ~                                                         
    ##     MM_EpD1 (dm21)    0.163    0.176    0.928  359.489    0.354   -0.183
    ##   MM_EpiDoc4 ~                                                          
    ##     Dp_EpD2           0.000                                        0.000
    ##     Dp_EpD1           0.028    0.025    1.121      Inf    0.262   -0.021
    ##   MM_EpiDoc2 ~                                                          
    ##     Dp_EpD1           0.016    0.015    1.087      Inf    0.277   -0.013
    ##   Dep_EpiDoc4 ~                                                         
    ##     Age               0.019    0.013    1.431  157.085    0.155   -0.007
    ##   MM_EpiDoc4 ~                                                          
    ##     Age               0.021    0.007    2.893      Inf    0.004    0.007
    ##   Dep_EpiDoc2 ~                                                         
    ##     Age               0.005    0.015    0.320   83.668    0.750   -0.025
    ##   MM_EpiDoc2 ~                                                          
    ##     Age              -0.000    0.002   -0.199      Inf    0.842   -0.005
    ##   Dep_EpiDoc1 ~                                                         
    ##     Age               0.041    0.019    2.200      Inf    0.028    0.004
    ##   MM_EpiDoc1 ~                                                          
    ##     Age               0.031    0.005    6.114      Inf    0.000    0.021
    ##   Dep_EpiDoc4 ~                                                         
    ##     Educatn           0.072    0.154    0.465      Inf    0.642   -0.231
    ##   Dep_EpiDoc2 ~                                                         
    ##     Educatn           0.541    0.144    3.758  491.648    0.000    0.258
    ##   Dep_EpiDoc1 ~                                                         
    ##     Educatn           0.495    0.164    3.021      Inf    0.003    0.174
    ##   MM_EpiDoc4 ~                                                          
    ##     Educatn           0.007    0.060    0.113      Inf    0.910   -0.111
    ##   MM_EpiDoc2 ~                                                          
    ##     Educatn          -0.045    0.030   -1.505      Inf    0.132   -0.103
    ##   MM_EpiDoc1 ~                                                          
    ##     Educatn           0.031    0.054    0.578      Inf    0.564   -0.075
    ##  ci.upper   Std.lv  Std.all
    ##                            
    ##     0.305    0.159    0.170
    ##     0.335    0.225    0.229
    ##                            
    ##     0.466    0.281    0.106
    ##     0.602    0.479    0.369
    ##                            
    ##     0.498    0.391    0.409
    ##                            
    ##     0.278    0.173    0.354
    ##                            
    ##     1.043    0.458    0.075
    ##     0.000    0.000    0.000
    ##                            
    ##     0.509    0.163    0.053
    ##                            
    ##     0.000    0.000    0.000
    ##     0.076    0.028    0.069
    ##                            
    ##     0.045    0.016    0.105
    ##                            
    ##     0.045    0.019    0.095
    ##                            
    ##     0.036    0.021    0.248
    ##                            
    ##     0.035    0.005    0.024
    ##                            
    ##     0.004   -0.000   -0.014
    ##                            
    ##     0.077    0.041    0.191
    ##                            
    ##     0.041    0.031    0.465
    ##                            
    ##     0.374    0.072    0.027
    ##                            
    ##     0.824    0.541    0.199
    ##                            
    ##     0.817    0.495    0.175
    ##                            
    ##     0.124    0.007    0.006
    ##                            
    ##     0.014   -0.045   -0.103
    ##                            
    ##     0.137    0.031    0.035
    ## 
    ## Covariances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##  .Dep_EpiDoc1 ~~                                                        
    ##    .MM_EpiDoc1        0.000                                        0.000
    ##  .Dep_EpiDoc2 ~~                                                        
    ##    .MM_EpiDoc2        0.000                                        0.000
    ##  .Dep_EpiDoc4 ~~                                                        
    ##    .MM_EpiDoc4        0.213    0.169    1.259  191.272    0.209   -0.120
    ##  ci.upper   Std.lv  Std.all
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.546    0.213    0.082
    ## 
    ## Intercepts:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##    .Dep_EpiDoc4       0.560    0.523    1.072  152.891    0.286   -0.472
    ##    .MM_EpiDoc4       -0.422    0.280   -1.504      Inf    0.133   -0.971
    ##    .Dep_EpiDoc2      -0.119    0.618   -0.192  105.346    0.848   -1.343
    ##    .MM_EpiDoc2        0.096    0.058    1.667      Inf    0.095   -0.017
    ##    .Dep_EpiDoc1      -0.575    0.718   -0.801      Inf    0.423   -1.982
    ##    .MM_EpiDoc1       -0.818    0.189   -4.320      Inf    0.000   -1.189
    ##  ci.upper   Std.lv  Std.all
    ##     1.592    0.560    0.195
    ##     0.128   -0.422   -0.341
    ##     1.106   -0.119   -0.041
    ##     0.209    0.096    0.206
    ##     0.832   -0.575   -0.188
    ##    -0.447   -0.818   -0.856
    ## 
    ## Variances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##    .Dep_EpiDoc4       6.821    0.968    7.044      Inf    0.000    4.922
    ##    .MM_EpiDoc4        0.989    0.108    9.159      Inf    0.000    0.777
    ##    .Dep_EpiDoc2       6.206    0.864    7.185  103.684    0.000    4.493
    ##    .MM_EpiDoc2        0.190    0.051    3.692      Inf    0.000    0.089
    ##    .Dep_EpiDoc1       8.427    1.156    7.290      Inf    0.000    6.161
    ##    .MM_EpiDoc1        0.700    0.086    8.150      Inf    0.000    0.532
    ##  ci.upper   Std.lv  Std.all
    ##     8.719    6.821    0.831
    ##     1.201    0.989    0.645
    ##     7.919    6.206    0.726
    ##     0.291    0.190    0.870
    ##    10.692    8.427    0.903
    ##     0.869    0.700    0.768
    ## 
    ## R-Square:
    ##                    Estimate
    ##     Dep_EpiDoc4       0.169
    ##     MM_EpiDoc4        0.355
    ##     Dep_EpiDoc2       0.274
    ##     MM_EpiDoc2        0.130
    ##     Dep_EpiDoc1       0.097
    ##     MM_EpiDoc1        0.232
    ## 
    ## 
    ## Group 2 [Female]:
    ## 
    ## Regressions:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##   Dep_EpiDoc4 ~                                                         
    ##     Dp_EpD1           0.322    0.070    4.604      Inf    0.000    0.185
    ##     Dp_EpD2 (dd13)    0.225    0.056    3.991  560.655    0.000    0.114
    ##   MM_EpiDoc4 ~                                                          
    ##     MM_EpD2 (mm14)    0.281    0.094    2.983      Inf    0.003    0.096
    ##     MM_EpD1 (mm04)    0.479    0.063    7.610      Inf    0.000    0.356
    ##   Dep_EpiDoc2 ~                                                         
    ##     Dp_EpD1 (dd01)    0.391    0.054    7.228  345.812    0.000    0.285
    ##   MM_EpiDoc2 ~                                                          
    ##     MM_EpD1 (mm01)    0.173    0.054    3.240      Inf    0.001    0.068
    ##   Dep_EpiDoc4 ~                                                         
    ##     MM_EpD2           0.497    0.303    1.643      Inf    0.100   -0.096
    ##     MM_EpD1           0.000                                        0.000
    ##   Dep_EpiDoc2 ~                                                         
    ##     MM_EpD1 (dm21)    0.163    0.176    0.928  359.489    0.354   -0.183
    ##   MM_EpiDoc4 ~                                                          
    ##     Dp_EpD2           0.000                                        0.000
    ##     Dp_EpD1           0.015    0.017    0.907      Inf    0.364   -0.018
    ##   MM_EpiDoc2 ~                                                          
    ##     Dp_EpD1           0.021    0.015    1.465      Inf    0.143   -0.007
    ##   Dep_EpiDoc4 ~                                                         
    ##     Age              -0.013    0.018   -0.729      Inf    0.466   -0.048
    ##   MM_EpiDoc4 ~                                                          
    ##     Age               0.024    0.006    3.978      Inf    0.000    0.012
    ##   Dep_EpiDoc2 ~                                                         
    ##     Age               0.027    0.020    1.383  505.162    0.167   -0.011
    ##   MM_EpiDoc2 ~                                                          
    ##     Age               0.003    0.004    0.763      Inf    0.446   -0.005
    ##   Dep_EpiDoc1 ~                                                         
    ##     Age               0.055    0.020    2.777      Inf    0.005    0.016
    ##   MM_EpiDoc1 ~                                                          
    ##     Age               0.049    0.005    9.625      Inf    0.000    0.039
    ##   Dep_EpiDoc4 ~                                                         
    ##     Educatn           0.174    0.177    0.980      Inf    0.327   -0.174
    ##   Dep_EpiDoc2 ~                                                         
    ##     Educatn           0.353    0.203    1.735  292.139    0.084   -0.047
    ##   Dep_EpiDoc1 ~                                                         
    ##     Educatn           0.777    0.212    3.665      Inf    0.000    0.361
    ##   MM_EpiDoc4 ~                                                          
    ##     Educatn           0.049    0.059    0.822      Inf    0.411   -0.068
    ##   MM_EpiDoc2 ~                                                          
    ##     Educatn           0.022    0.033    0.667      Inf    0.505   -0.042
    ##   MM_EpiDoc1 ~                                                          
    ##     Educatn           0.025    0.050    0.507      Inf    0.612   -0.073
    ##  ci.upper   Std.lv  Std.all
    ##                            
    ##     0.459    0.322    0.337
    ##     0.335    0.225    0.233
    ##                            
    ##     0.466    0.281    0.151
    ##     0.602    0.479    0.391
    ##                            
    ##     0.498    0.391    0.393
    ##                            
    ##     0.278    0.173    0.263
    ##                            
    ##     1.091    0.497    0.096
    ##     0.000    0.000    0.000
    ##                            
    ##     0.509    0.163    0.046
    ##                            
    ##     0.000    0.000    0.000
    ##     0.049    0.015    0.045
    ##                            
    ##     0.050    0.021    0.116
    ##                            
    ##     0.022   -0.013   -0.042
    ##                            
    ##     0.035    0.024    0.214
    ##                            
    ##     0.066    0.027    0.085
    ##                            
    ##     0.011    0.003    0.053
    ##                            
    ##     0.094    0.055    0.172
    ##                            
    ##     0.059    0.049    0.537
    ##                            
    ##     0.521    0.174    0.054
    ##                            
    ##     0.754    0.353    0.105
    ##                            
    ##     1.193    0.777    0.230
    ##                            
    ##     0.165    0.049    0.042
    ##                            
    ##     0.086    0.022    0.035
    ##                            
    ##     0.123    0.025    0.027
    ## 
    ## Covariances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##  .Dep_EpiDoc1 ~~                                                        
    ##    .MM_EpiDoc1        0.000                                        0.000
    ##  .Dep_EpiDoc2 ~~                                                        
    ##    .MM_EpiDoc2        0.000                                        0.000
    ##  .Dep_EpiDoc4 ~~                                                        
    ##    .MM_EpiDoc4       -0.003    0.191   -0.015      Inf    0.988   -0.377
    ##  ci.upper   Std.lv  Std.all
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.371   -0.003   -0.001
    ## 
    ## Intercepts:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##    .Dep_EpiDoc4       1.360    0.696    1.954  992.910    0.051   -0.006
    ##    .MM_EpiDoc4       -0.435    0.219   -1.987      Inf    0.047   -0.864
    ##    .Dep_EpiDoc2      -0.157    0.784   -0.200  253.073    0.841   -1.701
    ##    .MM_EpiDoc2       -0.128    0.171   -0.752      Inf    0.452   -0.463
    ##    .Dep_EpiDoc1      -0.568    0.795   -0.714      Inf    0.475   -2.126
    ##    .MM_EpiDoc1       -1.396    0.225   -6.204      Inf    0.000   -1.838
    ##  ci.upper   Std.lv  Std.all
    ##     2.727    1.360    0.378
    ##    -0.006   -0.435   -0.337
    ##     1.387   -0.157   -0.042
    ##     0.206   -0.128   -0.184
    ##     0.990   -0.568   -0.151
    ##    -0.955   -1.396   -1.325
    ## 
    ## Variances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##    .Dep_EpiDoc4       9.404    1.052    8.942      Inf    0.000    7.341
    ##    .MM_EpiDoc4        1.009    0.087   11.664      Inf    0.000    0.840
    ##    .Dep_EpiDoc2      10.736    1.338    8.021  342.074    0.000    8.104
    ##    .MM_EpiDoc2        0.427    0.073    5.806      Inf    0.000    0.283
    ##    .Dep_EpiDoc1      12.629    1.726    7.316      Inf    0.000    9.246
    ##    .MM_EpiDoc1        0.781    0.120    6.481      Inf    0.000    0.545
    ##  ci.upper   Std.lv  Std.all
    ##    11.467    9.404    0.728
    ##     1.179    1.009    0.605
    ##    13.369   10.736    0.770
    ##     0.571    0.427    0.881
    ##    16.013   12.629    0.895
    ##     1.017    0.781    0.703
    ## 
    ## R-Square:
    ##                    Estimate
    ##     Dep_EpiDoc4       0.272
    ##     MM_EpiDoc4        0.395
    ##     Dep_EpiDoc2       0.230
    ##     MM_EpiDoc2        0.119
    ##     Dep_EpiDoc1       0.105
    ##     MM_EpiDoc1        0.297

### testing if keep beign the same of sat

``` r
lavaan.mi::lavTestLRT.mi(resAR_BiD,resSatModel )
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
    ##             Df   Chisq Chisq diff Df diff Pr(>Chisq)     RIV     FMI
    ## resSatModel  2  2.3128                                              
    ## resAR_BiD   14 70.5426     42.433      12 2.8158e-05 0.20761 0.17192

``` r
lavaan.mi::lavTestLRT.mi(resAR_BiD,resAR,resRCov,resSatModel )
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
    ##             Df   Chisq Chisq diff Df diff Pr(>Chisq)     RIV      FMI
    ## resSatModel  2  2.3128                                               
    ## resRCov      5  5.5145      2.000       3    0.57247 0.21946 0.179963
    ## resAR       11 70.2953     39.513       6    0.00000 0.08258 0.076281
    ## resAR_BiD   14 70.5426      2.150       3    0.54187 0.44581 0.308346

``` r
partable<-parameterEstimates.mi(resAR_BiD,standardized =TRUE)
partable<-partable%>%mutate(
  label2=paste0(lhs,op,rhs),
  estimate_p=paste0(round(est,3),' (',round(pvalue,3),')'),
  CI=paste0('[',round(ci.lower,3),':',round(ci.upper,3),']')
)%>%select(-c(lhs,op,rhs))
gPartable<-reshape(partable%>%filter(!block==0), idvar = c('label2'),
                                                   timevar = "group", direction = "wide",sep='_')
gPartable%>%select(label2,estimate_p_1,estimate_p_2,CI_1,CI_2)
```

    ##                      label2   estimate_p_1   estimate_p_2              CI_1
    ## 1   Dep_EpiDoc4~Dep_EpiDoc1  0.159 (0.033)      0.322 (0)     [0.013:0.305]
    ## 2   Dep_EpiDoc4~Dep_EpiDoc2      0.225 (0)      0.225 (0)     [0.114:0.335]
    ## 3     MM_EpiDoc4~MM_EpiDoc2  0.281 (0.003)  0.281 (0.003)     [0.096:0.466]
    ## 4     MM_EpiDoc4~MM_EpiDoc1      0.479 (0)      0.479 (0)     [0.356:0.602]
    ## 5   Dep_EpiDoc2~Dep_EpiDoc1      0.391 (0)      0.391 (0)     [0.285:0.498]
    ## 6     MM_EpiDoc2~MM_EpiDoc1  0.173 (0.001)  0.173 (0.001)     [0.068:0.278]
    ## 7    Dep_EpiDoc4~MM_EpiDoc2  0.458 (0.125)    0.497 (0.1)    [-0.127:1.043]
    ## 8    Dep_EpiDoc4~MM_EpiDoc1         0 (NA)         0 (NA)             [0:0]
    ## 9    Dep_EpiDoc2~MM_EpiDoc1  0.163 (0.354)  0.163 (0.354)    [-0.183:0.509]
    ## 10   MM_EpiDoc4~Dep_EpiDoc2         0 (NA)         0 (NA)             [0:0]
    ## 11   MM_EpiDoc4~Dep_EpiDoc1  0.028 (0.262)  0.015 (0.364)    [-0.021:0.076]
    ## 12   MM_EpiDoc2~Dep_EpiDoc1  0.016 (0.277)  0.021 (0.143)    [-0.013:0.045]
    ## 13  Dep_EpiDoc1~~MM_EpiDoc1         0 (NA)         0 (NA)             [0:0]
    ## 14  Dep_EpiDoc2~~MM_EpiDoc2         0 (NA)         0 (NA)             [0:0]
    ## 15  Dep_EpiDoc4~~MM_EpiDoc4  0.213 (0.209) -0.003 (0.988)     [-0.12:0.546]
    ## 16          Dep_EpiDoc4~Age  0.019 (0.155) -0.013 (0.466)    [-0.007:0.045]
    ## 17           MM_EpiDoc4~Age  0.021 (0.004)      0.024 (0)     [0.007:0.036]
    ## 18          Dep_EpiDoc2~Age   0.005 (0.75)  0.027 (0.167)    [-0.025:0.035]
    ## 19           MM_EpiDoc2~Age      0 (0.842)  0.003 (0.446)    [-0.005:0.004]
    ## 20          Dep_EpiDoc1~Age  0.041 (0.028)  0.055 (0.005)     [0.004:0.077]
    ## 21           MM_EpiDoc1~Age      0.031 (0)      0.049 (0)     [0.021:0.041]
    ## 22    Dep_EpiDoc4~Education  0.072 (0.642)  0.174 (0.327)    [-0.231:0.374]
    ## 23    Dep_EpiDoc2~Education      0.541 (0)  0.353 (0.084)     [0.258:0.824]
    ## 24    Dep_EpiDoc1~Education  0.495 (0.003)      0.777 (0)     [0.174:0.817]
    ## 25     MM_EpiDoc4~Education   0.007 (0.91)  0.049 (0.411)    [-0.111:0.124]
    ## 26     MM_EpiDoc2~Education -0.045 (0.132)  0.022 (0.505)    [-0.103:0.014]
    ## 27     MM_EpiDoc1~Education  0.031 (0.564)  0.025 (0.612)    [-0.075:0.137]
    ## 28 Dep_EpiDoc4~~Dep_EpiDoc4      6.821 (0)      9.404 (0)     [4.922:8.719]
    ## 29   MM_EpiDoc4~~MM_EpiDoc4      0.989 (0)      1.009 (0)     [0.777:1.201]
    ## 30 Dep_EpiDoc2~~Dep_EpiDoc2      6.206 (0)     10.736 (0)     [4.493:7.919]
    ## 31   MM_EpiDoc2~~MM_EpiDoc2       0.19 (0)      0.427 (0)     [0.089:0.291]
    ## 32 Dep_EpiDoc1~~Dep_EpiDoc1      8.427 (0)     12.629 (0)    [6.161:10.692]
    ## 33   MM_EpiDoc1~~MM_EpiDoc1        0.7 (0)      0.781 (0)     [0.532:0.869]
    ## 34                 Age~~Age   204.021 (NA)   135.726 (NA) [204.021:204.021]
    ## 35           Age~~Education     6.794 (NA)      3.68 (NA)     [6.794:6.794]
    ## 36     Education~~Education     1.159 (NA)      1.24 (NA)     [1.159:1.159]
    ## 37            Dep_EpiDoc4~1   0.56 (0.286)   1.36 (0.051)    [-0.472:1.592]
    ## 38             MM_EpiDoc4~1 -0.422 (0.133) -0.435 (0.047)    [-0.971:0.128]
    ## 39            Dep_EpiDoc2~1 -0.119 (0.848) -0.157 (0.841)    [-1.343:1.106]
    ## 40             MM_EpiDoc2~1  0.096 (0.095) -0.128 (0.452)    [-0.017:0.209]
    ## 41            Dep_EpiDoc1~1 -0.575 (0.423) -0.568 (0.475)    [-1.982:0.832]
    ## 42             MM_EpiDoc1~1     -0.818 (0)     -1.396 (0)   [-1.189:-0.447]
    ## 43                    Age~1    44.188 (NA)    42.481 (NA)   [44.188:44.188]
    ## 44              Education~1     2.823 (NA)     2.409 (NA)     [2.823:2.823]
    ##                 CI_2
    ## 1      [0.185:0.459]
    ## 2      [0.114:0.335]
    ## 3      [0.096:0.466]
    ## 4      [0.356:0.602]
    ## 5      [0.285:0.498]
    ## 6      [0.068:0.278]
    ## 7     [-0.096:1.091]
    ## 8              [0:0]
    ## 9     [-0.183:0.509]
    ## 10             [0:0]
    ## 11    [-0.018:0.049]
    ## 12     [-0.007:0.05]
    ## 13             [0:0]
    ## 14             [0:0]
    ## 15    [-0.377:0.371]
    ## 16    [-0.048:0.022]
    ## 17     [0.012:0.035]
    ## 18    [-0.011:0.066]
    ## 19    [-0.005:0.011]
    ## 20     [0.016:0.094]
    ## 21     [0.039:0.059]
    ## 22    [-0.174:0.521]
    ## 23    [-0.047:0.754]
    ## 24     [0.361:1.193]
    ## 25    [-0.068:0.165]
    ## 26    [-0.042:0.086]
    ## 27    [-0.073:0.123]
    ## 28    [7.341:11.467]
    ## 29      [0.84:1.179]
    ## 30    [8.104:13.369]
    ## 31     [0.283:0.571]
    ## 32    [9.246:16.013]
    ## 33     [0.545:1.017]
    ## 34 [135.726:135.726]
    ## 35       [3.68:3.68]
    ## 36       [1.24:1.24]
    ## 37    [-0.006:2.727]
    ## 38   [-0.864:-0.006]
    ## 39    [-1.701:1.387]
    ## 40    [-0.463:0.206]
    ## 41     [-2.126:0.99]
    ## 42   [-1.838:-0.955]
    ## 43   [42.481:42.481]
    ## 44     [2.409:2.409]

## model with restriction on AR, bi-dir (part) and age

``` r
Model_ARAge<-'
# evolução
#Dep_EpiDoc4~dd03*Dep_EpiDoc1+dd13*Dep_EpiDoc2
Dep_EpiDoc4~Dep_EpiDoc1+dd13*Dep_EpiDoc2
MM_EpiDoc4~mm14*MM_EpiDoc2+mm04*MM_EpiDoc1
Dep_EpiDoc2~dd01*Dep_EpiDoc1
MM_EpiDoc2~mm01*MM_EpiDoc1 

# causalidade - c/ temporalidade/lag
Dep_EpiDoc4~MM_EpiDoc2+0*MM_EpiDoc1
Dep_EpiDoc2~dm21*MM_EpiDoc1
MM_EpiDoc4~0*Dep_EpiDoc2+Dep_EpiDoc1 
MM_EpiDoc2~Dep_EpiDoc1

# corrrelação mesmo t
MM_EpiDoc1~~0*Dep_EpiDoc1
MM_EpiDoc2~~0* Dep_EpiDoc2
MM_EpiDoc4~~Dep_EpiDoc4


# causalidade - idade
Dep_EpiDoc4~0*Age
MM_EpiDoc4~ma4*Age
Dep_EpiDoc2~da2*Age
MM_EpiDoc2~ma2*Age
Dep_EpiDoc1~da1*Age
MM_EpiDoc1~ma1*Age

# Escolaridade
Dep_EpiDoc4~Education
Dep_EpiDoc2~Education
Dep_EpiDoc1~Education
MM_EpiDoc4+MM_EpiDoc2+MM_EpiDoc1~Education

'
```

``` r
resARAge<-lavaan.mi::sem.mi(data=data.mi,
         model=Model_ARAge,
         estimator = "ML", test = "yuan.bentler.mplus",
         likelihood = "wishart",
           sampling.weights='ipw',group=GVAr,
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
summary(resARAge,fit.meas=TRUE,standardized=TRUE,rsquare=TRUE,ci=TRUE)
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
    ##   Number of equality constraints                    11
    ## 
    ##   Number of observations per group:                   
    ##     Male                                           355
    ##     Female                                         317
    ##   Sampling weights variable                        ipw
    ## 
    ## Model Test User Model:
    ## 
    ##                                                     Standard      Scaled
    ##   Test statistic                                      84.381      57.976
    ##   Degrees of freedom                                      21          21
    ##   P-value                                              0.000       0.000
    ##   Average scaling correction factor                                1.455
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
    ##   Comparative Fit Index (CFI)                    0.936       0.953
    ##   Tucker-Lewis Index (TLI)                       0.836       0.880
    ##                                                                   
    ##   Robust Comparative Fit Index (CFI)                         0.945
    ##   Robust Tucker-Lewis Index (TLI)                            0.859
    ## 
    ## Loglikelihood and Information Criteria:
    ## 
    ##   Loglikelihood user model (H0)              -7368.826   -7368.826
    ##   Scaling correction factor                                  1.540
    ##       for the MLR correction                                      
    ##   Loglikelihood unrestricted model (H1)      -7325.594   -7325.594
    ##   Scaling correction factor                                  1.717
    ##       for the MLR correction                                      
    ##                                                                   
    ##   Akaike (AIC)                               14851.652   14851.652
    ##   Bayesian (BIC)                             15108.737   15108.737
    ##   Sample-size adjusted Bayesian (SABIC)      14927.757   14927.757
    ## 
    ## Root Mean Square Error of Approximation:
    ## 
    ##   RMSEA                                          0.095       0.075
    ##   90 Percent confidence interval - lower         0.074       0.057
    ##   90 Percent confidence interval - upper         0.117       0.094
    ##   P-value H_0: RMSEA <= 0.050                    0.000       0.014
    ##   P-value H_0: RMSEA >= 0.080                    0.886       0.350
    ##                                                                   
    ##   Robust RMSEA                                               0.087
    ##   90 Percent confidence interval - lower                     0.061
    ##   90 Percent confidence interval - upper                     0.115
    ##   P-value H_0: Robust RMSEA <= 0.050                         0.012
    ##   P-value H_0: Robust RMSEA >= 0.080                         0.700
    ## 
    ## Standardized Root Mean Square Residual:
    ## 
    ##   SRMR                                           0.062       0.062
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
    ##     Dp_EpD1           0.172    0.072    2.409  670.260    0.016    0.032
    ##     Dp_EpD2 (dd13)    0.223    0.056    3.963  563.076    0.000    0.113
    ##   MM_EpiDoc4 ~                                                          
    ##     MM_EpD2 (mm14)    0.283    0.094    2.999      Inf    0.003    0.098
    ##     MM_EpD1 (mm04)    0.481    0.062    7.711      Inf    0.000    0.359
    ##   Dep_EpiDoc2 ~                                                         
    ##     Dp_EpD1 (dd01)    0.391    0.054    7.234  349.187    0.000    0.285
    ##   MM_EpiDoc2 ~                                                          
    ##     MM_EpD1 (mm01)    0.176    0.053    3.331      Inf    0.001    0.072
    ##   Dep_EpiDoc4 ~                                                         
    ##     MM_EpD2           0.528    0.289    1.825      Inf    0.068   -0.039
    ##     MM_EpD1           0.000                                        0.000
    ##   Dep_EpiDoc2 ~                                                         
    ##     MM_EpD1 (dm21)    0.179    0.175    1.022  308.896    0.307   -0.165
    ##   MM_EpiDoc4 ~                                                          
    ##     Dp_EpD2           0.000                                        0.000
    ##     Dp_EpD1           0.027    0.024    1.137      Inf    0.255   -0.020
    ##   MM_EpiDoc2 ~                                                          
    ##     Dp_EpD1           0.015    0.015    1.043      Inf    0.297   -0.013
    ##   Dep_EpiDoc4 ~                                                         
    ##     Age               0.000                                        0.000
    ##   MM_EpiDoc4 ~                                                          
    ##     Age      (ma4)    0.022    0.006    3.943      Inf    0.000    0.011
    ##   Dep_EpiDoc2 ~                                                         
    ##     Age      (da2)    0.010    0.013    0.757   91.611    0.451   -0.016
    ##   MM_EpiDoc2 ~                                                          
    ##     Age      (ma2)    0.000    0.002    0.110      Inf    0.912   -0.004
    ##   Dep_EpiDoc1 ~                                                         
    ##     Age      (da1)    0.045    0.015    3.102      Inf    0.002    0.017
    ##   MM_EpiDoc1 ~                                                          
    ##     Age      (ma1)    0.037    0.004    9.289      Inf    0.000    0.029
    ##   Dep_EpiDoc4 ~                                                         
    ##     Educatn           0.175    0.152    1.151  419.417    0.250   -0.124
    ##   Dep_EpiDoc2 ~                                                         
    ##     Educatn           0.509    0.141    3.611  344.794    0.000    0.232
    ##   Dep_EpiDoc1 ~                                                         
    ##     Educatn           0.472    0.158    2.981      Inf    0.003    0.162
    ##   MM_EpiDoc4 ~                                                          
    ##     Educatn           0.004    0.061    0.071      Inf    0.943   -0.115
    ##   MM_EpiDoc2 ~                                                          
    ##     Educatn          -0.049    0.028   -1.750      Inf    0.080   -0.104
    ##   MM_EpiDoc1 ~                                                          
    ##     Educatn          -0.004    0.051   -0.076      Inf    0.940   -0.104
    ##  ci.upper   Std.lv  Std.all
    ##                            
    ##     0.313    0.172    0.184
    ##     0.334    0.223    0.229
    ##                            
    ##     0.467    0.283    0.106
    ##     0.603    0.481    0.378
    ##                            
    ##     0.497    0.391    0.408
    ##                            
    ##     0.279    0.176    0.369
    ##                            
    ##     1.096    0.528    0.087
    ##     0.000    0.000    0.000
    ##                            
    ##     0.523    0.179    0.060
    ##                            
    ##     0.000    0.000    0.000
    ##     0.074    0.027    0.066
    ##                            
    ##     0.044    0.015    0.099
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.033    0.022    0.248
    ##                            
    ##     0.036    0.010    0.048
    ##                            
    ##     0.004    0.000    0.007
    ##                            
    ##     0.073    0.045    0.210
    ##                            
    ##     0.045    0.037    0.534
    ##                            
    ##     0.473    0.175    0.066
    ##                            
    ##     0.785    0.509    0.186
    ##                            
    ##     0.782    0.472    0.166
    ##                            
    ##     0.123    0.004    0.004
    ##                            
    ##     0.006   -0.049   -0.111
    ##                            
    ##     0.097   -0.004   -0.004
    ## 
    ## Covariances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##  .Dep_EpiDoc1 ~~                                                        
    ##    .MM_EpiDoc1        0.000                                        0.000
    ##  .Dep_EpiDoc2 ~~                                                        
    ##    .MM_EpiDoc2        0.000                                        0.000
    ##  .Dep_EpiDoc4 ~~                                                        
    ##    .MM_EpiDoc4        0.210    0.170    1.236  192.561    0.218   -0.125
    ##  ci.upper   Std.lv  Std.all
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.546    0.210    0.081
    ## 
    ## Intercepts:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##    .Dep_EpiDoc4       1.074    0.367    2.927  839.244    0.004    0.354
    ##    .MM_EpiDoc4       -0.433    0.234   -1.847      Inf    0.065   -0.893
    ##    .Dep_EpiDoc2      -0.258    0.567   -0.456  122.957    0.649   -1.380
    ##    .MM_EpiDoc2        0.077    0.057    1.363      Inf    0.173   -0.034
    ##    .Dep_EpiDoc1      -0.687    0.622   -1.104      Inf    0.270   -1.907
    ##    .MM_EpiDoc1       -0.983    0.174   -5.659      Inf    0.000   -1.324
    ##  ci.upper   Std.lv  Std.all
    ##     1.794    1.074    0.374
    ##     0.027   -0.433   -0.343
    ##     0.863   -0.258   -0.088
    ##     0.189    0.077    0.164
    ##     0.533   -0.687   -0.224
    ##    -0.643   -0.983   -0.991
    ## 
    ## Variances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##    .Dep_EpiDoc4       6.889    0.981    7.021      Inf    0.000    4.965
    ##    .MM_EpiDoc4        0.989    0.108    9.192      Inf    0.000    0.778
    ##    .Dep_EpiDoc2       6.210    0.862    7.205  104.192    0.000    4.501
    ##    .MM_EpiDoc2        0.190    0.051    3.721      Inf    0.000    0.090
    ##    .Dep_EpiDoc1       8.430    1.154    7.306      Inf    0.000    6.168
    ##    .MM_EpiDoc1        0.706    0.086    8.208      Inf    0.000    0.538
    ##  ci.upper   Std.lv  Std.all
    ##     8.813    6.889    0.838
    ##     1.200    0.989    0.622
    ##     7.919    6.210    0.719
    ##     0.291    0.190    0.853
    ##    10.691    8.430    0.898
    ##     0.875    0.706    0.717
    ## 
    ## R-Square:
    ##                    Estimate
    ##     Dep_EpiDoc4       0.162
    ##     MM_EpiDoc4        0.378
    ##     Dep_EpiDoc2       0.281
    ##     MM_EpiDoc2        0.147
    ##     Dep_EpiDoc1       0.102
    ##     MM_EpiDoc1        0.283
    ## 
    ## 
    ## Group 2 [Female]:
    ## 
    ## Regressions:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##   Dep_EpiDoc4 ~                                                         
    ##     Dp_EpD1           0.317    0.069    4.579      Inf    0.000    0.181
    ##     Dp_EpD2 (dd13)    0.223    0.056    3.963  563.076    0.000    0.113
    ##   MM_EpiDoc4 ~                                                          
    ##     MM_EpD2 (mm14)    0.283    0.094    2.999      Inf    0.003    0.098
    ##     MM_EpD1 (mm04)    0.481    0.062    7.711      Inf    0.000    0.359
    ##   Dep_EpiDoc2 ~                                                         
    ##     Dp_EpD1 (dd01)    0.391    0.054    7.234  349.187    0.000    0.285
    ##   MM_EpiDoc2 ~                                                          
    ##     MM_EpD1 (mm01)    0.176    0.053    3.331      Inf    0.001    0.072
    ##   Dep_EpiDoc4 ~                                                         
    ##     MM_EpD2           0.464    0.298    1.554      Inf    0.120   -0.121
    ##     MM_EpD1           0.000                                        0.000
    ##   Dep_EpiDoc2 ~                                                         
    ##     MM_EpD1 (dm21)    0.179    0.175    1.022  308.896    0.307   -0.165
    ##   MM_EpiDoc4 ~                                                          
    ##     Dp_EpD2           0.000                                        0.000
    ##     Dp_EpD1           0.016    0.017    0.967      Inf    0.334   -0.017
    ##   MM_EpiDoc2 ~                                                          
    ##     Dp_EpD1           0.023    0.015    1.554      Inf    0.120   -0.006
    ##   Dep_EpiDoc4 ~                                                         
    ##     Age               0.000                                        0.000
    ##   MM_EpiDoc4 ~                                                          
    ##     Age      (ma4)    0.022    0.006    3.943      Inf    0.000    0.011
    ##   Dep_EpiDoc2 ~                                                         
    ##     Age      (da2)    0.010    0.013    0.757   91.611    0.451   -0.016
    ##   MM_EpiDoc2 ~                                                          
    ##     Age      (ma2)    0.000    0.002    0.110      Inf    0.912   -0.004
    ##   Dep_EpiDoc1 ~                                                         
    ##     Age      (da1)    0.045    0.015    3.102      Inf    0.002    0.017
    ##   MM_EpiDoc1 ~                                                          
    ##     Age      (ma1)    0.037    0.004    9.289      Inf    0.000    0.029
    ##   Dep_EpiDoc4 ~                                                         
    ##     Educatn           0.144    0.167    0.859      Inf    0.390   -0.184
    ##   Dep_EpiDoc2 ~                                                         
    ##     Educatn           0.402    0.198    2.028  223.996    0.044    0.011
    ##   Dep_EpiDoc1 ~                                                         
    ##     Educatn           0.808    0.202    3.995      Inf    0.000    0.412
    ##   MM_EpiDoc4 ~                                                          
    ##     Educatn           0.053    0.059    0.900      Inf    0.368   -0.063
    ##   MM_EpiDoc2 ~                                                          
    ##     Educatn           0.029    0.035    0.827      Inf    0.408   -0.040
    ##   MM_EpiDoc1 ~                                                          
    ##     Educatn           0.060    0.056    1.069      Inf    0.285   -0.050
    ##  ci.upper   Std.lv  Std.all
    ##                            
    ##     0.453    0.317    0.331
    ##     0.334    0.223    0.231
    ##                            
    ##     0.467    0.283    0.155
    ##     0.603    0.481    0.384
    ##                            
    ##     0.497    0.391    0.395
    ##                            
    ##     0.279    0.176    0.256
    ##                            
    ##     1.048    0.464    0.089
    ##     0.000    0.000    0.000
    ##                            
    ##     0.523    0.179    0.048
    ##                            
    ##     0.000    0.000    0.000
    ##     0.049    0.016    0.048
    ##                            
    ##     0.051    0.023    0.124
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.033    0.022    0.203
    ##                            
    ##     0.036    0.010    0.031
    ##                            
    ##     0.004    0.000    0.004
    ##                            
    ##     0.073    0.045    0.140
    ##                            
    ##     0.045    0.037    0.431
    ##                            
    ##     0.471    0.144    0.045
    ##                            
    ##     0.792    0.402    0.121
    ##                            
    ##     1.204    0.808    0.240
    ##                            
    ##     0.169    0.053    0.047
    ##                            
    ##     0.098    0.029    0.047
    ##                            
    ##     0.169    0.060    0.066
    ## 
    ## Covariances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##  .Dep_EpiDoc1 ~~                                                        
    ##    .MM_EpiDoc1        0.000                                        0.000
    ##  .Dep_EpiDoc2 ~~                                                        
    ##    .MM_EpiDoc2        0.000                                        0.000
    ##  .Dep_EpiDoc4 ~~                                                        
    ##    .MM_EpiDoc4       -0.006    0.190   -0.033      Inf    0.973   -0.379
    ##  ci.upper   Std.lv  Std.all
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.366   -0.006   -0.002
    ## 
    ## Intercepts:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##    .Dep_EpiDoc4       0.911    0.374    2.437      Inf    0.015    0.178
    ##    .MM_EpiDoc4       -0.371    0.220   -1.687      Inf    0.092   -0.802
    ##    .Dep_EpiDoc2       0.450    0.614    0.733  171.822    0.465   -0.762
    ##    .MM_EpiDoc2       -0.026    0.103   -0.255      Inf    0.799   -0.228
    ##    .Dep_EpiDoc1      -0.201    0.686   -0.293      Inf    0.770   -1.545
    ##    .MM_EpiDoc1       -0.989    0.183   -5.416      Inf    0.000   -1.347
    ##  ci.upper   Std.lv  Std.all
    ##     1.643    0.911    0.254
    ##     0.060   -0.371   -0.296
    ##     1.662    0.450    0.121
    ##     0.176   -0.026   -0.038
    ##     1.143   -0.201   -0.054
    ##    -0.631   -0.989   -0.987
    ## 
    ## Variances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##    .Dep_EpiDoc4       9.428    1.063    8.868      Inf    0.000    7.343
    ##    .MM_EpiDoc4        1.010    0.086   11.708      Inf    0.000    0.841
    ##    .Dep_EpiDoc2      10.791    1.351    7.989  374.133    0.000    8.135
    ##    .MM_EpiDoc2        0.427    0.074    5.797      Inf    0.000    0.283
    ##    .Dep_EpiDoc1      12.643    1.726    7.325      Inf    0.000    9.260
    ##    .MM_EpiDoc1        0.797    0.125    6.358      Inf    0.000    0.552
    ##  ci.upper   Std.lv  Std.all
    ##    11.514    9.428    0.732
    ##     1.179    1.010    0.643
    ##    13.447   10.791    0.785
    ##     0.571    0.427    0.901
    ##    16.026   12.643    0.903
    ##     1.043    0.797    0.794
    ## 
    ## R-Square:
    ##                    Estimate
    ##     Dep_EpiDoc4       0.268
    ##     MM_EpiDoc4        0.357
    ##     Dep_EpiDoc2       0.215
    ##     MM_EpiDoc2        0.099
    ##     Dep_EpiDoc1       0.097
    ##     MM_EpiDoc1        0.206

``` r
lavaan.mi::lavTestLRT.mi(resARAge,resSatModel )
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
    ##             Df   Chisq Chisq diff Df diff Pr(>Chisq)    RIV     FMI
    ## resSatModel  2  2.3128                                             
    ## resARAge    21 84.3807     55.489      19 1.9558e-05 0.2007 0.16715

``` r
lavaan.mi::lavTestLRT.mi(resARAge,resAR_BiD,resAR,resRCov,resSatModel )
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
    ##             Df   Chisq Chisq diff Df diff Pr(>Chisq)     RIV      FMI
    ## resSatModel  2  2.3128                                               
    ## resRCov      5  5.5145      2.000       3    0.57247 0.21946 0.179963
    ## resAR       11 70.2953     39.513       6    0.00000 0.08258 0.076281
    ## resAR_BiD   14 70.5426      2.150       3    0.54187 0.44581 0.308346
    ## resARAge    21 84.3807     10.611       7    0.15652 0.18887 0.158865

``` r
partable<-parameterEstimates.mi(resARAge,standardized =TRUE)
partable<-partable%>%mutate(
  label2=paste0(lhs,op,rhs),
  estimate_p=paste0(round(est,3),' (',round(pvalue,3),')'),
  CI=paste0('[',round(ci.lower,3),':',round(ci.upper,3),']')
)%>%select(-c(lhs,op,rhs))
gPartable<-reshape(partable%>%filter(!block==0), idvar = c('label2'),
                                                   timevar = "group", direction = "wide",sep='_')
gPartable%>%select(label2,estimate_p_1,estimate_p_2,CI_1,CI_2)
```

    ##                      label2   estimate_p_1   estimate_p_2              CI_1
    ## 1   Dep_EpiDoc4~Dep_EpiDoc1  0.172 (0.016)      0.317 (0)     [0.032:0.313]
    ## 2   Dep_EpiDoc4~Dep_EpiDoc2      0.223 (0)      0.223 (0)     [0.113:0.334]
    ## 3     MM_EpiDoc4~MM_EpiDoc2  0.283 (0.003)  0.283 (0.003)     [0.098:0.467]
    ## 4     MM_EpiDoc4~MM_EpiDoc1      0.481 (0)      0.481 (0)     [0.359:0.603]
    ## 5   Dep_EpiDoc2~Dep_EpiDoc1      0.391 (0)      0.391 (0)     [0.285:0.497]
    ## 6     MM_EpiDoc2~MM_EpiDoc1  0.176 (0.001)  0.176 (0.001)     [0.072:0.279]
    ## 7    Dep_EpiDoc4~MM_EpiDoc2  0.528 (0.068)   0.464 (0.12)    [-0.039:1.096]
    ## 8    Dep_EpiDoc4~MM_EpiDoc1         0 (NA)         0 (NA)             [0:0]
    ## 9    Dep_EpiDoc2~MM_EpiDoc1  0.179 (0.307)  0.179 (0.307)    [-0.165:0.523]
    ## 10   MM_EpiDoc4~Dep_EpiDoc2         0 (NA)         0 (NA)             [0:0]
    ## 11   MM_EpiDoc4~Dep_EpiDoc1  0.027 (0.255)  0.016 (0.334)     [-0.02:0.074]
    ## 12   MM_EpiDoc2~Dep_EpiDoc1  0.015 (0.297)   0.023 (0.12)    [-0.013:0.044]
    ## 13  Dep_EpiDoc1~~MM_EpiDoc1         0 (NA)         0 (NA)             [0:0]
    ## 14  Dep_EpiDoc2~~MM_EpiDoc2         0 (NA)         0 (NA)             [0:0]
    ## 15  Dep_EpiDoc4~~MM_EpiDoc4   0.21 (0.218) -0.006 (0.973)    [-0.125:0.546]
    ## 16          Dep_EpiDoc4~Age         0 (NA)         0 (NA)             [0:0]
    ## 17           MM_EpiDoc4~Age      0.022 (0)      0.022 (0)     [0.011:0.033]
    ## 18          Dep_EpiDoc2~Age   0.01 (0.451)   0.01 (0.451)    [-0.016:0.036]
    ## 19           MM_EpiDoc2~Age      0 (0.912)      0 (0.912)    [-0.004:0.004]
    ## 20          Dep_EpiDoc1~Age  0.045 (0.002)  0.045 (0.002)     [0.017:0.073]
    ## 21           MM_EpiDoc1~Age      0.037 (0)      0.037 (0)     [0.029:0.045]
    ## 22    Dep_EpiDoc4~Education   0.175 (0.25)   0.144 (0.39)    [-0.124:0.473]
    ## 23    Dep_EpiDoc2~Education      0.509 (0)  0.402 (0.044)     [0.232:0.785]
    ## 24    Dep_EpiDoc1~Education  0.472 (0.003)      0.808 (0)     [0.162:0.782]
    ## 25     MM_EpiDoc4~Education  0.004 (0.943)  0.053 (0.368)    [-0.115:0.123]
    ## 26     MM_EpiDoc2~Education  -0.049 (0.08)  0.029 (0.408)    [-0.104:0.006]
    ## 27     MM_EpiDoc1~Education  -0.004 (0.94)   0.06 (0.285)    [-0.104:0.097]
    ## 28 Dep_EpiDoc4~~Dep_EpiDoc4      6.889 (0)      9.428 (0)     [4.965:8.813]
    ## 29   MM_EpiDoc4~~MM_EpiDoc4      0.989 (0)       1.01 (0)       [0.778:1.2]
    ## 30 Dep_EpiDoc2~~Dep_EpiDoc2       6.21 (0)     10.791 (0)     [4.501:7.919]
    ## 31   MM_EpiDoc2~~MM_EpiDoc2       0.19 (0)      0.427 (0)      [0.09:0.291]
    ## 32 Dep_EpiDoc1~~Dep_EpiDoc1       8.43 (0)     12.643 (0)    [6.168:10.691]
    ## 33   MM_EpiDoc1~~MM_EpiDoc1      0.706 (0)      0.797 (0)     [0.538:0.875]
    ## 34                 Age~~Age   204.021 (NA)   135.726 (NA) [204.021:204.021]
    ## 35           Age~~Education     6.794 (NA)      3.68 (NA)     [6.794:6.794]
    ## 36     Education~~Education     1.159 (NA)      1.24 (NA)     [1.159:1.159]
    ## 37            Dep_EpiDoc4~1  1.074 (0.004)  0.911 (0.015)     [0.354:1.794]
    ## 38             MM_EpiDoc4~1 -0.433 (0.065) -0.371 (0.092)    [-0.893:0.027]
    ## 39            Dep_EpiDoc2~1 -0.258 (0.649)   0.45 (0.465)     [-1.38:0.863]
    ## 40             MM_EpiDoc2~1  0.077 (0.173) -0.026 (0.799)    [-0.034:0.189]
    ## 41            Dep_EpiDoc1~1  -0.687 (0.27)  -0.201 (0.77)    [-1.907:0.533]
    ## 42             MM_EpiDoc1~1     -0.983 (0)     -0.989 (0)   [-1.324:-0.643]
    ## 43                    Age~1    44.188 (NA)    42.481 (NA)   [44.188:44.188]
    ## 44              Education~1     2.823 (NA)     2.409 (NA)     [2.823:2.823]
    ##                 CI_2
    ## 1      [0.181:0.453]
    ## 2      [0.113:0.334]
    ## 3      [0.098:0.467]
    ## 4      [0.359:0.603]
    ## 5      [0.285:0.497]
    ## 6      [0.072:0.279]
    ## 7     [-0.121:1.048]
    ## 8              [0:0]
    ## 9     [-0.165:0.523]
    ## 10             [0:0]
    ## 11    [-0.017:0.049]
    ## 12    [-0.006:0.051]
    ## 13             [0:0]
    ## 14             [0:0]
    ## 15    [-0.379:0.366]
    ## 16             [0:0]
    ## 17     [0.011:0.033]
    ## 18    [-0.016:0.036]
    ## 19    [-0.004:0.004]
    ## 20     [0.017:0.073]
    ## 21     [0.029:0.045]
    ## 22    [-0.184:0.471]
    ## 23     [0.011:0.792]
    ## 24     [0.412:1.204]
    ## 25    [-0.063:0.169]
    ## 26     [-0.04:0.098]
    ## 27     [-0.05:0.169]
    ## 28    [7.343:11.514]
    ## 29     [0.841:1.179]
    ## 30    [8.135:13.447]
    ## 31     [0.283:0.571]
    ## 32     [9.26:16.026]
    ## 33     [0.552:1.043]
    ## 34 [135.726:135.726]
    ## 35       [3.68:3.68]
    ## 36       [1.24:1.24]
    ## 37     [0.178:1.643]
    ## 38     [-0.802:0.06]
    ## 39    [-0.762:1.662]
    ## 40    [-0.228:0.176]
    ## 41    [-1.545:1.143]
    ## 42   [-1.347:-0.631]
    ## 43   [42.481:42.481]
    ## 44     [2.409:2.409]

## model with restriction on ED

``` r
Model_ED<-'
# evolução
#Dep_EpiDoc4~dd03*Dep_EpiDoc1+dd13*Dep_EpiDoc2
Dep_EpiDoc4~Dep_EpiDoc1+dd13*Dep_EpiDoc2
MM_EpiDoc4~mm14*MM_EpiDoc2+mm04*MM_EpiDoc1
Dep_EpiDoc2~dd01*Dep_EpiDoc1
MM_EpiDoc2~mm01*MM_EpiDoc1 

# causalidade - c/ temporalidade/lag
Dep_EpiDoc4~MM_EpiDoc2+0*MM_EpiDoc1
Dep_EpiDoc2~dm21*MM_EpiDoc1
MM_EpiDoc4~0*Dep_EpiDoc2+Dep_EpiDoc1 
MM_EpiDoc2~Dep_EpiDoc1

# corrrelação mesmo t
MM_EpiDoc1~~0*Dep_EpiDoc1
MM_EpiDoc2~~0* Dep_EpiDoc2
MM_EpiDoc4~~Dep_EpiDoc4


# causalidade - idade
Dep_EpiDoc4~0*Age
MM_EpiDoc4~ma4*Age
Dep_EpiDoc2~da2*Age
MM_EpiDoc2~ma2*Age
Dep_EpiDoc1~da1*Age
MM_EpiDoc1~ma1*Age

# Escolaridade
Dep_EpiDoc4~0*Education
Dep_EpiDoc2~de2*Education
Dep_EpiDoc1~Education
MM_EpiDoc4~0*Education
MM_EpiDoc2~0*Education
MM_EpiDoc1~Education

'
```

``` r
resED<-lavaan.mi::sem.mi(data=data.mi,
         model=Model_ED,
         estimator = "ML", test = "yuan.bentler.mplus",
         likelihood = "wishart",
           sampling.weights='ipw',group=GVAr,
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
summary(resED,fit.meas=TRUE,standardized=TRUE,rsquare=TRUE,ci=TRUE)
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
    ##   Number of model parameters                        62
    ##   Number of equality constraints                    12
    ## 
    ##   Number of observations per group:                   
    ##     Male                                           355
    ##     Female                                         317
    ##   Sampling weights variable                        ipw
    ## 
    ## Model Test User Model:
    ## 
    ##                                                     Standard      Scaled
    ##   Test statistic                                      93.912      69.559
    ##   Degrees of freedom                                      28          28
    ##   P-value                                              0.000       0.000
    ##   Average scaling correction factor                                1.350
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
    ##   Comparative Fit Index (CFI)                    0.934       0.948
    ##   Tucker-Lewis Index (TLI)                       0.872       0.899
    ##                                                                   
    ##   Robust Comparative Fit Index (CFI)                         0.943
    ##   Robust Tucker-Lewis Index (TLI)                            0.890
    ## 
    ## Loglikelihood and Information Criteria:
    ## 
    ##   Loglikelihood user model (H0)              -7373.108   -7373.108
    ##   Scaling correction factor                                  1.572
    ##       for the MLR correction                                      
    ##   Loglikelihood unrestricted model (H1)      -7325.594   -7325.594
    ##   Scaling correction factor                                  1.717
    ##       for the MLR correction                                      
    ##                                                                   
    ##   Akaike (AIC)                               14846.217   14846.217
    ##   Bayesian (BIC)                             15071.730   15071.730
    ##   Sample-size adjusted Bayesian (SABIC)      14912.976   14912.976
    ## 
    ## Root Mean Square Error of Approximation:
    ## 
    ##   RMSEA                                          0.084       0.069
    ##   90 Percent confidence interval - lower         0.065       0.052
    ##   90 Percent confidence interval - upper         0.103       0.086
    ##   P-value H_0: RMSEA <= 0.050                    0.002       0.037
    ##   P-value H_0: RMSEA >= 0.080                    0.651       0.142
    ##                                                                   
    ##   Robust RMSEA                                               0.077
    ##   90 Percent confidence interval - lower                     0.055
    ##   90 Percent confidence interval - upper                     0.100
    ##   P-value H_0: Robust RMSEA <= 0.050                         0.025
    ##   P-value H_0: Robust RMSEA >= 0.080                         0.448
    ## 
    ## Standardized Root Mean Square Residual:
    ## 
    ##   SRMR                                           0.066       0.066
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
    ##     Dp_EpD1           0.184    0.070    2.613  735.664    0.009    0.046
    ##     Dp_EpD2 (dd13)    0.235    0.055    4.284  649.850    0.000    0.127
    ##   MM_EpiDoc4 ~                                                          
    ##     MM_EpD2 (mm14)    0.285    0.094    3.024      Inf    0.002    0.100
    ##     MM_EpD1 (mm04)    0.480    0.062    7.704      Inf    0.000    0.358
    ##   Dep_EpiDoc2 ~                                                         
    ##     Dp_EpD1 (dd01)    0.390    0.054    7.219  333.311    0.000    0.284
    ##   MM_EpiDoc2 ~                                                          
    ##     MM_EpD1 (mm01)    0.177    0.053    3.319      Inf    0.001    0.072
    ##   Dep_EpiDoc4 ~                                                         
    ##     MM_EpD2           0.506    0.288    1.759      Inf    0.079   -0.058
    ##     MM_EpD1           0.000                                        0.000
    ##   Dep_EpiDoc2 ~                                                         
    ##     MM_EpD1 (dm21)    0.177    0.175    1.014  307.579    0.311   -0.167
    ##   MM_EpiDoc4 ~                                                          
    ##     Dp_EpD2           0.000                                        0.000
    ##     Dp_EpD1           0.027    0.023    1.177      Inf    0.239   -0.018
    ##   MM_EpiDoc2 ~                                                          
    ##     Dp_EpD1           0.012    0.015    0.827      Inf    0.408   -0.017
    ##   Dep_EpiDoc4 ~                                                         
    ##     Age               0.000                                        0.000
    ##   MM_EpiDoc4 ~                                                          
    ##     Age      (ma4)    0.022    0.006    3.925      Inf    0.000    0.011
    ##   Dep_EpiDoc2 ~                                                         
    ##     Age      (da2)    0.010    0.013    0.793   89.737    0.430   -0.016
    ##   MM_EpiDoc2 ~                                                          
    ##     Age      (ma2)   -0.001    0.002   -0.523      Inf    0.601   -0.004
    ##   Dep_EpiDoc1 ~                                                         
    ##     Age      (da1)    0.045    0.014    3.104      Inf    0.002    0.017
    ##   MM_EpiDoc1 ~                                                          
    ##     Age      (ma1)    0.037    0.004    9.295      Inf    0.000    0.029
    ##   Dep_EpiDoc4 ~                                                         
    ##     Educatn           0.000                                        0.000
    ##   Dep_EpiDoc2 ~                                                         
    ##     Educatn  (de2)    0.473    0.119    3.982  320.921    0.000    0.239
    ##   Dep_EpiDoc1 ~                                                         
    ##     Educatn           0.472    0.158    2.983      Inf    0.003    0.162
    ##   MM_EpiDoc4 ~                                                          
    ##     Educatn           0.000                                        0.000
    ##   MM_EpiDoc2 ~                                                          
    ##     Educatn           0.000                                        0.000
    ##   MM_EpiDoc1 ~                                                          
    ##     Educatn          -0.004    0.051   -0.076      Inf    0.940   -0.104
    ##  ci.upper   Std.lv  Std.all
    ##                            
    ##     0.321    0.184    0.196
    ##     0.343    0.235    0.240
    ##                            
    ##     0.470    0.285    0.107
    ##     0.603    0.480    0.377
    ##                            
    ##     0.497    0.390    0.408
    ##                            
    ##     0.281    0.177    0.371
    ##                            
    ##     1.070    0.506    0.084
    ##     0.000    0.000    0.000
    ##                            
    ##     0.522    0.177    0.060
    ##                            
    ##     0.000    0.000    0.000
    ##     0.072    0.027    0.066
    ##                            
    ##     0.041    0.012    0.078
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.033    0.022    0.252
    ##                            
    ##     0.037    0.010    0.051
    ##                            
    ##     0.002   -0.001   -0.027
    ##                            
    ##     0.073    0.045    0.210
    ##                            
    ##     0.045    0.037    0.534
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.707    0.473    0.174
    ##                            
    ##     0.782    0.472    0.166
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.097   -0.004   -0.004
    ## 
    ## Covariances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##  .Dep_EpiDoc1 ~~                                                        
    ##    .MM_EpiDoc1        0.000                                        0.000
    ##  .Dep_EpiDoc2 ~~                                                        
    ##    .MM_EpiDoc2        0.000                                        0.000
    ##  .Dep_EpiDoc4 ~~                                                        
    ##    .MM_EpiDoc4        0.207    0.172    1.207  193.697    0.229   -0.132
    ##  ci.upper   Std.lv  Std.all
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.546    0.207    0.079
    ## 
    ## Intercepts:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##    .Dep_EpiDoc4       1.509    0.194    7.758  399.443    0.000    1.126
    ##    .MM_EpiDoc4       -0.438    0.191   -2.296      Inf    0.022   -0.811
    ##    .Dep_EpiDoc2      -0.181    0.525   -0.344  129.504    0.731   -1.220
    ##    .MM_EpiDoc2       -0.004    0.060   -0.061      Inf    0.951   -0.122
    ##    .Dep_EpiDoc1      -0.687    0.622   -1.105      Inf    0.269   -1.906
    ##    .MM_EpiDoc1       -0.983    0.174   -5.663      Inf    0.000   -1.323
    ##  ci.upper   Std.lv  Std.all
    ##     1.891    1.509    0.527
    ##    -0.064   -0.438   -0.346
    ##     0.858   -0.181   -0.062
    ##     0.114   -0.004   -0.008
    ##     0.532   -0.687   -0.224
    ##    -0.643   -0.983   -0.991
    ## 
    ## Variances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##    .Dep_EpiDoc4       6.924    0.994    6.963      Inf    0.000    4.974
    ##    .MM_EpiDoc4        0.989    0.108    9.193      Inf    0.000    0.778
    ##    .Dep_EpiDoc2       6.214    0.863    7.199  104.827    0.000    4.503
    ##    .MM_EpiDoc2        0.193    0.052    3.697      Inf    0.000    0.090
    ##    .Dep_EpiDoc1       8.430    1.153    7.311      Inf    0.000    6.170
    ##    .MM_EpiDoc1        0.706    0.086    8.214      Inf    0.000    0.538
    ##  ci.upper   Std.lv  Std.all
    ##     8.875    6.924    0.844
    ##     1.200    0.989    0.619
    ##     7.926    6.214    0.725
    ##     0.295    0.193    0.859
    ##    10.689    8.430    0.898
    ##     0.875    0.706    0.717
    ## 
    ## R-Square:
    ##                    Estimate
    ##     Dep_EpiDoc4       0.156
    ##     MM_EpiDoc4        0.381
    ##     Dep_EpiDoc2       0.275
    ##     MM_EpiDoc2        0.141
    ##     Dep_EpiDoc1       0.102
    ##     MM_EpiDoc1        0.283
    ## 
    ## 
    ## Group 2 [Female]:
    ## 
    ## Regressions:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##   Dep_EpiDoc4 ~                                                         
    ##     Dp_EpD1           0.324    0.068    4.756      Inf    0.000    0.190
    ##     Dp_EpD2 (dd13)    0.235    0.055    4.284  649.850    0.000    0.127
    ##   MM_EpiDoc4 ~                                                          
    ##     MM_EpD2 (mm14)    0.285    0.094    3.024      Inf    0.002    0.100
    ##     MM_EpD1 (mm04)    0.480    0.062    7.704      Inf    0.000    0.358
    ##   Dep_EpiDoc2 ~                                                         
    ##     Dp_EpD1 (dd01)    0.390    0.054    7.219  333.311    0.000    0.284
    ##   MM_EpiDoc2 ~                                                          
    ##     MM_EpD1 (mm01)    0.177    0.053    3.319      Inf    0.001    0.072
    ##   Dep_EpiDoc4 ~                                                         
    ##     MM_EpD2           0.475    0.298    1.595      Inf    0.111   -0.109
    ##     MM_EpD1           0.000                                        0.000
    ##   Dep_EpiDoc2 ~                                                         
    ##     MM_EpD1 (dm21)    0.177    0.175    1.014  307.579    0.311   -0.167
    ##   MM_EpiDoc4 ~                                                          
    ##     Dp_EpD2           0.000                                        0.000
    ##     Dp_EpD1           0.020    0.016    1.254      Inf    0.210   -0.011
    ##   MM_EpiDoc2 ~                                                          
    ##     Dp_EpD1           0.026    0.014    1.808      Inf    0.071   -0.002
    ##   Dep_EpiDoc4 ~                                                         
    ##     Age               0.000                                        0.000
    ##   MM_EpiDoc4 ~                                                          
    ##     Age      (ma4)    0.022    0.006    3.925      Inf    0.000    0.011
    ##   Dep_EpiDoc2 ~                                                         
    ##     Age      (da2)    0.010    0.013    0.793   89.737    0.430   -0.016
    ##   MM_EpiDoc2 ~                                                          
    ##     Age      (ma2)   -0.001    0.002   -0.523      Inf    0.601   -0.004
    ##   Dep_EpiDoc1 ~                                                         
    ##     Age      (da1)    0.045    0.014    3.104      Inf    0.002    0.017
    ##   MM_EpiDoc1 ~                                                          
    ##     Age      (ma1)    0.037    0.004    9.295      Inf    0.000    0.029
    ##   Dep_EpiDoc4 ~                                                         
    ##     Educatn           0.000                                        0.000
    ##   Dep_EpiDoc2 ~                                                         
    ##     Educatn  (de2)    0.473    0.119    3.982  320.921    0.000    0.239
    ##   Dep_EpiDoc1 ~                                                         
    ##     Educatn           0.808    0.202    3.998      Inf    0.000    0.412
    ##   MM_EpiDoc4 ~                                                          
    ##     Educatn           0.000                                        0.000
    ##   MM_EpiDoc2 ~                                                          
    ##     Educatn           0.000                                        0.000
    ##   MM_EpiDoc1 ~                                                          
    ##     Educatn           0.060    0.056    1.070      Inf    0.285   -0.050
    ##  ci.upper   Std.lv  Std.all
    ##                            
    ##     0.457    0.324    0.337
    ##     0.343    0.235    0.244
    ##                            
    ##     0.470    0.285    0.157
    ##     0.603    0.480    0.386
    ##                            
    ##     0.497    0.390    0.391
    ##                            
    ##     0.281    0.177    0.258
    ##                            
    ##     1.059    0.475    0.091
    ##     0.000    0.000    0.000
    ##                            
    ##     0.522    0.177    0.048
    ##                            
    ##     0.000    0.000    0.000
    ##     0.052    0.020    0.061
    ##                            
    ##     0.054    0.026    0.141
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.033    0.022    0.208
    ##                            
    ##     0.037    0.010    0.033
    ##                            
    ##     0.002   -0.001   -0.015
    ##                            
    ##     0.073    0.045    0.140
    ##                            
    ##     0.045    0.037    0.431
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.707    0.473    0.141
    ##                            
    ##     1.204    0.808    0.240
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.169    0.060    0.066
    ## 
    ## Covariances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##  .Dep_EpiDoc1 ~~                                                        
    ##    .MM_EpiDoc1        0.000                                        0.000
    ##  .Dep_EpiDoc2 ~~                                                        
    ##    .MM_EpiDoc2        0.000                                        0.000
    ##  .Dep_EpiDoc4 ~~                                                        
    ##    .MM_EpiDoc4        0.003    0.190    0.014      Inf    0.989   -0.369
    ##  ci.upper   Std.lv  Std.all
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.000    0.000    0.000
    ##                            
    ##     0.375    0.003    0.001
    ## 
    ## Intercepts:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##    .Dep_EpiDoc4       1.191    0.235    5.059      Inf    0.000    0.729
    ##    .MM_EpiDoc4       -0.275    0.222   -1.239      Inf    0.215   -0.710
    ##    .Dep_EpiDoc2       0.258    0.554    0.465  143.878    0.643   -0.838
    ##    .MM_EpiDoc2        0.079    0.074    1.056      Inf    0.291   -0.067
    ##    .Dep_EpiDoc1      -0.201    0.685   -0.293      Inf    0.770   -1.544
    ##    .MM_EpiDoc1       -0.989    0.182   -5.419      Inf    0.000   -1.346
    ##  ci.upper   Std.lv  Std.all
    ##     1.652    1.191    0.331
    ##     0.160   -0.275   -0.221
    ##     1.354    0.258    0.069
    ##     0.225    0.079    0.115
    ##     1.142   -0.201   -0.054
    ##    -0.631   -0.989   -0.987
    ## 
    ## Variances:
    ##                    Estimate  Std.Err  t-value       df  P(>|t|) ci.lower
    ##    .Dep_EpiDoc4       9.455    1.071    8.825      Inf    0.000    7.353
    ##    .MM_EpiDoc4        1.013    0.087   11.591      Inf    0.000    0.842
    ##    .Dep_EpiDoc2      10.805    1.352    7.995  367.704    0.000    8.147
    ##    .MM_EpiDoc2        0.429    0.075    5.743      Inf    0.000    0.282
    ##    .Dep_EpiDoc1      12.643    1.725    7.329      Inf    0.000    9.262
    ##    .MM_EpiDoc1        0.797    0.125    6.363      Inf    0.000    0.552
    ##  ci.upper   Std.lv  Std.all
    ##    11.557    9.455    0.732
    ##     1.184    1.013    0.652
    ##    13.463   10.805    0.777
    ##     0.575    0.429    0.910
    ##    16.024   12.643    0.903
    ##     1.043    0.797    0.794
    ## 
    ## R-Square:
    ##                    Estimate
    ##     Dep_EpiDoc4       0.268
    ##     MM_EpiDoc4        0.348
    ##     Dep_EpiDoc2       0.223
    ##     MM_EpiDoc2        0.090
    ##     Dep_EpiDoc1       0.097
    ##     MM_EpiDoc1        0.206

``` r
lavaan.mi::lavTestLRT.mi(resED,resSatModel )
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
    ##             Df   Chisq Chisq diff Df diff Pr(>Chisq)     RIV     FMI
    ## resSatModel  2  2.3128                                              
    ## resED       28 93.9125     67.401      26 1.5737e-05 0.17391 0.14814

``` r
lavaan.mi::lavTestLRT.mi(resED, resARAge,resAR_BiD,resAR,resRCov, resSatModel )
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
    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.

    ## Test statistic(s) pooled using the D4 pooling method.
    ##   Pooled statistic: "standard"
    ##   Method to robustify pooled statistic:  "yuan.bentler.mplus"
    ## 
    ##             Df   Chisq Chisq diff Df diff Pr(>Chisq)     RIV      FMI
    ## resSatModel  2  2.3128                                               
    ## resRCov      5  5.5145      2.000       3    0.57247 0.21946 0.179963
    ## resAR       11 70.2953     39.513       6    0.00000 0.08258 0.076281
    ## resAR_BiD   14 70.5426      2.150       3    0.54187 0.44581 0.308346
    ## resARAge    21 84.3807     10.611       7    0.15652 0.18887 0.158865
    ## resED       28 93.9125      7.788       7    0.35162 0.10117 0.091874

## reported on fully adjusted

``` r
partable<-parameterEstimates.mi(resED,standardized =TRUE)
partable<-partable%>%mutate(
  label2=paste0(lhs,op,rhs),
  estimate_p=paste0(round(est,3),' (',round(pvalue,3),')'),
  CI=paste0('[',round(ci.lower,3),':',round(ci.upper,3),']')
)%>%select(-c(lhs,op,rhs))
gPartable<-reshape(partable%>%filter(!block==0), idvar = c('label2'),
                                                   timevar = "group", direction = "wide",sep='_')
gPartable%>%select(label2,estimate_p_1,estimate_p_2,CI_1,CI_2)
```

    ##                      label2   estimate_p_1   estimate_p_2              CI_1
    ## 1   Dep_EpiDoc4~Dep_EpiDoc1  0.184 (0.009)      0.324 (0)     [0.046:0.321]
    ## 2   Dep_EpiDoc4~Dep_EpiDoc2      0.235 (0)      0.235 (0)     [0.127:0.343]
    ## 3     MM_EpiDoc4~MM_EpiDoc2  0.285 (0.002)  0.285 (0.002)        [0.1:0.47]
    ## 4     MM_EpiDoc4~MM_EpiDoc1       0.48 (0)       0.48 (0)     [0.358:0.603]
    ## 5   Dep_EpiDoc2~Dep_EpiDoc1       0.39 (0)       0.39 (0)     [0.284:0.497]
    ## 6     MM_EpiDoc2~MM_EpiDoc1  0.177 (0.001)  0.177 (0.001)     [0.072:0.281]
    ## 7    Dep_EpiDoc4~MM_EpiDoc2  0.506 (0.079)  0.475 (0.111)     [-0.058:1.07]
    ## 8    Dep_EpiDoc4~MM_EpiDoc1         0 (NA)         0 (NA)             [0:0]
    ## 9    Dep_EpiDoc2~MM_EpiDoc1  0.177 (0.311)  0.177 (0.311)    [-0.167:0.522]
    ## 10   MM_EpiDoc4~Dep_EpiDoc2         0 (NA)         0 (NA)             [0:0]
    ## 11   MM_EpiDoc4~Dep_EpiDoc1  0.027 (0.239)    0.02 (0.21)    [-0.018:0.072]
    ## 12   MM_EpiDoc2~Dep_EpiDoc1  0.012 (0.408)  0.026 (0.071)    [-0.017:0.041]
    ## 13  Dep_EpiDoc1~~MM_EpiDoc1         0 (NA)         0 (NA)             [0:0]
    ## 14  Dep_EpiDoc2~~MM_EpiDoc2         0 (NA)         0 (NA)             [0:0]
    ## 15  Dep_EpiDoc4~~MM_EpiDoc4  0.207 (0.229)  0.003 (0.989)    [-0.132:0.546]
    ## 16          Dep_EpiDoc4~Age         0 (NA)         0 (NA)             [0:0]
    ## 17           MM_EpiDoc4~Age      0.022 (0)      0.022 (0)     [0.011:0.033]
    ## 18          Dep_EpiDoc2~Age    0.01 (0.43)    0.01 (0.43)    [-0.016:0.037]
    ## 19           MM_EpiDoc2~Age -0.001 (0.601) -0.001 (0.601)    [-0.004:0.002]
    ## 20          Dep_EpiDoc1~Age  0.045 (0.002)  0.045 (0.002)     [0.017:0.073]
    ## 21           MM_EpiDoc1~Age      0.037 (0)      0.037 (0)     [0.029:0.045]
    ## 22    Dep_EpiDoc4~Education         0 (NA)         0 (NA)             [0:0]
    ## 23    Dep_EpiDoc2~Education      0.473 (0)      0.473 (0)     [0.239:0.707]
    ## 24    Dep_EpiDoc1~Education  0.472 (0.003)      0.808 (0)     [0.162:0.782]
    ## 25     MM_EpiDoc4~Education         0 (NA)         0 (NA)             [0:0]
    ## 26     MM_EpiDoc2~Education         0 (NA)         0 (NA)             [0:0]
    ## 27     MM_EpiDoc1~Education  -0.004 (0.94)   0.06 (0.285)    [-0.104:0.097]
    ## 28 Dep_EpiDoc4~~Dep_EpiDoc4      6.924 (0)      9.455 (0)     [4.974:8.875]
    ## 29   MM_EpiDoc4~~MM_EpiDoc4      0.989 (0)      1.013 (0)       [0.778:1.2]
    ## 30 Dep_EpiDoc2~~Dep_EpiDoc2      6.214 (0)     10.805 (0)     [4.503:7.926]
    ## 31   MM_EpiDoc2~~MM_EpiDoc2      0.193 (0)      0.429 (0)      [0.09:0.295]
    ## 32 Dep_EpiDoc1~~Dep_EpiDoc1       8.43 (0)     12.643 (0)     [6.17:10.689]
    ## 33   MM_EpiDoc1~~MM_EpiDoc1      0.706 (0)      0.797 (0)     [0.538:0.875]
    ## 34                 Age~~Age   204.021 (NA)   135.726 (NA) [204.021:204.021]
    ## 35           Age~~Education     6.794 (NA)      3.68 (NA)     [6.794:6.794]
    ## 36     Education~~Education     1.159 (NA)      1.24 (NA)     [1.159:1.159]
    ## 37            Dep_EpiDoc4~1      1.509 (0)      1.191 (0)     [1.126:1.891]
    ## 38             MM_EpiDoc4~1 -0.438 (0.022) -0.275 (0.215)   [-0.811:-0.064]
    ## 39            Dep_EpiDoc2~1 -0.181 (0.731)  0.258 (0.643)     [-1.22:0.858]
    ## 40             MM_EpiDoc2~1 -0.004 (0.951)  0.079 (0.291)    [-0.122:0.114]
    ## 41            Dep_EpiDoc1~1 -0.687 (0.269)  -0.201 (0.77)    [-1.906:0.532]
    ## 42             MM_EpiDoc1~1     -0.983 (0)     -0.989 (0)   [-1.323:-0.643]
    ## 43                    Age~1    44.188 (NA)    42.481 (NA)   [44.188:44.188]
    ## 44              Education~1     2.823 (NA)     2.409 (NA)     [2.823:2.823]
    ##                 CI_2
    ## 1       [0.19:0.457]
    ## 2      [0.127:0.343]
    ## 3         [0.1:0.47]
    ## 4      [0.358:0.603]
    ## 5      [0.284:0.497]
    ## 6      [0.072:0.281]
    ## 7     [-0.109:1.059]
    ## 8              [0:0]
    ## 9     [-0.167:0.522]
    ## 10             [0:0]
    ## 11    [-0.011:0.052]
    ## 12    [-0.002:0.054]
    ## 13             [0:0]
    ## 14             [0:0]
    ## 15    [-0.369:0.375]
    ## 16             [0:0]
    ## 17     [0.011:0.033]
    ## 18    [-0.016:0.037]
    ## 19    [-0.004:0.002]
    ## 20     [0.017:0.073]
    ## 21     [0.029:0.045]
    ## 22             [0:0]
    ## 23     [0.239:0.707]
    ## 24     [0.412:1.204]
    ## 25             [0:0]
    ## 26             [0:0]
    ## 27     [-0.05:0.169]
    ## 28    [7.353:11.557]
    ## 29     [0.842:1.184]
    ## 30    [8.147:13.463]
    ## 31     [0.282:0.575]
    ## 32    [9.262:16.024]
    ## 33     [0.552:1.043]
    ## 34 [135.726:135.726]
    ## 35       [3.68:3.68]
    ## 36       [1.24:1.24]
    ## 37     [0.729:1.652]
    ## 38      [-0.71:0.16]
    ## 39    [-0.838:1.354]
    ## 40    [-0.067:0.225]
    ## 41    [-1.544:1.142]
    ## 42   [-1.346:-0.631]
    ## 43   [42.481:42.481]
    ## 44     [2.409:2.409]

``` r
source('./code/parTableToCSV.R')
parTableToCSV(resED,path=paste0('./Results/LaggedSEM_mi20_final/',folderName))
```

    ## Joining with `by = join_by(group)`

    ## [1] Group2     label      estimate_p CI         Group     
    ## <0 rows> (or 0-length row.names)

    ## $partable
    ##            lhs op         rhs group label block     est    se      t
    ## 1  Dep_EpiDoc4  ~ Dep_EpiDoc1     1           1   0.184 0.070  2.613
    ## 2  Dep_EpiDoc4  ~ Dep_EpiDoc2     1  dd13     1   0.235 0.055  4.284
    ## 3   MM_EpiDoc4  ~  MM_EpiDoc2     1  mm14     1   0.285 0.094  3.024
    ## 4   MM_EpiDoc4  ~  MM_EpiDoc1     1  mm04     1   0.480 0.062  7.704
    ## 5  Dep_EpiDoc2  ~ Dep_EpiDoc1     1  dd01     1   0.390 0.054  7.219
    ## 6   MM_EpiDoc2  ~  MM_EpiDoc1     1  mm01     1   0.177 0.053  3.319
    ## 7  Dep_EpiDoc4  ~  MM_EpiDoc2     1           1   0.506 0.288  1.759
    ## 8  Dep_EpiDoc4  ~  MM_EpiDoc1     1           1   0.000 0.000     NA
    ## 9  Dep_EpiDoc2  ~  MM_EpiDoc1     1  dm21     1   0.177 0.175  1.014
    ## 10  MM_EpiDoc4  ~ Dep_EpiDoc2     1           1   0.000 0.000     NA
    ## 11  MM_EpiDoc4  ~ Dep_EpiDoc1     1           1   0.027 0.023  1.177
    ## 12  MM_EpiDoc2  ~ Dep_EpiDoc1     1           1   0.012 0.015  0.827
    ## 13 Dep_EpiDoc1 ~~  MM_EpiDoc1     1           1   0.000 0.000     NA
    ## 14 Dep_EpiDoc2 ~~  MM_EpiDoc2     1           1   0.000 0.000     NA
    ## 15 Dep_EpiDoc4 ~~  MM_EpiDoc4     1           1   0.207 0.172  1.207
    ## 16 Dep_EpiDoc4  ~         Age     1           1   0.000 0.000     NA
    ## 17  MM_EpiDoc4  ~         Age     1   ma4     1   0.022 0.006  3.925
    ## 18 Dep_EpiDoc2  ~         Age     1   da2     1   0.010 0.013  0.793
    ## 19  MM_EpiDoc2  ~         Age     1   ma2     1  -0.001 0.002 -0.523
    ## 20 Dep_EpiDoc1  ~         Age     1   da1     1   0.045 0.014  3.104
    ## 21  MM_EpiDoc1  ~         Age     1   ma1     1   0.037 0.004  9.295
    ## 22 Dep_EpiDoc4  ~   Education     1           1   0.000 0.000     NA
    ## 23 Dep_EpiDoc2  ~   Education     1   de2     1   0.473 0.119  3.982
    ## 24 Dep_EpiDoc1  ~   Education     1           1   0.472 0.158  2.983
    ## 25  MM_EpiDoc4  ~   Education     1           1   0.000 0.000     NA
    ## 26  MM_EpiDoc2  ~   Education     1           1   0.000 0.000     NA
    ## 27  MM_EpiDoc1  ~   Education     1           1  -0.004 0.051 -0.076
    ## 28 Dep_EpiDoc4 ~~ Dep_EpiDoc4     1           1   6.924 0.994  6.963
    ## 29  MM_EpiDoc4 ~~  MM_EpiDoc4     1           1   0.989 0.108  9.193
    ## 30 Dep_EpiDoc2 ~~ Dep_EpiDoc2     1           1   6.214 0.863  7.199
    ## 31  MM_EpiDoc2 ~~  MM_EpiDoc2     1           1   0.193 0.052  3.697
    ## 32 Dep_EpiDoc1 ~~ Dep_EpiDoc1     1           1   8.430 1.153  7.311
    ## 33  MM_EpiDoc1 ~~  MM_EpiDoc1     1           1   0.706 0.086  8.214
    ## 34         Age ~~         Age     1           1 204.021 0.000     NA
    ## 35         Age ~~   Education     1           1   6.794 0.000     NA
    ## 36   Education ~~   Education     1           1   1.159 0.000     NA
    ## 37 Dep_EpiDoc4 ~1                 1           1   1.509 0.194  7.758
    ## 38  MM_EpiDoc4 ~1                 1           1  -0.438 0.191 -2.296
    ## 39 Dep_EpiDoc2 ~1                 1           1  -0.181 0.525 -0.344
    ## 40  MM_EpiDoc2 ~1                 1           1  -0.004 0.060 -0.061
    ## 41 Dep_EpiDoc1 ~1                 1           1  -0.687 0.622 -1.105
    ## 42  MM_EpiDoc1 ~1                 1           1  -0.983 0.174 -5.663
    ## 43         Age ~1                 1           1  44.188 0.000     NA
    ## 44   Education ~1                 1           1   2.823 0.000     NA
    ## 45 Dep_EpiDoc4  ~ Dep_EpiDoc1     2           2   0.324 0.068  4.756
    ## 46 Dep_EpiDoc4  ~ Dep_EpiDoc2     2  dd13     2   0.235 0.055  4.284
    ## 47  MM_EpiDoc4  ~  MM_EpiDoc2     2  mm14     2   0.285 0.094  3.024
    ## 48  MM_EpiDoc4  ~  MM_EpiDoc1     2  mm04     2   0.480 0.062  7.704
    ## 49 Dep_EpiDoc2  ~ Dep_EpiDoc1     2  dd01     2   0.390 0.054  7.219
    ## 50  MM_EpiDoc2  ~  MM_EpiDoc1     2  mm01     2   0.177 0.053  3.319
    ## 51 Dep_EpiDoc4  ~  MM_EpiDoc2     2           2   0.475 0.298  1.595
    ## 52 Dep_EpiDoc4  ~  MM_EpiDoc1     2           2   0.000 0.000     NA
    ## 53 Dep_EpiDoc2  ~  MM_EpiDoc1     2  dm21     2   0.177 0.175  1.014
    ## 54  MM_EpiDoc4  ~ Dep_EpiDoc2     2           2   0.000 0.000     NA
    ## 55  MM_EpiDoc4  ~ Dep_EpiDoc1     2           2   0.020 0.016  1.254
    ## 56  MM_EpiDoc2  ~ Dep_EpiDoc1     2           2   0.026 0.014  1.808
    ## 57 Dep_EpiDoc1 ~~  MM_EpiDoc1     2           2   0.000 0.000     NA
    ## 58 Dep_EpiDoc2 ~~  MM_EpiDoc2     2           2   0.000 0.000     NA
    ## 59 Dep_EpiDoc4 ~~  MM_EpiDoc4     2           2   0.003 0.190  0.014
    ## 60 Dep_EpiDoc4  ~         Age     2           2   0.000 0.000     NA
    ## 61  MM_EpiDoc4  ~         Age     2   ma4     2   0.022 0.006  3.925
    ## 62 Dep_EpiDoc2  ~         Age     2   da2     2   0.010 0.013  0.793
    ## 63  MM_EpiDoc2  ~         Age     2   ma2     2  -0.001 0.002 -0.523
    ## 64 Dep_EpiDoc1  ~         Age     2   da1     2   0.045 0.014  3.104
    ## 65  MM_EpiDoc1  ~         Age     2   ma1     2   0.037 0.004  9.295
    ## 66 Dep_EpiDoc4  ~   Education     2           2   0.000 0.000     NA
    ## 67 Dep_EpiDoc2  ~   Education     2   de2     2   0.473 0.119  3.982
    ## 68 Dep_EpiDoc1  ~   Education     2           2   0.808 0.202  3.998
    ## 69  MM_EpiDoc4  ~   Education     2           2   0.000 0.000     NA
    ## 70  MM_EpiDoc2  ~   Education     2           2   0.000 0.000     NA
    ## 71  MM_EpiDoc1  ~   Education     2           2   0.060 0.056  1.070
    ## 72 Dep_EpiDoc4 ~~ Dep_EpiDoc4     2           2   9.455 1.071  8.825
    ## 73  MM_EpiDoc4 ~~  MM_EpiDoc4     2           2   1.013 0.087 11.591
    ## 74 Dep_EpiDoc2 ~~ Dep_EpiDoc2     2           2  10.805 1.352  7.995
    ## 75  MM_EpiDoc2 ~~  MM_EpiDoc2     2           2   0.429 0.075  5.743
    ## 76 Dep_EpiDoc1 ~~ Dep_EpiDoc1     2           2  12.643 1.725  7.329
    ## 77  MM_EpiDoc1 ~~  MM_EpiDoc1     2           2   0.797 0.125  6.363
    ## 78         Age ~~         Age     2           2 135.726 0.000     NA
    ## 79         Age ~~   Education     2           2   3.680 0.000     NA
    ## 80   Education ~~   Education     2           2   1.240 0.000     NA
    ## 81 Dep_EpiDoc4 ~1                 2           2   1.191 0.235  5.059
    ## 82  MM_EpiDoc4 ~1                 2           2  -0.275 0.222 -1.239
    ## 83 Dep_EpiDoc2 ~1                 2           2   0.258 0.554  0.465
    ## 84  MM_EpiDoc2 ~1                 2           2   0.079 0.074  1.056
    ## 85 Dep_EpiDoc1 ~1                 2           2  -0.201 0.685 -0.293
    ## 86  MM_EpiDoc1 ~1                 2           2  -0.989 0.182 -5.419
    ## 87         Age ~1                 2           2  42.481 0.000     NA
    ## 88   Education ~1                 2           2   2.409 0.000     NA
    ##              df pvalue ci.lower ci.upper  std.lv std.all std.nox
    ## 1  7.356640e+02  0.009    0.046    0.321   0.184   0.196   0.196
    ## 2  6.498500e+02  0.000    0.127    0.343   0.235   0.240   0.240
    ## 3  3.226086e+09  0.002    0.100    0.470   0.285   0.107   0.107
    ## 4  1.553764e+08  0.000    0.358    0.603   0.480   0.377   0.377
    ## 5  3.333110e+02  0.000    0.284    0.497   0.390   0.408   0.408
    ## 6  1.545388e+26  0.001    0.072    0.281   0.177   0.371   0.371
    ## 7  3.804153e+03  0.079   -0.058    1.070   0.506   0.084   0.084
    ## 8            NA     NA    0.000    0.000   0.000   0.000   0.000
    ## 9  3.075790e+02  0.311   -0.167    0.522   0.177   0.060   0.060
    ## 10           NA     NA    0.000    0.000   0.000   0.000   0.000
    ## 11 4.717016e+08  0.239   -0.018    0.072   0.027   0.066   0.066
    ## 12 1.786306e+27  0.408   -0.017    0.041   0.012   0.078   0.078
    ## 13           NA     NA    0.000    0.000   0.000   0.000   0.000
    ## 14           NA     NA    0.000    0.000   0.000   0.000   0.000
    ## 15 1.936970e+02  0.229   -0.132    0.546   0.207   0.079   0.079
    ## 16           NA     NA    0.000    0.000   0.000   0.000   0.000
    ## 17 2.208014e+06  0.000    0.011    0.033   0.022   0.252   0.018
    ## 18 8.973700e+01  0.430   -0.016    0.037   0.010   0.051   0.004
    ## 19 1.072407e+23  0.601   -0.004    0.002  -0.001  -0.027  -0.002
    ## 20 8.621924e+19  0.002    0.017    0.073   0.045   0.210   0.015
    ## 21 2.273356e+19  0.000    0.029    0.045   0.037   0.534   0.037
    ## 22           NA     NA    0.000    0.000   0.000   0.000   0.000
    ## 23 3.209210e+02  0.000    0.239    0.707   0.473   0.174   0.162
    ## 24 2.209286e+19  0.003    0.162    0.782   0.472   0.166   0.154
    ## 25           NA     NA    0.000    0.000   0.000   0.000   0.000
    ## 26           NA     NA    0.000    0.000   0.000   0.000   0.000
    ## 27 1.094847e+20  0.940   -0.104    0.097  -0.004  -0.004  -0.004
    ## 28 1.647028e+03  0.000    4.974    8.875   6.924   0.844   0.844
    ## 29 1.631852e+12  0.000    0.778    1.200   0.989   0.619   0.619
    ## 30 1.048270e+02  0.000    4.503    7.926   6.214   0.725   0.725
    ## 31 2.117210e+26  0.000    0.090    0.295   0.193   0.859   0.859
    ## 32 6.351265e+22  0.000    6.170   10.689   8.430   0.898   0.898
    ## 33 5.001300e+20  0.000    0.538    0.875   0.706   0.717   0.717
    ## 34           NA     NA  204.021  204.021 204.021   1.000 204.021
    ## 35           NA     NA    6.794    6.794   6.794   0.442   6.794
    ## 36           NA     NA    1.159    1.159   1.159   1.000   1.159
    ## 37 3.994430e+02  0.000    1.126    1.891   1.509   0.527   0.527
    ## 38 1.082646e+06  0.022   -0.811   -0.064  -0.438  -0.346  -0.346
    ## 39 1.295040e+02  0.731   -1.220    0.858  -0.181  -0.062  -0.062
    ## 40 4.323210e+22  0.951   -0.122    0.114  -0.004  -0.008  -0.008
    ## 41 2.256199e+23  0.269   -1.906    0.532  -0.687  -0.224  -0.224
    ## 42 2.646125e+18  0.000   -1.323   -0.643  -0.983  -0.991  -0.991
    ## 43           NA     NA   44.188   44.188  44.188   3.094  44.188
    ## 44           NA     NA    2.823    2.823   2.823   2.622   2.823
    ## 45 3.329390e+03  0.000    0.190    0.457   0.324   0.337   0.337
    ## 46 6.498500e+02  0.000    0.127    0.343   0.235   0.244   0.244
    ## 47 3.226086e+09  0.002    0.100    0.470   0.285   0.157   0.157
    ## 48 1.553764e+08  0.000    0.358    0.603   0.480   0.386   0.386
    ## 49 3.333110e+02  0.000    0.284    0.497   0.390   0.391   0.391
    ## 50 1.545388e+26  0.001    0.072    0.281   0.177   0.258   0.258
    ## 51 3.209892e+03  0.111   -0.109    1.059   0.475   0.091   0.091
    ## 52           NA     NA    0.000    0.000   0.000   0.000   0.000
    ## 53 3.075790e+02  0.311   -0.167    0.522   0.177   0.048   0.048
    ## 54           NA     NA    0.000    0.000   0.000   0.000   0.000
    ## 55 1.711400e+09  0.210   -0.011    0.052   0.020   0.061   0.061
    ## 56 9.454530e+24  0.071   -0.002    0.054   0.026   0.141   0.141
    ## 57           NA     NA    0.000    0.000   0.000   0.000   0.000
    ## 58           NA     NA    0.000    0.000   0.000   0.000   0.000
    ## 59 3.844291e+03  0.989   -0.369    0.375   0.003   0.001   0.001
    ## 60           NA     NA    0.000    0.000   0.000   0.000   0.000
    ## 61 2.208014e+06  0.000    0.011    0.033   0.022   0.208   0.018
    ## 62 8.973700e+01  0.430   -0.016    0.037   0.010   0.033   0.003
    ## 63 1.072407e+23  0.601   -0.004    0.002  -0.001  -0.015  -0.001
    ## 64 8.621924e+19  0.002    0.017    0.073   0.045   0.140   0.012
    ## 65 2.273356e+19  0.000    0.029    0.045   0.037   0.431   0.037
    ## 66           NA     NA    0.000    0.000   0.000   0.000   0.000
    ## 67 3.209210e+02  0.000    0.239    0.707   0.473   0.141   0.127
    ## 68 2.482674e+19  0.000    0.412    1.204   0.808   0.240   0.216
    ## 69           NA     NA    0.000    0.000   0.000   0.000   0.000
    ## 70           NA     NA    0.000    0.000   0.000   0.000   0.000
    ## 71 1.709213e+22  0.285   -0.050    0.169   0.060   0.066   0.059
    ## 72 1.495549e+03  0.000    7.353   11.557   9.455   0.732   0.732
    ## 73 7.717348e+10  0.000    0.842    1.184   1.013   0.652   0.652
    ## 74 3.677040e+02  0.000    8.147   13.463  10.805   0.777   0.777
    ## 75 5.189000e+27  0.000    0.282    0.575   0.429   0.910   0.910
    ## 76 1.165974e+22  0.000    9.262   16.024  12.643   0.903   0.903
    ## 77 1.731922e+21  0.000    0.552    1.043   0.797   0.794   0.794
    ## 78           NA     NA  135.726  135.726 135.726   1.000 135.726
    ## 79           NA     NA    3.680    3.680   3.680   0.284   3.680
    ## 80           NA     NA    1.240    1.240   1.240   1.000   1.240
    ## 81 7.049060e+03  0.000    0.729    1.652   1.191   0.331   0.331
    ## 82 2.233103e+06  0.215   -0.710    0.160  -0.275  -0.221  -0.221
    ## 83 1.438780e+02  0.643   -0.838    1.354   0.258   0.069   0.069
    ## 84 3.856466e+22  0.291   -0.067    0.225   0.079   0.115   0.115
    ## 85 1.066903e+22  0.770   -1.544    1.142  -0.201  -0.054  -0.054
    ## 86 1.477186e+19  0.000   -1.346   -0.631  -0.989  -0.987  -0.987
    ## 87           NA     NA   42.481   42.481  42.481   3.646  42.481
    ## 88           NA     NA    2.409    2.409   2.409   2.163   2.409
    ##                      label2     estimate_p                CI  Group
    ## 1   Dep_EpiDoc4~Dep_EpiDoc1  0.184 (0.009)     [0.046:0.321]   Male
    ## 2   Dep_EpiDoc4~Dep_EpiDoc2      0.235 (0)     [0.127:0.343]   Male
    ## 3     MM_EpiDoc4~MM_EpiDoc2  0.285 (0.002)        [0.1:0.47]   Male
    ## 4     MM_EpiDoc4~MM_EpiDoc1       0.48 (0)     [0.358:0.603]   Male
    ## 5   Dep_EpiDoc2~Dep_EpiDoc1       0.39 (0)     [0.284:0.497]   Male
    ## 6     MM_EpiDoc2~MM_EpiDoc1  0.177 (0.001)     [0.072:0.281]   Male
    ## 7    Dep_EpiDoc4~MM_EpiDoc2  0.506 (0.079)     [-0.058:1.07]   Male
    ## 8    Dep_EpiDoc4~MM_EpiDoc1         0 (NA)             [0:0]   Male
    ## 9    Dep_EpiDoc2~MM_EpiDoc1  0.177 (0.311)    [-0.167:0.522]   Male
    ## 10   MM_EpiDoc4~Dep_EpiDoc2         0 (NA)             [0:0]   Male
    ## 11   MM_EpiDoc4~Dep_EpiDoc1  0.027 (0.239)    [-0.018:0.072]   Male
    ## 12   MM_EpiDoc2~Dep_EpiDoc1  0.012 (0.408)    [-0.017:0.041]   Male
    ## 13  Dep_EpiDoc1~~MM_EpiDoc1         0 (NA)             [0:0]   Male
    ## 14  Dep_EpiDoc2~~MM_EpiDoc2         0 (NA)             [0:0]   Male
    ## 15  Dep_EpiDoc4~~MM_EpiDoc4  0.207 (0.229)    [-0.132:0.546]   Male
    ## 16          Dep_EpiDoc4~Age         0 (NA)             [0:0]   Male
    ## 17           MM_EpiDoc4~Age      0.022 (0)     [0.011:0.033]   Male
    ## 18          Dep_EpiDoc2~Age    0.01 (0.43)    [-0.016:0.037]   Male
    ## 19           MM_EpiDoc2~Age -0.001 (0.601)    [-0.004:0.002]   Male
    ## 20          Dep_EpiDoc1~Age  0.045 (0.002)     [0.017:0.073]   Male
    ## 21           MM_EpiDoc1~Age      0.037 (0)     [0.029:0.045]   Male
    ## 22    Dep_EpiDoc4~Education         0 (NA)             [0:0]   Male
    ## 23    Dep_EpiDoc2~Education      0.473 (0)     [0.239:0.707]   Male
    ## 24    Dep_EpiDoc1~Education  0.472 (0.003)     [0.162:0.782]   Male
    ## 25     MM_EpiDoc4~Education         0 (NA)             [0:0]   Male
    ## 26     MM_EpiDoc2~Education         0 (NA)             [0:0]   Male
    ## 27     MM_EpiDoc1~Education  -0.004 (0.94)    [-0.104:0.097]   Male
    ## 28 Dep_EpiDoc4~~Dep_EpiDoc4      6.924 (0)     [4.974:8.875]   Male
    ## 29   MM_EpiDoc4~~MM_EpiDoc4      0.989 (0)       [0.778:1.2]   Male
    ## 30 Dep_EpiDoc2~~Dep_EpiDoc2      6.214 (0)     [4.503:7.926]   Male
    ## 31   MM_EpiDoc2~~MM_EpiDoc2      0.193 (0)      [0.09:0.295]   Male
    ## 32 Dep_EpiDoc1~~Dep_EpiDoc1       8.43 (0)     [6.17:10.689]   Male
    ## 33   MM_EpiDoc1~~MM_EpiDoc1      0.706 (0)     [0.538:0.875]   Male
    ## 34                 Age~~Age   204.021 (NA) [204.021:204.021]   Male
    ## 35           Age~~Education     6.794 (NA)     [6.794:6.794]   Male
    ## 36     Education~~Education     1.159 (NA)     [1.159:1.159]   Male
    ## 37            Dep_EpiDoc4~1      1.509 (0)     [1.126:1.891]   Male
    ## 38             MM_EpiDoc4~1 -0.438 (0.022)   [-0.811:-0.064]   Male
    ## 39            Dep_EpiDoc2~1 -0.181 (0.731)     [-1.22:0.858]   Male
    ## 40             MM_EpiDoc2~1 -0.004 (0.951)    [-0.122:0.114]   Male
    ## 41            Dep_EpiDoc1~1 -0.687 (0.269)    [-1.906:0.532]   Male
    ## 42             MM_EpiDoc1~1     -0.983 (0)   [-1.323:-0.643]   Male
    ## 43                    Age~1    44.188 (NA)   [44.188:44.188]   Male
    ## 44              Education~1     2.823 (NA)     [2.823:2.823]   Male
    ## 45  Dep_EpiDoc4~Dep_EpiDoc1      0.324 (0)      [0.19:0.457] Female
    ## 46  Dep_EpiDoc4~Dep_EpiDoc2      0.235 (0)     [0.127:0.343] Female
    ## 47    MM_EpiDoc4~MM_EpiDoc2  0.285 (0.002)        [0.1:0.47] Female
    ## 48    MM_EpiDoc4~MM_EpiDoc1       0.48 (0)     [0.358:0.603] Female
    ## 49  Dep_EpiDoc2~Dep_EpiDoc1       0.39 (0)     [0.284:0.497] Female
    ## 50    MM_EpiDoc2~MM_EpiDoc1  0.177 (0.001)     [0.072:0.281] Female
    ## 51   Dep_EpiDoc4~MM_EpiDoc2  0.475 (0.111)    [-0.109:1.059] Female
    ## 52   Dep_EpiDoc4~MM_EpiDoc1         0 (NA)             [0:0] Female
    ## 53   Dep_EpiDoc2~MM_EpiDoc1  0.177 (0.311)    [-0.167:0.522] Female
    ## 54   MM_EpiDoc4~Dep_EpiDoc2         0 (NA)             [0:0] Female
    ## 55   MM_EpiDoc4~Dep_EpiDoc1    0.02 (0.21)    [-0.011:0.052] Female
    ## 56   MM_EpiDoc2~Dep_EpiDoc1  0.026 (0.071)    [-0.002:0.054] Female
    ## 57  Dep_EpiDoc1~~MM_EpiDoc1         0 (NA)             [0:0] Female
    ## 58  Dep_EpiDoc2~~MM_EpiDoc2         0 (NA)             [0:0] Female
    ## 59  Dep_EpiDoc4~~MM_EpiDoc4  0.003 (0.989)    [-0.369:0.375] Female
    ## 60          Dep_EpiDoc4~Age         0 (NA)             [0:0] Female
    ## 61           MM_EpiDoc4~Age      0.022 (0)     [0.011:0.033] Female
    ## 62          Dep_EpiDoc2~Age    0.01 (0.43)    [-0.016:0.037] Female
    ## 63           MM_EpiDoc2~Age -0.001 (0.601)    [-0.004:0.002] Female
    ## 64          Dep_EpiDoc1~Age  0.045 (0.002)     [0.017:0.073] Female
    ## 65           MM_EpiDoc1~Age      0.037 (0)     [0.029:0.045] Female
    ## 66    Dep_EpiDoc4~Education         0 (NA)             [0:0] Female
    ## 67    Dep_EpiDoc2~Education      0.473 (0)     [0.239:0.707] Female
    ## 68    Dep_EpiDoc1~Education      0.808 (0)     [0.412:1.204] Female
    ## 69     MM_EpiDoc4~Education         0 (NA)             [0:0] Female
    ## 70     MM_EpiDoc2~Education         0 (NA)             [0:0] Female
    ## 71     MM_EpiDoc1~Education   0.06 (0.285)     [-0.05:0.169] Female
    ## 72 Dep_EpiDoc4~~Dep_EpiDoc4      9.455 (0)    [7.353:11.557] Female
    ## 73   MM_EpiDoc4~~MM_EpiDoc4      1.013 (0)     [0.842:1.184] Female
    ## 74 Dep_EpiDoc2~~Dep_EpiDoc2     10.805 (0)    [8.147:13.463] Female
    ## 75   MM_EpiDoc2~~MM_EpiDoc2      0.429 (0)     [0.282:0.575] Female
    ## 76 Dep_EpiDoc1~~Dep_EpiDoc1     12.643 (0)    [9.262:16.024] Female
    ## 77   MM_EpiDoc1~~MM_EpiDoc1      0.797 (0)     [0.552:1.043] Female
    ## 78                 Age~~Age   135.726 (NA) [135.726:135.726] Female
    ## 79           Age~~Education      3.68 (NA)       [3.68:3.68] Female
    ## 80     Education~~Education      1.24 (NA)       [1.24:1.24] Female
    ## 81            Dep_EpiDoc4~1      1.191 (0)     [0.729:1.652] Female
    ## 82             MM_EpiDoc4~1 -0.275 (0.215)      [-0.71:0.16] Female
    ## 83            Dep_EpiDoc2~1  0.258 (0.643)    [-0.838:1.354] Female
    ## 84             MM_EpiDoc2~1  0.079 (0.291)    [-0.067:0.225] Female
    ## 85            Dep_EpiDoc1~1  -0.201 (0.77)    [-1.544:1.142] Female
    ## 86             MM_EpiDoc1~1     -0.989 (0)   [-1.346:-0.631] Female
    ## 87                    Age~1    42.481 (NA)   [42.481:42.481] Female
    ## 88              Education~1     2.409 (NA)     [2.409:2.409] Female
    ## 
    ## $gPartable
    ##                      label2 estimate_p_Male estimate_p_Female           CI_Male
    ## 1   Dep_EpiDoc4~Dep_EpiDoc1   0.184 (0.009)         0.324 (0)     [0.046:0.321]
    ## 2   Dep_EpiDoc4~Dep_EpiDoc2       0.235 (0)         0.235 (0)     [0.127:0.343]
    ## 3     MM_EpiDoc4~MM_EpiDoc2   0.285 (0.002)     0.285 (0.002)        [0.1:0.47]
    ## 4     MM_EpiDoc4~MM_EpiDoc1        0.48 (0)          0.48 (0)     [0.358:0.603]
    ## 5   Dep_EpiDoc2~Dep_EpiDoc1        0.39 (0)          0.39 (0)     [0.284:0.497]
    ## 6     MM_EpiDoc2~MM_EpiDoc1   0.177 (0.001)     0.177 (0.001)     [0.072:0.281]
    ## 7    Dep_EpiDoc4~MM_EpiDoc2   0.506 (0.079)     0.475 (0.111)     [-0.058:1.07]
    ## 8    Dep_EpiDoc4~MM_EpiDoc1          0 (NA)            0 (NA)             [0:0]
    ## 9    Dep_EpiDoc2~MM_EpiDoc1   0.177 (0.311)     0.177 (0.311)    [-0.167:0.522]
    ## 10   MM_EpiDoc4~Dep_EpiDoc2          0 (NA)            0 (NA)             [0:0]
    ## 11   MM_EpiDoc4~Dep_EpiDoc1   0.027 (0.239)       0.02 (0.21)    [-0.018:0.072]
    ## 12   MM_EpiDoc2~Dep_EpiDoc1   0.012 (0.408)     0.026 (0.071)    [-0.017:0.041]
    ## 13  Dep_EpiDoc1~~MM_EpiDoc1          0 (NA)            0 (NA)             [0:0]
    ## 14  Dep_EpiDoc2~~MM_EpiDoc2          0 (NA)            0 (NA)             [0:0]
    ## 15  Dep_EpiDoc4~~MM_EpiDoc4   0.207 (0.229)     0.003 (0.989)    [-0.132:0.546]
    ## 16          Dep_EpiDoc4~Age          0 (NA)            0 (NA)             [0:0]
    ## 17           MM_EpiDoc4~Age       0.022 (0)         0.022 (0)     [0.011:0.033]
    ## 18          Dep_EpiDoc2~Age     0.01 (0.43)       0.01 (0.43)    [-0.016:0.037]
    ## 19           MM_EpiDoc2~Age  -0.001 (0.601)    -0.001 (0.601)    [-0.004:0.002]
    ## 20          Dep_EpiDoc1~Age   0.045 (0.002)     0.045 (0.002)     [0.017:0.073]
    ## 21           MM_EpiDoc1~Age       0.037 (0)         0.037 (0)     [0.029:0.045]
    ## 22    Dep_EpiDoc4~Education          0 (NA)            0 (NA)             [0:0]
    ## 23    Dep_EpiDoc2~Education       0.473 (0)         0.473 (0)     [0.239:0.707]
    ## 24    Dep_EpiDoc1~Education   0.472 (0.003)         0.808 (0)     [0.162:0.782]
    ## 25     MM_EpiDoc4~Education          0 (NA)            0 (NA)             [0:0]
    ## 26     MM_EpiDoc2~Education          0 (NA)            0 (NA)             [0:0]
    ## 27     MM_EpiDoc1~Education   -0.004 (0.94)      0.06 (0.285)    [-0.104:0.097]
    ## 28 Dep_EpiDoc4~~Dep_EpiDoc4       6.924 (0)         9.455 (0)     [4.974:8.875]
    ## 29   MM_EpiDoc4~~MM_EpiDoc4       0.989 (0)         1.013 (0)       [0.778:1.2]
    ## 30 Dep_EpiDoc2~~Dep_EpiDoc2       6.214 (0)        10.805 (0)     [4.503:7.926]
    ## 31   MM_EpiDoc2~~MM_EpiDoc2       0.193 (0)         0.429 (0)      [0.09:0.295]
    ## 32 Dep_EpiDoc1~~Dep_EpiDoc1        8.43 (0)        12.643 (0)     [6.17:10.689]
    ## 33   MM_EpiDoc1~~MM_EpiDoc1       0.706 (0)         0.797 (0)     [0.538:0.875]
    ## 34                 Age~~Age    204.021 (NA)      135.726 (NA) [204.021:204.021]
    ## 35           Age~~Education      6.794 (NA)         3.68 (NA)     [6.794:6.794]
    ## 36     Education~~Education      1.159 (NA)         1.24 (NA)     [1.159:1.159]
    ## 37            Dep_EpiDoc4~1       1.509 (0)         1.191 (0)     [1.126:1.891]
    ## 38             MM_EpiDoc4~1  -0.438 (0.022)    -0.275 (0.215)   [-0.811:-0.064]
    ## 39            Dep_EpiDoc2~1  -0.181 (0.731)     0.258 (0.643)     [-1.22:0.858]
    ## 40             MM_EpiDoc2~1  -0.004 (0.951)     0.079 (0.291)    [-0.122:0.114]
    ## 41            Dep_EpiDoc1~1  -0.687 (0.269)     -0.201 (0.77)    [-1.906:0.532]
    ## 42             MM_EpiDoc1~1      -0.983 (0)        -0.989 (0)   [-1.323:-0.643]
    ## 43                    Age~1     44.188 (NA)       42.481 (NA)   [44.188:44.188]
    ## 44              Education~1      2.823 (NA)        2.409 (NA)     [2.823:2.823]
    ##            CI_Female
    ## 1       [0.19:0.457]
    ## 2      [0.127:0.343]
    ## 3         [0.1:0.47]
    ## 4      [0.358:0.603]
    ## 5      [0.284:0.497]
    ## 6      [0.072:0.281]
    ## 7     [-0.109:1.059]
    ## 8              [0:0]
    ## 9     [-0.167:0.522]
    ## 10             [0:0]
    ## 11    [-0.011:0.052]
    ## 12    [-0.002:0.054]
    ## 13             [0:0]
    ## 14             [0:0]
    ## 15    [-0.369:0.375]
    ## 16             [0:0]
    ## 17     [0.011:0.033]
    ## 18    [-0.016:0.037]
    ## 19    [-0.004:0.002]
    ## 20     [0.017:0.073]
    ## 21     [0.029:0.045]
    ## 22             [0:0]
    ## 23     [0.239:0.707]
    ## 24     [0.412:1.204]
    ## 25             [0:0]
    ## 26             [0:0]
    ## 27     [-0.05:0.169]
    ## 28    [7.353:11.557]
    ## 29     [0.842:1.184]
    ## 30    [8.147:13.463]
    ## 31     [0.282:0.575]
    ## 32    [9.262:16.024]
    ## 33     [0.552:1.043]
    ## 34 [135.726:135.726]
    ## 35       [3.68:3.68]
    ## 36       [1.24:1.24]
    ## 37     [0.729:1.652]
    ## 38      [-0.71:0.16]
    ## 39    [-0.838:1.354]
    ## 40    [-0.067:0.225]
    ## 41    [-1.544:1.142]
    ## 42   [-1.346:-0.631]
    ## 43   [42.481:42.481]
    ## 44     [2.409:2.409]

# GENERATING PLOTs

## prepare for miPLOt

``` r
library(mice)
```

    ## 
    ## Attaching package: 'mice'

    ## The following object is masked from 'package:stats':
    ## 
    ##     filter

    ## The following objects are masked from 'package:base':
    ## 
    ##     cbind, rbind

``` r
data_bal=complete(data.mi,0)
res_dt0<-sem(data=data_bal,
                   model=Model_ED,
                    estimator = "ML", test = "yuan.bentler.mplus",
         likelihood = "wishart", 
           sampling.weights='ipw',
         orthogonal=TRUE
  )
```

## general plot

``` r
source('./code/plotLavaanMI.R')
lay=get_layout(
  "Age", "","", '',"","",'',"",  "","",
            "", "","", '',"","",'',"",  "","",
   "", "","", '',"","",'Dep_EpiDoc2',"",  "","Dep_EpiDoc4",
          "", '',"", "","", "", "", "","", "", 
            "",'',"","Dep_EpiDoc1", "", "","",  "","", "",
           "", "","", '',"","",'',"",  "","",
            "Education", "","", '',"","",'',"",  "","",
  "",'',"","", "","",'MM_EpiDoc2',"","",  "MM_EpiDoc4",
            "", "","", '',"","",'',"",  "","",
           "",'',"", 'MM_EpiDoc1',"","","","", "", '',
           rows = 10)
```

### automatic

## groups

``` r
dfPar<-get_edges(res_dt0)
p <- prepare_graph(res_dt0, layout = lay)
dfNPar<-p$nodes
partable<-parameterEstimates.mi(resED,standardized =TRUE)
dfPar2<-left_join(dfPar%>%select(rhs,lhs,op,arrow,est_sig_std),
                    partable)%>%mutate(
                      pvalue_temp=ifelse(pvalue<0.001,'<0.001',sprintf('=%.3f',pvalue)),
                      label2=sprintf('%.3f (p%s)',std.all,pvalue_temp),
                      label2=ifelse(op=='~~','',label2)
                    )
```

    ## Joining with `by = join_by(rhs, lhs, op)`

``` r
dfPar2<-dfPar2%>%mutate(
  path= paste0(rhs,op,lhs)
)
dfPar2<-dfPar2%>%group_by(path)%>%mutate(label=paste(paste(group,' ',est_sig_std),collapse=';\n'))
```

``` r
b<-prepare_graph(res_dt0, layout = lay)%>%
  edit_graph( {label = dfPar2$label},element='edges') %>%
  edit_graph({label_location = 0.3})%>%
  edit_graph({label_size = 2},element='edges')%>%
  color_var("white") %>%
  plot(size = c(15, 2))
```

    ## Warning in `[<-.data.frame`(`*tmp*`, nl, value = list(from = c("Dep_EpiDoc1", :
    ## replacement element 4 has 72 rows to replace 36 rows

``` r
partable<-parameterEstimates.mi(resED,standardized =TRUE)
partable<-partable%>%mutate(
  label= paste0(lhs,op,rhs)
)%>%filter(op=='~')
partable<-partable%>%mutate(
  index=as.numeric(rownames(partable)))

fitms<-fitMeasures(resED)
```

    ## Warning: lavaan->lav_options_set():  
    ##    observed.information for ALL test statistics is set to h1.

``` r
fits<-paste0('CFI: ',round(fitms[['cfi']],3),' TLI: ',round(fitms[['tli']],3),' RMSEA: ',round(fitms[['rmsea']],3),' SRMR: ',round(fitms[['srmr']],3))
print(fits)
```

    ## [1] "CFI: 0.934 TLI: 0.872 RMSEA: 0.084 SRMR: 0.066"

``` r
partable$group<-factor(partable$group,labels=c('male','female'))
gall<-ggplot(data=partable,aes(y = factor(group), x = est,color=factor(group)))+
  geom_point(shape = 18) +  
  geom_errorbarh(aes(xmin = ci.lower, xmax = ci.upper), height = 0.25) +
  geom_vline(xintercept =0, color = "red",  cex = 0.5, alpha = 0.5)+
    facet_wrap( ~ label,ncol=2)+
  scale_y_discrete(name = "",labels = unique(data_bal$depression)) +theme_bw() +xlab('estimate')+
  theme(strip.text.y = element_text(angle = 0))
```

    ## Warning: `geom_errobarh()` was deprecated in ggplot2 4.0.0.
    ## ℹ Please use the `orientation` argument of `geom_errorbar()` instead.
    ## This warning is displayed once per session.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

    ## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `linewidth` instead.
    ## This warning is displayed once per session.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

``` r
gall
```

    ## `height` was translated to `width`.

![](4Imp2_SEM_crossLagged_MGRoup_sex_git_files/figure-gfm/unnamed-chunk-34-1.png)<!-- -->

``` r
library(ggpubr)
ggarrange(b, gall, widths = c(5,3),          ncol = 2)
```

    ## `height` was translated to `width`.

![](4Imp2_SEM_crossLagged_MGRoup_sex_git_files/figure-gfm/unnamed-chunk-34-2.png)<!-- -->

``` r
ggsave(filename=paste0('./Results/LaggedSEM_mi20_final/',folderName,'/laggedSEM_group1.png'), width = 35, height = 15, units = "cm")
```

## general plot + forest plot for groups

## just ed and age

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
    ##  ! package         * version  date (UTC) lib source
    ##  P abind             1.4-8    2024-09-12 [?] RSPM
    ##  P backports         1.5.0    2024-05-23 [?] RSPM
    ##  P boot              1.3-31   2024-08-28 [?] RSPM
    ##  P broom             1.0.9    2025-07-28 [?] RSPM
    ##  P cachem            1.1.0    2024-05-16 [?] CRAN (R 4.5.2)
    ##  P car               3.1-3    2024-09-27 [?] RSPM
    ##  P carData           3.0-5    2022-01-06 [?] RSPM
    ##  P checkmate         2.3.3    2025-08-18 [?] RSPM
    ##  P cli               3.6.5    2025-04-23 [?] CRAN (R 4.5.2)
    ##  P coda              0.19-4.1 2024-01-31 [?] RSPM
    ##  P codetools         0.2-20   2024-03-31 [?] RSPM
    ##  P CompQuadForm      1.4.4    2025-07-13 [?] RSPM
    ##  P cowplot           1.2.0    2025-07-07 [?] RSPM
    ##  P data.table        1.17.8   2025-07-10 [?] RSPM
    ##  P dbscan            1.2.3    2025-08-20 [?] RSPM
    ##  P devtools          2.5.0    2026-03-14 [?] CRAN (R 4.5.2)
    ##  P digest            0.6.37   2024-08-19 [?] RSPM
    ##  P dplyr           * 1.1.4    2023-11-17 [?] CRAN (R 4.5.2)
    ##  P ellipsis          0.3.2    2021-04-29 [?] RSPM
    ##  P evaluate          1.0.4    2025-06-18 [?] RSPM
    ##  P farver            2.1.2    2024-05-13 [?] CRAN (R 4.5.2)
    ##  P fastDummies       1.7.5    2025-01-20 [?] RSPM
    ##  P fastmap           1.2.0    2024-05-15 [?] CRAN (R 4.5.2)
    ##  P foreach           1.5.2    2022-02-02 [?] RSPM
    ##  P Formula           1.2-5    2023-02-24 [?] RSPM
    ##  P fs                2.0.1    2026-03-24 [?] RSPM
    ##  P future            1.70.0   2026-03-14 [?] RSPM
    ##  P future.apply      1.20.0   2025-06-06 [?] RSPM
    ##  P generics          0.1.4    2025-05-09 [?] CRAN (R 4.5.2)
    ##  P ggplot2         * 4.0.0    2025-09-11 [?] RSPM
    ##  P ggpubr          * 0.6.2    2025-10-17 [?] CRAN (R 4.5.2)
    ##  P ggsignif          0.6.4    2022-10-13 [?] RSPM
    ##  P glmnet            4.1-10   2025-07-17 [?] RSPM
    ##  P globals           0.19.1   2026-03-13 [?] RSPM
    ##  P glue              1.8.0    2024-09-30 [?] CRAN (R 4.5.2)
    ##  P gridExtra         2.3      2017-09-09 [?] CRAN (R 4.5.2)
    ##  P gsubfn            0.7      2018-03-16 [?] RSPM
    ##  P gtable            0.3.6    2024-10-25 [?] CRAN (R 4.5.2)
    ##  P htmltools         0.5.8.1  2024-04-04 [?] RSPM
    ##  P httr              1.4.7    2023-08-15 [?] RSPM
    ##  P igraph            2.1.4    2025-01-23 [?] RSPM
    ##  P iterators         1.0.14   2022-02-05 [?] RSPM
    ##  P jomo              2.7-6    2023-04-15 [?] RSPM
    ##  P knitr             1.50     2025-03-16 [?] RSPM
    ##  P labeling          0.4.3    2023-08-29 [?] CRAN (R 4.5.2)
    ##  P lattice           0.22-6   2024-03-20 [?] RSPM
    ##  P lavaan          * 0.6-19   2024-09-26 [?] CRAN (R 4.5.2)
    ##  P lavaan.mi       * 0.1-0    2025-03-10 [?] RSPM
    ##  P lifecycle         1.0.5    2026-01-08 [?] CRAN (R 4.5.2)
    ##  P listenv           0.9.1    2024-01-29 [?] RSPM
    ##  P lme4              1.1-37   2025-03-26 [?] RSPM
    ##  P magrittr          2.0.3    2022-03-30 [?] RSPM
    ##  P MASS              7.3-61   2024-06-13 [?] RSPM
    ##  P Matrix            1.7-1    2024-10-18 [?] RSPM
    ##  P memoise           2.0.1    2021-11-26 [?] CRAN (R 4.5.2)
    ##  P mice            * 3.18.0   2025-05-27 [?] RSPM
    ##  P minqa             1.2.8    2024-08-17 [?] RSPM
    ##  P mitml             0.4-5    2023-03-08 [?] RSPM
    ##  P mnormt            2.1.1    2022-09-26 [?] RSPM
    ##  P MplusAutomation   1.2      2025-09-02 [?] RSPM
    ##  P mvtnorm           1.3-3    2025-01-10 [?] RSPM
    ##  P nlme              3.1-166  2024-08-14 [?] RSPM
    ##  P nloptr            2.2.1    2025-03-17 [?] RSPM
    ##  P nnet              7.3-19   2023-05-03 [?] RSPM
    ##  P nonnest2          0.5-8    2024-08-28 [?] RSPM
    ##  P pan               1.9      2023-12-07 [?] RSPM
    ##  P pander            0.6.6    2025-03-01 [?] RSPM
    ##  P parallelly        1.45.1   2025-07-24 [?] RSPM
    ##  P pbivnorm          0.6.0    2015-01-23 [?] RSPM
    ##  P pillar            1.11.0   2025-07-04 [?] RSPM
    ##  P pkgbuild          1.4.8    2025-05-26 [?] CRAN (R 4.5.2)
    ##  P pkgconfig         2.0.3    2019-09-22 [?] CRAN (R 4.5.2)
    ##  P pkgload           1.5.1    2026-04-01 [?] CRAN (R 4.5.2)
    ##  P plyr              1.8.9    2023-10-02 [?] RSPM
    ##  P progressr         0.15.1   2024-11-22 [?] RSPM
    ##  P proto             1.0.0    2016-10-29 [?] RSPM
    ##  P psych             2.5.6    2025-06-23 [?] RSPM
    ##  P purrr             1.2.1    2026-01-09 [?] CRAN (R 4.5.2)
    ##  P quadprog          1.5-8    2019-11-20 [?] RSPM
    ##  P R6                2.6.1    2025-02-15 [?] CRAN (R 4.5.2)
    ##  P ragg              1.5.2    2026-03-23 [?] CRAN (R 4.5.2)
    ##  P RANN              2.6.2    2024-08-25 [?] RSPM
    ##  P rbibutils         2.4.1    2026-01-21 [?] RSPM
    ##  P RColorBrewer      1.1-3    2022-04-03 [?] CRAN (R 4.5.2)
    ##  P Rcpp              1.1.0    2025-07-02 [?] RSPM
    ##  P Rdpack            2.6.4    2025-04-09 [?] RSPM
    ##  P reformulas        0.4.4    2026-02-02 [?] RSPM
    ##    renv              1.1.5    2025-07-24 [1] RSPM
    ##  P rlang             1.2.0    2026-04-06 [?] RSPM
    ##  P rmarkdown         2.29     2024-11-04 [?] RSPM
    ##  P rpart             4.1.23   2023-12-05 [?] RSPM
    ##  P rstatix           0.7.3    2025-10-18 [?] RSPM
    ##  P rstudioapi        0.17.1   2024-10-22 [?] RSPM
    ##  P S7                0.2.0    2024-11-07 [?] RSPM
    ##  P sandwich          3.1-1    2024-09-15 [?] CRAN (R 4.5.2)
    ##  P scales            1.4.0    2025-04-24 [?] CRAN (R 4.5.2)
    ##  P sessioninfo       1.2.3    2025-02-05 [?] CRAN (R 4.5.2)
    ##  P shape             1.4.6.1  2024-02-23 [?] RSPM
    ##  P survival          3.7-0    2024-06-05 [?] RSPM
    ##  P systemfonts       1.3.1    2025-10-01 [?] RSPM
    ##  P texreg            1.39.4   2024-07-24 [?] RSPM
    ##  P textshaping       1.0.4    2025-10-10 [?] RSPM
    ##  P tibble            3.3.0    2025-06-08 [?] RSPM
    ##  P tidyr             1.3.1    2024-01-24 [?] RSPM
    ##  P tidyselect        1.2.1    2024-03-11 [?] CRAN (R 4.5.2)
    ##  P tidySEM         * 0.2.9    2025-07-30 [?] RSPM
    ##  P usethis           3.2.1    2025-09-06 [?] CRAN (R 4.5.2)
    ##  P vctrs             0.7.2    2026-03-21 [?] CRAN (R 4.5.2)
    ##  P withr             3.0.2    2024-10-28 [?] CRAN (R 4.5.2)
    ##  P xfun              0.53     2025-08-19 [?] CRAN (R 4.5.2)
    ##  P xtable            1.8-4    2019-04-21 [?] RSPM
    ##  P yaml              2.3.10   2024-07-26 [?] RSPM
    ##  P zoo               1.8-14   2025-04-10 [?] RSPM
    ## 
    ##  [1] /mnt/Data/Projects/RutePortugal/renv/library/linux-ubuntu-noble/R-4.5/x86_64-pc-linux-gnu
    ##  [2] /home/tomoe/.cache/R/renv/sandbox/linux-ubuntu-noble/R-4.5/x86_64-pc-linux-gnu/9a444a72
    ## 
    ##  * ── Packages attached to the search path.
    ##  P ── Loaded and on-disk path mismatch.
    ## 
    ## ──────────────────────────────────────────────────────────────────────────────
