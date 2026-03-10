# epidoc_ageing_multimorb_depression
This repository presents the code for analysing the effect of ageing and multimorbidity in depression in the EpiDoC Cohort (Portugal) 

## The project
The project is a joint collaboration between Colab CRCH, NOVA and NOIS/PUC

## objective:
To conduct a Cross-Lagged Panel Model (CLPM) to capture the interplay among multimorbidity and depression in the EpiDoc Cohort data.
## paper
Not yet published
# data availability
The data (EpiDoc Chort data) is not included here, but can be made available upon reasonable request.
## Method:
The multivariate missing data were handled using multiple imputation in R’s mice package (version 3.18.0) with the CART (classification and regression trees) method. To address the possibility of non-randomness in missingness, the full dataset was used for imputation after excluding utility columns, such as storing weights or information on when and whether participants responded.
The algorithm was allowed to iterate over 40 iterations, with automatic variable predictor definition, in two parallel rounds with 10 imputations each, using parallel seeds of 500 and 700. It took over 3 weeks on Ubuntu 24.04.2 LTS, with an AMD Ryzen 7 7700 (16 threads), 32 GiB of RAM (and a 200GB swap file), and an NVIDIA GeForce RTX 3060. 

To address dropout cases naturally occurring in cohort studies, we decided to use only the data from participants who responded to all waves and to correct the non-response bias using inverse probability weighting.
Weight correction was made using inverted probability weighting, using R’s broom package. To compute the propensity weight, a Logit link model was used, considering depression (categorical levels: normal, mild, moderate, and severe), age categories, gender, multi-morbidity, family size, smoking, alcohol, exercise, employment status, and scholarship as covariates.


The imputed data, obtained with the CART method (multiple imputation using the mice package), were assessed for the structural equation model using lavaan.mi package (0.1-0) in R 4.3.3. Robust maximum likelihood estimation with Huber–White standard errors was used to account for non-normality and heteroskedasticity 

For the sequential adjustments applying constraints, the model differences were tested using the Yuan–Bentler scaled χ2 difference test, and, complementary, each model fit was assessed using the CFI (Comparative Fit Index), TFI (Tucker–Lewis Index), RMSEA (Root Mean Square Error of Approximation), and SRM (Standardised Root Mean Square Residual). 
Causal effects were estimated using standardised coefficients and their corresponding 95% confidence interval.  All the statistical t-tests for coefficients were two-sided, with p-values <0.05 considered significant. 

The longitudinal model applied is the Cross-Lagged Panel Model (CLPM). Formal testing for random-intercepts-CLPM, which is currently widely recommended, was not conducted. This decision was primarily due to the unequal intervals between waves, which were deemed potentially unsuitable for an RI-CLPM. Furthermore, an initial attempt to adjust for the equality of auto-regressive components was unsuccessful, leading to the exclusion of RI-CLPM as a viable option. 
