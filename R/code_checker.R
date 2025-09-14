#--------------------------- step by step check ------------------------------#
load("data/ExampleData.RData")
source("R/utils.R")

source("R/coxkl.R")
source("R/cv.coxkl.R")

source("R/coxkl_ridge.R")
source("R/cv.coxkl_ridge.R")

source("R/coxkl_highdim.R")
source("R/cv.coxkl_highdim.R")

Rcpp::sourceCpp("src/utils.cpp")
Rcpp::sourceCpp("src/KLCox.cpp")
Rcpp::sourceCpp("src/KLCox_highdim.cpp")


#-------------------- 1. check: compare with coxph under eta = 0 --------------------------# 
dat <- data.frame(time = ExampleData$time,
                  status = ExampleData$status,
                  stratum = ExampleData$stratum,
                  ExampleData$z)

library(survival)
coxph_res <- coxph(
  Surv(time, status) ~ Z1 + Z2 + Z3 + Z4 + Z5 + Z6 + strata(stratum),
  data = dat,
  control = coxph.control(eps = 1e-7, iter.max = 100)
)
coxph_res$coefficients
#        Z1         Z2         Z3         Z4         Z5         Z6 
# 0.1586184 -0.2285272  0.4641955 -0.3864544  0.3095823 -0.2555455 
coxph_res$loglik[2]
# [1] -1202.085



result1 <- coxkl(z = ExampleData$z,
                 delta = ExampleData$status,
                 time = ExampleData$time,
                 stratum = ExampleData$stratum,
                 RS = NULL,
                 beta = ExampleData$beta_external,
                 etas = 0,
                 Mstop = 50,
                 backtrack = F,
                 message = F)
result1$beta
#               0
# [1,]  0.1586184
# [2,] -0.2285272
# [3,]  0.4641955
# [4,] -0.3864544
# [5,]  0.3095823
# [6,] -0.2555455
result1$likelihood
#         0 
# -1202.085 


#-------------------- 2. check: compare with Lingfeng's original codes --------------------------# 
load("data/ExampleDataHighDim.RData")
result_new <- coxkl(z = ExampleData$z,
                    delta = ExampleData$status,
                    time = ExampleData$time,
                    RS = NULL,
                    beta = ExampleData$beta_external,
                    etas = seq(0, 5, 1),
                    message = T)


head(result_new$beta)  # the result is the same with that form Lingfeng's original code
#                0           1           2           3           4           5
# [1,]  0.33158984  0.30507909  0.29725498  0.29348346  0.29126065  0.28979449
# [2,] -0.45556536 -0.38588179 -0.36236166 -0.35044965 -0.34324199 -0.33840898
# [3,] -0.07571637 -0.04594781 -0.03470755 -0.02884734 -0.02525266 -0.02282281
# [4,] -0.07332686  0.00575570  0.03151854  0.04437497  0.05208957  0.05723476
# [5,]  0.62690796  0.54722031  0.52055560  0.50720223  0.49918453  0.49383736
# [6,]  0.19092080  0.09606307  0.06455724  0.04881102  0.03936478  0.03306787



#-------------------- 3. check: cv.coxkl --------------------------# 
load("data/ExampleData.RData")

etas <- generate_eta(method = "exponential", n = 10, max_eta = 5)
etas <- sample(etas) ## suffle eta orders
cv_res1 <- cv.coxkl(z = ExampleData$z,
                    delta = ExampleData$status,
                    time = ExampleData$time,
                    stratum = NULL,
                    RS = NULL,
                    beta = ExampleData$beta_external,
                    etas = etas,
                    nfolds = 5,
                    criteria = c("LinPred"),   #"V&VH", "LinPred", "CIndex_pooled", "CIndex_foldaverage"
                    message = T)
cv_res1
# $internal_stat
#           eta LinPred_Loss
# 1  0.00000000     2086.045
# 2  0.03374245     2077.520
# 3  0.09002825     2065.255
# 4  0.18391863     2049.000
# 5  0.34053721     2029.814
# 6  0.60179276     2010.365
# 7  1.03759328     1993.975
# 8  1.76455236     1982.737
# 9  2.97719318     1976.535
# 10 5.00000000     1973.824
# 
# $external_stat
# [1] 1973.572

plot(cv_res1$internal_stat[,1], cv_res1$internal_stat[,2], type = "b", xlab = "eta")



#---------------------------4. ridge (coxkl_ridge)------------------------------#
load("data/ExampleDataHighDim.RData")
result3 <- coxkl_ridge(z = ExampleData$z,
                       delta = ExampleData$status,
                       time = ExampleData$time,
                       stratum = NULL,
                       RS = NULL,
                       beta = ExampleData$beta_external,
                       message = T)
#   Cross-validation over lambda sequence:
#   |==============================| 100%

head(result3$beta[, 1:5])
#         138.1936    128.8798    120.1937    112.0931    104.5384
# [1,]  0.13888135  0.14450653  0.15017707  0.15588048  0.16160406
# [2,] -0.17547848 -0.18241304 -0.18944000 -0.19654779 -0.20372421
# [3,] -0.04192810 -0.04312603 -0.04429336 -0.04542704 -0.04652448
# [4,]  0.03213009  0.03176027  0.03127296  0.03066580  0.02993718
# [5,]  0.23628338  0.24513399  0.25408686  0.26312975  0.27225002
# [6,]  0.09007171  0.09307455  0.09608820  0.09910679  0.10212425


#---------------------------5. cv.coxkl_ridge ------------------------------#
set.seed(1)
etas <- generate_eta(method = "exponential", n = 10, max_eta = 5)
etas <- sample(etas) ## suffle eta orders
cv.result3 <- cv.coxkl_ridge(z = ExampleData$z,
                             delta = ExampleData$status,
                             time = ExampleData$time,
                             stratum = NULL,
                             RS = NULL,
                             beta = ExampleData$beta_external,
                             etas = etas,
                             nfolds = 5, 
                             cv.criteria = "CIndex_pooled",
                             message = T)
cv.result3$external_stat.best_per_eta
#           eta    lambda CIndex_pooled
# 1  0.00000000 64.143778     0.7174664
# 2  0.03374245 52.106135     0.7191440
# 3  0.09002825 52.224183     0.7204656
# 4  0.18391863 49.408333     0.7223465
# 5  0.34053721 27.022552     0.7260573
# 6  0.60179276  7.926052     0.7296157
# 7  1.03759328 10.082701     0.7330216
# 8  1.76455236 11.133753     0.7353599
# 9  2.97719318  2.826165     0.7371899
# 10 5.00000000  5.030473     0.7359699

cv.result3$external_stat
# [1] 0.7345466

plot(cv.result3$best_per_eta[,1], cv.result3$best_per_eta[,3], type = "b", xlab = "eta")

#---------------------------6.coxkl_highdim------------------------------#
load("data/ExampleDataHighDim.RData")
result4 <- coxkl_highdim(z = ExampleData$z,
                         delta = ExampleData$status,
                         time = ExampleData$time,
                         stratum = NULL,
                         RS = NULL,
                         beta = ExampleData$beta_external,
                         eta = 0,
                         alpha = 1.0,
                         message = T)


#---------------------------7. cv.coxkl_highdim ------------------------------#
set.seed(1)
etas <- generate_eta(method = "exponential", n = 10, max_eta = 5)
etas <- sample(etas)
result5 <- cv.coxkl_highdim(z = ExampleData$z,
                            delta = ExampleData$status,
                            time = ExampleData$time,
                            stratum = NULL,
                            RS = NULL,
                            beta = ExampleData$beta_external,
                            etas = etas,
                            alpha = 1.0,
                            nfolds = 5, 
                            cv.criteria = "V&VH",
                            message = T)
result5$best_per_eta
#           eta      lambda     Loss
# 1  0.00000000 0.028426151 2350.817
# 2  0.03374245 0.026549628 2349.435
# 3  0.09002825 0.024816367 2347.365
# 4  0.18391863 0.020393447 2344.343
# 5  0.34053721 0.016952548 2340.283
# 6  0.60179276 0.017459028 2336.117
# 7  1.03759328 0.014612411 2332.514
# 8  1.76455236 0.009900687 2329.408
# 9  2.97719318 0.007674849 2327.024
# 10 5.00000000 0.005514951 2325.576

result5$external_stat
# 2325.419



