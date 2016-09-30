### R code from vignette source 'demo.Rnw'

###################################################
### code chunk number 1: demo.Rnw:22-23 (eval = FALSE)
###################################################
## install.packages(c("Rcpp", "RcppArmadillo"))


###################################################
### code chunk number 2: demo.Rnw:27-28 (eval = FALSE)
###################################################
## install.packages("path/filename.tar.gz", repos = NULL, type = "source")


###################################################
### code chunk number 3: demo.Rnw:33-35 (eval = FALSE)
###################################################
## install.packages("devtools") # if you have not installed devools
## devtools::install_github("mphili/cdm")


###################################################
### code chunk number 4: demo.Rnw:40-41 (eval = FALSE)
###################################################
## install.packages(c("Rcpp", "Matrix", "limSolve", "gtools"))


###################################################
### code chunk number 5: demo.Rnw:45-46
###################################################
library("Rcdm")


###################################################
### code chunk number 6: demo.Rnw:50-67
###################################################
# install.packages("pks")
data("probability", package = "pks")

# reduce to 12 items from the first test
items <- sprintf("b1%.2i", 1:12)
resp <- probability[, items]
resp <- resp[complete.cases(resp),]

qmat <- t(read.table(header = FALSE, text = "
  0    1    0    0    1    1    0    0    0    1    1    0
  0    0    0    1    0    0    0    0    1    1    1    1
  1    0    0    0    1    1    1    1    1    0    1    1
  0    0    1    0    0    0    1    1    0    0    0    1
"))

colnames(qmat) <- c("cp", "id", "pb", "un")
rownames(qmat) <- colnames(resp)


###################################################
### code chunk number 7: demo.Rnw:71-72
###################################################
mGDINA <- gdina(resp, qmat)


###################################################
### code chunk number 8: demo.Rnw:76-77
###################################################
plot(mGDINA)


###################################################
### code chunk number 9: demo.Rnw:85-86
###################################################
mGDINA$dj


###################################################
### code chunk number 10: demo.Rnw:90-91
###################################################
mGDINA$pa


###################################################
### code chunk number 11: demo.Rnw:95-96
###################################################
mGDINA$pj


###################################################
### code chunk number 12: demo.Rnw:100-101
###################################################
coef(mGDINA)


###################################################
### code chunk number 13: demo.Rnw:108-109
###################################################
v0 <- vcov(mGDINA)


###################################################
### code chunk number 14: demo.Rnw:114-115
###################################################
v1 <- vcov(mGDINA, type = "partial")


###################################################
### code chunk number 15: demo.Rnw:120-121
###################################################
v2 <- vcov(mGDINA, type = "itemwise")


###################################################
### code chunk number 16: demo.Rnw:126-129
###################################################
cbind("full (correct)"       = sqrt(diag(v0))[seq(1, 2*nrow(qmat))],
      "partial (incorrect)"  = sqrt(diag(v1)),
      "itemwise (incorrect)" = sqrt(diag(v2)))


###################################################
### code chunk number 17: demo.Rnw:134-135
###################################################
sqrt(diag(v1))[-seq(1, 2*nrow(qmat))]


###################################################
### code chunk number 18: demo.Rnw:142-143
###################################################
confint(mGDINA, alpha = 0.05)


###################################################
### code chunk number 19: demo.Rnw:149-150
###################################################
gdina_wald(mGDINA)


###################################################
### code chunk number 20: demo.Rnw:155-158
###################################################
logLik(mGDINA)
BIC(mGDINA)
AIC(mGDINA)


###################################################
### code chunk number 21: demo.Rnw:164-165
###################################################
mDINA <- gdina(resp, qmat, rule = "DINA")


###################################################
### code chunk number 22: demo.Rnw:169-170
###################################################
mDINO <- gdina(resp, qmat, rule = "DINO")


###################################################
### code chunk number 23: demo.Rnw:174-175
###################################################
mACDM <- gdina(resp, qmat, rule = "ACDM")


###################################################
### code chunk number 24: demo.Rnw:181-183
###################################################
AIC(mDINA, mDINO, mACDM, mGDINA)
BIC(mDINA, mDINO, mACDM, mGDINA)


###################################################
### code chunk number 25: demo.Rnw:187-189
###################################################
anova(mDINA, mACDM, mGDINA)
anova(mDINO, mACDM, mGDINA)


