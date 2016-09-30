library("Rcdm")

## load data
data("probability", package = "pks")
# install.packages("pks")

items <- sprintf("b1%.2i", 1:12)
resp <- probability[, items]
resp <- resp[complete.cases(resp),]

qmat <- read.table(header = TRUE, text = "
  cp id pb un
   0  0  1  0
   1  0  0  0
   0  0  0  1
   0  1  0  0
   1  0  1  0
   1  0  1  0
   0  0  1  1
   0  0  1  1
   0  1  1  0
   1  1  0  0
   1  1  1  0
   0  1  1  1
")
rownames(qmat) <- colnames(resp)

## fit model
mGDINA <- gdina(resp, qmat)

## coefficients
coef(mGDINA)

## standard errors
# correct computation of the information matrix
vhat <- vcov(mGDINA)
sqrt(diag(vhat))

## standard errors
# incorrect computation of the information matrix
vhat <- vcov(mGDINA, type = "partial") # ignore skill parameters
sqrt(diag(vhat))

## standard errors
# incorrect computation of the information matrix
vhat <- vcov(mGDINA, type = "itemwise") # itemwise
sqrt(diag(vhat))

## confidente intervals for item parameters
confint(mGDINA, alpha = 0.05)

## plot conditional response probabilities per item
plot(mGDINA)

## some fit indices
logLik(mGDINA)
BIC(mGDINA)
AIC(mGDINA)

## fit DINA model
mDINA <- gdina(resp, qmat, rule = "DINA")
plot(mDINA)

## fit DINO model
mDINO <- gdina(resp, qmat, rule = "DINO")
plot(mDINO)

## fit A-CDM model
mACDM <- gdina(resp, qmat, rule = "ACDM")
plot(mACDM)

## compare (nested) models
anova(mDINA, mACDM, mGDINA)
anova(mDINO, mACDM, mGDINA)

## compare (non-nested) models
AIC(mDINA, mDINO, mACDM, mGDINA)
BIC(mDINA, mDINO, mACDM, mGDINA)

## Wald-Test
gdina_wald(mGDINA)
