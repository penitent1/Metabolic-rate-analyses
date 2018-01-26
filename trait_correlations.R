## Trait correlations
library(ape)


## Load the anolis tree and dataset
anolis_phy <- read.tree("anolis.tre.R")
## Tip: Setting `row.names=1` assigns the rownames to be equal to the first column (in this case the species names)
anolis_dat <- read.csv("anolis.csv", row.names=1)

## Using gls to test for an association of body size with hind limb (dumb, i know but bear with me)

## Regular OLS
ols <- gls(AVG_hl ~ AVG_SVL, data=anolis_dat, method="ML")

## PGLS assuming a Brownian correlation
pgls_bm <- gls(AVG_hl ~ AVG_SVL, data=anolis_dat, correlation = corBrownian(phy=anolis_phy), method="ML")

## PGLS assuming a Pagel's lambda correlation
pgls_lambda <- gls(AVG_hl ~ AVG_SVL, data=anolis_dat, 
                  correlation = corPagel(value=1, phy=anolis_phy, fixed=FALSE), method="ML")

pgls_lambda

## Compare models using AIC (lower is beter)
AIC(ols)
AIC(pgls_bm)
AIC(pgls_lambda)

## Using a likelihood ratio test to evaluate whether pgls is justified
anova(pgls_lambda, ols)

## Is lambda significantly better than Brownian motion
anova(pgls_bm, pgls_lambda)
