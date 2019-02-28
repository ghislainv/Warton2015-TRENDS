##############
## Starter code for fitting joint models:
## (1) A multivariate GLMM using the lme4 package
## (2) A LVM using the boral package
###############

## We will use the spider dataset available in mvabund, pitfall trap counts of 12 hunting spider species at 28 sites. First load it from mvabund:
library(mvabund)
data(spider)

## Cherry-pick a few variables (lme4 can't handle many responses, especially when we only have 28 sites)
spid.abund = spider$abund[,c(1:3,7,8,12)]

##########
## (1) Fitting a multivariate GLMM using the lme4 package ##
## Use lme4 to fit a multivariate GLMM like in equation 1 of main text, i.e. (random) site effect, and a quadratic effect of soil for each spp. 
## Using a Poisson model with a log link, other options are available via the family argument, as usual.
## Note this function is super-fussy about data - the model didn't converge even for these six species!
##########   

spid.vec = c(as.matrix(spid.abund))
n.site = dim(spid.abund)[1]
n.spp = dim(spid.abund)[2]

## construct a data frame with soil (standardized), species labels, and site labels:
X = data.frame(soil = rep(scale(spider$x[,1]), n.spp), spp = rep(dimnames(spid.abund)[[2]], each=n.site), site = rep(dimnames(spid.abund)[[1]], n.spp) )

## fit the GLMM using lme4 and look at results
library(lme4)
fit.glmm = glmer(spid.vec~0+spp+spp:soil+spp:I(soil^2)+(1|site)+(0+spp|site), data=X, family=poisson())
# Returns some warnings about non-convergence -- the best solution to this is more samples and less variables!!

## The key term in the above is (0+spp|site), which introduces a multivariate random intercept at each site.
## This technique will work for count data but for binary or ordinal data, the variance of the random effect
## needs to be fixed (e.g. at one) for identifiability, which lme4 doesn't do - instead something like MCMCglmm might work.

print(summary(fit.glmm),correlation=FALSE) # To look at estimated parameter values
confint(fit.glmm, method="Wald") # 95% confidence intervals for model parameters. Very approximate though, other more computationally intensive options are available.

## Use corrplot package to calculate and plot residual correlations between species, e.g. possibly due to species interaction etc...
library(corrplot)
vrcorrs=VarCorr(fit.glmm)
corrs=attr(vrcorrs$site.1,"corr")
corrplot(corrs, diag = F, type = "lower", title = "Residual correlations from GLMM", method = "color", tl.srt = 45)



##########
## (2) Fitting a Latent Variable Model using the boral package ##
## Use boral to fit a LVM like in equation 3 of main text, i.e. two latent variables, (random) site effect, and a quadratic effect of soil. 
## Also note that this can fit models to more data (e.g. try full dataset) but takes longer to fit it - for full spdier it will still be a few min, for full aravo dataset over half an hour.
##########   

## Covariates need to be stored as a matrix for boral:
covX <- cbind(scale(spider$x[,1]),scale(spider$x[,1])^2)

## fit the LVM using boral and look at results.  Need version 0.7 or later, available on CRAN.
## NOTE: boral uses JAGS so you need to first install JAGS-3.0.0.exe or higher from http://sourceforge.net/projects/mcmc-jags/files/JAGS/

library(boral)
fit.lvm <- boral(y = spid.abund, X = covX, num.lv = 2, family = "poisson", row.eff = "random", save.model = TRUE, calc.ics = F, hypparams = c(20,20,20,20))
summary(fit.lvm) # To look at estimated parameter values
fit.lvm$hpdintervals # 95% credible intervals for model parameters.

## compare model coefficients with GLMM:
fixef(fit.glmm)
cbind(fit.lvm$lv.coefs.mean,fit.lvm$X.coefs.mean)
## the slopes are in general agreement, whereas the intercepts differ by quite a bit between the two models.
## The main reason is that in boral, the random intercepts are sampled are a non-zero normal distribution. 

## Dunn-Smyth residual plots  to check model assumption, outliers etc...
plot(fit.lvm)

## Residual ordination plot of the sites (please see Figure 2b in the main text for relevant explanation)
## Please note that boral currently does not automatically produce biplots, although species loadings can be obtained from fit.lvm$lv.coefs.median[,2:3]
lvsplot(fit.lvm)

## Use corrplot package to plot residual correlations between species, e.g. possibly due to species interaction etc...
res.cors <- get.residual.cor(fit.lvm)
corrplot(res.cors$rescor, diag = F, type = "lower", title = "Residual correlations from LVM", method = "color", tl.srt = 45)
## Note correlations are similar to those from the GLMM, the main difference being stronger correlations for the
## two points at the bottom-right of the plot. Adding an extra latent variable to the model would fix this (num.lv=3). 


## Please see help file for the boral function to find further information into the functions available (?boral)
## One nice trick is the get.enviro.cor function, to produce correlations DUE TO environmental variables
