Robust backfitting
================
Matias Salibian
2020-07-04

## A robust backfitting algorithm

The `R` package `RBF` (available on CRAN
[here](https://cran.r-project.org/package=RBF)) implements the robust
back-fitting algorithm as proposed by Boente, Martinez and
Salibian-Barrera in

> Boente G, Martinez A, Salibian-Barrera M. (2017) Robust estimators for
> additive models using backfitting. Journal of Nonparametric
> Statistics. Taylor & Francis; 29, 744-767.
> [DOI: 10.1080/10485252.2017.1369077](https://doi.org/10.1080/10485252.2017.1369077)

This repository contains a development version of `RBF` which may differ
slightly from the one available on CRAN (until the CRAN version is
updated appropriately).

The package in this repository can be installed from within `R` by using
the following code (assuming the
[devtools](https://cran.r-project.org/package=devtools)) package is
available:

``` r
devtools::install_github("msalibian/RBF")
```

### An example

Here is a (longish) example on how `RBF` works. We use the Air Quality
data.

``` r
library(RBF)
data(airquality)
```

A scatter plot of the data

``` r
pairs(airquality[, c('Ozone', 'Solar.R', 'Wind', 'Temp')], 
      pch=19, col='gray30', cex=1.5)
```

![](README_files/figure-gfm/scatter-1.png)<!-- -->

The following bandwidths were obtained via a robust leave-one-out
cross-validation procedure (described in the paper). Here we just set
them to their optimal values:

``` r
bandw <- c(136.728453,   8.894283,   4.764985)
```

Now we use the robust backfitting algorithm to fit an additive model
using Tukey’s bisquare loss (the default tuning constant for this loss
function is 4.685). We remove cases with missing entries.

``` r
ccs <- complete.cases(airquality)
fit.full <- backf.rob(Ozone ~ Solar.R + Wind + Temp, data=airquality,
                subset=ccs, windows=bandw, degree=1, type='Tukey')
```

We display the 3 fits (one per additive component), being careful with
the axis limits (which are stored to use them later):

``` r
lim.cl <- lim.rob <- matrix(0, 2, 3)
x0 <- fit.full$Xp
for(j in 1:3) {
  re <- fit.full$y - fit.full$alpha - rowSums(fit.full$g.matrix[,-j])
  lim.rob[,j] <- c(min(re), max(re))
  plot(re ~ x0[,j], type='p', pch=19, col='gray30', 
       xlab=colnames(x0)[j], ylab='', cex=1.5)
  oo <- order(x0[,j])
  lines(x0[oo,j], fit.full$g.matrix[oo,j], lwd=5, col='blue')
}
```

<img src="README_files/figure-gfm/showfits-1.png" width="33%" /><img src="README_files/figure-gfm/showfits-2.png" width="33%" /><img src="README_files/figure-gfm/showfits-3.png" width="33%" />

NOTE: These plots could also be obtained using the `plot` method:

``` r
# Plot each component
plot(fit.full, which=1:3)
```

We now compute and display the classical backfitting fits, with
bandwidths chosen via leave-one-out CV. Below are plots of partial
residuals with the classical and robust fits overlaid:

``` r
aircomplete <- airquality[ complete.cases(airquality), ]
library(gam)
fit.gam <- gam(Ozone ~ lo(Solar.R, span=.7) + lo(Wind, span=.7)+
                 lo(Temp, span=.5), data=aircomplete)
# Plot both fits (robust and classical) 
x <- as.matrix( aircomplete[ , c('Solar.R', 'Wind', 'Temp')] )
y <- as.vector( aircomplete[ , 'Ozone'] )
fits <- predict(fit.gam, type='terms')
# alpha.gam <- attr(fits, 'constant')
for(j in 1:3) {
  re <- fit.full$yp - fit.full$alpha - rowSums(fit.full$g.matrix[,-j])
  plot(re ~ x[,j], type='p', pch=20, col='gray45', xlab=colnames(x)[j], ylab='')
  oo <- order(x[,j])
  lines(x[oo,j], fit.full$g.matrix[oo,j], lwd=5, col='blue')
  lines(x[oo,j], fits[oo,j], lwd=5, col='magenta')
}
```

<img src="README_files/figure-gfm/classicfits-1.png" width="33%" /><img src="README_files/figure-gfm/classicfits-2.png" width="33%" /><img src="README_files/figure-gfm/classicfits-3.png" width="33%" />

To identify potential outiers we look at the residuals from the robust
fit, and use the function `boxplot`:

``` r
re.ro <- residuals(fit.full)
ou.ro <- boxplot(re.ro, col='gray80')$out
in.ro <- (1:length(re.ro))[ re.ro %in% ou.ro ]
points(rep(1, length(in.ro)), re.ro[in.ro], pch=20, cex=1.5, col='red')
```

![](README_files/figure-gfm/outliers-1.png)<!-- -->

We highlight these suspicious observations on the scatter plot

``` r
cs <- rep('gray30', nrow(aircomplete))
cs[in.ro] <- 'red'
os <- 1:nrow(aircomplete)
os2 <- c(os[-in.ro], os[in.ro])
pairs(aircomplete[os2, c('Ozone', 'Solar.R', 'Wind', 'Temp')], 
      pch=19, col=cs[os2], cex=1.5)
```

![](README_files/figure-gfm/showouts-1.png)<!-- -->

and on the partial residuals plots

``` r
for(j in 1:3) {
  re <- fit.full$yp - fit.full$alpha - rowSums(fit.full$g.matrix[,-j])
  plot(re ~ x[,j], type='p', pch=20, col='gray45', xlab=colnames(x)[j], ylab='')
  points(re[in.ro] ~ x0[in.ro,j], pch=19, col='red', cex=1.5)
  oo <- order(x[,j])
  lines(x[oo,j], fit.full$g.matrix[oo,j], lwd=5, col='blue')
  lines(x[oo,j], fits[oo,j], lwd=5, col='magenta')
}
```

<img src="README_files/figure-gfm/showouts2-1.png" width="33%" /><img src="README_files/figure-gfm/showouts2-2.png" width="33%" /><img src="README_files/figure-gfm/showouts2-3.png" width="33%" />

We now compute the classical backfitting algorithm on the data without
the potential outliers identified by the robust fit (the optimal
smoothing parameters were computed using leave-one-out
cross-validation). Note that now both fits (robust and non-robust) are  
almost identical.

``` r
# Run the classical backfitting algorithm without outliers
airclean <- aircomplete[-in.ro, ]
fit.gam2 <- gam(Ozone ~ lo(Solar.R, span=.7) + lo(Wind, span=.8)+
                  lo(Temp, span=.3), data=airclean)
fits2 <- predict(fit.gam2, type='terms')
# alpha.gam2 <- attr(fits2, 'constant')
dd2 <- aircomplete[-in.ro, c('Solar.R', 'Wind', 'Temp')]
for(j in 1:3) {
  re <- fit.full$yp - fit.full$alpha - rowSums(fit.full$g.matrix[,-j])
  plot(re ~ x[,j], type='p', pch=20, col='gray45', xlab=colnames(x)[j], ylab='')
  points(re[in.ro] ~ x[in.ro,j], pch=20, col='red', cex=1.5)
  oo <- order(dd2[,j])
  lines(dd2[oo,j], fits2[oo,j], lwd=5, col='magenta')
  oo <- order(x[,j])
  lines(x[oo,j], fit.full$g.matrix[oo,j], lwd=5, col='blue')
}
```

<img src="README_files/figure-gfm/bothonclean-1.png" width="33%" /><img src="README_files/figure-gfm/bothonclean-2.png" width="33%" /><img src="README_files/figure-gfm/bothonclean-3.png" width="33%" />

Finally, we compare the prediction accuracy obtained with each of the
fits. Because we are not interested in predicting well any possible
outliers in the data, we evaluate the quality of the predictions using a
5%-trimmed mean squared prediction error (effectively measuring the
prediction accuracy on 95% of the data). We use this alpha-trimmed mean
squared function:

``` r
tms <- function(a, alpha=.1) {
  # alpha is the proportion to trim
  a2 <- sort(a^2)
  n0 <- floor( length(a) * (1 - alpha) )
  return( mean(a2[1:n0]) )
}
```

We use 100 runs of 5-fold CV to compare the 5%-trimmed mean squared
prediction error of the robust fit and the classical one.

``` r
dd <- airquality
dd <- dd[complete.cases(dd), c('Ozone', 'Solar.R', 'Wind', 'Temp')]
# 100 runs of K-fold CV
M <- 100
# 5-fold
K <- 5
n <- nrow(dd)
# store (trimmed) TMSPE for robust and gam, and also
tmspe.ro <- tmspe.gam <- vector('numeric', M)
set.seed(123)
ii <- (1:n)%%K + 1
for(runs in 1:M) {
  tmpro <- tmpgam <- vector('numeric', n)
  ii <- sample(ii)
  for(j in 1:K) {
    fit.full <- backf.rob(Ozone ~ Solar.R + Wind + Temp, 
                           point=dd[ii==j, -1], windows=bandw, 
                           epsilon=1e-6, degree=1, type='Tukey', 
                           subset = (ii!=j), data = dd)
    tmpro[ ii == j ] <- rowSums(fit.full$prediction) + fit.full$alpha
    fit.gam <- gam(Ozone ~ lo(Solar.R, span=.7) + lo(Wind, span=.7)+
                     lo(Temp, span=.5), data=dd[ii!=j, ])
    tmpgam[ ii == j ] <- predict(fit.gam, newdata=dd[ii==j, ], type='response')
  }
  tmspe.ro[runs] <- tms( dd$Ozone - tmpro, alpha=.05)
  tmspe.gam[runs] <- tms( dd$Ozone - tmpgam, alpha=.05)
}
```

These are the boxplots. We see that the robust fit consistently fits the
vast majority (95%) of the data better than the classical one.

``` r
boxplot(tmspe.ro, tmspe.gam, names=c('Robust', 'Classical'), 
        col=c('tomato3', 'gray80'), main='')
```

![](README_files/figure-gfm/pred3-1.png)<!-- -->
