
# Load data and set up response vector and design matrix
data(airquality)

# # Show data in pairwise plots
 png('ScatterPlot.png', bg='transparent')
pairs(airquality[, c('Ozone', 'Solar.R', 'Wind', 'Temp')], pch=19, col='gray30', cex=1.5)
 dev.off()

# Load RBF package
# devtools::install_github("msalibian/RBF")
library(RBF)

# Bandwidth selection with leave-one-out cross-validation
## Without outliers
# This takes a long time to compute (approx 380 minutes running
# R 3.6.1 on an Intel(R) Core(TM) i7-4790 CPU @ 3.60GHz)
ccs <- complete.cases(airquality)
aircomplete <- airquality[ccs, c('Ozone', 'Solar.R', 'Wind', 'Temp')]
a <- c(1/2, 1, 1.5, 2, 2.5, 3)
h1 <- a * sd(aircomplete[,2])
h2 <- a * sd(aircomplete[,3])
h3 <- a * sd(aircomplete[,4])
hh <- expand.grid(h1, h2, h3)
nh <- nrow(hh)
rmspe <- rep(NA, nh)
jbest <- 0
cvbest <- +Inf
n <- nrow(aircomplete)
# Optimal bandwidths
# Var1     Var2     Var3
# 33 136.7285 10.67314 4.764985

bandw <- c(136.7285, 10.67314, 4.764985)

# Compute the robust backfitting with the optimal bandwidths
fit.full <- backf.rob(Ozone ~ Solar.R + Wind + Temp, windows=bandw, 
                      degree=1, type='Tukey', subset = ccs, 
                      data=airquality)

 # Classical backfitting
 # Find optimal spans using leave-one-out cross validation
 library(gam)
 aircomplete <- airquality[ccs, c('Ozone', 'Solar.R', 'Wind', 'Temp')]
 a <- c(.3, .4, .5, .6, .7, .8, .9)
 hh <- expand.grid(a, a, a)
 nh <- nrow(hh)
 jbest <- 0
 cvbest <- +Inf
 n <- nrow(aircomplete)
 # Fit the backfitting algorithm with the optimal spans found above
 fit.gam <- gam(Ozone ~ lo(Solar.R, span=.7) + lo(Wind, span=.7)+
                  lo(Temp, span=.5), data=aircomplete)
 
 # Plot both fits (robust and classical) 
 x <- as.matrix( aircomplete[ , c('Solar.R', 'Wind', 'Temp')] )
 y <- as.vector( aircomplete[ , 'Ozone'] )
 fits <- predict(fit.gam, type='terms')
 par(mfrow=c(2,2))
 
 png('Figure-ozone-todos.png', bg='transparent', width = 720, height = 420, units="px") #480, 280
 par(mfrow=c(1,3))
 for(j in 1:3) {
   re <- fit.full$yp - fit.full$alpha - rowSums(fit.full$g.matrix[,-j])
   # png(paste0('Figure-ozone-res-h-g',j,'.png'), bg='transparent')
   plot(re ~ x[,j], type='p', pch=20, col='gray45', xlab=colnames(x)[j], ylab='')
   oo <- order(x[,j])
   lines(x[oo,j], fit.full$g.matrix[oo,j], lwd=5, col='blue', lty=1)
   lines(x[oo,j], fits[oo,j], lwd=5, col='magenta', lty=2)
   # dev.off()
 }
 
 dev.off()
 par(mfrow=c(1,1)) 
 
# Identify possible outliers
# first compute residuals
re.ro <- residuals(fit.full)
# use the function boxplot() to plot and identify potential outliers
 png('Figure-ozone-boxplot.png', bg='transparent')
ou.ro <- boxplot(re.ro, col='gray80')$out
# determine their indices
in.ro <- (1:length(re.ro))[ re.ro %in% ou.ro ]
# highlight potential outliers on the residual
# boxplot with red circles
points(rep(1, length(in.ro)), re.ro[in.ro], pch=20, cex=3, col='red')
 dev.off()
(in.ro)

# Identify potential outliers on the pairwise scatterplot
# using red circles
aircomplete <- airquality[ complete.cases(airquality), ]
cs <- rep('gray30', nrow(aircomplete))
cs[in.ro] <- 'red'
os <- 1:nrow(aircomplete)
os2 <- c(os[-in.ro], os[in.ro])
 #png('Figure-ozone-scat-h.png', bg='transparent')
pairs(aircomplete[os2, c('Ozone', 'Solar.R', 'Wind', 'Temp')], 
      pch=19, col=cs[os2], cex=1.5)
 #dev.off()


# Plot both fits (robust and classical) 
x <- as.matrix( aircomplete[ , c('Solar.R', 'Wind', 'Temp')] )
y <- as.vector( aircomplete[ , 'Ozone'] )
fits <- predict(fit.gam, type='terms')
par(mfrow=c(2,2))

#png('Figure-ozone-res-h-g-todos.png', bg='transparent', width = 720, height = 420, units="px") #480, 280
par(mfrow=c(1,3))
for(j in 1:3) {
  re <- fit.full$yp - fit.full$alpha - rowSums(fit.full$g.matrix[,-j])
  # png(paste0('Figure-ozone-res-h-g',j,'.png'), bg='transparent')
  plot(re ~ x[,j], type='p', pch=20, col='gray45', xlab=colnames(x)[j], ylab='')
  points(re[in.ro] ~ x[in.ro,j], pch=20, col='red', cex=2)
  oo <- order(x[,j])
  lines(x[oo,j], fit.full$g.matrix[oo,j], lwd=5, col='blue', lty=1)
  lines(x[oo,j], fits[oo,j], lwd=5, col='magenta', lty=2)
  # dev.off()
}
#dev.off()


# Fit gam without outliers
# Find optimal bandwidths
airclean <- aircomplete[-in.ro, c('Ozone', 'Solar.R', 'Wind', 'Temp')]
a <- c(.3, .4, .5, .6, .7, .8, .9)
hh <- expand.grid(a, a, a)
nh <- nrow(hh)
jbest <- 0
cvbest <- +Inf
n <- nrow(airclean)
# # Var1 Var2 Var3
# # 40  0.7  0.8  0.3

fit.gam2 <- gam(Ozone ~ lo(Solar.R, span=.7) + lo(Wind, span=.8)+
                  lo(Temp, span=.3), data=airclean) #aircomplete[-in.ro, ])

# Plot both fits (robust and classical computed w/o outliers) 
fits2 <- predict(fit.gam2, type='terms')
dd2 <- aircomplete[-in.ro, c('Solar.R', 'Wind', 'Temp')]
par(mfrow=c(2,2))
png('Figure-ozone-out-cla-rob.png', bg='transparent', width = 720, height = 420, units="px") #480, 280
par(mfrow=c(1,3))
for(j in 1:3) {
  re <- fit.full$yp - fit.full$alpha - rowSums(fit.full$g.matrix[,-j])
  # png(paste0('Figure-ozone-out-cla-rob-g',j,'.png'), bg='transparent')
  plot(re ~ x[,j], type='p', pch=20, col='gray45', xlab=colnames(x)[j], ylab='')
  points(re[in.ro] ~ x[in.ro,j], pch=20, col='red', cex=2)
  oo <- order(dd2[,j])
  lines(dd2[oo,j], fits2[oo,j], lwd=5, col='magenta', lty=2)
  oo <- order(x[,j])
  lines(x[oo,j], fit.full$g.matrix[oo,j], lwd=5, col='blue', lty=1)
   #dev.off()
}
dev.off()
par(mfrow=c(1,1))
