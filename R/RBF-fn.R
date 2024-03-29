#' Derivative of Tukey's bi-square loss function.
#'
#' This function evaluates the first derivative of Tukey's bi-square loss function.
#'
#' This function evaluates the first derivative of Tukey's bi-square loss function.
#'
#' @param r a vector of real numbers
#' @param k a positive tuning constant.
#'
#' @return A vector of the same length as \code{x}.
#'
#' @author Matias Salibian-Barrera, \email{matias@stat.ubc.ca}, Alejandra Martinez
#'
#' @examples
#' x <- seq(-2, 2, length=10)
#' psi.tukey(r=x, k = 1.5)
#'
#' @export
#' @import stats graphics
#' @useDynLib RBF, .registration = TRUE
psi.tukey <- function(r, k=4.685){
  u <- abs(r/k)
  w <- r*((1-u)*(1+u))^2
  w[u>1] <- 0
  return(w)
}


# Tukey's weight function "Psi(r)/r"
psi.w <- function(r, k= 4.685){
  u <- abs(r/k)
  w <- ((1 + u) * (1 - u))^2
  w[u > 1] <- 0
  return(w)
}

#Tukey's loss function
rho.tukey <- function(r, k=4.685){
  if(abs(r)<=k){
    return( 1-(1-(r/k)^2)^3 )
  }else{
    return(1)
  }
}

#' Derivative of Huber's loss function.
#'
#' This function evaluates the first derivative of Huber's loss function.
#'
#' This function evaluates the first derivative of Huber's loss function.
#'
#' @param r a vector of real numbers
#' @param k a positive tuning constant.
#'
#' @return A vector of the same length as \code{x}.
#'
#' @author Matias Salibian-Barrera, \email{matias@stat.ubc.ca}, Alejandra Martinez
#'
#' @examples
#' x <- seq(-2, 2, length=10)
#' psi.huber(r=x, k = 1.5)
#'
#' @export
psi.huber <- function(r, k=1.345)
  pmin(k, pmax(-k, r))

#Huber's weight function "Psi(r)/r"
psi.huber.w <- function(r, k=1.345)
  pmin(1, k/abs(r))


#Huber's loss function
rho.huber <- function(r, k=1.345){
  if(abs(r)<=k){
    return(r^2)
  }else{
    return(2*k*abs(r)-k^2)
  }
}

#' Epanechnikov kernel
#'
#' This function evaluates an Epanechnikov kernel
#'
#' This function evaluates an Epanechnikov kernel
#'
#' @param x a vector of real numbers
#'
#' @return A vector of the same length as \code{x} where each entry is
#' \code{0.75 * (1 - x^2)} if \code{x < 1} and 0 otherwise.
#'
#' @author Matias Salibian-Barrera, \email{matias@stat.ubc.ca}, Alejandra Martinez
#'
#' @examples
#' x <- seq(-2, 2, length=10)
#' k.epan(x)
#'
#' @export
k.epan<-function(x) {
  a <- 0.75*(1-x^2)
  tmp <- a*(abs(x)<=1)
  return(tmp)
}

# Euclidean norm
my.norm.2 <- function(x) sqrt(sum(x^2))


#' Classic Backfitting
#'
#' This function computes the standard backfitting algorithm for additive models.
#'
#' This function computes the standard backfitting algorithm for additive models,
#' using a squared loss function and local polynomial smoothers.
#'
#' @param formula an object of class \code{formula} (or one that can be coerced to 
#' that class): a symbolic description of the model to be fitted. 
#' @param data an optional data frame, list or environment (or object coercible 
#' by \link{as.data.frame} to a data frame) containing the variables in the model. 
#' If not found in \code{data}, the variables are taken from \code{environment(formula)}, 
#' typically the environment from which the function was called.
#' @param subset an optional vector specifying a subset of observations to be used in 
#' the fitting process.
#' @param point matrix of points where predictions will be computed and returned.
#' @param windows vector of bandwidths for the local polynomial smoother,
#' one per explanatory variable.
#' @param epsilon convergence criterion. Maximum allowed relative difference between
#' consecutive estimates
#' @param degree degree of the local polynomial smoother. Defaults to \code{0} (local constant).
#' @param prob vector of probabilities of observing each response (length n).
#' Defaults to \code{NULL} and in that case it is ignored.
#' @param max.it Maximum number of iterations for the algorithm.
#'
#' @return A list with the following components:
#' \item{alpha}{Estimate for the intercept.}
#' \item{g.matrix }{Matrix of estimated additive components (n by p).}
#' \item{prediction }{Matrix of estimated additive components for the points listed in
#' the argument \code{point}.}
#'
#' @references Hasie, TJ and Tibshirani, RJ. Generalized Additive Models, 1990. Chapman
#' and Hall, London.
#'
#' @author Matias Salibian-Barrera, \email{matias@stat.ubc.ca}, Alejandra Martinez
#'
#' @examples
#' data(airquality)
#' tmp <- backf.cl(Ozone ~ Solar.R + Wind + Temp, data=airquality, 
#' subset=complete.cases(airquality), windows=c(130, 9, 10), degree=1)
#'
#' @export
backf.cl <- function(formula, data, subset, point=NULL, windows, epsilon=1e-6, degree=0,
                     prob=NULL, max.it=100) {
  
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms") # allow model.frame to update it
  yp <- model.response(mf, "numeric")
  Xp <- model.matrix(mt, mf, NULL) 
  # above line typically is "model.matrix(mt, mf, contrasts)" and contrasts equals NULL

  # remove the default intercept, it is estimated separately  
  if( all( Xp[, 1] == 1 ) ) Xp <- Xp[, -1]
  

  # AUX <- get_all_vars(formula)
  # yp <- AUX[,1]
  # Xp <- AUX[,-1]
  
  n <- length(yp)
  Xp <- as.matrix(Xp)
  q <- dim(Xp)[2]
  corte <- 10*epsilon
  corte.bis <- 10*epsilon

  # Remove cases with missing responses
  yp <- yp[ tmp <- (!is.na(yp)) ]
  Xp <- Xp[tmp, , drop=FALSE]
  n.miss <- length(yp)
  if(is.null(prob)){
    prob <- rep(1,n.miss)
  } else {
    prob <- prob[tmp]
  }

  alpha <- mean(yp)

  # Start the backfitting algorithm.
  g.matriz <- matrix(0,n.miss,q)
  it <- 0
  while( (corte > epsilon) & (it < max.it)) {
    g.matriz.aux <- g.matriz
    for(j in 1:q) {
      y.tilde.bis <- yp - alpha - rowSums(g.matriz[,-j,drop=FALSE])
      for(i in 1:n.miss) {
        we <- k.epan( (Xp[,j] - Xp[i,j]) / as.numeric(windows[j]) )
        if(degree == 0)
          g.matriz[i,j] <- sum( we * y.tilde.bis ) / sum( we )
        if(degree > 0) {
          tmp <- outer( as.vector( Xp[,j] - Xp[i,j]) , 1:degree, "^" )
          tmp <- cbind(rep(1,n.miss), tmp)
          g.matriz[i,j] <- solve( t(tmp * we) %*% tmp, t(tmp) %*% (we * y.tilde.bis))[1]
        }
      }
    }
    aux1 <- sum( colMeans(g.matriz) )
    g.matriz <- scale(g.matriz,center=TRUE,scale=FALSE)
    corte <- my.norm.2(rowSums(g.matriz.aux)-rowSums(g.matriz))
    it <- it + 1
  }

  prediccion <- NULL

  if(!is.null(point)){
    
    #if(!is.matrix(point)) { #is.null(dim(point))) {
    #  prediccion <- mpunto <- as.matrix(point) # matrix(point, byrow=TRUE, ncol=length(point))
    #} else {
    #  prediccion <- mpunto <- point
    #}
    
    if(is.null(dim(point))){
      if(q==1){
        prediccion <- mpunto <- as.matrix(point)
      }else{
        prediccion <- mpunto <- t(as.matrix(point))
      }
    } else {
      prediccion <- mpunto <- point
    }
    
    np <- dim(mpunto)[1]
    for(k in 1:np){
      for(j in 1:q){
        y.tilde.bis <- yp - alpha - rowSums(g.matriz[,-j,drop=FALSE]) #- aux1
        we <- k.epan( (Xp[,j] - mpunto[k,j]) / as.numeric(windows[j]) )
        if(degree == 0)
          prediccion[k,j] <- sum( we * y.tilde.bis ) / sum( we )
        if(degree > 0) {
          tmp <- outer( as.vector( Xp[,j] - mpunto[k,j]) , 1:degree, "^" )
          tmp <- cbind(rep(1,n.miss), tmp)
          prediccion[k,j] <- solve( t(tmp * we) %*% tmp, t(tmp) %*% (we * y.tilde.bis))[1]
        }
      }
    }
  }
  object <- list(alpha=alpha, g.matrix=g.matriz, prediction=prediccion, 
                 Xp=Xp, yp=yp, formula=formula, call=cl)
  class(object) <- c("backf.cl", "backf", "list")
  return(object)
}


#' Robust Backfitting
#'
#' This function computes a robust backfitting algorithm for additive models
#'
#' This function computes a robust backfitting algorithm for additive models
#' using robust local polynomial smoothers.
#'
#' @param formula an object of class \code{formula} (or one that can be coerced to 
#' that class): a symbolic description of the model to be fitted. 
#' @param data an optional data frame, list or environment (or object coercible 
#' by \link{as.data.frame} to a data frame) containing the variables in the model. 
#' If not found in \code{data}, the variables are taken from \code{environment(formula)}, 
#' typically the environment from which the function was called.
#' @param subset an optional vector specifying a subset of observations to be used in 
#' the fitting process.
#' @param point matrix of points where predictions will be computed and returned.
#' @param windows vector of bandwidths for the local polynomial smoother,
#' one per explanatory variable.
#' @param epsilon convergence criterion. Maximum allowed relative difference between
#' consecutive estimates
#' @param degree degree of the local polynomial smoother. Defaults to \code{0} (local constant).
#' @param prob vector of probabilities of observing each response (length n).
#' Defaults to \code{NULL} and in that case it is ignored.
#' @param sigma.hat estimate of the residual standard error. If \code{NULL} (default) we use the
#' \link{mad} of the residuals obtained with local medians.
#' @param max.it Maximum number of iterations for the algorithm.
#' @param k.h tuning constant for a Huber-type loss function.
#' @param k.t tuning constant for a Tukey-type loss function.
#' @param type one of either \code{'Tukey'} or \code{'Huber'}.
#'
#' @return A list with the following components:
#' \item{alpha}{Estimate for the intercept.}
#' \item{g.matrix }{Matrix of estimated additive components (n by p).}
#' \item{prediction }{Matrix of estimated additive components for the points listed in
#' the argument \code{point}.}
#' \item{sigma.hat }{Estimate of the residual standard error.}
#'
#' @author Matias Salibian-Barrera, \email{matias@stat.ubc.ca}, Alejandra Martinez
#'
#' @references Boente G, Martinez A, Salibian-Barrera M. Robust estimators
#' for additive models using backfitting. Journal of Nonparametric Statistics,
#' 2017; 29:744-767. https://doi.org/10.1080/10485252.2017.1369077
#'
#' @examples
#' data(airquality)
#' tmp <- backf.rob(Ozone ~ Solar.R + Wind + Temp, data=airquality, 
#' subset=complete.cases(airquality), windows=c(136.7, 8.9, 4.8), degree=1)
#'
#' @export
backf.rob <- function(formula, data, subset, windows, point=NULL, epsilon=1e-6, degree=0,
                      sigma.hat=NULL, prob=NULL, max.it=50, k.h=1.345,
                      k.t = 4.685, type='Huber'){
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms") # allow model.frame to update it
  yp <- model.response(mf, "numeric")
  Xp <- model.matrix(mt, mf, NULL) 
  # above line typically is "model.matrix(mt, mf, contrasts)" and contrasts equals NULL
  
  # remove the default intercept, it is estimated separately  
  if( all( Xp[, 1] == 1 ) ) Xp <- Xp[, -1]
  
  # AUX <- get_all_vars(formula)
  # yp <- AUX[,1]
  # Xp <- AUX[,-1]
  
  Xp <- as.matrix(Xp)
  n <- length(yp)
  q <- dim(Xp)[2]
  corte <- 10*epsilon
  corte.bis <- 10*epsilon

  # Remove observations with missing responses
  yp <- yp[ tmp<-!is.na(yp) ]
  XX <- Xp
  Xp <- Xp[tmp, , drop=FALSE]
  n.miss <- length(yp)
  if(is.null(prob)){prob <- rep(1,n.miss)
  } else {
    prob <- prob[tmp]
  }


  # Estimate residual standard error
  if(is.null(sigma.hat)){
    ab <- rep(0,n.miss)
    for(i in 1:n.miss){
      xtildebis <- scale(Xp, center=Xp[i, ], scale=windows)
      a <- matrix(as.numeric(abs(xtildebis) < 1), n.miss, q)
      a <- apply(a, 1, prod)
      a[ a == 0 ] <- NA
      ab[i] <- median( a*yp, na.rm = TRUE)
    }
    sigma.hat <- mad(yp - ab)
    if( sigma.hat < 1e-10 ) sigma.hat <- 1e-10 # sigma.hat <- sd(yp-ab,na.rm=TRUE)
  }


  alpha <- 0

  #Additive components start with zeroes
  g.matriz <- matrix(0,n.miss,q)
  it <- 0

  # Main loop
  while( (corte > epsilon) & (it < max.it) ){
    g.matriz.aux <- g.matriz
    for(j in 1:q) {
      # partial residuals
      y.tilde.bis <- yp - alpha - rowSums(g.matriz[,-j,drop=FALSE])
      if( any(is.na(y.tilde.bis))) {
        it <- (max.it + 1)
        break()
      }
      for(i in 1:n.miss){
        if( degree == 0 ){
          if( type=='Huber') {
            mu.ini <- median( y.tilde.bis[ abs(Xp[,j] - Xp[i,j]) < windows[j]] )
            if(!is.na(mu.ini)) {
              g.matriz[i,j] <- .C("kernel_huber_pos", as.double(Xp[i,j]), as.double(Xp[,j]), as.integer(n.miss),
                                  as.double(y.tilde.bis), as.double(mu.ini), as.double(windows[j]),
                                  as.double(epsilon), as.double(sigma.hat),
                                  as.double(prob), as.double(k.h), as.integer(max.it), salida=as.double(0), PACKAGE='RBF')$salida
            } else {
              g.matriz[i,j] <- NA
            }
          }
          if( type=='Tukey') {
            mu.ini <- median( y.tilde.bis[ abs(Xp[,j] - Xp[i,j]) < windows[j]] )
            if(!is.na(mu.ini)) {
              g.matriz[i,j] <- .C("kernel_tukey_pos", as.double(Xp[i,j]), as.double(as.matrix(Xp[,j])), as.integer(n.miss),
                                  as.double(y.tilde.bis), as.double(mu.ini), as.double(windows[j]),
                                  as.double(epsilon), as.double(sigma.hat),
                                  as.double(prob), as.double(k.t), as.integer(max.it), salida=as.double(0), PACKAGE='RBF')$salida
             } else {
              g.matriz[i,j] <- NA
            }
          }
        }
        if( degree > 0 ) {
          tmp <- outer( as.vector( Xp[,j] - Xp[i,j]) , 1:degree, "^" )
          tmp <- cbind(rep(1,n.miss), tmp)
          beta.ini <- rep(0, degree+1)
          if( type == 'Huber') {
            g.matriz[i,j] <- .C("kernel_huber_lin", as.double(Xp[i,j]), as.double(Xp[,j]), as.integer(n.miss),
                                as.double(y.tilde.bis), as.double(tmp), as.integer(degree),
                                as.double(beta.ini), as.double(windows[j]), as.double(epsilon), as.double(sigma.hat),
                                as.double(prob), as.double(k.h), as.integer(max.it),
                                salida=as.double(rep(0, degree+1)), PACKAGE='RBF')$salida[1]
          }
          if( type == 'Tukey') {
            beta.ini <- .C("kernel_huber_lin", as.double(Xp[i,j]), as.double(Xp[,j]), as.integer(n.miss),
                           as.double(y.tilde.bis), as.double(tmp), as.integer(degree),
                           as.double(beta.ini), as.double(windows[j]), as.double(epsilon), as.double(sigma.hat),
                           as.double(prob), as.double(k.h), as.integer(max.it),
                           salida=as.double(rep(0, degree+1)), PACKAGE='RBF')$salida
            if(!any(is.na(beta.ini))) {
              g.matriz[i,j] <- .C("kernel_tukey_lin", as.double(Xp[i,j]), as.double(Xp[,j]), as.integer(n.miss),
                                  as.double(y.tilde.bis), as.double(tmp), as.integer(degree),
                                  as.double(beta.ini), as.double(windows[j]), as.double(epsilon), as.double(sigma.hat),
                                  as.double(prob), as.double(k.t), as.integer(max.it),
                                  salida=as.double(rep(0, degree+1)), PACKAGE='RBF')$salida[1]
            } else {
              g.matriz[i,j] <- NA
            }
          }
        }
      }
    }
    aux1 <- sum(colMeans(g.matriz))
    g.matriz <- scale(g.matriz,center=TRUE,scale=FALSE)
    corte <- my.norm.2(rowSums(g.matriz.aux)-rowSums(g.matriz))
    it <- it + 1
    # update intercept
    y.tilde.bis <- yp - rowSums(g.matriz)
    if(!any(is.na(y.tilde.bis))) {
      if( type=='Huber') {
        mu.ini <- median( y.tilde.bis )
        alpha <- .C("huber_pos", as.integer(n.miss),
                    as.double(y.tilde.bis), as.double(mu.ini),
                    as.double(epsilon), as.double(sigma.hat),
                    as.double(prob), as.double(k.h), as.integer(max.it),
                    salida=as.double(0), PACKAGE='RBF')$salida
      }
      if( type=='Tukey') {
        mu.ini <- median( y.tilde.bis )
        mu.ini <- .C("huber_pos", as.integer(n.miss),
                     as.double(y.tilde.bis), as.double(mu.ini),
                     as.double(epsilon), as.double(sigma.hat),
                     as.double(prob), as.double(k.h), as.integer(max.it),
                     salida=as.double(0), PACKAGE='RBF')$salida
        alpha <- .C("tukey_pos", as.integer(n.miss),
                    as.double(y.tilde.bis), as.double(mu.ini),
                    as.double(epsilon), as.double(sigma.hat),
                    as.double(prob), as.double(k.t), as.integer(max.it),
                    salida=as.double(0), PACKAGE='RBF')$salida
      }
    } else { stop('Error computing additive components'); break()
    }
  }

  prediccion <- NULL

  if(!is.null(point)){
    
    
    #if(!is.matrix(point)) { #is.null(dim(point))) {
    #  prediccion <- mpunto <-  as.matrix(point) #matrix(point, byrow=TRUE, ncol=length(point))
    #} else {
    #  prediccion <- mpunto <- point
    #}
    
    if(is.null(dim(point))){
      if(q==1){
        prediccion <- mpunto <- as.matrix(point)
      }else{
        prediccion <- mpunto <- t(as.matrix(point))
      }
    } else {
      prediccion <- mpunto <- point
    }
    
    np <- dim(mpunto)[1]
    for(k in 1:np){
      for(j in 1:q){
        y.tilde.bis <- yp - alpha - rowSums(g.matriz[,-j,drop=FALSE])
        if( degree == 0 ){
          mu.ini <- median( y.tilde.bis[ abs(Xp[,j] - mpunto[k,j]) < windows[j]] )
          if(!is.na(mu.ini)) {
            if( type=='Huber') {
              prediccion[k,j] <- .C("kernel_huber_pos", as.double(mpunto[k,j]), as.double(Xp[,j]), as.integer(n.miss),
                                    as.double(y.tilde.bis), as.double(mu.ini), as.double(windows[j]),
                                    as.double(epsilon), as.double(sigma.hat),
                                    as.double(prob), as.double(k.h), as.integer(max.it),salida=as.double(0), PACKAGE='RBF')$salida
            }
            if( type=='Tukey') {
              prediccion[k,j] <- .C("kernel_tukey_pos", as.double(mpunto[k,j]), as.double(as.matrix(Xp[,j])), as.integer(n.miss),
                                    as.double(y.tilde.bis), as.double(mu.ini), as.double(windows[j]),
                                    as.double(epsilon), as.double(sigma.hat),
                                    as.double(prob), as.double(k.t), as.integer(max.it),salida=as.double(0), PACKAGE='RBF')$salida
            }
          } else {
            prediccion[k,j] <- NA
          }
        }
        if( degree > 0 ) {
          tmp <- outer( as.vector( Xp[,j] - mpunto[k,j]) , 1:degree, "^" )
          tmp <- cbind(rep(1,n.miss), tmp)
          beta.ini <- rep(0, degree+1)
          if( type == 'Huber') {
            prediccion[k,j] <- .C("kernel_huber_lin", as.double(mpunto[k,j]), as.double(Xp[,j]), as.integer(n.miss),
                                  as.double(y.tilde.bis), as.double(tmp), as.integer(degree),
                                  as.double(beta.ini), as.double(windows[j]), as.double(epsilon), as.double(sigma.hat),
                                  as.double(prob), as.double(k.h), as.integer(max.it),
                                  salida=as.double(rep(0, degree+1)), PACKAGE='RBF')$salida[1]
          }
          if( type == 'Tukey') {
            beta.ini <- .C("kernel_huber_lin", as.double(mpunto[k,j]), as.double(Xp[,j]), as.integer(n.miss),
                           as.double(y.tilde.bis), as.double(tmp), as.integer(degree),
                           as.double(beta.ini), as.double(windows[j]), as.double(epsilon), as.double(sigma.hat),
                           as.double(prob), as.double(k.h), as.integer(max.it),
                           salida=as.double(rep(0, degree+1)), PACKAGE='RBF')$salida
            if(!any(is.na(beta.ini))) {
              prediccion[k,j] <- .C("kernel_tukey_lin", as.double(mpunto[k,j]), as.double(Xp[,j]), as.integer(n.miss),
                                    as.double(y.tilde.bis), as.double(tmp), as.integer(degree),
                                    as.double(beta.ini), as.double(windows[j]), as.double(epsilon), as.double(sigma.hat),
                                    as.double(prob), as.double(k.t), as.integer(max.it),
                                    salida=as.double(rep(0, degree+1)), PACKAGE='RBF')$salida[1]
            } else {
              prediccion[k,j] <- NA
            }
          }
        }
      }
    }
  }
  object <- list(alpha=alpha, g.matrix=g.matriz, sigma.hat=sigma.hat, prediction=prediccion, 
                 type=type, Xp=Xp, yp=yp, formula=formula, call=cl)
  class(object) <- c("backf.rob", "backf", "list")
  return(object)
  # return(list(alpha=alpha,g.matrix=g.matriz, sigma.hat=sigma.hat, prediction=prediccion))
}


# #' Cross-validation for Robust Backfitting
# #'
# #' This function performs one run of K-fold cross-validation using the
# #' robust backfitting algorithm.
# #'
# #' This function performs one run of K-fold cross-validation using the
# #' robust backfitting algorithm and returns a robust measure of the
# #' hold out prediction error.
# #'
# #' @param k a positive integer indicating the number of folds.
# #' @param Xp a matrix (n x p) containing the explanatory variables
# #' @param yp vector of responses (missing values are allowed)
# #' @param windows vector of bandwidths for the local polynomial smoother,
# #' one per explanatory variable.
# #' @param epsilon convergence criterion. Maximum allowed relative difference between
# #' consecutive estimates
# #' @param degree degree of the local polynomial smoother. Defaults to \code{0} (local constant).
# #' @param seed an integer used to set the seed of the pseudo-number generator that
# #' creates the \code{k} folds.
# #' @param max.it Maximum number of iterations for the algorithm.
# #' @param k.h tuning constant for a Huber-type loss function.
# #' @param k.t tuning constant for a Tukey-type loss function.
# #' @param type one of either \code{'Tukey'} or \code{'Huber'}.
# #'
# #' @return A real number with a robust measure of (hold-out) prediction error
# #'
# #' @author Matias Salibian-Barrera, \email{matias@stat.ubc.ca}, Alejandra Martinez
# #'
# #' @examples
# #' data(airquality)
# #' x <- airquality
# #' x <- x[complete.cases(x), c('Ozone', 'Solar.R', 'Wind', 'Temp')]
# #' y <- as.vector(x$Ozone)
# #' x <- as.matrix(x[, c('Solar.R', 'Wind', 'Temp')])
# #' backf.rob.cv(k=5, Xp = x, yp=y, windows=c(136.7, 8.9, 4.8), type='Tukey', degree=1)
# #'
# #' @export
# backf.rob.cv <- function(k=5, Xp, yp, windows, epsilon=1e-6, degree, type='Tukey', k.h=1.345,
#                          k.t = 4.685, seed=123, max.it=50) {
#   # does k-fold CV and returns "robust mean-squared prediction error"
#   n <- length(yp)
#   # k1 <- floor(n/k)
#   # ids <- rep(1:k, each=k1)
#   # if( length(ids) < n ) ids <- c(ids, 1:(n%%k))
#   # save existing random seed
#   if(exists(".Random.seed", where=.GlobalEnv)) old.seed <- .Random.seed
#   set.seed(seed)
#   ids <- sample( (1:n) %% k + 1 )
#   Xp <- as.matrix(Xp)
#   preds <- rep(NA, n)
#   for(j in 1:k) {
#     XX <- Xp[ids!=j, , drop=FALSE]
#     yy <- yp[ids!=j]
#     tmp <- try( backf.rob(yy ~ XX, point=Xp[ids==j,, drop=FALSE], windows=windows, epsilon=epsilon, degree=degree, type=type, max.it=max.it, k.h=k.h, k.t=k.t) )
#       #try( backf.rob(Xp=XX, yp=yy, point=Xp[ids==j,, drop=FALSE], windows=windows, epsilon=epsilon, degree=degree, type=type, max.it=max.it, k.h=k.h, k.t=k.t) )
#     if( class(tmp)[1] != 'try-error') {
#       preds[ids==j] <- rowSums(tmp$prediction) + tmp$alpha
#     }
#   }
#   # restore seed existing before call
#   if(exists('old.seed')) assign('.Random.seed', old.seed, envir=.GlobalEnv)
#   return( mad( (preds-yp), na.rm=TRUE )^2 + median( (preds-yp)^2, na.rm=TRUE ) )
# }
# 
# #' Cross-validation for the Classical Backfitting algorithm
# #'
# #' This function performs one run of K-fold cross-validation using the
# #' classical backfitting algorithm.
# #'
# #' This function performs one run of K-fold cross-validation using the
# #' classical backfitting algorithm and returns the mean squared
# #' hold out prediction error.
# #'
# #' @param k a positive integer indicating the number of folds.
# #' @param Xp a matrix (n x p) containing the explanatory variables
# #' @param yp vector of responses (missing values are allowed)
# #' @param windows vector of bandwidths for the local polynomial smoother,
# #' one per explanatory variable.
# #' @param epsilon convergence criterion. Maximum allowed relative difference between
# #' consecutive estimates
# #' @param degree degree of the local polynomial smoother. Defaults to \code{0} (local constant).
# #' @param seed an integer used to set the seed of the pseudo-number generator that
# #' creates the \code{k} folds.
# #' @param max.it Maximum number of iterations for the algorithm.
# #'
# #' @return A real number with the mean squared (hold-out) prediction error
# #'
# #' @author Matias Salibian-Barrera, \email{matias@stat.ubc.ca}, Alejandra Martinez
# #'
# #' @examples
# #' data(airquality)
# #' x <- airquality
# #' x <- x[complete.cases(x), c('Ozone', 'Solar.R', 'Wind', 'Temp')]
# #' y <- as.vector(x$Ozone)
# #' x <- as.matrix(x[, c('Solar.R', 'Wind', 'Temp')])
# #' backf.l2.cv(k=5, Xp = x, yp=y, windows=c(130, 9, 10), degree=1)
# #'
# #' @export
# backf.l2.cv <- function(k=5, Xp, yp, windows, epsilon=1e-6,
#                         degree, seed=123, max.it=50) {
#   # does k-fold CV and returns mean-squared prediction error
#   n <- length(yp)
#   # k1 <- floor(n/k)
#   # ids <- rep(1:k, each=k1)
#   # if( length(ids) < n ) ids <- c(ids, 1:(n%%k))
#   # save existing random seed
#   if(exists(".Random.seed", where=.GlobalEnv)) old.seed <- .Random.seed
#   set.seed(seed)
#   ids <- sample( (1:n) %% k + 1 )
#   Xp <- as.matrix(Xp)
#   preds <- rep(NA, n)
#   for(j in 1:k) {
#     XX <- Xp[ids!=j, , drop=FALSE]
#     yy <- yp[ids!=j]
#     tmp <- try( backf.cl(yy ~ XX, point=Xp[ids==j, , drop=FALSE], windows=windows, epsilon=epsilon, degree=degree, max.it=max.it) )
#       #try( backf.cl(Xp=XX, yp=yy, point=Xp[ids==j, , drop=FALSE], windows=windows, epsilon=epsilon, degree=degree, max.it=max.it) )
#     if( class(tmp)[1] != 'try-error') {
#       preds[ids==j] <- rowSums(tmp$prediction) + tmp$alpha
#     }
#   }
#   # restore seed existing before call
#   if(exists('old.seed')) assign('.Random.seed', old.seed, envir=.GlobalEnv)
#   return( mean( (preds-yp)^2, na.rm=TRUE ) )
# }


pos.est <- function(y, sigma.hat, typePhi, ini=NULL, epsilon=1e-6, iter.max=10){
  yp <- y[ tmp<-!is.na(y) ]
  if(is.null(ini)){
    ini <- median(yp)
  }
  corte <- 10
  iter <- 0
  n <- length(yp)
  prob <- rep(1,n)
  k.h <- 1.345
  k.t <- 4.685
  if(typePhi=='Huber'){
    beta <- .C("huber_pos", as.integer(n), as.double(yp), as.double(ini), as.double(epsilon), 
               as.double(sigma.hat), as.double(prob), as.double(k.h), as.integer(iter.max), salida=as.double(0) )$salida
  }
  if(typePhi=='Tukey'){
    beta <- .C("tukey_pos", as.integer(n), as.double(yp), as.double(ini), as.double(epsilon), 
               as.double(sigma.hat), as.double(prob), as.double(k.t), as.integer(iter.max), salida=as.double(0) )$salida
  }
  return(beta)
}



#' Residuals for objects of class \code{backf}
#'
#' This function returns the residuals of the fitted additive model using
#' the classical or robust backfitting estimators, as computed with \code{\link{backf.cl}} or
#' \code{\link{backf.rob}}.
#'
#' @param object an object of class \code{backf}, a result of a call to \code{\link{backf.cl}} or \code{\link{backf.rob}}.
#' @param ... additional other arguments. Currently ignored.
#'
#' @return A vector of residuals.
#'
#' @author Alejandra Mercedes Martinez \email{ale_m_martinez@hotmail.com}
#'
#' @export
residuals.backf <- function(object, ...){
  return( object$yp - rowSums(object$g.matrix) -object$alpha )
}

#' Fitted values for objects of class \code{backf}.
#'
#' This function returns the fitted values given the covariates of
#' the original sample under an additive model using the classical or
#' robust backfitting approach computed with \code{\link{backf.cl}} or
#' \code{\link{backf.rob}}.
#'
#' @param object an object of class \code{backf}, a result of a call to \code{\link{backf.cl}} or \code{\link{backf.rob}}.
#' @param ... additional other arguments. Currently ignored.
#'
#' @return A vector of fitted values.
#'
#' @author Alejandra Mercedes Martinez \email{ale_m_martinez@hotmail.com}
#'
#' @export
predict.backf <- function(object, ...){
  return( rowSums(object$g.matrix) + object$alpha )
}

#' Fitted values for objects of class \code{backf}
#'
#' This function returns the fitted values given the covariates of the original sample under an additive model using a classical or robust marginal integration procedure estimator computed with \code{backf.cl} or \code{backf.rob}.
#'
#' @param object an object of class \code{backf}, a result of a call to \code{\link{backf.cl}} or \code{\link{backf.rob}}.
#' @param ... additional other arguments. Currently ignored.
#'
#' @return A vector of fitted values.
#'
#' @author Alejandra Mercedes Martinez \email{ale_m_martinez@hotmail.com}
#'
#' @method fitted.values backf
#' @export
fitted.values.backf <- function(object,...){
  UseMethod("fitted")
}

#' @export
fitted.backf <- function(object,...){
  return(predict(object))
}


# plot.backf <- function(object, which=1:np, ask=FALSE,...){
#   Xp <- object$Xp
#   np <- dim(Xp)[2]
#   opar <- par(ask=ask)
#   on.exit(par(opar))
#   these <- rep(FALSE, np)
#   these[ which ] <- TRUE
#   if( np!= 1) {
#     for(i in 1:np) {
#       if(these[i]) {
#         ord <- order(Xp[,i])
#         x_name <- paste("x",i,sep="")
#         y_name <- bquote(paste(hat('g')[.(i)]))
#         if( !is.null(dim(object$g.matrix[,-i])) ){
#           res <- object$yp - rowSums(object$g.matrix[,-i])-object$alpha
#         } else {
#           res <- object$yp - object$g.matrix[,-i]-object$alpha
#         }
#         lim_cl <- c(min(res), max(res))
#         plot(Xp[ord,i],object$g.matrix[ord,i],type="l",lwd=3,main="",xlab=x_name,ylab=y_name, ylim=lim_cl)
#         points(Xp[,i], res, pch=20,col='gray45')
#       }
#     }
#   } else {
#     x_name <- "x"
#     y_name <- bquote(paste(hat('g')))
#     ord <- order(Xp)
#     res <- object$yp-object$alpha
#     lim_cl <- c(min(res), max(res))
#     plot( Xp[ord], object$g.matrix[ord], type='l', lwd=3, ylim=lim_cl, xlab=x_name, ylab=y_name)
#     points(Xp, res, pch=20, col='gray45')
#   }
# }


# plot.backf <- function(object, which=1:np, ask=FALSE, ...){
#   Xp <- object$Xp
#   np <- dim(Xp)[2]
#   opar <- par(ask=ask)
#   on.exit(par(opar))
#   these <- rep(FALSE, np)
#   these[ which ] <- TRUE
#     for(i in 1:np) {
#       if(these[i]) {
#         ord <- order(Xp[,i])
#         x_name <- paste("x",i,sep="")
#         y_name <- bquote(paste(hat('g')[.(i)]))
#         res <- object$yp - rowSums(object$g.matrix[,-i, drop=FALSE])-object$alpha
#         lim_cl <- c(min(res), max(res))
#         plot(Xp[ord,i], object$g.matrix[ord,i], type="l", lwd=3, main="", xlab=x_name, ylab=y_name, ylim=lim_cl)
#         points(Xp[,i], res, pch=20, col='gray45')
#       }
#     }
# }

#' Diagnostic plots for objects of class \code{backf}
#'
#' Plot method for objects of class \code{backf}.
#'
#' @param x an object of class \code{backf}, a result of a call to \code{\link{backf.cl}} or \code{\link{backf.rob}}.
#' @param which vector of indices of explanatory variables for which partial residuals plots will
#' be generaetd. Defaults to all available explanatory variables.
#' @param ask logical value. If \code{TRUE}, the graphical device will prompt for confirmation before
#' going to the next page/screen of output.
#' @param ... additional other arguments. Currently ignored.
#'
#' @author Alejandra Mercedes Martinez \email{ale_m_martinez@hotmail.com}
#'
#' @examples
#' tmp <- backf.rob(Ozone ~ Solar.R + Wind + Temp, data=airquality, 
#' subset=complete.cases(airquality), windows=c(136.7, 8.9, 4.8), degree=1)
#' plot(tmp, which=1:2)
#'
#' @export
plot.backf <- function(x, ask=FALSE, which=1:np, ...) {
  object <- x
  Xp <- object$Xp
  np <- dim(Xp)[2]
  opar <- par(ask=ask)
  on.exit(par(opar))
  these <- rep(FALSE, np)
  these[ which ] <- TRUE
  for(i in 1:np) {
    if(these[i]) {
      ord <- order(Xp[,i])
      if (is.null(colnames(Xp)) ){
        x_name <- bquote(paste('x')[.(i)])
      } else {
        x_name <- colnames(Xp)[i]
      }
      y_name <- bquote(paste(hat('g')[.(i)]))
      res <- object$yp - rowSums(object$g.matrix[,-i, drop=FALSE])-object$alpha
      lim_cl <- c(min(res), max(res))
      # plot(Xp[ord,i], object$g.matrix[ord,i], type="l", lwd=3, main="", xlab=x_name, ylab=y_name, ylim=lim_cl)
      # points(Xp[,i], res, pch=20, col='gray45')
      plot(Xp[,i], res, pch=20,col='gray45',main="",xlab=x_name,ylab=y_name, ylim=lim_cl,cex.lab=0.8)
      lines(Xp[ord,i],object$g.matrix[ord,i],lwd=3)
    }
  }
}


#Multiple R-squared 
R2 <- function(object,...){
  yp <- object$yp
  y <- yp[ tmp<-!is.na(yp) ]
  n <- length(y)
  res <- residuals(object)
  S02 <- sum((y-mean(y))^2)
  S2 <- sum(res^2)
  R2 <- (S02-S2)/S02
  return(R2)
}

#Robust multiple R-squared
R2.rob <- function(object,...){
  yp <- object$yp
  y <- yp[ tmp<-!is.na(yp) ]
  n <- length(y)
  S02 <- 0
  S2 <- 0
  res <- residuals(object)
  sigma.hat <- object$sigma.hat
  typePhi <- object$type
  pos <- pos.est(y, sigma.hat, typePhi=typePhi, ini=NULL, epsilon=1e-6, iter.max=50)
  if(typePhi=='Tukey'){
    for(i in 1:n){
      S02 <- S02 + rho.tukey((y[i]-pos)/sigma.hat)
      S2 <- S2 + rho.tukey(res[i]/sigma.hat)
    }
  }
  if(typePhi=='Huber'){
    for(i in 1:n){
      S02 <- S02 + rho.huber((y[i]-pos)/sigma.hat)
      S2 <- S2 + rho.huber(res[i]/sigma.hat)
    }
  }
  R2.rob <- (S02-S2)/S02
  return(R2.rob)
}


#' Summary for additive models fits using backfitting
#'
#' Summary method for class \code{backf}.
#'
#' This function returns the estimation of the intercept and also the
#' five-number summary and the mean of the residuals for both classical and
#' robust estimators. For the classical estimator, it also returns the R-squared.
#' For the robust estimator it returns a robust version of the R-squared and 
#' the estimate of the residual standard error.
#'
#' @param object an object of class \code{backf}, a result of a call to
#' \code{\link{backf.cl}} or \code{\link{backf.rob}}.
#' @param ... additional other arguments. Currently ignored.
#' 
#'
#' @author Alejandra Mercedes Martinez \email{ale_m_martinez@hotmail.com}
#'
#' @export
#' @aliases summary.backf summary.backf.cl summary.backf.rob
summary.backf <- function(object,...){
  NextMethod()
}

#' @export
summary.backf.cl <- function(object,...){
  #message("Estimate of the intercept: ", round(object$alpha,5))
  #message("Multiple R-squared: ", round(R2(object),5))
  #res <- residuals(object)
  #message("Residuals:")
  #summary(res)
  sal <- list(
    intercept=round(object$alpha,5),
    R2=round(R2(object),5),
    residuals=summary(residuals(object)),
    call = object$call
  )
  class(sal) <- c("summary.backf", "summary.backf.cl")
  return(sal)
}

#' @export
summary.backf.rob <- function(object,...){
  sal <- list(
    intercept=round(object$alpha,5),
    rse=round(object$sigma,5),
    R2=round(R2.rob(object),5),
    residuals=summary(residuals(object)),
    call = object$call
  )
  class(sal) <- c("summary.backf", "summary.backf.rob")
  return(sal)
}

#' @export
#' @aliases print.summary.backf print.summary.backf.cl print.summary.backf.rob
print.summary.backf <- function(x,...){
  NextMethod()
}

#' @export
print.summary.backf.cl <- function(x,...){
  cat("Estimate of the intercept: ", x$intercept, "\n")
  cat("Multiple R-squared: ", x$R2, "\n")
  #res <- residuals(object)
  cat("Residuals:\n", names(x$residuals), "\n", x$residuals, "\n")
  #summary(res)
}

#' @export
print.summary.backf.rob <- function(x,...){
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", 
                         collapse = "\n"), "\n\n", sep = "")
  cat("Estimate of the intercept: ", x$intercept, "\n")
  cat("Estimate of the residual standard error: ", x$rse, "\n")
  cat("Robust multiple R-squared: ", x$R2, "\n")
  #res <- residuals(object)
  cat("Residuals:\n", names(x$residuals), "\n", x$residuals, "\n")
  #summary(res)
}

#' Deviance for objects of class \code{backf}
#'
#' This function returns the deviance of the fitted additive model using one of the three
#' classical or robust marginal integration estimators, as computed with \code{\link{backf.cl}} or
#' \code{\link{backf.rob}}.
#'
#' @param object an object of class \code{backf}, a result of a call to \code{\link{backf.cl}} or \code{\link{backf.rob}}.
#' @param ... additional other arguments. Currently ignored.
#'
#' @return A real number.
#'
#' @author Alejandra Mercedes Martinez \email{ale_m_martinez@hotmail.com}
#'
#' @export
deviance.backf <- function(object, ...){
  NextMethod()
}

#' @export
deviance.backf.cl <- function(object, ...){
  return( sum( (residuals(object))^2) )
}

#' @export
deviance.backf.rob <- function(object, ...){
  yp <- object$yp
  y <- yp[ tmp<-!is.na(yp) ]
  n <- length(y)
  S2 <- 0
  res <- residuals(object)
  sigma.hat <- object$sigma.hat
  typePhi <- object$type
  if(typePhi=='Tukey'){
    for(i in 1:n){
      S2 <- S2 + rho.tukey(res[i]/sigma.hat)
    }
  }
  if(typePhi=='Huber'){
    for(i in 1:n){
      S2 <- S2 + rho.huber(res[i]/sigma.hat)
    }
  }
  return( S2 )
}



#' Additive model formula
#'
#' Description of the additive model formula extracted from an object of class \code{backf}.
#'
#' @param x an object of class \code{backf}, a result of a call to \code{\link{backf.cl}} or \code{\link{backf.rob}}.
#' @param ... additional other arguments. Currently ignored.
#'
#' @return A model formula.
#'
#' @author Alejandra Mercedes Martinez \email{ale_m_martinez@hotmail.com}
#'
#' @export
formula.backf <- function(x, ...){
  return(x$formula )
}


#' Print a Marginal Integration procedure
#'
#' The default print method for a \code{backf} object.
#'
#' @param x an object of class \code{backf}, a result of a call to \code{\link{backf.cl}} or \code{\link{backf.rob}}.
#' @param ... additional other arguments. Currently ignored.
#'
#' @return A real number.
#'
#' @author Alejandra Mercedes Martinez \email{ale_m_martinez@hotmail.com}
#'
#' @export
print.backf <- function(x, ...){
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", 
                         collapse = "\n"), "\n\n", sep = "")
  cat("\n")
  invisible(x)
}

