#' @useDynLib RBF, .registration = TRUE
# Tukey's Psi
#' @export
psi.tukey <- function(r, k=4.685){
  u <- abs(r/k)
  w <- r*((1-u)*(1+u))^2
  w[u>1] <- 0
  return(w)
}

#Tukey's weight function "Psi(r)/r"
#' @export
psi.w <- function(r, k= 4.685){
  u <- abs(r/k)
  w <- ((1 + u) * (1 - u))^2
  w[u > 1] <- 0
  return(w)
}

# Huber's Psi
#' @export
psi.huber <- function(r, k=1.345)
  pmin(k, pmax(-k, r))

#Huber's weight function "Psi(r)/r"
psi.huber.w <- function(r, k=1.345)
  pmin(1, k/abs(r))

#Nucleo de Epanechnikov
#' @export
k.epan<-function(x) {
  a <- 0.75*(1-x^2)
  tmp <- a*(abs(x)<=1)
  return(tmp)
}

#Norma 2 Euclidea
my.norm.2 <- function(x) sqrt(sum(x^2))


## Classic Backfitting
#' @export
backf.cl <- function(Xp, yp, point=NULL, windows, epsilon=1e-6, degree=0,
                     prob=NULL, max.it=100) {
  # Xp = covariance matrix (n x q)
  # yp = response vector (NA's are allowed).
  # point = vector of length q where prediction is computed.
  #      If missing, predictions are computed for each row of Xp.
  # windows = vector of length q with kernel windows
  # epsilon = convergence criterion
  # degree = degree of the local polynomial smoother
  #        A value of 0 means locally constant fits.
  # prob = probabilities of observing each response (n)
  # max.it = max number of iterations

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
        we <- k.epan( (Xp[,j] - Xp[i,j]) / windows[j] )
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
    if(is.null(dim(point))) {
      prediccion <- mpunto <- t(as.matrix(point))
    } else {
      prediccion <- mpunto <- point
    }
    np <- dim(mpunto)[1]
    for(k in 1:np){
      for(j in 1:q){
        y.tilde.bis <- yp - alpha - rowSums(g.matriz[,-j,drop=FALSE]) #- aux1
        we <- k.epan( (Xp[,j] - mpunto[k,j]) / windows[j] )
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
  object <- list(alpha=alpha, g.matrix=g.matriz, prediction=prediccion, Xp=Xp, yp=yp)
  class(object) <- c("backf.cl", "backf", "list")
  return(object)
  # return(list(alpha=alpha, g.matrix=g.matriz, prediction=prediccion))
}


# Robust Backfitting
#' @export
backf.rob <- function(Xp, yp, windows, point=NULL, epsilon=1e-6, degree=0,
                      sigma.hat=NULL, prob=NULL, max.it=20, k.h=1.345,
                      k.t = 4.685, type='Huber'){
  # Xp = covariance matrix (n x q)
  # yp = response vector (NA's are allowed).
  # point = vector of length q where prediction is computed.
  #      If missing, predictions are computed for each row of Xp.
  # windows = vector of length q with kernel windows
  # epsilon = convergence criterion
  # degree = degree of the local polynomial smoother
  #        A value of 0 means locally constant fits.
  # prob = probabilities of observing each response (n)
  # max.it = max number of iterations
  # sigma.hat = estimate of the residual standard error. If missing we use the
  # mad of the residuals obtained with local medians.
  # k.h = tuning constant for the Huber function
  # k.t = tuning constant for the Tukey function
  # type = 'Huber' or 'Tukey'

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
                                  as.double(prob), as.double(k.h), as.integer(max.it),salida=as.double(0) )$salida
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
                                  as.double(prob), as.double(k.t), as.integer(max.it),salida=as.double(0) )$salida
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
                                salida=as.double(rep(0, degree+1)) )$salida[1]
          }
          if( type == 'Tukey') {
            beta.ini <- .C("kernel_huber_lin", as.double(Xp[i,j]), as.double(Xp[,j]), as.integer(n.miss),
                           as.double(y.tilde.bis), as.double(tmp), as.integer(degree),
                           as.double(beta.ini), as.double(windows[j]), as.double(epsilon), as.double(sigma.hat),
                           as.double(prob), as.double(k.h), as.integer(max.it),
                           salida=as.double(rep(0, degree+1)) )$salida
            if(!any(is.na(beta.ini))) {
              g.matriz[i,j] <- .C("kernel_tukey_lin", as.double(Xp[i,j]), as.double(Xp[,j]), as.integer(n.miss),
                                  as.double(y.tilde.bis), as.double(tmp), as.integer(degree),
                                  as.double(beta.ini), as.double(windows[j]), as.double(epsilon), as.double(sigma.hat),
                                  as.double(prob), as.double(k.t), as.integer(max.it),
                                  salida=as.double(rep(0, degree+1)) )$salida[1]
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
                    salida=as.double(0) )$salida
      }
      if( type=='Tukey') {
        mu.ini <- median( y.tilde.bis )
        mu.ini <- .C("huber_pos", as.integer(n.miss),
                     as.double(y.tilde.bis), as.double(mu.ini),
                     as.double(epsilon), as.double(sigma.hat),
                     as.double(prob), as.double(k.h), as.integer(max.it),
                     salida=as.double(0) )$salida
        alpha <- .C("tukey_pos", as.integer(n.miss),
                    as.double(y.tilde.bis), as.double(mu.ini),
                    as.double(epsilon), as.double(sigma.hat),
                    as.double(prob), as.double(k.t), as.integer(max.it),
                    salida=as.double(0) )$salida
      }
    } else { stop('Error computing additive components'); break()
    }
  }

  prediccion <- NULL

  if(!is.null(point)){
    if(is.null(dim(point))) {
      prediccion <- mpunto <- t(as.matrix(point))
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
                                    as.double(prob), as.double(k.h), as.integer(max.it),salida=as.double(0) )$salida
            }
            if( type=='Tukey') {
              prediccion[k,j] <- .C("kernel_tukey_pos", as.double(mpunto[k,j]), as.double(as.matrix(Xp[,j])), as.integer(n.miss),
                                    as.double(y.tilde.bis), as.double(mu.ini), as.double(windows[j]),
                                    as.double(epsilon), as.double(sigma.hat),
                                    as.double(prob), as.double(k.t), as.integer(max.it),salida=as.double(0) )$salida
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
                                  salida=as.double(rep(0, degree+1)) )$salida[1]
          }
          if( type == 'Tukey') {
            beta.ini <- .C("kernel_huber_lin", as.double(mpunto[k,j]), as.double(Xp[,j]), as.integer(n.miss),
                           as.double(y.tilde.bis), as.double(tmp), as.integer(degree),
                           as.double(beta.ini), as.double(windows[j]), as.double(epsilon), as.double(sigma.hat),
                           as.double(prob), as.double(k.h), as.integer(max.it),
                           salida=as.double(rep(0, degree+1)) )$salida
            if(!any(is.na(beta.ini))) {
              prediccion[k,j] <- .C("kernel_tukey_lin", as.double(mpunto[k,j]), as.double(Xp[,j]), as.integer(n.miss),
                                    as.double(y.tilde.bis), as.double(tmp), as.integer(degree),
                                    as.double(beta.ini), as.double(windows[j]), as.double(epsilon), as.double(sigma.hat),
                                    as.double(prob), as.double(k.t), as.integer(max.it),
                                    salida=as.double(rep(0, degree+1)) )$salida[1]
            } else {
              prediccion[k,j] <- NA
            }
          }
        }
      }
    }
  }
  object <- list(alpha=alpha,g.matrix=g.matriz, sigma.hat=sigma.hat, prediction=prediccion, Xp=Xp, yp=yp)
  class(object) <- c("backf.rob", "backf", "list")
  return(object)
  # return(list(alpha=alpha,g.matrix=g.matriz, sigma.hat=sigma.hat, prediction=prediccion))
}


# CV functions
#' @export
backf.rob.cv <- function(k=5, Xp, yp, windows, epsilon,
                         degree, type, seed=123, max.it=50) {
  # does k-fold CV and returns "robust mean-squared prediction error"
  n <- length(yp)
  # k1 <- floor(n/k)
  # ids <- rep(1:k, each=k1)
  # if( length(ids) < n ) ids <- c(ids, 1:(n%%k))
  # save existing random seed
  if(exists(".Random.seed", where=.GlobalEnv)) old.seed <- .Random.seed
  set.seed(seed)
  ids <- sample( (1:n) %% k + 1 )
  Xp <- as.matrix(Xp)
  ids <- sample(ids)
  preds <- rep(NA, n)
  for(j in 1:k) {
    XX <- Xp[ids!=j, , drop=FALSE]
    yy <- yp[ids!=j]
    tmp <- try( backf.rob(Xp=XX, yp=yy, point=Xp[ids==j,], windows=windows, epsilon=epsilon,
                          degree=degree, type=type, max.it=max.it) )
    if( class(tmp)[1] != 'try-error') {
      preds[ids==j] <- rowSums(tmp$prediction) + tmp$alpha
    }
  }
  # restore seed existing before call
  if(exists('old.seed')) assign('.Random.seed', old.seed, envir=.GlobalEnv)
  return( mad( (preds-yp), na.rm=TRUE )^2 + median( (preds-yp), na.rm=TRUE )^2 )
}

#' @export
backf.l2.cv <- function(k=5, Xp, yp, windows, epsilon,
                        degree, seed=123, max.it=50) {
  # does k-fold CV and returns mean-squared prediction error
  n <- length(yp)
  # k1 <- floor(n/k)
  # ids <- rep(1:k, each=k1)
  # if( length(ids) < n ) ids <- c(ids, 1:(n%%k))
  # save existing random seed
  if(exists(".Random.seed", where=.GlobalEnv)) old.seed <- .Random.seed
  set.seed(seed)
  ids <- sample( (1:n) %% k + 1 )
  Xp <- as.matrix(Xp)
  ids <- sample(ids)
  preds <- rep(NA, n)
  for(j in 1:k) {
    XX <- Xp[ids!=j, , drop=FALSE]
    yy <- yp[ids!=j]
    tmp <- try( backf.cl(Xp=XX, yp=yy, point=Xp[ids==j,], windows=windows, epsilon=epsilon,
                         degree=degree, max.it=max.it) )
    if( class(tmp)[1] != 'try-error') {
      preds[ids==j] <- rowSums(tmp$prediction) + tmp$alpha
    }
  }
  # restore seed existing before call
  if(exists('old.seed')) assign('.Random.seed', old.seed, envir=.GlobalEnv)
  return( mean( (preds-yp)^2, na.rm=TRUE ) )
}


#' Residuals for objects of class \code{backf}
#'
#' This function returns the residuals of the fitted additive model using
#' the classical or robust backfitting estimators, as computed with \link{\code{backf.cl}} or
#' \link{\code{backf.rob}}.
#'
#' @param object an object of class \code{backf}, a result of a call to \code{\link{backf.cl}} or \code{\link{backf.rob}}.
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
#' robust backfitting approach computed with \link{\code{backf.cl}} or
#' \link{\code{backf.rob}}.
#'
#' @param object an object of class \code{backf}, a result of a call to \code{\link{backf.cl}} or \code{\link{backf.rob}}.
#'
#' @return A vector of fitted values.
#'
#' @author Alejandra Mercedes Martinez \email{ale_m_martinez@hotmail.com}
#'
#' @export
predict.backf <- function(object, ...){
  return( rowSums(object$g.matrix) + object$alpha )
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
#' @param object an object of class \code{backf}, a result of a call to \code{\link{backf.cl}} or \code{\link{backf.rob}}.
#' @param which vector of indices of explanatory variables for which partial residuals plots will
#' be generaetd. Defaults to all available explanatory variables.
#' @param ask logical value. If \code{TRUE}, the graphical device will prompt for confirmation before
#' going to the next page/screen of output.
#' @param ... additional other arguments. Currently ignored.
#'
#' @author Alejandra Mercedes Martinez \email{ale_m_martinez@hotmail.com}
#'
#' @export
plot.backf <- function(object, which=1:np, ask=FALSE, ...){
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
      plot(Xp[ord,i], object$g.matrix[ord,i], type="l", lwd=3, main="", xlab=x_name, ylab=y_name, ylim=lim_cl)
      points(Xp[,i], res, pch=20, col='gray45')
    }
  }
}

#' Summary for additive models fits using backfitting
#'
#' Summary method for class \code{backf}.
#'
#' This function returns the estimation of the intercept and also the
#' five-number summary and the mean of the residuals for both classical and
#' robust estimators. For the robust estimator it also returns the estimate of
#' the residual standard error.
#'
#' @param object an object of class \code{backf}, a result of a call to
#' \code{\link{backf.cl}} or \code{\link{backf.rob}}.
#' @param ... additional other arguments. Currently ignored.
#'
#' @author Alejandra Mercedes Martinez \email{ale_m_martinez@hotmail.com}
#'
#' @export
#' @aliases summary.backf summary.backf.cl summary.backf.rob
summary.backf <- function(object,...){
  NextMethod()
}

summary.backf.cl <- function(object,...){
  message("Estimate of the intercept: ", round(object$alpha,5))
  res <- residuals(object)
  message("Residuals:")
  summary(res)
}

summary.backf.rob <- function(object,...){
  message("Estimate of the intercept: ", round(object$alpha,5))
  message("Estimate of the residual standard error: ", round(object$sigma,5))
  res <- residuals(object)
  message("Residuals:")
  summary(res)
}

