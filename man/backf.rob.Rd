\name{backf.rob}
\alias{backf.rob}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A robust backfitting algorithm for additive models
%%  ~~function to do ... ~~
}
\description{A robust backfitting algorithm for additive models.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
backf.rob(Xp, yp, windows, point = NULL, epsilon = 1e-06, degree = 0, sigma.hat = NULL, prob = NULL, max.it = 20, k.h = 1.345, k.t = 4.685, type = "Huber")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Xp}{Matrix of explanatory variables (n by p).
%%     ~~Describe \code{Xp} here~~
}
  \item{yp}{Vector of responses (missing values are allowed).
%%     ~~Describe \code{yp} here~~
}
  \item{windows}{Vector of bandwidths for the local polynomial smoother, one
per explanatory variable.
%%     ~~Describe \code{windows} here~~
}
  \item{point}{Matrix of points where predictions will be computed and returned. 
%%     ~~Describe \code{point} here~~
}
  \item{epsilon}{Convergence criterion.
%%     ~~Describe \code{epsilon} here~~
}
  \item{degree}{Degree of the local polynomial smoother. Defaults to ``0'' (local 
constant).
%%     ~~Describe \code{degree} here~~
}
  \item{sigma.hat}{Estimate of the residual standard error. If NULL we use the mad of the residuals obtained with local medians.
%%     ~~Describe \code{sigma.hat} here~~
}
  \item{prob}{Probabilities of observing each response (n). Defaults to ``NULL'' and
in that case it is ignored.
%%     ~~Describe \code{prob} here~~
}
  \item{max.it}{Maximum number of iterations for the algorithm.
%%     ~~Describe \code{max.it} here~~
}
  \item{k.h}{Tuning constant for a Huber-type loss function. 
%%     ~~Describe \code{k.h} here~~
}
  \item{k.t}{Tuning constant for a Tukey-type loss function. 
%%     ~~Describe \code{k.t} here~~
}
  \item{type}{One of either ``Tukey'' or ``Huber''.
%%     ~~Describe \code{type} here~~
}
}
\details{A robust backfitting algorithm for additive models,
using local polynomial smoothers.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{alpha }{Estimate for the intercept.}
\item{g.matrix }{Matrix of estimated additive components (n by p).}
\item{sigma.hat }{Estimate of the residual standard error.}
\item{prediction }{Matrix of estimated additive components 
for the points listed in the argument point.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{Matias Salibian-Barrera, Alejandra Martinez
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\keyword{ Robust Backfitting }
