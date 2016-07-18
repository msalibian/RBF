\name{backf.cl}
\alias{backf.cl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Classic backfitting algorithm for additive models
%%  ~~function to do ... ~~
}
\description{The standard backfitting algorithm for additive models,
using a squared loss function and local polynomial smoothers.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
backf.cl(Xp, yp, point = NULL, windows, epsilon = 1e-06, degree = 0, prob = NULL, max.it = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Xp}{Matrix of explanatory variables (n by p).
%%     ~~Describe \code{Xp} here~~
}
  \item{yp}{Vector of responses (missing values are allowed).
%%     ~~Describe \code{yp} here~~
}
  \item{point}{Matrix of points where predictions will be computed and returned. 
%%     ~~Describe \code{point} here~~
}
  \item{windows}{Vector of bandwidths for the local polynomial smoother, one
per explanatory variable. 
%%     ~~Describe \code{windows} here~~
}
  \item{epsilon}{Convergence criterion.
%%     ~~Describe \code{epsilon} here~~
}
  \item{degree}{Degree of the local polynomial smoother. Defaults to ``0'' (local 
constant).
%%     ~~Describe \code{degree} here~~
}
  \item{prob}{Probabilities of observing each response (n). Defaults to ``NULL'' and
in that case it is ignored.
%%     ~~Describe \code{prob} here~~
}
  \item{max.it}{Maximum number of iterations for the algorithm.
%%     ~~Describe \code{max.it} here~~
}
}
\details{
The standard backfitting algorithm for additive models,
using a squared loss function and local polynomial smoothers.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{alpha }{Estimate for the intercept.}
\item{g.matrix }{Matrix of estimated additive components (n by p).}
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
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Additive Models }
