\name{sr_2_sal}
\alias{sr_2_sal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sr_2_sal(sr, srfw = 0.705264, srmar = 0.70918, confw = 74.6, conmar = 6819, salfw = 0.1, salmar = 31.8, sallim, fill = "NA")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sr}{
%%     ~~Describe \code{sr} here~~
}
  \item{srfw}{
%%     ~~Describe \code{srfw} here~~
}
  \item{srmar}{
%%     ~~Describe \code{srmar} here~~
}
  \item{confw}{
%%     ~~Describe \code{confw} here~~
}
  \item{conmar}{
%%     ~~Describe \code{conmar} here~~
}
  \item{salfw}{
%%     ~~Describe \code{salfw} here~~
}
  \item{salmar}{
%%     ~~Describe \code{salmar} here~~
}
  \item{sallim}{
%%     ~~Describe \code{sallim} here~~
}
  \item{fill}{
%%     ~~Describe \code{fill} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
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

## The function is currently defined as
function (sr, srfw = 0.705264, srmar = 0.70918, confw = 74.6, 
    conmar = 6819, salfw = 0.1, salmar = 31.8, sallim, fill = "NA") 
{
    if (any(sr < min(srfw, srmar) | sr > max(srfw, srmar), na.rm = T)) {
        warning("Some of your measured strontium ratio values are outside the bounds of your two endmembers, make sure that srfw and srmar are set correctly", 
            call. = F, immediate. = T)
    }
    if (is.na(sallim)) 
        stop("You have not set a high salinity limit (sallim argument)")
    sal = (((salfw * srmar * conmar) - (salfw * sr * conmar) - 
        (salmar * srmar * conmar) + (salmar * sr * conmar))/((sr * 
        confw) - (sr * conmar) - (srfw * confw) + (srmar * conmar))) + 
        salmar
    if (fill == "NA") {
        sal[sal > sallim] = sallim
    }
    else if (fill == "sallim") {
        sal[sal > sallim] = "NA"
    }
    return(sal)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
