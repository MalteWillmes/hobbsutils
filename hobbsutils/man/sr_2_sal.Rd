\name{sr_2_sal}
\alias{sr_2_sal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Convert measured 87Sr/86Sr isotope ratio to salinity

}
\description{
This function calculates salinity from measured Sr using a 2 end member mixing model. It requires the measured Sr ratio, Endmember 1 Strontium, Endmember 2 Strontium, concentrations for each, and salinity of the two endmembers. Default values only correct for Sacramento Delta
}
\usage{
sr_2_sal(sr, srfw = 0.705781, srmar = 0.70918, confw = 103.8, conmar = 7900)

}

\arguments{
  \item{sr}{Measured strontium. Should be giving as the ratio of Sr87/Sr86.}
  \item{srfw}{Strontium ratio of freshwater (or endmember 1).}
  \item{srmar}{Strontium ratio of the ocean (or endmember 2).}
  \item{confw}{Concentation (parts per million) of Strontium in freshwater (or endmember 1).}
  \item{conmar}{Concentration (parts per million) of Strontium in the ocean (or endmember 2).}
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

{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
