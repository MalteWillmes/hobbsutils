\name{colmatch}
\alias{colmatch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
colmatch(x, y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
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
function (x, y) 
{
    if (all(colnames(x) \%in\% colnames(y)) & all(colnames(y) \%in\% 
        colnames(x))) {
        return("All column names match between both objects")
    }
    else if (!all(colnames(x) \%in\% colnames(y)) & all(colnames(y) \%in\% 
        colnames(x))) {
        if (length(colnames(x)[!colnames(x) \%in\% colnames(y)]) > 
            1) {
            cat("\n\n", "columns ")
            cat(colnames(x)[!colnames(x) \%in\% colnames(y)], sep = " and ")
            cat(" are missing from second object\n")
        }
        else {
            cat("\n\n", "column ")
            cat(colnames(x)[!colnames(x) \%in\% colnames(y)], sep = " and ")
            cat(" is missing from second object\n")
        }
    }
    else if (all(colnames(x) \%in\% colnames(y)) & !all(colnames(y) \%in\% 
        colnames(x))) {
        if (length(colnames(y)[!colnames(y) \%in\% colnames(x)]) > 
            1) {
            cat("\n\n", "columns ")
            cat(colnames(y)[!colnames(y) \%in\% colnames(x)], sep = " and ")
            cat(" are missing from first object\n")
        }
        else {
            cat("\n\n", "column ")
            cat(colnames(y)[!colnames(y) \%in\% colnames(x)], sep = " and ")
            cat(" is missing from first object\n")
        }
    }
    else if (!all(colnames(x) \%in\% colnames(y)) & !all(colnames(y) \%in\% 
        colnames(x))) {
        if (length(colnames(y)[!colnames(y) \%in\% colnames(x)]) > 
            1 & length(colnames(x)[!colnames(x) \%in\% colnames(y)]) > 
            1) {
            cat("\n\n", "columns ")
            cat(colnames(y)[!colnames(y) \%in\% colnames(x)], sep = " and ")
            cat(" are missing from first object \n and")
            cat("\n", "columns ")
            cat(colnames(x)[!colnames(x) \%in\% colnames(y)], sep = " and ")
            cat(" are missing from second object\n")
        }
        else if (length(colnames(y)[!colnames(y) \%in\% colnames(x)]) > 
            1 & length(colnames(x)[!colnames(x) \%in\% colnames(y)]) == 
            1) {
            cat("\n\n", "columns ")
            cat(colnames(y)[!colnames(y) \%in\% colnames(x)], sep = " and ")
            cat(" are missing from first object \n and")
            cat("\n", "column ")
            cat(colnames(x)[!colnames(x) \%in\% colnames(y)], sep = " and ")
            cat(" is missing from second object\n")
        }
        else if (length(colnames(y)[!colnames(y) \%in\% colnames(x)]) == 
            1 & length(colnames(x)[!colnames(x) \%in\% colnames(y)]) > 
            1) {
            cat("\n\n", "column ")
            cat(colnames(y)[!colnames(y) \%in\% colnames(x)], sep = " and ")
            cat(" is missing from first object \n and")
            cat("\n", "columns ")
            cat(colnames(x)[!colnames(x) \%in\% colnames(y)], sep = " and ")
            cat(" are missing from second object\n")
        }
        else {
            cat("\n\n", "column ")
            cat(colnames(y)[!colnames(y) \%in\% colnames(x)], sep = " and ")
            cat(" is missing from first object \n and")
            cat("\n", "column ")
            cat(colnames(x)[!colnames(x) \%in\% colnames(y)], sep = " and ")
            cat(" is missing from second object\n")
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
