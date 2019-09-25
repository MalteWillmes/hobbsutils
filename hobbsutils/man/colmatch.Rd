\name{colnamefun}
\alias{colnamefun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Check if any column names are missing between two objects
}
\description{
This function will examine the column names of two provided objects and will return the names of any missing columns and which object is missing which columns
}
\usage{
colnamefun(x,y, match = F, join = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Name of first object with column names (generally a dataframe).
}
  \item{y}{Name of second object with column names (gnerally a dataframe).
}
  \item{match}{expects logical value either T/TRUE or F/FALSE
  }
  \item{join} {expects logical value either T/TRUE or F/FALSE
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{match }{Wether to return columns that match or return columns that are missing. Match = F/FALSE is default and will return missing columns. Match == T/TRUE will return matching columns.}
  \item{join }{  Wether or not to format return values for use in joining. Default is join = F/False. This will return the results formatted for legibility in the console and WILL NOT WORK FOR JOINING. Join = T/True will format the results more simply so that it is legible by the "by" argument in functions like left_join(). }
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
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
