\name{VNTR}
\alias{VNTR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Genotyping Variable Number Tandem Repeats (VNTR) for the genome sequence of monkeypox virus (MPXV)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
VNTR(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{}
  \item{VNTR}{}
  \item{match_s}{}
  \item{mismatch_s}{}
  \item{regionStart}{}
  \item{regionEnd}{}
  \item{baseonly}{}
  \item{VNTRoutput}{}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{ID}{sequence name}
\item{r}{the number of tandem repeats}
\item{match}{the number of matches}
\item{mismatch}{the number of mismatches}
\item{indel}{the number of indels}
\item{score}{alignment score for VNTR region}
\item{start_pos}{start position of the VNTR region}
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
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
