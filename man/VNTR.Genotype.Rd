\name{VNTR.Genotype}
\alias{VNTR.Genotype}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Genotyping Variable Number Tandem Repeats (VNTR) for the genome sequence of monkeypox virus (MPXV)
}
\description{
The funciton \code{VNTR.Genotype} computes the copy of the variable number tandem repeats.
}
\usage{
VNTR.Genotype(data, vntr,
      regionStart, regionEnd,
      match_s = match_s,
      mismatch_s = mismatch_s,
      baseonly = True,
      VNTRoutput = False,
      finder = False)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{sequences from a file in FASTA format}
  \item{vntr}{variable number tandem repeat}
  \item{match_s}{matching weight}
  \item{mismatch_s}{mismatching penalty}
  \item{regionStart}{start position of VNTR region}
  \item{regionEnd}{end position of VNTR region}
  \item{baseonly}{logical. If TRUE, only uses the letters in base alphabet i.e. A,C,G,T.}
  \item{VNTRoutput}{logical. If TRUE, export output to .csv file.}
  \item{finder}{logical. If TRUE, call function \code{\link{STR_finder}}. }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{ID}{name of sequence}
\item{r}{the copy of tandem repeats}
\item{match}{the number of matches}
\item{mismatch}{the number of mismatches}
\item{indel}{the number of indels}
\item{score}{alignment score of VNTR region for each query strain}
\item{start_pos}{start position of the VNTR region for each query strain}
}
\references{
%% ~put references to the literature/web site here ~

  Pagès H, Aboyoun P, Gentleman R, DebRoy S (2022)
  \emph{Biostrings: Efficient manipulation of biological strings}.
  R package version 2.64.0, \href{https://bioconductor.org/packages/Biostrings}{https://bioconductor.org/packages/Biostrings}.

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
## load example
data(example)


## VNTR
vntr <- c("T","TATGATGGA","AT","ATATACATT")
regionStart <- c(132436,150542,173240,178413)
regionEnd <- c(133216,151501,173320,179244)

## parameter settings
baseonly = T
match_s = 2
mismatch_s = -5
VNTRoutput = F
finder = F

## computes the copy of the variable number tandem repeats
out <- VNTR.Genotype(data = MPXVseq, vntr = vntr,
                     regionStart = regionStart, regionEnd = regionEnd,
                     match_s = match_s, mismatch_s=mismatch_s,
                     baseonly = baseonly,VNTRoutput = VNTRoutput,
                     finder = finder)

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
