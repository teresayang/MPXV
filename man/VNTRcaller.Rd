\name{VNTRcaller}
\alias{VNTRcaller}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calling Variable Number Tandem Repeats (VNTR) for the genome sequence of
monkeypox virus (MPXV)
%%  ~~function to do ... ~~
}
\description{The function \code{VNTRcaller} estimates the copy numbers of VNTRs.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
VNTRcaller(
  data,
  vntr = vntr,
  match_s = match_s,
  mismatch_s = mismatch_s,
  regionStart = regionStart,
  regionEnd = regionEnd,
  baseonly = T,
  fseqExtend = fseqExtend,
  VNTRoutput = F,
  tracker = F
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{MPXV sequences from a file in FASTA format}

\item{vntr}{sequence of the repetitive unit}

\item{match_s}{matching weight}

\item{mismatch_s}{mismatching penalty}

\item{regionStart}{starting position of a VNTR region}

\item{regionEnd}{ending position of a VNTR region}

\item{baseonly}{logical. If TRUE, only the letters of the nucleotide bases
(i.e. A,C,G,T) are considered. If FALSE, the degenerate codes are also
considered.}

\item{fseqExtend}{the length of extended flanking sequence}

\item{VNTRoutput}{logical. If TRUE, the output is written to .csv files.}

\item{tracker}{logical. If TRUE, call function \code{\link{VNTRtracker}}.}
}

\value{
\item{ID}{name of MPXV sequences} \item{r}{copy of tandem repeats}
\item{match}{number of matches} \item{mismatch}{number of mismatches}
\item{indel}{number of insertions and deletions (indels)}
\item{score}{alignment score of a VNTR region for a query strain}
\item{start_pos}{starting position of the VNTR region for a query strain}
\item{count_n}{number of N characters}
\item{count_deg}{number of non-ATCGN characters}
\item{count_len}{VNTR region length}
\item{vntr_match_norm}{proportion of number of N characters to VNTR region length}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
Yang, H.-C. et al (2022) Monkeypox genome contains variable number
tandem repeats enabling accurate virus tracking.

Pages H, Aboyoun P, Gentleman R, DebRoy S (2022) \emph{Biostrings: Efficient
manipulation of biological strings}.  R package version 2.64.0.
%%\href{https://bioconductor.org/packages/Biostringshttps://bioconductor.org/packages/Biostrings}.
}

\examples{
## Read MPXV example sequences from a file in FASTA format.
MPXVseq <- read.fasta(system.file("extdata/MPXV_seq_example.fasta.gz", package = "MPXV"), as.string = T)

## VNTR
vntr <- c("T","TATGATGGA","AT","ATATACATT")
regionStart <- c(132436,150542,173240,178413)
regionEnd <- c(133216,151501,173320,179244)

## parameter settings
match_s = 2
mismatch_s = -5
baseonly = TRUE
fseqExtend = 150
VNTRoutput = FALSE
tracker = FALSE

## computes the copy of the variable number tandem repeats
out <- VNTRcaller(data = MPXVseq, vntr = vntr,
                  regionStart = regionStart, regionEnd = regionEnd,
                  match_s = match_s, mismatch_s=mismatch_s,
                  baseonly = baseonly, fseqExtend = fseqExtend,
                  VNTRoutput = VNTRoutput,
                  tracker = tracker)

## For more details about the output of VNTRcaller, please vitsit https://github.com/teresayang/MPXV_VNTR.

}


