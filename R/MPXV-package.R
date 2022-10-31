

#' Genotyping Variable Number Tandem Repeats (VNTR) for the genome sequence of
#' monkeypox virus (MPXV)
#' 
#' The funciton \code{VNTR.Genotype} computes the copy of the variable number
#' tandem repeats.
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param data sequences from a file in FASTA format
#' @param VNTR variable number tandem repeat
#' @param match_s matching weight
#' @param mismatch_s mismatching penalty
#' @param regionStart start position of VNTR region
#' @param regionEnd end position of VNTR region
#' @param baseonly logical. If TRUE, only uses the letters in base alphabet
#' i.e. A,C,G,T.
#' @param VNTRoutput logical. If TRUE, export output to .csv file.
#' @param finder logical. If TRUE, call function \code{\link{STR_finder}}.
#' @return \item{ID}{name of sequence} \item{r}{the copy of tandem repeats}
#' \item{match}{the number of matches} \item{mismatch}{the number of
#' mismatches} \item{indel}{the number of indels} \item{score}{alignment score
#' of VNTR region for each query strain} \item{start_pos}{start position of the
#' VNTR region for each query strain}
#' @note %% ~~further notes~~
#' @author %% ~~who you are~~
#' @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#' @references %% ~put references to the literature/web site here ~
#' 
#' Pagès H, Aboyoun P, Gentleman R, DebRoy S (2022) \emph{Biostrings: Efficient
#' manipulation of biological strings}.  R package version 2.64.0,
#' \href{https://bioconductor.org/packages/Biostringshttps://bioconductor.org/packages/Biostrings}.
#' @examples
#' 
#' ## load example
#' data(example)
#' 
#' VNTR <- c("T","TATGATGGA","AT","ATATACATT")
#' regionStart <- c(132436,150542,173240,178413)
#' regionEnd <- c(133216,151501,173320,179244)
#' 
#' baseonly = T
#' match_s <- 2
#' mismatch_s <- -5
#' VNTRoutput = F
#' finder = F
#' 
#' 
#' out <- VNTR(data, VNTR=VNTR, match_s=match_s, mismatch_s=mismatch_s,
#'             regionStart=regionStart, regionEnd=regionEnd,baseonly = baseonly,VNTRoutput=VNTRoutput,finder=finder)
#' 
#' 
NULL



