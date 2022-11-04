\name{VNTRtracker}
\alias{VNTRtracker}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Tracking the origin of monkeypox virus (MPXV) by using Variable Number Tandem Repeats (VNTR) markers.
}
\description{
The function \code{VNTRtracker} calculates the MPXV strain-strain genetic distance to track the origin of a given MPXV by using VNTRs.
}
\usage{
VNTRtracker(r, 
	out_dir="/", 
	file_name="Tracking_output")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r}{A vector of copies of VNTR of MPXV.}
  \item{out_dir}{Directory name of the output files.}
  \item{file_name}{Filename of the output.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
The four VNTR markers nt133095[A25R−A26L] (T), nt150560[A47R−A49R] (TATGATGGA), nt173267[B14R−B15L] (AT), and nt179074[B18R−B19R] (ATATACATT) 
are proposed to track the original of MPXV strains. Based on these four VNTRs, we calculate three distances (PIC, entropy, and length of repetitive unit (see ref.)) between a query MPXV strain and the strains in our database.  
}
\value{
The output contains multiple columns, including Accession ID, Country name, lineage, Release Date, clade, lineage, Collection Date, and three distances mentioned above.

}
\references{
%% ~put references to the literature/web site here ~
Yang, H.-C. et al (2022)
\emph{Monkeypox genome contains variable number tandem repeats enabling accurate virus tracking}.
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
r=c(27,7,24,16) 

VNTRtracker(r, 
	out_dir="C:/", 
	file_name="Tracking_output")
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