\name{VNTRtracker}
\alias{VNTRtracker}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Tracking the origin of monkeypox virus (MPXV) by using Variable Number Tandem Repeats (VNTR) markers
}
\description{
The function \code{VNTRtracker} calculates the MPXV strain-strain genetic distance to track the origin of a given MPXV by using VNTRs.
}
\usage{
VNTRtracker(r, file_name = "Tracking_output")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r}{A vector of copies of VNTR of MPXV.}
  \item{file_name}{Filename of the output.}
}
\details{
 The four VNTR markers nt133095[A25R - A26L] (T), nt150560[A47R - A49R] (TATGATGGA), nt173267[B14R - B15L] (AT), and nt179074[B18R - B19R] (ATATACATT) are proposed to track the original of MPXV strains. Based on these four VNTRs, we calculate three distances (ploymorphism information content, entropy, and length of repetitive unit (see Yang et al. (2022)) between a query MPXV strain and the strains in our database.
}
\value{
The output contains multiple columns, including accession ID, country name, collection date, clade, lineage, and three distances mentioned above.
}
\references{
Yang, H.-C. et al (2022)
  Monkeypox genome contains variable number tandem repeats enabling accurate virus tracking.
}

\examples{
r=c(27,7,24,16)

## The output is a csv file with the file name given by the parameter "file_name".
## In this file, each row represents a strain in our database, and column variables contain accession ID, country name, collection date, clade, lineage, and three distances (see Yang et al. (2022)) between a query MPXV strain and a strains in our database. 
## The rank of strains are ordered by ploymorphism information content distance (see Yang et al. (2022)) from low to high. 
VNTRtracker(r, file_name="Tracking_output")

}

