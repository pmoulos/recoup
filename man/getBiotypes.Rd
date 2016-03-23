\name{getBiotypes}
\alias{getBiotypes}
\title{List default Ensembl biotypes}
\usage{
    getBiotypes(org)
}
\arguments{
    \item{org}{One of the supported \code{recoup} organisms
    See \code{\link{recoup}} for further information.}
}
\value{
    A character vector of biotypes.
}
\description{
    This function returns a character vector of Ensembl 
    biotypes for each supported organism. Mostly for 
    internal use, but can also be used to list the 
    biotypes and use some of them to subset initial 
    genomic regions to be profiled.
}
\examples{
    hg18.bt <- getBiotypes("hg18")
}
\author{
    Panagiotis Moulos
}

