\name{calcCoverage}
\alias{calcCoverage}
\title{Calculate coverages over a genomic region}
\usage{
    calcCoverage(input, mask, strand = NULL,
        ignore.strand = TRUE, rc = NULL)
}
\arguments{
    \item{input}{a \code{GRanges} object or a list of
    \code{GRanges} (not a \code{GRangesList}!) or the
    path to a BAM or BigWig file.}

    \item{mask}{a \code{GRanges} or \code{GRangesList}
    object.}
    
    \item{strand}{see the \code{strandedParams} in the 
    main \code{\link{recoup}} function.}
    
    \item{ignore.strand}{see the \code{strandedParams} in 
    the main \code{\link{recoup}} function.}
    
    \item{rc}{fraction (0-1) of cores to use in a multicore 
    system. It defaults to \code{NULL} (no parallelization).}
}
\value{
    A list of \code{Rle} objects representing the
    genomic coverages of interest.
}
\description{
    This function returns a coverage list for the genomic
    regions in \code{mask} argument. Generally it should not
    be used alone and is intended for internal use, although
    it is useful for calculating stand-alone coverages.
}
\details{
    \code{input} contains the short reads in one of the 
    formats described in the arguments section. When 
    \code{input} is a list, this list should contain one
    member per chromosome of the organism of interest.
    
    \code{mask} contains the genomic regions over which
    the coverage will be calculated from the input reads.
    When calculating RNA-Seq profiles, \code{mask} must be
    a named \code{GRangesList} where each member represents
    the exons of the respective gene.
}
\examples{
# Load some data
data("recoup_test_data",package="recoup")

# Calculate coverage Rle
mask <- makeGRangesFromDataFrame(df=test.genome,
    keep.extra.columns=TRUE)
small.cov <- calcCoverage(test.input[[1]]$ranges,mask)
}
\author{
    Panagiotis Moulos
}
