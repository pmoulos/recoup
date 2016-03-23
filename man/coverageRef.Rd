\name{coverageRef}
\alias{coverageRef}
\title{Calculate coverage in a set of reference genomic 
    regions (ChIP-Seq or unspliced mode)}
\usage{
    coverageRef(
        input, 
        genomeRanges, 
        region = c( "tss", "tes", "genebody", "custom"),
        flank = c(2000, 2000),
        strandedParams = list(strand=NULL, ignoreStrand=TRUE),
        bamParams = NULL,
        rc = NULL
    )
}
\arguments{
    \item{input}{an input list as in \code{\link{recoup}} but
    with the \code{ranges} field of each member filled 
    (e.g. after using \code{\link{preprocessRanges}}).}

    \item{genomeRanges}{the \code{genome} from 
    \code{\link{recoup}} as a \code{GRanges} object (e.g. 
    the output from 
    \code{\link{makeGRangesFromDataFrame}}).}
    
    \item{region}{see the \code{region} argument in the 
    main \code{\link{recoup}} function.}
    
    \item{flank}{see the \code{flank} argument in the 
    main \code{\link{recoup}} function.}
    
    \item{strandedParams}{see the \code{strandedParams} 
    argument in the main \code{\link{recoup}} function.}
    
    \item{bamParams}{see the \code{bamParams} 
    argument in the main \code{\link{recoup}} function.}
    
    \item{rc}{fraction (0-1) of cores to use in a multicore 
    system. It defaults to \code{NULL} (no parallelization).}
}
\value{
    Same as input with the \code{ranges} fields filled.
}
\description{
    This function fills the \code{coverage} field in the
    main \code{input} argument in \code{\link{recoup}} 
    function.
}
\examples{
# Load some data
data("recoup_test_data",package="recoup")

# Calculate coverages
testGenomeRanges <- makeGRangesFromDataFrame(df=test.genome,
    keep.extra.columns=TRUE)
test.input <- coverageRef(
    test.input, 
    genomeRanges=testGenomeRanges, 
    region="tss",
    flank=c(2000,2000)
)
}
\author{
    Panagiotis Moulos
}

