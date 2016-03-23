\name{profileMatrix}
\alias{profileMatrix}
\title{Calculate final profile matrices for plotting}
\usage{
    profileMatrix(input, flank, binParams, rc = NULL)
}
\arguments{
    \item{input}{an input list as in \code{\link{recoup}} but
    with the \code{ranges} the \code{coverage} fields of each 
    member filled (e.g. after using 
    \code{\link{preprocessRanges}} and 
    \code{\link{coverageRef}}).}

    \item{flank}{see the \code{flank} argument in the 
    main \code{\link{recoup}} function.}
    
    \item{binParams}{see the \code{binParams} argument in 
    the main \code{\link{recoup}} function.}

    \item{rc}{fraction (0-1) of cores to use in a multicore 
    system. It defaults to \code{NULL} (no parallelization).}
}
\value{
    Same as input with the \code{profile} fields filled.
}
\description{
    This function fills the \code{profile} field in the
    main \code{input} argument in \code{\link{recoup}} 
    function by calculating profile matrices from coverages
    which will be used for plotting.
}
\examples{
# Load some data
data("recoup_test_data",package="recoup")
# Do some work
testGenomeRanges <- makeGRangesFromDataFrame(df=test.genome,
    keep.extra.columns=TRUE)
test.input <- coverageRef(
    test.input, 
    genomeRanges=testGenomeRanges, 
    region="tss",
    flank=c(2000,2000)
)
test.input <- profileMatrix(
    test.input,
    flank=c(2000,2000),
    binParams=list(flankBinSize=50,regionBinSize=150,
        sumStat="mean",smooth=TRUE),
    rc=0.5
)
}
\author{
    Panagiotis Moulos
}

