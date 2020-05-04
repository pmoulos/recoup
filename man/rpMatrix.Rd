\name{rpMatrix}
\alias{rpMatrix}
\title{Reads per million profile matrices for plotting}
\usage{
    rpMatrix(input, mainRanges, flank, binParams,
        strandedParams = list(strand = NULL, ignoreStrand = TRUE),
        rc = NULL)
}
\arguments{
    \item{input}{an input list as in \code{\link{recoup}} but
    with the \code{ranges} field of each member filled (e.g. 
    after using \code{\link{preprocessRanges}}.}
    
    \item{mainRanges}{the \code{genome} from 
    \code{\link{recoup}} as a \code{GRanges} object (e.g. 
    the output from 
    \code{\link{makeGRangesFromDataFrame}}).}
    
    \item{flank}{see the \code{flank} argument in the 
    main \code{\link{recoup}} function.}
    
    \item{binParams}{see the \code{binParams} argument in 
    the main \code{\link{recoup}} function.}
    
    \item{strandedParams}{see the \code{strandedParams} 
    argument in the main \code{\link{recoup}} function.}

    \item{rc}{fraction (0-1) of cores to use in a multicore 
    system. It defaults to \code{NULL} (no parallelization).}
}
\value{
    Same as input with the \code{profile} fields filled.
}
\description{
    This function fills the \code{profile} field in the
    main \code{input} argument in \code{\link{recoup}} 
    function by calculating profile matrices using reads
    per million (rpm) or reads per kb per million reads
    (rpkm) over binned genomic areas of interest, instead
    of genomic coverage signals. The profile matrices are
    used for later plotting.
}
\details{
    Regarding the calculation of \code{rpm} and \code{rpkm}
    values, the calculations slightly differ from the 
    default defintions of these measurements in the sense
    that they are also corrected for the bin lengths so as
    to acquire human-friendly values for plotting.
    
    Note that the genomic ranges (BAM/BED files) must be 
    imported before using this function (as per the default
    \code{recoup} pipeline). We plan to support streaming 
    directly from BAM files in the future.
}
\examples{
# Load some data
data("recoup_test_data",package="recoup")
# Do some work
testGenomeRanges <- makeGRangesFromDataFrame(df=test.genome,
    keep.extra.columns=TRUE)
w <- width(testGenomeRanges)
testGenomeRanges <- promoters(testGenomeRanges,upstream=2000,downstream=0)
testGenomeRanges <- resize(testGenomeRanges,width=w+4000)
test.input <- rpMatrix(
    test.input, 
    mainRanges=testGenomeRanges,
    flank=c(2000,2000),
    binParams=list(flankBinSize=50,regionBinSize=150,binType="fixed",
        sumStat="mean",interpolation="spline"),
    rc=0.1
)
}
\author{
    Panagiotis Moulos
}
