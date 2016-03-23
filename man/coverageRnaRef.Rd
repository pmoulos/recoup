\name{coverageRnaRef}
\alias{coverageRnaRef}
\title{Calculate coverage in a set of reference genomic 
    regions (RNA-Seq or spliced mode)}
\usage{
    coverageRnaRef(
        input, 
        genomeRanges, 
        helperRanges,
        flank,
        strandedParams = list(strand=NULL, ignoreStrand=TRUE),
        bamParams = NULL,
        rc = NULL
    )
}
\arguments{
    \item{input}{an input list as in \code{\link{recoup}} but
    with the \code{ranges} field of each member filled 
    (e.g. after using \code{\link{preprocessRanges}}).}

    \item{genomeRanges}{a named \code{GRangesList} where
    list member names are genes and list members are 
    \code{GRanges} representing each gene's exons.}
    
    \item{helperRanges}{the \code{genome} from 
    \code{\link{recoup}} as a \code{GRanges} object (e.g. 
    the output from 
    \code{\link{makeGRangesFromDataFrame}}).}
        
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

# Note: the figures that will be produced will not look 
# realistic or pretty and will be "bumpy". This is because 
# package size limitations posed by Bioconductor guidelines 
# do not allow for a full test dataset. As a result, the input
# below is not an RNA-Seq dataset. Have a look at the
# vignette on how to test with more realistic data.

# Calculate coverages
testGenomeRanges <- makeGRangesFromDataFrame(df=test.genome,
    keep.extra.columns=TRUE)
test.input <- coverageRnaRef(
    test.input,
    genomeRanges=test.exons,
    helperRanges=testGenomeRanges,
    flank=c(2000,2000)
)
}
\author{
    Panagiotis Moulos
}

