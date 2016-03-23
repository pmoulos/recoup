\name{recoupCorrelation}
\alias{recoupCorrelation}
\title{Plot (faceted) average genomic coverage correlations}
\usage{
    recoupCorrelation(recoupObj, samples = NULL, rc = NULL)
}
\arguments{
    \item{recoupObj}{a list object created from 
    \code{\link{recoup}}.}
    
    \item{samples}{which samples to plot. Either numeric
    (denoting the sample indices) or sample ids. Defaults
    to \code{NULL} for all samples.}

    \item{rc}{fraction (0-1) of cores to use in a multicore 
    system. It defaults to \code{NULL} (no parallelization).}
}
\value{
    The function returns the \code{recoupObj} with the
    slot for the correlation plot filled. See also the 
    \code{recoupPlot}, \code{getr} and \code{setr}
    function.
}
\description{
    This function takes as input argument and output object
    from \code{\link{recoup}} and creates the average genomic
    curve correlations according to the options present in the 
    input object. It can be used with saved recoup outputs
    so as to recreate the plots without re-reading BAM/BED
    files and re-calculating coverages.
}
\examples{
# Load some data
data("recoup_test_data",package="recoup")

# Calculate coverages
test.tss <- recoup(
    test.input,
    design=NULL,
    region="tss",
    type="chipseq",
    genome=test.genome,
    flank=c(2000,2000),
    selector=NULL,
    plotParams=list(profile=FALSE,correlation=TRUE,
        heatmap=FALSE),
    rc=0.5
)

# Plot coverage correlations
recoupCorrelation(test.tss,rc=0.5)
}
\author{
    Panagiotis Moulos
}
