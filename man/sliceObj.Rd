\name{sliceObj}
\alias{sliceObj}
\title{Subset recoup output list objects}
\usage{
    sliceObj(obj, i = NULL, j = NULL, k = NULL, 
        dropPlots = FALSE, rc = NULL)
}
\arguments{
    \item{obj}{a list object created from 
    \code{\link{recoup}}.}

    \item{i}{vector of numeric or character indices, 
    corresponding to the index or rownames or names of
    reference genomic regions. The \code{design} object
    member will also be subset. If there is a 
    \code{selector} attached to the input object (see
    \code{\link{recoup}} arguments) it will be dropped.}
    
    \item{j}{vector of numeric indices corresponding to 
    the profile matrix vertical index (or base pair 
    position or bin of base pairs) so as to subset the
    profile. The function will do its best to "guess"
    new plotting x-axis labels.}
    
    \item{k}{vector of numeric or character indices
    corresponding to sample index or sample names. These
    will be returned.}
    
    \item{dropPlots}{if profile and/or heatmap plots are
    attached to the input object, they will be recalculated
    if \code{dropPlots=="TRUE"} (default) or dropped 
    otherwise}.
    
    \item{rc}{fraction (0-1) of cores to use in a multicore 
    system. It defaults to \code{NULL} (no parallelization).}
}
\value{
    A \code{recoup} list object, susbet according to 
    \code{i, j}.
}
\description{
    This function takes as input argument an output object
    from \code{\link{recoup}} and subsets it according to the
    inputs \code{i,j,k}. The attached plots may or may not be
    recalculated. Other input parameters stores in 
    \code{obj$callopts} are not changed apart from any
    \code{selector} option which is dropped. Note that when
    slicing vertically (by \code{j}), the \code{$coverage}
    member of the input data (if present) is \bold{not} 
    sliced, but remains as is. You can drop it using
    \code{\link{removeData}} as it is used to recalculate
    profile matrices only if bin sizes are changed in a 
    \code{\link{recoup}} call.
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
    plotParams=list(plot=FALSE,profile=TRUE,
        heatmap=TRUE,device="x11"),
    rc=0.1
)

# Plot coverage profiles
o <- sliceObj(test.tss,i=1:10,k=1)
}
\author{
    Panagiotis Moulos
}
