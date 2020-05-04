\name{kmeansDesign}
\alias{kmeansDesign}
\title{Apply k-means clustering to profile data}
\usage{
    kmeansDesign(input, design = NULL, kmParams)
}
\arguments{
    \item{input}{a list object created from 
    \code{\link{recoup}} or partially processed by
    \code{\link{recoup}} or its \code{data} member.
    See the main input to \code{\link{recoup}} for 
    further information.}

    \item{design}{See the respective argument in 
    \code{\link{recoup}} for further information}
    
    \item{kmParams}{Contains parameters for k-means
    clustering on profiles. See the respective argument 
    in \code{\link{recoup}} for further information.}
}
\value{
    The design data frame, either created from scratch or
    augmented by k-means clustering.
}
\description{
    This function performs k-means clustering on 
    \code{\link{recoup}} generated profile matrices and
    stores the result as a factor in the design element.
    If no design is present, then one is created from the
    k-means result.
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
    flank=c(1000,1000),
    selector=NULL,
    plotParams=list(plot=FALSE,profile=TRUE,
        heatmap=TRUE,device="x11"),
    rc=0.1
)

# Re-design based on k-means
kmParams=list(k=2,nstart=20,algorithm="MacQueen",iterMax=20,
    reference=NULL)
design <- kmeansDesign(test.tss$data,kmParams=kmParams)
}
\author{
    Panagiotis Moulos
}
