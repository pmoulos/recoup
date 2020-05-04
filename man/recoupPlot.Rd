\name{recoupPlot}
\alias{recoupPlot}
\title{Plot list objects returned by recoup}
\usage{
    recoupPlot(recoupObj, what = c("profile", "heatmap", 
        "correlation"), device = c("x11", "png", "jpg", "tiff", 
        "bmp", "pdf", "ps"), outputDir = ".",
        outputBase = paste(vapply(recoupObj,
        function(x) return(x$data$id), character(1)), 
        sep = "_"), mainh = 1, ...)
}
\arguments{
    \item{recoupObj}{a list object created from 
    \code{\link{recoup}}.}

    \item{what}{one or more of \code{"profile"}, 
    \code{"heatmap"} or \code{"correlation"}. See the 
    \code{plotParams} in the main \code{\link{recoup}} 
    function. A minimum valid version is provided for 
    default plotting.}
    
    \item{device}{a valid R graphics device. See the 
    \code{plotParams} in the main \code{\link{recoup}} 
    function.}
    
    \item{outputDir}{a valid directory when device is not
    \code{"x11"}. See the \code{plotParams} in the main 
    \code{\link{recoup}} function.}
    
    \item{outputBase}{a valid file name to be used as 
    basis when device is not \code{"x11"}. See the 
    \code{plotParams} in the main \code{\link{recoup}} 
    function. Defaults to a concatenation of sample ids.}
    
    \item{mainh}{the reference heatmap for ordering
    operations. Normally, calculated in 
    \code{\link{recoup}}. See also the \code{\link{draw}}
    function in the \code{ComplexHeatmap} package. Deafults
    to the first heatmap.}
    
   \item{...}{further parameters passed either to 
   \code{ggsave} or the base graphics devices of R.}
}
\value{
    This function does not returns anything, just plots the
    \code{recoup} plots.
}
\description{
    This function takes as input argument an output object
    from \code{\link{recoup}} and plots the \code{ggplot2} 
    and \code{ComplexHeatmap} objects stored there.
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
recoupPlot(test.tss)
}
\author{
    Panagiotis Moulos
}
