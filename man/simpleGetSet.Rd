\name{simpleGetSet}
\alias{getr}
\alias{setr}
\title{Get and set some reusable objects from a recoup 
    object}
\usage{
    getr(obj, key = c("design", "profile", "heatmap", 
        "correlation", "orderBy", "kmParams", "plotParams"))
    setr(obj, key, value = NULL)
}
\arguments{
    \item{obj}{a list object created from 
    \code{\link{recoup}}.}

    \item{key}{one of \code{"design"}, \code{"profile"},
    \code{"heatmap"}, \code{"correlation"}, \code{orderBy}, 
    \code{kmParams}, \code{plotParams}. For 
    \code{"profile"}, \code{"heatmap"}, the respective 
    plots are retrieved or changed according to which 
    function is called. For \code{setr} it can (and 
    preferably) be a named list of arguments to be changed
    in the \code{recoup} list object. The list names are
    the same as above. For the rest, see the main 
    \code{\link{recoup}} man page.}
    
    \item{value}{a valid \code{ggplot} or 
    \code{HeatmapList} object created from 
    \code{\link{recoupProfile}} or 
    \code{\link{recoupHeatmap}} 
    \code{\link{recoupCorrelation}} when changing plots. 
    Values for all other types are also checked for 
    validity.}
}
\value{
    For \code{getr}, the object asked to be retrieved. For 
    \code{setr}, the \code{obj} with the respective slots 
    filled or replaced with \code{value}.
}
\description{
    The \code{getr} and \code{setr} functions are used to 
    get several reusable/changeable objects of 
    \code{\link{recoup}} or replcace them (e.g. when the 
    user wishes to change some \code{ggplot} or 
    \code{ComplexHeatmap} parameters manually in a plot, or
    change the heatmap profile ordering mode).
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
    rc=0.5
)

# Plot coverage profiles

# Get the curve profile plot
pp <- getr(test.tss,"profile")

# Change some ggplot parameter
pp <- pp + 
    theme(axis.title.x=element_text(size=14))
    
# Store the new plot
test.tss <- setr(test.tss,"profile",pp)
## or even better
# test.tss <- setr(test.tss,list(profile=pp))
}
\author{
    Panagiotis Moulos
}
