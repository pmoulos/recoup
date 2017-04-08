\name{mergeRuns}
\alias{mergeRuns}
\title{Merge recoup outputs of same type}
\usage{
    mergeRuns(..., withDesign = c("auto", "drop"),
		dropPlots = TRUE)
}
\arguments{
    \item{...}{one or more \code{\link{recoup}} output 
    objects.}

    \item{withDesign}{one of \code{"auto"} (default) or 
	\code{"drop"}. Determines how to merge designs. See
	details for further information.}
	
	\item{dropPlots}{if profile and/or heatmap plots are
    attached to the input object(s), they will be 
    recalculated if \code{dropPlots=="TRUE"} (default) or 
    dropped otherwise}.
}
\value{
    A \code{\link{recoup}} output object with as many 
    samples as in \code{...}.
}
\description{
    This function accepts two or more \code{\link{recoup}}
    output objects holding single samples to a merged object
    so that all samples can be used together. This is useful
    when many coverages must be calulated/plotted and memory
    issues do not allow effective parallelization.
}
\details{
    The \code{withDesign} argument controls what should be 
    done if any input has an attached design. The default
    behavioir (\code{"auto"}) will try to do its best to 
    preserve compatible designs. If one or more inputs have
    the same design, it will be applied to the rest of the
    samples. If there is only one design, it will be applied
    to all samples (if you don't want this to happen, choose
    \code{"drop"}). If more than one sample has an attached
    design but these are incompatible (different numbers of
    rows/rownames, columns/columnnames), then all designs
    are dropped. Obviously, \code{withDesign="drop"} drops
    all attached designs and the output object is free of
    a design data frame.
}
\examples{
# Load some data
data("recoup_test_data",package="recoup")

test.input.shift <- test.input
names(test.input.shift) <- paste(names(test.input.shift),"_1",sep="")
test.input.shift[[1]]$id <- paste0(test.input.shift[[1]]$id,"_1")
test.input.shift[[1]]$ranges <- 
	shift(test.input.shift[[1]]$ranges,100)
test.input.shift[[2]]$id <- paste0(test.input.shift[[2]]$id,"_1")
test.input.shift[[2]]$ranges <- 
	shift(test.input.shift[[2]]$ranges,100)

test.tss.1 <- recoup(
    test.input,
    design=NULL,
    region="tss",
    type="chipseq",
    genome=test.genome,
    flank=c(2000,2000),
    selector=NULL,
    rc=0.5
)

test.tss.2 <- recoup(
    test.input.shift,
    design=NULL,
    region="tss",
    type="chipseq",
    genome=test.genome,
    flank=c(2000,2000),
    selector=NULL,
    rc=0.5
)

test.tss <- mergeRuns(test.tss.1,test.tss.2)
}

\author{
    Panagiotis Moulos
}
