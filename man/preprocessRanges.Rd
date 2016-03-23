\name{preprocessRanges}
\alias{preprocessRanges}
\title{Read and preprocess BAM/BED files to GRanges}
\usage{
    preprocessRanges(input, preprocessParams, 
        bamParams = NULL, rc = NULL)
}
\arguments{
    \item{input}{an input list as in \code{\link{recoup}} but
    with the \code{ranges} field of each member filled 
    (e.g. after using \code{\link{preprocessRanges}}).}

    \item{preprocessParams}{see the \code{preprocessParams} 
    argument in the main \code{\link{recoup}} function.}
    
    \item{bamParams}{see the \code{bamParams} argument in 
    the main \code{\link{recoup}} function.}

    \item{rc}{fraction (0-1) of cores to use in a multicore 
    system. It defaults to \code{NULL} (no parallelization).}
}
\value{
    This function fills the \code{ranges} field in the
    main \code{input} argument in \code{\link{recoup}} 
    function.
}
\description{
    This function reads the BAM/BED files present in the 
    input list object and fills the \code{ranges} field of 
    the latter. At the same time it takes care of certain
    preprocessing steps like normalization.
}
\examples{
    # This example only demonstrates the usage of the
    # preprocessRanges function. The input BAM files
    # included with the package will not produce
    # realistic plots as they contain only a very small
    # subset of the original data presented in the 
    # vignettes (50k reads). Please see recoup vignettes 
    # for further demonstrations.
    test.in <- list(
        WT_H4K20me1=list(
            id="WT_H4K20me1",
            name="WT H4K20me1",
            file=system.file("extdata",
                "WT_H4K20me1_50kr.bam", 
                package="recoup"),
            format="bam",
            color="#EE0000"
        ),
        Set8KO_H4K20me1=list(
            id="Set8KO_H4K20me1",
            name="Set8KO H4K20me1",
            file=system.file("extdata", 
                "Set8KO_H4K20me1_50kr.bam", 
                package="recoup"),
            format="bam",
            color="#00BB00"
        )
    )
    
    pp=list(
        normalize="none",
        spliceAction="split",
        spliceRemoveQ=0.75
    )
    
    test.in <- preprocessRanges(test.in,pp)
}
\author{
    Panagiotis Moulos
}

