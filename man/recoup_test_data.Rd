\docType{data}
\name{recoup_test_data}
\alias{recoup_test_data}
\alias{test.design}
\alias{test.exons}
\alias{test.genome}
\alias{test.input}
\title{Reference and genomic sample regions for recoup testing}
\format{\code{data.frame} and \code{list} objects whose format 
is accepted by recoup.}
\source{
    Personal communication with the Talianids lab at BSRC 
    'Alexander Fleming'. Unpublished data.
}
\description{
    The testing data package containes a small gene set, a 
    design data frame, some genomic regions and an input 
    object for testing of recoup with ChIP-Seq and RNA-Seq 
    data. Specifically:
    \itemize{
        \item \code{test.input}: A small data set which 
        contains 10000 reads from H4K20me1 ChIP-Seq 
        data from WT adult mice and Set8 (Pr-Set7) KO 
        mice. The tissue is liver.
        \item \code{test.genome}: A small gene set (100 
        genes) and their coordinates from mouse mm9 
        chromosome 12.
        \item{test.design}: A data frame containing the 100 
        above genes categorized according to expression and 
        strand.
        \item{test.exons}: A \code{GRangesList} containing 
        the exons of the 100 above genes for use with 
        \code{recoup} RNA-Seq mode.
    }
}
\author{
    Panagiotis Moulos
}
\keyword{datasets}
