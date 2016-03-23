\name{removeData}
\alias{removeData}
\title{Remove data from recoup list object}
\usage{
    removeData(input, type = c("ranges", "coverage", 
        "profile"))
}
\arguments{
    \item{input}{a list object created from 
    \code{\link{recoup}} or its \code{data} member.}

    \item{type}{one of \code{"ranges"}, \code{"coverage"},
    \code{"profile"}.}
}
\value{
    A list which is normally the output of 
    \code{\link{recoup}} without the members that have 
    been removed from it.
}
\description{
    This function clears members of the \code{\link{recoup}}
    output object that must be cleared in order to apply
    a new set of parameters without completely rerunning
    \code{\link{recoup}}.
}
\details{
    This function clears members of the \code{\link{recoup}}
    output object which typically take some time to be 
    calculated but it is necessary to clean them if the 
    user wants to change input parameters that cause 
    recalculations of these members. For example, if the
    user changes the \code{binParams}, the profile matrices
    (\code{"profile"} object member) have to be 
    recalculated.
    
    \code{type} controls what data will be removed. 
    \code{"ranges"} removes the reads imported from BAM/BED
    files. This is useful when for example the normalization
    method is changed. \code{"coverage"} removes the 
    calculated coverages over the reference genomic regions.
    This is required again when the normalization method 
    changes. \code{"profile"} removes the profile matrices
    derived from coverages. This is required for example 
    when the \code{binParams} main argument changes.
}
\examples{
# Load some data
data("recoup_test_data",package="recoup")

# Before removing
names(test.input)

# Remove a member
test.input <- removeData(test.input,"ranges")

# Removed
names(test.input)
}
\author{
    Panagiotis Moulos
}
