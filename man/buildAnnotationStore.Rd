\name{buildAnnotationStore}
\alias{buildAnnotationStore}
\title{Build a local annotation database for recoup}
\usage{
    buildAnnotationStore(organisms, sources,
    home = file.path(path.expand("~"), ".recoup"),
    forceDownload = TRUE, rc = NULL)
}
\arguments{
    \item{organisms}{a character vector of organisms
    for which to download and build annotations. Check 
    the main \code{\link{recoup}} help page for details 
    on supported organisms.}

    \item{sources}{a character vector of public sources
    from which to download and build annotations. Check 
    the main \code{\link{recoup}} help page for details 
    on supported annotation sources.}
    
    \item{home}{a valid path (accessible at least by the
    current user) where the annotation database will be 
    set up. It defaults to \code{".recoup"} inside the
    current user's home directory.}
    
    \item{forceDownload}{by default, 
    \code{buildAnnotationStore} will not download an
    existing annotation again (\code{FALSE}). Set to 
    \code{TRUE} if you wish to update the annotation 
    database.}

    \item{rc}{fraction (0-1) of cores to use in a multicore 
    system. It defaults to \code{NULL} (no parallelization).
    It is used in the case of \code{type="exon"} to process
    the return value of the query to the UCSC Genome 
    Browser database.}
}
\value{
    The function does not return anything. Only the annotation
    directory and contents are created.
}
\description{
    This function creates a local annotation database to be
    used with recoup so as to avoid long time on the fly 
    annotation downloads and formatting.
}
\examples{
\donttest{
    buildAnnotationStore("mm9","ensembl")
}
}
\author{
    Panagiotis Moulos
}

