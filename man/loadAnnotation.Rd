\name{loadAnnotation}
\alias{loadAnnotation}
\title{Load a recoup annotation element}
\usage{
    loadAnnotation(genome, refdb, 
        type = c("gene", "exon", "utr"), version="auto",
        db = file.path(system.file(package = "recoup"),
            "annotation.sqlite"), summarized = FALSE, 
            asdf = FALSE, rc = NULL)
}
\arguments{
    \item{genome}{a \code{\link{recoup}} supported
    organisms or a custom, imported by the user, name. See 
    also the main \code{\link{recoup}} man page.}

    \item{refdb}{a \code{\link{recoup}} supported
    annotation source or a custom, imported by the user, name.
    See also the main \code{\link{recoup}} man page.}
    
    \item{type}{one of the \code{"gene"}, \code{"exon"} 
    or \code{"utr"}.}
    
    \item{version}{same as the \code{version} in 
    \code{\link{recoup}}.}
    
    \item{db}{same as the \code{db} in 
    \code{\link{buildAnnotationDatabase}}.}
    
    \item{summarized}{if \code{TRUE}, retrieve summarized,
    non-overlaping elements where appropriate (e.g. exons).}
    
    \item{asdf}{return the result as a \code{\link{data.frame}}
    (default \code{FALSE}).}
    
    \item{rc}{same as the \code{rc} in 
    \code{\link{buildAnnotationDatabase}}.}
}
\value{
    The function returns a \code{GenomicRanges} object with
    the requested annotation.
}
\description{
    This function creates loads an annotation element from
    the local annotation database to be used with recoup.
    If the annotation is not found and the organism is 
    supported, the annotation is created on the fly but not
    imported in the local database. Use
    \code{buildAnnotationDatabase} for this purpose.
    
}
\examples{
db <- file.path(system.file(package="recoup"),
  "annotation.sqlite")
if (file.exists(db))
  gr <- loadAnnotation(genome="hg19",refdb="ensembl",
    type="gene",db=db)
}
\author{
    Panagiotis Moulos
}
