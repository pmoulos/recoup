\name{getAnnotation}
\alias{getAnnotation}
\title{Annotation downloader}
\usage{
    getAnnotation(org, type, refdb = "ensembl", ver = NULL,
        rc = NULL)
}
\arguments{
    \item{org}{the organism for which to download
    annotation. Check the main \code{\link{recoup}}
    help page for details on supported organisms.}

    \item{type}{either \code{"gene"} or \code{"exon"}.}

    \item{refdb}{the online source to use to fetch 
    annotation. It can be \code{"ensembl"} (default), 
    \code{"ucsc"} or \code{"refseq"}. In the later two
    cases, an SQL connection is opened with the UCSC 
    public databases.}
    
    \item{ver}{the version of \code{org} to use as
    related to \code{refdb} or \code{NULL} for latest
    versions. See also the main \code{recoup} help page.}

    \item{rc}{fraction (0-1) of cores to use in a multicore 
    system. It defaults to \code{NULL} (no parallelization).
    It is used in the case of \code{type="exon"} to process
    the return value of the query to the UCSC Genome 
    Browser database.}
}
\value{
    A data frame with the canonical (not isoforms!) genes or
    exons of the requested organism. When
    \code{type="genes"}, the data frame has the following
    columns: chromosome, start, end, gene_id, gc_content,
    strand, gene_name, biotype. When \code{type="exon"} the
    data frame has the following columns: chromosome, start,
    end, exon_id, gene_id, strand, gene_name, biotype. The
    gene_id and exon_id correspond to Ensembl gene and exon
    accessions respectively. The gene_name corresponds to
    HUGO nomenclature gene names.
}
\description{
    This function connects to the EBI's Biomart service using
    the package biomaRt and downloads annotation elements
    (gene co-ordinates, exon co-ordinates, gene
    identifications, biotypes etc.) for each of the supported
    organisms. See the help page of \code{\link{recoup}}
    for a list of supported organisms. The function downloads
    annotation for an organism genes or exons. It also uses
    the UCSC public database connection API to download UCSC
    and RefSeq annotations.
}
\note{
    The data frame that is returned contains only "canonical"
    chromosomes for each organism. It does not contain
    haplotypes or random locations and does not contain
    chromosome M.
}
\examples{
    mm9.genes <- getAnnotation("mm9","gene")
}
\author{
    Panagiotis Moulos
}

