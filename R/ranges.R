preprocessRanges <- function(input,preprocessParams,bamParams=NULL,rc=NULL) {
    hasRanges <- sapply(input,function(x) is.null(x$ranges))
    if (!any(hasRanges))
        return(input)
    # Else, BAM/BED files must be read, so we check if they exist
    if (!all(sapply(input,function(x) {
        file.exists(x$file)
    })))
        stop("One or more input files cannot be found! Check the validity of ",
            "the file paths.")
    switch(preprocessParams$normalize,
        none = {
            ranges <- cmclapply(input,function(x,pp,p) {
                message("Reading sample ",x$name)
                return(readRanges(x$file,x$format,pp$spliceAction,
                    pp$spliceRemoveQ,pp$bedGenome,params=p))
            },preprocessParams,bamParams,rc=rc)
            names(ranges) <- names(input)
            for (i in 1:length(input))
                input[[i]]$ranges <- ranges[[i]]
        },
        linear = { # Same as none but will process after coverage calculation
            ranges <- cmclapply(input,function(x,pp,p) {
                message("Reading sample ",x$name)
                return(readRanges(x$file,x$format,pp$spliceAction,
                    pp$spliceRemoveQ,pp$bedGenome,params=p))
            },preprocessParams,bamParams,rc=rc)
            names(ranges) <- names(input)
            for (i in 1:length(input))
                input[[i]]$ranges <- ranges[[i]]
        },
        downsample = {
            ranges <- cmclapply(input,function(x,pp,p) {
                message("Reading sample ",x$name)
                return(readRanges(x$file,x$format,pp$spliceAction,
                    pp$spliceRemoveQ,pp$bedGenome,params=p))
            },preprocessParams,bamParams,rc=rc)
            libSizes <- lengths(ranges)
            downto = min(libSizes)
            set.seed(preprocessParams$seed)
            downsampleIndex <- lapply(libSizes,function(x,s) {
                return(sort(sample(x,s)))
            },downto)
            names(downsampleIndex) <- names(input)
            for (i in 1:length(input))
                input[[i]]$ranges <- ranges[[i]][downsampleIndex[[i]]]
        },
        sampleto = {
            ranges <- cmclapply(input,function(x,pp,p) {
                message("Reading sample ",x$name)
                return(readRanges(x$file,x$format,pp$spliceAction,
                    pp$spliceRemoveQ,pp$bedGenome,params=p))
            },preprocessParams,bamParams,rc=rc)
            set.seed(preprocessParams$seed)
            libSizes <- lengths(ranges)
            downsampleIndex <- lapply(libSizes,function(x,s) {
                return(sort(sample(x,s)))
            },preprocessParams$sampleTo)
            names(downsampleIndex) <- names(input)
            for (i in 1:length(input))
                input[[i]]$ranges <- ranges[[i]][downsampleIndex[[i]]]
        }
    )
    return(input)
}

getRegionalRanges <- function(ranges,region,flank) {
    switch(region,
        genebody = {
            w <- width(ranges)
            ranges <- promoters(ranges,upstream=flank[1],downstream=0)
            return(resize(ranges,width=w+flank[1]+flank[2]))
        },
        tss = {
            return(promoters(ranges,upstream=flank[1],downstream=flank[2]))
        },
        tes = {
            tmp <- resize(ranges,width=1,fix="end")
            return(promoters(tmp,upstream=flank[1],downstream=flank[2]))
        },
        custom = {
            if (all(width(ranges)==1))
                return(promoters(ranges,upstream=flank[1],downstream=flank[2]))
            else {
                w <- width(ranges)
                ranges <- promoters(ranges,upstream=flank[1],downstream=0)
                return(resize(ranges,width=w+flank[1]+flank[2]))
            }
        }
    )
}

getFlankingRanges <- function(ranges,flank,dir=c("upstream","downstream")) {
    dir = dir[1]
    if (dir=="upstream")
        return(promoters(ranges,upstream=flank,downstream=0))
    else if (dir=="downstream") {
        return(flank(ranges,width=flank,start=FALSE,both=FALSE))
    }
}

readRanges <- function(input,format,sa,sq,bg,params=NULL) {
    if (format=="bam")
        return(readBam(input,sa,sq,params))
    else if (format=="bed")
        return(readBed(input,bg))
}

readBam <- function(bam,sa=c("keep","remove","split"),sq=0.75,params=NULL) {
    sa = sa[1]
    checkTextArgs("sa",sa,c("keep","remove","split"),multiarg=FALSE)
    checkNumArgs("sq",sq,"numeric",c(0,1),"botheq")
    switch(sa,
        keep = {
            return(trim(as(readGAlignments(file=bam),"GRanges")))
        },
        split = {
            return(trim(unlist(grglist(readGAlignments(file=bam)))))
        },
        remove = {
            reads <- trim(as(readGAlignments(file=bam),"GRanges"))
            qu <- quantile(width(reads),sq)
            rem <- which(width(reads)>qu)
            if (length(rem)>0)
                reads <- reads[-rem]
            message("  Excluded ",length(rem)," reads")
            return(reads)
        }
    )
}

readBed <- function(bed,bg) {
    if (!requireNamespace("GenomeInfoDb"))
        stop("R package GenomeInfoDb is required to retrieve chromosome ",
            "lengths when importing reads in bed files!")
    bed <- trim(import.bed(bed,trackLine=FALSE))
    sf <- GenomeInfoDb::fetchExtendedChromInfoFromUCSC(getUcscOrganism(bg))
    rownames(sf) <- as.character(sf[,1])
    sf <- sf[seqlevels(bed),]
    sf <- Seqinfo(seqnames=sf[,1],seqlengths=sf[,2],
        isCircular=sf[,3],genome=getUcscOrganism(bg))
    seqinfo(bed) <- sf
    return(bed)
}
