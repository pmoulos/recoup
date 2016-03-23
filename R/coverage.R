coverageRef <- function(input,genomeRanges,region=c("tss","tes","genebody",
    "custom"),flank=c(2000,2000),strandedParams=list(strand=NULL,
    ignoreStrand=TRUE),bamParams=NULL,rc=NULL) {
    hasCoverage <- sapply(input,function(x) is.null(x$coverage))
    if (!any(hasCoverage))
        return(input)
    if (region %in% c("tss","tes","custom")) {
        if (region %in% c("tss","tes"))
            input <- coverageBaseRef(input,genomeRanges,region,flank,
                strandedParams,rc=rc)#,bamParams)
        else if (region=="custom") {
            if (all(width(genomeRanges)==1))
                input <- coverageBaseRef(input,genomeRanges,region,flank,
                    strandedParams,rc=rc)#,bamParams)
            else
                input <- coverageAreaRef(input,genomeRanges,region,flank,
                    strandedParams,rc=rc)#,bamParams)
        }
    }
    else # is genebody
        input <- coverageAreaRef(input,genomeRanges,region,flank,strandedParams,
            rc=rc)#,bamParams)
    return(input)
}

coverageBaseRef <- function(input,genomeRanges,region,flank,strandedParams,
    rc=NULL) {
    mainRanges <- getRegionalRanges(genomeRanges,region,flank)
    names(input) <- sapply(input,function(x) return(x$id))
    for (n in names(input)) {
        message("Calculating ",region," coverage for ",input[[n]]$name)
        if (!is.null(input[[n]]$ranges))
            input[[n]]$coverage <- calcCoverage(input[[n]]$ranges,mainRanges,
                strand=strandedParams$strand,
                ignore.strand=strandedParams$ignoreStrand,rc=rc)
        else
            input[[n]]$coverage <- calcCoverage(input[[n]]$file,mainRanges,
                strand=strandedParams$strand,
                ignore.strand=strandedParams$ignoreStrand,rc=rc)
    }
    return(input)
}

coverageAreaRef <- function(input,genomeRanges,region,flank,strandedParams,
    bamParams=NULL,rc=NULL) {
    mainRanges <- getRegionalRanges(genomeRanges,region,flank)
    #leftRanges <- getFlankingRanges(mainRanges,flank[1],"upstream")
    #rightRanges <- getFlankingRanges(mainRanges,flank[2],"downstream")
    
    names(input) <- sapply(input,function(x) return(x$id))
    for (n in names(input)) {
        theRanges <- splitBySeqname(input[[n]]$ranges,rc=rc)
        message("Calculating ",region," coverage for ",input[[n]]$name)
        input[[n]]$coverage <- calcCoverage(theRanges,mainRanges,
            strand=strandedParams$strand,
            ignore.strand=strandedParams$ignoreStrand,rc=rc)
        #message(" center...")
        #input[[n]]$coverage$center <- calcCoverage(theRanges,mainRanges,
        #    strand=strandedParams$strand,
        #    ignore.strand=strandedParams$ignoreStrand,rc=rc)
        #message(" upstream...")
        #input[[n]]$coverage$upstream <- calcCoverage(theRanges,leftRanges,
        #    strand=strandedParams$strand,
        #    ignore.strand=strandedParams$ignoreStrand,rc=rc)
        #message(" downstream...")
        #input[[n]]$coverage$downstream <- calcCoverage(theRanges,rightRanges,
        #    strand=strandedParams$strand,
        #    ignore.strand=strandedParams$ignoreStrand,rc=rc)
    }
    return(input)
}

coverageRnaRef <- function(input,genomeRanges,helperRanges,flank,
    strandedParams=list(strand=NULL,ignoreStrand=TRUE),bamParams=NULL,rc=NULL) {
    hasCoverage <- sapply(input,function(x) is.null(x$coverage))
    if (!any(hasCoverage))
        return(input)
    if (flank[1]==0)
        leftRanges <- getFlankingRanges(helperRanges,1,"upstream")
    else
        leftRanges <- getFlankingRanges(helperRanges,flank[1],"upstream")
    if (flank[1]==0)
        rightRanges <- getFlankingRanges(helperRanges,1,"downstream")
    else
        rightRanges <- getFlankingRanges(helperRanges,flank[2],"downstream")
    names(input) <- sapply(input,function(x) return(x$id))
    for (n in names(input)) {
        theRanges <- splitBySeqname(input[[n]]$ranges)
        message("Calculating genebody coverage for ",input[[n]]$name)
        message(" center...")
        #input[[n]]$coverage$center <- 
        center <- calcCoverage(theRanges,genomeRanges,
            strand=strandedParams$strand,
            ignore.strand=strandedParams$ignoreStrand,rc=rc)
        message(" upstream...")
        #input[[n]]$coverage$upstream <- 
        left <- calcCoverage(theRanges,leftRanges,
            strand=strandedParams$strand,
            ignore.strand=strandedParams$ignoreStrand,rc=rc)
        message(" downstream...")
        #input[[n]]$coverage$downstream <- 
        right <- calcCoverage(theRanges,rightRanges,
            strand=strandedParams$strand,
            ignore.strand=strandedParams$ignoreStrand,rc=rc)
        message("   merging...")
        input[[n]]$coverage <- cmclapply(1:length(center),function(i,le,ce,ri) {
            if (any(is.null(le[[i]]),is.null(ce[[i]]),is.null(ri[[i]])))
                return(NULL)
            else
                return(c(le[[i]],ce[[i]],ri[[i]]))
        },left,center,right,rc=rc)
        names(input[[n]]$coverage) <- names(genomeRanges)
    }
    return(input)
}

calcCoverage <- function(input,mask,strand=NULL,ignore.strand=TRUE,rc=NULL) {
    if (!is(input,"GRanges") && !is.list(input) && is.character(input)
        && !file.exists(input))
        stop("The input argument must be a GenomicRanges object or a valid ",
            "BAM file or a list of GenomicRanges")
    if (!is(mask,"GRanges") && !is(mask,"GRangesList"))
        stop("The mask argument must be a GRanges or GRangesList object")
    isBam <- FALSE
    if (is.character(input) && file.exists(input))
        isBam <- TRUE
    if (!is.null(strand) && !is.list(strand) && !isBam) {
        message("Retrieving ",strand," reads...")
        input <- input[strand(input)==strand]
    }
    if (!is.list(input) && !isBam)
        input <- splitBySeqname(input)
    index <- 1:length(mask)
    if (isBam)
        cov <- cmclapply(index,coverageFromBam,mask,input,ignore.strand,rc=rc)
    else
        cov <- cmclapply(index,coverageFromRanges,mask,input,ignore.strand,
            rc=rc)
    names(cov) <- names(mask)
    
    # Failed effort to chunk the reference genomic ranges, takes more time
    #if (length(index)<=1000)
    #    ii <- splitVector(index,10)
    #else
    #    ii <- splitVector(index,100)
    #if (isBam)
    #    cov <- unlist(lapply(ii,function(i,mask,input,ignore.strand,rc=rc) {
    #        return(cmclapply(i,coverageFromBam,mask,input,ignore.strand,rc=rc))
    #    },mask,input,ignore.strand,rc=rc),recursive=FALSE)
    #    
    #else
    #    cov <- unlist(lapply(1:length(ii),function(i,mask,input,ignore.strand,rc=rc) {
    #        return(cmclapply(i,coverageFromRanges,mask,input,ignore.strand,rc=rc))
    #    },mask,input,ignore.strand,rc=rc),recursive=FALSE)
    
    names(cov) <- names(mask)
    gc(verbose=FALSE)
    return(cov) # Rle
}

coverageFromRanges <- function(i,mask,input,ignore.strand) {
    if (is(mask,"GRangesList"))
        x <- mask[[i]]
    else
        x <- mask[i]
    y<-list(
        chromosome=as.character(seqnames(x))[1],
        start=start(x),
        end=end(x),
        strand=as.character(strand(x))[1],
        reads=NULL,
        coverage=NULL
    )
    if (!is.null(input[[y$chromosome]])) {
        y$reads <- input[[y$chromosome]][
            subjectHits(findOverlaps(x,input[[y$chromosome]],
                ignore.strand=ignore.strand))]
    }
    else {
        message(y$chromosome," not found!")
        y$reads <- NULL
    }
    if (length(y$reads)>0) {
        tryCatch({
            cc <- as.character(seqnames(y$reads))[1]
            y$coverage <- coverage(y$reads)
            if (length(y$start)>1) { # List of exons, RNA, merge exons
                i2k <- unlist(lapply(1:length(y$start),function(j,s,e) {
                    return(s[j]:e[j])
                },y$start,y$end))
                y$coverage <- y$coverage[[cc]][i2k]
            }
            else
                y$coverage <- y$coverage[[cc]][y$start:y$end]
            if (y$strand=="+")
                return(y$coverage)
            else if (y$strand=="-")
                return(rev(y$coverage))
            else
                return(y$coverage)
        },
        error=function(e) {
            message("Caught invalid genomic area!")
            print(mask[i])
            message("Will return zero coverage")
            return(NULL)
        },finally={})
    }
    else
        return(NULL)
}

coverageFromBam <- function(i,mask,input,ignore.strand,pp) {
    if (is(mask,"GRangesList"))
        x <- mask[[i]]
    else
        x <- mask[i]
    y<-list(
        chromosome=as.character(seqnames(x)),
        start=start(x),
        end=end(x),
        strand=as.character(strand(x)),
        reads=NULL,
        coverage=NULL
    )
    bam.file <- input
    bam.index <- paste(input,"bai",sep=".")
    bp <- ScanBamParam(which=x)
    
    switch(pp$spliceAction,
        keep = {
            y$reads <- as(readGAlignments(file=bam.file,index=bam.index,
                param=bp,with.which_label=TRUE),"GRanges")
        },
        split = {
            y$reads <- unlist(grglist(readGAlignments(file=bam.file,
                index=bam.index,param=bp,with.which_label=TRUE)))
        },
        remove = {
            y$reads <- as(readGAlignments(file=bam.file,index=bam.index,
                param=bp,with.which_label=TRUE),"GRanges")
            qu <- quantile(width(y$reads),pp$spliceRemoveQ)
            rem <- which(width(y$reads)>qu)
            if (length(rem)>0)
                y$reads <- y$reads[-rem]
            message("  Excluded ",length(rem)," reads")
        }
    )
    
    if (length(y$reads)>0) {
        tryCatch({
            seqlevels(y$reads) <- as.character(seqnames(x))
            cc <- as.character(seqnames(y$reads))[1]
            y$coverage <- coverage(y$reads)
            if (length(y$start)>1) { # List of exons, RNA, merge exons
                i2k <- unlist(lapply(1:length(y$start),function(j,s,e) {
                    return(s[j]:e[j])
                },y$start,y$end))
                y$coverage <- y$coverage[[cc]][i2k]
            }
            else
                y$coverage <- y$coverage[[cc]][y$start:y$end]
            if (y$strand=="+")
                return(y$coverage)
            else if (y$strand=="-")
                return(rev(y$coverage))
            else
                return(y$coverage)
        },
        error=function(e) {
            message("Caught invalid genomic area!")
            print(mask[i])
            message("Will return zero coverage")
            return(NULL)
        },finally={})
    }
    else
        return(NULL)
}
