preprocessRanges <- function(input,preprocessParams,genome,bamRanges=NULL,
    bamParams=NULL,rc=NULL) {
    hasRanges <- sapply(input,function(x) is.null(x$ranges))
    if (!any(hasRanges))
        return(input)
    # Else, BAM/BED files must be read, so we check if they exist
    if (!all(sapply(input,function(x) {
        file.exists(x$file)
    })))
        stop("One or more input files cannot be found! Check the validity of ",
            "the file paths.")
    
    # Since this function is exported, some check for preprocessParams is
    # required
    if (is.null(preprocessParams$fragLen))
        preprocessParams$fragLen <- NA
    if (is.null(preprocessParams$cleanLevel))
        preprocessParams$cleanLevel <- 0
    if (is.null(preprocessParams$normalize))
        preprocessParams$normalize <- "none"
    
    # Define the BAM ranges to read if present. If yes, files need to be sorted,
    # indexed etc.
    if (!is.null(bamRanges)) {
        # BAM ranges must be extended by fragLen if fragLen, otherwise we lose
        # reads falling in this small area which causes a signal drop artifact
        # on the edges of the areas to plot coverage for
        if (!is.na(preprocessParams$fragLen)) {
            w <- width(bamRanges)
            bamRanges <- promoters(bamRanges,upstream=preprocessParams$fragLen,
                downstream=0)
            bamRanges <- resize(bamRanges,width=w+2*preprocessParams$fragLen)
        }
        bi <- all(sapply(input,function(x) {
            file.exists(paste0(x$file,".bai"))
        }))
        if (!all(bi)) {
            nbi <- which(!bi)
            cmclapply(nbi,prepareBam,input,rc=rc)
        }
        if (is.null(bamParams))
            bamParams <- ScanBamParam(which=bamRanges)
        else
            bamWhich(bamParams) <- bamRanges
    }
    
    # Read the ranges
    ranges <- cmclapply(input,function(x,pp,p) {
        message("Reading sample ",x$name)
        return(readRanges(x$file,x$format,pp$spliceAction,
            pp$spliceRemoveQ,pp$bedGenome,params=p))
    },preprocessParams,bamParams,rc=rc)
    names(ranges) <- names(input)
    
    # Resize if requested
    if (!is.na(preprocessParams$fragLen))
        ranges <- resizeRanges(ranges,preprocessParams$fragLen,rc=rc)
    
    # With 0 we do nothing
    # With 1, remove unanchored regions, keep chrM
    if (preprocessParams$cleanLevel==1) {
        message("Removing unanchored reads from all samples")
        ranges <- cmclapply(ranges,cleanRanges,1,genome,rc=rc)
    }
    # With 2, remove unanchored regions and chrM
    else if (preprocessParams$cleanLevel==2) {
        message("Removing unanchored and mitochondrial reads from all samples")
        ranges <- cmclapply(ranges,cleanRanges,2,genome,rc=rc)
    }
    # With 3, remove unanchored regions, chrM and unique reads
    else if (preprocessParams$cleanLevel==3) {
        message("Removing unanchored, mitochondrial and duplicate reads ",
                "from all samples")
        ranges <- cmclapply(ranges,cleanRanges,3,genome,rc=rc)
    }
        
    switch(preprocessParams$normalize,
        none = {
            for (i in 1:length(input))
                input[[i]]$ranges <- ranges[[i]]
        },
        linear = { # Same as none but will process after coverage calculation
            for (i in 1:length(input))
                input[[i]]$ranges <- ranges[[i]]
        },
        downsample = {
            libSizes <- lengths(ranges)
            downto = min(libSizes)
            set.seed(preprocessParams$seed)
            downsampleIndex <- lapply(libSizes,function(x,s) {
                return(sort(sample(x,s)))
            },downto)
            names(downsampleIndex) <- names(input)
            for (i in 1:length(input)) {
                message("Downsampling sample ",input[[i]]$name," to ",downto,
                    " reads")
                input[[i]]$ranges <- ranges[[i]][downsampleIndex[[i]]]
            }
        },
        sampleto = {
            set.seed(preprocessParams$seed)
            libSizes <- lengths(ranges)
            downsampleIndex <- lapply(libSizes,function(x,s) {
                if (s>x) s <- x
                return(sort(sample(x,s)))
            },preprocessParams$sampleTo)
            names(downsampleIndex) <- names(input)
            for (i in 1:length(input)) {
                message("Sampling sample ",input[[i]]$name," to ",
                    preprocessParams$sampleTo," reads")
                input[[i]]$ranges <- ranges[[i]][downsampleIndex[[i]]]
            }
        }
    )
    
    return(input)
}

getMainRanges <- function(genomeRanges,helperRanges=NULL,type,region,flank,
    rc=NULL) {
    if (type=="rnaseq" && is.null(helperRanges))
        stop("helperRanges must be supplied when type is \"rnaseq\"")
        
    message("Getting main ranges for measurements")
    message("  measurement type: ", type)
    message("  genomic region type: ", region)
    
    if (type=="chipseq") {
        mainRanges <- getRegionalRanges(genomeRanges,region,flank)
        return(list(mainRanges=mainRanges,bamRanges=mainRanges))
    }
    else if (type=="rnaseq") {
        bamRanges <- getRegionalRanges(helperRanges,region,flank)
        return(list(mainRanges=genomeRanges,bamRanges=bamRanges))
    }
}

getMainRnaRangesOnTheFly <- function(genomeRanges,flank,rc=NULL) {
    message("Creating summarized exon flanking region ",flank[1]," bps and ",
        flank[2]," bps")
    if (is.null(rc)) {
        # For some progress recording...
        flankedSexon <- lapply(genomeRanges,flankFirstLast,
            f[1],f[2],rc=rc)
    }
    else
        flankedSexon <- cmclapply(genomeRanges,flankFirstLast,
            flank[1],flank[2],rc=rc)
    message("Creating a GRangesList with the flanking regions")
    return(GRangesList(flankedSexon))
}

#getMainRnaRangesOnTheFly <- function(helperRanges,flank,rc=NULL) {
#    message("Creating summarized exon flanking region ",flank[1]," bps and ",
#        flank[2]," bps")
#    leftRanges <- getFlankingRanges(helperRanges,flank[1],"upstream")
#    elementMetadata(leftRanges) <- elementMetadata(leftRanges)[,c(2,1,3,4)]
#    names(elementMetadata(leftRanges))[1] <- "exon_id"
#    leftRanges <- as(leftRanges,"GRangesList")
#    rightRanges <- getFlankingRanges(helperRanges,flank[2],"downstream")
#    elementMetadata(rightRanges) <- elementMetadata(rightRanges)[,c(2,1,3,4)]
#    names(elementMetadata(rightRanges))[1] <- "exon_id"
#    rightRanges <- as(rightRanges,"GRangesList")
#    if (is.null(rc)) {
#        # For some progress recording...
#        flankedSexon <- GRangesList()
#        for (i in 1:length(genomeRanges)) {
#            if (i%%1000 == 0)
#                message("  processed ",i," ranges")
#            flankedSexon[[i]] <- c(
#                leftRanges[[i]],
#                genomeRanges[[i]],
#                rightRanges[[i]]
#            )
#        }
#    }
#    else {
#        # Dangerous... must issue warning
#        warning("Parallel creation of GRangesList object may cause a memory ",
#            "leak. Please monitor your system.",immediate.=TRUE)
#        flankedSexon <- cmcmapply(c,leftRanges,helperRanges,
#            rightRanges,rc=rc)
#    }
#    return(flankedSexon)
#}

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

resizeRanges <- function(ranges,fragLen,rc=NULL) {
    if (!is.null(names(ranges)))
        return(cmclapply(names(ranges),function(n,dat,fl) {
            message("Resizing ",n," to ",fl," bases")
            return(trim(resize(dat[[n]],width=fl,fix="start")))
        },ranges,fragLen,rc=rc))
    else
        return(cmclapply(1:length(ranges),function(i,dat,fl) {
            message("Resizing rangeset",i," to ",fl," bases")
            return(trim(resize(dat[[i]],width=fl,fix="start")))
        },ranges,fragLen,rc=rc))
}

flankFirstLast <- function(x,u,d) {
    s <- as.character(strand(x[1]))
    n <- length(x)
    w1 <- width(x[1])
    w2 <- width(x[n])
    if (s=="+") {
        x[1] <- promoters(x[1],upstream=u,downstream=0)
        x[1] <- resize(x[1],width=w1+u)
        x[n] <- resize(x[n],width=w2+d)
    }
    else if (s=="-") {
        x[1] <- resize(x[1],width=w1+u)
        x[n] <- promoters(x[n],upstream=d,downstream=0)
        x[n] <- resize(x[n],width=w2+d)
    }
    else {
        x[1] <- promoters(x[1],upstream=u,downstream=0)
        x[1] <- resize(x[1],width=w1+u)
        x[n] <- resize(x[n],width=w2+d)
    }
    return(x)
}

cleanRanges <- function(aRange,level,org) {
    if (level==1) {
        chrs <- getValidChrsWithMit(org)
        aRange <- aRange[seqnames(aRange) %in% chrs]
        newsi <- which(seqlevels(aRange) %in% chrs)
        seqlevels(aRange) <- seqlevels(aRange)[newsi]
        seqinfo(aRange) <- seqinfo(aRange)[seqlevels(aRange)]
    }
    else if (level==2) {
        chrs <- getValidChrs(org)
        aRange <- aRange[seqnames(aRange) %in% chrs]
        newsi <- which(seqlevels(aRange) %in% chrs)
        seqlevels(aRange) <- seqlevels(aRange)[newsi]
        seqinfo(aRange) <- seqinfo(aRange)[seqlevels(aRange)]
    }
    else if (level==3) {
        chrs <- getValidChrs(org)
        aRange <- aRange[seqnames(aRange) %in% chrs]
        newsi <- which(seqlevels(aRange) %in% chrs)
        seqlevels(aRange) <- seqlevels(aRange)[newsi]
        seqinfo(aRange) <- seqinfo(aRange)[seqlevels(aRange)]
        aRange <- unique(aRange)
    }
    return(aRange)
}

readRanges <- function(input,format,sa,sq,bg,params=NULL) {
    if (format=="bam")
        return(readBam(input,sa,sq,params))
    else if (format=="bed")
        return(readBed(input,bg))
    else if (format=="bigwig")
        return(NULL)
}

readBam <- function(bam,sa=c("keep","remove","split"),sq=0.75,params=NULL) {
    sa = sa[1]
    checkTextArgs("sa",sa,c("keep","remove","split"),multiarg=FALSE)
    checkNumArgs("sq",sq,"numeric",c(0,1),"botheq")
    switch(sa,
        keep = {
            return(trim(as(readGAlignments(file=bam,param=params),"GRanges")))
        },
        split = {
            return(trim(unlist(grglist(readGAlignments(file=bam,
                param=params)))))
        },
        remove = {
            reads <- trim(as(readGAlignments(file=bam,param=params),"GRanges"))
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

prepareBam <- function(i,input) {
    tryCatch({
        # Will fail if BAM unsorted
        message("Indexing BAM file ",input[[i]]$file)
        indexBam(input[[i]]$file)
    },error=function(e) {
        warning("Caught error ",e," while indexing BAM file ",
            input[[i]]$file,"! Will try to sort now...",
            immediate.=TRUE)
        message("Sorting BAM file ",input[[i]]$file)
        file.rename(input[[i]]$file,paste0(input[[i]]$file,".uns"))
        ff <- sub(pattern="(.*)\\..*$",replacement="\\1",input[[i]]$file)
        sortBam(paste0(input[[i]]$file,".uns"),ff)
        file.remove(paste0(input[[i]]$file,".uns"))
        message("Indexing BAM file ",input[[i]]$file)
        indexBam(input[[i]]$file)
    },finally="")
}
