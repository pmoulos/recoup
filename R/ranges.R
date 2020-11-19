preprocessRanges <- function(input,preprocessParams,genome,bamRanges=NULL,
    bamParams=NULL,rc=NULL) {
    hasRanges <- vapply(input,function(x) is.null(x$ranges),logical(1))
    if (!any(hasRanges))
        return(input)
    # Else, BAM/BED files must be read, so we check if they exist
    if (!all(vapply(input,function(x) {
        file.exists(x$file)
    },logical(1))))
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
        # Find the first available BAM file as input may be mixed
        bfsl <- vapply(input,function(x) {
            if (!is.null(x$format)) {
                if (x$format == "bam")
                    return(TRUE)
            }
            else {
                if (grepl("\\.bam$",x$file,ignore.case=TRUE,perl=TRUE))
                    return(TRUE)
            }
            return(FALSE)
        },logical(1))
        bfsi <- which(bfsl)
        if (length(bfsi) > 0) {
            fbi <- bfsi[1]
            fb <- input[[fbi]]$file
            # Subset the bamRanges according to the BAM header
            bamRanges <- .subsetGRangesByBamHeader(bamRanges,fb)
        }
        
        # BAM ranges must be extended by fragLen if fragLen, otherwise we lose
        # reads falling in this small area which causes a signal drop artifact
        # on the edges of the areas to plot coverage for
        if (!is.na(preprocessParams$fragLen)) {
            w <- width(bamRanges)
            bamRanges <- promoters(bamRanges,upstream=preprocessParams$fragLen,
                downstream=0)
            bamRanges <- resize(bamRanges,width=w+2*preprocessParams$fragLen)
        }
        bi <- vapply(input,function(x) {
            if (!is.null(x$format)) {
                if (x$format == "bigwig")
                    return(TRUE)
            }
            else {
                if (grepl("\\.(bigwig|bw|wig|bg)$",x$file,ignore.case=TRUE,
                    perl=TRUE))
                    return(TRUE)
            }
            return(file.exists(paste0(x$file,".bai")))
        },logical(1))
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
    #if (is.character(genome)) {
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
    #}
        
    switch(preprocessParams$normalize,
        none = {
            for (i in seq_len(length(input)))
                input[[i]]$ranges <- ranges[[i]]
        },
        linear = { # Same as none but will process after coverage calculation
            for (i in seq_len(length(input)))
                input[[i]]$ranges <- ranges[[i]]
        },
        downsample = {
            libSizes <- lengths(ranges)
            downto = min(libSizes)
            #set.seed(preprocessParams$seed)
            downsampleIndex <- lapply(libSizes,function(x,s) {
                return(sort(sample(x,s)))
            },downto)
            names(downsampleIndex) <- names(input)
            for (i in seq_len(length(input))) {
                message("Downsampling sample ",input[[i]]$name," to ",downto,
                    " reads")
                input[[i]]$ranges <- ranges[[i]][downsampleIndex[[i]]]
            }
        },
        sampleto = {
            #set.seed(preprocessParams$seed)
            libSizes <- lengths(ranges)
            downsampleIndex <- lapply(libSizes,function(x,s) {
                if (s>x) s <- x
                return(sort(sample(x,s)))
            },preprocessParams$sampleTo)
            names(downsampleIndex) <- names(input)
            for (i in seq_len(length(input))) {
                message("Sampling sample ",input[[i]]$name," to ",
                    preprocessParams$sampleTo," reads")
                input[[i]]$ranges <- ranges[[i]][downsampleIndex[[i]]]
            }
        }
    )
    
    return(input)
}

getMainRanges <- function(genomeRanges,helperRanges=NULL,type,region,flank,
    onFlankFail,rc=NULL) {
    if (type=="rnaseq" && is.null(helperRanges))
        stop("helperRanges must be supplied when type is \"rnaseq\"")
        
    message("Getting main ranges for measurements")
    message("  measurement type: ", type)
    message("  genomic region type: ", region)
    
    offB <- FALSE
    if (type=="chipseq") {
        suppressWarnings(mainRanges <- 
            getRegionalRanges(genomeRanges,region,flank))
        if (any(start(mainRanges) < 1)) {
            offB <- TRUE
            mainRanges <- .dropOrTrimGR(mainRanges,onFlankFail)
        }
        return(list(mainRanges=mainRanges,bamRanges=mainRanges,offBound=offB))
    }
    else if (type=="rnaseq") {
        suppressWarnings(bamRanges <- 
            getRegionalRanges(helperRanges,region,flank))
        if (any(start(bamRanges) < 1) || which(end(bamRanges)>seqlengths(
            bamRanges)[as.character(seqnames(bamRanges))])) {
            offB <- TRUE
            bamRanges <- .dropOrTrimGR(bamRanges,onFlankFail)
            # mainRanges will be offended as well
            suppressWarnings(mainRanges <- 
                getMainRnaRanges(genomeRanges,flank))
            mainRanges <- .dropOrTrimGRL(mainRanges,onFlankFail)
        }
        else
            mainRanges <- getMainRnaRanges(genomeRanges,flank)
        return(list(mainRanges=mainRanges,bamRanges=bamRanges,offBound=offB))
    }
}

.dropOrTrimGR <- function(gr,opt) {
    if (opt == "drop") {
        offStart <- which(start(gr)<1)
        offEnd <- which(end(gr)>seqlengths(gr)[as.character(seqnames(gr))])
        bad <- union(offStart,offEnd)
        gr <- gr[-bad]
    }
    else if (opt == "trim")
        gr <- trim(gr)
    return(gr)
}

.dropOrTrimGRL <- function(grl,opt) {
    if (opt == "drop") {
        offStart <- which(lengths(which(start(grl)<1))>0)
        tmp <- unlist(grl)
        tf <- end(tmp)>seqlengths(tmp)[as.character(seqnames(tmp))]
        tfs <- split(tf,tmp$gene_id)
        offEnd <- which(vapply(tfs,function(x) any(x),logical(1)))
        bad <- union(offStart,offEnd)
        grl <- grl[-bad]
    }
    else if (opt == "trim")
        grl <- trim(grl)
    return(grl)
}

getMainRnaRanges <- function(genomeRanges,flank) {
    message("Creating summarized exon flanking region ",flank[1]," bps and ",
        flank[2]," bps")
    #if (is.null(rc)) {
    #    # For some progress recording...
    #    flankedSexon <- lapply(genomeRanges,flankFirstLastOld,
    #        f[1],f[2],rc=rc)
    #}
    #else
    #    flankedSexon <- cmclapply(genomeRanges,flankFirstLastOld,
    #        flank[1],flank[2],rc=rc)
    #message("Creating a GRangesList with the flanking regions")
    #return(GRangesList(flankedSexon))
    if (!is(genomeRanges,"CompressedRangesList"))
        genomeRanges <- updateObject(genomeRanges)
    return(flankFirstLast(genomeRanges,flank[1],flank[2]))
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
        utr3 = {
            w <- width(ranges)
            ranges <- promoters(ranges,upstream=flank[1],downstream=0)
            return(resize(ranges,width=w+flank[1]+flank[2]))
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
        return(cmclapply(seq_len(length(ranges)),function(i,dat,fl) {
            message("Resizing rangeset ",i," to ",fl," bases")
            return(trim(resize(dat[[i]],width=fl,fix="start")))
        },ranges,fragLen,rc=rc))
}

flankFirstLast <- function(x,u,d) {
    # Names of x for later respliting
    n <- names(x)
    f <- rep(n,lengths(x))
    
    # Collapse initial gene ranges list
    y <- unlist(x)
    
    # Find the index of first exons
    fst <- grep("_MEX_1$",names(y))
    # How many do we have in each gene
    dif <- diff(fst) - 1
    # Based on how many, find the index of last exons
    lst <- c(fst[seq_len(length(fst)-1)] + dif,length(y))
    
    # Define upstream flanking based on strand of collapsed list
    up <- rep(u,length(fst))
    up[which(as.character(strand(y[fst])) == "-")] <- d
    do <- rep(d,length(lst))
    do[which(as.character(strand(y[lst])) == "-")] <- u
    
    # Resize
    start(y[fst]) <- start(y[fst]) - up
    end(y[lst]) <- end(y[lst]) + do
    
    # Resplit
    names(y) <- as.character(y$exon_id)
    return(split(y,f))
}

cleanRanges <- function(aRange,level,org) {
    # Need to convert seqnames to character because %in% does not work...
    if (level==1) {
        chrs <- getValidChrsWithMit(org)
        aRange <- aRange[as.character(seqnames(aRange)) %in% chrs]
        newsi <- which(as.character(seqlevels(aRange)) %in% chrs)
        seqlevels(aRange) <- seqlevels(aRange)[newsi]
        seqinfo(aRange) <- seqinfo(aRange)[seqlevels(aRange)]
    }
    else if (level==2) {
        chrs <- getValidChrs(org)
        aRange <- aRange[as.character(seqnames(aRange)) %in% chrs]
        newsi <- which(as.character(seqlevels(aRange)) %in% chrs)
        seqlevels(aRange) <- seqlevels(aRange)[newsi]
        seqinfo(aRange) <- seqinfo(aRange)[seqlevels(aRange)]
    }
    else if (level==3) {
        chrs <- getValidChrs(org)
        aRange <- aRange[as.character(seqnames(aRange)) %in% chrs]
        newsi <- which(as.character(seqlevels(aRange)) %in% chrs)
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

readBamIntervals <- function(bam,gr,sa=c("keep","remove","split"),sq=0.75,
    params=NULL) {
    sa = sa[1]
    checkTextArgs("sa",sa,c("keep","remove","split"),multiarg=FALSE)
    checkNumArgs("sq",sq,"numeric",c(0,1),"botheq")
    bamIndex <- paste(bam,"bai",sep=".")
    if (is.null(params))
        params <- ScanBamParam(which=gr)
    else
        bamWhich(params) <- gr
    switch(sa,
        keep = {
            return(as(readGAlignments(file=bam,index=bamIndex,
                param=params,with.which_label=TRUE),"GRanges"))
        },
        split = {
            return(unlist(grglist(readGAlignments(file=bam,
                index=bamIndex,param=params,with.which_label=TRUE))))
        },
        remove = {
            reads <- as(readGAlignments(file=bam.file,index=bamIndex,
                param=params,with.which_label=TRUE),"GRanges")
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

splitRanges <- function(x,n,type=c("variable","fixed"),flank=NULL,where=NULL,
    rc=rc) {
    if (!is(x,"GRanges"))
        stop("x must be a GRanges object!")
    
    type <- type[1]
    if (!is.null(flank)) {
        switch(where,
            center = {
                x <- narrow(x,start=flank[1]+1,width=width(x)-flank[1]-flank[2])
            },
            upstream = {
                x <- resize(x,width=flank[1],fix="start")
            },
            downstream = {
                x <- resize(x,width=flank[2],fix="end")
            }
        )
    }
    
    # There is the possibility that some reference ranges are trimmed, where the
    # flanking regions will be slighlty smaller. Still, since binning is 
    # mandatory with the non-coverage approach, this will have practically no
    # effect to general profiles.
    message("  getting areas for profiling")
    if (type == "fixed") {
        windows <- tileRanges(x,n)
        windows <- windows[names(x)]
    }
    else if (type == "variable") { # Dynamic
        # sqrt function provides better width distribution than log and of
        # course better than natural scale. For better results we also exlcude
        # width extremes
        message("  creating dynamic bins")
        w <- width(x)
        qs <- quantile(w,c(0.01,0.99))
        w[w<=qs[1]] <- qs[1]
        w[w>=qs[2]] <- qs[2]
        cutFac <- cut(sqrt(w),breaks=n,labels=seq_len(n))
        #cutFac <- cut(w,breaks=n,labels=seq_len(n))
        #cutFac <- cut(log(w),breaks=n,labels=seq_len(n))
        cutFacFreq <- table(cutFac)
        zeroFreq <- which(cutFacFreq==0)
        if (length(zeroFreq) > 0)
            cutFacFreq <- cutFacFreq[-zeroFreq]
        binSizes <- as.numeric(names(cutFacFreq))
        windows <- cmclapply(binSizes,function(i,B,Y) {
            return(tileRanges(Y[which(B==i)],i))
        },cutFac,x,rc=rc)
        message("  getting final areas for profiling")
        windows <- do.call("c",windows)
        windows <- windows[names(x)]
    }
    return(windows)
}

splitRangesList <- function(x,n,type=c("variable","fixed"),flank=NULL,
    where=NULL,rc=NULL) {
    if (!is(x,"GRangesList"))
        stop("x must be a GRangesList object!")
    
    type <- type[1]
    # Store the order, will be needed later
    theOrder <- names(x)
    
    # lapply over large GRangesList is insanely slow
    message("  getting elements to bin")
    if (!is.null(flank)) {
        # Collapse the GRangesList
        y <- suppressWarnings(trim(unlist(x)))
        # Since the list contents are fixed recoup annotation elements, we can
        # use a dirty solution of grepping specific strings to get indices
        # corresponding to areas we want
        switch(where,
            center = {
                centerIndex <- grep("MEX_",names(y))
                x <- y[centerIndex]
            },
            upstream = {
                flankIndex <- grep("MEX_",names(y),invert=TRUE)
                ii <- seq(from=1,to=length(flankIndex),by=2)
                x <- y[flankIndex[ii]]
                theTiles <- splitRanges(x,n)
                # They have wrong names
                names(theTiles) <- 
                    vapply(strsplit(names(x),".",fixed=TRUE),function(x) { 
                        return(x[1]) 
                    },character(1))
                return(theTiles)
            },
            downstream = {
                flankIndex <- grep("MEX_",names(y),invert=TRUE)
                ii <- seq(from=2,to=length(flankIndex),by=2)
                x <- y[flankIndex[ii]]
                theTiles <- splitRanges(x,n)
                # They have wrong names
                names(theTiles) <- 
                    vapply(strsplit(names(x),".",fixed=TRUE),function(x) { 
                        return(x[1]) 
                    },character(1))
                return(theTiles)
            }
        )
    }
    
    # Resplit and reorder
    geme <- strsplit(names(x),".",fixed=TRUE)
    genes <- vapply(geme,function(x) {
        return(x[1])
    },character(1))
    names(x) <- vapply(geme,function(x) {
        return(x[2])
    },character(1))
    tmp <- split(x,genes)
    tmp <- tmp[theOrder]
    
    # Now, in the case of exons, we must bin each exon according to its length
    # and the number of bins must add to n, so we must construct a binning
    # scheme which is proportional to the exon lengths
    message("  creating proportional bin sizes")
    wList <- width(tmp) # IntegerList
    if (type == "fixed")
        w <- lapply(wList,fixedRangedBinSizes,n)
    else {
        prew <- sum(wList)
        qs <- quantile(prew,c(0.01,0.99))
        prew[prew<=qs[1]] <- qs[1]
        prew[prew>=qs[2]] <- qs[2]
        N <- as.numeric(cut(sqrt(prew),breaks=n,labels=seq_len(n)))
        w <- cmclapply(seq_len(length(N)),variableRangedBinSizes,wList,N,rc=rc)
    }
    
    # Correct the contents of tmp (GRangesList) by getting the offending genes
    # from w above through the flag list member and remove the exons that can't
    # be binned because of the allowed number of bins
    message("  veryfiying the bin sizes (may take a few minutes)")
    flagInd <- which(vapply(w,function(a) return(a$flag),numeric(1)))
    for (i in flagInd) { # How many can they be...
        message("    veryfying ",i)
        gr <- GRangesList(tmp[[i]][w[[i]]$exonsToUse])
        names(gr) <- names(tmp[i])
        tmp[i] <- gr
        # Never attempt again a usage of the '[[' operator in GRangesList...
        # This also means that we never should use lapply methods... Better
        # unlist and re-list using split and some factor like we do here.
        # https://support.bioconductor.org/p/77300/#77447
    }
    # Now unlist tmp...
    Y <- unlist(tmp)
    # ...and get the bin size vectors
    W <- unlist(lapply(w,function(a) {
        return(a$binSize)
    }),use.names=FALSE)
    names(W) <- names(Y)
    
    # Create a hash of bin sizes and then apply the tileRanges function with
    # n=member_of_table. It should runs faster and in the end combine the 
    # GRangesLists with "c" and reorder according to original names
    message("  creating the actual bins")
    hash <- table(W)
    binSizes <- as.numeric(names(hash))
    # Consider cmclapply
    theTiles <- cmclapply(binSizes,function(i,W,Y) {
        return(tileRanges(Y[which(W==i)],i))
    },W,Y,rc=rc)
    theTiles <- do.call("c",theTiles)
    theTiles <- theTiles[names(Y)]
    
    # Now unlist the tiles and recombine to a new GRangesList based on gene ids
    message("  fixing bin and bin collection names")
    theTiles <- unlist(theTiles)
    # Took a little bit more than expected, until changed the pattern to "." (we
    # don't need a regex) and fixed=TRUE
    # FIXME: The split approach won't work properly with UCSC accessions as
    # they contain dots (.). We need to figure this out.
    tileGenes <- vapply(strsplit(names(theTiles),".",fixed=TRUE),function(x) {
        return(x[1])
    },character(1))
    # We need only the gene names in the end, the tiles must be anonymous
    theTiles <- unname(theTiles)
    # ...and this should conclude
    message("  getting final areas for profiling")
    theTiles <- split(theTiles,tileGenes)
    theTiles <- theTiles[theOrder]
    
    #    len <- length(gr)
    #    if (n - len < len) # Sample w/o replacement
    #        ts <- sample(gr,n-len)
    #    else
    #        ts <- sample(gr,n-len,replace=TRUE)
    
    return(theTiles)
    
    # This is the final first version of this function. After this, for the
    # genes that have <n or >n bins, the following heuristics will be applied
    # after after calculating the overlaps, rpm etc.
    # <n: Interpolate on the measurement up to n
    # >n: Interpolate on the measurement down to n
}

tileRanges <- function(x,n) {
    # Perform the tiling operation
    # 1. Find widths of potentially resized ranges
    w <- width(x)
    # 2. Find those potentially smaller than number of bins
    offending <- which(w<n)
    rest <- which(w>=n)
    # 3. Tile the rest
    if (length(rest) > 0)
        restWindows <- tile(x[rest],n=n)
    else
        restWindows <- GRangesList()
    # 4. Tile the offending by the minimum of their width minus 1bp so as to 
    #    avoid problems with later interpolation
    if (length(offending) > 0)
        offendingWindows <- tile(x[offending],n=min(w))
    else
        offendingWindows <- GRangesList()
    # 5. Merge and restore original order
    windows <- c(restWindows,offendingWindows)
    ord <- sort(c(rest,offending),index.return=TRUE)
    windows <- windows[ord$ix]
    names(windows) <- names(x)
    return(windows)
}

fixedRangedBinSizes <- function(x,n) {
    # Firstly, if a gene has >n exons, we should downsample to n, by sorting
    # acording to length, selecting the first n longest ones and assign a 
    # unit bin size to each of them. Otherwise we end up with a lot of 
    # 0-sized bins which cannot be handled by the heuristics below. We also
    # store a flag so as to determine which indices are offending, get the
    # names and apply the coresponding indices directly to the initial 
    # GRangesList, otherwise very slow...
    if (length(x) == 1 && x < n)
        return(list(binSize=x,exonsToUse=1,flag=FALSE))
    if (length(x) > n) {
         s <- sort(x,decreasing=TRUE,index.return=TRUE)
         j <- sort(s$ix[seq_len(n)])
         return(list(binSize=rep(1,n),exonsToUse=j,flag=TRUE))
    }
    sumLen <- sum(x)
    binSize <- round(x*n/sumLen)
    if (any(binSize==0))
        binSize[binSize==0] <- 1
    check <- sum(binSize)
    if (check > n) {
        # Randomly subtract, but protect agains zero by removing from the
        # sampling space, bins of size 1
        d <- check - n
        nonOneInds <- seq_len(length(binSize))
        sizeOne <- which(binSize==1)
        if (length(sizeOne) > 0)
            nonOneInds <- which(binSize > 1)
        # nonOneInds is the protected sampling space, however it may end 
        # smaller than the number of bins to become smaller... In that case
        # we remove 2 from the remaining largest bins, unless the remaining
        # largest bins are of length 2 where in this case we have to remove
        # 1... aargh... or we simply eradicate remaining exons in this case
        # after all this as they will not be informative anyway and flag.
        lnoi <- length(nonOneInds)
        if (lnoi < d) {
            dd <- d - lnoi
            ss <- sort(binSize[nonOneInds],decreasing=TRUE,
                index.return=TRUE)
            removeSize <- rep(1,lnoi)
            removeSize[ss$ix[seq_len(dd)]] <- 2
            binSize[nonOneInds] <- binSize[nonOneInds] - removeSize
            if (any(binSize[nonOneInds]==0)) {
                ze <- which(binSize==0)
                exonsToUse <- seq_len(length(binSize))
                binSize <- binSize[-ze]
                exonsToUse <- exonsToUse[-ze]
                return(list(binSize=binSize,exonsToUse=exonsToUse,
                    flag=TRUE))
            }
        }
        else {
            ii <- sample(nonOneInds,d)
            binSize[ii] <- binSize[ii] - 1
        }
    }
    else if (check < n) {
        # Randomly add, in this case there is no risk of letting a 0-size
        # bin go, as the minimum bin size is one from above check
        d <- n - check
        ii <- sample(length(binSize),d)
        binSize[ii] <- binSize[ii] + 1
    }
    return(list(binSize=binSize,exonsToUse=seq_len(length(binSize)),flag=FALSE))
}

variableRangedBinSizes <- function(i,w,n) {
    x <- w[[i]]
    if (length(x) == 1 && x < n[i])
        return(list(binSize=x,exonsToUse=1,flag=FALSE))
    if (length(x) > n[i]) {
         s <- sort(x,decreasing=TRUE,index.return=TRUE)
         j <- sort(s$ix[seq_len(n[i])])
         return(list(binSize=rep(1,n[i]),exonsToUse=j,flag=TRUE))
    }
    sumLen <- sum(x)
    binSize <- round(x*n[i]/sumLen)
    if (any(binSize==0))
        binSize[binSize==0] <- 1
    check <- sum(binSize)
    if (check > n[i]) {
        d <- check - n[i]
        nonOneInds <- seq_len(length(binSize))
        sizeOne <- which(binSize==1)
        if (length(sizeOne) > 0)
            nonOneInds <- which(binSize > 1)
        lnoi <- length(nonOneInds)
        if (lnoi < d) {
            dd <- d - lnoi
            ss <- sort(binSize[nonOneInds],decreasing=TRUE,
                index.return=TRUE)
            removeSize <- rep(1,lnoi)
            removeSize[ss$ix[seq_len(dd)]] <- 2
            binSize[nonOneInds] <- binSize[nonOneInds] - removeSize
            if (any(binSize[nonOneInds]==0)) {
                print("WTF")
                ze <- which(binSize==0)
                exonsToUse <- seq_len(length(binSize))
                binSize <- binSize[-ze]
                exonsToUse <- exonsToUse[-ze]
                return(list(binSize=binSize,exonsToUse=exonsToUse,
                    flag=TRUE))
            }
        }
        else {
            ii <- sample(nonOneInds,d)
            binSize[ii] <- binSize[ii] - 1
        }
    }
    else if (check < n[i]) {
        d <- n[i] - check
        ii <- sample(length(binSize),d)
        binSize[ii] <- binSize[ii] + 1
    }
    return(list(binSize=binSize,exonsToUse=seq_len(length(binSize)),flag=FALSE))
}

.subsetGRangesByBamHeader <- function(gr,bam) {
    ci <- .chromInfoFromBAM(bam)
    
    if (!any(rownames(ci) %in% seqlevels(gr)))
        # Trying to subset by non-overlaping sequences??? Return empty
        return(gr[-seq_along(gr)])
    
    if (nrow(ci) >= length(seqlevels(gr)))
        # Nothing to do or cannot subset
        return(gr)
    
    # Create a new Seqinfo object for the subset
    sf <- seqinfo(gr)
    circ <- isCircular(sf)[rownames(ci)]
    gen <- genome(sf)[rownames(ci)]
    subSf <- Seqinfo(seqnames=rownames(ci),seqlengths=ci[,1],isCircular=circ,
        genome=gen)
    
    # Do the actual subsetting
    subGr <- gr[seqnames(gr) %in% rownames(ci)]
    seqlevels(subGr) <- seqlevels(subSf)
    seqinfo(subGr) <- subSf
    if (is(gr,"GRangesList")) {
        subi <- which(lengths(subGr)>0)
        if (length(subi) > 0)
            subGr <- subGr[subi]
    }
    
    return(subGr)
}

.subsetGRangesByChrs <- function(gr,chrs) {
    if (!any(chrs %in% seqlevels(gr)))
        # Trying to subset by non-overlaping sequences??? Return empty
        return(gr[-seq_along(gr)])
    
    if (length(chrs) >= length(seqlevels(gr)))
        # Nothing to do or cannot subset
        return(gr)
    
    # Create a new Seqinfo object for the subset
    sf <- seqinfo(gr)
    lens <- seqlengths(sf)[chrs]
    circ <- isCircular(sf)[chrs]
    gen <- genome(sf)[chrs]
    subSf <- Seqinfo(seqnames=chrs,seqlengths=lens,isCircular=circ,genome=gen)
    
    # Do the actual subsetting
    subGr <- gr[seqnames(gr) %in% chrs]
    seqlevels(subGr) <- seqlevels(subSf)
    seqinfo(subGr) <- subSf
    if (is(gr,"GRangesList")) {
        subi <- which(lengths(subGr)>0)
        if (length(subi) > 0)
            subGr <- subGr[subi]
    }
    
    return(subGr)
}

################################################################################

#flankFirstLastOld <- function(x,u,d) {
#    s <- as.character(strand(x[1]))
#    n <- length(x)
#    w1 <- width(x[1])
#    w2 <- width(x[n])
#    if (s=="+") {
#        x[1] <- promoters(x[1],upstream=u,downstream=0)
#        x[1] <- resize(x[1],width=w1+u)
#        x[n] <- resize(x[n],width=w2+d)
#    }
#    else if (s=="-") {
#        x[1] <- resize(x[1],width=w1+u)
#        x[n] <- promoters(x[n],upstream=d,downstream=0)
#        x[n] <- resize(x[n],width=w2+d)
#    }
#    else {
#        x[1] <- promoters(x[1],upstream=u,downstream=0)
#        x[1] <- resize(x[1],width=w1+u)
#       x[n] <- resize(x[n],width=w2+d)
#    }
#    return(x)
#}
#
