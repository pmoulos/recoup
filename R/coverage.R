coverageRef <- function(input,mainRanges,
    strandedParams=list(strand=NULL,ignoreStrand=TRUE),rc=NULL) {
    hasCoverage <- sapply(input,function(x) is.null(x$coverage))
    if (!any(hasCoverage))
        return(input)
    names(input) <- sapply(input,function(x) return(x$id))
    for (n in names(input)) {
        message("Calculating requested regions coverage for ",input[[n]]$name)
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

coverageRnaRef <- function(input,mainRanges,
    strandedParams=list(strand=NULL,ignoreStrand=TRUE),rc=NULL) {
    .Deprecated("coverageRef")
    return(coverageRef(input,mainRanges,strandedParams,rc=NULL))
    #hasCoverage <- sapply(input,function(x) is.null(x$coverage))
    #if (!any(hasCoverage))
    #    return(input)
    #names(input) <- sapply(input,function(x) return(x$id))
    #for (n in names(input)) {
    #    message("Calculating genebody exon coverage for ",input[[n]]$name)
    #    if (!is.null(input[[n]]$ranges))
    #        theRanges <- splitBySeqname(input[[n]]$ranges)
    #    else
    #        theRanges <- input[[n]]$file
    #    input[[n]]$coverage <- calcCoverage(theRanges,mainRanges,
    #        strand=strandedParams$strand,
    #        ignore.strand=strandedParams$ignoreStrand,rc=rc)
    #    names(input[[n]]$coverage) <- names(mainRanges)
    #}
    #return(input)
}

calcCoverage <- function(input,mask,strand=NULL,ignore.strand=TRUE,rc=NULL) {
    if (!is(input,"GRanges") && !is.list(input) && is.character(input)
        && !file.exists(input))
        stop("The input argument must be a GenomicRanges object or a valid ",
            "BAM/BigWig file or a list of GenomicRanges")
    if (!is(mask,"GRanges") && !is(mask,"GRangesList"))
        stop("The mask argument must be a GRanges or GRangesList object")
    isBam <- isBigWig <- FALSE
    if (is.character(input) && file.exists(input)) {
        if (length(grep("\\.bam$",input,ignore.case=TRUE,perl=TRUE))>0)
            isBam <- TRUE
        else if (length(grep("\\.(bigwig|bw|wig|bg)$",input,ignore.case=TRUE,
            perl=TRUE))>0)
            isBigWig <- TRUE
    }
    if (!is.null(strand) && !is.list(strand) && !isBam && !isBigWig) {
        message("Retrieving ",strand," reads...")
        input <- input[strand(input)==strand]
    }
    #if (!is.list(input) && !isBam && !isBigWig)
    #    input <- splitBySeqname(input)
    #if (is(mask,"GRanges"))
    #    #index <- 1:length(mask)
    #    index <- split(1:length(mask),as.character(seqnames(mask)))
    #else if (is(mask,"GRangesList"))
    #    index <- split(1:length(mask),as.character(runValue(seqnames(mask))))
    if (isBam) #{
        #index2 <- 1:length(mask)
        #cov <- cmclapply(index2,coverageFromBam,mask,input,ignore.strand,rc=rc)
        cov <- coverageFromBam(input,mask,ignore.strand,rc=rc)
    #}
    else if (isBigWig)
        #cov <- cmclapply(index,coverageFromBigWig,mask,input,rc=rc)
        #cov <- unlist(lapply(index,coverageFromBigWig,mask,input,rc=rc),
        #    use.names=FALSE)
        cov <- coverageFromBigWig(input,mask,rc=rc)
    else
        #cov <- cmclapply(index,coverageFromRanges,mask,input,ignore.strand,
        #    rc=rc)
        #cov <- unlist(lapply(index,coverageFromRanges,mask,input,
        #    ignore.strand,rc=rc))
        cov <- coverageFromRanges(input,mask,ignore.strand,rc=rc)
    
    # cov returning already named in this implementation
    ## Naming is problematic with this implementation if a chromosome is missing
    ## from one of the samples so is not filtered in the main recoup function
    ## object preparation
    ##names(cov) <- names(mask) # !Will not work in the case described above
    #theNames <- 
    #    unlist(lapply(strsplit(names(cov),"\\."),
    #        function(x) return(x[length(x)])))
    #names(cov) <- theNames
    gc(verbose=FALSE)
    return(cov) # Rle
}

# In the code below, x is a set of GRanges for each sequence, let's say the set
# of genes in chr1. Many of them overlap (depending on annotation set) so the
# final coverage has to be "normalized" by assigning to each range the coverage
# proportionally to overlap. So if 2 ranges overlap by 5 bases, then the 
# coverage over those 5 bases must be divided by 2. This produces smoother
# plots but not sure if it's completely biologically relevant. The older version
# considered each range individually but lacked on speed. If we do not normalize
# the potentially overlapping areas (e.g. TSSs) will maybe present artifacts.
# Maybe the coverage normalization should be added as an option in the future.
coverageFromRanges <- function(input,mask,ignore.strand,rc=NULL) {
    if (is(mask,"GRanges")) {
        chrs <- as.character(unique(seqnames(input)))
        message("  calculating total coverage")
        
        # If mask is coming from custom genome and not from recoup annotation,
        # then seqlengths are non-existent or do not match
        if (all(is.na(seqlengths(mask))))
            seqlengths(mask) <- seqlengths(input)
        
        preCov <- coverage(input)
        preCov <- preCov[chrs]
        maskList <- split(mask,seqnames(mask))
        normFactor <- coverage(maskList)
        normFactor <- normFactor[chrs]
        covs <- cmclapply(names(maskList),function(x,maskList,preCov,
            normFactor) {
            return(lazyRangesCoverage(x,maskList,preCov,normFactor))
            #cs <- lazyRangesCoverage(x,maskList,preCov,normFactor)
            #haveEqualLengths <- covLengthsEqual(cs)
            #ps <- profileMatrixSample(cs,flank,bp,haveEqualLengths,rc=NULL)
            #return(list(coverage=cs,profile=ps))
        },maskList,preCov,normFactor,rc=rc)
        covs <- unlist(covs)
        names(covs) <- names(mask)
        return(covs)
    }
    else if (is(mask,"GRangesList")) {
        # Coverage of split ranges (RNA-Seq)
        chrs <- as.character(unique(seqnames(input)))
        preCov <- coverage(input)
        preCov <- preCov[chrs]
        # The normalization factor, not really required with summarized exons
        #normFactor <- coverage(mask)
        #preCov <- preCov[names(normFactor)]
        #normCov <- preCov/normFactor
        #normCov[is.na(normCov)] <- 0
        #normCov[normCov==Inf] <- 1
        return(lazyRangesListCoverage(mask,preCov,chrs,rc=rc))
        #cs <- lazyRangesListCoverage(mask,preCov,chrs,rc=rc)
        #haveEqualLengths <- covLengthsEqual(cs)
        #ps <- profileMatrixSample(cs,flank,bp,haveEqualLengths,rc=rc)
        #return(list(coverage=cs,profile=ps))
    }
}

coverageFromBigWig <- function(input,mask,rc=NULL) {
    if (!requireNamespace("rtracklayer"))
        stop("R package rtracklayer is required to calculate coverage from ",
            "BigWig files!")
    preCov <- import.bw(input,selection=BigWigSelection(mask),as="RleList")
    chrs <- names(preCov)
    #
    if (all(is.na(seqlengths(mask))))
        seqlengths(mask) <- seqlengths(preCov)
    #
    if (is(mask,"GRanges")) {
        maskList <- split(mask,seqnames(mask))
        covs <- cmclapply(chrs,function(x,maskList,preCov) {
            return(lazyRangesCoverage(x,maskList,preCov))
        },maskList,preCov,rc=rc)
        covs <- unlist(covs)
        names(covs) <- names(mask)
        return(covs)
    }
    else if (is(x,"GRangesList"))
        return(lazyRangesListCoverage(mask,preCov,chrs,rc=rc))
}

coverageFromBam <- function(input,mask,ignore.strand,
    pp=list(spliceAction="keep",spliceRemoveQ=0.75),rc=NULL) {
    bamFile <- input
    if (is(mask,"GRanges")) {
        chrs <- as.character(unique(seqnames(mask)))
        maskList <- split(mask,seqnames(mask))
        covs <- cmclapply(chrs,function(x,maskList) {
            message("  Reading sequence ",x)
            M <- maskList[[x]]
            if (!is.null(M)) {
                reads <- readBamIntervals(bamFile,M,pp$spliceAction,
                    pp$spliceRemoveQ)
                if (length(reads)>0) {
                    seqlevels(reads) <- as.character(seqnames(M)[1])
                    tmp <- coverage(reads)
                    tmp <- tmp[[x]]
                    if (!is.null(cov)) {
                        V <- Views(tmp,ranges(M))
                        cot <- unlist(viewApply(V,function(x) x))
                        names(cot) <- names(M)
                        inv <- which(strand(M)=="-")
                        cot[inv] <- lapply(cot[inv],rev)
                        return(cot)
                    }
                }
                else {
                    message("    Sequence ",x," not found!")
                    return(Rle(NA))
                }
            }
        },maskList,rc=rc)
        covs <- unlist(covs)
        names(covs) <- names(mask)
        return(covs)
    }
    else if (is(mask,"GRangesList")) {
        Mraw <- trim(unlist(mask))
        chrs <- as.character(unique(seqnames(Mraw)))
        maskList <- split(Mraw,seqnames(Mraw))
        maskList <- maskList[chrs]
        covs <- lapply(chrs,function(x,maskList) {
            message("  Reading sequence ",x)
            M <- maskList[[x]]
            if (!is.null(M)) {
                reads <- readBamIntervals(bamFile,M,pp$spliceAction,
                    pp$spliceRemoveQ)
                if (length(reads)>0) {
                    seqlevels(reads) <- as.character(seqnames(M)[1])
                    tmp <- coverage(reads)
                    tmp <- tmp[[x]]
                    if (!is.null(cov)) {
                        V <- Views(tmp,ranges(M))
                        cot <- unlist(viewApply(V,function(x) x))
                        names(cot) <- names(M)
                        inv <- which(strand(M)=="-")
                        cot[inv] <- lapply(cot[inv],rev)
                        genes <- sapply(strsplit(names(cot),".",fixed=TRUE),
                            function(x) { return(x[2]) })
                        names(cot) <- genes
                        ugenes <- unique(genes)
                        cot <- cmclapply(ugenes,function(x,cot) {
                            return(Reduce("c",cot[names(cot)==x]))
                        },cot,rc=rc)
                        names(cot) <- ugenes
                        return(cot)
                    }
                }
                else {
                    message("    Sequence ",x," not found!")
                    return(Rle(NA))
                }
            }
        })
        covs <- unlist(covs)
        names(covs) <- names(mask)
        return(covs)
    }
}

lazyRangesCoverage <- function(x,maskList,preCov,normFactor=NULL) {
    message("  processing ",x)
    m <- maskList[[x]]
    pre <- preCov[[x]]
    if (!is.null(m) && !is.null(pre)) { # Sanity...
        co <- pre
        if (!is.null(normFactor)) {
            co <- pre/normFactor[[x]]
            co[is.na(co)] <- 0
            co[co==Inf] <- 1
        }
        V <- Views(co,ranges(m))
        cot <- unlist(viewApply(V,function(x) x))
        names(cot) <- names(m)
        
        #inv <- which(strand(m)=="-")
        inv <- which(as.logical(strand(m)=="-"))
        if (length(inv) > 0)
            cot[inv] <- lapply(cot[inv],rev)
        return(cot)
    }
    else {
        message("    ",x," not found!")
        return(Rle(NA))
    }
}

lazyRangesListCoverage <- function(mask,preCov,chrs,rc=NULL) {
    # Collapse the summarized exons structure
    Mraw <- trim(unlist(mask))
    # ...and split per chromosome
    maskList <- split(Mraw,seqnames(Mraw))
    maskList <- maskList[chrs]
    # ...finally construct Views like the simple ChIP-Seq case
    # Profiling showed that the RNA-Seq case does not benefit from 
    # parallelization of the outer loop as in the ChIP-Seq case
    V <- Views(preCov,ranges(maskList))
    # Coverage and handling strand
    covs <- unlist(viewApply(V,function(x) x))
    Msimple <- unlist(maskList)
    #inv <- which(strand(Msimple)=="-")
    inv <- which(as.logical(strand(Msimple)=="-"))
    if (length(inv) > 0)
        covs[inv] <- cmclapply(covs[inv],rev,rc=rc)
    # Now we must reconstruct the initial summarized exon structure
    names(covs) <- names(Msimple)
    # Get unique gene names from composite Rle list names
    genes <- sapply(strsplit(names(covs),".",fixed=TRUE),function(x) {
        return(x[2])
    })
    names(covs) <- genes
    ugenes <- unique(genes)
    # Finally, reconstruct a gene coverage Rle
    covs <- cmclapply(ugenes,function(x,covs) {
        return(Reduce("c",covs[names(covs)==x]))
    },covs,rc=rc)
    names(covs) <- ugenes
    return(covs)
}

################################################################################

#~ # Legacy functions

#~ coverageFromRangesOld <- function(i,mask,input,ignore.strand,rc=NULL) {
#~     x <- mask[i]
#~     if (is(x,"GRanges")) {
#~         y<-list(
#~             chromosome=as.character(seqnames(x))[1], 
#~             start=start(x),
#~             end=end(x),
#~             strand=as.character(strand(x)),
#~             reads=NULL,
#~             coverage=NULL
#~         )
#~         message("  processing ",y$chromosome)
#~         if (!is.null(input[[y$chromosome]])) {
#~             #S <- unique(subjectHits(findOverlaps(x,input[[y$chromosome]],
#~             #    ignore.strand=ignore.strand)))
#~             #y$reads <- input[[y$chromosome]][S]
#~             y$reads <- input[[y$chromosome]][
#~                 unique(subjectHits(findOverlaps(x,input[[y$chromosome]],
#~                     ignore.strand=ignore.strand)))]
#~         }
#~         else {
#~             message(y$chromosome," not found!")
#~             y$reads <- NULL
#~         }
#~         if (length(y$reads)>0) {
#~             xCov <- coverage(x)
#~             cc <- as.character(seqnames(y$reads))[1]
#~             #tt <- table(S)
#~             #w <- rep(1,length(S))
#~             #for (ii in unique(tt)) 
#~             #    w[tt==ii] <- w[tt==ii]/ii
#~             #y$coverage <- coverage(y$reads,weight=w)
#~             y$coverage <- coverage(y$reads)
#~             y$coverage <- y$coverage[[cc]]
#~             xCov <- xCov[[cc]]
#~             covs <- cmclapply(1:length(y$start),function(j,s,e,d,r) {
#~                 tryCatch({
#~                     if (d[j] == "+")
#~                         return(r[s[j]:e[j]]/xCov[s[j]:e[j]])
#~                     else if (d[j] == "-")
#~                         return(rev(r[s[j]:e[j]]/xCov[s[j]:e[j]]))
#~                     else
#~                         return(r[s[j]:e[j]]/xCov[s[j]:e[j]])
#~                 },
#~                 error=function(e) {
#~                     message("Caught invalid genomic area!")
#~                     print(mask[i][j])
#~                     message("Will return zero coverage")
#~                     return(Rle(NA))
#~                     #return(NA)
#~                 },finally={})
#~             },y$start,y$end,y$strand,y$coverage,rc=rc)
#~             names(covs) <- names(x)
#~             return(covs)
#~         }
#~         else
#~             return(Rle(NA))
#~             #return(NA)
#~     }
#~     else if (is(x,"GRangesList")) {
#~         y<-list(
#~             # Chrom and strand should always be the same from higher 
#~             # functions
#~             chromosome=as.character(seqnames(x)[[1]])[1],
#~             start=start(x), # Integer list
#~             end=end(x), # Integer list
#~             strand=as.character(strand(x)[[1]])[1],
#~             reads=NULL,
#~             coverage=NULL
#~         )
#~         message("  processing ",y$chromosome)
#~         if (!is.null(input[[y$chromosome]])) {
#~             y$reads <- input[[y$chromosome]][
#~                 subjectHits(findOverlaps(x,input[[y$chromosome]],
#~                     ignore.strand=ignore.strand))]
#~         }
#~         else {
#~             message(y$chromosome," not found!")
#~             y$reads <- NULL
#~         }
#~         if (length(y$reads)>0) {
#~             #xCov <- coverage(x)
#~             cc <- as.character(seqnames(y$reads))[1]
#~             y$coverage <- coverage(y$reads)
#~             y$coverage <- y$coverage[[cc]]
#~             #xCov <- xCov[[cc]]
#~             covs <- cmclapply(1:length(y$start),function(j,s,e,d,r) {
#~                 tryCatch({
#~                     subinds <- unlist(lapply(1:length(s[[j]]),
#~                     function(k,ss,ee) {
#~                         return(ss[[k]]:ee[[k]])
#~                     },s[[j]],e[[j]]))
#~                     if (d == "+")
#~                         #return(r[subinds]/xCov[subinds])
#~                         return(r[subinds])
#~                     else if (d == "-")
#~                         #return(rev(r[subinds]/xCov[subinds]))
#~                         return(rev(r[subinds]))
#~                     else
#~                         #return(r[subinds]/xCov[subinds])
#~                         return(r[subinds])
#~                 },error=function(e) {
#~                         print(e)
#~                         message("Caught invalid genomic area!")
#~                         print(mask[i][j])
#~                         message("Will return zero coverage")
#~                         return(Rle(NA))
#~                     },finally={})
#~                 },y$start,y$end,y$strand,y$coverage,rc=NULL)
#~             names(covs) <- names(x)
#~             return(covs)
#~         }
#~         else
#~             return(Rle(NA))
#~     }
#~ }

#~ coverageFromRangesOlder <- function(i,mask,input,ignore.strand) {
#~     if (is(mask,"GRangesList"))
#~         x <- mask[[i]]
#~     else
#~         x <- mask[i]
#~     y<-list(
#~         chromosome=as.character(seqnames(x))[1],
#~         start=start(x),
#~         end=end(x),
#~         strand=as.character(strand(x))[1],
#~         reads=NULL,
#~         coverage=NULL
#~     )
#~     if (!is.null(input[[y$chromosome]])) {
#~         y$reads <- input[[y$chromosome]][
#~             unique(subjectHits(findOverlaps(x,input[[y$chromosome]],
#~                 ignore.strand=ignore.strand)))]
#~     }
#~     else {
#~         message(y$chromosome," not found!")
#~         y$reads <- NULL
#~     }
#~     if (length(y$reads)>0) {
#~         tryCatch({
#~             cc <- as.character(seqnames(y$reads))[1]
#~             y$coverage <- coverage(y$reads)
#~             if (length(y$start)>1) { # List of exons, RNA, merge exons
#~                 i2k <- unlist(lapply(1:length(y$start),function(j,s,e) {
#~                     return(s[j]:e[j])
#~                 },y$start,y$end))
#~                 y$coverage <- y$coverage[[cc]][i2k]
#~             }
#~             else
#~                 y$coverage <- y$coverage[[cc]][y$start:y$end]
#~             if (y$strand=="+")
#~                 return(y$coverage)
#~             else if (y$strand=="-")
#~                 return(rev(y$coverage))
#~             else
#~                 return(y$coverage)
#~         },
#~         error=function(e) {
#~             message("Caught invalid genomic area!")
#~             print(mask[i])
#~             message("Will return zero coverage")
#~             return(NULL)
#~         },finally={})
#~     }
#~     else
#~         return(NULL)
#~ }

#~ coverageFromBigWigOld <- function(i,mask,input,rc=NULL) {
#~     if (!requireNamespace("rtracklayer"))
#~         stop("R package rtracklayer is required to calculate coverage from ",
#~             "BigWig files!")
#~     x <- mask[i]
#~     if (is(mask,"GRanges")) {
#~         chr <- as.character(seqnames(x))[1]
#~         bwrle <- import.bw(input,selection=BigWigSelection(x),as="RleList")
#~         if (chr %in% names(bwrle)) {
#~             covs <- cmclapply(x,function(y) {
#~                 tryCatch(
#~                     return(bwrle[[chr]][start(y):end(y)]),
#~                     error=function(e) {
#~                         message("Caught invalid genomic area!")
#~                         print(y)
#~                         message("Will return zero coverage")
#~                         return(Rle(NA))
#~                     },
#~                     finally={}
#~                 )
#~             },rc=rc)
#~             names(covs) <- names(x)
#~             return(covs)
#~         }
#~         else {
#~             message(chr,"not found!")
#~             return(Rle(NA))
#~         }
#~     }
#~     else if (is(x,"GRangesList")) {
#~         chr <- as.character(seqnames(x)[[1]])[1]
#~         message("  processing ",chr)
#~         covs <- cmclapply(x,function(y) {
#~             bwrle <- import.bw(input,selection=BigWigSelection(y),
#~              as="RleList")
#~             if (chr %in% names(bwrle)) {
#~                 tmp <- bwrle[[chr]]
#~                 i2k <- unlist(lapply(1:length(start(y)),function(j,s,e) {
#~                     return(s[j]:e[j])
#~                 },start(y),end(y)))
#~                 tryCatch(
#~                     return(tmp[i2k]),
#~                     error=function(e) {
#~                         message("Caught invalid genomic area!")
#~                         print(y)
#~                         message("Will return zero coverage")
#~                         return(Rle(NA))
#~                     },
#~                     finally={}
#~                 )
#~             }
#~             else {
#~                 message(chr,"not found!")
#~                 return(Rle(NA))
#~             }
#~         },rc=rc)
#~         names(covs) <- names(x)
#~         return(covs)
#~     }
#~ }

#~ coverageFromBigWigOlder <- function(i,mask,input) {
#~     if (!requireNamespace("rtracklayer")) {
#~         stop("R package rtracklayer is required to calculate coverage from ",
#~             "BigWig files!")
#~     }
#~     if (is(mask,"GRangesList"))
#~         x <- mask[[i]]
#~     else
#~         x <- mask[i]
#~     chr <- as.character(seqnames(x))[1]
#~     bwrle <- import.bw(input,selection=BigWigSelection(x),as="RleList")
#~     bwcov <- NULL
#~     if (chr %in% names(bwrle)) {
#~         bwcov <- tryCatch(
#~             bwrle[[chr]][start(x):end(x)],
#~             error=function(e) {
#~                 message("Caught invalid genomic area!")
#~                 print(mask[i])
#~                 message("Will return zero coverage")
#~                 return(NULL)
#~             },
#~             finally={}
#~         )
#~         return(bwcov)
#~     }
#~     else {
#~         message(chr,"not found!")
#~         return(NULL)
#~     }
#~ }

#~ coverageFromBamOld <- function(i,mask,input,ignore.strand,
#~     pp=list(spliceAction="keep")) {
#~     if (is(mask,"GRangesList"))
#~         x <- mask[[i]]
#~     else
#~         x <- mask[i]
#~     y<-list(
#~         chromosome=as.character(seqnames(x)),
#~         start=start(x),
#~         end=end(x),
#~         strand=as.character(strand(x)),
#~         reads=NULL,
#~         coverage=NULL
#~     )
#~     bam.file <- input
#~     bam.index <- paste(input,"bai",sep=".")
#~     bp <- ScanBamParam(which=x)
    
#~     switch(pp$spliceAction,
#~         keep = {
#~             y$reads <- as(readGAlignments(file=bam.file,index=bam.index,
#~                 param=bp,with.which_label=TRUE),"GRanges")
#~         },
#~         split = {
#~             y$reads <- unlist(grglist(readGAlignments(file=bam.file,
#~                 index=bam.index,param=bp,with.which_label=TRUE)))
#~         },
#~         remove = {
#~             y$reads <- as(readGAlignments(file=bam.file,index=bam.index,
#~                 param=bp,with.which_label=TRUE),"GRanges")
#~             qu <- quantile(width(y$reads),pp$spliceRemoveQ)
#~             rem <- which(width(y$reads)>qu)
#~             if (length(rem)>0)
#~                 y$reads <- y$reads[-rem]
#~             message("  Excluded ",length(rem)," reads")
#~         }
#~     )
    
#~     if (length(y$reads)>0) {
#~         tryCatch({
#~             seqlevels(y$reads) <- as.character(seqnames(x))
#~             cc <- as.character(seqnames(y$reads))[1]
#~             y$coverage <- coverage(y$reads)
#~             if (length(y$start)>1) { # List of exons, RNA, merge exons
#~                 i2k <- unlist(lapply(1:length(y$start),function(j,s,e) {
#~                     return(s[j]:e[j])
#~                 },y$start,y$end))
#~                 y$coverage <- y$coverage[[cc]][i2k]
#~             }
#~             else
#~                 y$coverage <- y$coverage[[cc]][y$start:y$end]
#~             if (y$strand=="+")
#~                 return(y$coverage)
#~             else if (y$strand=="-")
#~                 return(rev(y$coverage))
#~             else
#~                 return(y$coverage)
#~         },
#~         error=function(e) {
#~             message("Caught invalid genomic area!")
#~             print(mask[i])
#~             message("Will return zero coverage")
#~             return(NULL)
#~         },finally={})
#~     }
#~     else
#~         return(NULL)
#~ }
