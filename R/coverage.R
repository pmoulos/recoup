coverageRef <- function(input,mainRanges,
    strandedParams=list(strand=NULL,ignoreStrand=TRUE),rc=NULL) {
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
    hasCoverage <- sapply(input,function(x) is.null(x$coverage))
    if (!any(hasCoverage))
        return(input)
    names(input) <- sapply(input,function(x) return(x$id))
    for (n in names(input)) {
        message("Calculating genebody exon coverage for ",input[[n]]$name)
        if (!is.null(input[[n]]$ranges))
            theRanges <- splitBySeqname(input[[n]]$ranges)
        else
            theRanges <- input[[n]]$file
        input[[n]]$coverage <- calcCoverage(theRanges,mainRanges,
            strand=strandedParams$strand,
            ignore.strand=strandedParams$ignoreStrand,rc=rc)
        names(input[[n]]$coverage) <- names(mainRanges)
    }
    return(input)
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
    if (!is.list(input) && !isBam && !isBigWig)
        input <- splitBySeqname(input)
    if (is(mask,"GRanges"))
        #index <- 1:length(mask)
        index <- split(1:length(mask),as.character(seqnames(mask)))
    else if (is(mask,"GRangesList"))
        index <- split(1:length(mask),as.character(runValue(seqnames(mask))))
    if (isBam) {
        index2 <- 1:length(mask)
        cov <- cmclapply(index2,coverageFromBam,mask,input,ignore.strand,rc=rc)
    }
    else if (isBigWig)
        #cov <- cmclapply(index,coverageFromBigWig,mask,input,rc=rc)
        cov <- unlist(lapply(index,coverageFromBigWig,mask,input,rc=rc),
            use.names=FALSE)
    else
        #cov <- cmclapply(index,coverageFromRanges,mask,input,ignore.strand,
        #    rc=rc)
        cov <- unlist(lapply(index,coverageFromRanges,mask,input,
            ignore.strand,rc=rc))
    
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
    
    # Naming is problematic with this implementation if a chromosome is missing
    # from one of the samples so is not filtered in the main recoup function
    # object preparation
    #names(cov) <- names(mask) # !Will not work in the case described above
    theNames <- 
        unlist(lapply(strsplit(names(cov),"\\."),
            function(x) return(x[length(x)])))
    names(cov) <- theNames
    gc(verbose=FALSE)
    return(cov) # Rle
}

coverageFromRanges <- function(i,mask,input,ignore.strand,rc=NULL) {
    x <- mask[i]
    if (is(x,"GRanges")) {
        y<-list(
            chromosome=as.character(seqnames(x))[1], 
            start=start(x),
            end=end(x),
            strand=as.character(strand(x)),
            reads=NULL,
            coverage=NULL
        )
        message("  processing ",y$chromosome)
        if (!is.null(input[[y$chromosome]])) {
            #S <- unique(subjectHits(findOverlaps(x,input[[y$chromosome]],
            #    ignore.strand=ignore.strand)))
            #y$reads <- input[[y$chromosome]][S]
            y$reads <- input[[y$chromosome]][
                unique(subjectHits(findOverlaps(x,input[[y$chromosome]],
                    ignore.strand=ignore.strand)))]
        }
        else {
            message(y$chromosome," not found!")
            y$reads <- NULL
        }
        if (length(y$reads)>0) {
            xCov <- coverage(x)
            cc <- as.character(seqnames(y$reads))[1]
            #tt <- table(S)
            #w <- rep(1,length(S))
            #for (ii in unique(tt)) 
            #    w[tt==ii] <- w[tt==ii]/ii
            #y$coverage <- coverage(y$reads,weight=w)
            y$coverage <- coverage(y$reads)
            y$coverage <- y$coverage[[cc]]
            xCov <- xCov[[cc]]
            covs <- cmclapply(1:length(y$start),function(j,s,e,d,r) {
                tryCatch({
                    if (d[j] == "+")
                        return(r[s[j]:e[j]]/xCov[s[j]:e[j]])
                    else if (d[j] == "-")
                        return(rev(r[s[j]:e[j]]/xCov[s[j]:e[j]]))
                    else
                        return(r[s[j]:e[j]]/xCov[s[j]:e[j]])
                },
                error=function(e) {
                    message("Caught invalid genomic area!")
                    print(mask[i][j])
                    message("Will return zero coverage")
                    return(Rle(NA))
                    #return(NA)
                },finally={})
            },y$start,y$end,y$strand,y$coverage,rc=rc)
            names(covs) <- names(x)
            return(covs)
        }
        else
            return(Rle(NA))
            #return(NA)
    }
    else if (is(x,"GRangesList")) {
        y<-list(
            # Chrom and strand should always be the same from higher functions
            chromosome=as.character(seqnames(x)[[1]])[1],
            start=start(x), # Integer list
            end=end(x), # Integer list
            strand=as.character(strand(x)[[1]])[1],
            reads=NULL,
            coverage=NULL
        )
        message("  processing ",y$chromosome)
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
            #xCov <- coverage(x)
            cc <- as.character(seqnames(y$reads))[1]
            y$coverage <- coverage(y$reads)
            y$coverage <- y$coverage[[cc]]
            #xCov <- xCov[[cc]]
            covs <- cmclapply(1:length(y$start),function(j,s,e,d,r) {
                tryCatch({
                    subinds <- unlist(lapply(1:length(s[[j]]),function(k,ss,ee) {
                        return(ss[[k]]:ee[[k]])
                    },s[[j]],e[[j]]))
                    if (d == "+")
                        #return(r[subinds]/xCov[subinds])
                        return(r[subinds])
                    else if (d == "-")
                        #return(rev(r[subinds]/xCov[subinds]))
                        return(rev(r[subinds]))
                    else
                        #return(r[subinds]/xCov[subinds])
                        return(r[subinds])
                },error=function(e) {
                        print(e)
                        message("Caught invalid genomic area!")
                        print(mask[i][j])
                        message("Will return zero coverage")
                        return(Rle(NA))
                    },finally={})
                },y$start,y$end,y$strand,y$coverage,rc=NULL)
            names(covs) <- names(x)
            return(covs)
        }
        else
            return(Rle(NA))
    }
}

coverageFromRangesOld <- function(i,mask,input,ignore.strand) {
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
            unique(subjectHits(findOverlaps(x,input[[y$chromosome]],
                ignore.strand=ignore.strand)))]
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

coverageFromBigWig <- function(i,mask,input,rc=NULL) {
    if (!requireNamespace("rtracklayer"))
        stop("R package rtracklayer is required to calculate coverage from ",
            "BigWig files!")
    x <- mask[i]
    if (is(mask,"GRanges")) {
        chr <- as.character(seqnames(x))[1]
        bwrle <- import.bw(input,selection=BigWigSelection(x),as="RleList")
        if (chr %in% names(bwrle)) {
            covs <- cmclapply(x,function(y) {
                tryCatch(
                    return(bwrle[[chr]][start(y):end(y)]),
                    error=function(e) {
                        message("Caught invalid genomic area!")
                        print(y)
                        message("Will return zero coverage")
                        return(Rle(NA))
                    },
                    finally={}
                )
            },rc=rc)
            names(covs) <- names(x)
            return(covs)
        }
        else {
            message(chr,"not found!")
            return(Rle(NA))
        }
    }
    else if (is(x,"GRangesList")) {
        chr <- as.character(seqnames(x)[[1]])[1]
        message("  processing ",chr)
        covs <- cmclapply(x,function(y) {
            bwrle <- import.bw(input,selection=BigWigSelection(y),as="RleList")
            if (chr %in% names(bwrle)) {
                tmp <- bwrle[[chr]]
                i2k <- unlist(lapply(1:length(start(y)),function(j,s,e) {
                    return(s[j]:e[j])
                },start(y),end(y)))
                tryCatch(
                    return(tmp[i2k]),
                    error=function(e) {
                        message("Caught invalid genomic area!")
                        print(y)
                        message("Will return zero coverage")
                        return(Rle(NA))
                    },
                    finally={}
                )
            }
            else {
                message(chr,"not found!")
                return(Rle(NA))
            }
        },rc=rc)
        names(covs) <- names(x)
        return(covs)
    }
}

coverageFromBigWigOld <- function(i,mask,input) {
    if (!requireNamespace("rtracklayer")) {
        stop("R package rtracklayer is required to calculate coverage from ",
            "BigWig files!")
    }
    if (is(mask,"GRangesList"))
        x <- mask[[i]]
    else
        x <- mask[i]
    chr <- as.character(seqnames(x))[1]
    bwrle <- import.bw(input,selection=BigWigSelection(x),as="RleList")
    bwcov <- NULL
    if (chr %in% names(bwrle)) {
        bwcov <- tryCatch(
            bwrle[[chr]][start(x):end(x)],
            error=function(e) {
                message("Caught invalid genomic area!")
                print(mask[i])
                message("Will return zero coverage")
                return(NULL)
            },
            finally={}
        )
        return(bwcov)
    }
    else {
        message(chr,"not found!")
        return(NULL)
    }
}

coverageFromBam <- function(i,mask,input,ignore.strand,
    pp=list(spliceAction="keep")) {
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

#coverageFromRanges <- function(i,mask,input,ignore.strand,flank,binParams,
#   rc=NULL) {
#    x <- mask[i]
#    if (is(x,"GRanges")) {
#        y<-list(
#            chromosome=as.character(seqnames(x))[1], 
#            start=start(x),
#            end=end(x),
#            strand=as.character(strand(x)),
#            reads=NULL,
#            coverage=NULL
#        )
#        message("  processing ",y$chromosome)
#        if (!is.null(input[[y$chromosome]])) {
#            #S <- unique(subjectHits(findOverlaps(x,input[[y$chromosome]],
#            #    ignore.strand=ignore.strand)))
#            #y$reads <- input[[y$chromosome]][S]
#            y$reads <- input[[y$chromosome]][
#                unique(subjectHits(findOverlaps(x,input[[y$chromosome]],
#                    ignore.strand=ignore.strand)))]
#        }
#        else {
#            message(y$chromosome," not found!")
#            y$reads <- NULL
#        }
#        if (length(y$reads)>0) {
#            xCov <- coverage(x)
#            cc <- as.character(seqnames(y$reads))[1]
#            #tt <- table(S)
#            #w <- rep(1,length(S))
#            #for (ii in unique(tt)) 
#            #    w[tt==ii] <- w[tt==ii]/ii
#            #y$coverage <- coverage(y$reads,weight=w)
#            y$coverage <- coverage(y$reads)
#            y$coverage <- y$coverage[[cc]]
#            xCov <- xCov[[cc]]
#            covs <- cmclapply(1:length(y$start),function(j,s,e,d,r) {
#                tryCatch({
#                    if (d[j] == "+")
#                        return(r[s[j]:e[j]]/xCov[s[j]:e[j]])
#                    else if (d[j] == "-")
#                        return(rev(r[s[j]:e[j]]/xCov[s[j]:e[j]]))
#                    else
#                        return(r[s[j]:e[j]]/xCov[s[j]:e[j]])
#                },
#               error=function(e) {
#                    message("Caught invalid genomic area!")
#                    print(mask[i][j])
#                    message("Will return zero coverage")
#                    return(Rle(NA))
#                    #return(NA)
#                },finally={})
#            },y$start,y$end,y$strand,y$coverage,rc=rc)
#            profs <- profileMatrixSeg(covs,flank,binParams,rc=rc)
#            return(list(covs=covs,profs=profs))
#        }
#        else
#            return(Rle(NA))
#            #return(NA)
#    }
#    else if (is(x,"GRangesList")) {
#        y<-list(
#            # Chrom and strand should always be the same from higher functions
#            chromosome=as.character(seqnames(x)[[1]])[1],
#            start=start(x), # Integer list
#            end=end(x), # Integer list
#            strand=as.character(strand(x)[[1]])[1],
#            reads=NULL,
#            coverage=NULL
#        )
#        message("  processing ",y$chromosome)
#        if (!is.null(input[[y$chromosome]])) {
#            y$reads <- input[[y$chromosome]][
#                subjectHits(findOverlaps(x,input[[y$chromosome]],
#                    ignore.strand=ignore.strand))]
#        }
#        else {
#            message(y$chromosome," not found!")
#            y$reads <- NULL
#        }
#        if (length(y$reads)>0) {
#            #xCov <- coverage(x)
#            cc <- as.character(seqnames(y$reads))[1]
#            y$coverage <- coverage(y$reads)
#            y$coverage <- y$coverage[[cc]]
#            #xCov <- xCov[[cc]]
#            covs <- cmclapply(1:length(y$start),function(j,s,e,d,r) {
#                tryCatch({
#                    subinds <- unlist(lapply(1:length(s[[j]]),function(k,ss,ee) {
#                        return(ss[[k]]:ee[[k]])
#                    },s[[j]],e[[j]]))
#                    if (d == "+")
#                        #return(r[subinds]/xCov[subinds])
#                        return(r[subinds])
#                    else if (d == "-")
#                        #return(rev(r[subinds]/xCov[subinds]))
#                        return(rev(r[subinds]))
#                    else
#                        #return(r[subinds]/xCov[subinds])
#                        return(r[subinds])
#                },error=function(e) {
#                        print(e)
#                        message("Caught invalid genomic area!")
#                        print(mask[i][j])
#                        message("Will return zero coverage")
#                        return(Rle(NA))
#                },finally={})
#            },y$start,y$end,y$strand,y$coverage,rc=NULL)
#            profs <- profileMatrixSeg(covs,flank,binParams,rc=rc)
#            return(list(covs=covs,profs=profs))
#        }
#        else
#            return(Rle(NA))
#    }
#}
