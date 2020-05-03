rpMatrix <- function(input,mainRanges,flank,binParams,
    strandedParams=list(strand=NULL,ignoreStrand=TRUE),rc=NULL) {
    hasProfile <- sapply(input,function(x) is.null(x$profile))
    if (!any(hasProfile))
        return(input)
    
    # Nullify the coverage slot
    input <- lapply(input,function(x) {
        x$coverage <- NULL
        return(x)
    })
    
    haveEqualLengths <- FALSE
    if (is(mainRanges,"GRanges")) # May be from equal length regions
        haveEqualLengths <- covLengthsEqual(mainRanges)
        
    # Keep a log of which entities are reverse stranded (if any) for correcting
    # the final profile matrix (if required!)
    revInd <- NULL
    if (is(mainRanges,"GRanges"))
        #revInd <- which(strand(mainRanges) == "-")
        revInd <- which(as.logical(strand(mainRanges) == "-"))
    else if (is(mainRanges,"GRangesList")) {
        s <- unlist(runValue(strand(mainRanges)))
        revInd <- which(s=="-")
    }
    
    logS <- TRUE
    
    message("Binning reference regions to create profile")
    if (!haveEqualLengths) {
        # First get the areas to calculate overlaps ONCE
        message(" center")
        if (is(mainRanges,"GRanges"))
            areasCenter <- splitRanges(mainRanges,binParams$regionBinSize,
                binParams$binType,flank,where="center",rc=rc)
        else if (is(mainRanges,"GRangesList"))
            areasCenter <- splitRangesList(mainRanges,binParams$regionBinSize,
                binParams$binType,flank,where="center",#seed=binParams$seed,
                rc=rc)
        
        if (flank[1]==0)
            areasUpstream <- NULL
        else {
            message(" upstream")
            if (is(mainRanges,"GRanges"))
                areasUpstream <- splitRanges(mainRanges,binParams$flankBinSize,
                    binParams$binType,flank,where="upstream",rc=rc)
            else if (is(mainRanges,"GRangesList"))
                areasUpstream <- splitRangesList(mainRanges,
                    binParams$flankBinSize,binParams$binType,flank,
                    #where="upstream",seed=binParams$seed,rc=rc)
                    where="upstream",rc=rc)
        }
        if (flank[2]==0)
            areasDownstream <- NULL
        else {
            message(" downstream")
            if (is(mainRanges,"GRanges"))
                areasDownstream <- splitRanges(mainRanges,
                    binParams$flankBinSize,binParams$binType,flank,
                    where="downstream",rc=rc)
            else if (is(mainRanges,"GRangesList"))
                areasDownstream <- splitRangesList(mainRanges,
                    binParams$flankBinSize,binParams$binType,flank,
                    #where="downstream",seed=binParams$seed,rc=rc)
                    where="downstream",rc=rc)
        }
        
        # The calculate the profiles
        names(input) <- sapply(input,function(x) return(x$id))
        for (n in names(input)) {
            message("Calculating requested regions rpm for ",input[[n]]$name)
            if (!is.null(input[[n]]$ranges)) {
                message(" center")
                center <- rpMatrixSample(input[[n]]$ranges,areasCenter,
                    binParams$regionBinSize,binParams$interpolation,
                    strandedParams$strand,strandedParams$ignoreStrand,
                    #logScale=logS,binParams$seed,rc=rc)
                    logScale=logS,rc=rc)
                
                r <- flank/sum(flank)
                message(" upstream")
                if (is.null(areasUpstream))
                    left <- NULL
                else
                    left <- rpMatrixSample(input[[n]]$ranges,areasUpstream,
                        round(2*binParams$flankBinSize*r[1]),
                        binParams$interpolation,strandedParams$strand,
                        strandedParams$ignoreStrand,logScale=logS,rc=rc)
                        #binParams$seed,rc=rc)
                
                message(" downstream")
                if (is.null(areasDownstream))
                    right <- NULL
                right <- rpMatrixSample(input[[n]]$ranges,areasDownstream,
                    round(2*binParams$flankBinSize*r[2]),
                    binParams$interpolation,strandedParams$strand,
                    strandedParams$ignoreStrand,logScale=logS,#binParams$seed,
                    rc=rc)
                
                input[[n]]$profile <- cbind(left,center,right)
                if (length(revInd) > 0)
                    input[[n]]$profile[revInd,] <- 
                        input[[n]]$profile[revInd,ncol(input[[n]]$profile):1]
            }
            else {
                warning("rpm signal type not yet supported for BAM file ",
                    "streaming!",immediate.=TRUE)
                input[[n]]$profile <- NULL
            }
        }
    }
    else {
        if (is(mainRanges,"GRanges"))
            areas <- splitRanges(mainRanges,binParams$regionBinSize,rc=rc)
        else if (is(mainRanges,"GRangesList"))
            areas <- splitRangesList(mainRanges,binParams$regionBinSize,
                binParams$binType,flank,
                #where="center",seed=binParams$seed,rc=rc)
                where="center",rc=rc)
                
        # The calculate the profiles
        names(input) <- sapply(input,function(x) return(x$id))
        for (n in names(input)) {
            message("Calculating requested regions for ",input[[n]]$name)
            if (!is.null(input[[n]]$ranges))
                input[[n]]$profile <- rpMatrixSample(input[[n]]$ranges,
                    areas,binParams$regionBinSize,binParams$interpolation,
                    strandedParams$strand,strandedParams$ignoreStrand,
                    #logScale=logS,binParams$seed,rc=rc)
                    logScale=logS,rc=rc)
                if (length(revInd) > 0)
                    input[[n]]$profile[revInd,] <- 
                        input[[n]]$profile[revInd,ncol(input[[n]]$profile):1]
            else {
                warning("rpm signal type not yet supported for BAM file ",
                    "streaming!",immediate.=TRUE)
                input[[n]]$profile <- NULL
            }
        }
    }
    
    return(input)
}

# A read should not be assigned to a bin multiple times... In the same manner
# as overlapping genes
rpMatrixSample <- function(input,mask,binSize=200,
    interpolation=c("auto","spline","linear","neighborhood"),
    strand=NULL,ignoreStrand=TRUE,logScale=TRUE,rc=NULL) {
    interpolation=interpolation[1]
    
    # Maintain order (in the galaxy)
    theOrder <- names(mask)
    # Unlist the input because GRangesList sums the sub GRanges overlaps
    message("  joining")
    grs <- unlist(mask)
    # Count and re-split
    message("  counting")
    rl <- width(input)[1]
    rawCounts <- countOverlaps(grs,input,ignore.strand=ignoreStrand)#,
        #minoverlap=0.9*rl)
    message("  splitting")
    countList <- split(unname(rawCounts),names(rawCounts))
    countList <- countList[theOrder]
    # Find those with length < or > binSize so as to interpolate
    L <- lengths(countList)
    # ...and a failsafe for NULLs caused maybe by trimming
    ze <- which(L == 0)
    if (length(ze) > 0)
        countList[ze] <- lapply(countList[ze],function(x,b) {
            return(rep(0,b))
        },binSize)
    offending <- which(L != binSize)
    # Interpolate
    if (length(offending) > 0)
        countList[offending] <- lapply(countList[offending],function(x,b,p) {
            interpolated <- interpolateSignal(x,b,p)
            interpolated[interpolated < 0] <- 0
            return(interpolated)
        },binSize,interpolation)
    # Assemble the matrix
    countMatrix <- do.call("rbind",countList)
    # Get the measurement
    normFactor <- length(input)/binSize
    if (logScale)
        statMatrix <- log(sweep(countMatrix+1,2,1e+6/normFactor,"*"))
    else
        statMatrix <- sweep(countMatrix,2,1e+6/normFactor,"*")
    
    return(statMatrix)
}

# There is a problem with the bin creation... Forcing the bin size to the size
# specified, enforces very small bins, where the same read is counted an insane
# number of times. The number of bins should differ for each gene/region 
# according to its length and then inter/extrapolate... We need a dynamic
# bin creation procedure which should not be that hard...
# The number of bins should be dependent on i) the total number of bins
# (starting from 1), ii) the read length, iii) the gene length
# - The minimum bin length should be equal to the read length, so as to be able
#   to assess where a read belongs
# - One read is assigned to the area where >50% of the read belongs to
#   (minoverlap = round(0.5*read length)
#
# So it goes like: 
# number of bins = gene length / read length
# the minimum is 2, if we have one we rep
# there is no maximum as we finally interpolate to binSize
