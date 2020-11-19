#splitBySeqname <- function(gr,rc=NULL) {
#    #message("Splitting input regions by seqname...")
#    grList <- cmclapply(levels(seqnames(gr)),function(x,lib) {
#        #message("  ",x)
#        tmp <- lib[seqnames(lib)==x]
#        if (length(tmp)>0) return(tmp) else return(NULL)
#    },gr,rc=rc)
#    names(grList) <- levels(seqnames(gr))
#    null <- which(sapply(grList,is.null))
#    if (length(null)>0)
#        grList <- grList[-null]
#    return(grList)
#}

contVector <- function(x,size=NULL,flank=NULL,where=NULL) {
    if (!is.null(flank)) {
        switch(where,
            upstream = {
                size <- flank[1]
                if (is(x,"Rle") && !is.na(runValue(x)))
                    x <- x[seq_len(flank[1])]
                else
                    return(rep(0,size))
            },
            downstream = {
                size <- flank[2]
                if (is(x,"Rle") && !is.na(runValue(x)))
                    x <- x[(length(x)-flank[2]+1):length(x)]
                else
                    return(rep(0,size))
            }
        )
    }
    else {
        if (!is(x,"Rle") || all(is.na(runValue(x))))
            return(rep(0,size))
    }
    return(as.numeric(x))
}

splitVector <- function(x,n,flank,where,interp,stat) {
    isRle <- ifelse(is(x,"Rle"),TRUE,FALSE)
    if (!is.null(flank)) {
        switch(where,
            center = {
                if (is(x,"Rle") && all(!is.na(runValue(x))))
                    x <- x[(flank[1]+1):(length(x)-flank[2])]
                else
                    return(rep(0,n))
            },
            upstream = {
                if (is(x,"Rle") && all(!is.na(runValue(x))))
                    x <- x[seq_len(flank[1])]
                else
                    return(rep(0,n))
            },
            downstream = {
                if (is(x,"Rle") && all(!is.na(runValue(x))))
                    x <- x[(length(x)-flank[2]+1):length(x)]
                else
                    return(rep(0,n))
            },
            locus = {
                if (is(x,"Rle") && all(!is.na(runValue(x))))
                    x <- x[seq_len(length(x))]
                else
                    return(rep(0,n))
            }
        )
    }
    else {
        if (!is(x,"Rle") || all(is.na(runValue(x))))
            return(rep(0,n))
    }
    if (length(x)<n) {
        if (isRle)
            x <- as.numeric(x)
        x <- interpolateSignal(x,n,interp)
        if (isRle)
            x <- Rle(x)
    }
    binSize <- floor(length(x)/n)
    dif <- length(x) - binSize*n 
    binFac <- rep(binSize,n)
    # Random bin increase size to avoid problems
    add <- sample(seq_len(n),dif)
    binFac[add] <- binFac[add]+1
    f <- factor(rep(seq_len(n),binFac))
    #S <- split(x,f)
    #return(llply(S,stat))
    S <- split(seq_len(length(x)),f)
    s <- vapply(S,function(x) x[1],integer(1))
    return(aggregate(x,FUN=stat,start=s,width=binFac))
}

interpolateSignal <- function(x,n,interp) {
    switch(interp,
        auto = {
            d <- (n-length(x))/n
            if (d < 0.2) { # Then quite safe for neighborhood method
                y <- rep(NA,n)
                y[c(1,2)] <- x[c(1,2)]
                y[(n-1):n] <- x[(length(x)-1):length(x)]
                orig.pos <- sort(sample(3:(n-2),length(x)-4))
                y[orig.pos] <- x[3:(length(x)-2)]
                na <- which(is.na(y))
                avinds <- lapply(na,function(z) {
                    return(c(z-2,z-1,z+1,z+2))
                })
                xx <- unlist(lapply(avinds,function(ii,yy) {
                    return(mean(yy[ii],na.rm=TRUE))
                },y))
                y[na] <- xx
                x <- y
            }
            else { # Spline is the safest
                x <- spline(x,n=n)$y
                x[x<0] <- 0
            }
        },
        spline = {
            x <- spline(x,n=n)$y
            x[x<0] <- 0
        },
        linear = {
            x <- approx(x,n=n)$y
            x[x<0] <- 0
        },
        neighborhood = {
            y <- rep(NA,n)
            y[c(1,2)] <- x[c(1,2)]
            y[(n-1):n] <- x[(length(x)-1):length(x)]
            orig.pos <- sort(sample(3:(n-2),length(x)-4))
            y[orig.pos] <- x[3:(length(x)-2)]
            na <- which(is.na(y))
            avinds <- lapply(na,function(z) {
                return(c(z-2,z-1,z+1,z+2))
            })
            xx <- unlist(lapply(avinds,function(ii,yy) {
                return(mean(yy[ii],na.rm=TRUE))
            },y))
            y[na] <- xx
            x <- y
        }
    )
    return(x)
}

covLengthsEqual <- function(rl) {
    if (is(rl,"RleList") || is(rl,"list"))
        len <- lengths(rl)
    else if (is(rl,"GRanges"))
        len <- width(rl)
    else if (is(rl,"GRangesList"))
        len <- vapply(width(rl),sum,integer(1))
    z <- which(len==1)
    if (length(z)>0)
        len <- len[-z]
    return(ifelse(all(len==len[1]),TRUE,FALSE))
}

readConfig <- function(input) {
    if (missing(input) || !file.exists(input))
        stop("File to read sample info from should be a valid existing text ",
            "file!")
    tab <- read.delim(input)
    if (is.null(tab$id))
        stop("Sample id column not found in ",input,"!")
    if (is.null(tab$file))
        stop("Sample file path column not found in ",input,"!")
    if (is.null(tab$format))
        stop("Sample file format column not found in ",input,"!")
    samples <- as.character(tab$id)
    if (length(samples) != length(unique(samples)))
        stop("Sample identifiers must be unique for each sample!")
    if (length(grep(" ",samples))>0)
        stop("White space is not allowed in sample ids!")
    if (is.null(tab$name))
        nams <- as.character(tab$name)
    else
        nams <- samples
    files <- as.character(tab$file)
    if (any(!file.exists(files))) {
        bi <- which(!file.exists(files))
        stop("Input file ",files[bi]," does not exist! Please check paths...")
    }
    formats <- as.character(tab$format)
    if (!all(formats %in% c("bam","bed","bigwig")))
        stop("Input formats must be one of \"bam\", \"bed\", \"bigwig\"")
    cls <- NULL
    if (!is.null(tab$color)) {
        cls <- as.character(tab$color)
        chkcls <- areColors(cls)
        if (!all(chkcls)) {
            warning("Invalid colors found in color column in config file ",
                input,"! Will use automatic colors...")
            cls <- NULL
        }
    }
    output <- vector("list",nrow(tab))
    for (i in seq_len(nrow(tab))) {
        output[[i]]$id <- samples[i]
        output[[i]]$name <- nams[i]
        output[[i]]$file <- files[i]
        output[[i]]$format <- formats[i]
        output[[i]]$color <- cls[i]
        output[[i]]$ranges <- NULL
        output[[i]]$coverage <- NULL
        output[[i]]$profile <- NULL
    }
    names(output) <- samples
    return(output)
}

kmeansDesign <- function(input,design=NULL,kmParams) {
    if (missing(kmParams)) {
        if (!is.null(input$data))
            kmParams <- getr(input,"kmParams")
        else
            kmParams <- getDefaultListArgs("kmParams")
    }
    else {
        kmParamsDefault <- getDefaultListArgs("kmParams")
        kmParams <- setArg(kmParamsDefault,kmParams)
        kmParams <- validateListArgs("kmParams",kmParams)
    }
    if (!is.null(input$data)) # Fed with recoup object
        input <- input$data
    hasProfile <- vapply(input,function(x) is.null(x$profile),logical(1))
    if (any(hasProfile))
        stop("Profile matrices for k-means clustering are missing from the ",
            "input object. Have you called the profileMatrix function?")
    if (kmParams$k>0) {
        if (is.null(kmParams$reference)) {
            # Merge matrices to one and perform k-means. As normally coverages 
            # are normalized (the user is responsible to tell recoup how to do 
            # this and has been done at this point), we are legalized to do that
            message("Performing k-means (k=",kmParams$k,") clustering on ",
                "total profiles")
            theBigMatrix <- do.call("cbind",lapply(input,function(x) {
                return(as.matrix(x$profile))
            }))
            kcl <- kmeans(theBigMatrix,centers=kmParams$k,
                iter.max=kmParams$iterMax,nstart=kmParams$nstart,
                algorithm=kmParams$algorithm)
        }
        else {
            message("Performing k-means (k=",kmParams$k,") clustering using ",
                "the ",input[[kmParams$reference]]$name," sample profile as ",
                "reference")
            #set.seed(kmParams$seed)
            kcl <- kmeans(input[[kmParams$reference]]$profile,
                centers=kmParams$k,iter.max=kmParams$iterMax,
                nstart=kmParams$nstart,algorithm=kmParams$algorithm)
        }
        kmorder <- kcl$cluster
        #names(kmorder) <- rownames(input[[1]]$profile)
        if (!is.null(design)) {
            kmorder <- kmorder[rownames(design)]
            card <- table(kmorder)[kmorder]
            design$kcluster <- as.factor(paste("Cluster ",kmorder," (",card,
                ")",sep=""))
        }
        else {
            card <- table(kmorder)[kmorder]
            design <- data.frame(kcluster=paste("Cluster ",kmorder," (",card,
                ")",sep=""))
            rownames(design) <- names(kmorder)
        }
    }
    return(design)
}

# Experimental! Not working properly at the moment.
# In the final version, the design variable should be replaced by an "obj"
# variable. The obj can be a recoup list object (with design stored) or a design
# data frame
reorderClusters <- function(design,newOrder) {
    # Check if newOrder is numeric
    if (!is.numeric(newOrder))
        stop("newOrder must be numeric")
    # Check if design is a data frame
    if (!is.data.frame(design))
        stop("design must be a recoup design data frame!")
    # Check if design is a proper clutering design
    if (!("kcluster" %in% names(design)))
        stop("design must hold clustering information! Are you sure you run ",
            "kmeansDesign first?")
    # Check correct number of clusters in new order
    if (length(newOrder) != length(unique(as.character(design$kcluster))))
        stop("The number of clusters in the new order vector does not match ",
            "with the existing one!")
    
    # Do job
    cur <- as.character(design$kcluster)
    # They are of the form Cluster X (YYY)
    spl <- strsplit(cur," ")
    num <- as.numeric(vapply(spl,function(x) return(x[2]),character(1)))
    
    # Before going furher, check that newOrder contains the same clusters
    unum <- unique(num)
    poff <- intersect(unum,newOrder)
    if (length(poff) != length(unique(num)))
        stop("The provided clusters are not the same as in design! Cluster ",
            setdiff(newOrder,unum)," not found in design.")
    
    # Proceed
    car <- vapply(spl,function(x) return(x[3]),character(1))
    
    names(newOrder) <- seq_len(length(newOrder))
    
    # Finally, reorder
    newClus <- num
    for (i in unum)
        newClus[num==i] <- newOrder[as.character(i)]
    
    design$kcluster <- paste("Cluster",newClus,car)
    
    return(design)
}

sliceObj <- function(obj,i=NULL,j=NULL,k=NULL,dropPlots=FALSE,rc=NULL) {
    if (is.null(obj$data))
        stop("No data slot found in obj! Are you sure it's an output from",
            "recoup?")
    if (!is.null(obj$callopts$selector))
        obj$callopts$selector <- NULL
    if (!is.null(i)) {
        if (!is.numeric(i) && !is.character(i))
            stop("Horizontal indexing must be numeric or character!")
        for (s in seq_len(length(obj$data))) {
            if (!is.null(obj$data[[s]]$coverage))
                obj$data[[s]]$coverage <- obj$data[[s]]$coverage[i]
            if (!is.null(obj$data[[s]]$profile))
                obj$data[[s]]$profile <- obj$data[[s]]$profile[i,,drop=FALSE]
        }
        if (!is.null(obj$design))
            obj$design <- obj$design[i,,drop=FALSE]
    }
    if (!is.null(j)) {
        if (!is.numeric(j))
            stop("Vertical indexing must be numeric!")
        n <- ncol(obj$data[[1]]$profile)
        if (obj$callopts$region %in% c("tss","tes")
            || obj$callopts$customIsBase) {
            if (obj$callopts$binParams$regionBinSize==0)
                obj$callopts$flank <- c(min(j),n-max(j))
            else {
                f <- obj$callopts$flank/obj$callopts$binParams$regionBinSize
                obj$callopts$flank <- c(round(min(j)*f[1]),
                    round((n-max(j))*f[2]))
            }
        }
        else {
            if (obj$callopts$binParams$flankBinSize==0) {
                if (min(j)<obj$callopts$flank[1])
                    obj$callopts$flank[1] <- min(j)
                else
                    obj$callopts$flank[1] <- 0
                if ((n-max(j))<obj$callopts$flank[2])
                    obj$callopts$flank[2] <- n-max(j)
                else
                    obj$callopts$flank[2] <- 0
            }
            else {
                flankSpace <- 2*obj$callopts$binParams$flankBinSize
                oFlank <- obj$callopts$flank
                r <- obj$callopts$flank/sum(obj$callopts$flank)
                original.start.index <- round(r[1]*flankSpace)
                original.end.index <- round(n - r[2]*flankSpace)
                nonZeroFlank <- FALSE
                if (min(j)<original.start.index) {
                    obj$callopts$flank[1] <- 
                        round((original.start.index - min(j))*
                            obj$callopts$flank[1]/original.start.index)
                    nonZeroFlank <- TRUE
                }
                else
                    obj$callopts$flank[1] <- 0
                if (max(j)>original.end.index) {
                    obj$callopts$flank[2] <- 
                        round((max(j) - original.end.index)*
                            obj$callopts$flank[2]/(n - original.end.index))
                    nonZeroFlank <- TRUE
                }
                else
                    obj$callopts$flank[2] <- 0
                if (nonZeroFlank)
                    obj$callopts$binParams$flankBinSize <- round(0.5*flankSpace*
                            sum(obj$callopts$flank)/sum(oFlank))
                obj$callopts$binParams$forcedBinSize[1] <- 
                    obj$callopts$binParams$flankBinSize
                if (all(obj$callopts$flank==0))
                    obj$callopts$region <- "custom"
            }
        }
        for (s in seq_len(length(obj$data))) {
            if (!is.null(obj$data[[s]]$profile))
                obj$data[[s]]$profile <- obj$data[[s]]$profile[,j,drop=FALSE]
        }
    }
    if (!is.null(k)) {
        if (!is.numeric(k) && !is.character(k))
            stop("Sample indexing must be numeric or character!")
        obj$data <- obj$data[k]
    }
    if (dropPlots)
        obj$plots <- list(profile=NULL,heatmap=NULL,correlation=NULL)
    else {
        message("Recalculating profile plots after slicing...")
        if (obj$callopts$plotParams$profile) {
            message("  profile")
            obj <- recoupProfile(obj,rc=rc)
        }
        if (obj$callopts$plotParams$heatmap) {
            message("  heatmap")
            obj <- recoupHeatmap(obj,rc=rc)
        }
        if (obj$callopts$plotParams$correlation) {
            message("  correlation")
            obj <- recoupCorrelation(obj,rc=rc)
        }
    }
    return(obj)
}

# TODO: Check if the inputs have no names and work with indices
# TODO: Check the callopts one by one. For example, there is the risk of 
#       different organisms if we don't check for "organism".
# TODO: Recalculate plots after merge
mergeRuns <- function(...,withDesign=c("auto","drop"),dropPlots=TRUE) {
    tmp <- list(...)
    withDesign <- tolower(withDesign[1])
    
    # Initialize new object
    merged <- list()
    
    # Check design
    if (withDesign == "auto") {
        allHaveDesign <- 
            all(vapply(tmp,function(x) !is.null(x$design),logical(1)))
        whichDesign <- vapply(tmp,function(x) is.null(x$design),logical(1))
        oneHasDesign <- ifelse(length(which(!whichDesign))==1,TRUE,FALSE)
        noneHasDesign <- all(whichDesign)
        if (noneHasDesign)
            merged$design <- NULL
        if (allHaveDesign) {
            tryCatch({
                # Will fail if some design has different number of columns
                test1 <- 
                    do.call("rbind",lapply(tmp,function(x) return(x$design)))
                # Will fail if some design has different number of rows
                test2 <- 
                    do.call("cbind",lapply(tmp,function(x) return(x$design)))
                # Will fail if some design has different row names
                n1 <- sort(Reduce("intersect",
                    lapply(tmp,function(x) return(rownames(x$design)))))
                n2 <- sort(rownames(tmp[[1]]$design))
                if (!identical(n1,n2))
                    stop("Incompatible row names")
                # Will fail if some design has different column names
                n1 <- sort(Reduce("intersect",
                    lapply(tmp,function(x) return(colnames(x$design)))))
                n2 <- sort(colnames(tmp[[1]]$design))
                if (!identical(n1,n2))
                    stop("Incompatible column names")
                
                # If nothing fails, then assign the design to the new object
                merged$design <- tmp[[1]]$design
            },error=function(e) {
                warning("Caught incompatible designs! Dropping design...",
                    immediate.=TRUE)
                message("----------")
                print(e)
                message("----------")
                merged$design <- NULL
            },
            finally="")
        }
        if (oneHasDesign) # Apply to all
            merged$design <- tmp[[which(whichDesign)]]$design
        # Some have design...
        if (!noneHasDesign && !allHaveDesign && !oneHasDesign) {
            # Gather them
            designs <- lapply(tmp,function(x) {
                if (!is.null(x$design))
                    return(x$design)
                else
                    return(NULL)
            })
            
            # ...and repeat the previous pattern
            tryCatch({
                # Will fail if some design has different number of columns
                test1 <- do.call("rbind",designs)
                # Will fail if some design has different number of rows
                test2 <- do.call("cbind",designs)
                # Will fail if some design has different row names
                n1 <- sort(Reduce("intersect",
                    lapply(designs,function(x) return(rownames(x)))))
                n2 <- sort(rownames(designs[[1]]))
                if (!identical(n1,n2))
                    stop("Incompatible row names")
                # Will fail if some design has different column names
                n1 <- sort(Reduce("intersect",
                    lapply(designs,function(x) return(colnames(x)))))
                n2 <- sort(colnames(designs[[1]]))
                if (!identical(n1,n2))
                    stop("Incompatible column names")
                
                # If nothing fails, then assign the design to the new object
                merged$design <- designs[[1]]
            },error=function(e) {
                warning("Caught incompatible designs! Dropping design...",
                    immediate.=TRUE)
                message("----------")
                print(e)
                message("----------")
                merged$design <- NULL
            },
            finally="")
        }
    }
    else if (withDesign=="drop")
        merged$design <- NULL
    
    # Quick for now
    merged$callopts <- tmp[[1]]$callopts
    
    # Try to merge actual data
    newLen <- sum(vapply(tmp,function(x) return(length(x$data)),integer(1)))
    merged$data <- vector("list",newLen)
    theNames <- character(0)
    for (i in seq_len(length(tmp)))
        theNames <- c(theNames,names(tmp[[i]]$data))
    names(merged$data) <- theNames

    for (m in tmp) {
        for (n in names(m$data)) {
            merged$data[[n]] <- m$data[[n]]
            if (!is.null(merged$design)) {
                if (!is.null(merged$data[[n]]$coverage))
                    merged$data[[n]]$coverage <- 
                        merged$data[[n]]$coverage[rownames(merged$design)]
                if (!is.null(merged$data[[n]]$profile))
                    merged$data[[n]]$profile <- 
                        merged$data[[n]]$profile[rownames(merged$design),]
            }
        }
    }
    
    # Drop all plots for now, later recalculate
    if (dropPlots)
        merged$plots <- list(profile=NULL,heatmap=NULL,correlation=NULL)
    else
        merged$plots <- NULL
    
    return(merged)
}

decideChanges <- function(input,currCall,prevCall) {
    if (is.null(prevCall))
        return(input)
    # Check region and flank
    if (currCall$region != prevCall$region 
        || !all(currCall$flank == prevCall$flank))
        input <- removeData(input,c("coverage","profile"))
    # Check binParams
    if (currCall$binParams$flankBinSize != prevCall$binParams$flankBinSize
        || currCall$binParams$regionBinSize != prevCall$binParams$regionBinSize
        || currCall$binParams$sumStat != prevCall$binParams$sumStat
        || currCall$binParams$interpolation != prevCall$binParams$interpolation
        || currCall$binParams$forceHeatmapBinning != 
            prevCall$binParams$forceHeatmapBinning
        || currCall$binParams$forcedBinSize != prevCall$binParams$forcedBinSize
        || currCall$binParams$chunking != prevCall$binParams$chunking)
        input <- removeData(input,"profile")
    if (currCall$preprocessParams$normalize != 
        prevCall$preprocessParams$normalize
        || currCall$preprocessParams$sampleTo != 
        prevCall$preprocessParams$sampleTo
        || currCall$preprocessParams$spliceAction != 
        prevCall$preprocessParams$spliceAction
        || currCall$preprocessParams$spliceRemoveQ != 
        prevCall$preprocessParams$spliceRemoveQ
        || currCall$preprocessParams$spliceRemoveQ != 
        prevCall$preprocessParams$spliceRemoveQ)
        #|| currCall$preprocessParams$seed != prevCall$preprocessParams$seed)
        input <- removeData(input,c("ranges","coverage","profile"))
        
        return(input)
}

removeData <- function(input,type=c("ranges","coverage","profile")) {
    type <- tolower(type)
    checkTextArgs("type",type,c("ranges","coverage","profile"),multiarg=TRUE)
    if (!is.null(input$data)) # Gave recoup output object
        for (i in seq_len(length(input$data)))
            input$data[[i]][type] <- NULL
    else
        for (i in seq_len(length(input)))
            input[[i]][type] <- NULL
    return(input)
}

calcLinearFactors <- function(input,preprocessParams) {
    hasRanges <- vapply(input,function(x) is.null(x$ranges),logical(1))
    if (any(hasRanges))
        stop("Please provide input reads before calculation normalization ",
            "factors")
    libSizes <- vapply(input,function(x) return(length(x$ranges)),integer(1))
    if (preprocessParams$normalize=="linear" 
        || preprocessParams$normalize=="downsample")
        linFac <- min(libSizes)/libSizes
    else if (preprocessParams$normalize=="sampleto")
        linFac <- preprocessParams$sampleTo/libSizes
    names(linFac) <- names(input)
    return(linFac)
}

cmclapply <- function(...,rc) {
    if (suppressWarnings(!requireNamespace("parallel")) 
        || .Platform$OS.type!="unix")
        m <- FALSE
    else {
        m <- TRUE
        ncores <- parallel::detectCores()
        if (ncores==1) 
            m <- FALSE
        else {
            if (!missing(rc) && !is.null(rc))
                ncores <- ceiling(rc*ncores)
            else 
                m <- FALSE
        }
    }
    if (m)
        return(mclapply(...,mc.cores=ncores,mc.set.seed=FALSE))
    else
        return(lapply(...))
}

ssCI <- function(fit) {
    res <- (fit$yin - fit$y)/(1-fit$lev)
    sigma <- sqrt(var(res))
    upper <- fit$y + 3.0*sigma*sqrt(fit$lev)
    lower <- fit$y - 3.0*sigma*sqrt(fit$lev)
    return(list(lower=lower,upper=upper))
}

getDefaultListArgs <- function(arg,design=NULL,genome=NULL) {
    switch(arg,
        orderBy = {
            return(list(
                what=c("none","suma","sumn","maxa","maxn","hc"),
                order=c("descending","ascending"),
                custom=NULL
            ))
        },
        binParams = {
            return(list(
                flankBinSize=0,
                regionBinSize=0,
                sumStat=c("mean","median"),
                interpolation=c("auto","spline","linear","neighborhood"),
                binType=c("variable","fixed"),
                forceHeatmapBinning=TRUE,
                forcedBinSize=c(50,200),
                chunking=FALSE#,
                #seed=42
            ))
        },
        preprocessParams = {
            return(list(
                fragLen=NA,
                cleanLevel=c(0,1,2,3),
                normalize=c("none","linear","downsample","sampleto"),
                sampleTo=1e+6,
                spliceAction=c("split","keep","remove"),
                spliceRemoveQ=0.75,
                #bedGenome=ifelse(!is.null(genome) && genome %in% c("hg18",
                #    "hg19","hg38","mm9","mm10","rn5","dm3","danrer7","pantro4",
                #    "susscr3"),genome,NA),
                bedGenome=NA#,
                #seed=42
            ))
        },
        selector = {
            return(list(
                id=NULL,
                biotype=NULL,
                exonType=NULL
            ))
        },
        strandedParams = {
            return(list(
                strand=NULL,
                ignoreStrand=TRUE
            ))
        },
        plotParams = {
            return(list(
                plot=TRUE,
                profile=TRUE,
                heatmap=TRUE,
                correlation=TRUE,
                signalScale=c("natural","log2"),
                heatmapScale=c("each","common"),
                heatmapFactor=1,
                corrScale=c("normalized","each"),
                sumStat=c("mean","median"),
                smooth=TRUE,
                corrSmoothPar=ifelse(is.null(design),0.1,0.5),
                singleFacet=c("none","wrap","grid"),
                multiFacet=c("wrap","grid"),
                conf=TRUE,
                device=c("x11","png","jpg","tiff","bmp","pdf","ps"),
                outputDir=".",
                outputBase=NULL
            ))
        },
        saveParams = {
            return(list(
                ranges=TRUE,
                coverage=TRUE,
                profile=TRUE,
                profilePlot=TRUE,
                heatmapPlot=TRUE,
                correlationPlot=TRUE
            ))
        },
        kmParams = {
            return(list(
                k=0,
                nstart=20,
                algorithm=c("Hartigan-Wong","Lloyd","Forgy","MacQueen"),
                iterMax=20,
                reference=NULL#,
                #seed=42
            ))
        }
    )
}

#setr <- function(obj,key=c("design","profile","heatmap","correlation",
#   "orderBy","kmParams","plotParams"),value) {
setr <- function(obj,key,value=NULL) {
    if (is.character(key)) { # Property name, property value pairs
        if (is.null(value))
            stop("When the key argument is a list of property names, the ",
                "value argument must not be NULL.")
        if (length(key)>1) { # Value must be a list of lists
            if (length(key)!=length(lengths(value)))
                stop("The number values to set must be equal to the number of ",
                    "properties to change")
        }
        else
            value <- list(value)
        valid <- key %in% c("design","profile","heatmap","correlation",
            "orderBy","kmParams","plotParams")
        not.valid <- which(!valid)
        if (length(not.valid)>0) {
            warning("The following parameters to set are are invalid and will ",
                "be ignored: ",paste(key[not.valid],collapse=", "),
                    immediate.=TRUE)
            key[not.valid] <- NULL
            value[not.valid] <- NULL
        }
        if (length(key)==0)
            return(obj)
        names(value) <- key
        key <- value
    }
    else if (is.list(key)) { # Only a named list with property names and values
        valid <- names(key) %in% c("design","profile","heatmap","correlation",
            "orderBy","kmParams","plotParams")
        not.valid <- which(!valid)
        if (length(not.valid)>0) {
            warning("The following parameters to set are are invalid and will ",
                "be ignored: ",paste(names(key)[not.valid],collapse=", "),
                    immediate.=TRUE)
            key[not.valid] <- NULL
        }
        if (length(key)==0)
            return(obj)
    }
    for (n in names(key)) {
        switch(n,
            design = {
                if (!is.null(key[[n]])) {
                    if (!is.data.frame(key[[n]]))
                        key[[n]] <- read.delim(key[[n]],row.names=1)
                    nfac <- ncol(key[[n]])
                    if (length(obj$data)>1 && nfac>2)
                        stop("When more than one files are provided for ",
                            "coverage generation, the maximum number of ",
                            "allowed value factors is 2")
                    if (length(obj$data)>1 && nfac>1 
                        && obj$callopts$kmParams$k>0)
                        stop("When more than one files are provided for ",
                            "coverage generation and k-means clustering is ",
                            "also requested, the maximum number of allowed ",
                            "design factors is 1")
                    if (length(obj$data)==1 && nfac>3)
                        stop("The maximum number of allowed design factors is",
                            " 3")
                    if (length(obj$data)==1 && nfac>2 
                        && obj$callopts$kmParams$k>0)
                        stop("The maximum number of allowed design factors ",
                            "when k-means clustering is requested is 2")
                }
                obj$design <- key[[n]]
            },
            profile = {
                if (!is(key[[n]],"gg") && !is(key[[n]],"ggplot")) {
                    warning("The supplied profile plot is not a ggplot ",
                        "object! Ignoring...",immediate.=TRUE)
                    return(obj)
                }
                obj$plots$profile <- key[[n]]
            },
            heatmap = {
                if (!is(key[[n]],"HeatmapList")) {
                    warning("The supplied heatmap plot is not a Heatmap ",
                        "object! Ignoring...",immediate.=TRUE)
                    return(obj)
                }
                obj$plots$heatmap <- key[[n]]
            },
            correlation = {
                if (!is(key[[n]],"gg") && !is(key[[n]],"ggplot")) {
                    warning("The supplied correlation plot is not a ggplot ",
                        " object! Ignoring...",immediate.=TRUE)
                    return(obj)
                }
                obj$plots$correlation <- key[[n]]
            },
            orderBy = {
                orderByDefault <- getDefaultListArgs("orderBy")
                orderBy <- setArg(orderByDefault,key[[n]])
                orderBy <- validateListArgs("orderBy",orderBy)
                obj$callopts$orderBy <- orderBy
                # Should change also complexHeatmapParams
                if (length(grep("hc",orderBy$what))>0
                    && !(obj$callopts$complexHeatmapParams$main$cluster_rows 
                    || obj$callopts$complexHeatmapParams$group$cluster_rows)) {
                    message("Changing also hierarchical clustering parameters ",
                        "to comply with new orderBy settings.")
                    obj$callopts$complexHeatmapParams$main$cluster_rows <- TRUE
                    obj$callopts$complexHeatmapParams$group$cluster_rows <- TRUE
                }
                if ((obj$callopts$complexHeatmapParams$main$cluster_rows 
                    || obj$callopts$complexHeatmapParams$group$cluster_rows)
                    && length(grep("hc",orderBy$what))==0) {
                    message("Changing also hierarchical clustering parameters ",
                        "to comply with new orderBy settings.")
                    obj$callopts$complexHeatmapParams$main$cluster_rows <- FALSE
                    obj$callopts$complexHeatmapParams$group$cluster_rows <- 
                        FALSE
                }
            },
            kmParams = {
                kmParamsDefault <- getDefaultListArgs("kmParams")
                kmParams <- setArg(kmParamsDefault,key[[n]])
                kmParams <- validateListArgs("kmParams",kmParams)
                obj$callopts$kmParams <- kmParams
            },
            plotParams = {
                plotParamsDefault <- getDefaultListArgs("plotParams")
                plotParams <- setArg(plotParamsDefault,key[[n]])
                plotParams <- validateListArgs("plotParams",plotParams)
                obj$callopts$plotParams <- plotParams
            }
        )
    }
    return(obj)
}

getr <- function(obj,key=c("design","profile","heatmap","correlation","orderBy",
    "kmParams","plotParams")) {
    checkTextArgs("key",key,c("design","profile","heatmap","correlation",
        "orderBy","kmParams","plotParams"))
    switch(key,
        design = {
            return(obj$design)
        },
        profile = {
            return(obj$plots$profile)
        },
        heatmap = {
            return(obj$plots$heatmap)
        },
        correlation = {
            return(obj$plots$correlation)
        },
        orderBy = {
            return(obj$callopts$orderBy)
        },
        kmParams = {
            return(obj$callopts$kmParams)
        },
        plotParams = {
            return(obj$callopts$plotParams)
        }
    )
}

getArg <- function(arg.list,arg.name) {
    return(arg.list[arg.name])
}

setArg <- function(arg.list,arg.name,arg.value=NULL) {
    if (is.list(arg.name))
        arg.list[names(arg.name)] <- arg.name
    else if (is.character(arg.name)) {
        tmp <- vector("list",length(arg.name))
        names(tmp) <- arg.name
        i <- 0
        for (n in arg.name) {
            i <- i + 1
            tmp[[n]] <- arg.value[i]
        }
        arg.list[arg.name] <- tmp
    }
    return(arg.list)
}

areColors <- function(x) {
    return(vapply(x,function(y) {
        tryCatch(is.matrix(col2rgb(y)),
            error=function(e) FALSE)
    },logical(1)))
}
