profileMatrix <- function(input,flank,binParams,rc=NULL) {
    hasProfile <- sapply(input,function(x) is.null(x$profile))
    if (!any(hasProfile))
        return(input)
        
    len <- lengths(input[[1]]$coverage)
    z <- which(len==0)
    if (length(z)>0)
        len <- len[-z]
    haveEqualLengths <- ifelse(all(len==len[1]),TRUE,FALSE)
    
    #if (!is.null(input[[1]]$coverage$center)) {
    if (!haveEqualLengths) {
        for (n in names(input)) {
            message("Calculating profile for ",input[[n]]$name)
            message(" center")
            #center <- binCoverageMatrix(input[[n]]$coverage$center,
            #    binSize=binParams$regionBinSize,stat=binParams$sumStat,
            #    interpolation=binParams$interpolation,rc=rc)
            center <- binCoverageMatrix(
                input[[n]]$coverage,binSize=binParams$regionBinSize,
                stat=binParams$sumStat,interpolation=binParams$interpolation,
                flank=flank,where="center",rc=rc
            )
            if (binParams$flankBinSize!=0) {
                r <- flank/sum(flank)
                if (flank[1]==0)
                    left <- NULL
                else {
                    message(" upstream")
                    #left <- binCoverageMatrix(input[[n]]$coverage$upstream,
                    #    binSize=binParams$flankBinSize,stat=binParams$sumStat,
                    #    interpolation=binParams$interpolation,rc=rc)
                    left <- binCoverageMatrix(
                        input[[n]]$coverage,
                        binSize=round(2*binParams$flankBinSize*r[1]),
                        stat=binParams$sumStat,
                        interpolation=binParams$interpolation,flank=flank,
                        where="upstream",rc=rc
                    )
                }
                if (flank[2]==0)
                    right <- NULL
                else {
                    message(" downstream")
                    #right <- binCoverageMatrix(input[[n]]$coverage$downstream,
                    #    binSize=binParams$flankBinSize,stat=binParams$sumStat,
                    #    interpolation=binParams$interpolation,rc=rc)
                    right <- binCoverageMatrix(
                        input[[n]]$coverage,
                        binSize=round(2*binParams$flankBinSize*r[2]),
                        stat=binParams$sumStat,
                        interpolation=binParams$interpolation,flank=flank,
                        where="downstream",rc=rc
                    )
                }
            }
            else {
                if (flank[1]==0)
                    left <- NULL
                else {
                    message(" upstream")
                    #left <- baseCoverageMatrix(input[[n]]$coverage$upstream,
                    #   rc=rc)
                    left <- baseCoverageMatrix(input[[n]]$coverage,flank=flank,
                        where="upstream",rc=rc)
                }
                if (flank[2]==0)
                    right <- NULL
                else {
                    message(" downstream")
                    #right <- baseCoverageMatrix(input[[n]]$coverage$downstream,
                    #    rc=rc)
                    right <- baseCoverageMatrix(input[[n]]$coverage,flank=flank,
                        where="downstream",rc=rc)
                }
            }
            input[[n]]$profile <- cbind(left,center,right)
            rownames(input[[n]]$profile) <- names(input[[n]]$coverage)
            input[[n]]$profile <- input[[n]]$profile
        }
    }
    else {
        for (n in names(input)) {
            message("Calculating profile for ",input[[n]]$name)
            if (binParams$regionBinSize!=0)
                input[[n]]$profile <- 
                    binCoverageMatrix(input[[n]]$coverage,
                        binSize=binParams$regionBinSize,stat=binParams$sumStat,
                            rc=rc)
            else
                input[[n]]$profile <- 
                    baseCoverageMatrix(input[[n]]$coverage,rc=rc)
            input[[n]]$profile <- input[[n]]$profile
        }
    }
    return(input)
}

baseCoverageMatrix <- function(cvrg,flank=NULL,
        where=c("upstream","downstream"),rc=NULL) {
    if (is.null(flank)) {
        size <- length(cvrg[[1]])
        if (size==0) {
            for (i in 2:length(cvrg)) {
                if (!is.null(cvrg[[i]])) {
                    size <- length(cvrg[[i]])
                    break
                }
            }
        }
        tmp <- cmclapply(cvrg,function(x) {
            if (class(x)=="Rle")
                return(as.numeric(x))
        },rc=rc)
        null <- which(sapply(tmp,is.null))
        if (length(null)>0) {
            for (j in null) {
                fill <- rep(0,size)
                tmp[[j]] <- fill
            }
        }
    }
    else {
        nr <- lengths(cvrg)
        switch(where,
            upstream = {
                size <- flank[1]
                tmp <- cmclapply(1:length(cvrg),function(i,cvrg,flank) {
                    if (class(cvrg[[i]])=="Rle")
                        return(as.numeric(cvrg[[i]][1:flank[1]]))
                },cvrg,flank,rc=rc)
            },
            downstream = {
                size <- flank[2]
                tmp <- cmclapply(1:length(cvrg),function(i,cvrg,flank,nr) {
                    if (class(cvrg[[i]])=="Rle")
                        return(as.numeric(cvrg[[i]][(nr[i]-flank[2]+1):nr[i]]))
                },cvrg,flank,nr,rc=rc)
            }
        )
        null <- which(sapply(tmp,is.null))
        if (length(null)>0) {
            for (j in null) {
                fill <- rep(0,size)
                tmp[[j]] <- fill
            }
        }
    }
    return(do.call("rbind",tmp))
}

binCoverageMatrix <- function(cvrg,binSize=1000,stat=c("mean","median"),
    interpolation=c("auto","spline","linear","neighborhood"),flank=NULL,
    where=c("center","upstream","downstream"),rc=NULL) {
    stat <- stat[1]
    interpolation=interpolation[1]
    if (is.null(flank))
        tmp <- cmclapply(cvrg,function(x) {
            if (class(x)=="Rle")
                return(as.numeric(x))
        },rc=rc)
    else {
        nr <- lengths(cvrg)
        switch(where,
            center = {
                tmp <- cmclapply(1:length(cvrg),function(i,cvrg,flank,nr) {
                    if (class(cvrg[[i]])=="Rle") {
                        return(
                            as.numeric(cvrg[[i]][(flank[1]+1):(nr[i]-flank[2])])
                        )
                    }
                },cvrg,flank,nr,rc=rc)
            },
            upstream = {
                tmp <- cmclapply(1:length(cvrg),function(i,cvrg,flank) {
                    if (class(cvrg[[i]])=="Rle") {
                        return(as.numeric(cvrg[[i]][1:flank[1]]))
                    }
                },cvrg,flank,rc=rc)
            },
            downstream = {
                tmp <- cmclapply(1:length(cvrg),function(i,cvrg,flank,nr) {
                    if (class(cvrg[[i]])=="Rle") {
                        return(as.numeric(cvrg[[i]][(nr[i]-flank[2]+1):nr[i]]))
                    }
                },cvrg,flank,nr,rc=rc)
            }
        )
    }
    null <- which(sapply(tmp,is.null))
    if (length(null)>0) {
        for (j in null) {
            fill <- rep(0,binSize)
            tmp[[j]] <- fill
        }
    }
    tmp <- cmclapply(tmp,function(x) splitVector(x,binSize,interpolation,stat),
        rc=rc)
    ############################################################################
    # Good try but took too much time...
    #tmp <- cmclapply(cvrg,function(x) {
    #    if (is.null(x))
    #         x <- Rle(rep(0,binSize))
    #    splitVector(x,binSize,interpolation)
    #},rc=rc)
    ############################################################################
    statMatrix <- do.call("rbind",cmclapply(tmp,function(x) unlist(x)))
    #statMatrix <- do.call("rbind",cmclapply(tmp,function(x) 
    #    unlist(lapply(x,stat)),rc=rc))
    return(statMatrix)
}
