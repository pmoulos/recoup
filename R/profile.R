profileMatrix <- function(input,flank,binParams,rc=NULL) {
    hasProfile <- vapply(input,function(x) is.null(x$profile),logical(1))
    if (!any(hasProfile))
        return(input)
    haveEqualLengths <- covLengthsEqual(input[[1]]$coverage)
    for (n in names(input)) {
        message("Calculating profile for ",input[[n]]$name)
        input[[n]]$profile <- profileMatrixSample(input[[n]]$coverage,flank,
            binParams,haveEqualLengths,rc=rc)
    }
    return(input)
}

profileMatrixSample <- function(x,flank,binParams,haveEqualLengths,rc=NULL) {
    if (!haveEqualLengths) {
        message(" center")
        center <- binCoverageMatrix(
            x,binSize=binParams$regionBinSize,stat=binParams$sumStat,
            interpolation=binParams$interpolation,flank=flank,where="center",
            chunks=binParams$chunks,rc=rc
        )
        if (binParams$flankBinSize!=0) {
            r <- flank/sum(flank)
            if (flank[1]==0)
                left <- NULL
            else {
                message(" upstream")
                left <- binCoverageMatrix(
                    x,binSize=round(2*binParams$flankBinSize*r[1]),
                    stat=binParams$sumStat,
                    interpolation=binParams$interpolation,flank=flank,
                    where="upstream",chunks=binParams$chunks,rc=rc
                )
            }
            if (flank[2]==0)
                right <- NULL
            else {
                message(" downstream")
                right <- binCoverageMatrix(
                    x,binSize=round(2*binParams$flankBinSize*r[2]),
                    stat=binParams$sumStat,
                    interpolation=binParams$interpolation,flank=flank,
                    where="downstream",chunks=binParams$chunks,rc=rc
                )
            }
        }
        else {
            if (flank[1]==0)
                left <- NULL
            else {
                message(" upstream")
                left <- baseCoverageMatrix(x,flank=flank,where="upstream",
                    chunks=binParams$chunks,rc=rc)
            }
            if (flank[2]==0)
                right <- NULL
            else {
                message(" downstream")
                right <- baseCoverageMatrix(x,flank=flank,where="downstream",
                    chunks=binParams$chunks,rc=rc)
            }
        }
        prof <- cbind(left,center,right)
        rownames(prof) <- names(x)
    }
    else {
        if (binParams$regionBinSize!=0)
            prof <- binCoverageMatrix(x,binSize=binParams$regionBinSize,
                stat=binParams$sumStat,chunks=binParams$chunks,rc=rc)
        else
            prof <- baseCoverageMatrix(x,chunks=binParams$chunks,rc=rc)
    }
    return(prof)
}

baseCoverageMatrix <- function(cvrg,flank=NULL,
        where=c("upstream","downstream"),chunks=NULL,rc=NULL) {
    if (is.null(flank)) {
        if (is.null(chunks)) {
            size <- length(cvrg[[1]])
            if (size==0) {
                for (i in 2:length(cvrg)) {
                    if (!is.na(runValue(cvrg[[i]]))) {
                        size <- length(cvrg[[i]])
                        break
                    }
                }
            }
            tmp <- cmclapply(cvrg,contVector,size,NULL,NULL,rc=rc)
        }
        else {
            tmpc <- vector("list",length(chunks))
            names(tmpc) <- names(chunks)
            for (n in names(chunks)) {
                cvrgc <- cvrg[chunks[[n]]]
                size <- length(cvrgc[[1]])
                if (size==0) {
                    for (i in 2:length(cvrgc)) {
                        if (!is.na(runValue(cvrgc[[i]]))) {
                            size <- length(cvrgc[[i]])
                            break
                        }
                    }
                }
                tmpc[[n]] <- cmclapply(cvrgc,contVector,size,NULL,NULL,rc=rc)
            }
            tmp <- do.call("c",tmpc)
        }
    }
    else {
        if (is.null(chunks))
            tmp <- cmclapply(cvrg,contVector,NULL,flank,where)
        else {
            tmpc <- vector("list",length(chunks))
            names(tmpc) <- names(chunks)
            for (n in names(chunks)) {
                cvrgc <- cvrg[chunks[[n]]]
                tmpc[[n]] <- cmclapply(cvrgc,contVector,NULL,flank,where,rc=rc)
            }
            tmp <- do.call("c",tmpc)
        }
    }
    return(do.call("rbind",tmp))
}

binCoverageMatrix <- function(cvrg,binSize=1000,stat=c("mean","median"),
    interpolation=c("auto","spline","linear","neighborhood"),flank=NULL,
    where=c("center","upstream","downstream"),chunks=NULL,rc=NULL) {
    where <- where[1]
    stat <- stat[1]
    interpolation=interpolation[1]
    if (is.null(chunks)) {
        tmp <- cmclapply(cvrg,function(x) splitVector(x,binSize,flank,where,
            interpolation,stat),rc=rc)
    }
    else {
        tmpc <- vector("list",length(chunks))
        names(tmpc) <- names(chunks)
        for (n in names(chunks)) {
            message("    processing chunk ",n)
            cvrgc <- cvrg[chunks[[n]]]
            tmpc[[n]] <- cmclapply(cvrgc,function(x) splitVector(x,binSize,
                flank,where,interpolation,stat),rc=rc)
        }
        tmp <- do.call("c",tmpc)
        names(tmp) <- names(cvrg)
    }
    statMatrix <- do.call("rbind",cmclapply(tmp,function(x) unlist(x)))
    return(statMatrix)
}

imputeProfile <- function(input,method=c("knn","simple"),rc=NULL) {
    hasMissing <- vapply(input,function(x) any(is.na(x$profile)),logical(1))
    if (!any(hasMissing))
        return(input)
    
    method <- method[1]
    if (method == "knn" && !requireNamespace("impute"))
        stop("R package impute is required to perform the imputation!")
    
    ii <- which(hasMissing)
    mn <- vapply(input[ii],function(x) x$name,character(1))
    message("Missing values detected in profiles of ",paste0(mn,collapse=", "),
        "! Imputing...")
    for (n in names(ii)) {
        message("Imputing missing values for ",input[[n]]$name)
        # NaNs are extremely rare and sparse, so few neighbors
        if (method == "simple") {
            # Impute based on the average of its 4 neighboring values
            miss <- which(is.na(input[[n]]$profile),arr.ind=TRUE)
            for (r in seq_len(nrow(miss))) {
                input[[n]]$profile[miss[r,1],miss[r,2]] <- 
                    mean(input[[n]]$profile[miss[r,1],c(miss[r,2]-2,miss[r,2]-1,
                        miss[r,2]+1,miss[r,2]+2)])
            }
        }
        else if (method == "knn") {
            tmp <- impute.knn(input[[n]]$profile,k=3)
            input[[n]]$profile <- tmp$data
        }
    }
    return(input)
}

################################################################################

# Legacy functions

baseCoverageMatrixOld <- function(cvrg,flank=NULL,
        where=c("upstream","downstream"),chunks=NULL,rc=NULL) {
    if (is.null(flank)) {
        if (is.null(chunks)) {
            size <- length(cvrg[[1]])
            if (size==0) {
                for (i in 2:length(cvrg)) {
                    #if (!is.null(cvrg[[i]])) {
                    if (!is.na(runValue(cvrg[[i]]))) {
                        size <- length(cvrg[[i]])
                        break
                    }
                }
            }
            tmp <- cmclapply(cvrg,function(x) {
                if (is(x,"Rle") && !is.na(runValue(x)))
                    return(as.numeric(x))
                else
                    return(NULL)
            },rc=rc)
        }
        else {
            tmpc <- vector("list",length(chunks))
            names(tmpc) <- names(chunks)
            for (n in names(chunks)) {
                cvrgc <- cvrg[chunks[[n]]]
                size <- length(cvrgc[[1]])
                if (size==0) {
                    for (i in 2:length(cvrgc)) {
                        if (!is.na(runValue(cvrgc[[i]]))) {
                            size <- length(cvrgc[[i]])
                            break
                        }
                    }
                }
                tmpc[[n]] <- cmclapply(cvrgc,function(x) {
                    if (is(x,"Rle") && !is.na(runValue(x)))
                        return(as.numeric(x))
                    else
                        return(NULL)
                },rc=rc)
            }
            tmp <- do.call("c",tmpc)
        }
        null <- which(vapply(tmp,is.null,logical(1)))
        if (length(null)>0) {
            for (j in null) {
                fill <- rep(0,size)
                tmp[[j]] <- fill
            }
        }
    }
    else {
        if (is.null(chunks)) {
            nr <- lengths(cvrg)
            switch(where,
                upstream = {
                    size <- flank[1]
                    tmp <- cmclapply(1:length(cvrg),function(i,cvrg,flank) {
                        if (is(cvrg[[i]],"Rle") && !is.na(runValue(cvrg[[i]])))
                            return(as.numeric(cvrg[[i]][1:flank[1]]))
                         else
                            return(NULL)
                    },cvrg,flank,rc=rc)
                },
                downstream = {
                    size <- flank[2]
                    tmp <- cmclapply(1:length(cvrg),function(i,cvrg,flank,nr) {
                        if (is(cvrg[[i]],"Rle") && !is.na(runValue(cvrg[[i]])))
                            return(as.numeric(
                                cvrg[[i]][(nr[i]-flank[2]+1):nr[i]]))
                        else
                            return(NULL)
                    },cvrg,flank,nr,rc=rc)
                }
            )
        }
        else {
            tmpc <- vector("list",length(chunks))
            names(tmpc) <- names(chunks)
            for (n in names(chunks)) {
                cvrgc <- cvrg[chunks[[n]]]
                nr <- lengths(cvrgc)
                switch(where,
                    upstream = {
                        size <- flank[1]
                        tmpc[[n]] <- 
                            cmclapply(1:length(cvrgc),function(i,cvrgc,flank) {
                            if (is(cvrgc[[i]],"Rle") 
                                && !is.na(runValue(cvrgc[[i]])))
                                return(as.numeric(cvrgc[[i]][1:flank[1]]))
                             else
                                return(NULL)
                        },cvrgc,flank,rc=rc)
                    },
                    downstream = {
                        size <- flank[2]
                        tmpc[[n]] <- 
                            cmclapply(1:length(cvrgc),
                            function(i,cvrgc,flank,nr) {
                            if (is(cvrgc[[i]],"Rle") 
                                && !is.na(runValue(cvrgc[[i]])))
                                return(as.numeric(
                                    cvrgc[[i]][(nr[i]-flank[2]+1):nr[i]]))
                            else
                                return(NULL)
                        },cvrgc,flank,nr,rc=rc)
                    }
                )
            }
            tmp <- do.call("c",tmpc)
        }
        null <- which(vapply(tmp,is.null,logical(1)))
        if (length(null)>0) {
            for (j in null) {
                fill <- rep(0,size)
                tmp[[j]] <- fill
            }
        }
    }
    return(do.call("rbind",tmp))
}



################################################################################

# Legacy functions

#binCoverageMatrixOld <- function(cvrg,binSize=1000,stat=c("mean","median"),
#    interpolation=c("auto","spline","linear","neighborhood"),flank=NULL,
#    where=c("center","upstream","downstream"),chunks=NULL,rc=NULL) {
#    stat <- stat[1]
#    interpolation=interpolation[1]
#    if (is.null(flank)) {
#        if (is.null(chunks)) {
#            tmp <- cmclapply(cvrg,function(x) {
#                if (is(x,"Rle") && !is.na(runValue(x)))
#                    return(as.numeric(x))
#                else
#                    return(NULL)
#            },rc=rc)
#        }
#        else {
#            tmpc <- vector("list",length(chunks))
#            names(tmpc) <- names(chunks)
#            for (n in names(chunks)) {
#                message("    processing chunk ",n)
#                tmpc[[n]] <- cmclapply(cvrg[chunks[[n]]],function(x) {
#                    if (is(x,"Rle") && !is.na(runValue(x)))
#                        return(as.numeric(x))
#                    else
#                        return(NULL)
#                },rc=rc)
#            }
#            tmp <- do.call("c",tmpc)
#            names(tmp) <- names(cvrg)
#        }
#    }
#    else {
#        if (is.null(chunks)) {
#            nr <- lengths(cvrg)
#            switch(where,
#                center = {
#                    tmp <- cmclapply(1:length(cvrg),function(i,cvrg,flank,nr) {
#                        if (is(cvrg[[i]],"Rle") 
#                            && !is.na(runValue(cvrg[[i]]))) {
#                            return(
#                                as.numeric(
#                                    cvrg[[i]][(flank[1]+1):(nr[i]-flank[2])])
#                            )
#                        }
#                        else
#                            return(NULL)
#                    },cvrg,flank,nr,rc=rc)
#                },
#                upstream = {
#                    tmp <- cmclapply(1:length(cvrg),function(i,cvrg,flank) {
#                        if (is(cvrg[[i]],"Rle") 
#                            && !is.na(runValue(cvrg[[i]]))) {
#                            return(as.numeric(cvrg[[i]][1:flank[1]]))
#                        }
#                        else
#                            return(NULL)
#                    },cvrg,flank,rc=rc)
#                },
#                downstream = {
#                    tmp <- cmclapply(1:length(cvrg),function(i,cvrg,flank,nr) {
#                        if (is(cvrg[[i]],"Rle") 
#                            && !is.na(runValue(cvrg[[i]]))) {
#                            return(as.numeric(
#                                cvrg[[i]][(nr[i]-flank[2]+1):nr[i]]))
#                        }
#                        else
#                            return(NULL)
#                    },cvrg,flank,nr,rc=rc)
#                }
#            )
#        }
#        else {
#            tmpc <- vector("list",length(chunks))
#            names(tmpc) <- names(chunks)
#            for (n in names(chunks)) {
#                message("    processing chunk ",n)
#                cvrgc <- cvrg[chunks[[n]]]
#                nr <- lengths(cvrgc)
#                switch(where,
#                    center = {
#                        tmpc[[n]] <- 
#                            cmclapply(1:length(cvrgc),
#                            function(i,cvrgc,flank,nr) {
#                            if (is(cvrgc[[i]],"Rle") 
#                                && !is.na(runValue(cvrgc[[i]]))) {
#                                return(
#                                    as.numeric(cvrgc[[i]][
#                                        (flank[1]+1):(nr[i]-flank[2])])
#                                )
#                            }
#                            else
#                                return(NULL)
#                        },cvrgc,flank,nr,rc=rc)
#                    },
#                    upstream = {
#                        tmpc[[n]] <- 
#                            cmclapply(1:length(cvrgc),function(i,cvrgc,flank) {
#                            if (is(cvrgc[[i]],"Rle") 
#                                && !is.na(runValue(cvrgc[[i]]))) {
#                                return(as.numeric(cvrgc[[i]][1:flank[1]]))
#                            }
#                            else
#                                return(NULL)
#                        },cvrgc,flank,rc=rc)
#                    },
#                    downstream = {
#                        tmpc[[n]] <- 
#                            cmclapply(1:length(cvrgc),
#                            function(i,cvrgc,flank,nr) {
#                            if (is(cvrgc[[i]],"Rle") 
#                                && !is.na(runValue(cvrgc[[i]]))) {
#                                return(as.numeric(
#                                    cvrgc[[i]][(nr[i]-flank[2]+1):nr[i]]))
#                            }
#                            else
#                                return(NULL)
#                        },cvrgc,flank,nr,rc=rc)
#                    }
#                )
#            }
#            tmp <- do.call("c",tmpc)
#            names(tmp) <- names(cvrg)
#        }
#    }
#    null <- which(sapply(tmp,is.null))
#    if (length(null)>0) {
#        for (j in null) {
#            fill <- rep(0,binSize)
#            tmp[[j]] <- fill
#        }
#    }
#    if (is.null(chunks))
#        tmp <- cmclapply(tmp,function(x) splitVector(x,binSize,
#            interpolation,stat),rc=rc)
#    else {
#        tmpc <- vector("list",length(chunks))
#        names(tmpc) <- names(chunks)
#        for (n in names(chunks)) {
#            message("    processing chunk ",n)
#            y <- tmp[chunks[[n]]]
#            tryCatch({
#            tmpc[[n]] <- 
#                cmclapply(y,function(x) splitVector(x,binSize,interpolation,
#                    stat),rc=rc)
#            },error=function(e) print(e),finally="")
#        }
#        tmp <- do.call("c",tmpc)
#        names(tmp) <- names(cvrg)
#    }
#    ###########################################################################
#    # Good try but took too much time...
#    #tmp <- cmclapply(cvrg,function(x) {
#    #    if (is.null(x))
#    #         x <- Rle(rep(0,binSize))
#    #    splitVector(x,binSize,interpolation)
#    #},rc=rc)
#    ###########################################################################
#    statMatrix <- do.call("rbind",cmclapply(tmp,function(x) unlist(x)))
#    #statMatrix <- do.call("rbind",cmclapply(tmp,function(x) 
#    #    unlist(lapply(x,stat)),rc=rc))
#    return(statMatrix)
#}
