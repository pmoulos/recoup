profileMatrix <- function(input,flank,binParams,rc=NULL) {
    hasProfile <- sapply(input,function(x) is.null(x$profile))
    if (!any(hasProfile))
        return(input)
        
    len <- lengths(input[[1]]$coverage)
    z <- which(len==1)
    if (length(z)>0)
        len <- len[-z]
    haveEqualLengths <- ifelse(all(len==len[1]),TRUE,FALSE)
    
    #if (!is.null(input[[1]]$coverage$center)) {
    if (!haveEqualLengths) {
        for (n in names(input)) {
            message("Calculating profile for ",input[[n]]$name)
            message(" center")
            center <- binCoverageMatrix(
                input[[n]]$coverage,binSize=binParams$regionBinSize,
                stat=binParams$sumStat,interpolation=binParams$interpolation,
                flank=flank,where="center",chunks=binParams$chunks,rc=rc
            )
            if (binParams$flankBinSize!=0) {
                r <- flank/sum(flank)
                if (flank[1]==0)
                    left <- NULL
                else {
                    message(" upstream")
                    left <- binCoverageMatrix(
                        input[[n]]$coverage,
                        binSize=round(2*binParams$flankBinSize*r[1]),
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
                        input[[n]]$coverage,
                        binSize=round(2*binParams$flankBinSize*r[2]),
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
                    left <- baseCoverageMatrix(input[[n]]$coverage,flank=flank,
                        where="upstream",chunks=binParams$chunks,rc=rc)
                }
                if (flank[2]==0)
                    right <- NULL
                else {
                    message(" downstream")
                    right <- baseCoverageMatrix(input[[n]]$coverage,flank=flank,
                        where="downstream",chunks=binParams$chunks,rc=rc)
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
                            chunks=binParams$chunks,rc=rc)
            else
                input[[n]]$profile <- 
                    baseCoverageMatrix(input[[n]]$coverage,
                        chunks=binParams$chunks,rc=rc)
            input[[n]]$profile <- input[[n]]$profile
        }
    }
    return(input)
}

baseCoverageMatrix <- function(cvrg,flank=NULL,
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
                #if (class(x)=="Rle")
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
        null <- which(sapply(tmp,is.null))
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
                        #if (class(cvrg[[i]])=="Rle")
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
                            return(as.numeric(cvrg[[i]][(nr[i]-flank[2]+1):nr[i]]))
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
                        tmpc[[n]] <- cmclapply(1:length(cvrgc),function(i,cvrgc,flank) {
                            if (is(cvrgc[[i]],"Rle") && !is.na(runValue(cvrgc[[i]])))
                                return(as.numeric(cvrgc[[i]][1:flank[1]]))
                             else
                                return(NULL)
                        },cvrgc,flank,rc=rc)
                    },
                    downstream = {
                        size <- flank[2]
                        tmpc[[n]] <- cmclapply(1:length(cvrgc),function(i,cvrgc,flank,nr) {
                            if (is(cvrgc[[i]],"Rle") && !is.na(runValue(cvrgc[[i]])))
                                return(as.numeric(cvrgc[[i]][(nr[i]-flank[2]+1):nr[i]]))
                            else
                                return(NULL)
                        },cvrgc,flank,nr,rc=rc)
                    }
                )
            }
            tmp <- do.call("c",tmpc)
        }
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
    where=c("center","upstream","downstream"),chunks=NULL,rc=NULL) {
    stat <- stat[1]
    interpolation=interpolation[1]
    if (is.null(flank)) {
        if (is.null(chunks)) {
            tmp <- cmclapply(cvrg,function(x) {
                #if (class(x)=="Rle")
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
                message("    processing chunk ",n)
                tmpc[[n]] <- cmclapply(cvrg[chunks[[n]]],function(x) {
                    if (is(x,"Rle") && !is.na(runValue(x)))
                        return(as.numeric(x))
                    else
                        return(NULL)
                },rc=rc)
            }
            tmp <- do.call("c",tmpc)
            names(tmp) <- names(cvrg)
        }
    }
    else {
        if (is.null(chunks)) {
            nr <- lengths(cvrg)
            switch(where,
                center = {
                    tmp <- cmclapply(1:length(cvrg),function(i,cvrg,flank,nr) {
                        #if (class(cvrg[[i]])=="Rle") {
                        if (is(cvrg[[i]],"Rle") && !is.na(runValue(cvrg[[i]]))) {
                            return(
                                as.numeric(cvrg[[i]][(flank[1]+1):(nr[i]-flank[2])])
                            )
                        }
                        else
                            return(NULL)
                    },cvrg,flank,nr,rc=rc)
                },
                upstream = {
                    tmp <- cmclapply(1:length(cvrg),function(i,cvrg,flank) {
                        if (is(cvrg[[i]],"Rle") && !is.na(runValue(cvrg[[i]]))) {
                            return(as.numeric(cvrg[[i]][1:flank[1]]))
                        }
                        else
                            return(NULL)
                    },cvrg,flank,rc=rc)
                },
                downstream = {
                    tmp <- cmclapply(1:length(cvrg),function(i,cvrg,flank,nr) {
                        if (is(cvrg[[i]],"Rle") && !is.na(runValue(cvrg[[i]]))) {
                            return(as.numeric(cvrg[[i]][(nr[i]-flank[2]+1):nr[i]]))
                        }
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
                message("    processing chunk ",n)
                cvrgc <- cvrg[chunks[[n]]]
                nr <- lengths(cvrgc)
                switch(where,
                    center = {
                        tmpc[[n]] <- cmclapply(1:length(cvrgc),function(i,cvrgc,flank,nr) {
                            #if (class(cvrg[[i]])=="Rle") {
                            if (is(cvrgc[[i]],"Rle") && !is.na(runValue(cvrgc[[i]]))) {
                                return(
                                    as.numeric(cvrgc[[i]][(flank[1]+1):(nr[i]-flank[2])])
                                )
                            }
                            else
                                return(NULL)
                        },cvrgc,flank,nr,rc=rc)
                    },
                    upstream = {
                        tmpc[[n]] <- cmclapply(1:length(cvrgc),function(i,cvrgc,flank) {
                            if (is(cvrgc[[i]],"Rle") && !is.na(runValue(cvrgc[[i]]))) {
                                return(as.numeric(cvrgc[[i]][1:flank[1]]))
                            }
                            else
                                return(NULL)
                        },cvrgc,flank,rc=rc)
                    },
                    downstream = {
                        tmpc[[n]] <- cmclapply(1:length(cvrgc),function(i,cvrgc,flank,nr) {
                            if (is(cvrgc[[i]],"Rle") && !is.na(runValue(cvrgc[[i]]))) {
                                return(as.numeric(cvrgc[[i]][(nr[i]-flank[2]+1):nr[i]]))
                            }
                            else
                                return(NULL)
                        },cvrgc,flank,nr,rc=rc)
                    }
                )
            }
            tmp <- do.call("c",tmpc)
            names(tmp) <- names(cvrg)
        }
    }
    null <- which(sapply(tmp,is.null))
    if (length(null)>0) {
        for (j in null) {
            fill <- rep(0,binSize)
            tmp[[j]] <- fill
        }
    }
    if (is.null(chunks))
        tmp <- cmclapply(tmp,function(x) splitVector(x,binSize,interpolation,stat),
            rc=rc)
    else {
        tmpc <- vector("list",length(chunks))
        names(tmpc) <- names(chunks)
        for (n in names(chunks)) {
            message("    processing chunk ",n)
            y <- tmp[chunks[[n]]]
            tryCatch({
            tmpc[[n]] <- 
                cmclapply(y,function(x) splitVector(x,binSize,interpolation,stat),
                    rc=rc)
            },error=function(e) print(e),finally="")
        }
        tmp <- do.call("c",tmpc)
        names(tmp) <- names(cvrg)
    }
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

#profileMatrixSeg <- function(covs,flank,binParams,rc=NULL) {
#    hasProfile <- sapply(input,function(x) is.null(x$profile))
#    if (!any(hasProfile))
#        return(input)
#        
#    len <- lengths(covs)
#    z <- which(len==0)
#    if (length(z)>0)
#        len <- len[-z]
#    haveEqualLengths <- ifelse(all(len==len[1]),TRUE,FALSE)
#    
#    if (!haveEqualLengths) {
#       message(" center")
#       center <- binCoverageMatrix(
#           covs,binSize=binParams$regionBinSize,
#           stat=binParams$sumStat,interpolation=binParams$interpolation,
#           flank=flank,where="center",rc=rc
#       )
#       if (binParams$flankBinSize!=0) {
#           r <- flank/sum(flank)
#           if (flank[1]==0)
#               left <- NULL
#           else {
#               message(" upstream")
#               left <- binCoverageMatrix(
#                   covs,binSize=round(2*binParams$flankBinSize*r[1]),
#                   stat=binParams$sumStat,flank=flank,where="upstream",
#                   interpolation=binParams$interpolation,rc=rc
#               )
#           }
#           if (flank[2]==0)
#               right <- NULL
#           else {
#               message(" downstream")
#               right <- binCoverageMatrix(
#                   covs,binSize=round(2*binParams$flankBinSize*r[2]),
#                   stat=binParams$sumStat,flank=flank,where="downstream",
#                   interpolation=binParams$interpolation,rc=rc
#               )
#           }
#       }
#       else {
#           if (flank[1]==0)
#               left <- NULL
#           else {
#               message(" upstream")
#               left <- baseCoverageMatrix(covs,flank=flank,where="upstream",
#                   rc=rc)
#           }
#           if (flank[2]==0)
#               right <- NULL
#           else {
#               message(" downstream")
#               right <- baseCoverageMatrix(covs,flank=flank,where="downstream",
#                   rc=rc)
#           }
#       }
#       profile <- cbind(left,center,right)
#       rownames(profile) <- names(covs)
#   }
#    else {
#       for (n in names(input)) {
#            if (binParams$regionBinSize!=0)
#                profile <- 
#                    binCoverageMatrix(covs,binSize=binParams$regionBinSize,
#                       stat=binParams$sumStat,rc=rc)
#            else
#                profile <- baseCoverageMatrix(covs,rc=rc)
#        }
#    }
#    return(profile)
#}
