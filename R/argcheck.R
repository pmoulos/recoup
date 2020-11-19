checkMainArgs <- function(main.args) {
    in.args <- names(main.args)[-1] # 1st member name of calling function
    valid.args <- c("input","design","region","type","signal","genome","refdb",
        "flank","onFlankFail","fraction","orderBy","binParams","selector",
        "preprocessParams","plotParams","saveParams","kmParams",
        "strandedParams","ggplotParams","complexHeatmapParams","bamParams",
        "onTheFly","localDb","rc")
    invalid <- setdiff(in.args,valid.args)
    if (length(invalid) > 0) {
        for (i in seq_len(length(invalid)))
            warning("Unknown input argument to recoup function: ",invalid[i],
                " ...Ignoring...",immediate.=TRUE)
    }
}

checkTextArgs <- function(arg.name,arg.value,arg.list,multiarg=FALSE) {
    if (multiarg) {
        arg.value <- arg.value
        if (!all(arg.value %in% arg.list))
            stop(arg.name," parameter must be one or more of ",
                paste(paste("\"",arg.list,sep=""),collapse="\", "),"\"!")
    }
    else {
        arg.value <- arg.value[1]
        if (!(arg.value %in% arg.list))
            stop(arg.name," parameter must be one of ",
                paste(paste("\"",arg.list,sep=""),collapse="\", "),"\"!")
    }
}

checkNumArgs <- function(arg.name,arg.value,arg.type,arg.bounds,direction) {
    switch(arg.type,
        numeric = {
            if (!is.numeric(arg.value)) {
                # Can it be converted to numeric?
                arg.value <- suppressWarnings(as.numeric(arg.value))
                if (is.na(arg.value)) # No
                    stop(arg.name," parameter must be a numeric value!")
            }
            if (!missing(arg.bounds)) {
                switch(direction,
                    both = {
                        if (arg.value<=arg.bounds[1] ||
                            arg.value>=arg.bounds[2])
                            stop(arg.name," parameter must be a numeric value ",
                                "larger than ",arg.bounds[1]," and smaller ",
                                "than ",arg.bounds[2],"!")
                    },
                    botheq = {
                        if (arg.value<arg.bounds[1] || arg.value>arg.bounds[2])
                            stop(arg.name," parameter must be a numeric value ",
                                "larger than ",arg.bounds[1]," and smaller ",
                                "than ",arg.bounds[2],"!")
                    },
                    gt = {
                        if (arg.value<=arg.bounds[1])
                            stop(arg.name," parameter must be a numeric value ",
                                "greater than ",arg.bounds[1],"!")
                    },
                    lt = {
                        if (arg.value>=arg.bounds[1])
                            stop(arg.name," parameter must be a numeric value ",
                                "lower than ",arg.bounds[1],"!")
                    },
                    gte = {
                        if (arg.value<arg.bounds[1])
                            stop(arg.name," parameter must be a numeric value ",
                                "greater than or equal to ",arg.bounds[1],"!")
                    },
                    lte = {
                        if (arg.value>arg.bounds[1])
                            stop(arg.name," parameter must be a numeric value ",
                                "lower than or equal to ",arg.bounds[1],"!")
                    }
                )
            }
        },
        integer = {
            if (!is.integer(arg.value)) {
                # Can it be converted to integer?
                arg.value <- suppressWarnings(as.integer(arg.value))
                if (is.na(arg.value)) # No
                    stop(arg.name," parameter must be an integer value!")
            }
            if (!missing(arg.bounds)) {
                switch(direction,
                    both = {
                        if (arg.value<=arg.bounds[1] ||
                            arg.value>=arg.bounds[2])
                            stop(arg.name," parameter must be an integer ",
                                "larger than ",arg.bounds[1]," and smaller ",
                                "than ",arg.bounds[2],"!")
                    },
                    botheq = {
                        if (arg.value<arg.bounds[1] || arg.value>arg.bounds[2])
                            stop(arg.name," parameter must be an integer ",
                                "larger than or equal to ",arg.bounds[1],
                                " and smaller than or equal to ",arg.bounds[2],
                                "!")
                    },
                    gt = {
                        if (arg.value<=arg.bounds[1])
                            stop(arg.name," parameter must be an integer ",
                                "greater than ",arg.bounds[1],"!")
                    },
                    lt = {
                        if (arg.value>=arg.bounds[1])
                            stop(arg.name," parameter must be an integer ",
                                "lower than ",arg.bounds[1],"!")
                    },
                    gte = {
                        if (arg.value<arg.bounds[1])
                            stop(arg.name," parameter must be an integer ",
                                "greater than or equal to ",arg.bounds[1],"!")
                    },
                    lte = {
                        if (arg.value>arg.bounds[1])
                            stop(arg.name," parameter must be an integer ",
                                "lower than or equal to ",arg.bounds[1],"!")
                    }
                )
            }
        }
    )
}

checkFileArgs <- function(arg.name,arg.value) {
    if (!file.exists(arg.value))
        stop(arg.name," parameter must be an existing file!")
}

validateListArgs <- function(what,arg.list) {
    if (!(what %in% c("binParams","preprocessParams","selector",
        "strandedParams","saveParams","plotParams","kmParams","orderBy")))
        stop("Input list type for validation must be one of \"binParams\", ",
            "\"preprocessParams\", \"selector\", \"strandedParams\", ",
            "\"saveParams\", \"plotParams\", \"kmParams\", \"orderBy\"")
    if (!is.list(arg.list))
        stop(what," argument must be a list!")
    switch(what,
        orderBy = {
            valid.1 <- names(arg.list) %in% c("what","order","custom")
            not.valid.1 <- which(!valid.1)
            if (length(not.valid.1)>0) {
                warning("The following ",what," argument names are invalid ",
                    "and will be ignored: ",paste(names(arg.list)[not.valid.1],
                        collapse=", "),immediate.=TRUE)
                arg.list[not.valid.1] <- NULL
            }
            if (length(arg.list)>0) {
                for (n in names(arg.list)) {
                    switch(n,
                        what = {
                            arg.list$what <- tolower(arg.list$what[1])
                            if (length(grep("^(none|sum|max|avg|hc)",
                                arg.list$what,perl=TRUE))==0)
                                stop("The what option of orderBy parameter ",
                                    "must be one of \"none\", \"suma\", ",
                                    "\"sumn\", \"maxa\", \"maxn\", \"avga\",",
                                    "\"avgn\", \"hcn\".")
                        },
                        order = {
                            arg.list$order <- tolower(arg.list$order[1])
                            checkTextArgs("The order option of orderBy",
                                arg.list$order,c("descending","ascending"))
                        },
                        custom = {
                            if (!is.null(arg.list$custom) 
                                && !is.numeric(arg.list$custom))
                                stop("The custom option of orderBy parameter ",
                                    "must be a numeric vector or NULL!")
                        }
                    )
                }
            }
        },
        binParams = {
            valid.1 <- names(arg.list) %in% c("flankBinSize","regionBinSize",
                "sumStat","interpolation","binType","forceHeatmapBinning",
                "forcedBinSize","chunking")#,"seed")
            not.valid.1 <- which(!valid.1)
            if (length(not.valid.1)>0) {
                warning("The following ",what," argument names are invalid ",
                    "and will be ignored: ",paste(names(arg.list)[not.valid.1],
                        collapse=", "),immediate.=TRUE)
                arg.list[not.valid.1] <- NULL
            }
            if (length(arg.list)>0) {
                for (n in names(arg.list)) {
                    switch(n,
                        flankBinSize = {
                            checkNumArgs("The flankBinSize option of binParams",
                                arg.list$flankBinSize,"numeric",0,"gte")
                            
                        },
                        regionBinSize = {
                            checkNumArgs(
                                "The regionBinSize option of binParams",
                                arg.list$regionBinSize,"numeric",0,"gte"
                            )
                        },
                        sumStat = {
                            arg.list$sumStat <- tolower(arg.list$sumStat[1])
                            checkTextArgs("The sumStat option of binParams",
                                arg.list$sumStat,c("mean","median"))
                        },
                        interpolation = {
                            arg.list$interpolation <- 
                                tolower(arg.list$interpolation[1])
                            checkTextArgs(
                                "The interpoloation option of binParams",
                                arg.list$interpolation,
                                c("auto","spline","linear","neighborhood")
                            )
                        },
                        binType = {
                            arg.list$binType <- 
                                tolower(arg.list$binType[1])
                            checkTextArgs(
                                "The binType option of binParams",
                                arg.list$binType,c("variable","fixed")
                            )
                        },
                        forceHeatmapBinning = {
                            if (!is.logical(arg.list$forceHeatmapBinning))
                                stop("The forceHeatmapBinning option of ",
                                    "binParams parameter must be TRUE or ",
                                    "FALSE!")
                        },
                        forcedBinSize = {
                            checkNumArgs(
                                "The first forcedBinSize option of binParams",
                                arg.list$forcedBinSize[1],"numeric",0,"gte"
                            )
                            checkNumArgs(
                                "The second forcedBinSize option of binParams",
                                arg.list$forcedBinSize[2],"numeric",0,"gte"
                            )
                        },
                        chunking = {
                            if (!is.logical(arg.list$chunking))
                                stop("The chunking option of binParams ",
                                    "parameter must be TRUE or FALSE!")
                        }#,
                        #seed = {
                        #    checkNumArgs(
                        #        "The seed option of binParams",arg.list$seed,
                        #        "numeric"
                        #    )
                        #}
                    )
                }
            }
        },
        selector = {
            valid.1 <- names(arg.list) %in% c("id","biotype","exonType")
            not.valid.1 <- which(!valid.1)
            if (length(not.valid.1)>0) {
                warning("The following ",what," argument names are invalid ",
                    "and will be ignored: ",paste(names(arg.list)[not.valid.1],
                        collapse=", "),immediate.=TRUE)
                arg.list[not.valid.1] <- NULL
            }
            # selector members cannot be really checked.In absolute custom cases
            # id can be whatever, text, number. Also biotype requires additional
            # argument organism to this function and also exonType. Thus, the
            # program might crash elsewhere.
        },
        preprocessParams = {
            valid.1 <- names(arg.list) %in% c("fragLen","cleanLevel",
                "normalize","sampleTo","spliceAction","spliceRemoveQ",
                "bedGenome")#,"seed")
            not.valid.1 <- which(!valid.1)
            if (length(not.valid.1)>0) {
                warning("The following ",what," argument names are invalid ",
                    "and will be ignored: ",paste(names(arg.list)[not.valid.1],
                        collapse=", "),immediate.=TRUE)
                arg.list[not.valid.1] <- NULL
            }
            if (length(arg.list)>0) {
                for (n in names(arg.list)) {
                    switch(n,
                        fragLen = {
                            if (!is.na(arg.list$fragLen))
                                checkNumArgs(
                                    "The fragLen option of preprocessParams",
                                    arg.list$fragLen,"numeric",0,"gt"
                                )
                        },
                        cleanLevel = {
                            arg.list$cleanLevel <- 
                                as.integer(arg.list$cleanLevel[1])
                            checkNumArgs(
                                "The cleanLevel option of preprocessParams",
                                arg.list$cleanLevel,"integer",c(0,3),"botheq"
                            )
                        },
                        normalize = {
                            arg.list$normalize <- tolower(arg.list$normalize[1])
                            checkTextArgs(
                                "The normalize option of preprocessParams",
                                arg.list$normalize,
                                c("none","linear","downsample","sampleto")
                            )
                            
                        },
                        sampleTo = {
                            checkNumArgs(
                                "The sampleTo option of preprocessParams",
                                arg.list$sampleTo,"numeric",1e+5,"gte"
                            )
                        },
                        spliceAction = {
                            arg.list$spliceAction <- 
                                tolower(arg.list$spliceAction[1])
                            checkTextArgs(
                                "The spliceAction option of preprocessParams",
                                arg.list$spliceAction,c("keep","remove","split")
                            )
                        },
                        spliceRemoveQ = {
                            checkNumArgs(
                                "The spliceRemoveQ option of preprocessParams",
                                arg.list$spliceRemoveQ,"numeric",c(0,1),"both"
                            )
                        },
                        bedGenome = {
                            if (!is.na(arg.list$bedGenome))
                                checkTextArgs(
                                    "The bedGenome option of preprocessParams",
                                    arg.list$bedGenome,c("hg18","hg19","hg38",
                                        "mm9","mm10","rn5","dm3","danrer7",
                                        "pantro4","susscr3")
                                )
                        }#,
                        #seed = {
                        #    if (!is.numeric(arg.list$seed))
                        #        stop("The seed option of preprocessParams ",
                        #            "parameter must be numeric!")
                        #}
                    )
                }
            }
        },
        plotParams = {
            valid.1 <- names(arg.list) %in% c("plot","profile","heatmap",
                "correlation","device","signalScale","heatmapScale",
                "heatmapFactor","corrScale","sumStat","smooth","corrSmoothPar",
                "singleFacet","multiFacet","singleFacetDirection","conf",
                "outputDir","outputBase")
            not.valid.1 <- which(!valid.1)
            if (length(not.valid.1)>0) {
                warning("The following ",what," argument names are invalid ",
                    "and will be ignored: ",paste(names(arg.list)[not.valid.1],
                        collapse=", "),immediate.=TRUE)
                arg.list[not.valid.1] <- NULL
            }
            if (length(arg.list)>0) {
                for (n in names(arg.list)) {
                    switch(n,
                        plot = {
                            if (!is.logical(arg.list$plot))
                                stop("The plot option of plotParams parameter ",
                                    "must be TRUE or FALSE!")
                        },
                        profile = {
                            if (!is.logical(arg.list$profile))
                                stop("The profile option of plotParams ",
                                    "parameter must be TRUE or FALSE!")
                        },
                        heatmap = {
                            if (!is.logical(arg.list$heatmap))
                                stop("The heatmap option of plotParams ",
                                    "parameter must be TRUE or FALSE!")
                        },
                        correlation = {
                            if (!is.logical(arg.list$correlation))
                                stop("The correlation option of plotParams ",
                                    "parameter must be TRUE or FALSE!")
                        },
                        conf = {
                            if (!is.logical(arg.list$conf))
                                stop("The conf option of plotParams parameter ",
                                    "must be TRUE or FALSE!")
                        },
                        device = {
                            arg.list$device <- tolower(arg.list$device[1])
                            checkTextArgs(
                                "The device option of plotParams",
                                arg.list$device,
                                c("x11","png","jpg","tiff","bmp","pdf","ps")
                            )
                        },
                        signalScale = {
                            arg.list$signalScale <- 
                                tolower(arg.list$signalScale[1])
                            checkTextArgs(
                                "The signalScale option of plotParams",
                                arg.list$signalScale,
                                c("natural","log2")
                            )
                        },
                        heatmapScale = {
                            arg.list$heatmapScale <- 
                                tolower(arg.list$heatmapScale[1])
                            checkTextArgs(
                                "The signalScale option of plotParams",
                                arg.list$heatmapScale,
                                c("each","common")
                            )
                        },
                        heatmapFactor = {
                            checkNumArgs(
                                "The heatmapFactor option of plotParams",
                                arg.list$heatmapFactor,"numeric",0,"gt"
                            )
                        },
                        corrScale = {
                            arg.list$corrScale <- 
                                tolower(arg.list$corrScale[1])
                            checkTextArgs(
                                "The corrScale option of plotParams",
                                arg.list$corrScale,
                                c("normalized","each")
                            )
                        },
                        sumStat = {
                            arg.list$sumStat <- tolower(arg.list$sumStat[1])
                            checkTextArgs("The sumStat option of plotParams",
                                arg.list$sumStat,c("mean","median"))
                        },
                        smooth = {
                            if (!is.logical(arg.list$smooth))
                                stop("The smooth option of plotParams ",
                                    "parameter must be TRUE or FALSE!")
                        },
                        corrSmoothPar = {
                            checkNumArgs(
                                "The corrSmoothPar option of plotParams",
                                arg.list$corrSmoothPar,"numeric",c(0,1),"both"
                            )
                        },
                        singleFacet = {
                            arg.list$singleFacet <- 
                                tolower(arg.list$singleFacet[1])
                            checkTextArgs(
                                "The facet option of plotParams",
                                arg.list$singleFacet,
                                c("none","wrap","grid")
                            )
                        },
                        multiFacet = {
                            arg.list$multiFacet <- 
                                tolower(arg.list$multiFacet[1])
                            checkTextArgs(
                                "The facet option of plotParams",
                                arg.list$multiFacet,
                                c("wrap","grid")
                            )
                        },
                        singleFacetDirection = {
                            arg.list$singleFacetDirection <- 
                                tolower(arg.list$singleFacetDirection[1])
                            checkTextArgs(
                                "The facet option of plotParams",
                                arg.list$singleFacetDirection,
                                c("horizontal","vertical")
                            )
                        },
                        outputDir = {
                            if (!dir.exists(arg.list$outputDir))
                                stop("The outputDir option of plotParams must ",
                                    "be an existing directory!")
                        },
                        outputBase = {
                            if (!is.null(arg.list$outputBase)
                                && !is.character(arg.list$outputBase))
                                stop("The outputBase option of plotParams ",
                                    "parameter must be either NULL for ",
                                    "automatic filename generation or a ",
                                    "string!")
                        }
                    )
                }
            }
        },
        saveParams = {
            valid.1 <- names(arg.list) %in% c("ranges","coverage","profile",
                "profilePlot","heatmapPlot","correlationPlot")
            not.valid.1 <- which(!valid.1)
            if (length(not.valid.1)>0) {
                warning("The following ",what," argument names are invalid ",
                    "and will be ignored: ",paste(names(arg.list)[not.valid.1],
                        collapse=", "),immediate.=TRUE)
                arg.list[not.valid.1] <- NULL
            }
            if (length(arg.list)>0) {
                for (n in names(arg.list)) {
                    switch(n,
                        ranges = {
                            if (!is.logical(arg.list$ranges))
                                stop("The ranges option of saveParams ",
                                    "parameter must be TRUE or FALSE!")
                        },
                        coverage = {
                            if (!is.logical(arg.list$coverage))
                                stop("The coverage option of saveParams ",
                                    "parameter must be TRUE or FALSE!")
                        },
                        profile = {
                            if (!is.logical(arg.list$profile))
                                stop("The profile option of saveParams ",
                                    "parameter must be TRUE or FALSE!")
                        },
                        profilePlot = {
                            if (!is.logical(arg.list$profilePlot))
                                stop("The profilePlot option of saveParams ",
                                    "parameter must be TRUE or FALSE!")
                        },
                        heatmapPlot = {
                            if (!is.logical(arg.list$heatmapPlot))
                                stop("The heatmapPlot option of saveParams ",
                                    "parameter must be TRUE or FALSE!")
                        },
                        correlationPlot = {
                            if (!is.logical(arg.list$correlationPlot))
                                stop("The correlationPlot option of ",
                                    "saveParams parameter must be TRUE or ",
                                    "FALSE!")
                        }
                    )
                }
            }
        },
        kmParams = {
            valid.1 <- names(arg.list) %in% c("k","nstart","algorithm",
                "reference","iterMax")#,"seed")
            not.valid.1 <- which(!valid.1)
            if (length(not.valid.1)>0) {
                warning("The following ",what," argument names are invalid ",
                    "and will be ignored: ",paste(names(arg.list)[not.valid.1],
                        collapse=", "),immediate.=TRUE)
                arg.list[not.valid.1] <- NULL
            }
            if (length(arg.list)>0) {
                for (n in names(arg.list)) {
                    switch(n,
                        k = {
                            checkNumArgs("The k option of kmParams",arg.list$k,
                            "numeric",0,"gte")
                        },
                        nstart = {
                            checkNumArgs("The nstart option of kmParams",
                                arg.list$nstart,"numeric",0,"gt")
                        },
                        algorithm = {
                            arg.list$algorithm <- arg.list$algorithm[1]
                            checkTextArgs(
                                "The algorithm option of kmParams",
                                arg.list$algorithm,
                                c("Hartigan-Wong","Lloyd","Forgy","MacQueen")
                            )
                        },
                        reference = {
                            if (!is.null(arg.list$reference) 
                                && !is.character(arg.list$reference))
                                stop("The reference option of kmParams must ",
                                    "be either NULL or on of the sample ids!")
                        },
                        iterMax = {
                            checkNumArgs("The iterMax option of kmParams",
                                arg.list$iterMax,"numeric",0,"gt")
                        }#,
                        #seed = {
                        #    if (!is.numeric(arg.list$seed))
                        #        stop("The seed option of kmParams parameter ",
                        #            "must be numeric!")
                        #}
                    )
                }
            }
        },
        strandedParams = {
            valid.1 <- names(arg.list) %in% c("strand","ignoreStrand")
            not.valid.1 <- which(!valid.1)
            if (length(not.valid.1)>0) {
                warning("The following ",what," argument names are invalid ",
                    "and will be ignored: ",paste(names(arg.list)[not.valid.1],
                        collapse=", "),immediate.=TRUE)
                arg.list[not.valid.1] <- NULL
            }
            if (length(arg.list)>0) {
                for (n in names(arg.list)) {
                    switch(n,
                        strand = {
                            if (!is.null(arg.list$strand))
                                checkTextArgs(paste0("The strand option of ",
                                    "strandParams "),arg.list$strand,c("+","-"))
                        },
                        ignoreStrand = {
                            if (!is.logical(arg.list$ignoreStrand))
                                stop("The ignoreStrand option of strandParams ",
                                    "parameter must be TRUE or FALSE!")
                        }
                    )
                }
            }
        }
    )
    return(arg.list)
}

checkInput <- function(input) {
    # Input must have  id, file and format
    for (i in seq_len(length(input))) {
        if (is.null(input[[i]]$id))
            stop("All input list members must have an id field! Member ",i,
                " does not have one.")
        if (is.null(input[[i]]$file))
            stop("All input list members must have a file field! Member ",i,
                " does not have one.")
        if (is.null(input[[i]]$format))
            stop("All input list members must have a format field! Member ",i,
                " does not have one.")
    }
}

validateOrganismArg <- function(org) {
    org <- tolower(org[1])
    checkTextArgs("org",org,getSupportedOrganisms(),multiarg=FALSE)
}
