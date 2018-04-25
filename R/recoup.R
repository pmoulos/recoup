recoup <- function(
    input,
    design=NULL,
    region=c("genebody","tss","tes","custom"),
    type=c("chipseq","rnaseq"),
    genome=c("hg18","hg19","hg38","mm9","mm10","rn5","rn6","dm3","dm6",
        "danrer7","danrer10","pantro4","pantro5","susscr3","susscr11",
        "ecucab2","tair10"),
    version="auto",
    refdb=c("ensembl","ucsc","refseq"),
    flank=c(2000,2000),
    fraction=1,
    orderBy=list(
        what=c("none","suma","sumn","maxa","maxn","avga","avgn","hcn"),
        order=c("descending","ascending"),
        custom=NULL
    ),
    #orderBy=getDefaultListArgs("orderBy"),
    binParams=list(
        flankBinSize=0,
        regionBinSize=0,
        sumStat=c("mean","median"),
        interpolation=c("auto","spline","linear","neighborhood"),
        forceHeatmapBinning=TRUE,
        forcedBinSize=c(50,200),
        chunking=FALSE
    ),
    #binParams=getDefaultListArgs("binParams"),
    selector=NULL,
    #selector=getDefaultListArgs("selector"),
    preprocessParams=list(
        fragLen=NA,
        cleanLevel=c(0,1,2,3),
        normalize=c("none","linear","downsample","sampleto"),
        sampleTo=1e+6,
        spliceAction=c("split","keep","remove"),
        spliceRemoveQ=0.75,
        #bedGenome=ifelse(genome %in% c("hg18","hg19","hg38","mm9","mm10","rn5",
        #    "dm3","danrer7","pantro4","susscr3"),genome,NA),
        bedGenome=NA,
        seed=42
    ),
    #preprocessParams=getDefaultListArgs("preprocessParams"),
    plotParams=list(
        plot=TRUE,
        profile=TRUE,
        heatmap=TRUE,
        correlation=TRUE,
        signalScale=c("natural","log2"),
        heatmapScale=c("common","each"),
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
    ),
    #plotParams=getDefaultListArgs("plotParams",design),
    saveParams=list(
        ranges=TRUE,
        coverage=TRUE,
        profile=TRUE,
        profilePlot=TRUE,
        heatmapPlot=TRUE,
        correlationPlot=TRUE
    ),
    #saveParams=getDefaultListArgs("saveParams"),
    kmParams=list(
        k=0, # Do not perform kmeans
        nstart=20,
        algorithm=c("Hartigan-Wong","Lloyd","Forgy","MacQueen"),
        iterMax=20,
        reference=NULL,
        seed=42
    ),
    #kmParams=getDefaultListArgs("kmParams"),
    strandedParams=list(
        strand=NULL,
        ignoreStrand=TRUE
    ),
    #strandedParams=getDefaultListArgs("strandedParams"),
    ggplotParams=list(
        title=element_text(size=12),
        axis.title.x=element_text(size=10,face="bold"),
        axis.title.y=element_text(size=10,face="bold"),
        axis.text.x=element_text(size=9,face="bold"),
        axis.text.y=element_text(size=10,face="bold"),
        strip.text.x=element_text(size=10,face="bold"),
        strip.text.y=element_text(size=10,face="bold"),
        legend.position="bottom",
        panel.margin=grid::unit(1,"lines")
    ),
    #ggplotParams=getDefaultListArgs("ggplotParams"),
    complexHeatmapParams=list(
        main=list(
            cluster_rows=ifelse(length(grep("hc",orderBy$what))>0,TRUE,FALSE),
            cluster_columns=FALSE,
            column_title_gp=grid::gpar(fontsize=10,font=2),
            show_row_names=FALSE,
            show_column_names=FALSE,
            heatmap_legend_param=list(
                color_bar="continuous"
            )
        ),
        group=list(
            cluster_rows=ifelse(length(grep("hc",orderBy$what))>0,TRUE,FALSE),
            cluster_columns=FALSE,
            column_title_gp=grid::gpar(fontsize=10,font=2),
            show_row_names=FALSE,
            show_column_names=FALSE,
            row_title_gp=grid::gpar(fontsize=8,font=2),
            gap=unit(5,"mm"),
            heatmap_legend_param=list(
                color_bar="continuous"
            )
        )
    ),
    #complexHeatmapParams=getDefaultListArgs("complexHeatmapParams"),
    bamParams=NULL,
    onTheFly=FALSE, # Directly from BAM w/o Ranges storing, also N/A with BED,
    localDbHome=file.path(path.expand("~"),".recoup"),
    rc=NULL
) {
    ############################################################################
    # Begin of simple parameter checking for a new or a restored object
    if (is.list(input) && !is.null(input$data)) { # Refeeding recoup object
        message("Detected a previous recoup output as input. Any existing ",
            "data requiring\nlengthy calculations (short reads, coverages and ",
            "profile matrices will be\nreused depending on other parameters. ",
            "If you want a complete recalculation,\neither input a fresh list ",
            "of arguments or remove any of the above with the\nremoveData ",
            "function.\n")
        prevCallParams <- input$callopts
        input <- input$data
    }
    else
        prevCallParams <- NULL
    if (!is.list(input) && file.exists(input))
        input <- readConfig(input)
    checkInput(input)
    if (is.null(names(input)))
        names(input) <- sapply(input,function(x) return(x$id))
    
    # Check if there are any mispelled or invalid parameters and throw a warning
    checkMainArgs(as.list(match.call()))
    
    # Check rest of arguments
    region <- tolower(region[1])
    refdb <- tolower(refdb[1])
    type <- tolower(type[1])
    if (!is.null(design) && !is.data.frame(design))
        checkFileArgs("design",design)
    if (!is.data.frame(genome) && file.exists(genome))
        checkFileArgs("genome",genome)
    else if (is.character(genome))
        checkTextArgs("genome",genome,getSupportedOrganisms(),multiarg=FALSE)
    checkTextArgs("refdb",refdb,getSupportedRefDbs(),multiarg=FALSE)
    checkTextArgs("type",type,c("chipseq","rnaseq"),multiarg=FALSE)
    checkNumArgs("fraction",fraction,"numeric",c(0,1),"botheq")
    if (any(flank<0))
        stop("The minimum flanking allowed is 0 bp")
    if (any(flank>50000))
        stop("The maximum flanking allowed is 50000 bp")
    # Check the version argument
    version <- version[1]
    if (is.character(version)) {
        version <- tolower(version)
        checkTextArgs("version",version,c("auto"),multiarg=FALSE)
    }
    else
        checkNumArgs("version",version,"numeric")
    # If type is rnaseq, only genebody plots are valid
    if (type=="rnaseq" && region!="genebody") {
        warning("When type is \"rnaseq\", plots can be created only on ",
            "genebodies! Switching to genebody regions...",immediate.=TRUE)
        region <- "genebody"
    }
    
    # If type is rnaseq, read extension is illegal because of potential splicing
    if (type=="rnaseq" && !is.null(preprocessParams$fragLen)
        && !is.na(preprocessParams$fragLen)) {
        warning("When type is \"rnaseq\", read extension/reduction should not ",
            "happen because of potential splicing! Ignoring...",immediate.=TRUE)
        preprocessParams$fragLen <- NA
    }
    
    # and list arguments
    orderByDefault <- getDefaultListArgs("orderBy")
    binParamsDefault <- getDefaultListArgs("binParams")
    selectorDefault <- getDefaultListArgs("selector")
    preprocessParamsDefault <- getDefaultListArgs("preprocessParams",
        genome=genome)
    plotParamsDefault <- getDefaultListArgs("plotParams",design=design)
    saveParamsDefault <- getDefaultListArgs("saveParams")
    kmParamsDefault <- getDefaultListArgs("kmParams")
    strandedParamsDefault <- getDefaultListArgs("strandedParams")
    
    orderBy <- setArg(orderByDefault,orderBy)
    binParams <- setArg(binParamsDefault,binParams)
    if (!is.null(selector))
        selector <- setArg(selectorDefault,selector)
    preprocessParams <- setArg(preprocessParamsDefault,preprocessParams)
    plotParams <- setArg(plotParamsDefault,plotParams)
    saveParams <- setArg(saveParamsDefault,saveParams)
    kmParams <- setArg(kmParamsDefault,kmParams)
    strandedParams <- setArg(strandedParamsDefault,strandedParams)
    
    orderBy <- validateListArgs("orderBy",orderBy)
    binParams <- validateListArgs("binParams",binParams)
    if (!is.null(selector))
        selector <- validateListArgs("selector",selector)
    preprocessParams <- validateListArgs("preprocessParams",preprocessParams)
    plotParams <- validateListArgs("plotParams",plotParams)
    saveParams <- validateListArgs("saveParams",saveParams)
    kmParams <- validateListArgs("kmParams",kmParams)
    strandedParams <- validateListArgs("strandedParams",strandedParams)
    
    if (is.null(plotParams$outputBase))
        plotParams$outputBase <- paste(sapply(input,function(x) return(x$id)),
            collapse="-")
    
    if (any(sapply(input,function(x) return(tolower(x$format)=="bed")))
        && !(preprocessParams$bedGenome %in% c("hg18","hg19","hg38","mm9",
            "mm10","rn5","dm3","danrer7","pantro4","susscr3")))
        stop("When short read files are in BED format, either the genome ",
            "parameter should be one\nof the supported organisms, or ",
            "preprocessParams$bedGenome must be specified as one of them.")
   
    # End of simple parameter checking for a new or a restored object
    ############################################################################
    
    ############################################################################
    # Begin of complex parameter storage procedure and parameter recall from a
    # previous call
    
    # Store all parameters to an obect for later reference to a next call, after
    # checking if something has changed (e.g. using setr)...
    thisCall <- as.list(match.call())[-1]
    if (!is.null(prevCallParams)) {
        callParams <- list(
            region=if (is.null(thisCall$region)) prevCallParams$region else
                region,
            type=if (is.null(thisCall$type)) prevCallParams$type else type,
            genome=if (is.null(thisCall$genome)) prevCallParams$genome else
                genome,
            refdb=if (is.null(thisCall$refdb)) prevCallParams$refdb else
                refdb,
            version=if (is.null(thisCall$version)) prevCallParams$version else
                version,
            flank=if (is.null(thisCall$flank)) prevCallParams$flank else flank,
            fraction=if (is.null(thisCall$fraction)) prevCallParams$fraction
                else fraction,
            orderBy=if (is.null(thisCall$orderBy)) prevCallParams$orderBy else
                orderBy,
            binParams=if (is.null(thisCall$binParams)) prevCallParams$binParams
                else binParams,
            selector=selector, # selector is a special case
            preprocessParams=if (is.null(thisCall$preprocessParams))
                prevCallParams$preprocessParams else preprocessParams,
            plotParams=if (is.null(thisCall$plotParams)) 
                prevCallParams$plotParams else plotParams,
            saveParams=if (is.null(thisCall$saveParams))
                prevCallParams$saveParams else saveParams,
            kmParams=if (is.null(thisCall$kmParams))
                prevCallParams$kmParams else kmParams,
            strandedParams=if (is.null(thisCall$strandedParams))
                prevCallParams$strandedParams else strandedParams,
            ggplotParams=if (is.null(thisCall$ggplotParams))
                prevCallParams$ggplotParams else ggplotParams,
            complexHeatmapParams=if (is.null(thisCall$complexHeatmapParams))
                prevCallParams$complexHeatmapParams else complexHeatmapParams,
            #bamParams=if (is.null(thisCall$bamParams),
            #    prevCallParams$bamParams,bamParams), ## Unused
            onTheFly=if (is.null(thisCall$onTheFly)) prevCallParams$onTheFly 
                else onTheFly,
            localDbHome=if (is.null(thisCall$localDbHome)) 
                prevCallParams$localDbHome else localDbHome,
            rc=if (is.null(thisCall$rc)) prevCallParams$rc else rc,
            customIsBase=NULL # Additional placeholder
        )
    }
    else {
        callParams <- list(
            region=region,
            type=type,
            genome=genome,
            version=version,
            refdb=refdb,
            flank=flank,
            fraction=fraction,
            orderBy=orderBy,
            binParams=binParams,
            selector=selector,
            preprocessParams=preprocessParams,
            plotParams=plotParams,
            saveParams=saveParams,
            kmParams=kmParams,
            strandedParams=strandedParams,
            ggplotParams=ggplotParams,
            complexHeatmapParams=complexHeatmapParams,
            bamParams=bamParams,
            onTheFly=onTheFly,
            localDbHome=localDbHome,
            rc=rc,
            customIsBase=NULL # Additional placeholder
        )
    }
    
    # ...and check if there has been a previous call and decide on big
    # recalculations...
    input <- decideChanges(input,callParams,prevCallParams)
    
    # Redefine all final arguments for this call
    region=callParams$region
    type=callParams$type
    genome=callParams$genome
    version=callParams$version
    refdb=callParams$refdb
    flank=callParams$flank
    fraction=callParams$fraction
    orderBy=callParams$orderBy
    binParams=callParams$binParams
    selector=callParams$selector
    preprocessParams=callParams$preprocessParams
    plotParams=callParams$plotParams
    saveParams=callParams$saveParams
    kmParams=callParams$kmParams
    strandedParams=callParams$strandedParams
    ggplotParams=callParams$ggplotParams
    complexHeatmapParams=callParams$complexHeatmapParams
    bamParams=callParams$bamParams
    onTheFly=callParams$onTheFly
    localDbHome=callParams$localDbHome
    rc=callParams$rc
    # End of complex parameter storage procedure and parameter recall from a
    # previous call
    ############################################################################
    
    # Continue with actual work
    if (!is.data.frame(genome)) {
        if (file.exists(genome)) {
            genome <- read.delim(genome) # Must be bed like
            rownames(genome) <- as.character(genome[,4])
            genomeRanges <- makeGRangesFromDataFrame(
                df=genome,
                keep.extra.columns=TRUE
            )
        }
        else {
            # Check if local storage has been set
            storageHome <- file.path(localDbHome,refdb,genome)
            if (dir.exists(storageHome)) {
                # Some compatibility with previous annoation stores
                if ("gene.rda" %in% dir(storageHome)) { # Older version
                    warning("Older annotation storage detected! Please ",
                        "rebuild your annotation storage using the ",
                        "buildAnnotationStore function\nso as to have the ",
                        "ability of versioning genomic coordinate annotations.",
                        "\nIn the future, older versions may become unusable.",
                        immediate.=TRUE)
                    storageHomeVersion <- storageHome
                }
                else { # Newer version
                    # Decide on version
                    if (version != "auto") {
                        storageHomeVersion <- file.path(storageHome,version)
                        if (!dir.exists(storageHomeVersion)) {
                            warning("The annotation directory ",
                                storageHomeVersion,
                                " does not seem to exist! Have you run ",
                                "buildAnnotationStorage? Will use newest ",
                                "existing version.",immediate.=TRUE)
                            version <- "auto"
                        }
                    }
                    if (version == "auto") {
                        vers <- suppressWarnings(as.numeric(dir(storageHome)))
                        if (any(is.na(vers))) {
                            of <- vers[which(is.na(vers))]
                            stop("Corrupted annotation storage directory ",
                                "contents! ->",of,"<- Either remove offending ",
                                "files/directories or rebuild.")
                        }
                        vers <- sort(vers,decreasing=TRUE)
                        version <- vers[1]
                        storageHomeVersion <- file.path(storageHome,version)
                    }
                }
                
                if (type=="chipseq") {
                    g <- load(file.path(storageHomeVersion,"gene.rda"))
                    genomeRanges <- gene
                    helperRanges <- NULL
                }
                else if (type=="rnaseq") {
                    # Load the helper ranges anyway
                    g <- load(file.path(storageHomeVersion,"gene.rda"))
                    helperRanges <- gene
                    if (all(flank==0)) {
                        g <- load(file.path(storageHomeVersion,
                            "summarized_exon.rda"))
                        genomeRanges <- sexon
                    }
                    else if (all(flank==500)) {
                        g <- load(file.path(storageHomeVersion,
                            "summarized_exon_F500.rda"))
                        genomeRanges <- flankedSexon
                    }
                    else if (all(flank==1000)) {
                        g <- load(file.path(storageHomeVersion,
                            "summarized_exon_F1000.rda"))
                        genomeRanges <- flankedSexon
                    }
                    else if (all(flank==2000)) {
                        g <- load(file.path(storageHomeVersion,
                            "summarized_exon_F2000.rda"))
                        genomeRanges <- flankedSexon
                    }
                    else if (all(flank==5000)) {
                        g <- load(file.path(storageHomeVersion,
                            "summarized_exon_F2000.rda"))
                        genomeRanges <- flankedSexon
                    }
                    else {
                        warning("When using recoup in RNA-Seq mode, it is ",
                            "much faster to use a precalculated set of ",
                            "regions and their flanks (0, 500, 1000, 2000 and ",
                            "5000 bps supported) and subset this one for a ",
                            "custom set of genes. Otherwise calculations may ",
                            "take too long to complete.",immediate.=TRUE)
                        genomeRanges <- 
                            getMainRnaRangesOnTheFly(helperRanges,flank,rc=rc)
                    }
                }
            }
            else { # On the fly
                message("Getting annotation on the fly for ",genome," from ",
                    refdb)
                if (type=="chipseq") {
                    genome <- getAnnotation(genome,"gene",refdb=refdb,rc=rc)
                    helperRanges <- NULL
                }
                else if (type=="rnaseq") {
                    helper <- getAnnotation(genome,"gene",refdb=refdb,rc=rc)
                    helperRanges <- makeGRangesFromDataFrame(
                        df=helper,
                        keep.extra.columns=TRUE,
                        seqnames.field="chromosome"
                    )
                    names(helperRanges) <- as.character(helperRanges$gene_id)
                    genome <- getAnnotation(genome,"exon",refdb=refdb,rc=rc)
                    annGr <- makeGRangesFromDataFrame(
                        df=genome,
                        keep.extra.columns=TRUE,
                        seqnames.field="chromosome"
                    )
                    message("Merging exons")
                    annGr <- reduceExons(annGr,rc=rc)
                    names(annGr) <- as.character(annGr$exon_id)
                    genomeRanges <- split(annGr,annGr$gene_id)
                }
            }
        }
    }
    else {
        rownames(genome) <- as.character(genome[,4])
        genomeRanges <- makeGRangesFromDataFrame(
            df=genome,
            keep.extra.columns=TRUE
        )
    }
    
    # After genome read, check if we have a valid custom orderer
    if (!is.null(orderBy$custom)) {
        if (length(orderBy$custom) != nrow(genome)) {
            warning("The custom orderer provide with orderBy parameter does ",
                "not have length equal to the number of elements in genome ",
                "argument and will be ignored!",immediate.=TRUE)
            orderBy$custom <- NULL
        }
    }
    
    # Read and check design compatibilities. Check if k-means is requested and
    # message accordingly. If k-means is requested it will be added to the 
    # design data frame
    hasDesign <- FALSE
    if (!is.null(design)) {
        if (!is.data.frame(design))
            design <- read.delim(design,row.names=1)
        nfac <- ncol(design)
        if (length(input)>1 && nfac>2)
            stop("When more than one files are provided for coverage ",
                "generation, the maximum number of allowed design factors is 2")
        if (length(input)>1 && nfac>1 && kmParams$k>0)
            stop("When more than one files are provided for coverage ",
                "generation and k-means clustering is also requested, the ",
                "maximum number of allowed design factors is 1")
        if (length(input)==1 && nfac>3)
            stop("The maximum number of allowed design factors is 3")
        if (length(input)==1 && nfac>2 && kmParams$k>0)
            stop("The maximum number of allowed design factors when k-means ",
                "clustering is requested is 2")
        if (length(input)==1 && nfac>2 && plotParams$singleFacet!="none")
            stop("The maximum number of allowed design factors whith one ",
                "sample but wishing a gridded profile layout is 2")
        if (length(input)==1 && nfac>1 && kmParams$k>0 
            && plotParams$singleFacet!="none")
            stop("The maximum number of allowed design factors whith one ",
                "sample but wishing for k-means clustering and gridded ",
                "profile layout is 1")
        # Reduce the genomeRanges according to design or the other way
        if (nrow(design)>length(genomeRanges))
            design <- tryCatch({
                design[names(genomeRanges),,drop=FALSE]
            },error=function(e) {
                stop("Unexpected error occured! Are you sure that element ",
                    "(row) names in the design file are of the same type as ",
                    "the genome file?")
            },finally={})
        else if (nrow(design)<=length(genomeRanges)) {
            genomeRanges <- tryCatch({
                genomeRanges[rownames(design)]
            },error=function(e) {
                stop("Unexpected error occured! Are you sure that element ",
                    "(row) names in the design file are of the same type as ",
                    "the genome file?")
            },finally={})
            if (type=="rnaseq")
                helperRanges <- tryCatch({
                    helperRanges[rownames(design)]
                },error=function(e) {
                    stop("Unexpected error occured! Are you sure that element ",
                        "(row) names in the design file are of the same type as ",
                        "the genome file?")
                },finally={})
        }
        # ...but maybe the names are incompatible
        if (length(genomeRanges)==0)
            stop("No ranges left after using the identifiers provided with ",
                "the design file. Are you sure the identifiers between the ",
                "two files are compatible?")
        if (nrow(design)==0)
            stop("No design elements left after using the identifiers ",
                "provided with the genome file. Are you sure the identifiers ",
                "between the two files are compatible?")
    }
    
    # Apply the rest of the filters if any to reduce later computational burden
    if (!is.null(selector)) {
        if (type=="chipseq")
            genomeRanges <- applySelectors(genomeRanges,selector,rc=rc)
        if (type=="rnaseq") {
            helperRanges <- applySelectors(helperRanges,selector)
            # Filter and align names if we have helperRanges around
            genomeRanges <- genomeRanges[names(helperRanges)]
        }
    }
    
    # Now we must follow two paths according to region type, genebody and custom
    # areas with equal/unequal widths, or tss, tes and 1-width custom areas.
    callParams$customIsBase <- FALSE
    if (region=="custom" && all(width(genomeRanges)==1)) {
        if (all(flank==0)) {
            warning("Flanking cannot be zero bp in both directions when the ",
                "reference region is only 1bp! Setting to default ",
                "(2000,2000)...",immediate.=TRUE)
            flank <- c(2000,2000)
            callParams$flank <- flank
        }
        callParams$customIsBase <- TRUE
    }
    if (region %in% c("tss","tes")) {
        if (all(flank==0)) {
            warning("Flanking cannot be zero bp in both directions when the ",
                "reference region is \"tss\" or \"tes\"! Setting to default ",
                "(2000,2000)...", immediate.=TRUE)
            flank <- c(2000,2000)
            callParams$flank <- flank
        }
    }
    
    # Here we must determine the final ranges to pass to preprocessRanges and
    # then work with these later on
    intermRanges <- getMainRanges(genomeRanges,helperRanges=helperRanges,type,
        region,flank,rc=rc)
    mainRanges <- intermRanges$mainRanges
    bamRanges <- intermRanges$bamRanges

    # Here we must write code for the reading and normalization of bam files
    # The preprocessRanges function looks if there is a valid (not null) ranges
    # field in input
    #if (!onTheFly)
    #input <- preprocessRanges(input,preprocessParams,bamParams=bamParams,rc=rc)
    input <- preprocessRanges(input,preprocessParams,genome,bamRanges,
        bamParams=bamParams,rc=rc)
    
    # At this point we must apply the fraction parameter if <1. We choose this
    # point in order not to restrict later usage of the read ranges and since it
    # does not take much time to load into memory.
    #if (fraction<1) {
    #    newSize <- round(fraction*length(genomeRanges))
    #    set.seed(preprocessParams$seed)
    #    refIndex <- sort(sample(length(genomeRanges),newSize))
    #    genomeRanges <- genomeRanges[refIndex]
    #    if (type=="rnaseq")
    #        helperRanges <- helperRanges[refIndex]
    #    for (i in 1:length(input)) {
    #        if (!is.null(input[[i]]$ranges)) {
    #            newSize <- round(fraction*length(input[[i]]$ranges))
    #            set.seed(preprocessParams$seed)
    #            fracIndex <- sort(sample(length(input[[i]]$ranges),newSize))
    #            input[[i]]$ranges <- input[[i]]$ranges[fracIndex]
    #        }
    #        if (!is.null(input[[i]]$coverage)) {
    #            #if (!is.null(input[[i]]$coverage$center)) {
    #            #    input[[i]]$coverage$center <- 
    #            #        input[[i]]$coverage$center[refIndex]
    #            #}
    #            #else
    #            input[[i]]$coverage <- input[[i]]$coverage[refIndex]
    #        }
    #        if (!is.null(input[[i]]$profile))
    #            input[[i]]$profile <- input[[i]]$profile[refIndex,]
    #    }
    #}
    if (fraction<1) {
        newSize <- round(fraction*length(mainRanges))
        set.seed(preprocessParams$seed)
        refIndex <- sort(sample(length(mainRanges),newSize))
        mainRanges <- mainRanges[refIndex]
        for (i in 1:length(input)) {
            if (!is.null(input[[i]]$ranges)) {
                newSize <- round(fraction*length(input[[i]]$ranges))
                set.seed(preprocessParams$seed)
                fracIndex <- sort(sample(length(input[[i]]$ranges),newSize))
                input[[i]]$ranges <- input[[i]]$ranges[fracIndex]
            }
            if (!is.null(input[[i]]$coverage))
                input[[i]]$coverage <- input[[i]]$coverage[refIndex]
            if (!is.null(input[[i]]$profile))
                input[[i]]$profile <- input[[i]]$profile[refIndex,]
        }
    }

    # Remove unwanted seqnames from reference ranges
    chrs <- unique(unlist(lapply(input,function(x) {
        if (x$format %in% c("bam","bed"))
            return(as.character(runValue(seqnames(x$ranges))))
        else if (x$format=="bigwig") {
            if (!requireNamespace("rtracklayer"))
                stop("R package rtracklayer is required to read and import ",
                    "BigWig files!")
            return(as.character(seqnames(seqinfo(BigWigFile(x$file)))))
        }
    })))
    #if (type=="chipseq") {
    #    keep <- which(as.character(seqnames(genomeRanges)) %in% chrs)
    #    genomeRanges <- genomeRanges[keep]
    #}
    #else if (type=="rnaseq") {
    #    keeph <- which(as.character(seqnames(helperRanges)) %in% chrs)
    #    helperRanges <- helperRanges[keeph]
    #    genomeRanges <- genomeRanges[names(helperRanges)]
    #    ########################################################################
    #    ## There must be an R bug with `lengths` here as although it runs in 
    #    ## Rcmd, it does not pass package building or vignette kniting... But 
    #    ## for the time being it seems that it is not needed as the name 
    #    ## filtering works
    #    #lens <- which(lengths(genomeRanges)==0)
    #    #if (length(lens)>0)
    #    #    genomeRanges[lens] <- NULL
    #    ########################################################################
    #}
    if (type=="chipseq") {
        keep <- which(as.character(seqnames(mainRanges)) %in% chrs)
        mainRanges <- mainRanges[keep]
    }
    else if (type=="rnaseq") {
        keeph <- which(as.character(seqnames(helperRanges)) %in% chrs)
        helperRanges <- helperRanges[keeph]
        mainRanges <- mainRanges[names(helperRanges)]
        ########################################################################
        ## There must be an R bug with `lengths` here as although it runs in 
        ## Rcmd, it does not pass package building or vignette kniting... But 
        ## for the time being it seems that it is not needed as the name 
        ## filtering works
        #lens <- which(lengths(genomeRanges)==0)
        #if (length(lens)>0)
        #    genomeRanges[lens] <- NULL
        ########################################################################
    }
    
    #if (type=="chipseq")
    #    input <- coverageRef(input,genomeRanges,region,flank,strandedParams,
    #        rc=rc)#,bamParams)
    #else if (type=="rnaseq")
    #    input <- coverageRnaRef(input,genomeRanges,helperRanges,flank,
    #        strandedParams,rc=rc)#,bamParams)
    if (type=="chipseq")
        input <- coverageRef(input,mainRanges,strandedParams,rc=rc)
    else if (type=="rnaseq") #{
        #if (flank[1]==0 && flank[2]==0)
        #    mainRnaRanges <- genomeRanges
        #else {
        #    mainRnaRanges <- tryCatch(
        #        loadMainRnaRanges(genome,refdb,flank),
        #        error=function(e) { # Insanely slow fallback switch
        #            warning("The requested flanking regions are not ",
        #                "compatible with any of the precalculated ones. They ",
        #                "will be calculated on the fly which is a lenghty ",
        #                "process. If you wish to use precalculated flanking ",
        #                "regions, stop the execution and adjust the flanking ",
        #                "parameter.",immediate.=TRUE)
        #            getMainRnaRangesOnTheFly(helperRanges,flank,rc=rc)
        #        },finally=""
        #    )
        #}
        input <- coverageRnaRef(input,mainRanges,strandedParams,rc=rc)
    #}
    # If normalization method is linear, we must adjust the coverages
    if (preprocessParams$normalize=="linear") {
        linFac <- calcLinearFactors(input,preprocessParams)
        for (n in names(input)) {
            if (linFac[n]==1)
                next
            #if (is.null(input[[n]]$coverage$center))
            input[[n]]$coverage <- cmclapply(input[[n]]$coverage,
                function(x,f) {
                    return(x*f)
                },linFac[n]
            )
            #else {
            #    input[[n]]$coverage$center <- 
            #        cmclapply(input[[n]]$coverage$center,function(x,f) {
            #            return(x*f)
            #        },linFac[n])
            #}
        }
    }
    
    # Now we must summarize and create the matrices. If genebody or unequal 
    # custom lengths, bin is a must, else we leave to user
    mustBin <- FALSE
    if (region=="genebody")
        mustBin <- TRUE
    if (region=="custom") {
        #w <- width(genomeRanges)
        w <- width(mainRanges)
        if (any(w!=w[1]))
            mustBin <- TRUE
    }
    if (mustBin) {
        if (binParams$regionBinSize==0) {
            warning("Central region bin size not set for a region that must ",
                "be binned! Setting to 1000...",immediate.=TRUE)
            binParams$regionBinSize <- 1000
        }
    }
    
    # Chunking?
    if (binParams$chunking) {
        if (type=="chipseq")
            binParams$chunks <- split(1:length(genomeRanges),
                as.character(seqnames(genomeRanges)))
        else if (type=="rnaseq")
            binParams$chunks <- split(1:length(helperRanges),
                as.character(seqnames(helperRanges)))
    }
    input <- profileMatrix(input,flank,binParams,rc)
    
    # Perform the k-means clustering if requested and append to design (which
    # has been checked, if we are allowed to do so)
    if (kmParams$k>0)
        design <- kmeansDesign(input,design,kmParams)
    
    # Coverages and profiles calculated... Now depending on plot option, we go 
    # further or return the enriched input object for saving
    if (!plotParams$profile && !plotParams$heatmap) {
        recoupObj <- toOutput(input,design,saveParams,callParams=callParams)
        return(recoupObj)
    }
    else {
        recoupObj <- toOutput(input,design,list(ranges=TRUE,coverage=TRUE,
            profile=TRUE),callParams=callParams)
    }
            
    ## Our plot objects
    recoupPlots <- list()
    
    # We must pass the matrices to plotting function
    if (plotParams$profile) {
        message("Constructing genomic coverage profile curve(s)")
        #theProfile <- recoupProfile(recoupObj,rc=rc)
        recoupObj <- recoupProfile(recoupObj,rc=rc)
        theProfile <- getr(recoupObj,"profile")
        recoupPlots$profilePlot <- theProfile
    }
    
    # Some default binning MUST be applied for the heatmap... Otherwise it could
    # take insanely long time and space to draw/store
    if (plotParams$heatmap) {
        # Inform the user about enforced binning (or not)
        if (region %in% c("tss","tes") || callParams$customIsBase) {
            if (binParams$regionBinSize==0 && binParams$forceHeatmapBinning) {
                message("The resolution of the requested profiles will be ",
                    "lowered to avoid\nincreased computation time and/or ",
                    "storage space for heatmap profiles...")
                
            }
            else if (binParams$regionBinSize==0
                && !binParams$forceHeatmapBinning)
                warning("forceHeatmapBinning is turned off for high ",
                    "resolution plotting. Be prepared for\nlong computational ",
                    "times and big figures!",immediate.=TRUE)
        }
        else {
            if ((binParams$regionBinSize==0 || binParams$flankBinSize==0)
                && binParams$forceHeatmapBinning) {
                message("The resolution of the requested profiles will be ",
                    "lowered to avoid\nincreased computation time and/or ",
                    "storage space for heatmap profiles...")
                
            }
            else if ((binParams$regionBinSize==0 || binParams$flankBinSize==0)
                && !binParams$forceHeatmapBinning)
                warning("forceHeatmapBinning is turned off for high ",
                    "resolution plotting. Be prepared for\nlong computational ",
                    "times and big figures!",immediate.=TRUE)
        }
        
        if (binParams$forceHeatmapBinning 
            && (binParams$regionBinSize==0 || binParams$flankBinSize==0)) {
            helpInput <- recoupObj$data
            if (region %in% c("tss","tes") || callParams$customIsBase) {
                for (n in names(helpInput)) {
                    message("Calculating ",region," profile for ",
                        helpInput[[n]]$name)
                    helpInput[[n]]$profile <- 
                        binCoverageMatrix(helpInput[[n]]$coverage,
                            binSize=binParams$forcedBinSize[2],
                            stat=binParams$sumStat,rc=rc)
                    helpInput[[n]]$profile <- helpInput[[n]]$profile
                }
            }
            else {
                for (n in names(helpInput)) {
                    message("Calculating ",region," profile for ",
                        helpInput[[n]]$name)
                    message(" center")
                    #center <- binCoverageMatrix(helpInput[[n]]$coverage$center,
                    #    binSize=binParams$forcedBinSize[2],
                    #    stat=binParams$sumStat,rc=rc)
                    center <- binCoverageMatrix(
                        input[[n]]$coverage,binSize=binParams$forcedBinSize[2],
                        stat=binParams$sumStat,
                        interpolation=binParams$interpolation,
                        flank=flank,where="center",rc=rc
                    )
                    message(" upstream")
                    #left <- binCoverageMatrix(helpInput[[n]]$coverage$upstream,
                    #    binSize=binParams$forcedBinSize[1],
                    #    stat=binParams$sumStat,rc=rc)
                    left <- binCoverageMatrix(
                        input[[n]]$coverage,binSize=binParams$forcedBinSize[1],
                        stat=binParams$sumStat,
                        interpolation=binParams$interpolation,flank=flank,
                        where="upstream",rc=rc
                    )
                    message(" downstream")
                    #right <- binCoverageMatrix(
                    #    helpInput[[n]]$coverage$downstream,
                    #    binSize=binParams$forcedBinSize[1],
                    #    stat=binParams$sumStat,rc=rc)
                    right <- binCoverageMatrix(
                        input[[n]]$coverage,binSize=binParams$forcedBinSize[1],
                        stat=binParams$sumStat,
                        interpolation=binParams$interpolation,flank=flank,
                        where="downstream",rc=rc
                    )
                    helpInput[[n]]$profile <- cbind(left,center,right)
                    rownames(helpInput[[n]]$profile) <-
                        names(input[[n]]$coverage)
                    helpInput[[n]]$profile <- helpInput[[n]]$profile
                }
            }
        }
        else
            helpInput <- recoupObj$data
        
        helpObj <- recoupObj
        helpObj$data <- helpInput
        message("Constructing genomic coverage heatmap(s)")
        #theHeatmap <- recoupHeatmap(helpObj,rc=rc)
        helpObj <- recoupHeatmap(helpObj,rc=rc)
        theHeatmap <- getr(helpObj,"heatmap")
        recoupObj <- setr(recoupObj,"heatmap",theHeatmap)
        recoupPlots$heatmapPlot <- theHeatmap
        
        # Derive the main heatmap in case of hierarchical clustering        
        mainh <- 1
        if (length(grep("hc",orderBy$what))>0) {
            nc <- nchar(orderBy$what)
            mh <- suppressWarnings(as.numeric(substr(orderBy$what,nc,nc)))
            if (is.na(mh))
                warning("Reference profile for hierarchical clustering order ",
                    "not recognized! Using the 1st...",immediate.=TRUE)
            else if (mh > length(input)) {
                warning("Reference profile (",mh,") for hierarchical ",
                    "clustering order does not exist (the input has only ",
                    length(input)," sources! Using the last...",
                    immediate.=TRUE)
                    mainh <- length(input)
            }
            else
                mainh <- mh
        }
    }
    
    # We must pass the matrices to plotting function
    if (plotParams$correlation) {
        message("Constructing coverage correlation profile curve(s)")
        recoupObj <- recoupCorrelation(recoupObj,rc=rc)
        theCorrelation <- getr(recoupObj,"correlation")
        recoupPlots$correlationPlot <- theCorrelation
    }
    
    # Overwrite final object so as to return it
    recoupObj <- toOutput(input,design,saveParams,recoupPlots,callParams)
    
    # Make any plots asked
    if (plotParams$plot) {
        what <- character(0)
        if (plotParams$profile)
            what <- c(what,"profile")
        if (plotParams$heatmap)
            what <- c(what,"heatmap")
        if (plotParams$correlation)
            what <- c(what,"correlation")
        if (length(what)>0)
            recoupPlot(recoupObj,what,plotParams$device,plotParams$outputDir,
                plotParams$outputBase,mainh)
    }
    
    # Return the enriched input object according to save options
    return(recoupObj)
}

applySelectors <- function(ranges,selector,rc=NULL) {
    if (!is.null(selector$id)) {
        ranges <- ranges[selector$ids]
        if (length(ranges)==0)
            stop("No ranges left after using the identifiers provided with ",
                "the selector field. Are you sure the identifiers between the ",
                "two files are compatible?")
    }
    if (!is.null(selector$biotype) && !is.null(ranges$biotype)) {
        good <- which(ranges$biotype %in% selector$biotype)
        ranges <- ranges[good]
        if (length(ranges)==0)
            stop("No ranges left after using the biotypes provided with the ",
                "selector field. Are you sure the identifiers between the two ",
                "files are compatible?")
    }
    if (!is.null(selector$exonType) && !is.null(ranges$exon_type)) {
        good <- which(ranges$exonType %in% selector$exonType)
        ranges <- ranges[good]
        if (length(ranges)==0)
            stop("No ranges left after using the exon types provided with ",
                "the selector field. Are you sure the identifiers between the ",
                "two files are compatible?")
    }
    return(ranges)
}

toOutput <- function(input,design=NULL,saveParams,plotObjs=NULL,
    callParams=NULL) {
    if (!saveParams$ranges)
        input <- removeData(input,"ranges")
    if (!saveParams$coverage)
        input <- removeData(input,"coverage")
    if (!saveParams$profile)
        input <- removeData(input,"profile")
    plots <- list(profile=NULL,heatmap=NULL,correlation=NULL)
    if (!is.null(plotObjs) && saveParams$profilePlot)
        plots$profile <- plotObjs$profilePlot
    if (!is.null(plotObjs) && saveParams$heatmapPlot)
        plots$heatmap <- plotObjs$heatmapPlot
    if (!is.null(plotObjs) && saveParams$correlationPlot)
        plots$correlation <- plotObjs$correlationPlot
    return(list(
        data=input,
        design=design,
        plots=plots,
        callopts=callParams
    ))
}
