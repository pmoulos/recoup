recoupPlot <- function(recoupObj,what=c("profile","heatmap","correlation"),
    device=c("x11","png","jpg","tiff","bmp","pdf","ps"),outputDir=".",
    outputBase=paste(sapply(recoupObj,function(x) return(x$data$id)),sep="_"),
    mainh=1,...) {
    what <- tolower(what)
    device <- tolower(device[1])
    checkTextArgs("what",what,c("profile","heatmap","correlation"),
        multiarg=TRUE)
    checkTextArgs("device",device,c("x11","png","jpg","tiff","bmp","pdf","ps"))
    
    if ("profile" %in% what) {
        if (!is.null(recoupObj$plots$profile)) {
            theProfile <- recoupObj$plots$profile
            if (device=="x11") {
                dev.new()
                plot(theProfile)
            }
            else
                ggsave(filename=paste(outputBase,"_profile.",device,sep=""),
                    plot=theProfile,path=outputDir,device=device,...)
        }
        else
            message("Profile to plot not found! Are you sure you called ",
                "recoupProfile first?")
    }
    if ("heatmap" %in% what) {
        if (!is.null(recoupObj$plots$heatmap)) {
            theHeatmap <- recoupObj$plots$heatmap
            if (device=="x11") {
                dev.new()
                draw(theHeatmap,gap=grid::unit(1,"cm"),main_heatmap=mainh)
            }
            else {
                graphicsOpen(device,paste(outputBase,"_heatmap.",device,
                    sep=""),...)
                #if (device == "pdf") {
                #   # Starting from width=4, we add 2 inches for each heatmap
                #   iw <- 4 + (length(recoupObj$data)-1)*2
                #    graphicsOpen(device,paste(outputBase,"_heatmap.",
                #        device,sep=""),width=iw,...)
                #}
                #else {
                #   # Starting from width=400, we add 200 pixels for each 
                #   # heatmap
                #   iw <- 400 + (length(recoupObj$data)-1)*200
                #   graphicsOpen(device,paste(outputBase,"_heatmap.",
                #       device,sep=""),width=iw,...)
                #    graphicsOpen(device,paste(outputBase,"_heatmap.",device,
                #       sep=""),...)
                draw(theHeatmap,gap=grid::unit(1,"cm"),main_heatmap=mainh)
                graphicsClose(device)
            }
        }
        else
            message("Heatmap to plot not found! Are you sure you called ",
                "recoupHeatmap first?")
    }
    if ("correlation" %in% what) {
        if (!is.null(recoupObj$plots$correlation)) {
            theCorrelation <- recoupObj$plots$correlation
            if (device=="x11") {
                dev.new()
                plot(theCorrelation)
            }
            else
                ggsave(filename=paste(outputBase,"_correlation.",device,sep=""),
                    plot=theCorrelation,path=outputDir,device=device,...)
        }
        else
            message("Correlation to plot not found! Are you sure you called ",
                "recoupCorrelation first?")
    }
}

recoupProfile <- function(recoupObj,samples=NULL,rc=NULL) {
    # Retrieve data
    input <- recoupObj$data
    design <- recoupObj$design
    #opts <- recoupObj$plotopts
    
    # Attach some config options for profile and heatmap. irrespectively of 
    # subsequent plotting
    ggplotParams <- recoupObj$callopts$ggplotParams
    lineSize <- ifelse(
        recoupObj$callopts$plotParams$device %in% c("x11","png","jpg","bmp"),
        0.7,0.6
    )
    ggplotParams$lineSize <- lineSize
    ggplotParams$singleFacet <- recoupObj$callopts$plotParams$singleFacet
    ggplotParams$multiFacet <- recoupObj$callopts$plotParams$multiFacet
    opts <- list(
        xAxisParams=list(
            region=recoupObj$callopts$region,
            flank=recoupObj$callopts$flank,
            customIsBase=recoupObj$callopts$customIsBase
        ),
        yAxisParams=list(
            signalScale=recoupObj$callopts$plotParams$signalScale,
            heatmapScale=recoupObj$callopts$plotParams$heatmapScale,
            heatmapFactor=recoupObj$callopts$plotParams$heatmapFactor,
            conf=recoupObj$callopts$plotParams$conf
        ),
        binParams=recoupObj$callopts$binParams,
        ggplotParams=ggplotParams
    )
    
    if (is.null(input[[1]]$profile))
        stop("Profile matrix not found in the input object! Are you sure you ",
            "saved it while running the main recoup function? This may occur ",
            "when using recoupProfile and/or recoupHeatmap on an object ",
            "returned by recoup having changed the default saveParams ",
            "parameter")
    
    if (!requireNamespace("grid"))
        stop("R package grid is required to create average profile plots!")
    
    # Create the x-axis breaks and labels
    width <- ncol(input[[1]]$profile)
    lb <- makeHorizontalAnnotation(width,opts,"profile")
    breaks <- lb$breaks
    labels <- lb$labels
    
    # Filter the profiles here
    if (!is.null(samples)) {
        if (is.numeric(samples))
            ta <- 1:length(input)
        else 
            ta <- names(input)
        if (!all(samples %in% ta) || length(samples)==0)
            warning("Wrong indexing for profile plot subsetting! ",
                "Ignoring",immediate.=TRUE)
        else
            input <- input[samples]
    }
    
    profileColors <- unlist(sapply(input,function(x) return(x$color)))
    if (!is.null(profileColors))
        names(profileColors) <- unlist(sapply(input,function(x) return(x$name)))
    
    ggParams <- opts$ggplotParams
    
    if (is.null(design)) {
        # Create ggplot data
        profiles <- calcPlotProfiles(input,opts,2,rc)
        index <- 1:length(profiles[[1]][[1]])
        names <- sapply(input,function(x) return(x$name))
        position <- rep(index,length(input))
        signal <- unlist(lapply(profiles,function(x) return(x$profile)))
        ymin <- unlist(lapply(profiles,function(x) return(x$lower)))
        ymax <- unlist(lapply(profiles,function(x) return(x$upper)))
        condition <- rep(names,each=length(index))
        
        # Create ggplot data frame
        ggplot.data <- data.frame(
            Position=position,
            Signal=signal,
            Condition=factor(condition,levels=unique(condition)),
            ymin=ymin,
            ymax=ymax
        )
        
        # Create ggplot plot
        ggplot.plot <-
            ggplot(ggplot.data,mapping=aes(x=Position,y=Signal,
                colour=Condition)) + 
            geom_line(size=ggParams$lineSize)
        
        if (opts$yAxisParams$conf)
            ggplot.plot <- ggplot.plot +
                geom_ribbon(aes(x=Position,ymin=ymin,ymax=ymax,colour=Condition,
                    fill=Condition),alpha=0.3,size=0) 
        
        ggplot.plot <- ggplot.plot +
            theme_bw() +
            xlab("\nPosition in bp") +
            ylab("Average signal\n") +
            theme(title=ggParams$title,
                axis.title.x=ggParams$axis.title.x,
                axis.title.y=ggParams$axis.title.y,
                axis.text.x=ggParams$axis.text.x,
                axis.text.y=ggParams$axis.text.y,
                legend.position=ggParams$legend.position) +
             scale_x_continuous(breaks=breaks,labels=labels)
                
        if (!is.null(profileColors))
            ggplot.plot <- ggplot.plot + 
                scale_fill_manual(values=profileColors) +
                scale_color_manual(values=profileColors)
    }
    else {
        message("Using provided design to facet the coverage profiles")
        subcovmat <- lapply(input,function(x,d) {
            splitter <- split(rownames(x$profile),as.list(d),sep="|",drop=TRUE)
            out <- vector("list",length(splitter))
            names(out) <- names(splitter)
            for (n in names(splitter))
                out[[n]] <- x$profile[splitter[[n]],,drop=FALSE]
            return(out)
        },design)
        subProfiles <- lapply(subcovmat,function(x,opts,rc) {
            return(calcDesignPlotProfiles(x,opts,2,rc))
        },opts,rc)
        
        designProfiles <- lapply(subProfiles,function(x) {
            d <- names(x) # Should be the "|" separated factor names
            dsplit <- strsplit(d,split="|",fixed=TRUE)
            # Replication factors
            pop <- sapply(x,function(x) return(length(x$profile)))
            tmp <- vector("list",length(x))
            for (i in 1:length(dsplit)) {
                if (length(dsplit[[i]])>1)
                    tmp[[i]] <- t(replicate(pop[i],dsplit[[i]]))
                else
                    tmp[[i]] <- as.matrix(rep(dsplit[[i]],pop[i]))
            }
            o <- list()
            o$profile <- unlist(lapply(x,function(y) {
                return(y$profile)
            }),use.names=FALSE)
            o$upper <- unlist(lapply(x,function(y) {
                return(y$upper)
            }),use.names=FALSE)
            o$lower <- unlist(lapply(x,function(y) {
                return(y$lower)
            }),use.names=FALSE)
            o$design <- do.call("rbind",tmp)
            return(o)
        })
        
        index <- 1:ncol(input[[1]]$profile)
        names <- sapply(input,function(x) return(x$name))
        faceter <- do.call("rbind",lapply(designProfiles,function(x) 
                return(x$design)))
        m <- length(which(!duplicated(faceter)))
        position <- rep(index,length(input)*m)
        condition <- rep(names,each=length(index)*m)
        signal <- unlist(lapply(designProfiles,function(x)
            return(x$profile)),use.names=FALSE)
        ymin <- unlist(lapply(designProfiles,function(x) 
            return(x$lower)),use.names=FALSE)
        ymax <- unlist(lapply(designProfiles,function(x) 
            return(x$upper)),use.names=FALSE)
        
        if (length(input)>1) { # Case where two factors max
            ggplot.data <- data.frame(
                Position=position,
                Signal=signal,
                Condition=factor(condition,levels=unique(condition)),
                ymin=ymin,
                ymax=ymax
            )
            
            if (ncol(design)==1)
                ggplot.data$fac1 <- factor(as.character(faceter[,1]),
                    levels=unique(as.character(faceter[,1])))
            if (ncol(design)==2) {
                ggplot.data$fac1 <- factor(as.character(faceter[,1]),
                    levels=unique(as.character(faceter[,1])))
                ggplot.data$fac2 <- factor(as.character(faceter[,2]),
                    levels=unique(as.character(faceter[,2])))
            }
            
            ggplot.plot <-
                ggplot(ggplot.data,mapping=aes(x=Position,y=Signal,
                    colour=Condition)) + 
                geom_line(size=ggParams$lineSize) 
                
            if (opts$yAxisParams$conf)
                ggplot.plot <- ggplot.plot +
                    geom_ribbon(aes(x=Position,ymin=ymin,ymax=ymax,
                        colour=Condition,fill=Condition),alpha=0.3,size=0) 
                        
            ggplot.plot <- ggplot.plot +
                theme_bw() +
                xlab("\nPosition in bp") +
                ylab("Average signal\n") +
                theme(title=ggParams$title,
                    axis.title.x=ggParams$axis.title.x,
                    axis.title.y=ggParams$axis.title.y,
                    axis.text.x=ggParams$axis.text.x,
                    axis.text.y=ggParams$axis.text.y,
                    strip.text.x=ggParams$strip.text.x,
                    strip.text.y=ggParams$strip.text.y,
                    legend.position=ggParams$legend.position,
                    panel.margin=ggParams$panel.margin) +
                scale_x_continuous(breaks=breaks,labels=labels)

            if (!is.null(profileColors))
                ggplot.plot <- ggplot.plot + 
                    scale_fill_manual(values=profileColors) +
                    scale_color_manual(values=profileColors)
            
            if (ncol(design)==1) {
                if (ggParams$multiFacet=="wrap")
                    ggplot.plot <- ggplot.plot + facet_wrap(~ fac1)
                else if (ggParams$multiFacet=="grid")
                    ggplot.plot <- ggplot.plot + facet_grid(fac1~.)
            }
            if (ncol(design)==2)
                ggplot.plot <- ggplot.plot + facet_grid(fac1~fac2)
        }
        else {
            ggplot.data <- data.frame(
                Position=position,
                Signal=signal,
                Condition=condition,
                ymin=ymin,
                ymax=ymax
            )
            
            if (ncol(design)==1) {
                ggplot.data$Design <- factor(as.character(faceter[,1]),
                    levels=unique(as.character(faceter[,1])))
            }
            if (ncol(design)==2) {
                ggplot.data$Design <- factor(as.character(faceter[,1]),
                    levels=unique(as.character(faceter[,1])))
                ggplot.data$fac2 <- factor(as.character(faceter[,2]),
                    levels=unique(as.character(faceter[,2])))
            }
            if (ncol(design)==3) {
                ggplot.data$Design <- factor(as.character(faceter[,1]),
                    levels=unique(as.character(faceter[,1])))
                ggplot.data$fac2 <- factor(as.character(faceter[,2]),
                    levels=unique(as.character(faceter[,2])))
                ggplot.data$fac3 <- factor(as.character(faceter[,3]),
                    levels=unique(as.character(faceter[,3])))
            }
            
            if (opts$yAxisParams$conf) {
                if (ggParams$singleFacet=="none")
                    ggplot.plot <- ggplot(ggplot.data,mapping=aes(x=Position,
                        y=Signal,colour=Design)) +
                        geom_line(size=ggParams$lineSize) +
                        geom_ribbon(aes(x=Position,ymin=ymin,ymax=ymax,
                            colour=Design,fill=Design),alpha=0.3,size=0)
                else
                    ggplot.plot <- ggplot(ggplot.data,mapping=aes(x=Position,
                        y=Signal,colour=Condition)) +
                        geom_line(size=ggParams$lineSize) +
                        geom_ribbon(aes(x=Position,ymin=ymin,ymax=ymax,
                            colour=Condition,fill=Condition),alpha=0.3,size=0)
            }
            else {
                if (ggParams$singleFacet=="none")
                    ggplot.plot <- ggplot(ggplot.data,mapping=aes(x=Position,
                        y=Signal,colour=Design)) +
                        geom_line(size=ggParams$lineSize)
                else
                    ggplot.plot <- ggplot(ggplot.data,mapping=aes(x=Position,
                        y=Signal,colour=Condition)) +
                        geom_line(size=ggParams$lineSize)
            }
            ggplot.plot <- ggplot.plot +
                theme_bw() +
                xlab("\nPosition in bp") +
                ylab("Average signal\n") +
                theme(title=ggParams$title,
                    axis.title.x=ggParams$axis.title.x,
                    axis.title.y=ggParams$axis.title.y,
                    axis.text.x=ggParams$axis.text.x,
                    axis.text.y=ggParams$axis.text.y,
                    strip.text.x=ggParams$strip.text.x,
                    strip.text.y=ggParams$strip.text.y,
                    legend.position=ggParams$legend.position,
                    panel.margin=ggParams$panel.margin) +
                scale_x_continuous(breaks=breaks,labels=labels)
            
             if (!is.null(profileColors) && ggParams$singleFacet!="none") 
                ggplot.plot <- ggplot.plot + 
                    scale_fill_manual(values=profileColors) +
                    scale_color_manual(values=profileColors)
            
            if (ncol(design)==1) {
                if (ggParams$singleFacet=="wrap")
                    ggplot.plot <- ggplot.plot + facet_wrap(~ Design)
                else if (ggParams$singleFacet=="grid")
                    ggplot.plot <- ggplot.plot + facet_grid(Design~.)
            }
            if (ncol(design)==2) {
                if (ggParams$multiFacet=="wrap")
                    ggplot.plot <- ggplot.plot + facet_wrap(~ fac2)
                else if (ggParams$multiFacet=="grid")
                    ggplot.plot <- ggplot.plot + facet_grid(fac2~.)
            }
            if (ncol(design)==3)
                ggplot.plot <- ggplot.plot + facet_grid(fac2~fac3)
        }
   }
   #return(ggplot.plot)
   recoupObj <- setr(recoupObj,"profile",ggplot.plot)
}

recoupHeatmap <- function(recoupObj,samples=NULL,rc=NULL) {
    input <- recoupObj$data
    design <- recoupObj$design
    #opts <- recoupObj$plotopts
    
    # Attach some config options for profile and heatmap, irrespectively of 
    # subsequent plotting
    opts <- list(
        xAxisParams=list(
            region=recoupObj$callopts$region,
            flank=recoupObj$callopts$flank,
            customIsBase=recoupObj$callopts$customIsBase
        ),
        yAxisParams=list(
            signalScale=recoupObj$callopts$plotParams$signalScale,
            heatmapScale=recoupObj$callopts$plotParams$heatmapScale,
            heatmapFactor=recoupObj$callopts$plotParams$heatmapFactor
        ),
        binParams=recoupObj$callopts$binParams,
        orderBy=recoupObj$callopts$orderBy,
        complexHeatmapParams=recoupObj$callopts$complexHeatmapParams
    )
    
    if (is.null(input[[1]]$profile))
        stop("Profile matrix not found in the input object! Are you sure you ",
            "saved it while running the main recoup function? This may occur ",
            "when using recoupProfile and/or recoupHeatmap on an object ",
            "returned by recoup having changed the default saveParams ",
            "parameter")
    
    # Check compatibility of orderBy argument and ComplexHeatmap parameters
    # Hierarchical clustering asked in orderBy but otherwise in the heatmap
    # parameters. Clustering is performed.
    if (length(grep("hc",opts$orderBy$what))>0
        && !(opts$complexHeatmapParams$main$cluster_rows 
        || opts$complexHeatmapParams$group$cluster_rows)) {
        warning("Hierarchical clustering asked in the orderBy parameter but ",
            "is set to\nFALSE in complexHeatmapParams! Will auto-correct to ",
            "perform hierarchical\nclustering ordering...",immediate.=TRUE)
        opts$complexHeatmapParams$main$cluster_rows <- TRUE
        opts$complexHeatmapParams$group$cluster_rows <- TRUE
    }
    # Hierarchical clustering asked in heatmap parameters but not in the 
    # orderBy directives. Clustering is not performed.
    if ((opts$complexHeatmapParams$main$cluster_rows 
        || opts$complexHeatmapParams$group$cluster_rows)
        && length(grep("hc",opts$orderBy$what))==0) {
        warning("Hierarchical clustering asked in the complexHeatmapParams ",
            "parameter but\nnot in orderBy parameter! Hierarchical clustering ",
            "in the heatmap\nprofile will be turned off",immediate.=TRUE)
        opts$complexHeatmapParams$main$cluster_rows <- FALSE
        opts$complexHeatmapParams$group$cluster_rows <- FALSE
    }
    
    width <- ncol(input[[1]]$profile)
    labCol <- rep("",width)
    lb <- makeHorizontalAnnotation(width,opts,"heatmap")
    labCol[round(lb$breaks)] <- lb$labels
    
    # Filter the profiles here
    if (!is.null(samples)) {
        if (is.numeric(samples))
            ta <- 1:length(input)
        else 
            ta <- names(input)
        if (!all(samples %in% ta) || length(samples)==0)
            warning("Wrong indexing for profile plot subsetting! ",
                "Ignoring",immediate.=TRUE)
        else
            input <- input[samples]
    }
    
    if (opts$yAxisParams$signalScale=="log2") {
        for (n in names(input)) {
            input[[n]]$profile <- input[[n]]$profile + 1
            input[[n]]$profile <- log2(input[[n]]$profile)
        }
    }
    
    haCol <- HeatmapAnnotation(cn=function(index) {
        width <- ncol(input[[1]]$profile)
        labCol <- rep("",width)
        lb <- makeHorizontalAnnotation(width,opts,"heatmap")
        labCol[round(lb$breaks)] <- lb$labels
        grid.text(labCol,(1:width)/width,1,vjust=1,
            gp=grid::gpar(cex=0.7))
    })
    
    # Here we need to apply the orderBy directives. If orderBy$what is 
    # something else that hierarchical clustering, parameter control in main
    # has already fixed incompatibilities by turning off clustering. If it
    # is clustering, it's turned on, so we need to decide which is the main
    # profile.
    # The same ordering ideas must be applied to the design elements... We must
    # reorder the rownames of the design matrix in such a way as as to be ready 
    # when fed to Heatmap function.
    if (is.null(design))
        sorter <- orderProfiles(input,opts,rc=rc)
    else
        sorter <- orderProfilesByDesign(input,design,opts,rc=rc)
    
    colorFuns <- vector("list",length(input))
    names(colorFuns) <- names(input)
    profileColors <- unlist(sapply(input,function(x) return(x$color)))
    if (opts$yAxisParams$heatmapScale=="each") {
        for (n in names(colorFuns)) {
            qs <- c(0.95,0.96,0.97,0.98,0.99,0.995,0.999)
            pos <- 1
            sup <- quantile(input[[n]]$profile,qs[pos])
            while(sup==0) {
                pos <- pos + 1
                sup <- quantile(input[[n]]$profile,qs[pos])
                if (sup!=0)
                    break
            }
            supp <- opts$yAxisParams$heatmapFactor * sup
            if (!is.null(profileColors))
                colorFuns[[n]] <- colorRamp2(c(0,supp),c("white",
                    input[[n]]$color))
            else
                colorFuns[[n]] <- colorRamp2(c(0,supp),c("white","red2"))
        }
    }
    else if (opts$yAxisParams$heatmapScale=="common") {
        sups <- unlist(sapply(input,function(x) {
            return(quantile(x$profile,0.95))
        }))
        sup <- max(sups)
        supp <- opts$yAxisParams$heatmapFactor * sup
        for (n in names(colorFuns)) {
            if (!is.null(profileColors))
                colorFuns[[n]] <- colorRamp2(c(0,supp), c("white",
                    input[[n]]$color))
            else
                colorFuns[[n]] <- colorRamp2(c(0,supp), c("white","red2"))
        }
    }
    
    if (is.null(design)) {
        hParams <- opts$complexHeatmapParams$main
        hmList <- NULL
        for (n in names(input)) {
            input[[n]]$profile <- as.matrix(input[[n]]$profile)
            colnames(input[[n]]$profile) <- labCol
            hmList <- hmList + 
                Heatmap(
                    input[[n]]$profile[sorter$ix,],
                    name=input[[n]]$name,
                    cluster_rows=hParams$cluster_rows,
                    cluster_columns=hParams$cluster_columns,
                    column_title_gp=hParams$column_title_gp,
                    show_row_names=hParams$show_row_names,
                    show_column_names=hParams$show_column_names,
                    col=colorFuns[[n]],
                    column_title=paste(input[[n]]$name,"signal"),
                    heatmap_legend_param=hParams$heatmap_legend_param,
                    bottom_annotation=haCol
                )
        }
    }
    else {
        message("Using provided design to facet the coverage profiles")
        hParams <- opts$complexHeatmapParams$group
        hmList <- NULL
        for (n in names(input)) {
            input[[n]]$profile <- as.matrix(input[[n]]$profile)
            colnames(input[[n]]$profile) <- labCol
            hmList <- hmList + 
                Heatmap(
                    input[[n]]$profile,
                    name=input[[n]]$name,
                    cluster_rows=hParams$cluster_rows,
                    cluster_columns=hParams$cluster_columns,
                    column_title_gp=hParams$column_title_gp,
                    show_row_names=hParams$show_row_names,
                    show_column_names=hParams$show_column_names,
                    col=colorFuns[[n]],
                    column_title=paste(input[[n]]$name,"signal"),
                    heatmap_legend_param=hParams$heatmap_legend_param,
                    bottom_annotation=haCol,
                    split=design,
                    row_order=sorter,
                    row_title_gp=hParams$row_title_gp,
                    gap=hParams$gap
                )
       }
   }
   #return(hmList)
   recoupObj <- setr(recoupObj,"heatmap",hmList)
}

recoupCorrelation <- function(recoupObj,samples=NULL,rc=NULL) {
    # Retrieve data
    input <- recoupObj$data
    design <- recoupObj$design
    
    # Attach some config options for profile and heatmap. irrespectively of 
    # subsequent plotting
    ggplotParams <- recoupObj$callopts$ggplotParams
    lineSize <- ifelse(
        recoupObj$callopts$plotParams$device %in% c("x11","png","jpg","bmp"),
        0.7,0.6
    )
    ggplotParams$lineSize <- lineSize
    ggplotParams$singleFacet <- recoupObj$callopts$plotParams$singleFacet
    ggplotParams$multiFacet <- recoupObj$callopts$plotParams$multiFacet
    opts <- list(
        yAxisParams=list(
            signalScale=recoupObj$callopts$plotParams$signalScale,
            conf=recoupObj$callopts$plotParams$conf,
            corrScale=recoupObj$callopts$plotParams$corrScale,
            corrSmoothPar=recoupObj$callopts$plotParams$corrSmoothPar
        ),
        orderBy=recoupObj$callopts$orderBy,
        binParams=recoupObj$callopts$binParams,
        ggplotParams=ggplotParams
    )
    
    if (is.null(input[[1]]$profile))
        stop("Profile matrix not found in the input object! Are you sure you ",
            "saved it while running the main recoup function? This may occur ",
            "when using recoupProfile and/or recoupHeatmap on an object ",
            "returned by recoup having changed the default saveParams ",
            "parameter")
    
    if (!requireNamespace("grid"))
        stop("R package grid is required to create average profile plots!")
    
    # Filter the profiles here
    if (!is.null(samples)) {
        if (is.numeric(samples))
            ta <- 1:length(input)
        else 
            ta <- names(input)
        if (!all(samples %in% ta) || length(samples)==0)
            warning("Wrong indexing for profile plot subsetting! ",
                "Ignoring",immediate.=TRUE)
        else
            input <- input[samples]
    }
    
    # Normalize the profile if requested
    if (opts$yAxisParams$corrScale=="normalized") {
        for (n in names(input))
            input[[n]]$profile <- input[[n]]$profile/max(input[[n]]$profile)
    }
    
    profileColors <- unlist(sapply(input,function(x) return(x$color)))
    if (!is.null(profileColors))
        names(profileColors) <- unlist(sapply(input,function(x) return(x$name)))
    
    ggParams <- opts$ggplotParams
    
    if (is.null(design)) {
        # Create ggplot data
        profiles <- calcPlotProfiles(input,opts,1,rc)
        sorter <- orderSignals(profiles,opts)
        index <- 1:length(profiles[[1]][[1]])
        names <- sapply(input,function(x) return(x$name))
        position <- rep(index,length(input))
        signal <- unlist(lapply(profiles,function(x,s,p) {
            return(lowess(x$profile[s],f=p)$y)
            #return(smooth.spline(x$profile[s])$y)
        },sorter$ix,opts$yAxisParams$corrSmoothPar))
        ymin <- unlist(lapply(profiles,function(x,s,p) {
            return(lowess(x$lower[s],f=p)$y)
            #return(smooth.spline(x$lower[s])$y)
        },sorter$ix,opts$yAxisParams$corrSmoothPar))
        ymax <- unlist(lapply(profiles,function(x,s,p) {
            return(lowess(x$upper[s],f=p)$y)
            #return(smooth.spline(x$upper[s])$y)
        },sorter$ix,opts$yAxisParams$corrSmoothPar))
        condition <- rep(names,each=length(index))
        
        # Create ggplot data frame
        ggplot.data <- data.frame(
            Index=position,
            Coverage=signal,
            Condition=factor(condition,levels=unique(condition)),
            ymin=ymin,
            ymax=ymax
        )
        
        # Create ggplot plot
        ggplot.plot <-
            ggplot(ggplot.data,mapping=aes(x=Index,y=Coverage,
                colour=Condition)) + 
            geom_line(size=ggParams$lineSize)
       
        if (opts$yAxisParams$conf)
            ggplot.plot <- ggplot.plot +
                geom_ribbon(aes(x=Index,ymin=ymin,ymax=ymax,colour=Condition,
                    fill=Condition),alpha=0.3,size=0)
        
        ggplot.plot <- ggplot.plot +
            theme_bw() +
            xlab("\nIndex") +
            ylab("Average coverage\n") +
            theme(title=ggParams$title,
                axis.title.x=ggParams$axis.title.x,
                axis.title.y=ggParams$axis.title.y,
                axis.text.x=ggParams$axis.text.x,
                axis.text.y=ggParams$axis.text.y,
                legend.position=ggParams$legend.position)
                
        if (!is.null(profileColors))
            ggplot.plot <- ggplot.plot + 
                scale_fill_manual(values=profileColors) +
                scale_color_manual(values=profileColors)
    }
    else {
        message("Using provided design to facet the coverage profiles")
        subcovmat <- lapply(input,function(x,d) {
            splitter <- split(rownames(x$profile),as.list(d),sep="|",drop=TRUE)
            out <- vector("list",length(splitter))
            names(out) <- names(splitter)
            for (n in names(splitter))
                out[[n]] <- x$profile[splitter[[n]],,drop=FALSE]
            return(out)
        },design)
        subProfiles <- lapply(subcovmat,function(x,opts,rc) {
            return(calcDesignPlotProfiles(x,opts,1,rc))
        },opts,rc)
        
        designProfiles <- lapply(subProfiles,function(x) {
            d <- names(x) # Should be the "|" separated factor names
            dsplit <- strsplit(d,split="|",fixed=TRUE)
            # Replication factors
            pop <- sapply(x,function(x) return(length(x$profile)))
            tmp <- vector("list",length(x))
            for (i in 1:length(dsplit)) {
                if (length(dsplit[[i]])>1)
                    tmp[[i]] <- t(replicate(pop[i],dsplit[[i]]))
                else
                    tmp[[i]] <- as.matrix(rep(dsplit[[i]],pop[i]))
            }
            o <- list()
            o$profile <- unlist(lapply(x,function(y) {
                if (any(is.na(y$profile)))
                    y$profile[is.na(y$profile)] <- 0
                return(y$profile)
            }),use.names=FALSE)
            o$upper <- unlist(lapply(x,function(y) {
                if (any(is.na(y$upper)))
                    y$upper[is.na(y$upper)] <- 0
                return(y$upper)
            }),use.names=FALSE)
            o$lower <- unlist(lapply(x,function(y) {
                if (any(is.na(y$lower)))
                    y$lower[is.na(y$lower)] <- 0
                return(y$lower)
            }),use.names=FALSE)
            o$design <- do.call("rbind",tmp)
            return(o)
        })
        
        sorter <- orderDesignSignals(designProfiles,design,opts)
        names <- sapply(input,function(x) return(x$name))
        faceter <- do.call("rbind",lapply(designProfiles,function(x) 
                return(x$design)))
        position <- unlist(lapply(subProfiles,function(x) {
            return(unlist(lapply(x,function(y) return(1:length(y$profile)))))
        }))
        condition <- rep(names,sapply(designProfiles,function(x) {
            return(length(x$profile))
        }))
        signal <- unlist(lapply(designProfiles,function(x,s,p) {
            #return(lowess(x$profile[s],f=0.05)$y),sorter),
            if (any(is.na(x$lower[s])) || length(x$profile[s])<4)
                return(x$profile[s])
            else
                return(smooth.spline(x$profile[s],spar=p)$y)
        },sorter,opts$yAxisParams$corrSmoothPar),use.names=FALSE)
        ymin <- unlist(lapply(designProfiles,function(x,s,p) {
            #return(lowess(x$lower[s],f=0.05)$y),sorter),
            if (any(is.na(x$lower[s])) || length(x$lower[s])<4)
                return(x$lower[s])
            else
                return(smooth.spline(x$lower[s],spar=p)$y)
        },sorter,opts$yAxisParams$corrSmoothPar),use.names=FALSE)
        ymax <- unlist(lapply(designProfiles,function(x,s,p) {
            #return(lowess(x$upper[s],f=0.05)$y),sorter),
            if (any(is.na(x$lower[s])) || length(x$lower[s])<4)
                return(x$lower[s])
            else
                return(smooth.spline(x$upper[s],spar=p)$y)
        },sorter,opts$yAxisParams$corrSmoothPar),use.names=FALSE)
            
        if (length(input)>1) { # Case where two factors max
            ggplot.data <- data.frame(
                Index=position,
                Coverage=signal,
                Condition=factor(condition,levels=unique(condition)),
                ymin=ymin,
                ymax=ymax
            )
            
            if (ncol(design)==1)
                ggplot.data$fac1 <- factor(as.character(faceter[,1]),
                    levels=unique(as.character(faceter[,1])))
            if (ncol(design)==2) {
                ggplot.data$fac1 <- factor(as.character(faceter[,1]),
                    levels=unique(as.character(faceter[,1])))
                ggplot.data$fac2 <- factor(as.character(faceter[,2]),
                    levels=unique(as.character(faceter[,2])))
            }
            
            ggplot.plot <-
                ggplot(ggplot.data,mapping=aes(x=Index,y=Coverage,
                    colour=Condition)) + 
                geom_line(size=ggParams$lineSize)
            
            if (opts$yAxisParams$conf)
                ggplot.plot <- ggplot.plot +
                    geom_ribbon(aes(x=Index,ymin=ymin,ymax=ymax,
                        colour=Condition,fill=Condition),alpha=0.3,size=0)
                        
            ggplot.plot <- ggplot.plot +
                theme_bw() +
                xlab("\nIndex") +
                ylab("Average coverage\n") +
                theme(title=ggParams$title,
                    axis.title.x=ggParams$axis.title.x,
                    axis.title.y=ggParams$axis.title.y,
                    axis.text.x=ggParams$axis.text.x,
                    axis.text.y=ggParams$axis.text.y,
                    strip.text.x=ggParams$strip.text.x,
                    strip.text.y=ggParams$strip.text.y,
                    legend.position=ggParams$legend.position,
                    panel.margin=ggParams$panel.margin)

            if (!is.null(profileColors))
                ggplot.plot <- ggplot.plot + 
                    scale_fill_manual(values=profileColors) +
                    scale_color_manual(values=profileColors)
            
            if (ncol(design)==1) {
                if (ggParams$multiFacet=="wrap")
                    ggplot.plot <- ggplot.plot + facet_wrap(~ fac1,
                        scales="free_x")
                else if (ggParams$multiFacet=="grid")
                    ggplot.plot <- ggplot.plot + facet_grid(fac1~.,
                        scales="free_x")
            }
            if (ncol(design)==2)
                ggplot.plot <- ggplot.plot + facet_grid(fac1~fac2,
                        scales="free_x")
        }
        else {
            ggplot.data <- data.frame(
                Index=position,
                Coverage=signal,
                Condition=condition,
                ymin=ymin,
                ymax=ymax
            )
            
            if (ncol(design)==1) {
                ggplot.data$Design <- factor(as.character(faceter[,1]),
                    levels=unique(as.character(faceter[,1])))
            }
            if (ncol(design)==2) {
                ggplot.data$Design <- factor(as.character(faceter[,1]),
                    levels=unique(as.character(faceter[,1])))
                ggplot.data$fac2 <- factor(as.character(faceter[,2]),
                    levels=unique(as.character(faceter[,2])))
            }
            if (ncol(design)==3) {
                ggplot.data$Design <- factor(as.character(faceter[,1]),
                    levels=unique(as.character(faceter[,1])))
                ggplot.data$fac2 <- factor(as.character(faceter[,2]),
                    levels=unique(as.character(faceter[,2])))
                ggplot.data$fac3 <- factor(as.character(faceter[,3]),
                    levels=unique(as.character(faceter[,3])))
            }
            
            if (ggParams$singleFacet=="none") {
                ggplot.plot <- ggplot(ggplot.data,mapping=aes(x=Index,
                    y=Coverage,colour=Design)) +
                    geom_line(size=ggParams$lineSize)
                if (opts$yAxisParams$conf)
                    ggplot.plot <- ggplot.plot +
                        geom_ribbon(aes(x=Index,ymin=ymin,ymax=ymax,
                            colour=Design,fill=Design),alpha=0.3,size=0)
            }
            else {
                ggplot.plot <- ggplot(ggplot.data,mapping=aes(x=Index,
                    y=Coverage,colour=Condition)) +
                    geom_line(size=ggParams$lineSize)
                if (opts$yAxisParams$conf)
                    ggplot.plot <- ggplot.plot +
                        geom_ribbon(aes(x=Index,ymin=ymin,ymax=ymax),
                            alpha=0.3,size=0)
            }
            ggplot.plot <- ggplot.plot +
                theme_bw() +
                xlab("\nIndex") +
                ylab("Average coverage\n") +
                theme(title=ggParams$title,
                    axis.title.x=ggParams$axis.title.x,
                    axis.title.y=ggParams$axis.title.y,
                    axis.text.x=ggParams$axis.text.x,
                    axis.text.y=ggParams$axis.text.y,
                    strip.text.x=ggParams$strip.text.x,
                    strip.text.y=ggParams$strip.text.y,
                    legend.position=ggParams$legend.position,
                    panel.margin=ggParams$panel.margin)
            
             if (!is.null(profileColors) && ggParams$singleFacet!="none") 
                ggplot.plot <- ggplot.plot + 
                    scale_fill_manual(values=profileColors) +
                    scale_color_manual(values=profileColors)
            
            if (ncol(design)==1) {
                if (ggParams$singleFacet=="wrap")
                    ggplot.plot <- ggplot.plot + facet_wrap(~ Design)
                else if (ggParams$singleFacet=="grid")
                    ggplot.plot <- ggplot.plot + facet_grid(Design~.)
            }
            if (ncol(design)==2) {
                if (ggParams$multiFacet=="wrap")
                    ggplot.plot <- ggplot.plot + facet_wrap(~ fac2)
                else if (ggParams$multiFacet=="grid")
                    ggplot.plot <- ggplot.plot + facet_grid(fac2~.)
            }
            if (ncol(design)==3)
                ggplot.plot <- ggplot.plot + facet_grid(fac2~fac3)
        }
   }
   #return(ggplot.plot)
   recoupObj <- setr(recoupObj,"correlation",ggplot.plot)
}

calcPlotProfiles <- function(input,opts,sdim=c(2,1),rc) {
    sdim <- sdim[1]
    if (opts$binParams$smooth)
        profiles <- cmclapply(input,function(x,avgfun,scale) {
            if (scale=="log2") {
                x$profile <- x$profile + 1
                x$profile <- log2(x$profile)
            }
            o <- list()
            tryCatch({
                fit <- smooth.spline(apply(x$profile,sdim,avgfun))
                ci <- ssCI(fit)
                o$profile <- fit$y
                o$upper <- ci$upper
                o$lower <- ci$lower
                return(o)
            },error=function(e) {
                message("Caught splines error: ",e)
                o$profile <- apply(x$profile,sdim,avgfun)
                varfun <- ifelse(avgfun=="mean","sd","mad")
                va <- apply(x,2,varfun)
                o$upper <- o$profile + va
                o$lower <- o$profile - va
                return(o)
            },finally="")
        },opts$binParams$sumStat,opts$yAxisParams$signalScale,rc=rc)
    else
        profiles <- cmclapply(input,function(x,avgfun,scale) {
            if (scale=="log2") {
                x$profile <- x$profile + 1
                x$profile <- log2(x$profile)
            }
            o <- list()
            o$profile <- apply(x$profile,sdim,avgfun)
            varfun <- ifelse(avgfun=="mean","sd","mad")
            va <- apply(x,2,varfun)
            o$upper <- o$profile + va
            o$lower <- o$profile - va
            return(o)
        },opts$binParams$sumStat,opts$yAxisParams$signalScale,rc=rc)
    return(profiles)
}

calcDesignPlotProfiles <- function(covmat,opts,sdim=c(2,1),rc) {
    sdim <- sdim[1]
    if (opts$binParams$smooth)
        profiles <- cmclapply(covmat,function(x,avgfun,scale) {
            if (scale=="log2") {
                x <- x + 1
                x <- log2(x)
            }
            o <- list()
            tryCatch({
                fit <- smooth.spline(apply(x,sdim,avgfun))
                ci <- ssCI(fit)
                o$profile <- fit$y
                o$upper <- ci$upper
                o$lower <- ci$lower
                return(o)
            },error=function(e) {
                message("Caught splines error: ",e)
                o$profile <- apply(x,sdim,avgfun)
                varfun <- ifelse(avgfun=="mean","sd","mad")
                va <- apply(x,2,varfun)
                o$upper <- o$profile + va
                o$lower <- o$profile - va
                return(o)
            },finally="")
        },opts$binParams$sumStat,opts$yAxisParams$signalScale,rc=rc)
    else
        profiles <- cmclapply(covmat,function(x,avgfun,scale) {
            if (scale=="log2") {
                x <- x + 1
                x <- log2(x)
            }
            o <- list()
            o$profile <- apply(x,sdim,avgfun)
            varfun <- ifelse(avgfun=="mean","sd","mad")
            va <- apply(x,2,varfun)
            o$upper <- o$profile + va
            o$lower <- o$profile - va
            return(o)
        },opts$binParams$sumStat,opts$yAxisParams$signalScale,rc=rc)
    return(profiles)
}

orderProfiles <- function(input,opts,rc=NULL) {
    if (!is.null(opts$orderBy$custom)) {
        if (opts$orderBy$order=="descending")
            sorter <- sort(opts$orderBy$custom,decreasing=TRUE,
                index.return=TRUE)
        else
            sorter <- sort(opts$orderBy$custom,index.return=TRUE)
        return(sorter)
    }
    refh <- 1
    if (length(grep("^(sum|max|avg)",opts$orderBy$what,perl=TRUE))>0) {
        nc <- nchar(opts$orderBy$what)
        rh <- substr(opts$orderBy$what,nc,nc)
        if (rh!="a") {
            rh <- suppressWarnings(as.numeric(rh))
            if (is.na(rh))
                warning("Reference profile for heatmap ordering not ",
                    "recognized! Using the 1st...",immediate.=TRUE)
            else
                refh <- rh
        }
        else
            refh <- 0 # Flag to indicate order by sum/max of all profiles
    }
    byMax <- bySum <- byAvg <- FALSE
    sorter <- list(ix=1:nrow(input[[1]]$profile))
    if (length(grep("^sum",opts$orderBy$what,perl=TRUE))>0)
        bySum <- TRUE
    if (length(grep("^max",opts$orderBy$what,perl=TRUE))>0)
        byMax <- TRUE
    if (length(grep("^avg",opts$orderBy$what,perl=TRUE))>0)
        byAvg <- TRUE
    
    if (bySum) {
        if (refh==0) {
            #tmp <- do.call("cbind",lapply(input,function(x,rc) {
            #   #if (is.null(x$coverage$center))
            #    theCov <- x$coverage
            #    #else
            #    #    theCov <- x$coverage$center
            #    s <- cmclapply(theCov,function(y) {
            #        if (is.null(y))
            #            return(0)
            #        return(sum(y))
            #    },rc=rc)
            #    return(unlist(s))
            #},rc=rc))
            tmp <- do.call("cbind",lapply(input,function(x) {
                return(apply(x$profile,1,sum))
            }))
            theVal <- apply(tmp,1,sum)
            names(theVal) <- rownames(input[[1]]$profile)
        }
        else #{
            #if (is.null(input[[refh]]$coverage$center))
            #theCov <- input[[refh]]$coverage
            #else
            #    theCov <- input[[refh]]$coverage$center
            #theVal <- unlist(cmclapply(theCov,function(y) {
            #    if (is.null(y))
            #        return(0)
            #    return(sum(y))
            #},rc=rc))
            #names(theVal) <- rownames(input[[refh]]$profile)
            theVal <- apply(input[[refh]]$profile,1,sum)
        #}
        if (opts$orderBy$order=="descending")
            sorter <- sort(theVal,decreasing=TRUE,index.return=TRUE)
        else
            sorter <- sort(theVal,index.return=TRUE)
    }
    if (byMax) {
        if (refh==0) {
            #tmp <- do.call("cbind",lapply(input,function(x,rc) {
            #    #if (is.null(x$coverage$center))
            #    theCov <- x$coverage
            #    #else
            #    #    theCov <- x$coverage$center
            #    s <- cmclapply(theCov,function(y) {
            #        if (is.null(y))
            #            y <- 0
            #        m <- max(y)
            #        mp <- which(y==m)
            #        if (length(mp)>1)
            #            return(as.numeric(y[sample(mp,1)]))
            #        else
            #            return(as.numeric(y[mp]))
            #    },rc=rc)
            #    return(unlist(s))
            #},rc=rc))
            tmp <- do.call("cbind",lapply(input,function(x) {
                return(apply(x$profile,1,function(y) {
                    m <- max(y)
                    mp <- which(y==m)
                    if (length(mp)>1)
                        return(y[sample(mp,1)])
                    else
                        return(y[mp])
                }))
            }))
            theVal <- apply(tmp,1,max)
            names(theVal) <- rownames(input[[1]]$profile)
        }
        else #{
            #if (is.null(input[[refh]]$coverage$center))
            #theCov <- input[[refh]]$coverage
            #else
            #    theCov <- input[[refh]]$coverage$center
            #theVal <- unlist(cmclapply(theCov,function(y) {
            #    if (is.null(y))
            #        y <- 0
            #    m <- max(y)
            #    mp <- which(y==m)
            #    if (length(mp)>1)
            #        return(as.numeric(y[sample(mp,1)]))
            #    else
            #        return(as.numeric(y[mp]))
            #},rc=rc))
            #names(theVal) <- rownames(input[[refh]]$profile)
            theVal <- apply(input[[refh]]$profile,1,function(y) {
                m <- max(y)
                mp <- which(y==m)
                if (length(mp)>1)
                    return(y[sample(mp,1)])
                else
                    return(y[mp])
            })
        #}
        if (opts$orderBy$order=="descending")
            sorter <- sort(theVal,decreasing=TRUE,index.return=TRUE)
        else
            sorter <- sort(theVal,index.return=TRUE)
    }
    if (byAvg) {
        if (refh==0) {
            tmp <- do.call("cbind",lapply(input,function(x) {
                return(apply(x$profile,1,mean))
            }))
            theVal <- apply(tmp,1,mean)
            names(theVal) <- rownames(input[[1]]$profile)
        }
        else 
            theVal <- apply(input[[refh]]$profile,1,mean)
        if (opts$orderBy$order=="descending")
            sorter <- sort(theVal,decreasing=TRUE,index.return=TRUE)
        else
            sorter <- sort(theVal,index.return=TRUE)
    }
    return(sorter)
}

orderProfilesByDesign <- function(input,design,opts,rc=NULL) {
    splitter <- split(1:nrow(design),design,drop=TRUE)
    sortlist <- sorter <- vector("list",length(splitter))
    names(sortlist) <- names(sorter) <- names(splitter)
    
    if (!is.null(opts$orderBy$custom)) {
        for (n in names(splitter)) {
            S <- splitter[[n]]
            if (opts$orderBy$order=="descending")
                sortlist[[n]] <- 
                    sort(opts$orderBy$custom[S],decreasing=TRUE,
                        index.return=TRUE)$ix
            else
                sortlist[[n]] <- sort(opts$orderBy$custom[S],
                    index.return=TRUE)$ix
            sorter[[n]] <- S[sortlist[[n]]]
        }
        return(unlist(sorter))
    }
    
    for (n in names(splitter)) {
        S <- splitter[[n]]
        subcov <- lapply(input,function(x,s) {
            #if (!is.null(x$coverage$center))
            #    return(x$coverage$center[s])
            #else
                #return(x$coverage[s])
                return(x$profile[s,])
                
        },S)
            
        refh <- 1
        if (length(grep("^(sum|max|avg)",opts$orderBy$what,perl=TRUE))>0) {
            nc <- nchar(opts$orderBy$what)
            rh <- substr(opts$orderBy$what,nc,nc)
            if (rh!="a") {
                rh <- suppressWarnings(as.numeric(rh))
                if (is.na(rh))
                    warning("Reference profile for heatmap ordering not ",
                        "recognized! Using the 1st...",immediate.=TRUE)
                else
                    refh <- rh
            }
            else
                refh <- 0 # Flag to indicate order by sum/max of all profiles
        }
        byMax <- bySum <- byAvg <- FALSE
        sortlist[[n]] <- 1:length(S)
        if (length(grep("^sum",opts$orderBy$what,perl=TRUE))>0)
            bySum <- TRUE
        if (length(grep("^max",opts$orderBy$what,perl=TRUE))>0)
            byMax <- TRUE
        if (length(grep("^avg",opts$orderBy$what,perl=TRUE))>0)
            byAvg <- TRUE
        
        if (bySum) {
            if (refh==0) {
                #tmp <- do.call("cbind",lapply(subcov,function(x,rc) {
                #    s <- cmclapply(x,function(y) {
                #        if (is.null(y))
                #            return(0)
                #        return(sum(y))
                #    },rc=rc)
                #    return(unlist(s))
                #},rc=rc))
                #theVal <- apply(tmp,1,sum)
                tmp <- do.call("cbind",lapply(subcov,function(x) {
                    return(apply(x,1,sum))
                }))
                theVal <- apply(tmp,1,sum)
            }
            else #{
                #theVal <- unlist(cmclapply(subcov[[refh]],function(y) {
                #    if (is.null(y))
                #        return(0)
                #    return(sum(y))
                #},rc=rc))
                theVal <- apply(subcov[[refh]],1,sum)
            #}
            if (opts$orderBy$order=="descending")
                sortlist[[n]] <- 
                    sort(theVal,decreasing=TRUE,index.return=TRUE)$ix
            else
                sortlist[[n]] <- sort(theVal,index.return=TRUE)$ix
        }
        if (byMax) {
            if (refh==0) {
                #tmp <- do.call("cbind",lapply(subcov,function(x,rc) {
                #    s <- cmclapply(x,function(y) {
                #        if (is.null(y))
                #            y <- 0
                #        m <- max(y)
                #        mp <- which(y==m)
                #        if (length(mp)>1)
                #            return(as.numeric(y[sample(mp,1)]))
                #        else
                #            return(as.numeric(y[mp]))
                #    },rc=rc)
                #    return(unlist(s))
                #},rc=rc))
                tmp <- do.call("cbind",lapply(subcov,function(x) {
                    return(apply(x,1,function(y) {
                        m <- max(y)
                        mp <- which(y==m)
                        if (length(mp)>1)
                            return(y[sample(mp,1)])
                        else
                            return(y[mp])
                    }))
                }))
                theVal <- apply(tmp,1,max)
            }
            else {
                #theVal <- unlist(cmclapply(subcov[[refh]],function(y) {
                #    if (is.null(y))
                #        y <- 0
                #    m <- max(y)
                #    mp <- which(y==m)
                #    if (length(mp)>1)
                #        return(as.numeric(y[sample(mp,1)]))
                #    else
                #        return(as.numeric(y[mp]))
                #},rc=rc))
                theVal <- apply(subcov[[refh]],1,function(y) {
                    m <- max(y)
                    mp <- which(y==m)
                    if (length(mp)>1)
                        return(y[sample(mp,1)])
                    else
                        return(y[mp])
                })
            }
            if (opts$orderBy$order=="descending")
                sortlist[[n]] <- 
                    sort(theVal,decreasing=TRUE,index.return=TRUE)$ix
            else
                sortlist[[n]] <- sort(theVal,index.return=TRUE)$ix
        }
        if (byAvg) {
            if (refh==0) {
                tmp <- do.call("cbind",lapply(subcov,function(x) {
                    return(apply(x,1,mean))
                }))
                theVal <- apply(tmp,1,mean)
            }
            else
                theVal <- apply(subcov[[refh]],1,mean)
            if (opts$orderBy$order=="descending")
                sortlist[[n]] <- 
                    sort(theVal,decreasing=TRUE,index.return=TRUE)$ix
            else
                sortlist[[n]] <- sort(theVal,index.return=TRUE)$ix
        }
        sorter[[n]] <- S[sortlist[[n]]]
   }
   return(unlist(sorter))
}

orderSignals <- function(input,opts) {
    if (!is.null(opts$orderBy$custom)) {
        if (opts$orderBy$order=="descending")
            sorter <- sort(opts$orderBy$custom,decreasing=TRUE,
                index.return=TRUE)
        else
            sorter <- sort(opts$orderBy$custom,index.return=TRUE)
        return(sorter)
    }
    refh <- 1
    if (length(grep("^(sum|max|avg)",opts$orderBy$what,perl=TRUE))>0) {
        nc <- nchar(opts$orderBy$what)
        rh <- substr(opts$orderBy$what,nc,nc)
        rh <- suppressWarnings(as.numeric(rh))
        if (is.na(rh))
            warning("Reference profile for heatmap ordering not ",
                "recognized! Using the 1st...",immediate.=TRUE)
        else
            refh <- rh
    }
    sorter <- list(ix=1:length(input[[1]]$profile))
    if (opts$orderBy$order=="descending")
        sorter <- sort(input[[refh]]$profile,decreasing=TRUE,index.return=TRUE)
    else
        sorter <- sort(input[[refh]]$profile,index.return=TRUE)
    return(sorter)
}

orderDesignSignals <- function(input,design,opts) {
    splitter <- split(1:nrow(design),design,drop=TRUE)
    sortlist <- sorter <- vector("list",length(splitter))
    names(sortlist) <- names(sorter) <- names(splitter)
    
    if (!is.null(opts$orderBy$custom)) {
        for (n in names(splitter)) {
            S <- splitter[[n]]
            if (opts$orderBy$order=="descending")
                sortlist[[n]] <- 
                    sort(input[[refh]]$profile[S],decreasing=TRUE,
                        index.return=TRUE)$ix
            else
                sortlist[[n]] <- sort(input[[refh]]$profile[S],
                        index.return=TRUE)$ix
            sorter[[n]] <- S[sortlist[[n]]]
        }
        return(unlist(sorter))
    }

    for (n in names(splitter)) {
        S <- splitter[[n]]
        refh <- 1
        if (length(grep("^(sum|max|avg)",opts$orderBy$what,perl=TRUE))>0) {
            nc <- nchar(opts$orderBy$what)
            rh <- substr(opts$orderBy$what,nc,nc)
            rh <- suppressWarnings(as.numeric(rh))
            if (is.na(rh))
                warning("Reference profile for heatmap ordering not ",
                    "recognized! Using the 1st...",immediate.=TRUE)
            else
                refh <- rh
        }
        sortlist[[n]] <- list(ix=1:length(input[[1]]$profile[S]))
        if (opts$orderBy$order=="descending")
            sortlist[[n]] <- sort(input[[refh]]$profile[S],decreasing=TRUE,
                index.return=TRUE)$ix
        else
            sortlist[[n]] <- sort(input[[refh]]$profile[S],index.return=TRUE)$ix
        sorter[[n]] <- S[sortlist[[n]]]
    }
    return(unlist(sorter))
}

makeHorizontalAnnotation <- function(width,opts,type=c("profile","heatmap")) {
    fl <- opts$xAxisParams$flank
    flb <- fl
    fl[1] <- -fl[1]
    if (type=="heatmap" && opts$binParams$forceHeatmapBinning && !all(fl==0))
        opts$binParams$flankBinSize <- opts$binParams$forcedBinSize[1]
    if (!any(fl==0))
        edgeLabels <- paste(round(fl/1000,1),"kb",sep="")
    else {
        if (fl[1]==0 && fl[2]!=0)
            edgeLabels <- paste(round(fl[2]/1000,1),"kb",sep="")
        else if (fl[1]!=0 && fl[2]==0)
            edgeLabels <- paste(round(fl[1]/1000,1),"kb",sep="")
        else
            edgeLabels <- NULL
    }
    switch(opts$xAxisParams$region,
        tss = {
            midLabels <- "TSS"
            breaks <- c(1,round(width/2,1),width)
        },
        tes = {
            midLabels <- "TES"
            breaks <- c(1,round(width/2,1),width)
        },
        genebody = {
            midLabels <- c("TSS","TES")
            if (opts$binParams$flankBinSize==0)
                breaks <- c(
                    1,
                    abs(opts$xAxisParams$flank[1]),
                    width-opts$xAxisParams$flank[2],
                    width
                )
            else {
                f <- flb/max(flb)
                r <- flb/sum(flb)
                rdiff <- round(abs(opts$binParams$regionBinSize - 
                    (width - opts$binParams$flankBinSize*f[1] - 
                    opts$binParams$flankBinSize*f[2])))
                breaks <- c(
                    1,
                    round(opts$binParams$flankBinSize*f[1] + rdiff*r[1]),
                    round(width-(opts$binParams$flankBinSize*f[2] + 
                        rdiff*r[2])),
                    #round(opts$binParams$flankBinSize*f[1]),
                    #round(width-(opts$binParams$flankBinSize*f[2])),
                    width
                )
            }
        },
        custom = {
            if (opts$xAxisParams$customIsBase) {
                midLabels <- "Center"
                breaks <- round(c(width/8,width/2,width-width/8),1)
            }
            else {
                midLabels <- c("Start","End")
                if (opts$binParams$flankBinSize==0)
                    breaks <- c(
                        1,
                        abs(opts$xAxisParams$flank[1]),
                        width-opts$xAxisParams$flank[2],
                        width
                    )
                else {
                    f <- flb/max(flb)
                    r <- flb/sum(flb)
                    rdiff <- round(abs(opts$binParams$regionBinSize - 
                        (width - opts$binParams$flankBinSize*f[1] - 
                        opts$binParams$flankBinSize*f[2])))
                    breaks <- c(
                        1,
                        round(opts$binParams$flankBinSize*f[1] + rdiff*r[1]),
                        round(width-(opts$binParams$flankBinSize*f[2] + 
                            rdiff*r[2])),
                        width
                    )
                }
            }
        }
    )
    if (!any(fl==0))
        labels <- c(edgeLabels[1],midLabels,edgeLabels[2])
    else {
        if (fl[1]==0 && fl[2]!=0) {
            labels <- c(midLabels,edgeLabels)
            breaks <- breaks[c(1,3,4)]
        }
        else if (fl[1]!=0 && fl[2]==0) {
            labels <- c(edgeLabels,midLabels)
            breaks <- breaks[c(1,2,4)]
        }
        else {
            labels <- midLabels
            breaks <- breaks[c(1,4)]
        }
    }
    return(list(breaks=breaks,labels=labels))
}

graphicsOpen <- function(o,f,...) {
    if(o!="x11" && is.null(f))
        stop("Please specify an output file name for your plot")
    switch(o,
        x11 = { dev.new(...) },
        pdf = { pdf(file=f,pointsize=10,...) },
        ps = { postscript(file=f,pointsize=10,...) },
        png = { png(filename=f,pointsize=12,...) },
        jpg = { jpeg(filename=f,pointsize=12,quality=100,...) },
        bmp = { bmp(filename=f,pointsize=12,...) },
        tiff = { tiff(filename=f,pointsize=12,...) }
    )
}

graphicsClose <- function(o) {
    if (!is.element(o,c("x11","png","jpg","tiff","bmp","pdf","ps")))
        return(FALSE)
    if (o!="x11")
        dev.off()
}
