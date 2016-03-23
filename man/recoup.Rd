\name{recoup}
\alias{recoup}
\title{Create genomic signal profiles in predefined or 
    custom areas using short sequence reads}
\usage{
    recoup(
        input,
        design = NULL,
        region = c("genebody", "tss", "tes", "custom"),
        type = c("chipseq", "rnaseq"),
        genome = c("hg18", "hg19", "hg38", "mm9" ,"mm10",
            "rn5", "dm3", "danrer7", "pantro4", "susscr3"),
        refdb = c("ensembl", "ucsc", "refseq"),
        flank = c(2000, 2000),
        fraction = 1,
        orderBy = list(
            what = c("none", "suma", "sumn", 
            "maxa","maxn", "avga", "avgn", "hcn"), 
            order = c("descending", "ascending"),
            custom = NULL
        ),
        binParams = list(
            flankBinSize = 0,
            regionBinSize = 0,
            sumStat = c("mean", "median"),
            smooth = TRUE,
            interpolation = c("auto", "spline", "linear", 
                "neighborhood"),
            forceHeatmapBinning = TRUE,
            forcedBinSize = c(50, 200)
        ),
        selector = NULL,
        preprocessParams = list(
            normalize = c("none", "linear", 
                "downsample","sampleto"),
            sampleTo = 1e+6,
            spliceAction = c("split", "keep", "remove"),
            spliceRemoveQ = 0.75,
            bedGenome = ifelse(genome \%in\% c("hg18", "hg19",
                "hg38", "mm9", "mm10", "rn5", "dm3", 
                "danrer7", "pantro4", "susscr3"), genome, 
                NA),
            seed = 42
        ),
        plotParams = list(
            plot = TRUE,
            profile = TRUE,
            heatmap = TRUE,
            correlation = TRUE,
            signalScale = c("natural", "log2"),
            heatmapScale = c("common", "each" ),
            heatmapFactor = 1,
            corrScale = c("normalized", "each"),
            corrSmoothPar = ifelse(is.null(design), 0.1, 
                0.5),
            singleFacet = c("none", "wrap", "grid"),
            multiFacet = c("wrap", "grid"),
            conf = TRUE,
            device = c("x11", "png", "jpg", "tiff", "bmp",
                "pdf", "ps"),
            outputDir = ".",
            outputBase = NULL
        ),
        saveParams = list(
            ranges = TRUE,
            coverage = TRUE,
            profile = TRUE,
            profilePlot = TRUE,
            heatmapPlot = TRUE,
            correlationPlot = TRUE
        ),
        kmParams = list(
            k = 0,
            nstart = 20,
            algorithm = c("Hartigan-Wong", 
                "Lloyd", "Forgy", "MacQueen"),
            iterMax = 20,
            reference = NULL,
            seed = 42
        ),
        strandedParams = list(
            strand = NULL, 
            ignoreStrand = TRUE
        ),
        ggplotParams = list(
            title = element_text(size = 12),
            axis.title.x = element_text(size = 10, 
                face = "bold"),
            axis.title.y = element_text(size = 10, 
                face = "bold"),
            axis.text.x = element_text(size = 9, 
                face = "bold"),
            axis.text.y = element_text(size = 10, 
                face = "bold"),
            strip.text.x = element_text(size = 10, 
                face = "bold"),
            strip.text.y = element_text(size = 10, 
                face = "bold"),
            legend.position = "bottom",
            panel.margin = grid::unit(1, "lines")
        ),
        complexHeatmapParams = list(
            main = list(
                cluster_rows = ifelse(length(grep(
                    "hc", orderBy$what)) > 0, TRUE, FALSE),
                cluster_columns = FALSE,
                column_title_gp = grid::gpar(fontsize = 10, 
                    font = 2),
                show_row_names = FALSE,
                show_column_names = FALSE,
                heatmap_legend_param = list(
                    color_bar = "continuous"
                )
            ),
            group=list(
                cluster_rows = ifelse(length(grep(
                    "hc", orderBy$what)) > 0, TRUE, FALSE),
                cluster_columns = FALSE,
                column_title_gp = grid::gpar(fontsize = 10, 
                    font = 2),
                show_row_names = FALSE,
                show_column_names = FALSE,
                row_title_gp = grid::gpar(fontsize = 8,
                    font = 2),
                gap = unit(5, "mm"),
                heatmap_legend_param = list(
                    color_bar = "continuous"
                )
            )
        ),
        bamParams = NULL,
        onTheFly = FALSE,
        localDbHome = file.path(path.expand("~"), ".recoup"), 
        rc = NULL
    )
}
\arguments{
    \item{input}{the main input to \code{recoup} can be 
    either a list or a configuration file (with essentially
    the same contents as the list). In case of list input,
    it is a list of n lists, where n the number of samples.
    See Details for the inner list contents. Alternatively,
    \code{input} can be a text tab delimited file with a
    specific header (the same fields as each inner list
    when \code{input} is a list) and one row for each
    sample. Again, see Details section for the field
    specifications.}

    \item{design}{either a data frame with grouping factors 
    as columns (e.g. two grouping factors can be strand,  
    and Ensembl biotype) or a tab delimited text file with  
    the same content (grouping factors in columns). If a  
    data frame, the \code{row.names} attribute must  
    correspond to the names (e.g. rownames) of the 
    \code{genome} argument, or be a superset or subset of
    them. If a file, the first column must correspond to 
    the names (e.g. rownames) of the \code{genome} 
    argument or be a superset or subset of them.}

    \item{region}{one of \code{"tss"}, \code{"tes"}, 
    \code{"genebody"}, \code{"custom"}.}
    
    \item{type}{one of \code{"chipseq"}, \code{"rnaseq"}.}
    
    \item{genome}{when \code{region} is \code{"tss"}, 
    \code{"tes"} or \code{"genebody"}, \code{genome} can be 
    one of \code{"hg38"}, \code{"hg19"}, \code{"hg18"}, 
    \code{"mm10"}, \code{"mm9"}, \code{"dm3"}, 
    \code{"rn5"}, \code{"danrer7"}, \code{"pantro4"}, 
    \code{"susscr3"} for human, mouse, fruitfly, rat, 
    zebrafish and chimpanzee genomes respectively. When 
    when \code{region} is one of the above or 
    \code{"custom"}, \code{genome} can be a tab delimited 
    BED-like text file or a data frame.}
    
    \item{refdb}{one of \code{"ensembl"}, \code{"ucsc"} or 
    \code{"refseq"}. It will be used to retrieve genomic 
    reference regions when \code{genome} argument is one of
    the supported organisms.}

    \item{flank}{a vector of length two with the number of
    base pairs to flank upstream and downstream the 
    \code{region}. Minimum flank is 0bp and maximum is
    50kb. It is always expressed in bp.}
    
    \item{fraction}{a number from 0 to 1 (default) denoting
    the fraction of total data to be used. See Details for
    further information.}

    \item{orderBy}{a named list whose members control the
    order of the genomic regions (related to the 
    \code{genome} and \code{region} arguments as they 
    appear in heatmap profiles. The list has the following
    fields:
    \itemize{
        \item \code{what}: one of \code{"none"} (default),
        \code{"suma"}, \code{"sumn"}, \code{"maxa"},
        \code{"maxn"}, \code{"avga"}, \code{"avgn"}, 
        \code{"hcn"}, where \code{n} in \code{"sumn"},
        \code{"maxn"}, \code{"hcn"} is the index of the 
        profile which sould be used as reference. See 
        Details for further information.
        \item \code{order}: either \code{"descending"}
        (default) for ordering coverages from highest to 
        lowest, or \code{"ascending"} for the opposite.
        \item \code{custom}: a numeric vector of custom
        values (e.g. RNA abundance) that will be used to
        sort all the profiles. If provided, \code{what}
        will be ignored. Defaults to \code{NULL}.
    }}
    
    \item{binParams}{a named list whose members control the
    resolution of the coverage profiles. The list has the 
    following fields:
    \itemize{
        \item \code{flankBinSize}: the number of intervals
        (bins) into which the upstream and downstream 
        regions are split and the per-base coverage is 
        averaged across. If \code{0} (default), no binning
        is performed and the profiles are calculated based
        at the base-pair level (the highest possible 
        resolution).
        \item \code{regionBinSize}: the number of intervals
        (bins) into which the main region is split and 
        the per-base coverage is averaged across. 
        If \code{0} (default), no binning is performed 
        and the profiles are calculated based
        at the base-pair level (the highest possible 
        resolution).
        \item \code{sumStat}: the statistic which is used
        to summarize the bin coverage. Can be \code{"mean"}
        (default) or \code{"median"}.
        \item \code{smooth}: if \code{TRUE} (default), the
        final curve profiles are smoothed using splines.
        Set to \code{FALSE} for no smoothing. If the 
        reference genomic regions are many, the differences
        are minimal.
        \item \code{interpolation}: the interpolation 
        method to be used for coverage interpolation 
        when the reference regions are of unequal lengths
        (e.g. gene bodies) and the \code{regionBinSize}
        is larger than some of the former. Can be 
        \code{"auto"} (default), \code{"spline"}, 
        \code{"linear"} or \code{"neighborhood"}. See 
        Details for further explanations of each option.
        \item \code{forceHeatmapBinning}: if \code{TRUE}
        (default) and the profile resolution is very high
        (see \code{flankBinSize} and \code{regionBinSize}
        above), binning is applied prior to heatmap profile
        generation, otherwise the heatmaps will be oversized 
        and will take a lot of time to render. Set to 
        \code{FALSE} if both \code{flankBinSize} and 
        \code{regionBinSize} are not zero so as to avoid
        unecessary profile recalculations.
        \item \code{forcedBinSize}: a vector with two
        integers representing the \code{flankBinSize} and 
        \code{regionBinSize} to be used with 
        \code{forceHeatmapBinning} above.
    }
    See Details for a few further notes in the usage of 
    \code{binParams}.}
    
    \item{selector}{a named list whose members control 
    some subsetting abilities regarding the input reference
    genomic regions (\code{genome} argument) or NULL 
    (default) when the \code{genome} argument is/may be 
    custom. When list, it has the following fields:
    \itemize{
        \item \code{id}: a vector of ids of the same type
        as those present in the \code{genome} file or
        organism type/version.
        \item \code{bioype}: a vector of Ensembl biotypes
        that will be used to filter the \code{genome} when
        the latter is one of the supported organisms. Not
        used when \code{genome} is a custom file.
        \item \code{exonType}: currently not used.
    }}
    
    \item{preprocessParams}{a named list whose members 
    control certain preprocessing steps applied to the 
    \code{\link{GRanges}} objects obtained while or after
    reading the input BAM/BED files with short reads. The 
    list has the following fields:
    \itemize{
        \item \code{normalize}: one of \code{"none"}
        (default), \code{"linear", "downsample", 
        "sampleto"}. Controls how the coverages are 
        normalized across samples. See Details for 
        explanation of these options.
        \item \code{sampleTo}: a fixed library size for
        downsampling to be used with \code{"sampleto"}
        option above.
        \item \code{spliceAction}: one of \code{"split"}
        (default), \code{"keep", "remove"}. Controls the 
        action to be performed with spliced reads in the 
        case of RNA-Seq samples. See Details for 
        explanation of these options.
        \item \code{spliceRemoveQ}: the quantile of
        putative joint spliced read length to be used
        for read filtering when \code{spliceAction} is
        \code{"remove"}. See Details for further 
        explanations.
        \item \code{bedGenome}: one of the supported
        genomes, as when reading from bed files, 
        chromosomal lenghts are not available and must
        be retrieved with another way.
        \item \code{seed}: seed for random number 
        generation to be used when \code{spliceAction} is 
        \code{"downsample"} or \code{"sampleto"}. See 
        Details for further explanations. Defaults to 42.
    }}

    \item{plotParams}{a named list whose members control
    profile (curve and heatmap) plotting parameters. The 
    list has the following fields:
    \itemize{
        \item \code{plot}: if set to \code{TRUE} 
        (default), the plots created with the calculated
        profiles with \code{recoup} are displayed. Set to
        \code{FALSE} to plot later using the output object.
        \item \code{profile}: if set to \code{TRUE} 
        (default), the average coverage profile across
        the genomic regions of preference is calculated.
        Set to \code{FALSE} to suppress this.
        \item \code{heatmap}: if set to \code{TRUE} 
        (default), the coverage heatmap profile across
        the genomic regions of preference is calculated.
        Set to \code{FALSE} to suppress this.
        \item \code{correlation}: if set to \code{TRUE} 
        (default), the plots created with the calculated
        coverage correlations are displayed. Set to
        \code{FALSE} to plot later using the output object.
        \item \code{signalScale}: one of \code{"natural"}
        (default) or \code{"log2"} to control the signal
        scale of the final coverage plots. Hint: use
        \code{log2} scale for RNA-Seq profiles as it 
        produces much smoother plots.
        \item \code{heatmapScale}: one of \code{"common"}
        (default) or \code{"each"}. When \code{"common"},
        a common heatmap color scale is calculated for 
        all samples. When \code{"each"}, each heatmap has
        its own color scale.
        \item \code{heatmapFactor}: a positive numeric 
        value by which the upper color scale limit of
        the heatmap profile is multiplied. Defaults to 
        \code{1}. See Details for further information.
        \item \code{corrScale}: either \code{"normalized"}
        (default) or \code{"each"}. Controls the scale
        display in coverage correlation plots. See Details 
        for further information.
        \item \code{corrSmoothPar}: a numeric value between
        0 and 1 which controls the smoothing of correlation
        plots. Its default value is controlled by the 
        presence of \code{design}. See Details for further 
        information.
        \item \code{corrScale}: either \code{"normalized"}
        (default) or \code{"each"}. See Details for further
        information.
        \item \code{singleFacet}: how should \code{ggplot2}
        should facet the profiles with 1-factor design and
        only one sample whose profile will be plotted. Can
        be \code{"none"} (default), \code{"wrap"} or 
        \code{"grid"}. When \code{"none"}, no gridding is
        applied and design factors are distinguished by 
        colour. With more than one design factors, the 
        \code{multiFacet} option below is used.
        \item \code{multiFacet}: how should \code{ggplot2}
        should facet the profiles with 1-factor design and
        more than one samples whose profile will be plotted.
        Can be \code{"wrap"} (default) or \code{"grid"}.
        2 or 3 (3rd would be colour) factor designs are 
        faceted with \code{"grid"}.
        \item \code{conf}: plot also confidence intervals
        using \code{geom_ribbon} in profile or correlation
        plots.
        \item \code{device}: the R plotting device to 
        redirect the plots to. Can be 
        \code{"x11"} (default), \code{"png", "jpg", 
        "tiff", "bmp", "pdf", "ps"}.
        \item \code{outputDir}: the directory to place
        profiles when the plotting device is not 
        \code{"x11"}. Defaults to \code{"."}.
        \item \code{outputBase}: the naming template
        for output files when the plotting device is 
        not \code{"x11"}. The extensions \code{"_profile"}
        and \code{"_heatmap"} will be appended to 
        distinguish each plot type. Leave \code{NULL}
        (default) for automatic filename generation.
    }}

    \item{saveParams}{a named list which controls the
    information to be stored in the \code{recoup} output
    list object. The list has the following fields:
    \itemize{
        \item \code{ranges}: set to \code{TRUE} (default)
        to store the \code{GRanges} object obtained from 
        the BAM/BED files. Set to \code{FALSE} for not
        saving.
        \item \code{coverage}: set to \code{TRUE} (default)
        to store the \code{Rle} list object obtained from 
        the coverage calculations. Set to \code{FALSE} for 
        not saving.
        \item \code{profile}: set to \code{TRUE} (default)
        to store the profile matrices exracted from 
        coverage summarizations. Set to \code{FALSE} for 
        not saving. It must be present when using the 
        \code{recoup} output in the plotting functions
        \code{\link{recoupProfile}} and 
        \code{\link{recoupHeatmap}}.
        \item \code{profilePlot}: set to \code{TRUE} 
        (default) to store \code{ggplot} object containing
        the average coverage plot. Set to \code{FALSE} for 
        not saving. Must be \code{TRUE} if you wish to use
        the \code{recoup} output later with 
        \code{\link{recoupPlot}}.
        \item \code{heatmapPlot}: set to \code{TRUE} 
        (default) to store \code{ComplexHeatmap} object 
        containing the coverage heatmap plot. Set to 
        \code{FALSE} for not saving. Must be \code{TRUE} 
        if you wish to use the \code{recoup} output later 
        with \code{\link{recoupPlot}}.
        \item \code{correlationPlot}: set to \code{TRUE} 
        (default) to store \code{ggplot} object containing
        coverage correlation plot. Set to \code{FALSE} for 
        not saving. Must be \code{TRUE} if you wish to use 
        the \code{recoup} output later with 
        \code{\link{recoupPlot}}.
    }
    See the Details section for some additional 
    information.}
    
    \item{kmParams}{a named list which controls the 
    execution of k-means clustering using standard R base
    function \code{\link{kmeans}} otherwise. The list has
    the following fields:
    \itemize{
        \item \code{k}: the number of clusters for k-means
        clustering. When \code{0} (default), no k-means
        clustering is performed.
        \item \code{nstart}: See \code{\link{kmeans}}.
        \item \code{algorithm}: See \code{\link{kmeans}}.
        \item \code{iterMax}: See the \code{iter.max} 
        parameter in \code{\link{kmeans}}.
        \item \code{reference}: which profile to use as 
        reference for the determination of clusters and 
        ordering. The rest of the heatmaps will be ordered
        according to the reference clustering. It can be 
        either a sample id or \code{NULL} (default). If 
        the latter, all the profile matrices are merged 
        into one big matrix and k-means clustering is 
        performed on that matrix.
        \item \code{seed}: See \code{\link{kmeans}}.
    }
    The result of k-means clustering will be appended to
    \code{design} as an additional field. If \code{design}
    is \code{NULL}, it will be created and passed to the
    plotting functions.}
    
    \item{strandedParams}{a named list which controls how
    strand information will be treated (if present). The 
    list has the following fields:
    \itemize{
        \item \code{strand}: if set to \code{NULL} (default)
        then reads from both strands are used from the input
        BAM/BED files. If \code{"+"} or \code{"-"}, then
        only the respective strands are used.
        \item \code{ignoreStrand}: \code{TRUE} (default)
        or \code{FALSE}. Passed to the \code{ignore.strand}
        argument in the \code{\link{findOverlaps}} function
        used during coverage calculations.
    }}
    
    \item{ggplotParams}{a named list with theme parameters
    passed to the \code{\link{ggplot}} function of the
    \code{ggplot2} package. See the documentation of 
    \code{ggplot2} for further details. Only the parameters
    mentioned in the function call are used.}
    
    \item{complexHeatmapParams}{a named list with groups
    of parameters passed to the \code{\link{Heatmap}}
    function of the \code{ComplexHeatmap} package. The
    list has the following fields:
    \itemize{
        \item \code{main}: \code{\link{Heatmap}} parameters
        applied to each non-split (according to 
        \code{design}) heatmap. See the \code{recoup} 
        function call for supported parameters and 
        \code{\link{Heatmap}} for further details.
        \item \code{main}: \code{\link{Heatmap}} parameters
        applied to each split (according to \code{design}) 
        heatmap. See the \code{recoup} function call for 
        supported parameters and \code{\link{Heatmap}} 
        for further details.
    }}

    \item{bamParams}{BAM file read parameters passed to 
    \code{BamFile}. See the related function.
    Currently this is not used.}

    \item{onTheFly}{Read short reads directly from BAM files
    when input contains paths to BAM files. In this case the
    storage of short reads in the output list object as a 
    GRanges object is not possible. This argument is not
    used in this release and is a placeholder for the next
    release.}
    
    \item{localDbHome}{local path to retrieve annotation 
    from when \code{genome} is one of the supported 
    organisms. See also \code{\link{buildAnnotationStore}}.}
    
    \item{rc}{fraction (0-1) of cores to use in a multicore 
    system. It defaults to \code{NULL} (no parallelization).}
}
\value{
    a named list with five members:
    \itemize{
        \item \code{data}: the \code{input} argument if it
        was a list or the resulting list from the unexported
        internal \code{readConfig} function, with the 
        \code{ranges}, \code{coverage} and \code{profile}
        fields filled according to \code{saveParams}. This
        \code{data} member can be used again as an argument
        to \code{recoup}. The \code{coverage} and 
        \code{profile} fields will be recalculated according
        to \code{recoup} parameters but the \code{ranges} 
        will be resued if the input files are not changed.
        \item \code{design}: the design data frame which is
        used to facet the profiles.
        %\item \code{plotopts}: plotting options read from
        %\code{\link{recoupProfile}} and 
        %\code{\link{recoupHeatmap}} functions in order to
        %create the respective plots. It is a simple list 
        %object and users can modify the parameters which 
        %are mostly self-explanatory.
        \item \code{plots}: the \code{ggplot2} and/or
        \code{Heatmap} objects created by \code{recoup}.
        \item \code{callopts}: the majority of \code{recoup}
        call parameters. Their storage serves the reuse of
        a \code{recoup} list object so that only certain 
        elements of plots are recalculated.
    }
}
\description{
    This function calculates and plots signal profiles 
    created from short sequence reads derived from Next 
    Generation Sequencing technologies. The profiles 
    provided are either sumarized curve profiles or 
    heatmap profiles. Currently, \code{recoup} supports 
    genomic profile plots for reads derived from ChIP-Seq 
    and RNA-Seq experiments. The function uses ggplot2 and
    ComplexHeatmap graphics facilities for curve and 
    heatmap coverage profiles respectively. The output
    list object can be reused as input to this function
    which will automatically recognize which profile 
    elements needto be recalculated, saving time.
}
\details{
    When \code{input} is a list, is should contain as many
    sublists as the number of samples. Each sublist must 
    have at least the following fields:
    \itemize{ 
        \item \code{id}: a unique identifier for each 
        sample which should not contain whitespaces
        and preferably no special characters.
        \item \code{file}: the full path to the BAM or
        BED file. If the path to the BAM is a hyperlink,
        the BAM file must be indexed.
        \item \code{format}: either \code{"bam"} or
        \code{"bed"}.
    }
    Additionally, each sublist may also contain the 
    following fields:
    \itemize{
        \item \code{name}: a sample name which will 
        appear in plots.
        \item \code{color}: either an R color (see the 
        \code{\link{colors}} function) or a hexadecimal
        color (e.g. \code{"#FF0000"}).
    }
    When \code{input} is a text file, this should be 
    strictly tab-delimited (no other delimiters like 
    comma), should contain a header with the same names
    (case sensitive) as the sublist fields above (\code{
    id, file, format} are mandatory and \code{name, color}
    are optional).
    
    When \code{genome} is not one of the supported 
    organisms, it should be a text tab delimited file
    (only tabs supported) with a header line, or a 
    data frame, where the basic BED field must be 
    present, that means that at least \code{"chromosome"}, 
    \code{"start"}, \code{"end"}, a unique identifier 
    column and \code{"strand"} must be present, \
    preferably in this order. This file is read in a 
    data frame and then passed to the 
    \code{\link{makeGRangesFromDataFrame}} function 
    from the \code{GenomicRanges} package which takes 
    care of the rest. See also the 
    \code{\link{makeGRangesFromDataFrame}}'s documentation.
    When \code{genome} is one of the supported organisms,
    \code{recoup} takes care of the rest.
    
    When \code{region} is \code{"tss"}, the curve and 
    heatmap profiles are centered around the TSS of the
    (gene) regions provided with the \code{genome} 
    argument, flanked accordind to the \code{flank} 
    argument. The same applied for \code{region="tes"}
    where the plots are centered around the transcription
    end site. When \code{region} is \code{"genebody"},
    the profiles consist of two flanking parts (upstream
    of the TSS and downstream of the TES) and a middle
    part consisting of the gene body coverage profile.
    The latter is constructed by creating a fixed number
    of intervals (bins) along each gene and averaging
    the coverage of each interval. In some extreme cases
    (e.g. for small genes), the number of bins may be
    larger than the gene length. In these cases, a few
    zeros are distributed randomly across the number of
    bins to reach the predefined number of gene body 
    intervals. When \code{region} is \code{"custom"}
    the behavior depends on the custom regions length.
    If it contains single-base intervals (e.g. ChIP-Seq 
    peak centers), then the behavior is similar to the
    TSS behavior above. If it contains genomic intervals
    of equal or unequal size, the behavior is similar
    to the gene body case.
    
    The \code{fraction} parameter controls the total 
    fraction of both total reads and genomic regions
    to be used for profile creation. This means that
    the total reads for each sample are randomly
    downsampled to \code{fraction*100\%} of the 
    original reads and the same applied to the input
    genomic regions. This practice is followed by 
    similar packages (like ngs.plot) and serves the
    purpose of a quick overview of how the actual 
    profiles look before profiling the whole genome.
    
    Regarding the \code{orderBy} parameters, for the
    options of the \code{what} parameter \code{"sum"} 
    type of options order profiles according to i) the 
    sum of coverages of all samples in each genomic
    region when \code{orderBy$what="suma"} or ii) the
    sum of coverages of sample n (e.g. 2) in each 
    genomic region when \code{orderBy$what="sumn"}
    (e.g. \code{orderBy$what="sum1"}). The same apply
    for the \code{"max"} type of options but this time
    the ordering is performed according to the position
    of the highest coverage in each genomic profile. Ties
    in the position of highest coverage are broken 
    randomly and sorting is performed with the default R 
    \code{\link{sort}}. Similarly for \code{"avg"} type
    of options, the ordering is performed according to
    the average total coverage of a reference region.
    For the \code{"hc"} type of options, hierarchical 
    clustering is performed on the selected (n) 
    reference profile (e.g. \code{orderBy$what="hc1"}) 
    and this ordering is applied to the rest of the 
    sample profiles. When \code{what="none"}, no 
    ordering is performed and the input order is used 
    (\code{genome} argument). If any design is present 
    through the \code{design} argument or k-means 
    clustering is also performed (through the 
    \code{kmParams} argument), the orderBy directives 
    are applied to each sub-profile created by design 
    or k-means clustering.
    
    Regarding the \code{flankBinSize} field of 
    \code{binParams}, it is used only when 
    \code{region="genebody"} or \code{region="custom"}
    and the custom regions are not single-base regions.
    This happens as when the genomic regions to be 
    profiled are single-base regions (e.g. TSSs or
    ChIP-Seq peak centers), these regions are merged
    with the flanking areas and alltogether form the
    main genomic region. In these cases, only the 
    \code{regionBinSize} field value is used. Note that
    when \code{type="rnaseq"} or \code{region="genebody"} 
    or \code{region="custom"} with non single-base regions 
    the values of \code{flankBinSize} and 
    \code{regionBinSize} offer a fine control over how
    the flanking and the main regions are presented in the
    profiles. For example, when \code{flankBinSize=100} 
    and \code{regionBinSize=100} with a gene body profile
    plot, the outcome will look kind of "unrealistic" as 
    the e.g. 2kb flanking regions will look very similar
    to the usually larger gene bodies. On the other hand
    if \code{flankBinSize=50} and \code{regionBinSize=200},
    this setting will create more "realistic" gene body
    profiles as the flanking regions will be squished and
    the gene body area will look expanded. Within the same
    parameter group is also \code{interpolation}. When 
    working with reference regions of different lengths
    (e.g. gene bodies), it happens very often that their 
    lengths are a little to a lot smaller than the number 
    of bins into which they should be split and averaged 
    in order to be able to create the average curve and 
    heatmap profiles. \code{recoup} allows for dynamic
    resolutions by permitting to the user to set the 
    number of bins into which genomic areas will be 
    binned or by allowing a per-base resolution where 
    possible. The \code{interpolation} parameter controls
    what happens in such cases. When \code{"spline"}, the
    R function \code{\link{spline}} is used, with the
    default method, to produce a spline interpolation
    of the same size as the \code{regionBinSize} option
    and is used as the coverage for that region. When 
    \code{"linear"}, the procedure is the same as above
    but using \code{\link{approx}}. When 
    \code{"neighborhood"}, a number of \code{NA} values 
    are distributed randomly across the small area 
    coverage vector, excluding the first and the last
    two positions, in order to reach \code{regionBinSize}.
    Then, each \code{NA} position is filled with the mean
    value of the two values before and the two values
    after the \code{NA} position, with \code{na.rm=TRUE}.
    This method should be avoided when >20\% of the values
    of the extended vector are \code{NA}'s as it may cause
    a crash. However, it should be the most accurate one
    in the opposite case (few \code{NA}'s). When 
    \code{"auto"} (the default), a hybrid of 
    \code{"spline"} and \code{"neighborhood"} is applied.
    If the \code{NA}'s constitute more than 20\% of the 
    extended vector, \code{"spline"} is used, otherwise
    \code{"neighborhood"}. None of the above is applied 
    to regions of equal length as there is no need for 
    that.
    
    Regarding the usage of \code{selector$id} field, this
    requires some careful usage, as if the ids present 
    there and the ids of the \code{genome} areas do not
    match, there will be no genomic regions left to 
    calculate coverage profiles on and the program will
    crash.
    
    Regarding the usage of the \code{preprocessParams}
    argument, the \code{normalize} field controls how
    the \code{GRanges} representing the reads extracted
    from BAM/BED files will be normalized. When 
    \code{"none"}, no normalization is applied and external
    normalization is assumed. When \code{"linear"}, all
    the library sizes are divided by the maximum one and
    a normalization factor is calculated for each sample.
    The coverage of this sample across the input genomic
    regions is then multiplied by this factor. When 
    \code{"downsample"}, all libraries are downsampled to
    the minimum library size among samples. When 
    \code{"sampleto"}, all libraries are downsampled to
    a fixed number of reads. The \code{"seed"} field is
    used to define a seed for the generation of the random
    indices of reads to be removed for downsampling. It
    is required for reproducible results. The 
    \code{sampleTo} field of \code{preprocessParams}
    tells recoup the fixed number of reads to downsample
    all libraries when 
    \code{preprocessParams$normalize="sampleTo"}. It 
    defaults to 1 million reads (\code{1e+6}). The 
    \code{spliceAction} field of \code{preprocessParams} is
    used to control the action to be taken in the presence
    of RNA-Seq spliced reads (implies \code{type="rnaseq"}).
    When \code{"keep"}, no action is performed regarding
    the spliced reads (represented as very long reads 
    spanning intronic regions in the \code{GRanges} object).
    When \code{"remove"}, these reads are excluded from 
    coverage calculations according to their length as 
    follows: firstly the length distribution of all reads
    lengths (using the \code{width} function for 
    \code{GRanges}) is calculated. Then the quantile 
    defined by the field \code{spliceRemoveQ} of 
    \code{preprocessParams} is calculated and reads above
    the length corresponding to this quantile are excluded.
    When \code{"splice"}, then splice junction information
    inferred from CIGAR strings (if) present in the BAM 
    files is used to splice the longer reads and calculate
    real coverages. This option is not available with
    BED files, however, BED files can contain pre-spliced
    reads using for example BEDTools for conversion.
    
    Regarding the \code{heatmapFactor} option of 
    \code{plotParams}, it controls the color scale of the
    heatmap as follows: the default value (1) causes the
    extremes of the heatmap colors to be linearly and
    equally distributed across the actual coverage profile
    values. If set smaller than 1, the the upper extreme of
    the coverage values (which by default maps to the upper
    color point) is multiplied by this factor and this new
    value is set as the upper color break (limit). This has
    the effect of decreasing the brightness of the heatmap
    as color is saturated before reaching the maximum
    coverage value. If set greater than 1, then the heatmap
    brightness is increased. Regarding the \code{correlation}
    option of \code{plotParams}, if \code{TRUE} then 
    \code{recoup} calculates average coverage values for
    each reference region (row-wise in the profile matrices)
    instead of the average coverage in each base of the
    reference regions (column-wise in the profile matrices).
    This is particularly useful for checking whether total
    genome profiles for some biological factor/condition
    correlate with each other. This potential correlation 
    is becoming even clearer when \code{orderBy$what} is 
    not \code{"none"}. Regarding the \code{corrScale} option
    of \code{plotParams}, it controls whether the average
    coverage curves over the set of reference genomic 
    regions (one average coverage vale per genomic region,
    note the difference with the profiles where the coverage
    is calculated over the genomic locations themselves)
    should be normalized to a 0-1 scale or not. This is 
    particularly useful when plotting data from different
    libraries (e.g. PolII and H3K27me1 occupancy over gene
    bodies) where other types of normalization (e.g. read
    downsampling cannot be applied). Regarding the 
    \code{corrSmoothPar} option of \code{plotParams}, it 
    controls the smoothing parameter for coverage 
    correlation curves. If \code{design} is present, 
    spline smoothing is applied 
    (\code{\link{smooth.spline}}) with \code{spar=0.5}
    else lowess smoothing is applied (\code{\link{lowess}})
    with \code{f=0.1}. \code{corrSmoothPar} controls the
    \code{spar} and \code{f} respectively.
    
    Regarding the usage of \code{saveParams} argument, 
    this is useful for several purposes: one is for re-using
    recoup without re-reading BAM/BED files. If the ranges 
    are present in the input object to recoup, they are 
    not re-calculated. If not stored, the memory/storage 
    usage is reduced but the object can be used only for
    simply replotting the profiles using 
    \code{\link{recoupProfile}} and/or 
    \code{\link{recoupHeatmap}} functions.
    
    As a note regarding parallel calculations, the number 
    of cores assigned to \code{recoup} depends both on the
    number of cores and the available RAM in your system. 
    The most RAM expensive part of \code{recoup} is 
    currently the construction of binned profile matrices.
    If you have a lot of cores (e.g. 16) but less than 
    128Gb of RAM for this number of cores, you should avoid
    using all cores, especially with large BAM files. Half of
    them would be more appropriate.
    
    Finally, the output list of \code{recoup} can be provided
    as input again to \code{recoup} with some input parameters
    changed. \code{recoup} will then automatically recognize 
    what has been changed and recalculate some, all or none
    of the genomic region profiles, depending on what input 
    parameters have changed. For example, if any of the 
    ordering options change (e.g. from no profile ordering to
    k-means clustering), then no recalculations are performed 
    and the process is very fast. If region binning is changed
    (\code{binParams$flankBinSize} or 
    \code{binParams$regionBinSize}), then only profile matrices
    are recalculated and coverages are maintained. If any of the
    \code{preprocessParams} changes, this causes all object 
    including the short reads to be reimported and profiles
    recalculated from the beginning.
}
\examples{
# Load some sample data
data("recoup_test_data",package="recoup")

# Note: the figures that will be produced will not look 
# realistic and will be "bumpy". This is because package
# size limitations posed by Bioconductor guidelines do not
# allow for a full test dataset. Have a look at the
# vignette on how to test with more realistic data.

# TSS high resolution profile with no design
test.tss <- recoup(
    test.input,
    design=NULL,
    region="tss",
    type="chipseq",
    genome=test.genome,
    flank=c(2000,2000),
    selector=NULL,
    rc=0.5
)
\donttest{
# Genebody low resolution profile with 2-factor design, 
# wide genebody and more narrow flanking
test.gb <- recoup(
    test.input,
    design=test.design,
    region="genebody",
    type="chipseq",
    genome=test.genome,
    flank=c(2000,2000),
    binParams=list(flankBinSize=50,regionBinSize=150),
    orderBy=list(what="hc1"),
    selector=NULL,
    rc=0.5
)
}}
\author{
    Panagiotis Moulos
}
