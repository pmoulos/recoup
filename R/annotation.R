buildAnnotationStore <- function(organisms,sources,
    home=file.path(path.expand("~"),".recoup"),forceDownload=TRUE,rc=NULL) {

    if (missing(organisms))
        organisms <- c("hg18","hg19","hg38","mm9","mm10","rn5","dm3","danrer7",
            "pantro4","susscr3")
    if (missing(sources))
        sources <- c("ensembl","ucsc","refseq")
    
    checkTextArgs("organisms",organisms,c("hg18","hg19","hg38","mm9","mm10",
        "rn5","dm3","danrer7","pantro4","susscr3"),multiarg=FALSE)
    checkTextArgs("sources",sources,c("ensembl","ucsc","refseq"),multiarg=FALSE)
    
    if (!requireNamespace("GenomeInfoDb"))
        stop("R package GenomeInfoDb is required to construct annotation ",
            "stores!")
    
    for (s in sources) {
        for (o in organisms) {
            # Retrieving genome info
            message("Retrieving genome information for ",o)
            sf <- GenomeInfoDb::fetchExtendedChromInfoFromUCSC(
                getUcscOrganism(o))
            rownames(sf) <- as.character(sf[,1])
            sf <- sf[getValidChrs(o),]
            sf <- Seqinfo(seqnames=sf[,1],seqlengths=sf[,2],
                isCircular=sf[,3],genome=getUcscOrganism(o))
        
            # Retrieve gene annotations
            store.path <- file.path(home,s,o)
            if (!dir.exists(store.path))
                dir.create(store.path,recursive=TRUE,mode="0755")
            if (file.exists(file.path(store.path,"gene.rda")) && !forceDownload)
                message("Gene annotation for ",o," from ",s," has already ",
                    "created and will be skipped. If you wish to recreate it ",
                    "choose forceDownload = TRUE.")
            else {
                message("Retrieving gene annotation for ",o," from ",s)
                ann <- getAnnotation(o,"gene",refdb=s,rc=rc)
                gene <- makeGRangesFromDataFrame(
                    df=ann,
                    seqinfo=sf,
                    keep.extra.columns=TRUE,
                    seqnames.field="chromosome"
                )
                save(gene,file=file.path(store.path,"gene.rda"),compress=TRUE)
            }
            
            # Retrieve exon annotations
            store.path <- file.path(home,s,o)
            if (!dir.exists(store.path))
                dir.create(store.path,recursive=TRUE,mode="0755")
            
            if (file.exists(file.path(store.path,"exon.rda")) && !forceDownload)
                message("Exon annotation for ",o," from ",s," has already ",
                    "created and will be skipped. If you wish to recreate it ",
                    "choose forceDownload = TRUE.")
            else {
                message("Retrieving exon annotation for ",o," from ",s)
                ann <- getAnnotation(o,"exon",refdb=s,rc=rc)
                ann.gr <- makeGRangesFromDataFrame(
                    df=ann,
                    seqinfo=sf,
                    keep.extra.columns=TRUE,
                    seqnames.field="chromosome"
                )
                exon <- split(ann.gr,ann.gr$gene_id)                
                save(exon,file=file.path(store.path,"exon.rda"),compress=TRUE)
            }
            
            # Then summarize the exons and write again with type sum_exon
            if (file.exists(file.path(store.path,"summarized_exon.rda")) 
                && !forceDownload)
                message("Summarized exon annotation for ",o," from ",s,
                    " has already been created and will be skipped. If you ",
                    "wish to recreate it choose forceDownload = TRUE.")
            else {
                message("Merging exons for ",o," from ",s)
                ann.list <- reduceExons(ann.gr,rc=rc)
                ann.gr <- ann.list$model
                names(ann.gr) <- as.character(ann.gr$exon_id)
                sexon <- split(ann.gr,ann.gr$gene_id)
                activeLength <- ann.list$length
                names(activeLength) <- names(sexon)
                #names(ann.gr) <- as.character(ann.gr$exon_id)
                #sexon <- split(ann.gr,ann.gr$gene_id)
                save(sexon,activeLength,
                    file=file.path(store.path,"summarized_exon.rda"),
                    compress=TRUE)
            }
        }
    }
}

reduceExons <- function(gr,rc=NULL) {
    gene <- unique(as.character(gr$gene_id))
    if (!is.null(gr$gene_name))
        gn <- gr$gene_name
    else
        gn <- NULL
    if (!is.null(gr$biotype))
        bt <- gr$biotype   
    else
        bt <- NULL
    red.list <- cmclapply(gene,function(x,a,g,b) {
        tmp <- a[a$gene_id==x]
        if (!is.null(g))
            gena <- as.character(tmp$gene_name[1])
        if (!is.null(b))
            btty <- as.character(tmp$biotype[1])
        merged <- reduce(tmp)
        n <- length(merged)
        meta <- DataFrame(
            exon_id=paste(x,"MEX",1:n,sep="_"),
            gene_id=rep(x,n)
        )
        if (!is.null(g))
            meta$gene_name <- rep(gena,n)
        if (!is.null(b))
            meta$biotype <- rep(btty,n)
        mcols(merged) <- meta
        return(merged)
    },gr,gn,bt,rc=rc)
    len <- unlist(cmclapply(red.list,function(x) {
        return(sum(width(x)))
    },rc=rc))
    names(len) <- names(red.list)
    #return(list(model=GRangesList(red.list),length=len))
    #return(do.call("c",red.list))
    return(list(model=do.call("c",red.list),length=len))
}

getAnnotation <- function(org,type,refdb="ensembl",rc=NULL) {
    org <- tolower(org)
    switch(refdb,
        ensembl = { return(getEnsemblAnnotation(org,type)) },
        ucsc = { return(getUcscAnnotation(org,type,refdb,rc=rc)) },
        refseq = { return(getUcscAnnotation(org,type,refdb,rc=rc)) }
    )
}

getEnsemblAnnotation <- function(org,type) {
    dat <- "ENSEMBL_MART_ENSEMBL"
    mart <- tryCatch({
            useMart(biomart=dat,host=getHost(org),dataset=getDataset(org))
        },
        error=function(e) {
            useMart(biomart=dat,host=getAltHost(org),
            dataset=getDataset(org))
        },
        finally={})
    chrs.exp <- paste(getValidChrs(org),collapse="|")
    if (type=="gene") {
        bm <- getBM(attributes=getGeneAttributes(org),mart=mart)
        ann <- data.frame(
            chromosome=paste("chr",bm$chromosome_name,sep=""),
            start=bm$start_position,
            end=bm$end_position,
            gene_id=bm$ensembl_gene_id,
            gc_content=bm$percentage_gc_content,
            strand=ifelse(bm$strand==1,"+","-"),
            gene_name=if (org %in% c("hg18","mm9","tair10")) bm$external_gene_id 
                else bm$external_gene_name,
            biotype=bm$gene_biotype
        )
        rownames(ann) <- ann$gene_id
    }
    else if (type=="exon") {
        bm <- getBM(attributes=getExonAttributes(org),mart=mart)
        if (org == "hg19") {
            message("  Bypassing problem with hg19 Ensembl combined gene-exon ",
                "annotation... Will take slightly longer...")
            bmg <- getBM(attributes=getGeneAttributes(org),mart=mart)
            gene_name <- bmg$external_gene_name
            names(gene_name) <- bmg$ensembl_gene_id
            ann <- data.frame(
                chromosome=paste("chr",bm$chromosome_name,sep=""),
                start=bm$exon_chrom_start,
                end=bm$exon_chrom_end,
                exon_id=bm$ensembl_exon_id,
                gene_id=bm$ensembl_gene_id,
                strand=ifelse(bm$strand==1,"+","-"),
                gene_name=gene_name[bm$ensembl_gene_id],
                biotype=bm$gene_biotype
            )
            rownames(ann) <- ann$exon_id
        }
        else
            ann <- data.frame(
                chromosome=paste("chr",bm$chromosome_name,sep=""),
                start=bm$exon_chrom_start,
                end=bm$exon_chrom_end,
                exon_id=bm$ensembl_exon_id,
                gene_id=bm$ensembl_gene_id,
                strand=ifelse(bm$strand==1,"+","-"),
                gene_name=if (org %in% c("hg18","mm9","tair10")) 
                    bm$external_gene_id else bm$external_gene_name,
                biotype=bm$gene_biotype
            )
            rownames(ann) <- ann$exon_id
    }
    ann <- ann[order(ann$chromosome,ann$start),]
    ann <- ann[grep(chrs.exp,ann$chromosome),]
    ann$chromosome <- as.character(ann$chromosome)
    return(ann)
}

getUcscAnnotation <- function(org,type,refdb="ucsc",rc=NULL) {
    if (!requireNamespace("RMySQL")) {
        rmysql.present <- FALSE
        warning("R package RMySQL is not present! Annotation will be ",
            "retrieved by downloading temporary files from UCSC and the usage
            of a temporary SQLite database...",immediate.=TRUE)
    }
    else
        rmysql.present <- TRUE
    if (!requireNamespace("RSQLite"))
        stop("R package RSQLite is required to use annotation from UCSC!")

    valid.chrs <- getValidChrs(org)
    chrs.exp <- paste("^",paste(valid.chrs,collapse="$|^"),"$",sep="")

    db.org <- getUcscOrganism(org)
    if (rmysql.present) {
        db.creds <- getUcscCredentials()
        drv <- dbDriver("MySQL")
        con <- dbConnect(drv,user=db.creds[2],password=NULL,dbname=db.org,
            host=db.creds[1])
        query <- getUcscQuery(org,type,refdb)
        raw.ann <- dbGetQuery(con,query)
        dbDisconnect(con)
    }
    else {
        # This should return the same data frame as the db query
        tmp.sqlite <- getUcscDbl(org,type,refdb)
        drv <- dbDriver("SQLite")
        con <- dbConnect(drv,dbname=tmp.sqlite)
        query <- getUcscQuery(org,type,refdb)
        raw.ann <- dbGetQuery(con,query)
        dbDisconnect(con)
    }
    if (type=="gene") {
        ann <- raw.ann
        ann <- ann[grep(chrs.exp,ann$chromosome,perl=TRUE),]
        ann$chromosome <- as.character(ann$chromosome)
        rownames(ann) <- ann$gene_id
    }
    else if (type=="exon") {
        raw.ann <- raw.ann[grep(chrs.exp,raw.ann$chromosome,perl=TRUE),]
        ex.list <- cmclapply(as.list(1:nrow(raw.ann)),function(x,d,s) {
            r <- d[x,]
            starts <- as.numeric(strsplit(r[,"start"],",")[[1]])
            ends <- as.numeric(strsplit(r[,"end"],",")[[1]])
            nexons <- length(starts)
            ret <- data.frame(
                rep(r[,"chromosome"],nexons),
                starts,ends,
                paste(r[,"exon_id"],"_e",1:nexons,sep=""),
                rep(r[,"strand"],nexons),
                rep(r[,"gene_id"],nexons),
                rep(r[,"gene_name"],nexons),
                rep(r[,"biotype"],nexons)
            )
            names(ret) <- names(r)
            rownames(ret) <- ret$exon_id
            ret <- makeGRangesFromDataFrame(
                df=ret,
                keep.extra.columns=TRUE,
                seqnames.field="chromosome",
                seqinfo=s
            )
            return(ret)
        },raw.ann,valid.chrs,rc=rc)
        tmp.ann <- do.call("c",ex.list)
        ann <- data.frame(
            chromosome=as.character(seqnames(tmp.ann)),
            start=start(tmp.ann),
            end=end(tmp.ann),
            exon_id=as.character(tmp.ann$exon_id),
            gene_id=as.character(tmp.ann$gene_id),
            strand=as.character(strand(tmp.ann)),
            gene_name=as.character(tmp.ann$gene_name),
            biotype=tmp.ann$biotype
        )
        rownames(ann) <- ann$exon_id
    }
    
    gc.content <- getGcContent(ann,org)
    ann$gc_content <- gc.content
    ann <- ann[order(ann$chromosome,ann$start),]
    return(ann)
}

getGcContent <- function(ann,org) {
    if (missing(ann))
        stop("A valid annotation data frame must be provided in order to ",
            "retrieve GC-content.")
    org <- tolower(org[1])
    checkTextArgs("org",org,c("hg18","hg19","hg38","mm9","mm10","rn5","dm3",
        "danrer7","pantro4","susscr3"),multiarg=FALSE)
    # Convert annotation to GRanges
    message("Converting annotation to GenomicRanges object...")
    if (packageVersion("GenomicRanges")<1.14)
        ann.gr <- GRanges(
            seqnames=Rle(ann[,1]),
            ranges=IRanges(start=ann[,2],end=ann[,3]),
            strand=Rle(ann[,6]),
            name=as.character(ann[,4])
        )
    else
        ann.gr <- makeGRangesFromDataFrame(
            df=ann,
            keep.extra.columns=TRUE,
            seqnames.field="chromosome"
        )
    bsg <- loadBsGenome(org)
    message("Getting DNA sequences...")
    seqs <- getSeq(bsg,names=ann.gr)
    message("Getting GC content...")
    freq.matrix <- alphabetFrequency(seqs,as.prob=TRUE,baseOnly=TRUE)
    gc.content <- apply(freq.matrix,1,function(x) round(100*sum(x[2:3]),
        digits=2))
    names(gc.content) <- as.character(ann[,4])
    return(gc.content)
}

getUcscOrganism <- function(org) {
    switch(org,
        hg18 = { return("hg18") },
        hg19 = { return("hg19") },
        hg38 = { return("hg38") },
        mm9 = { return("mm9") },
        mm10 = { return("mm10") },
        rn5 = { return("rn5") },
        dm3 = { return("dm3") },
        danrer7 = { return("danRer7") },
        pantro4 = { return("panTro4") },
        susscr3 = { return("susScr3") }
    )
}

getBsOrganism <- function(org) {
    switch(org,
        hg18 = {
            return("BSgenome.Hsapiens.UCSC.hg18")
        },
        hg19 = {
            return("BSgenome.Hsapiens.UCSC.hg19")
        },
        hg38 = {
            return("BSgenome.Hsapiens.UCSC.hg38")
        },
        mm9 = {
            return("BSgenome.Mmusculus.UCSC.mm9")
        },
        mm10 = {
            return("BSgenome.Mmusculus.UCSC.mm10")
        },
        rn5 = {
            return("BSgenome.Rnorvegicus.UCSC.rn5")
        },
        dm3 = {
            return("BSgenome.Dmelanogaster.UCSC.dm3")
        },
        danrer7 = {
            return("BSgenome.Drerio.UCSC.danRer7")
        },
        pantro4 = {
            stop("panTro4 is not yet supported by BSgenome! Please use ",
                "Ensembl as annoation source.")
        },
        susscr3 = {
            return("BSgenome.Sscrofa.UCSC.susScr3")
        }
    )
}

loadBsGenome <- function(org) {
    if (!requireNamespace("BiocInstaller"))
        stop("The Bioconductor package BiocInstaller is required to ",
            "proceed!")
    if (!requireNamespace("BSgenome"))
        stop("The Bioconductor package BSgenome is required to ",
            "proceed!")
    bs.org <- getBsOrganism(org)
    if (bs.org %in% installed.genomes())
        bs.obj <- getBSgenome(getUcscOrganism(org))
    else {
        biocLite(bs.org)
        bs.obj <- getBSgenome(getUcscOrganism(org))
    }
    return(bs.obj)
}

getBiotypes <- function(a) {
    return(as.character(unique(a$biotype)))
}

getHost <- function(org) {
    switch(org,
        hg18 = { return("may2009.archive.ensembl.org") },
        hg19 = { return("grch37.ensembl.org") },
        hg38 = { return("www.ensembl.org") },
        mm9 = { return("may2012.archive.ensembl.org") },
        mm10 = { return("www.ensembl.org") },
        rn5 = { return("www.ensembl.org") },
        dm3 = { return("www.ensembl.org") },
        danrer7 = { return("www.ensembl.org") },
        pantro4 = { return("www.ensembl.org") },
        susscr3 = { return("www.ensembl.org") }
    )
}

getAltHost <- function(org) {
    switch(org,
        hg18 = { return("may2009.archive.ensembl.org") },
        hg19 = { return("grch37.ensembl.org") },
        hg38 = { return("uswest.ensembl.org") },
        mm9 = { return("may2012.archive.ensembl.org") },
        mm10 = { return("uswest.ensembl.org") },
        rn5 = { return("uswest.ensembl.org") },
        dm3 = { return("uswest.ensembl.org") },
        danrer7 = { return("uswest.ensembl.org") },
        pantro4 = { return("uswest.ensembl.org") },
        susscr3 = { return("uswest.ensembl.org") }
    )
}

getDataset <- function(org) {
    switch(org,
        hg18 = { return("hsapiens_gene_ensembl") },
        hg19 = { return("hsapiens_gene_ensembl") },
        hg38 = { return("hsapiens_gene_ensembl") },
        mm9 = { return("mmusculus_gene_ensembl") },
        mm10 = { return("mmusculus_gene_ensembl") },
        rn5 = { return("rnorvegicus_gene_ensembl") },
        dm3 = { return("dmelanogaster_gene_ensembl") },
        danrer7 = { return("drerio_gene_ensembl") },
        pantro4 = { return("ptroglodytes_gene_ensembl") },
        susscr3 = { return("sscrofa_gene_ensembl") }
    )
}

getValidChrs <- function(org)
{
    switch(org,
        hg18 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr20","chr21","chr22","chr3",
                "chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY"
            ))
        },
        hg19 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr20","chr21","chr22","chr3",
                "chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY"
            ))
        },
        hg38 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr20","chr21","chr22","chr3",
                "chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY"
            ))
        },
        mm9 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr3","chr4","chr5","chr6",
                "chr7","chr8","chr9","chrX","chrY"
            ))
        },
        mm10 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr3","chr4","chr5","chr6",
                "chr7","chr8","chr9","chrX","chrY"
            ))
        },
        rn5 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr3","chr4","chr5","chr6",
                "chr7","chr8","chr9","chrX"
            ))
        },
        dm3 = {
            return(c(
                "chr2L","chr2LHet","chr2R","chr2RHet","chr3L","chr3LHet",
                "chr3R","chr3RHet","chr4","chrU","chrUextra","chrX","chrXHet",
                "chrYHet"
            ))
        },
        danrer7 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr20","chr21","chr22","chr23",
                "chr24","chr25","chr3","chr4","chr5","chr6","chr7","chr8","chr9"
            ))
        },
        pantro4 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr20","chr21","chr22","chr2A","chr2B",
                "chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY"
            ))
        },
        susscr3 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr2","chr3","chr4","chr5","chr6","chr7",
                "chr8","chr9","chrX","chrY"
            ))
        }
    )
}

getGeneAttributes <- function(org) {
    if (org %in% c("hg18","mm9"))
        return(c(
            "chromosome_name",
            "start_position",
            "end_position",
            "ensembl_gene_id",
            "percentage_gc_content",
            "strand",
            "external_gene_id",
            "gene_biotype"
        ))
    else
        return(c(
            "chromosome_name",
            "start_position",
            "end_position",
            "ensembl_gene_id",
            "percentage_gc_content",
            "strand",
            "external_gene_name",
            "gene_biotype"
        ))
}

getExonAttributes <- function(org) {
    if (org %in% c("hg18","mm9"))
        return(c(
            "chromosome_name",
            "exon_chrom_start",
            "exon_chrom_end",
            "ensembl_exon_id",
            "strand",
            "ensembl_gene_id",
            "external_gene_id",
            "gene_biotype"
        ))
    else if (org == "hg19")
        return(c(
            "chromosome_name",
            "exon_chrom_start",
            "exon_chrom_end",
            "ensembl_exon_id",
            "strand",
            "ensembl_gene_id",
            "gene_biotype"
        ))
    else
        return(c(
            "chromosome_name",
            "exon_chrom_start",
            "exon_chrom_end",
            "ensembl_exon_id",
            "strand",
            "ensembl_gene_id",
            "external_gene_name",
            "gene_biotype"
        ))
}

getUcscDbl <- function(org,type,refdb="ucsc") {
    type <- tolower(type[1])
    org <- tolower(org[1])
    refdb <- tolower(refdb[1])
    checkTextArgs("type",type,c("gene","exon"))
    checkTextArgs("org",org,c("hg18","hg19","hg38","mm9","mm10","rn5","dm3",
        "danrer7","pantro4","susscr3"),multiarg=FALSE)
    checkTextArgs("refdb",refdb,c("ucsc","refseq"))
    
    if (!requireNamespace("RSQLite"))
        stop("R package RSQLite is required to use annotation from UCSC!")

    http.base <- paste("http://hgdownload.soe.ucsc.edu/goldenPath/",
        getUcscOrganism(org),"/database/",sep="")
    table.defs <- getUcscTabledef(org,type,refdb,"fields")
    file.list <- vector("list",length(table.defs))
    names(file.list) <- names(table.defs)
    for (n in names(file.list))
        file.list[[n]] <- paste(http.base,n,".txt.gz",sep="")
        
    # Fill the fields for each table
    drv <- dbDriver("SQLite")
    db.tmp <- tempfile()
    con <- dbConnect(drv,dbname=db.tmp)
    message("  Retrieving tables for temporary SQLite ",refdb," ",org," ",type,
        " subset database")
    for (n in names(file.list)) {
        message("    Retrieving table ",n)
        download.file(file.list[[n]],file.path(tempdir(),
            paste(n,".txt.gz",sep="")),quiet=TRUE)
        if (.Platform$OS.type == "unix")
            system(paste("gzip -df",file.path(tempdir(),
                paste(n,".txt.gz",sep=""))))
        else
            unzip(file.path(tempdir(),paste(n,".txt.gz",sep="")))
        sql.df <- read.delim(file.path(tempdir(),paste(n,".txt",sep="")),
            row.names=NULL,header=FALSE,strip.white=TRUE)
        names(sql.df) <- table.defs[[n]]
        dbWriteTable(con,n,sql.df,row.names=FALSE)
    }
    dbDisconnect(con)
    return(db.tmp)
}

getUcscTabledef <- function(org,type,refdb="ucsc",what="queries") {
    type <- tolower(type[1])
    org <- tolower(org[1])
    refdb <- tolower(refdb[1])
    what <- tolower(what[1])
    checkTextArgs("type",type,c("gene","exon"))
    checkTextArgs("org",org,c("hg18","hg19","hg38","mm9","mm10","rn5","dm3",
        "danrer7","pantro4","susscr3"),multiarg=FALSE)
    checkTextArgs("refdb",refdb,c("ucsc","refseq"))
    checkTextArgs("what",what,c("queries","fields"))
    switch(type,
        gene = {
            switch(refdb,
                ucsc = {
                    switch(org,
                        hg18 = {
                            return(list(
                                knownCanonical=
                                    getUcscTblTpl("knownCanonical",what),
                                knownGene=getUcscTblTpl("knownGene",what),
                                knownToRefSeq=
                                    getUcscTblTpl("knownToRefSeq",what),
                                refFlat=getUcscTblTpl("refFlat",what)
                            ))
                        },
                        hg19 = {
                            return(list(
                                knownCanonical=
                                    getUcscTblTpl("knownCanonical",what),
                                knownGene=getUcscTblTpl("knownGene",what),
                                knownToRefSeq=
                                    getUcscTblTpl("knownToRefSeq",what),
                                knownToEnsembl=
                                    getUcscTblTpl("knownToEnsembl",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what),
                                refFlat=getUcscTblTpl("refFlat",what)
                            ))
                        },
                        hg38 = {
                            return(list(
                                knownCanonical=
                                    getUcscTblTpl("knownCanonical",what),
                                knownGene=getUcscTblTpl("knownGene",what),
                                knownToRefSeq=
                                    getUcscTblTpl("knownToRefSeq",what),
                                refFlat=getUcscTblTpl("refFlat",what)
                            ))
                        },
                        mm9 = {
                            return(list(
                                knownCanonical=
                                    getUcscTblTpl("knownCanonical",what),
                                knownGene=getUcscTblTpl("knownGene",what),
                                knownToRefSeq=
                                    getUcscTblTpl("knownToRefSeq",what),
                                knownToEnsembl=
                                    getUcscTblTpl("knownToEnsembl",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what),
                                refFlat=getUcscTblTpl("refFlat",what)
                            ))
                        },
                        mm10 = {
                            return(list(
                                knownCanonical=
                                    getUcscTblTpl("knownCanonical",what),
                                knownGene=getUcscTblTpl("knownGene",what),
                                knownToRefSeq=
                                    getUcscTblTpl("knownToRefSeq",what),
                                knownToEnsembl=
                                    getUcscTblTpl("knownToEnsembl",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what),
                                refFlat=getUcscTblTpl("refFlat",what)
                            ))
                        },
                        rn5 = {
                            return(list(
                                mgcGenes=getUcscTblTpl("mgcGenes",what),
                                ensemblToGeneName=
                                    getUcscTblTpl("ensemblToGeneName",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what)
                            ))
                        },
                        dm3 = {
                            return(list(
                                flyBaseCanonical=
                                    getUcscTblTpl("flyBaseCanonical",what),
                                flyBaseGene=
                                    getUcscTblTpl("flyBaseGene",what),
                                flyBaseToRefSeq=
                                    getUcscTblTpl("flyBaseToRefSeq",what),
                                ensemblToGeneName=
                                    getUcscTblTpl("ensemblToGeneName",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what)
                            ))
                        },
                        danrer7 = {
                            return(list(
                                mgcGenes=getUcscTblTpl("mgcGenes",what),
                                ensemblToGeneName=
                                    getUcscTblTpl("ensemblToGeneName",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what)
                            ))
                        },
                        pantro4 = {
                            warning("No UCSC Genome annotation for Pan ",
                                "troglodytes! Will use RefSeq instead...",
                                immediate.=TRUE)
                            return(list(
                                refFlat=getUcscTblTpl("refFlat",what),
                                ensemblToGeneName=
                                    getUcscTblTpl("ensemblToGeneName",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what)
                            ))
                        },
                        susscr3 = {
                            warning("No UCSC Genome annotation for Sus ",
                                "scrofa! Will use RefSeq instead...",
                                immediate.=TRUE)
                            return(list(
                                refFlat=getUcscTblTpl("refFlat",what),
                                ensemblToGeneName=
                                    getUcscTblTpl("ensemblToGeneName",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what)
                            ))
                        }
                    )
                },
                refseq = {
                    switch(org,
                        hg18 = {
                            return(list(
                                refFlat=getUcscTblTpl("refFlat",what),
                                knownToRefSeq=
                                    getUcscTblTpl("knownToRefSeq",what),
                                knownCanonical=
                                    getUcscTblTpl("knownCanonical",what)
                            ))
                        },
                        hg19 = {
                            return(list(
                                refFlat=getUcscTblTpl("refFlat",what),
                                knownToRefSeq=
                                    getUcscTblTpl("knownToRefSeq",what),
                                knownCanonical=
                                    getUcscTblTpl("knownCanonical",what),
                                knownToEnsembl=
                                    getUcscTblTpl("knownToEnsembl",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what)
                            ))
                        },
                        hg38 = {
                            return(list(
                                refFlat=getUcscTblTpl("refFlat",what),
                                knownToRefSeq=
                                    getUcscTblTpl("knownToRefSeq",what),
                                knownCanonical=
                                    getUcscTblTpl("knownCanonical",what)
                            ))
                        },
                        mm9 = {
                            return(list(
                                refFlat=getUcscTblTpl("refFlat",what),
                                knownToRefSeq=
                                    getUcscTblTpl("knownToRefSeq",what),
                                knownCanonical=
                                    getUcscTblTpl("knownCanonical",what),
                                knownToEnsembl=
                                    getUcscTblTpl("knownToEnsembl",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what)
                            ))
                        },
                        mm10 = {
                            return(list(
                                refFlat=getUcscTblTpl("refFlat",what),
                                knownToRefSeq=
                                    getUcscTblTpl("knownToRefSeq",what),
                                knownCanonical=
                                    getUcscTblTpl("knownCanonical",what),
                                knownToEnsembl=
                                    getUcscTblTpl("knownToEnsembl",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what)
                            ))
                        },
                        rn5 = {
                            return(list(
                                refFlat=getUcscTblTpl("refFlat",what),
                                ensemblToGeneName=
                                    getUcscTblTpl("ensemblToGeneName",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what)
                            ))
                        },
                        dm3 = {
                            return(list(
                                refFlat=getUcscTblTpl("refFlat",what),
                                ensemblToGeneName=
                                    getUcscTblTpl("ensemblToGeneName",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what)
                            ))
                        },
                        danrer7 = {
                            return(list(
                                refFlat=getUcscTblTpl("refFlat",what),
                                ensemblToGeneName=
                                    getUcscTblTpl("ensemblToGeneName",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what)
                            ))
                        },
                        pantro4 = {
                            return(list(
                                refFlat=getUcscTblTpl("refFlat",what),
                                ensemblToGeneName=
                                    getUcscTblTpl("ensemblToGeneName",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what)
                            ))
                        },
                        susscr3 = {
                            return(list(
                                refFlat=getUcscTblTpl("refFlat",what),
                                ensemblToGeneName=
                                    getUcscTblTpl("ensemblToGeneName",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what)
                            ))
                        }
                    )
                }
            )
        },
        exon = {
            switch(refdb,
                ucsc = {
                    switch(org,
                        hg18 = {
                            return(list(
                                knownCanonical=
                                    getUcscTblTpl("knownCanonical",what),
                                knownGene=getUcscTblTpl("knownGene",what),
                                knownToRefSeq=
                                    getUcscTblTpl("knownToRefSeq",what),
                                refFlat=getUcscTblTpl("refFlat",what)
                            ))
                        },
                        hg19 = {
                            return(list(
                                knownCanonical=
                                    getUcscTblTpl("knownCanonical",what),
                                knownGene=getUcscTblTpl("knownGene",what),
                                knownToRefSeq=
                                    getUcscTblTpl("knownToRefSeq",what),
                                knownToEnsembl=
                                    getUcscTblTpl("knownToEnsembl",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what),
                                refFlat=getUcscTblTpl("refFlat",what)
                            ))
                        },
                        hg38 = {
                            return(list(
                                knownCanonical=
                                    getUcscTblTpl("knownCanonical",what),
                                knownGene=getUcscTblTpl("knownGene",what),
                                knownToRefSeq=
                                    getUcscTblTpl("knownToRefSeq",what),
                                refFlat=getUcscTblTpl("refFlat",what)
                            ))
                        },
                        mm9 = {
                            return(list(
                                knownCanonical=
                                    getUcscTblTpl("knownCanonical",what),
                                knownGene=getUcscTblTpl("knownGene",what),
                                knownToRefSeq=
                                    getUcscTblTpl("knownToRefSeq",what),
                                knownToEnsembl=
                                    getUcscTblTpl("knownToEnsembl",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what),
                                refFlat=getUcscTblTpl("refFlat",what)
                            ))
                        },
                        mm10 = {
                            return(list(
                                knownCanonical=
                                    getUcscTblTpl("knownCanonical",what),
                                knownGene=getUcscTblTpl("knownGene",what),
                                knownToRefSeq=
                                    getUcscTblTpl("knownToRefSeq",what),
                                knownToEnsembl=
                                    getUcscTblTpl("knownToEnsembl",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what),
                                refFlat=getUcscTblTpl("refFlat",what)
                            ))
                        },
                        rn5 = {
                            return(list(
                                mgcGenes=getUcscTblTpl("mgcGenes",what),
                                ensemblToGeneName=
                                    getUcscTblTpl("ensemblToGeneName",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what)
                            ))
                        },
                        dm3 = {
                            return(list(
                                flyBaseCanonical=
                                    getUcscTblTpl("flyBaseCanonical",what),
                                flyBaseGene=
                                    getUcscTblTpl("flyBaseGene",what),
                                flyBaseToRefSeq=
                                    getUcscTblTpl("flyBaseToRefSeq",what),
                                ensemblToGeneName=
                                    getUcscTblTpl("ensemblToGeneName",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what)
                            ))
                        },
                        danrer7 = {
                            return(list(
                                mgcGenes=getUcscTblTpl("mgcGenes",what),
                                ensemblToGeneName=
                                    getUcscTblTpl("ensemblToGeneName",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what)
                            ))
                        },
                        pantro4 = {
                            warning("No UCSC Genome annotation for Pan ",
                                "troglodytes! Will use RefSeq instead...",
                                immediate.=TRUE)
                            return(list(
                                refFlat=getUcscTblTpl("refFlat",what),
                                ensemblToGeneName=
                                    getUcscTblTpl("ensemblToGeneName",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what)
                            ))
                        },
                        susscr3 = {
                            warning("No UCSC Genome annotation for Sus ",
                                "scrofa! Will use RefSeq instead...",
                                immediate.=TRUE)
                            return(list(
                                refFlat=getUcscTblTpl("refFlat",what),
                                ensemblToGeneName=
                                    getUcscTblTpl("ensemblToGeneName",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what)
                            ))
                        }
                    )
                },
                refseq = {
                    switch(org,
                        hg18 = {
                            return(list(
                                refFlat=getUcscTblTpl("refFlat",what),
                                knownToRefSeq=
                                    getUcscTblTpl("knownToRefSeq",what),
                                knownCanonical=
                                    getUcscTblTpl("knownCanonical",what)
                            ))
                        },
                        hg19 = {
                            return(list(
                                refFlat=getUcscTblTpl("refFlat",what),
                                knownToRefSeq=
                                    getUcscTblTpl("knownToRefSeq",what),
                                knownCanonical=
                                    getUcscTblTpl("knownCanonical",what),
                                knownToEnsembl=
                                    getUcscTblTpl("knownToEnsembl",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what)
                            ))
                        },
                        hg38 = {
                            return(list(
                                refFlat=getUcscTblTpl("refFlat",what),
                                knownToRefSeq=
                                    getUcscTblTpl("knownToRefSeq",what),
                                knownCanonical=
                                    getUcscTblTpl("knownCanonical",what)
                            ))
                        },
                        mm9 = {
                            return(list(
                                refFlat=getUcscTblTpl("refFlat",what),
                                knownToRefSeq=
                                    getUcscTblTpl("knownToRefSeq",what),
                                knownCanonical=
                                    getUcscTblTpl("knownCanonical",what),
                                knownToEnsembl=
                                    getUcscTblTpl("knownToEnsembl",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what)
                            ))
                        },
                        mm10 = {
                            return(list(
                                refFlat=getUcscTblTpl("refFlat",what),
                                knownToRefSeq=
                                    getUcscTblTpl("knownToRefSeq",what),
                                knownCanonical=
                                    getUcscTblTpl("knownCanonical",what),
                                knownToEnsembl=
                                    getUcscTblTpl("knownToEnsembl",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what)
                            ))
                        },
                        rn5 = {
                            return(list(
                                refFlat=getUcscTblTpl("refFlat",what),
                                ensemblToGeneName=
                                    getUcscTblTpl("ensemblToGeneName",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what)
                            ))
                        },
                        dm3 = {
                            return(list(
                                refFlat=getUcscTblTpl("refFlat",what),
                                ensemblToGeneName=
                                    getUcscTblTpl("ensemblToGeneName",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what)
                            ))
                        },
                        danrer7 = {
                            return(list(
                                refFlat=getUcscTblTpl("refFlat",what),
                                ensemblToGeneName=
                                    getUcscTblTpl("ensemblToGeneName",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what)
                            ))
                        },
                        pantro4 = {
                            return(list(
                                refFlat=getUcscTblTpl("refFlat",what),
                                ensemblToGeneName=
                                    getUcscTblTpl("ensemblToGeneName",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what)
                            ))
                        },
                        susscr3 = {
                            return(list(
                                refFlat=getUcscTblTpl("refFlat",what),
                                ensemblToGeneName=
                                    getUcscTblTpl("ensemblToGeneName",what),
                                ensemblSource=
                                    getUcscTblTpl("ensemblSource",what)
                            ))
                        }
                    )
                }
            )
        }
    )
}

getUcscTblTpl <- function(tab,what="queries") {
    if (what=="queries") {
        switch(tab,
            knownCanonical = {
                return(paste(
                    "CREATE TABLE",
                    "`knownCanonical` (",
                    "`chrom` TEXT NOT NULL DEFAULT '',",
                    "`chromStart` INTEGER NOT NULL DEFAULT '0',",
                    "`chromEnd` INTEGER NOT NULL DEFAULT '0',",
                    "`clusterId` INTEGER NOT NULL DEFAULT '0',",
                    "`transcript` TEXT NOT NULL DEFAULT '',",
                    "`protein` TEXT NOT NULL DEFAULT ''",
                    ")",collapse=" "
                ))
            },
            knownGene = {
                return(paste(
                    "CREATE TABLE",
                    "`knownGene` (",
                    "`name` TEXT NOT NULL DEFAULT '',",
                    "`chrom` TEXT NOT NULL DEFAULT '',",
                    "`strand` TEXT NOT NULL DEFAULT '',",
                    "`txStart` UNSIGNED INTEGER NOT NULL DEFAULT '0',",
                    "`txEnd` UNSIGNED INTEGER NOT NULL DEFAULT '0',",
                    "`cdsStart` UNSIGNED INTEGER NOT NULL DEFAULT '0',",
                    "`cdsEnd` UNSIGNED INTEGER NOT NULL DEFAULT '0',",
                    "`exonCount` UNSIGNED INTEGER NOT NULL DEFAULT '0',",
                    "`exonStarts` TEXT NOT NULL,",
                    "`exonEnds` TEXT NOT NULL,",
                    "`proteinID` TEXT NOT NULL DEFAULT '',",
                    "`alignID` TEXT NOT NULL DEFAULT ''",
                    ")",collapse=" "
                ))
            },
            knownToRefSeq = {
                return(paste(
                    "CREATE TABLE",
                    "`knownToRefSeq` (",
                    "`name` TEXT NOT NULL DEFAULT '',",
                    "`value` TEXT NOT NULL DEFAULT ''",
                    ")",collapse=" "
                ))
            },
            refFlat = {
                return(paste("CREATE TABLE",
                    "`refFlat` (",
                    "`geneName` TEXT NOT NULL,",
                    "`name` TEXT NOT NULL,",
                    "`chrom` TEXT NOT NULL,",
                    "`strand` TEXT NOT NULL,",
                    "`txStart` UNSIGNED INTEGER NOT NULL,",
                    "`txEnd` UNSIGNED INTEGER NOT NULL,",
                    "`cdsStart` UNSIGNED INTEGER NOT NULL,",
                    "`cdsEnd` UNSIGNED INTEGER NOT NULL,",
                    "`exonCount` UNSIGNED INTEGER NOT NULL,",
                    "`exonStarts` TEXT NOT NULL,",
                    "`exonEnds` TEXT NOT NULL",
                    ")",collapse=" "
                ))
            },
            knownToEnsembl = {
                return(paste(
                    "CREATE TABLE",
                    "`knownToEnsembl` (",
                    "`name` TEXT NOT NULL DEFAULT '',",
                    "`value` TEXT NOT NULL DEFAULT ''",
                    ")",collapse=" "
                ))
            },
            ensemblSource = {
                return(paste(
                    "CREATE TABLE",
                    "`ensemblSource` (",
                    "`name` TEXT NOT NULL DEFAULT '',",
                    "`source` TEXT NOT NULL DEFAULT ''",
                    ")",collapse=" "
                ))
            },
            mgcGenes = {
                return(paste(
                    "CREATE TABLE `mgcGenes` (",
                    "`bin` UNSIGNED INTEGER NOT NULL,",
                    "`name` TEXT NOT NULL,",
                    "`chrom` TEXT NOT NULL,",
                    "`strand` TEXT NOT NULL,",
                    "`txStart` UNSIGNED INTEGER NOT NULL,",
                    "`txEnd` UNSIGNED INTEGER NOT NULL,",
                    "`cdsStart` UNSIGNED INTEGER NOT NULL,",
                    "`cdsEnd` UNSIGNED INTEGER NOT NULL,",
                    "`exonCount` UNSIGNED INTEGER NOT NULL,",
                    "`exonStarts` TEXT NOT NULL,",
                    "`exonEnds` TEXT NOT NULL,",
                    "`score` INTEGER DEFAULT NULL,",
                    "`name2` TEXT NOT NULL,",
                    "`cdsStartStat` TEXT NOT NULL,",
                    "`cdsEndStat` TEXT NOT NULL,",
                    "`exonFrames` TEXT NOT NULL",
                    ")",collapse=" "
                ))
            },
            ensemblToGeneName = {
                return(paste(
                    "CREATE TABLE",
                    "`knownToGeneName` (",
                    "`name` TEXT NOT NULL,",
                    "`value` TEXT NOT NULL",
                    ")",collapse=" "
                ))
            },
            flyBaseCanonical = {
                return(paste(
                    "CREATE TABLE",
                    "`flyBaseCanonical` (",
                    "`chrom` TEXT NOT NULL DEFAULT '',",
                    "`chromStart` INTEGER NOT NULL DEFAULT '0',",
                    "`chromEnd` INTEGER NOT NULL DEFAULT '0',",
                    "`clusterId` INTEGER unsigned NOT NULL DEFAULT '0',",
                    "`transcript` TEXT NOT NULL DEFAULT '',",
                    "`protein` TEXT NOT NULL DEFAULT ''",
                    ")",collapse=" "
                ))
            },
            flyBaseGene = {
                return(paste(
                    "CREATE TABLE",
                    "`flyBaseGene` (",
                    "`bin` UNSIGNED INTEGER NOT NULL DEFAULT '0',",
                    "`name` TEXT NOT NULL DEFAULT '',",
                    "`chrom` TEXT NOT NULL DEFAULT '',",
                    "`strand` TEXT NOT NULL DEFAULT '',",
                    "`txStart` UNSIGNED INTEGER NOT NULL DEFAULT '0',",
                    "`txEnd` UNSIGNED INTEGER NOT NULL DEFAULT '0',",
                    "`cdsStart` UNSIGNED INTEGER NOT NULL DEFAULT '0',",
                    "`cdsEnd` UNSIGNED INTEGER NOT NULL DEFAULT '0',",
                    "`exonCount` UNSIGNED INTEGER NOT NULL DEFAULT '0',",
                    "`exonStarts` TEXT NOT NULL,",
                    "`exonEnds` TEXT NOT NULL",
                    ")",collapse=" "
                ))
            },
            flyBaseToRefSeq = {
                return(paste(
                    "CREATE TABLE",
                    "`flyBaseToRefSeq` (",
                    "`name` TEXT NOT NULL DEFAULT '',",
                    "`value` TEXT NOT NULL DEFAULT ''",
                    ")",collapse=" "
                ))
            }
        )
    }
    else if (what=="fields") {
        switch(tab,
            knownCanonical = {
                return(c("chrom","chromStart","chromEnd","clusterId",
                "transcript","protein"))
            },
            knownGene = {
                return(c("name","chrom","strand","txStart","txEnd","cdsStart",
                    "cdsEnd","exonCount","exonStarts","exonEnds","proteinID",
                    "alignID"))
            },
            knownToRefSeq = {
                return(c("name","value"))
            },
            refFlat = {
                return(c("geneName","name","chrom","strand","txStart","txEnd",
                    "cdsStart","cdsEnd","exonCount","exonStarts","exonEnds"))
            },
            knownToEnsembl = {
                return(c("name","value"))
            },
            ensemblSource = {
                return(c("name","source"))
            },
            mgcGenes = {
                return(c("name","chrom","strand","txStart","txEnd","cdsStart",
                    "cdsEnd","exonCount","exonStarts","exonEnds","score",
                    "name2","cdsStartStat","cdsEndStat","exonFrames"
                ))
            },
            ensemblToGeneName = {
                return(c("name","value"))
            },
            flyBaseCanonical = {
                return(c("chrom","chromStart","chromEnd","clusterId",
                    "transcript","protein"))
            },
            flyBaseGene = {
                return(c("bin","name","chrom","strand","txStart","txEnd",
                    "cdsStart","cdsEnd","exonCount","exonStarts","exonEnds"))
            },
            flyBaseToRefSeq = {
                return(c("name","value"))
            }
        )
    }
}

getUcscQuery <- function(org,type,refdb="ucsc") {
    type <- tolower(type[1])
    org <- tolower(org[1])
    refdb <- tolower(refdb[1])
    checkTextArgs("type",type,c("gene","exon"))
    checkTextArgs("org",org,c("hg18","hg19","hg38","mm9","mm10","rn5","dm3",
        "danrer7","pantro4","susscr3"),multiarg=FALSE)
    checkTextArgs("refdb",refdb,c("ucsc","refseq"))
    switch(type,
        gene = {
            switch(refdb,
                ucsc = {
                    switch(org,
                        hg18 = {
                            return(paste("SELECT knownCanonical.chrom AS ",
                                "`chromosome`,`chromStart` AS `start`,",
                                "`chromEnd` AS `end`,`transcript` AS ",
                                "`gene_id`,0 AS `gc_content`,knownGene.strand ",
                                "AS `strand`,`geneName` AS `gene_name`,'NA' ",
                                "AS `biotype` FROM `knownCanonical` INNER ",
                                "JOIN `knownGene` ON ",
                                "knownCanonical.transcript=knownGene.name ",
                                "INNER JOIN `knownToRefSeq` ON ",
                                "knownCanonical.transcript=knownToRefSeq.name ",                                
                                "INNER JOIN `refFlat` ON ",
                                "knownToRefSeq.value=refFlat.name GROUP BY ",
                                "`gene_id` ORDER BY `chromosome`, `start`",
                                sep=""))
                        },
                        hg19 = {
                            return(paste("SELECT knownCanonical.chrom AS ",
                                "`chromosome`,`chromStart` AS `start`,",
                                "`chromEnd` AS `end`,`transcript` AS ",
                                "`gene_id`,0 AS `gc_content`,",
                                "knownGene.strand AS `strand`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`knownCanonical` INNER JOIN `knownGene` ON ",
                                "knownCanonical.transcript=knownGene.name ",
                                "INNER JOIN `knownToRefSeq` ON ",
                                "knownCanonical.transcript=knownToRefSeq.name ",
                                "INNER JOIN `knownToEnsembl` ON ",
                                "knownCanonical.transcript=knownToEnsembl.name",
                                " INNER JOIN `ensemblSource` ON ",
                                "knownToEnsembl.value=ensemblSource.name ",
                                "INNER JOIN `refFlat` ON ",
                                "knownToRefSeq.value=refFlat.name GROUP BY ",
                                "`gene_id` ORDER BY `chromosome`, `start`",
                                sep=""))
                        },
                        hg38 = {
                            return(paste("SELECT knownCanonical.chrom AS ",
                                "`chromosome`,`chromStart` AS `start`,",
                                "`chromEnd` AS `end`,`transcript` AS ",
                                "`gene_id`,0 AS `gc_content`,knownGene.strand ",
                                "AS `strand`,`geneName` AS `gene_name`,'NA' ",
                                "AS `biotype` FROM `knownCanonical` INNER ",
                                "JOIN `knownGene` ON ",
                                "knownCanonical.transcript=knownGene.name ",
                                "INNER JOIN `knownToRefSeq` ON ",
                                "knownCanonical.transcript=knownToRefSeq.name ",                                
                                "INNER JOIN `refFlat` ON ",
                                "knownToRefSeq.value=refFlat.name GROUP BY ",
                                "`gene_id` ORDER BY `chromosome`, `start`",
                                sep=""))
                            # Should be the same as hg19 but is like hg18
                        },
                        mm9 = {
                            return(paste("SELECT knownCanonical.chrom AS ",
                                "`chromosome`,`chromStart` AS `start`,",
                                "`chromEnd` AS `end`,`transcript` AS ",
                                "`gene_id`,0 AS `gc_content`,",
                                "knownGene.strand AS `strand`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`knownCanonical` INNER JOIN `knownGene` ON ",
                                "knownCanonical.transcript=knownGene.name ",
                                "INNER JOIN `knownToRefSeq` ON ",
                                "knownCanonical.transcript=knownToRefSeq.name ",
                                "INNER JOIN `knownToEnsembl` ON ",
                                "knownCanonical.transcript=knownToEnsembl.name",
                                " INNER JOIN `ensemblSource` ON ",
                                "knownToEnsembl.value=ensemblSource.name ",
                                "INNER JOIN `refFlat` ON ",
                                "knownToRefSeq.value=refFlat.name GROUP BY ",
                                "`gene_id` ORDER BY `chromosome`, `start`",
                                sep=""))
                        },
                        mm10 = {
                            return(paste("SELECT knownCanonical.chrom AS ",
                                "`chromosome`,`chromStart` AS `start`,",
                                "`chromEnd` AS `end`,`transcript` AS ",
                                "`gene_id`,0 AS `gc_content`,",
                                "knownGene.strand AS `strand`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`knownCanonical` INNER JOIN `knownGene` ON ",
                                "knownCanonical.transcript=knownGene.name ",
                                "INNER JOIN `knownToRefSeq` ON ",
                                "knownCanonical.transcript=knownToRefSeq.name ",
                                "INNER JOIN `knownToEnsembl` ON ",
                                "knownCanonical.transcript=knownToEnsembl.name",
                                " INNER JOIN `ensemblSource` ON ",
                                "knownToEnsembl.value=ensemblSource.name ",
                                "INNER JOIN `refFlat` ON ",
                                "knownToRefSeq.value=refFlat.name GROUP BY ",
                                "`gene_id` ORDER BY `chromosome`, `start`",
                                sep=""))
                        },
                        rn5 = {
                            return(paste("SELECT mgcGenes.chrom AS ",
                                "`chromosome`,`txStart` AS `start`,`txEnd` ",
                                "AS `end`,mgcGenes.name AS `gene_id`,0 AS ",
                                "`gc_content`,mgcGenes.strand AS `strand`,",
                                "`name2` AS `gene_name`,`source` AS `biotype` ",
                                "FROM `mgcGenes` INNER JOIN ",
                                "`ensemblToGeneName` ON ",
                                "mgcGenes.name2=ensemblToGeneName.value INNER ",
                                "JOIN `ensemblSource` ON ",
                                "ensemblToGeneName.name=ensemblSource.name ",
                                "GROUP BY `gene_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        dm3 = {
                            return(paste("SELECT flyBaseCanonical.chrom AS ",
                                "`chromosome`,`chromStart` AS `start`,",
                                "`chromEnd` AS `end`,`transcript` AS ",
                                "`gene_id`,0 AS `gc_content`,",
                                "flyBaseGene.strand AS `strand`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`flyBaseCanonical` INNER JOIN `flyBaseGene` ",
                                "ON flyBaseCanonical.transcript=",
                                "flyBaseGene.name INNER JOIN ",
                                "`flyBaseToRefSeq` ON ",
                                "flyBaseCanonical.transcript=",
                                "flyBaseToRefSeq.name INNER JOIN `refFlat` ON ",
                                "flyBaseToRefSeq.value=refFlat.name INNER ",
                                "JOIN `ensemblToGeneName` ON ",
                                "ensemblToGeneName.value=refFlat.geneName ",
                                "INNER JOIN `ensemblSource` ON ",
                                "ensemblToGeneName.name=ensemblSource.name ",
                                "GROUP BY `gene_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        danrer7 = {
                            return(paste("SELECT mgcGenes.chrom AS ",
                                "`chromosome`,`txStart` AS `start`,`txEnd` ",
                                "AS `end`,mgcGenes.name AS `gene_id`,0 AS ",
                                "`gc_content`,mgcGenes.strand AS `strand`,",
                                "`name2` AS `gene_name`,`source` AS `biotype` ",
                                "FROM `mgcGenes` INNER JOIN ",
                                "`ensemblToGeneName` ON ",
                                "mgcGenes.name2=ensemblToGeneName.value INNER ",
                                "JOIN `ensemblSource` ON ",
                                "ensemblToGeneName.name=ensemblSource.name ",
                                "GROUP BY `gene_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        pantro4 = {
                            warning("No UCSC Genome annotation for Pan ",
                                "troglodytes! Will use RefSeq instead...",
                                immediate.=TRUE)
                            return(paste("SELECT refFlat.chrom AS ",
                                "`chromosome`,refFlat.txStart AS `start`,",
                                "refFlat.txEnd AS `end`,refFlat.name AS ",
                                "`gene_id`,0 AS `gc_content`,",
                                "refFlat.strand AS `strand`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`refFlat` INNER JOIN `ensemblToGeneName` ON ",
                                "refFlat.geneName=ensemblToGeneName.value ",
                                "INNER JOIN `ensemblSource` ON ",
                                "ensemblToGeneName.name=ensemblSource.name ",
                                "GROUP BY `gene_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        susscr3 = {
                            warning("No UCSC Genome annotation for Sus ",
                                "scrofa! Will use RefSeq instead...",
                                immediate.=TRUE)
                            return(paste(
                                "SELECT refFlat.chrom AS `chromosome`,",
                                "refFlat.txStart AS `start`, refFlat.txEnd AS ",
                                "`end`, refFlat.name AS `gene_id`, 0 AS ",
                                "`gc_content`, refFlat.strand AS `strand`,",
                                "`geneName` AS `gene_name`, `source` AS ",
                                "`biotype` FROM `refFlat` INNER JOIN ",
                                "`ensemblToGeneName` ON ",
                                "refFlat.geneName=ensemblToGeneName.value ",
                                "INNER JOIN `ensemblSource` ON ",
                                "ensemblToGeneName.name=ensemblSource.name ",
                                "GROUP BY `gene_id` ORDER BY `chromosome`,",
                                "`start`",
                                sep=""
                            ))
                        }
                    )
                },
                refseq = {
                    switch(org,
                        hg18 = {
                            return(paste("SELECT  refFlat.chrom AS ",
                                "`chromosome`,refFlat.txStart AS `start`,",
                                "refFlat.txEnd AS `end`,refFlat.name AS ",
                                "`gene_id`,0 AS `gc_content`,refFlat.strand ",
                                "AS `strand`,`geneName` AS `gene_name`,'NA' ",
                                "AS `biotype` FROM `refFlat` INNER JOIN ",
                                "`knownToRefSeq` ON ",
                                "refFlat.name=knownToRefSeq.value INNER JOIN ",
                                "`knownCanonical` ON ",
                                "knownToRefSeq.name=knownCanonical.transcript ",
                                "GROUP BY refFlat.name ORDER BY `chromosome`,",
                                " `start`",
                                sep=""))
                        },
                        hg19 = {
                            return(paste("SELECT  refFlat.chrom AS ",
                                "`chromosome`,refFlat.txStart AS `start`,",
                                "refFlat.txEnd AS `end`,refFlat.name AS ",
                                "`gene_id`,0 AS `gc_content`,",
                                "refFlat.strand AS `strand`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`refFlat` INNER JOIN `knownToRefSeq` ON ",
                                "refFlat.name=knownToRefSeq.value INNER ",
                                "JOIN `knownCanonical` ON ",
                                "knownToRefSeq.name=knownCanonical.transcript ",
                                "INNER JOIN `knownToEnsembl` ON ",
                                "knownCanonical.transcript=knownToEnsembl.name",
                                " INNER JOIN `ensemblSource` ON ",
                                "knownToEnsembl.value=ensemblSource.name ",
                                "GROUP BY refFlat.name ORDER BY `chromosome`,",
                                " `start`",
                                sep=""))
                        },
                        hg38 = {
                            return(paste("SELECT  refFlat.chrom AS ",
                                "`chromosome`,refFlat.txStart AS `start`,",
                                "refFlat.txEnd AS `end`,refFlat.name AS ",
                                "`gene_id`,0 AS `gc_content`,refFlat.strand ",
                                "AS `strand`,`geneName` AS `gene_name`,'NA' ",
                                "AS `biotype` FROM `refFlat` INNER JOIN ",
                                "`knownToRefSeq` ON ",
                                "refFlat.name=knownToRefSeq.value INNER JOIN ",
                                "`knownCanonical` ON ",
                                "knownToRefSeq.name=knownCanonical.transcript ",
                                "GROUP BY refFlat.name ORDER BY `chromosome`,",
                                " `start`",
                                sep=""))
                            # Should be the same as hg19 but is as hg18
                        },
                        mm9 = {
                            return(paste("SELECT  refFlat.chrom AS ",
                                "`chromosome`,refFlat.txStart AS `start`,",
                                "refFlat.txEnd AS `end`,refFlat.name AS ",
                                "`gene_id`,0 AS `gc_content`,",
                                "refFlat.strand AS `strand`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`refFlat` INNER JOIN `knownToRefSeq` ON ",
                                "refFlat.name=knownToRefSeq.value INNER ",
                                "JOIN `knownCanonical` ON ",
                                "knownToRefSeq.name=knownCanonical.transcript ",
                                "INNER JOIN `knownToEnsembl` ON ",
                                "knownCanonical.transcript=knownToEnsembl.name",
                                " INNER JOIN `ensemblSource` ON ",
                                "knownToEnsembl.value=ensemblSource.name ",
                                "GROUP BY refFlat.name ORDER BY `chromosome`,",
                                " `start`",
                                sep=""))
                        },
                        mm10 = {
                            return(paste("SELECT  refFlat.chrom AS ",
                                "`chromosome`,refFlat.txStart AS `start`,",
                                "refFlat.txEnd AS `end`,refFlat.name AS ",
                                "`gene_id`,0 AS `gc_content`,",
                                "refFlat.strand AS `strand`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`refFlat` INNER JOIN `knownToRefSeq` ON ",
                                "refFlat.name=knownToRefSeq.value INNER ",
                                "JOIN `knownCanonical` ON ",
                                "knownToRefSeq.name=knownCanonical.transcript ",
                                "INNER JOIN `knownToEnsembl` ON ",
                                "knownCanonical.transcript=knownToEnsembl.name",
                                " INNER JOIN `ensemblSource` ON ",
                                "knownToEnsembl.value=ensemblSource.name ",
                                "GROUP BY refFlat.name ORDER BY `chromosome`,",
                                " `start`",
                                sep=""))
                        },
                        rn5 = {
                            return(paste("SELECT refFlat.chrom AS ",
                                "`chromosome`,refFlat.txStart AS `start`,",
                                "refFlat.txEnd AS `end`,refFlat.name AS ",
                                "`gene_id`,0 AS `gc_content`,",
                                "refFlat.strand AS `strand`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`refFlat` INNER JOIN `ensemblToGeneName` ON ",
                                "refFlat.geneName=ensemblToGeneName.value ",
                                "INNER JOIN `ensemblSource` ON ",
                                "ensemblToGeneName.name=ensemblSource.name ",
                                "GROUP BY `gene_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        dm3 = {
                            return(paste("SELECT refFlat.chrom AS ",
                                "`chromosome`,refFlat.txStart AS `start`,",
                                "refFlat.txEnd AS `end`,refFlat.name AS ",
                                "`gene_id`,0 AS `gc_content`,refFlat.strand ",
                                "AS `strand`,`geneName` AS `gene_name`,",
                                "`source` AS `biotype` FROM `refFlat` INNER ",
                                "JOIN `ensemblToGeneName` ON ",
                                "ensemblToGeneName.value=refFlat.geneName ",
                                "INNER JOIN `ensemblSource` ON ",
                                "ensemblToGeneName.name=ensemblSource.name ",
                                "GROUP BY `gene_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        danrer7 = {
                            return(paste("SELECT refFlat.chrom AS ",
                                "`chromosome`,refFlat.txStart AS `start`,",
                                "refFlat.txEnd AS `end`,refFlat.name AS ",
                                "`gene_id`,0 AS `gc_content`,",
                                "refFlat.strand AS `strand`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`refFlat` INNER JOIN `ensemblToGeneName` ON ",
                                "refFlat.geneName=ensemblToGeneName.value ",
                                "INNER JOIN `ensemblSource` ON ",
                                "ensemblToGeneName.name=ensemblSource.name ",
                                "GROUP BY `gene_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        pantro4 = {
                            return(paste("SELECT refFlat.chrom AS ",
                                "`chromosome`,refFlat.txStart AS `start`,",
                                "refFlat.txEnd AS `end`,refFlat.name AS ",
                                "`gene_id`,0 AS `gc_content`,",
                                "refFlat.strand AS `strand`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`refFlat` INNER JOIN `ensemblToGeneName` ON ",
                                "refFlat.geneName=ensemblToGeneName.value ",
                                "INNER JOIN `ensemblSource` ON ",
                                "ensemblToGeneName.name=ensemblSource.name ",
                                "GROUP BY `gene_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        susscr3 = {
                            return(paste(
                                "SELECT refFlat.chrom AS `chromosome`,",
                                "refFlat.txStart AS `start`, refFlat.txEnd AS ",
                                "`end`, refFlat.name AS `gene_id`, 0 AS ",
                                "`gc_content`, refFlat.strand AS `strand`,",
                                "`geneName` AS `gene_name`, `source` AS ",
                                "`biotype` FROM `refFlat` INNER JOIN ",
                                "`ensemblToGeneName` ON ",
                                "refFlat.geneName=ensemblToGeneName.value ",
                                "INNER JOIN `ensemblSource` ON ",
                                "ensemblToGeneName.name=ensemblSource.name ",
                                "GROUP BY `gene_id` ORDER BY `chromosome`,",
                                "`start`",
                                sep=""
                            ))
                        }
                    )
                }
            )
        },
        exon = {
            switch(refdb,
                ucsc = {
                    switch(org,
                        hg18 = {
                            return(paste("SELECT knownGene.chrom AS ",
                                "`chromosome`,knownGene.exonStarts AS `start`,",
                                "knownGene.exonEnds AS `end`,knownGene.name ",
                                "AS `exon_id`,knownGene.strand AS `strand`,",
                                "`transcript` AS `gene_id`,`geneName` AS ",
                                "`gene_name`,'NA' AS `biotype` FROM ",
                                "`knownGene` INNER JOIN `knownCanonical` ON ",
                                "knownGene.name=knownCanonical.transcript ",
                                "INNER JOIN `knownToRefSeq` ON ",
                                "knownCanonical.transcript=knownToRefSeq.name ",
                                "INNER JOIN `refFlat` ON ",
                                "knownToRefSeq.value=refFlat.name GROUP BY ",
                                "knownGene.name ORDER BY `chromosome`, `start`",
                                sep=""))
                        },
                        hg19 = {
                            return(paste("SELECT knownGene.chrom AS ",
                                "`chromosome`,knownGene.exonStarts AS `start`,",
                                "knownGene.exonEnds AS `end`,knownGene.name ",
                                "AS `exon_id`,knownGene.strand AS `strand`,",
                                "`transcript` AS `gene_id`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`knownGene` INNER JOIN `knownCanonical` ON ",
                                "knownGene.name=knownCanonical.transcript ",
                                "INNER JOIN `knownToRefSeq` ON ",
                                "knownCanonical.transcript=knownToRefSeq.name ",
                                "INNER JOIN `knownToEnsembl` ON ",
                                "knownCanonical.transcript=knownToEnsembl.name",
                                " INNER JOIN `ensemblSource` ON ",
                                "knownToEnsembl.value=ensemblSource.name ",
                                "INNER JOIN `refFlat` ON ",
                                "knownToRefSeq.value=refFlat.name GROUP BY ",
                                "knownGene.name ORDER BY `chromosome`, `start`",
                                sep=""))
                        },
                        hg38 = {
                            return(paste("SELECT knownGene.chrom AS ",
                                "`chromosome`,knownGene.exonStarts AS `start`,",
                                "knownGene.exonEnds AS `end`,knownGene.name ",
                                "AS `exon_id`,knownGene.strand AS `strand`,",
                                "`transcript` AS `gene_id`,`geneName` AS ",
                                "`gene_name`,'NA' AS `biotype` FROM ",
                                "`knownGene` INNER JOIN `knownCanonical` ON ",
                                "knownGene.name=knownCanonical.transcript ",
                                "INNER JOIN `knownToRefSeq` ON ",
                                "knownCanonical.transcript=knownToRefSeq.name ",
                                "INNER JOIN `refFlat` ON ",
                                "knownToRefSeq.value=refFlat.name GROUP BY ",
                                "knownGene.name ORDER BY `chromosome`, `start`",
                                sep=""))
                            # Should be the same as hg19 but is as hg18
                        },
                        mm9 = {
                            return(paste("SELECT knownGene.chrom AS ",
                                "`chromosome`,knownGene.exonStarts AS `start`,",
                                "knownGene.exonEnds AS `end`,knownGene.name ",
                                "AS `exon_id`,knownGene.strand AS `strand`,",
                                "`transcript` AS `gene_id`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`knownGene` INNER JOIN `knownCanonical` ON ",
                                "knownGene.name=knownCanonical.transcript ",
                                "INNER JOIN `knownToRefSeq` ON ",
                                "knownCanonical.transcript=knownToRefSeq.name ",
                                "INNER JOIN `knownToEnsembl` ON ",
                                "knownCanonical.transcript=knownToEnsembl.name",
                                " INNER JOIN `ensemblSource` ON ",
                                "knownToEnsembl.value=ensemblSource.name ",
                                "INNER JOIN `refFlat` ON ",
                                "knownToRefSeq.value=refFlat.name GROUP BY ",
                                "knownGene.name ORDER BY `chromosome`, `start`",
                                sep=""))
                        },
                        mm10 = {
                            return(paste("SELECT knownGene.chrom AS ",
                                "`chromosome`,knownGene.exonStarts AS `start`,",
                                "knownGene.exonEnds AS `end`,knownGene.name ",
                                "AS `exon_id`,knownGene.strand AS `strand`,",
                                "`transcript` AS `gene_id`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`knownGene` INNER JOIN `knownCanonical` ON ",
                                "knownGene.name=knownCanonical.transcript ",
                                "INNER JOIN `knownToRefSeq` ON ",
                                "knownCanonical.transcript=knownToRefSeq.name ",
                                "INNER JOIN `knownToEnsembl` ON ",
                                "knownCanonical.transcript=knownToEnsembl.name",
                                " INNER JOIN `ensemblSource` ON ",
                                "knownToEnsembl.value=ensemblSource.name ",
                                "INNER JOIN `refFlat` ON ",
                                "knownToRefSeq.value=refFlat.name GROUP BY ",
                                "knownGene.name ORDER BY `chromosome`, `start`",
                                sep=""))
                        },
                        rn5 = {
                            return(paste("SELECT mgcGenes.chrom AS ",
                                "`chromosome`,`exonStarts` AS `start`,",
                                "`exonEnds` AS `end`,mgcGenes.name AS ",
                                "`exon_id`,mgcGenes.strand AS `strand`,",
                                "mgcGenes.name AS `gene_id`,`name2` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`mgcGenes` INNER JOIN `ensemblToGeneName` ON ",
                                "mgcGenes.name2=ensemblToGeneName.value INNER ",
                                "JOIN `ensemblSource` ON ",
                                "ensemblToGeneName.name=ensemblSource.name ",
                                "GROUP BY `gene_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        dm3 = {
                            return(paste("SELECT flyBaseCanonical.chrom AS ",
                                "`chromosome`,flyBaseGene.exonStarts AS ",
                                "`start`,flyBaseGene.exonEnds AS `end`,",
                                "`transcript` AS `exon_id`,flyBaseGene.strand ",
                                "AS `strand`,`transcript` AS `gene_id`,",
                                "`geneName` AS `gene_name`,`source` AS ",
                                "`biotype` FROM `flyBaseCanonical` INNER JOIN ",
                                "`flyBaseGene` ON ",
                                "flyBaseCanonical.transcript=flyBaseGene.name ",
                                "INNER JOIN `flyBaseToRefSeq` ON ",
                                "flyBaseCanonical.transcript=",
                                "flyBaseToRefSeq.name INNER JOIN `refFlat` ON ",
                                "flyBaseToRefSeq.value=refFlat.name ",
                                "INNER JOIN `ensemblToGeneName` ON ",
                                "ensemblToGeneName.value=refFlat.geneName ",
                                "INNER JOIN `ensemblSource` ON ",
                                "ensemblToGeneName.name=ensemblSource.name ",
                                "GROUP BY `gene_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        danrer7 = {
                            return(paste("SELECT mgcGenes.chrom AS ",
                                "`chromosome`,`exonStarts` AS `start`,",
                                "`exonEnds` AS `end`,mgcGenes.name AS ",
                                "`exon_id`,mgcGenes.strand AS `strand`,",
                                "mgcGenes.name AS `gene_id`,`name2` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`mgcGenes` INNER JOIN `ensemblToGeneName` ON ",
                                "mgcGenes.name2=ensemblToGeneName.value INNER ",
                                "JOIN `ensemblSource` ON ",
                                "ensemblToGeneName.name=ensemblSource.name ",
                                "GROUP BY `gene_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        pantro4 = {
                            warning("No UCSC Genome annotation for Pan ",
                                "troglodytes! Will use RefSeq instead...",
                                immediate.=TRUE)
                            return(paste("SELECT refFlat.chrom AS ",
                                "`chromosome`,refFlat.exonStarts AS `start`,",
                                "refFlat.exonEnds AS `end`,refFlat.name AS ",
                                "`exon_id`,refFlat.strand AS `strand`,",
                                "refFlat.name AS `gene_id`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`refFlat` INNER JOIN `ensemblToGeneName` ON ",
                                "refFlat.geneName=ensemblToGeneName.value ",
                                "INNER JOIN `ensemblSource` ON ",
                                "ensemblToGeneName.name=ensemblSource.name ",
                                "GROUP BY `exon_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        susscr3 = {
                            warning("No UCSC Genome annotation for Sus ",
                                "scrofa! Will use RefSeq instead...",
                                immediate.=TRUE)
                            return(paste("SELECT refFlat.chrom AS ",
                                "`chromosome`,refFlat.exonStarts AS `start`,",
                                "refFlat.exonEnds AS `end`,refFlat.name AS ",
                                "`exon_id`,refFlat.strand AS `strand`,",
                                "refFlat.name AS `gene_id`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`refFlat` INNER JOIN `ensemblToGeneName` ON ",
                                "refFlat.geneName=ensemblToGeneName.value ",
                                "INNER JOIN `ensemblSource` ON ",
                                "ensemblToGeneName.name=ensemblSource.name ",
                                "GROUP BY `exon_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        }
                    )
                },
                refseq = {
                    switch(org,
                        hg18 = {
                            return(paste("SELECT refFlat.chrom AS ",
                                "`chromosome`,refFlat.exonStarts AS `start`,",
                                "refFlat.exonEnds AS `end`,refFlat.name AS ",
                                "`exon_id`,refFlat.strand AS `strand`,",
                                "refFlat.name AS `gene_id`,`geneName` AS ",
                                "`gene_name`,'NA' AS `biotype` FROM `refFlat` ",
                                "INNER JOIN `knownToRefSeq` ON ",
                                "refFlat.name=knownToRefSeq.value INNER JOIN ",
                                "`knownCanonical` ON ",
                                "knownToRefSeq.name=knownCanonical.transcript ",
                                "GROUP BY `exon_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        hg19 = {
                            return(paste("SELECT refFlat.chrom AS ",
                                "`chromosome`,refFlat.exonStarts AS `start`,",
                                "refFlat.exonEnds AS `end`,refFlat.name AS ",
                                "`exon_id`,refFlat.strand AS `strand`,",
                                "refFlat.name AS `gene_id`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`refFlat` INNER JOIN `knownToRefSeq` ON ",
                                "refFlat.name=knownToRefSeq.value INNER JOIN ",
                                "`knownCanonical` ON ",
                                "knownToRefSeq.name=knownCanonical.transcript ",
                                "INNER JOIN `knownToEnsembl` ON ",
                                "knownCanonical.transcript=knownToEnsembl.name",
                                " INNER JOIN `ensemblSource` ON ",
                                "knownToEnsembl.value=ensemblSource.name ",
                                "GROUP BY `exon_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        hg38 = {
                            return(paste("SELECT refFlat.chrom AS ",
                                "`chromosome`,refFlat.exonStarts AS `start`,",
                                "refFlat.exonEnds AS `end`,refFlat.name AS ",
                                "`exon_id`,refFlat.strand AS `strand`,",
                                "refFlat.name AS `gene_id`,`geneName` AS ",
                                "`gene_name`,'NA' AS `biotype` FROM `refFlat` ",
                                "INNER JOIN `knownToRefSeq` ON ",
                                "refFlat.name=knownToRefSeq.value INNER JOIN ",
                                "`knownCanonical` ON ",
                                "knownToRefSeq.name=knownCanonical.transcript ",
                                "GROUP BY `exon_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                            # Should be the same as hg19 but is as hg18
                        },
                        mm9 = {
                            return(paste("SELECT refFlat.chrom AS ",
                                "`chromosome`,refFlat.exonStarts AS `start`,",
                                "refFlat.exonEnds AS `end`,refFlat.name AS ",
                                "`exon_id`,refFlat.strand AS `strand`,",
                                "refFlat.name AS `gene_id`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`refFlat` INNER JOIN `knownToRefSeq` ON ",
                                "refFlat.name=knownToRefSeq.value INNER JOIN ",
                                "`knownCanonical` ON ",
                                "knownToRefSeq.name=knownCanonical.transcript ",
                                "INNER JOIN `knownToEnsembl` ON ",
                                "knownCanonical.transcript=knownToEnsembl.name",
                                " INNER JOIN `ensemblSource` ON ",
                                "knownToEnsembl.value=ensemblSource.name ",
                                "GROUP BY `exon_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        mm10 = {
                            return(paste("SELECT refFlat.chrom AS ",
                                "`chromosome`,refFlat.exonStarts AS `start`,",
                                "refFlat.exonEnds AS `end`,refFlat.name AS ",
                                "`exon_id`,refFlat.strand AS `strand`,",
                                "refFlat.name AS `gene_id`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`refFlat` INNER JOIN `knownToRefSeq` ON ",
                                "refFlat.name=knownToRefSeq.value INNER JOIN ",
                                "`knownCanonical` ON ",
                                "knownToRefSeq.name=knownCanonical.transcript ",
                                "INNER JOIN `knownToEnsembl` ON ",
                                "knownCanonical.transcript=knownToEnsembl.name",
                                " INNER JOIN `ensemblSource` ON ",
                                "knownToEnsembl.value=ensemblSource.name ",
                                "GROUP BY `exon_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        rn5 = {
                            return(paste("SELECT refFlat.chrom AS ",
                                "`chromosome`,refFlat.exonStarts AS `start`,",
                                "refFlat.exonEnds AS `end`,refFlat.name AS ",
                                "`exon_id`,refFlat.strand AS `strand`,",
                                "refFlat.name AS `gene_id`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`refFlat` INNER JOIN `ensemblToGeneName` ON ",
                                "refFlat.geneName=ensemblToGeneName.value ",
                                "INNER JOIN `ensemblSource` ON ",
                                "ensemblToGeneName.name=ensemblSource.name ",
                                "GROUP BY `exon_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        dm3 = {
                            return(paste("SELECT refFlat.chrom AS ",
                                "`chromosome`,refFlat.exonStarts AS `start`,",
                                "refFlat.exonEnds AS `end`,refFlat.name AS ",
                                "`exon_id`,refFlat.strand AS `strand`,",
                                "refFlat.name AS `gene_id`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`refFlat` INNER JOIN `ensemblToGeneName` ON ",
                                "ensemblToGeneName.value=refFlat.geneName ",
                                "INNER JOIN `ensemblSource` ON ",
                                "ensemblToGeneName.name=ensemblSource.name ",
                                "GROUP BY `gene_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        danrer7 = {
                            return(paste("SELECT refFlat.chrom AS ",
                                "`chromosome`,refFlat.exonStarts AS `start`,",
                                "refFlat.exonEnds AS `end`,refFlat.name AS ",
                                "`exon_id`,refFlat.strand AS `strand`,",
                                "refFlat.name AS `gene_id`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`refFlat` INNER JOIN `ensemblToGeneName` ON ",
                                "refFlat.geneName=ensemblToGeneName.value ",
                                "INNER JOIN `ensemblSource` ON ",
                                "ensemblToGeneName.name=ensemblSource.name ",
                                "GROUP BY `exon_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        pantro4 = {
                            return(paste("SELECT refFlat.chrom AS ",
                                "`chromosome`,refFlat.exonStarts AS `start`,",
                                "refFlat.exonEnds AS `end`,refFlat.name AS ",
                                "`exon_id`,refFlat.strand AS `strand`,",
                                "refFlat.name AS `gene_id`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`refFlat` INNER JOIN `ensemblToGeneName` ON ",
                                "refFlat.geneName=ensemblToGeneName.value ",
                                "INNER JOIN `ensemblSource` ON ",
                                "ensemblToGeneName.name=ensemblSource.name ",
                                "GROUP BY `exon_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        susscr3 = {
                            return(paste("SELECT refFlat.chrom AS ",
                                "`chromosome`,refFlat.exonStarts AS `start`,",
                                "refFlat.exonEnds AS `end`,refFlat.name AS ",
                                "`exon_id`,refFlat.strand AS `strand`,",
                                "refFlat.name AS `gene_id`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`refFlat` INNER JOIN `ensemblToGeneName` ON ",
                                "refFlat.geneName=ensemblToGeneName.value ",
                                "INNER JOIN `ensemblSource` ON ",
                                "ensemblToGeneName.name=ensemblSource.name ",
                                "GROUP BY `exon_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        }
                    )
                }
            )
        }
    )
}

getUcscCredentials <- function() {
    return(c(
        host="genome-mysql.cse.ucsc.edu",
        user="genome",
        password=""
    ))
}
