<!-- badges: start -->
  [![Travis build status](https://travis-ci.org/pmoulos/recoup.svg?branch=master)](https://travis-ci.org/pmoulos/recoup)
  ![GitHub](https://img.shields.io/github/license/pmoulos/recoup)
  ![GitHub repo size](https://img.shields.io/github/repo-size/pmoulos/recoup)
  ![GitHub issues](https://img.shields.io/github/issues/pmoulos/recoup)
<!-- badges: end -->

Genomic coverages remastered! The recoup package.
==================================================

The explosion of the usage of Next Generation Sequencing techniques during the 
past few years due to the seemingly endless portfolio of applications has 
created the need for new NGS data analytics tools which are able to offer 
comprehensive and at the same time flexible visualizations under several 
experimental settings and factors. An established visualization tool in NGS
experiments is the visualization of the signal created by short reads after the
application of every NGS protocol. Genome Browsers (e.g. the UCSC Genome 
Browser) serve very well this purpose considering single genomic areas. They are
very good when it comes to the visualization of the abundance of a single or a 
few genes or the strength of a few protein-DNA interaction sites. However, when 
it comes to the visualization of average signal profiles over multiple genomic 
locations (gene regions or others like DNA methylation sites or transcription
factor binding sites), Genome Browsers fail to depict such information in a 
comprehensive way. Furthermore, they cannot visualize average signal profiles of
factored data, for example a set of genes categorized in high, medium and low
expression or even by strand and they cannot visualize all signals of interest
mapped on all the genomic regions of interest in a compact way, something that
can be done using for example a heatmap.

In such cases, bioinformaticians use several toolkits like BEDTools and 
facilities from R/Bioconductor to read/import short reads and overlap them with 
genomic regions of interest, summarize them by binning/averaging overlaps to 
control the resolution of final graphics etc. This procedure often requires the
usage of multiple tools and custom scripts with multiple steps. One of the most
comprehensive and easy-to-use tools up to date is 
[ngs.plot](https://github.com/shenlab-sinai/ngsplot). It is sufficiently fast 
for most applications and has a low memory footprint which allows users to run 
it in low end machines too. It is command line based and users can run it by 
using a simple configuration file most of the times, has a rich database of 
genome annotation and genomic features and uses R/Bioconductor for underlying 
calculations and plotting of profiles. However, ngs.plot is not up to date with 
modern R graphics systems like `ggplot2` and `ComplexHeatmap`. As a result, 
among others, it is impossible to create faceted genomic profiles using a 
statistical design and in such cases, a lot of additional manual work and 
computational time is required in order to reach the desired outcomes. The same 
applies to heatmap profiles. Furthermore, the resolution of genomic profiles 
(e.g. per base coverage or per bin of base-pairs coverage) cannot be controlled 
and this can cause problems in cases where extreme levels of resolution (e.g. 
DNAse-Seq experiments) is required to reach meaningful biological conclusions. 
Last but not least, ngs.plot requires a not so straightforward setup in order 
to run, does not run in a unified working environment (e.g. R environment with
its graphics displaying mechanisms) and in some cases produces oversized and 
complex output.

The recoup package comes to fill such gaps by stepping on the shoulders of 
giants. It uses the now standardized and stable Bioconductor facilities to read
and import short reads from BAM/BED files and also modern R graphics systems,
namely [ggplot2](http://ggplot2.org/) and
[ComplexHeatmap](https://github.com/jokergoo/ComplexHeatmap) in order to create 
comprehensive averaged genomic profiles and genomic profile heatmaps. In 
addition it offers a lot of (easy to use) customization options and automation 
at various levels. Inexperienced users can gather their data in a simple text 
file and just choose one of the supported organisms and recoup does the rest 
for them. More experienced users can play with several options and provide more 
flexible input so as to produce optimal results. This vignette, although it 
covers basic usage of the package, it offers the basis for more sophisticated 
usage. recoup is not as fast as ngs.plot but we are working on this! Also, 
recoup is not here to replace other more mature packages. It is here to offer 
more options to users that need more sophisticated genomic profile 
visualizations. Finally, it offers a very flexible way to reuse genomic profiles
whose calculations may be computationally and time expensive by offering a smart
way to recognize which parameters have changed and acting based on these.

## 1. Getting started

Detailed instructions on how to run the recoup genomic profile creation 
pipeline can be found under the main documentation of the package:

```
library(recoup)
help(recoup)
```

Briefly, to run recoup you need:

* A set of BAM or BED files that contain the aligned short reads from a protein-
DNA interaction experiment (ChIP-Seq from transcrption factors, DNA methylation
signals, DNAse-Seq etc.) or gene expression experiments (RNA-Seq, spliced or
unspliced). Theoretically, alignments from other protocols can be used (e.g. to
measure the coverage of Exome-Seq).
* A set of reference genomic regions to calculate average profiles and heatmaps
over. These can be provided as a BED-like text file with a header (see the data
attached to the package, look at `recoup` man page). They can also be provided
as an organism version keyword, e.g. `hg19` or `mm9` and the respective regions
will be either retrieved from a local `recoup` annotation database setup, or
downloaded on the fly (takes significantly more time as some extra operations 
are required).
* A design file in the case that you wish to categorize your profiles (e.g. a 
set of H3K27me1 or H4K20me1 profiles categorized by gene transcription or 
expression levels -high, medium, low, or different levels of binding strength
of a transcription factor close to the TSS). The design file should have in the
first column unique identifiers which should be the same, a subset or a superset
of those in the reference genomic regions.

The package contains a small dataset which serves only for package building and
testing purposes as well as complying with Bioconductor guidelines. These data
are useful for the user to check how the input data to recoup should look like.
For a more complete test dataset (a small one) have a look and download from
[here](https://drive.google.com/file/d/0BxxrqIl3Nb0NSVNqdGNPa3M3cnc/view?usp=sharing)

## 2. Getting some data

In order to run smoothly the rest of the examples in these vignettes and produce
some realistic results, you need to download a set of example BAM files, genomic
regions and design files from 
[here](https://drive.google.com/file/d/0BxxrqIl3Nb0NSVNqdGNPa3M3cnc/view?usp=sharing). 
Following, a description of each file in the archive (the tissue is always 
liver):

* For the ChIP-Seq example
    + _WT_H4K20me1.bam_: H4K20me1 signals from wild type adult mice, chromosome
    12 only
    + _WT_H4K20me1.bam.bai_: The index of the above
    + _Set8KO_H4K20me1.bam_: H4K20me1 signals from adult mice where the Set8 
    (Pr-Set7) gene is knocked out, chromosome 12 only
    + _Set8KO_H4K20me1.bam.bai_: The index of the above
    + _mm9_custom_chr12.txt_: Chromosome 12 mouse genes in a BED-like format
    + _design.txt_: A design file for profile plot faceting (or separation of
    heatmap profiles) with categorization of genes according to their expression
    (high, medium, low) and their direction of transcrpiption (strand).
* For the RNA-Seq example
    + _WT.bam_: Gene expression RNA signals from wild type adult mice, 
    chromosome 12 only
    + _WT.bam.bai_: The index of the above
    + _Set8KO.bam_: Gene expression RNA signals from adult mice where the Set8 
    (Pr-Set7) gene is knocked out, chromosome 12 only
    + _Set8KO.bam.bai_: The index of the above
    + _design_mm9_rna.txt_: A design file for profile plot faceting (or 
    separation of heatmap profiles) with categorization of genes according to 
    their Ensembl biotype and their direction of transcrpiption (strand). This
    type of categorization is also present in the fixed `recoup` annotation
    database
    
In order to run the vignette examples, you should download and extract the 
archive to a path of your preference, e.g. `/home/me/recoup_tutorial`. In the 
rest of this tutorial, we assume that the path where the test data are placed is
`/home/me/recoup_tutorial`.

## 3. Building a local annotation store

Apart from a user specified file, the reference genomic regions used by recoup 
to construct average profiles over, can be predefined gene set from a few common
organisms supported by recoup. See the `recoup` man page for a list of these
organisms. In order to use this "database" of predefined genomic areas, you 
should run the function `buildAnnotationStore` with a list of organisms, a list
of annotation sources (Ensembl, RefSeq and UCSC supported) and a desired path to
store the annotations (defaults to `/home/me/.recoup`). For example:

```
buildAnnotationStore(c("hg19","mm9","rn4"),c("ensembl","refseq"))
```

See the man page of `buildAnnotationStore` for more details. This step is not
necessary for recoup to run as these annotations can be also downloaded on the
fly. However, if subsets of the supported organisms are to be used often, it is
much more preferrable to spend some time building the local store as it can save
a lot of running time.

## 4. Running recoup

The `recoup` function can be used to create coverage profiles from ChIP-Seq 
like experiments (signals over continuous genomic regions) or from RNA-Seq 
experiments (signals over non-continuous genomic regions). More details 
regarding each type can be found in

* The ChIP-Seq coverage profile vignette (see the package vignettes)
* The RNA-Seq coverage profile vignette (see the package vignettes)

## Quick examples

```
library(recoup)
test.path <- "/home/me/recoup_tutorial/chipseq"

chip.input <- list(
    WT_H4K20me1=list(
        id="WT_H4K20me1",
        name="WT H4K20me1",
        file=file.path(test.path,"WT_H4K20me1.bam"),
        format="bam",
        color="#EE0000"
    ),
    Set8KO_H4K20me1=list(
        id="Set8KO_H4K20me1",
        name="Set8KO H4K20me1",
        file=file.path(test.path,"Set8KO_H4K20me1.bam"),
        format="bam",
        color="#00BB00"
    )
)
```
### ChIP-Seq TSS profiles

```
genome <- file.path(test.path,"mm9_custom_chr12.txt")

test <- recoup(
    input=chip.input,
    region="tss",
    type="chipseq",
    genome=genome,
    flank=c(2000,2000),
    binParams=list(flankBinSize=0,regionBinSize=100),
    selector=NULL,
    plotParams=list(heatmapScale="common"),
    saveParams=list(ranges=TRUE,coverage=TRUE,profile=TRUE),
    rc=0.5
)
```

### RNA-Seq gene body profiles

```
library(recoup)
test.path <- "/home/me/recoup_tutorial/rnaseq"

rna.input <- list(
    list(
        id="WT",
        name="WT",
        file=file.path(test.path,"WT.bam"),
        format="bam"
    ),
    list(
        id="Set8KO",
        name="Set8KO",
        file=file.path(test.path,"Set8KO.bam"),
        format="bam"
    )
)

test <- recoup(
    input=rna.input,
    type="rnaseq",
    genome="mm9",
    flank=c(1000,1000),
    binParams=list(flankBinSize=50,regionBinSize=100),
    selector=NULL,
    preprocessParams=list(normalize="linear"),
    plotParams=list(signalScale="log2"),
    rc=0.5
)
```

