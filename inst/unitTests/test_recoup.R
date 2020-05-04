test_recoup <- function() {
    data("recoup_test_data",package="recoup")

    test.tss <- recoup(
        test.input,
        design=NULL,
        region="tss",
        type="chipseq",
        genome=test.genome,
        flank=c(2000,2000),
        selector=NULL,
        rc=0.1
    )
    
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
        rc=0.1
    )
    
    checkTrue(is.list(test.tss))
    checkTrue(!is.null(test.tss$data))
    checkTrue(is.list(test.gb))
    checkTrue(!is.null(test.gb$data))
}
