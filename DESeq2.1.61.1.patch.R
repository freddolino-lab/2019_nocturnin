lfcShrink_patched <- function(dds, coef, contrast, res, type="normal") {
  if (is.null(dispersions(dds))) {
    stop("lfcShrink requires dispersion estimates, first call estimateDispersions()")
  }

  # match the shrinkage type
  type <- match.arg(type, choices=c("normal"))

  # fit MLE coefficients... TODO skip this step
  dds <- estimateMLEForBetaPriorVar(dds)

  stopifnot(missing(coef) | missing(contrast))
  if (missing(contrast)) {
    modelMatrixType <- "standard"
  } else {
    modelMatrixType <- "expanded"
  }
  attr(dds,"modelMatrixType") <- modelMatrixType
  betaPriorVar <- estimateBetaPriorVar(dds)

  dds.shr <- nbinomWaldTest(dds,
                            betaPrior=TRUE,
                            betaPriorVar=betaPriorVar,
                            modelMatrixType=modelMatrixType,
                            quiet=TRUE)

  if (missing(contrast)) {
    rn <- resultsNames(dds.shr)
    res.shr <- results(dds.shr, name=rn[coef])
  } else {
    res.shr <- results(dds.shr, contrast=contrast)
  }
  
  if (!missing(res)) {
    ##############################
    # This is where the patch is #
    ##############################
    # original: 
    # res <- res[,c("baseMean","log2FoldChange","stat","pvalue","padj")]
    # edited: 
    
    res <- res[,c("baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
    res$log2FoldChange <- res.shr$log2FoldChange
    res$lfcSE <- res.shr$lfcSE
    # res$stat <- res.shr$stat
    mcols(res)$description[2] <- mcols(res.shr)$description[2]
    return(res)
  } else {
    return(res.shr$log2FoldChange)
  }
}
