if (("xtable" %in% rownames(installed.packages()) == FALSE) || (!require("DESeq2", quietly = TRUE))) {
    if (!require("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
        }
    BiocManager::install("DESeq2",update=TRUE, ask=FALSE, force=TRUE)
}
