if (("limma" %in% rownames(installed.packages()) == FALSE) || (!require("limma", quietly = TRUE))) {
    if (!require("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
        }
    BiocManager::install("limma",update=TRUE, ask=FALSE, force=TRUE)
}
