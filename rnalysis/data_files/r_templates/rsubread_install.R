if (("Rsubread" %in% rownames(installed.packages()) == FALSE) || (!require("Rsubread", quietly = TRUE))) {
    if (!require("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
        }
    BiocManager::install("Rsubread",update=TRUE, ask=FALSE, force=TRUE)
}
