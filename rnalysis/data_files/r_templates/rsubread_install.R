if (("Rsubread" %in% rownames(installed.packages()) == FALSE) || (!require("Rsubread", quietly = TRUE))) {
    options(repos = c(CRAN="https://cloud.r-project.org/"))
    if (!require("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
        }
    BiocManager::install("Rsubread",update=TRUE, ask=FALSE, force=TRUE)
}
