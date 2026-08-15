options(timeout = max(300, getOption("timeout")))

if (("Rsubread" %in% rownames(installed.packages()) == FALSE) || (!require("Rsubread", quietly = TRUE))) {
    options(repos = c(CRAN="https://cloud.r-project.org/"))
    if (!require("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
    }
    BiocManager::install("Rsubread",update=TRUE, ask=FALSE, force=TRUE)
    if (!requireNamespace("Rsubread", quietly = TRUE)) {
        stop("Failed to install Rsubread. The package is still unavailable.")
    }
}
