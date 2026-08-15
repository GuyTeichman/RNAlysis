options(timeout = max(300, getOption("timeout")))

pgks <- list("statmod")
for (pkg in pgks) {
    tryCatch(
        tryCatch(
        {install.packages(pkg)
        require(pkg)},
        warning = function(e) {
        install.packages(pkg, type = "binary")},
        error = function(e) {
        install.packages(pkg, type = "binary")}),
    error = function(e) {}
    )
}
if (("limma" %in% rownames(installed.packages()) == FALSE) || (!require("limma", quietly = TRUE))) {
    options(repos = c(CRAN="https://cloud.r-project.org/"))
    if (!require("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
    }
    BiocManager::install("limma",update=TRUE, ask=FALSE, force=TRUE)
    if (!requireNamespace("limma", quietly = TRUE)) {
        stop("Failed to install limma. The package is still unavailable.")
    }
}
