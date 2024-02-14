if (("DESeq2" %in% rownames(installed.packages()) == FALSE) || (!require("DESeq2", quietly = TRUE))) {
    options(repos = c(CRAN="https://cloud.r-project.org/"))
    install.packages("png")
    pgks <- list("XML", "vctrs", "RCurl")
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
    if (!require("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
        }
    BiocManager::install("DelayedArray",update=TRUE, ask=FALSE, force=TRUE)
    BiocManager::install("RSQLite",update=TRUE, ask=FALSE, force=TRUE)
    BiocManager::install("DESeq2",update=TRUE, ask=FALSE, force=TRUE)
}
