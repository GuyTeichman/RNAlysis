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
    if (!require("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
        }
    BiocManager::install("limma",update=TRUE, ask=FALSE, force=TRUE)
}

