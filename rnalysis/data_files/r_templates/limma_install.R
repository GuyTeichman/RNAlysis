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
    # Bump R's default download timeout (60s) so slow/degraded connections to CRAN/Bioconductor
    # mirrors have more room to finish the dependency download below.
    options(timeout = max(300, getOption("timeout")))
    if (!require("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
        }
    # See deseq2_install.R for why this retries: a transient network/SSL blip on one of limma's
    # dependencies leaves BiocManager::install() non-fatally warning instead of erroring, so limma
    # can silently fail to load later. Retry a few times, verifying limma actually loads, before
    # giving up (mirrors the dl() retry wrapper used for the kallisto/bowtie2 downloads in
    # build_ci.yml).
    limma_install_attempts <- 5
    for (attempt in seq_len(limma_install_attempts)) {
        BiocManager::install("limma",update=TRUE, ask=FALSE, force=TRUE)
        if (require("limma", quietly = TRUE)) {
            break
        }
        if (attempt < limma_install_attempts) {
            message(sprintf("limma failed to install/load (attempt %d/%d); retrying in 15s...",
                             attempt, limma_install_attempts))
            Sys.sleep(15)
        } else {
            stop(sprintf("Failed to install limma after %d attempts (see download errors above).",
                          limma_install_attempts))
        }
    }
}

