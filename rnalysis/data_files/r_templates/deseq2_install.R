if (("DESeq2" %in% rownames(installed.packages()) == FALSE) || (!require("DESeq2", quietly = TRUE))) {
    options(repos = c(CRAN="https://cloud.r-project.org/"))
    # Bump R's default download timeout (60s) so slow/degraded connections to CRAN/Bioconductor
    # mirrors have more room to finish the many-dependency DESeq2 download below.
    options(timeout = max(300, getOption("timeout")))
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
    # A single transient network/SSL blip while downloading one of DESeq2's ~30 Bioconductor/CRAN
    # dependencies (e.g. Biobase) does not make BiocManager::install() raise an error -- it just
    # warns and leaves that one dependency (and anything that needs it, like SummarizedExperiment)
    # uninstalled, so DESeq2 silently fails to load later with a confusing "could not find function"
    # error. Retry the install a few times, verifying DESeq2 actually loads afterwards, before giving
    # up (mirrors the dl() retry wrapper used for the kallisto/bowtie2 downloads in build_ci.yml).
    deseq2_install_attempts <- 5
    for (attempt in seq_len(deseq2_install_attempts)) {
        BiocManager::install("DelayedArray",update=TRUE, ask=FALSE, force=TRUE)
        BiocManager::install("RSQLite",update=TRUE, ask=FALSE, force=TRUE)
        BiocManager::install("DESeq2",update=TRUE, ask=FALSE, force=TRUE)
        if (require("DESeq2", quietly = TRUE)) {
            break
        }
        if (attempt < deseq2_install_attempts) {
            message(sprintf("DESeq2 failed to install/load (attempt %d/%d); retrying in 15s...",
                             attempt, deseq2_install_attempts))
            Sys.sleep(15)
        } else {
            stop(sprintf("Failed to install DESeq2 after %d attempts (see download errors above).",
                          deseq2_install_attempts))
        }
    }
}
