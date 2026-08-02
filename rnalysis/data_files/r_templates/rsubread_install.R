if (("Rsubread" %in% rownames(installed.packages()) == FALSE) || (!require("Rsubread", quietly = TRUE))) {
    options(repos = c(CRAN="https://cloud.r-project.org/"))
    # Bump R's default download timeout (60s) so slow/degraded connections to CRAN/Bioconductor
    # mirrors have more room to finish the dependency download below.
    options(timeout = max(300, getOption("timeout")))
    if (!require("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
        }
    # See deseq2_install.R for why this retries: a transient network/SSL blip on one of Rsubread's
    # dependencies leaves BiocManager::install() non-fatally warning instead of erroring, so
    # Rsubread can silently fail to load later. Retry a few times, verifying Rsubread actually
    # loads, before giving up (mirrors the dl() retry wrapper used for the kallisto/bowtie2
    # downloads in build_ci.yml).
    rsubread_install_attempts <- 5
    for (attempt in seq_len(rsubread_install_attempts)) {
        BiocManager::install("Rsubread",update=TRUE, ask=FALSE, force=TRUE)
        if (require("Rsubread", quietly = TRUE)) {
            break
        }
        if (attempt < rsubread_install_attempts) {
            message(sprintf("Rsubread failed to install/load (attempt %d/%d); retrying in 15s...",
                             attempt, rsubread_install_attempts))
            Sys.sleep(15)
        } else {
            stop(sprintf("Failed to install Rsubread after %d attempts (see download errors above).",
                          rsubread_install_attempts))
        }
    }
}
