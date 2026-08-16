=======
History
=======

4.3.0 (unreleased)
-------------------

Added
******
* You can now filter a table by **any** attribute in a GTF/GFF annotation file, not just biotype (``Filter.filter_by_gtf_attribute``) — for example, keep only the genes on a particular chromosome or strand, from a particular annotation source, or of a particular biotype. The reserved attribute names ``chromosome``, ``source`` and ``strand`` read the corresponding columns of the annotation file, while any other name is looked up as a standard attribute (such as ``gene_biotype`` or ``gene_name``). Works with both GTF and GFF3 files.
* You can now annotate a table with a feature attribute drawn from a GTF/GFF annotation file (``Filter.annotate_from_gtf``), adding a new column that labels each gene (or transcript) with, for example, its biotype, chromosome, strand, or source. Features that are absent from the annotation file are left blank. Works with both GTF and GFF3 files.
* The Principal Component Analysis functions (``pca``, ``sort_by_principal_component``, ``split_by_principal_components``) and every clustering function (``split_kmeans``, ``split_kmedoids``, ``split_hierarchical``, ``split_hdbscan``, ``split_clicom``) now let you **choose which transform** is applied to your data before it is standardized and analyzed, instead of only switching the Box-Cox power transform on and off. The new ``'log'`` option applies log2(x+1) to every value — the transform most biologists already use for count data. It is simpler than Box-Cox, applies one identical formula to every gene, and can never become numerically unstable, which makes it the transform to reach for when Box-Cox reports that it cannot handle your table (see Fixed, below).

Changed
*******
* The ``power_transform`` parameter of the PCA and clustering functions is now a drop-down menu of named transforms (``'box-cox'``, ``'log'``, ``'none'``) instead of a True/False switch. ``'box-cox'`` is the default and is exactly what ``True`` always did, and ``'none'`` is exactly what ``False`` always did, so **results are unchanged**. The old ``True``/``False`` values remain valid indefinitely and are mapped to the corresponding transforms, so scripts, saved Pipelines, and exported parameter files written by earlier versions keep running and keep producing identical results.
* Renamed several vague GUI button labels to clearer, action-oriented text: multi-value input buttons ("Set input") now read "Choose values", the load-screen file-picker button now reads "Choose table", the load-screen confirm button now reads "Load", and the Pipeline-creation confirm button now reads "Create Pipeline". This is a labeling change only; behavior is unchanged.
* The single main action of each GUI screen (for example "Load", "Apply", "Create Pipeline", "Run", and external tools' "Start <function>") is now visually distinguished as a filled, accented primary button, while preparatory/secondary buttons keep the default flat style. This is a visual change only; behavior is unchanged.
* Recolored the GUI's boolean toggle switch: the off/False state is now a neutral inactive gray instead of a red-orange that read as an error, and the on/True state is now a warm gold instead of green. Because the off state is a neutral gray, the gold/gray pairing stays off the red-green confusion axis, and the switch still encodes state redundantly through knob position and a "True"/"False" text label, so it remains distinguishable without relying on color.
* Sped up GO and KEGG enrichment analysis. P-value computation is now memoized across the many ontology terms that share the same gene counts (typically over 90% of terms in a GO run), and the ``elim`` propagation method no longer deep-copies the full annotation set on every run. Enrichment results are unchanged; the hypergeometric test and ``elim`` propagation benefit the most (up to ~20× and ~2-9× faster respectively on large ontologies).
* Reduced the time GO enrichment spends fetching and preparing annotations: when the fetched annotations come from several databases, the per-source gene-ID translation (one UniProt round-trip per source) now runs concurrently across sources instead of one after another, and the Gene Ontology graph is downloaded in parallel with organism/gene-ID-type detection during setup rather than before it. Enrichment results are unchanged.
* KEGG enrichment now fetches its pathway-to-gene annotations with a single compact request (KEGG's ``link`` endpoint) instead of downloading a full flat-file record for every pathway one chunk at a time. On *C. elegans* this cut the annotation-fetch step from roughly 63 seconds to about 3 seconds. The gene set of every pathway that was already included is unchanged; see the Fixed section for an accompanying correction to *which* pathways are included.
* Updated the required ``polars`` version to 1.43.x (from 1.41.x). Analysis results are unchanged; this was verified against the RNAlysis test suite.
* The minimum required version of ``aiohttp`` is now 3.12 (from 3.11), whose HTTP client resolves hosts faster and lowers request latency on Python 3.12 and newer — this benefits every RNAlysis feature that talks to a web service (gene-ID translation, ortholog mapping, GO/KEGG enrichment). Analysis results are unchanged.
* The minimum required version of the optional ``numba`` package (used by the randomization tests) is now 0.63 (from 0.61), which keeps the ``randomization`` extra installable alongside current NumPy versions and adds Python 3.14 support. Results are unchanged.
* CLICOM ensemble clustering (``CountFilter.split_clicom``) now runs faster, by computing each distinct power-transform/standardization of the data only once per run and reusing it across the many clustering setups that share it, instead of recomputing the identical transform for every setup. Clustering results are bit-for-bit identical.
* CLICOM's clique-finding step (the second-slowest part of ``CountFilter.split_clicom``) now runs up to tens of times faster, by replacing its cubic pure-Python loop over sets with a vectorized NumPy bitset implementation. Clustering results are bit-for-bit identical.
* The Principal Component Analysis functions (``pca``, ``sort_by_principal_component``, and ``split_by_principal_components``) now run substantially faster on large datasets, by parallelizing the per-gene Box-Cox power transform across CPU cores. Results are unchanged up to floating-point precision (verified against the test suite's reference outputs).
* Ensembl ortholog and paralog mapping now caches each gene's homology results in the daily cache, so repeating a mapping (or mapping gene sets that overlap a previous one) no longer re-requests the same data from Ensembl. This reduces load on the Ensembl REST API and speeds up repeated mappings.
* Mapping genes to orthologs through **PhylomeDB** now runs substantially faster. Each mapping previously loaded PhylomeDB's database-wide ID-conversion table (which spans every species, millions of rows) into two whole-table Python dictionaries on every call, even though only the queried genes are looked up; the ID-conversion and per-organism ortholog tables are now filtered to just the required identifiers before being materialized. Mapping results are unchanged.
* Mapping genes to orthologs through **OrthoInspector** in automatic-database mode now runs much faster. Automatic mode tries OrthoInspector's databases largest-first, but its largest database (``Eukaryota2023``) currently accepts the connection and then never answers ortholog queries — so every automatic mapping wasted a full read-timeout stalling on it before falling back to the next database. That known-stalling database is now skipped during automatic selection, and the per-database species lists used to pick a database are fetched concurrently instead of one at a time. On a typical worm→briggsae mapping this cut the cold-cache time from ~45s to ~10s. Mapping results are unchanged (a database that never responds never contributed a result); if OrthoInspector restores that database, remove it from the skip-list in ``io.py``.
* RNAlysis can now read **GFF3** genome-annotation files (in addition to GTF) everywhere it uses annotations — filtering and summarizing features by biotype, RPKM/TPM normalization, and summing transcript-level counts to genes. It also accepts the ``.gff`` extension and **gzip-compressed** (``.gz``) annotation files, and detects whether a file is GTF or GFF3 from its contents rather than relying solely on the file extension. For GFF3 files, the SO-term prefixes that some sources (e.g. Ensembl) add to identifiers — such as ``gene:`` and ``transcript:`` — are stripped so the resulting gene/transcript IDs match those in your count table.
* When a GTF/GFF annotation file contains malformed lines (not the expected 9 tab-separated columns), RNAlysis now consistently skips them and emits a single warning reporting how many lines were skipped. Previously the different annotation-based features handled such lines inconsistently — some skipped them silently, while summing transcript counts to genes raised an error and aborted.
* When filtering or summarizing features by biotype from a GTF/GFF file, the biotype attribute name now automatically resolves across common naming conventions (for example Ensembl's ``gene_biotype`` versus GENCODE's ``gene_type``): the default option finds the right attribute regardless of which convention your file uses, while an explicitly chosen attribute name still takes precedence.
* Gene-ID translation through UniProt (used by ``map_gene_ids`` and by GO/KEGG enrichment analysis) now runs faster: the ID-mapping job status is polled with an adaptive backoff that starts almost immediately instead of waiting a fixed few seconds between polls, and results for UniProtKB-target mappings are retrieved in a single streaming request instead of walking the result pages one at a time. Translation results are unchanged, with one intentional improvement to reproducibility: when a gene maps to several targets that share the *same* highest UniProt annotation score, the winning target is now chosen deterministically (by target ID) instead of depending on the order the results happened to arrive, so translation — and any enrichment that relies on it — is now fully reproducible. Verified against the test suite, including live UniProt round-trips.
* Gene-ID translation through UniProt (used by ``map_gene_ids`` and by GO and KEGG enrichment analysis) now caches each gene's mapping in the daily cache, so re-running an analysis on the same or an overlapping set of genes only sends the genes that haven't been translated yet to UniProt, instead of re-translating the whole set. This speeds up the common workflow of repeatedly re-filtering or re-clustering and re-running enrichment on overlapping gene sets, and reduces load on the UniProt service. Translation results are unchanged.
* Several ``filtering.CountFilter``/``Filter`` operations now run in a single fused Polars pass instead of scanning the table multiple times: the count-based filters (``filter_low_reads``, ``split_by_reads``, ``filter_by_row_sum``), all per-sample and per-gene normalizations, ``fold_change``, replicate averaging (``average_replicate_samples``), and the ``opposite=True`` path of every filter. Results are bit-for-bit identical, verified by a new eager-vs-lazy equivalence test suite.
* On the "New table" screen, RNAlysis now auto-detects the most likely table type from the chosen file and pre-selects it in the table-type menu (count matrix, differential expression, or fold change), instead of always defaulting to "Other table". The detection is conservative — when the content is ambiguous it falls back to "Other table" rather than guess wrong — and it only sets the default; you can still change the table type manually. This only affects the pre-selected type, not how tables are actually loaded.
* In the graphical interface, numeric input boxes for decimal (``float``) parameters now adapt to the magnitude of the value: an integer-valued default such as a read-count threshold of ``5`` is shown as ``5`` (rather than ``5.00``) and steps by whole units, while small fractional parameters keep a finer step and more decimals. This is a display-only change — the exact value you enter is still passed to the analysis unchanged, and arbitrary precise values can still be typed.
* Starting RNAlysis, and importing the Python API, is now faster and no longer requires an internet connection. The lists of supported organisms and gene-ID types that fill the drop-down menus of the ortholog-mapping, paralog-finding and gene-ID translation functions were previously fetched from UniProtKB, PantherDB, Ensembl and PhylomeDB every time RNAlysis started — four live web requests on the first startup of each day, and close to a second of waiting even when they were already cached. They are now read from a list packaged with RNAlysis and refreshed with every release. On an offline or firewalled computer (a compute cluster, a locked-down lab PC) these menus are now fully populated instead of empty, and the available values are identical on every machine running the same version of RNAlysis, instead of depending on the day the services were contacted. Actual analyses — gene-ID translation, ortholog and paralog mapping, and enrichment analysis — still query these services live, exactly as before, so results are unchanged. An organism or ID type that a service added after the release can still be typed into these fields by hand, as before.
* *RNAlysis* now starts up faster. The heaviest third-party packages (seaborn, scikit-learn, pandas, statsmodels, matplotlib-venn, UpSetPlot, kmedoids and hdbscan) are no longer imported when you import *RNAlysis* — they are loaded on demand, the first time a function that actually needs them runs. Importing ``rnalysis.filtering`` or ``rnalysis.enrichment`` is about 12-20% faster, and every script, notebook, and parallel worker process pays the smaller cost. Analysis results are unchanged.
* The randomization tests (``FoldChangeFilter.randomization_test`` and the randomization enrichment test) now cache their compiled ``numba`` code on disk, so only the very first run after installing or upgrading pays the compilation cost, instead of every new Python process (including every parallel worker process) paying it again. On this machine the first randomization call in a fresh process dropped from 0.83s to 0.15s. Results are unchanged. Caching is switched off in the standalone (frozen) app, where ``numba``'s cache location cannot yet be verified.
* The minimum required version of ``scikit-learn`` is now 1.6.0 (from 1.5.0) — the first version that can load its submodules on demand, which is what lets *RNAlysis* defer importing it.
* *RNAlysis* now resolves its settings and cache directories with the actively-maintained ``platformdirs`` package instead of the unmaintained ``appdirs`` (archived, last released in 2020). On every supported operating system the directories resolve to the same locations as before, so existing settings and caches are picked up unchanged. The one exception is macOS machines where the ``XDG_DATA_HOME``/``XDG_CACHE_HOME`` environment variables are set: there ``platformdirs`` honors them (``appdirs`` did not), so *RNAlysis* automatically moves its existing data folder — including your settings — to the new location the first time it runs, and starts a fresh daily web cache there.
* *RNAlysis* no longer uses bare ``assert`` statements to validate your input. Every input check now raises a real exception, so validation can no longer be silently removed: previously, running Python with the ``-O`` flag (or with ``PYTHONOPTIMIZE`` set in your environment) stripped all 341 of them, and invalid inputs flowed straight into the analysis, either failing later with a cryptic internal error or — worse — producing a wrong result without complaint. This affected the API/script workflow only; the standalone app was never affected.
* The exceptions *RNAlysis* raises for invalid input are now typed, so scripts can tell the different failure modes apart: a wrong-type argument raises ``rnalysis.exceptions.InvalidTypeError`` (which is also a ``TypeError``), a bad value — out of range, not one of the legal choices, wrong shape — raises ``rnalysis.exceptions.InvalidValueError`` (also a ``ValueError``), and a violated internal invariant (which always indicates a bug in *RNAlysis*) raises ``rnalysis.exceptions.InternalError``. All three share the common base ``rnalysis.exceptions.RNAlysisError``, and the first two additionally share ``RNAlysisInputError``. **If your scripts catch ``AssertionError`` around *RNAlysis* calls, change them to catch ``RNAlysisInputError`` (or the built-in ``TypeError``/``ValueError``), which is what those calls now raise.** Which inputs are accepted or rejected is unchanged, and every input-error message is word-for-word what it was before; internal-invariant messages additionally gain a standardized line pointing you at the issue tracker. Analysis results are unaffected.

Changed
*******
* Exporting an analysis report now streams each table straight from parquet to CSV instead of loading the whole table into memory first, reducing peak memory use when exporting reports with large tables.

Fixed
******
* Fixed a confusing crash (``ValueError: Input X contains NaN``, on a table containing no missing values at all) in the Principal Component Analysis functions and in the power-transform clustering functions. The Box-Cox power transform fits a separate exponent to every gene; on a gene whose values are near-constant at a high magnitude and measured across only a handful of samples — a reporter transgene such as GFP, for example — that exponent explodes and the transform overflows. *RNAlysis* now detects this (and the related case of a transform that stays finite but returns absurdly large values, which used to distort the analysis without crashing) and stops with a clear message **naming the offending genes**, so you can either filter them out or switch to the new ``'log'`` transform, which is immune to the problem. Tables that already worked are entirely unaffected, and their results are unchanged.
* Fixed a bug in the graphical interface where importing a parameter file into a function window could silently ignore a value that is not one of the listed choices: the drop-down kept displaying (and the analysis kept using) the previously-selected choice instead of switching to the imported value.
* Fixed a bug where reading gene/transcript attributes from a GTF file (used when filtering or summarizing features by biotype) with the "use version" option enabled raised an error on annotation files that have no version fields — such as GENCODE files, whose versions are embedded directly in the IDs. The version suffix is now simply omitted in that case instead of failing.
* Fixed a bug where mapping transcripts to an attribute *by feature name* could create a spurious, empty-named entry for a transcript that had a gene name but no transcript name; such transcripts are now correctly skipped.
* Fixed a bug where averaging replicate samples (``CountFilter.average_replicate_samples``) crashed when the averaging function was set to ``'median'`` or ``'geometric_mean'`` (both selectable in the GUI), because the underlying Polars calls no longer existed in current Polars versions. Both options now work correctly, computing the row-wise median / geometric mean across each group of replicate columns.
* Restored compatibility with current Polars versions in two internal code paths that relied on since-removed Polars behavior: loading a table while skipping comment lines, and parsing annotation scores during gene-ID translation.
* Fixed a bug in ``CountFilter.pairplot`` where the Spearman correlation shown for the first row/column of the plot was computed against the gene-index column instead of a sample, producing an incorrect value (and, with recent NumPy versions, an error). The correlations are now always computed between the correct pair of samples.
* Fixed a bug where mapping orthologs through the OrthoInspector service could hang indefinitely when one of its databases stopped responding (OrthoInspector relocated its API, and its largest database can stall). OrthoInspector requests now use a timeout, target OrthoInspector's current API host directly, and automatically fall back to the next available database.
* Fixed a bug where gene-ID translation and enrichment analyses (which rely on UniProt) could intermittently fail with an ``HTTP 400`` error while waiting for a UniProt ID-mapping job, especially when several analyses ran at once. Both steps that query a running job — checking its status and fetching its result link — now retry such transient errors a bounded number of times with jittered backoff, instead of failing on the first one.
* Fixed a bug where gene-ID translation and GO/KEGG enrichment (which rely on UniProt's asynchronous ID-mapping service) could hang indefinitely if a UniProt mapping job never became ready — for example when UniProt wedged a job under heavy load. RNAlysis now stops waiting for a single job after a bounded time, warns, and continues with a partial mapping: the genes it could not translate this run are retried on the next run rather than being cached as unmappable, and any genes that *were* translated are unaffected. The bound is deliberately set well above how long a legitimately slow job takes, so it never cuts short a job that would have finished, and is measured against wall-clock time so it also covers slow network requests, not just the wait between status checks. Every request this feature makes to UniProt now also has its own connection/read timeout, so a request that hangs at the network level (no response ever sent) is bounded too, instead of blocking forever before RNAlysis even gets a chance to apply the bound above. A UniProt job that fails outright (rather than merely running long) now also degrades the same way instead of crashing the analysis. As part of this fix, the per-poll message shown while waiting for a slow UniProt job was reworded from the misleading ``Retrying in Ns`` — which read like a failed request being retried — to make clear that the request succeeded and RNAlysis is simply waiting for UniProt to finish the job, and the degradation warning now correctly distinguishes "no results" from "partial results" instead of always claiming a partial result.
* Fixed a bug where closing the RNAlysis window could occasionally crash the application, because its background worker and console threads were signalled to stop but not waited for before the window was destroyed. Closing now shuts those threads down cleanly first.
* Hardened the gene-set visualization and set-operation windows against a rare crash that could occur while their preview plot was being refreshed: the previous plot is now fully detached from the window before it is scheduled for deletion, so a pending redraw can no longer target a plot that is being torn down.
* Fixed a bug where enrichment analyses (GO, KEGG, user-defined, and single-set) could crash with a confusing ``not enough values to unpack`` error when no results could be produced — either because the enrichment gene set ended up empty (for example when a web service temporarily mapped none of the genes to the required ID type, or when every gene was filtered out), or because the chosen statistical test could not be run (for example when the optional ``xlmhglite`` package is missing for single-set analysis). Such analyses now finish gracefully with an empty result instead of crashing.
* Fixed a bug where running KEGG enrichment for two different organisms on the same day could make the second analysis use the **first** organism's pathways. The per-organism KEGG pathway list was cached under an organism-agnostic filename, so the second organism loaded the first's cached list instead of fetching its own. The pathway-list cache is now organism-specific. (The compound/glycan caches are organism-agnostic by nature and were unaffected.)
* Fixed a bug where KEGG enrichment silently omitted every pathway whose KEGG record has no ``COMPOUND`` section — including important non-metabolic pathways such as Ribosome, Spliceosome, Proteasome, RNA polymerase, and several glycan-metabolism pathways — because the flat-file parser stopped collecting genes at the (absent) ``COMPOUND`` line and so read none. KEGG annotations are now read from KEGG's ``link`` endpoint, which lists all pathways that have genes, so these pathways are included: on *C. elegans* this recovered 36 previously-missing pathways. KEGG's global/overview maps (e.g. "Metabolic pathways", "Carbon metabolism") remain excluded, as they were before, because these ~thousand-gene superset maps are not biologically specific enrichment terms. **Because this changes which pathways are tested, KEGG enrichment results — and their multiple-testing correction — may differ from previous versions.** The gene sets of all pathways that were already being tested are unchanged (verified against the previous implementation on live KEGG data).
* Fixed a bug where automatically detecting the organism or gene-ID type from a set of gene IDs (via Ensembl) failed with an ``HTTP 400`` error, because the request body sent to Ensembl was double-encoded. The request is now formatted correctly.
* Fixed a bug where a malformed or unreachable external service (PantherDB, UniProt, Ensembl, or PhylomeDB) could crash RNAlysis on startup while it loaded the lists of supported organisms and gene-ID types. These lookups now degrade gracefully to an empty list with a warning, and their network requests use timeouts so a stalled service can't freeze startup.
* GUI table caching now runs off the UI thread and is guaranteed to finish before a session is saved, a report is exported, or the cache is cleared. Previously the asynchronous write introduced in 4.2.0 was "fire-and-forget", which could race with those operations and produce incomplete cached files.
* Fixed a bug where automatic common-name generation (used by the "smart" paired-end sample-naming option in the FASTQ functions) could produce an incorrect name when input file pairs differed in length, because a filtering step reused a loop variable left over from the last pair instead of each pair's own values.
* Fixed a crash (``IndexError``) in automatic common-name generation for inputs of uneven length.
* Fixed a bug where mapping orthologs or paralogs through the PantherDB service could crash the whole analysis when PantherDB intermittently returned an empty response for a gene (an empty ``HTTP 200`` body, which its retry mechanism does not cover). Such requests are now retried, and if the response stays empty that single gene is skipped with a warning instead of aborting the entire mapping.
* Fixed a bug where the ``random_seed`` parameter had no effect on the randomization (permutation) enrichment test, so its p-values were not reproducible even when a seed was set. The seed was applied to NumPy's random generator, but the permutation sampling runs under numba's separate generator, which was never seeded; the seed is now applied where the sampling actually happens. As a result, a given ``random_seed`` now reproduces the same randomization p-values exactly — and identically across all parallel backends (sequential, multiprocessing, loky, threading). Randomization p-values for a given seed are now deterministic; because they were previously non-deterministic, no earlier analysis could have relied on specific values.
* Fixed a crash on startup that prevented the application from opening for users whose ``settings.yaml`` held a stored Attribute or Biotype Reference Table path of ``null`` (for example after a previously-set reference table had been cleared). The Settings window, which is built during startup, fed that ``None`` value into a file-path field and raised a ``TypeError``. Such a value is now treated the same as "No file chosen", and the underlying path-legality helpers no longer error when given ``None``.
* Fixed a bug where choosing an invalid enrichment bar-plot style crashed with a confusing ``AttributeError`` instead of telling you which values are allowed — the error message tried to quote the plot style from an attribute that had not been set yet. An invalid plot style now reports the value you passed and the legal options.
* Fixed a bug where a failed installation of an R package (DESeq2, limma, or Rsubread) showed only R's raw exit status instead of *RNAlysis*'s guidance. The most common cause on fresh machines and computing clusters — no write permission to R's library folder — now surfaces that explanation along with a suggestion to install the package manually, while R's own output remains visible as the underlying cause. The guidance existed but was attached to the wrong error and could never actually be reached.
* Failed Bioconductor installs now stop immediately if DESeq2, limma, or Rsubread is still unavailable after installation, instead of exiting successfully and surfacing a misleading missing-function error during the later analysis. The download timeout for these installs is also raised from R's 60-second default to at least 300 seconds.
* **Saving a session over an existing session file can no longer destroy it.** Until now, saving assembled the new session directly on top of the old one — the existing ``.rnal`` file was deleted *before* anything replaced it, and the replacement archive was only produced at the very end. If the save failed partway through (a full disk, a permissions problem, a crash, or closing the laptop lid), the previous session was already gone and the new one had never been written, losing the entire analysis. The new session is now assembled under a temporary name alongside the target file and moved into place in a single atomic step, so your previous session file stays untouched until its replacement is complete on disk; if a save fails, you keep the previous session and RNAlysis cleans up after itself. Saving also now recovers gracefully if an earlier crash left a half-built session *folder* where the session file should be.
* Fixed a bug where loading a session modified the session file itself: RNAlysis renamed your ``.rnal`` file to ``.rnal.zip`` while reading it and renamed it back afterwards. Loading a read-only session file (for example one received as an email attachment or stored on write-protected media), or one that was momentarily locked by antivirus software or a cloud-sync client, could therefore fail outright — and a crash at exactly the wrong moment could leave your session renamed to ``.rnal.zip``. Sessions are now read without renaming or otherwise touching the file.
* Fixed a bug where opening a corrupt, truncated, or partially-downloaded session file reported an *internal RNAlysis bug* and asked you to file a bug report, instead of explaining what was actually wrong. Such a session now produces a clear message stating that the session file is corrupt or incomplete, naming the missing piece and — where it can still be read from the file — the version of RNAlysis that saved the session. A failed load also no longer leaves a half-extracted copy of the session behind in the cache.
* Fixed a bug where exporting a Pipeline crashed with a cryptic ``RepresenterError: cannot represent an object`` if any of its function parameters was a NumPy value (for example a threshold taken straight from a results table), a ``pathlib.Path``, or an ``Enum`` value. Worse, when you were overwriting an existing Pipeline file, that crash happened *after* the file had already been emptied — so the failed export destroyed the Pipeline you were replacing. Such values are now converted to their plain equivalents when the Pipeline is saved; the Pipeline is fully prepared before the target file is opened at all, so a failed export leaves your existing Pipeline file exactly as it was; and a value that genuinely cannot be saved to a Pipeline file now produces a clear message naming the function and the parameter at fault. Note that YAML has no tuple type: a tuple nested inside a parameter is saved — and re-loaded — as a list. This is now documented in ``export_pipeline`` and ``import_pipeline``.
* Fixed a bug where mistyping the name of a Pipeline file made RNAlysis try to interpret the file *path* as the contents of a Pipeline, failing with an unrelated ``TypeError: string indices must be integers``. A path that does not exist now raises a ``FileNotFoundError`` naming the file, while passing the contents of a Pipeline file directly (as a string) keeps working as before.
* Pipelines record the version of RNAlysis that exported them, but nothing ever read that stamp. When a Pipeline cannot be loaded — because it names a function that no longer exists, or because the file is malformed — the error now states which version of RNAlysis exported the Pipeline and which version you are running, instead of failing with a bare ``AttributeError``.
* Fixed a bug where the organism drop-down menu of "Map genes to nearest orthologs (using PhylomeDB)" was always empty, so an organism could only be entered as free text. RNAlysis fetched PhylomeDB's list of supported species correctly, but discarded every entry while filtering that list, so not a single species ever reached the menu. The menu now offers all species supported by PhylomeDB (6,284 of them), sorted alphabetically; entering an organism name or taxon ID by hand still works exactly as before. This changes only which values the menu offers - ortholog mapping results are unchanged.
* Fixed a bug where importing a gene set from a ``.csv`` or ``.tsv`` file kept stray leading/trailing whitespace around the gene identifiers (for example a trailing space after an ID), producing broken identifiers that match nothing in reference or count tables. Identifiers are now trimmed on import — and cells that contain only whitespace are ignored instead of importing as an empty identifier — matching how gene IDs are read when loading tables.
* Fixed a bug where passing a ``NaN`` read-count threshold to ``CountFilter.filter_low_reads``, ``CountFilter.split_by_reads``, or ``CountFilter.filter_by_row_sum`` silently discarded **every** gene and returned an empty table instead of reporting the bad value. Since no number compares as greater than or equal to ``NaN``, the threshold check let it through and the filter then matched nothing. A ``NaN`` threshold is now rejected like any other invalid value. (A ``NaN`` could reach these functions from a spreadsheet cell that failed to parse, or from an upstream calculation that divided by zero.)
* Fixed a bug where importing a saved parameter file that contains covariates or Likelihood Ratio Test factors into the **simplified** differential-expression window reported an internal error and asked you to file a bug report. It now explains that the simplified window does not support those settings, and tells you to open the full version of the window or remove them from the file.
* Fixed a bug where a single OrthoInspector database reporting a non-success status made the whole ortholog-mapping lookup fail with a "please report this bug" message, even though the fault lay with the external service. A database in that state is now skipped with a warning naming it and the status it reported, and the remaining databases stay usable.


4.2.0 (2026-05-30)
-------------------

Changed
*******
* Added Python 3.14 support
* Performance improvements, especially when using automatic report generation.

Fixed
******
* Fixed bug that caused RNAlysis to crash when PantherDB is unavailable
* Adapted RNAlysis to the updated KEGG taxon mapping
* Fixed bug that caused Ensembl ortholog mapping to fail in some environments
* Fixed bug that caused session loading to fail in some environments


4.1.2 (2025-06-04)
-------------------

Changed
*******
* Added Python 3.13 support

Fixed
******
* Fixed bug that caused Windows standalone version to sometimes fail to launch
* Fixed bug that caused RNAlysis to crash when trying to clear app cache without write access

4.1.1 (2025-01-11)
-------------------

Changed
*******
* Updated dependency versions
* Made jdk an optional dependency

Fixed
******
* Fixed bug where kallisto quantification would sometimes fail to sum transcripts to genes properly.
* Fixed bug where automatic mapping of gene IDs would sometimes raise an exception.

4.1.0 (2024-09-16)
-------------------
Version 4.1.0 of RNAlysis features improved performance and stability across the board, thanks to a switch to Qt6 and Numpy 2.


Changed
*******
* RNAlysis now runs on Qt6 instead of Qt5. This change should improve both performance and stability across RNAlysis.
* When importing parameters for differential expression analysis, RNAlysis will now automatically load the design matrix and chosen comparisons. This should work even with parameters that were exported in previous versions of RNAlysis.
* Made small improvements to the RNAlysis graphical interface.
* RNAlysis and its dependencies now run on Numpy 2 instead of Numpy 1.
* RNAlysis now uses a different implementation of the K-Medoids clustering algorithm, which should be more stable and faster than the previous implementation. However, note that the two implementations may give slightly different results.
* When running differential expression or feature counting, RNAlysis session reports will automatically include a logfile with R session info.
* Added optional parameters to all differential expression functions, allowing users to return a path to a logfile with R session info.

Fixed
******
* Fixed bug where table previews in the graphical interface would sometimes fail to update when applying functions inplace.

4.0.0 (2024-06-29)
-------------------
Version 4.0.0 of RNAlysis introduces two major improvements:
First, the software has switched from Pandas to Polars, which should substantially improve overall performance, particularly when generating automatic analysis reports.

Second, significant enhancements have been made to the automatic analysis reports.
These include new interactive features such as the ability to highlight analysis paths and filter the report by node types,
improved hierarchical layout for better readability, performance optimizations, and additional customization options.
The reports now offer a more comprehensive and user-friendly way to navigate and understand analysis workflows.

This update also includes various bug fixes and minor improvements to other functions within the software.

Happy analysis!

Added
******
* Added additional optional parameters to automatic report generation.
* Added additional parameters to 'Hierarchical clustergram plot' (CountFilter.clustergram).

Changed
*******
* RNAlysis now uses Polars instead of Pandas. This change should improve the performance of RNAlysis by an order of magnitude, especially when generating automatic analysis reports.
* Automatic analysis reports are now laid out in a hierarchical structure, making it easier to navigate through the report.
* When clicking on a node in an automatic analysis report, the analysis path leading to that node will be highlighted.
* When clicking on a legend node in an automatic analysis report, all nodes of that type will be highlighted.
* The 'Hierarchical clustergram plot' function should now run faster on large datasets.
* Added support for running tasks using 'sequential' backend on standalone releases of RNAlysis.

Fixed
******
* Fixed bug where the visualization functions 'Plot histogram of p-values' and 'Histogram of a column' would not display a graph immediately when using the stand-alone app.
* Fixed bug re-loading a saved session report would sometimes lead to missing connections in the analysis flow graph.
* Fixed bug where applying not-inplace operations to differential expression tables would change the expected name of the p-value column.
* Fixed bug where table previews in automatic analysis reports would display an extra empty line.

Backwards-incompatible changes
*******************************
* Since RNAlysis now uses Polars instead of Pandas, all functions that previously returned Pandas DataFrames will now return Polars DataFrames. To keep your code compatible with the new version, you may need to adjust the way you interact with the returned DataFrames, or convert them back into Pandas DataFrames.
* Sessions saved with previous versions of RNAlysis are mostly compatible with RNAlysis 4.0.0, but some of the new features of automatic analysis reports (such as highlighting the report by node type) may not work in older sessions.

3.12.0 (2024-05-14)
-------------------
I'm happy to announce the latest release of RNAlysis, which brings a variety of new features and improvements.

One of the key additions in this version is expanded support for advanced differential expression analysis using DESeq2 and Limma-Voom.
You can now test continuous covariates, perform likelihood ratio tests for factors, interactions, and polynomials, and take advantage of sample quality weights in Limma-Voom analysis.
We've also added several new visualization options, such as a p-value histogram plot and more customization choices for gene expression plots.

This release also includes various usability enhancements and bug fixes to provide a smoother analysis experience. We've clarified error messages, improved parameter annotations, and resolved issues that could lead to crashes or incorrect results in certain scenarios. RNAlysis now integrates with the latest versions of the OrthoInspector and ShortStack APIs as well.

As always, your support and feedback are greatly appreciated.
If you have any questions, encounter issues, or have suggestions for future development, please don't hesitate to reach out.
Happy analysis!
Guy

Added
******
* Added support for advanced differential expression analysis with DESeq2/Limma-Voom, including testing continuous covariates, as well as likelihood ratio tests for factors, interactions, and polynomials.
* Added support for sample quality weights in Limma-Voom differential expression analysis.
* Added support for the 'cooksCutoff' parameter in DESq2 differential expression analysis.
* Added support for user-provided scaling factors in DESeq2 differential expression analysis.
* bowtie2 and kallisto paired-end modes now support the new file naming method 'smart', that attempts to automatically determine the common name of each pair of paired-end fastq files.
* Added a p-value histogram plot (DESeqFilter.pval_histogram) that displays the distribution of p-values in a differential expression table.
* Added new visualization options for 'plot expression of specific genes' function (CountFilter.plot_expression).
* Added new function 'Histogram of a column' (Filter.histogram) that generates a histogram of a column in a table.
* Added new function 'Concatenate two tables' (Filter.concatenate) that concatenates two tables along their rows.
* The filtering function 'Filter with a number filter' now supports filtering based on absolute values.

Changed
********
* When running differential expression analysis, RNAlysis will automatically ensure that the order of samples in your design matrix matches the order of samples in your count matrix, avoiding erronious results.
* The 'Filter genes with low expression in all columns' function (CountFilter.filter_low_reads) now supports the 'n_samples' parameter, allowing users to filter genes with a minimal expression threshold in a specific number of samples.
* The 'Plot expressino of specific genes' function (CountFilter.plot_expression) now supports the 'jitter', 'group_names' and 'log_scale' parameters, allowing users to further customize the plot.
* The 'Scatter plot - sample VS sample' function (CountFilter.scatter_sample_vs_sample) now always displays highlighted points on top of the plot, making it easier to see which points are highlighted.
* Kallisto quantification (kallisto_quantify_single_end and kallisto_quantify_paired_end) now supports the 'summation_method' parameter, allowing users to choose between 'raw' and 'scaled_tpm' transcript summation methods. The default behavior of the functions did not change (it corresponds to 'scaled_tpm').
* Enrichment bar plots now have optional parameters that control font sizes for titles and labels.
* Moved the enrichment analysis and enrichment graphs actions to the "Enrichment" menu in the graphical interface to make the actions easier to find.
* Improved the clarity of error messages when attempting to read an invalid GTF file.
* RNAlysis now supports the latest version of OrthoInspector API.
* RNAlysis now supports the latest version of Ensembl Ortholog API.
* Improved annotation for the 'metric' parameter of the 'Hierarchical clustergram plot' function (CountFilter.clustergram).
* Improved performance of RNAlysis when generating automatic session reports.
* RNAlysis now offers default values for differential expression tables' column names.
* Functions that average replicates now display clearer group names by default.
* The RNAlysis interface to ShortStack now uses the most recent API (replaced 'knownRNAs' with 'known_miRNAs').
* When running differential expression, RNAlysis session reports will automatically include the auto-generated R script, as well as the sanitized design matrix used.
* Added optional parameters to all differential expression functions, allowing users to return a path to the auto-generated R script and data sanitized design matrix used.

Fixed
*******
* Fixed bug where some FASTQ/SAM functions could not be added to a FASTQ pipeline.
* Fixed bug where bowtie2 could not be run in/from directories with spaces in their names.
* Fixed bug where RNAlysis would crash when launched without an internet connection.
* Fixed bug that cause ID-mapping functions to raise an error when called from the MacOS stand-alone app (thanks to `Mitchzw <https://github.com/Mitchzw>`_ in `#34 <https://github.com/GuyTeichman/RNAlysis/issues/34>`_).
* Fixed bug that caused R package installations (DESeq2, limma, etc) to fail on some computers (thanks to `Celine-075 <https://github.com/Celine-075>`_ in `#35 <https://github.com/GuyTeichman/RNAlysis/issues/35>`_).
* Fixed bug where Limma-Voom differential expression would fit numerical covariates as categorical factors.
* Fixed bug where generating enrichment bar plots with ylim='auto' would cause bars with 100% depletion (log2FC=-inf) to disappear.
* Fixed bug where defining 10 or more sample groups in the graphical interface would cause the groups to be ordered incorrectly in graphs.
* Fixed bug where the 'return_scaling_factors' argument would not return the normalization scaling factors on the graphical interface.
* Fixed various visual issues in the graphical interface
* Fixed bug where the 'filter_by_row_sum' function would raise an error when the table contains non-numerical columns
* Fixed bug where running enrichment on an empty gene set would raise an error.
* Fixed bug where RNAlysis would suggest resuming an auto-report from loaded session even when auto-report is turned off.
* Fixed bug where disabling auto-report in the middle of the session would raise errors when trying to create new graphs.
* Fixed bug where generating multiple gene expression plots (split_plots=True) with auto-generated report would only add the last graph to the session report.
* Fixed bug where the function 'normalize_to_quantile' generated unclear table names.
* Fixed bug where sometimes the first operation performed in a session would not display correctly in the automatic session report.

New Contributors
*****************
* `Mitchzw`_ in `#34`_
* `Celine-075`_ in `#35`_

3.11.0 (2024-01-05)
-------------------
This release brings several exciting new features.
Notably, these inclue the ability to run functions from the Picardtools suite, and the ability to automatically save interactive session reports and later resume them from any saved session.
In addition, this release includes several visual upgrades, bug fixes, and quality-of-life improvements.
Happy analysis!

Added
*******
* *RNAlysis* can now run functions from the Picardtools suite, including conversion functions (BAM to SAM, SAM to FASTQ, FASTQ to SAM, etc), quality control (validate SAM), and post-processing functions (remove PCR duplicates, sort SAM, create BAM index).
* *RNAlysis* interactive session reports can now be resumed from any saved session, instead of having to start a new report from scratch. When loading a session created by *RNAlysis* 3.11 and beyond, you will have the option to resume the interactive report from the last saved state.

Changed
********
* CutAdapt adapter trimming functions can now receive an optional "new_filenames" parameter, which allows users to specify the names of the output files.
* Hierarchical clustergram plot (CountFilter.clustergram) now supports the 'colormap' parameter, which allows users to specify a custom color map for the plot.
* Hierarchical clustergram plot (CountFilter.clustergram) now displays continuous values on the color bar, instead of discrete values.
* Generally improved the looks of Hierarchical clustergram plot (CountFilter.clustergram).
* Previously-added functions (such as ortholog/paralog mapping, biotype summary by GTF file, etc) can now be applied to gene sets. Previously, some of these functions could only be applied to data tables.

Fixed
*******
* Function "Summarize feature biotypes (based on a reference table)" (biotypes_from_ref_table) now treats rows with missing values as "_missing_from_biotype_reference" instead of ignoring them entirely.
* Fixed bug where the Ensembl paralog-finding function would appear under the wrong tab in the graphical interface.
* Fixed bug where the description of the MA Plot function and parameters would not display correctly in the graphical interface.

3.10.1 (2023-11-22)
-------------------
Version 3.10.1 introduces several bug fixes, as well as well as support for random effect analysis in Limma-Voom differential expression.

Added
*******
* Limma-Voom differential expression can now fit mixed linear models containing a random effect (e.g. nested design).


Fixed
*******
* Fixed bug where trying to load sessions created with RNAlysis version 3.10.0 would result in an error.
* Fixed bug where using the OrthoInspector ortholog mapping function with database='auto' would sometimes fail to find an appropriate mapping database, even when one exists.
* Fixed bug where kallisto paired-end quantification window would display the 'read1' and 'read2' parameters twice.
* Fixed bug where empty sub-menus would appear under the FASTQ menu. These sub-menus will be implemented in future versions of RNAlysis.

3.10.0 (2023-10-31)
-------------------
I'm thrilled to introduce RNAlysis version 3.10.0.
This version includes features that were requested by users for a while, alongside quality-of-life improvements and bug fixes.
Here is a brief highlight of the most important additions:

**Ortholog Mapping:** *RNAlysis* can now map genes to their closest orthologs in different organisms.
You can map genes to their orthologs using four different databases - Ensembl, Panther, PhylomeDB, and OrthoInspector - extracting both one-to-one and one-to-many ortholog relationships and filtering them based on their reliability.

**Discovering Paralogs:** In the same vein, *RNAlysis* now facilitates the discovery of paralogs within a specific organism, using either the Ensembl or Panther databases.

**New visualization and analysis options for Principal Component Analysis (PCA):** We've introduced new functions and parameters to allow you to get more out of your principal component analysis.

I would also like to extend my personal apology for the delay in bringing you this update.
Due to personal reasons, this release, originally scheduled for the end of August, took longer than expected.
Your patience and support have been invaluable, and I'm eager to share these exciting additions with you.
Thank you for being a part of the RNAlysis community, and stay tuned for more updates in the near future!

Added
*******
* Added new functions to the filtering module that map genes to their closest orthologs in a different organism, using four different databases: Ensembl, Panther, PhylomeDB, and OrthoInspector.
* Added new functions to the filtering module that find paralogs of genes in a given organism, using two different databases: Ensembl and Panther.
* Added new function 'Sort table by contribution to a Principal Component (PCA)' (CountFilter.sort_by_principal_component), which allows sorting of genes in a count matrix by their contribution (gene loadings) to a principal component.
* Added a new parameter called 'legend' to 'Principal Component Analysis (PCA) plot' (CountFilter.pca), which allows users to display a legend on the PCA plot with a name for each sample group/color.

Changed
********
* RNAlysis now issues a warning when users run PCA or PCA-based functions on an unnormalized count matrix.
* The 'seek_fusion_genes' and 'learn_bias' arguments for kallisto quantification (fastq.kallisto_quant_single_end and fastq.kallisto_quant_paired_end), which were depracated in kallisto versions >0.48, are no longer displayed on the graphical interface. Old Pipelines that contain these arguments will still run, but new Pipelines will not contain them.
* Long-running functions now run in the background even when the 'inplace' parameter is set to True, instead of freezing the entire graphical interface.

Fixed
*******
* Fixed bug where functions would sometimes fail to run without displaying an error message.
* Fixed bug where progress bars in the graphical interface would sometimes not disappear after reaching 100% completion.
* RNAlysis should no longer display warning messages about graph layout when graphs are scaled down.
* Fixed bug where the clustergram function (CountFilter.clustergram) would raise an error with specific sets of dependecy versions.
* Loading tables no longer raises a depracation warning when using newer versions of Pandas.

3.9.2 (2023-06-23)
------------------
This patch contains bug fixes and improved functionality for enrichment lollipop plots,
as well as bug fixes for issues with the stand-alone version.

Changed
********
* Single-set enrichment result tables now contain observed/expected values based on the XL-mHG test cutoff.
* When loading a differential expression design matrix, RNAlysis now issues an error if the design matrix column names contain invalid characters.
* Updated the scaling of enrichment lollipop plots to make small 'observed' values easier to discern.

Fixed
*******
* Fixed bug where an error message would sometimes appear after RNAlysis finishes generating an automatic session report on the stand-alone app.
* Fixed bug where enrichment lollipop plots in horizontal mode would display the observed/expected values in reverse order.
* Fixed bug where enrichment lollipop plots and the 'show_exp' parameter would not work on single-set enrichment data.
* Fixed bug where sometimes the tab label of clustering/differential expression output tables would not match the name of the generated table.

3.9.1 (2023-06-19)
------------------

Version 3.9.1 of RNAlysis introduces several improvements and fixes to further improve your analysis experience.
The release includes new optional parameters for single-set enrichment functions, compatibility improvements with newer Python versions,
improved error messaging for R scripts, and adresses minor issues related to enrichment analysis, documentation, plotting parameters, and Pipeline saving.

Added
*******
* Added new optional parameters to single-set enrichment functions, allowing users to determine the top and bottom cutoffs for the XL-mHG test ("X" and "L").

Changed
********
* RNAlysis single-set enrichment analysis using the XL-mHG test now supports Python versions >= 3.8.
* RNAlysis stand-alone app is now built on Python 3.11, improving overall performance.
* Error messages caused by running R tools such as DESeq2, Limma-Voom, and FeatureCounts will now clearly state the reason the script failed, making it easier to understand what went wrong.

Fixed
*******
* Fixed bug where enrichment analysis would raise an error when running enrichment analysis on a gene set with no relevant annotations, or a gene set that does not intersect at all with the background gene set.
* Added missing documentation for plotting parameters in some enrichment functions.
* Depracation Warning should no longer appear when generating a box-plot or enhanced box-plot with scatter=True (CountFilter.box_plot, CountFilter.enhanced_box_plot)
* Fixed bug in featureCounts single-end mode where the 'output_folder' parameter could appear as disabled.
* Fixed bug where RNAlysis would raise an error message after saving a Pipeline, even when the Pipeline was saved successfully.

3.9.0 (2023-06-09)
------------------
Version 3.9.0 of *RNAlysis* introduces several enhancements and fixes to improve your experience.
The release includes additional enrichment plot styles, a new option for PCA plots,
the ability to load and save data tables in Parquet format, and new actions in the Help menu for reporting issues and suggesting improvements.
The update also improves the performance of various functions, ensures consistency in font and theme settings,
and addresses multiple bug fixes, including issues with automatic session reports and visualization functions.

Added
*******
* Added additional parameters to enrichment bar plots (enrichment.enrichment_bar_plot), including a new plot style ('lollipop') and observed/expected labels on the graph.
* Added a new parameter to Principal Component Analysis plots (CountFilter.pca) 'plot_grid', which can enable or disable adding a grid to PCA plots.
* RNAlysis can now load and save data tables in Parquet format (.parquet)
* Added new actions to the Help menu, allowing users to report issues, suggest issues, or open discussions.

Changed
********
* Functions in the FASTQ model are now added to automatic session reports.
* Many of the functions in RNAlysis should now run faster.
* Font type, size, and color for help tooltips should now match the global font settings.
* True/False toggle switches now scale with font size.
* When loading data tables into RNAlysis, you will now see only supported file formats by default.
* Clustering PCA plots are now plotted in proportion to the % variance explained by each PC.
* The legend in clustering PCA plots is now draggable.

Fixed
*******
* Fixed bug where data tables generated through the FASTQ model would not display properly in automatic session reports.
* Fixed bug where graphs generated through the Visualize Gene Sets window would not be added to automatic analysis report.
* When saving a file through the graphical interface, automatically-suggested filenames no longer contain illegal characters.
* Improved clarity of error message when R installation folder is not found.
* Fixed bug where some input parameter widgets in the RNAlysis graphical interface would not display properly.
* RNAlysis now provides a clearer warning message when attempting to run HDBSCAN clustering, if the hdbscan package is not installed.
* Label text in PCA plots and hierarchical clustergrams should no longer be cropped outside of the visible region of the plot.
* Fixed bug where some visualization functions, such as pair-plot (CountFilter.pairplot) would not display properly due to version mismatches between pandas and seaborn.
* Improved clarity of error messages when external apps' (kallisto, bowtie2, etc) installation folders are not found.
* Fixed bug where running the RNAlysis graphical interface on a new computer would sometimes raise an error (thanks to `NeuroRookie <https://github.com/NeuroRookie>`_ in `#25 <https://github.com/GuyTeichman/RNAlysis/issues/25>`_).
* Fixed a bug where the 'min_samples' parameter in HDBSCAN clustering could not be disabled.
* Fixed a bug where applying a function to a gene set with inplace=False would cause the new gene set to be called 'New Table'.
* Fixed a bug where RNAlysis would display the message "Pipeline saved successfully", even when the user cancels the save operation.

New Contributors
*****************
* `NeuroRookie`_ in `#25`_

3.8.0 (2023-05-07)
------------------
Version 3.8.0 of *RNAlysis* comes with several exciting new features, including the ability to generate interactive analysis reports automatically.
This feature allows users to create an interactive graph of all the datasets they loaded into RNAlysis, the functions applied, and the outputs generated during the analysis.
You can read more about this feature in the `Tutorial chapter <https://guyteichman.github.io/RNAlysis/build/tutorial.html#create-analysis-report>`_ and the `User Guide chapter <https://guyteichman.github.io/RNAlysis/build/user_guide_gui.html#rnalysis-interactive-analysis-reports>`_.

The new release also includes Pipelines for FASTQ functions, the ability to export normalization scaling factors, and other changes to improve the software's performance.
RNAlysis now supports Python 3.11, and many functions should now run faster. The software's graphic interface has also improved significantly, and users will now see a clearer error message when attempting to load unsupported file formats.
Lastly, the release also fixes several bugs.

Please note that, since Python 3.7 will be reaching end of life as of June 2023, new versions of *RNAlysis* will no longer support it.

Added
*******
* You can now generate interactive analysis reports automatically using the RNAlysis graphical interface. Read more about this feature `here <https://guyteichman.github.io/RNAlysis/build/user_guide_gui.html>`_.
* Added Pipelines for the FASTQ module (SingleEndPipeline and PairedEndPipeline), allowing you to apply a series of functions (adapter trimming, alignment, quantification, etc) to a batch of sequence/alignment files.
* Added new parameter 'return_scaling_factors' to normalization functions, that allows you to access the scaling factors calculated by RNAlysis.
* Added new parameter 'gzip_output' to CutAdapt adapter trimming (fastq.trim_adapters_single_end and fastq.trim_adapters_paired_end), allowing users to determine whether or not the output files will be gzipped.

Changed
*******
* RNAlysis now supports Python 3.11, and no longer supports Python 3.7.
* Many of the functions in RNAlysis should now run faster.
* The RNAlysis graphic interface should now boot up significantly faster.
* RNAlysis now shows an easier to understand error message when users attempt to load a table in an unsupported format (e.g. Excel files).
* CutAdapt adapter trimming (fastq.trim_adapters_single_end and fastq.trim_adapters_paired_end) now outputs non-gzipped files by default.
* Standardized all plotting functions in the filtering module to return a matplotlib Figure object, which can be further modified by users.

Fixed
*******
* RNAlysis failing to map gene IDs during GO enrichment analysis should no longer raise an error (thanks to `clockgene <https://github.com/clockgene>`_ in `#16 <https://github.com/GuyTeichman/RNAlysis/issues/16>`_).
* Fixed bug where the Command History window would not display history of the current tab immediately after clearing the current session.
* Fixed bug where adapter trimming would fail to run when using CutAdapt version >= 4.4.0.
* Fixed bug where 'Filter rows with duplicate names/IDs' (Filter.filter_duplicate_ids) would raise an error when applied to some tables.

New Contributors
*****************
* `clockgene`_ in `#16`_


3.7.0 (2023-04-07)
------------------
This version introduces small RNA read alignment using ShortStack, new filtering functions, a new optional parameter for Principal Component Analysis, improvements to the graphical interface, and bug fixes.

Added
*******
* Added small RNA read alignment using ShortStack (fastq.shortstack_align_smallrna).
* Added new filtering function 'Filter specific rows by name' (Filter.filter_by_row_name).
* Added new filtering function 'Filter rows with duplicate names/IDs' (Filter.filter_duplicate_ids).
* Function parameters in pop-up windows in the graphical interface can now be imported and exported.
* Added new parameter to Principal Component Analysis (CountFilter.pca) 'proportional_axes', that allows you to make the PCA projection axes proportional to the percentage of variance each PC explains.
* Improved clarity of error messages in the graphical interface.

Changed
*******
* Tables loaded into RNAlysis that use integer-indices will now be converted to use string-indices.
* Refactored CountFilter.from_folder into CountFilter.from_folder_htseqcount, and added a new CountFilter.from_folder method that accepts a folder of count files in any format.
* In the RNAlysis graphical interface, optional parameters that can be disabled will now display the hint "disable this parammeter?" instead of "disable input?".
* Added optional parameter 'ylim' to 'create enrichment bar-plot' function (enrichment.enrichment_bar_plot), allowing users to determine the Y-axis limits of the bar plot.
* Updated function signatures of 'Visualize gene ontology' and 'Visualize KEGG pathway' (enrichment.gene_ontology_graph and enrichment.kegg_pathway_graph) to make more sense.
* Removed parameter 'report_read_assignment_path' from featureCounts feature counting (fastq.featurecounts_single_end and fastq.featurecounts_paired_end).
* The RNAlysis graphical interface should now load more quickly.
* Progress bars in the graphical interface should now reflect elapsed/remaining time for tasks more accurately.

Fixed
*******
* Fixed bug in the function 'Split into Higly and Lowly expressed genes' (Filter.split_by_reads) where the two resulting tables would be named incorrectly (highly-expressed genes would be labeled 'belowXreads' and vice-versa).
* Fixed bug where the 'column' parameter of some functions ('Filter by percentile', 'Split by percentile', 'Filter with a number filter', 'Filter with a text filter') would not automatically detect column names in the graphical interface.
* Fixed bug where the 'numerator' and 'denominator' parameters of of the function 'Calculate fold change' would not automatically detect column names in the graphical interface.
* Fixed bug where tables with integer-indices could not be visualized properly through the graphical interface.
* Fixed bug in the function 'featureCounts single-end' (fastq.featurecounts_single_end) where setting parameter 'report_read_assignment' to any value other than None would raise an error.
* Functions that take column name/names as input (transform, filter_missing_values, filter_percentile, split_percentile) can now be applied to fold change tables (FoldCangeFilter objects).
* Fixed bug where the description for the 'n_bars' parameter of the 'create enrichment bar-plot' function (enrichment.enrichment_bar_plot) would not display properly.
* 'Visualize gene ontology' and 'Visualize KEGG pathway' (enrichment.gene_ontology_graph and enrichment.kegg_pathway_graph) now have proper parameter descriptions.
* Fixed bug where in-place intersection and difference in the filtering module would fail when using recent versions of pandas.
* Fixed bug where graphs generated through the Visualize Gene Sets window would not immediately display when using the RNAlysis stand-alone app.

3.6.2 (2023-03-25)
------------------
This version introduces quality-of-life changes to the graphical interface, as well as bug fixes.

Added
*******
* Sample groupings in functions such as PCA, Pair-plot, etc., can now be imported and exported through the graphical interface.

Fixed
*******
* Fixed bug where the stand-alone Mac version of RNAlysis would sometimes fail to map gene IDs (directly or in enrichment analysis).

3.6.1 (2023-03-22)
------------------
This version introduces minor bug fixes.

Changed
********
* DESeq2 automatic installation should now work more reliably.

Fixed
******
* Fixed bug where PCA plots would not display on the stand-alone app unless another visualization function was applied afterwards.
* Fixed bug where Pipelines that contain specific functions (such as translating gene IDs/filtering biotypes from GTF file) would fail to run through the graphical interface.
* GO Annotations annotated by ComplexPortal are now supported by RNAlysis.

3.6.0 (2023-03-07)
------------------
This version introduces improvements to the usability and clarity of the graphic interface,
new methods for automatic estimation of the number of clusters in a dataset,
and various bug fixes.

Added
******
* Added three new methods for automatic estimation of the number of clusters in a dataset: Callinski-Harabasz, Davies-Bouldin, and Bayesian Information Criterieon.
* Added a 'Close all Figures' actions to the 'View' menu of the *RNAlysis* graphic interface.
* Added an 'interactive' parameter to Volcano Plots (DESeqFilter.volcano_plot) and 'Scatter Sample Vs Sample' (CountFilter.scatter_sample_vs_sample), allowing user to label data points interactively by clicking on them.
* Added more optional plotting parameters to Volcano Plots (DESeqFilter.volcano_plot) and 'Scatter Sample Vs Sample' (CountFilter.scatter_sample_vs_sample).

Changed
********
* Progress bars are now integrated into the main *RNAlysis* window instead of opening as a dialog box.
* Information about running proccesses and functions is now displayed in the main *RNAlysis* window.
* It is now possible to cancel queued jobs through the *RNAlysis* graphic interface.
* When loading multiple data tables at the same time, it is now possible to change the table type of all data tables at once, instead of one-by-one.

Fixed
******
* RNAlysis KEGG enrichment should now match the new KEGG annotation format from March 1st 2023.
* Fixed bug where importing *RNAlysis* would raise ImportError when cutadapt is not installed.
* Fixed bug where the 'Run' button in the Enrichment Analysis window would grey out whenever the enrichment dataset is changed.
* Fixed bug where the *RNAlysis* stand-alone versions were unable to export Figures in specific formats (e.g. PDF, SVG).
* Fixed bug where functions that depend on R scripts (such as DESeq2 and limma) would sometimes fail to run on MacOS (thanks to Matthias Wilm and `sandyl27 <https://github.com/sandyl27>`_ in `#12 <https://github.com/GuyTeichman/RNAlysis/issues/12>`_).
* Fixed bug where running limma-voom with a design matrix whose column names contained spaces or special characterse would raise an error.
* Fixed bug where the 'highlight' parameter of CountFilter.scatter_sample_vs_sample would not work when used through the graphic interface.
* Fixed bug where enrichment analysis would sometimes fail to run when 'exclude_unannotated_genes' is set to False.
* Fixed bug where translate_gene_ids() would fail for RankedSet objects.
* Fixed bug where filtering gene sets by user-defined attributes (FeatureSet.filter_by_attribute()) would occasionally fail to run.

New Contributors
*****************
* `sandyl27`_ in `#12`_

3.5.2 (2023-02-23)
------------------
This version includes bug fixes for a few quality-of-life issues which were introduced in version RNAlysis 3.5.0.

Changed
********
* It is now possible to change the parallel backend of performance-intensive functions such as clustering an enrichment analysis in non-standalone versions of RNAlysis.
* Expanded the Frequently Asked Questions page.
* Added Perl as a dependency for Windows users on the How to Install page.
* Automatic row colours in column-picking tables should no longer mismatch with font colors in a way that obscures visibility.

Fixed
*****
* Fixed bug where occasionally newly-created tabs would open twice instead of once.
* Fixed bug where the 'Add Function' button of the Pipeline window would appear in the wrong location.
* Fixed bug where RNAlysis sessions saved after version 3.5.0 which contain gene sets would raise an error when loaded.
* Fixed typos in the RNAlysis tutorial.

3.5.1 (2023-02-22)
------------------
This version introduces minor bug fixes.

Fixed
*****
* Fixed bug where the *RNAlysis* stand-alone app would sometimes fail to run CutAdapt (thanks to Matthias Wilm).
* Fixed bug where the User Guide action in the graphical interface would point to the Tutorial, and vice versa.
* The X and Y axis labels on volcano plots should now format the 'log' in the label correctly.

3.5.0 (2023-02-08)
------------------
This version introduces new features such as differential expression analysis with the Limma-Voom pipeline,
customizable databases for the quick-search function, basic filtering and procrssing functions for gene sets,
improved progammatic API for FeatureSet and RankedSet objects, and RPKM and TPM read count normalization.
Several changes have been made to improve the user experience, including updated documentation,
improved clarity of function tooltips, clearer output formats, and faster download speeds for tutorial videos.

Added
*******
* Added differential expression analysis with the Limma-Voom pipelines (CountFilter.differential_expression_limma_voom)
* You can now select which databases to display in the right-click quick-search menu, using the settings menu.
* Gene sets now support some basic operations (filtering, gene ID translating, etc) through the graphical interface.
* enrichment.FeatureSet and enrichment.RankedSet now support some filtering operations from the filtering module (such as filtering by user-defined attributes, GO terms, or KEGG pathways).
* Added reads-per-kilobase-million (RPKM) and transcripts-per-million (TPM) normalization methods (CountFilter.normalize_to_rpkm() and CountFilter.normalize_to_tpm()).

Changed
********
* Classes enrichment.FeatureSet and enrichment.RankedSet now inherit from Python base-class set, and can be interacted with like other Python sets. The old API and attributes of these classes were maintained as they were.
* Improved documentation for some functions.
* Function selection tooltips should now display information more clearly.
* Pipelines that contain consecutive clustering/splitting functions will now return their outputs in a clearer format.
* Enrichment bar-plots should now adjust the x-axis limits more tightly to fit the displayed data.
* Improved clarity of automatic titles in enrichment plots.
* Download/update speed of tutorial videos has improved significantly.

Fixed
******
* Fixed bug where Pipelines would not always properly run 'translate_gene_ids'

3.4.2 (2023-02-01)
------------------
This version introduces minor bug fixes.

Fixed
******
* Fixed bug where updating RNAlysis through the graphical interface would not update some of the optional dependencies.
* Fixed typos in the documentation.

3.4.0 (2023-02-01)
------------------
From this release forward, *RNAlysis* is made available as a stand-alone app for Windows and MacOS. You can download these stand-alone versions from the GitHub Releases page.
In addition, new features were added, including new plots, filtering functions, integration of the external tools bowtie2 and featureCounts, and the ability to generate Gene Ontology Graphs and KEGG Pathway Graphs without running enrichment analysis from scratch.

Added
******
* Added a Scree Plot (explained variance per PC) to Principle Component Analysis
* Added CountFilter.split_by_principal_component(), allowing users to filter genes based on their contributions (loadings) to PCA Principal Components.
* Added Filter.drop_columns
* Added support for the Sharpened Cosine distance metric in clustering analyses
* KEGG enrichment can now generate KEGG pathway graphs for pathways that were found to be statistically significant
* Added functions to the enrichment module that can generate KEGG Pathway or Gene Ontology plots based on previously-generated enrichment results
* You can now clear the *RNAlysis* cache from the *RNAlysis* GUI, or through the general module.
* Added bowtie2 alignment to the fastq module.
* Added FeatureCounts feature-counting to the fastq module.
* You can now choose whether or not to discard genes from enrichment analysis if they have no annotations associated with them.
* When right-clicking on a specific cell in a table or line in a gene set view, a context menu will open, allowing you to copy the associated value, or look it up in one of few common biology databases.
* Added sections to the programmatic user guide about the `fastq` module.

Changed
********
* Replaced the 'parallel' parameter in enrichment functions with the 'parallel_backend' parameter, allowing users to choose which parallel backend (if any) will be used in the function.
* Added 'parallel_backend' parameter to all clustering functions under the filtering module.
* When generating Gene Ontology/KEGG Pathway graphs, users can choose whether or not to generate the figure in an additional separate file.
* Updated type annotations of some functions to be more precise and helpful (for example, setting a lower bound on some int function parameters).
* The colorbar ticks in enrichment bar plots now match the axis ticks on the main axis.
* Slight improvements in GUI performance, stability, and looks.
* Slight improvements in performance of enrichment analysis when examining a small number of attributes.
* enrichment.plot_enrichment() was replaced by enrichment.enrichment_bar_plot().
* CountFilter.differential_expression() has new optional parameter `output_folder`, which allows users to save the generated data tables and the R script that generated them into a specified folder.

Fixed
******
* In CountFilter.differential_expression_deseq2(), fixed a bug where design matrix files with non-comma delimiters would cause an error (thanks to `Mintxoklet <https://github.com/Mintxoklet>`_ in `#7 <https://github.com/GuyTeichman/RNAlysis/issues/7>`_)
* Fixed bug where setup.py would install a directory named tests into site-packages folder (thanks to `Bipin Kumar <https://github.com/kbipinkumar>`_ in `#9 <https://github.com/GuyTeichman/RNAlysis/issues/9>`_)
* Fixed bug where the windows of some functions (differential expression, adapter trimming, etc) did not show a link to the function's documentation page.
* Fixed typos in some parts of the *RNAlysis* documentation
* When filtering a table by a single user-defined attribute, the automatic table name will now be more informative about the operation applied.
* Fixed bug where occasionally a Pipeline or Function would generate multiple tables of the same name, but only one of them will appear in the GUI.
* Fixed bug where occasionally significance asterisks on enrichment bar-plots would appear in the wrong location.
* Fixed bug where fastq.create_kallisto_index() (Create Kallisto Index) would not make use of the `make_unique` parameter (thanks to Matthias Wilm)

Removed
********
* Removed the previously-deprecated functions `enrichment.enrich_randomization()` and `enrichment.enrich_hypergeometric()`.



New Contributors
*****************
* `Mintxoklet`_ in `#7`_
* `Bipin Kumar`_ in `#9`_
* Matthias Wilm

3.3.0 (2022-12-02)
------------------
* This version introduced quality-of-life improvements to the graphical user interface.

Added
******
* Added a Frequently Asked Questions page, and linked all *RNAlysis* help material inside the graphical interface Help menu.
* Pipelines can now be edited and deleted through the Pipeline menu of the graphical interface.
* Added Contributing page to the documentation, with basic guidelines on how you can contribute to the *RNAlysis* project!

Changed
*******
* All open tabs are now always visible in the main menu screen. Tab names are now shortened with ellipsis if nessecary.
* The right-click context menu of the main menu tabs now allows users to open a new tab at a specific position, or close a specific tab/tabs to the right/tabs to the left/all other tabs.
* *RNAlysis* documentation is now split into GUI documentation (quick-start video guide, tutorial, GUI user guide), and programmatic documentation (programmatic user guide)
* Improved readability of *RNAlysis* logs
* Pipelines are now exported with additional metadata - the version of *RNAlysis* they were exported from, and the date and time it was exported. This metadata should not affect Pipelines that were created in older versions, and does not affect the way Pipelines are applied to data tables.

Fixed
******
* *RNAlysis* now warns users if they attempt to overwrite an existing Pipeline.
* Fixed an incorrect keyboard shortcut for Export Pipeline action

3.2.2 (2022-11-25)
------------------


Fixed
******
* Fixed bug with DESeq2 automatic installation on Windows computers.

3.2.1 (2022-11-25)
------------------

Changed
*******
* Updated citation information for *RNAlysis*

Fixed
******
* Fixed typos in the *RNAlysis* tutorial

3.2.0 (2022-11-23)
------------------
* This version introduces quality-of-life changes to the graphical user interface, functions for translating gene IDs and running differential expression analysis, and extends RNAlysis to support Python versions 3.9 and 3.10.

Added
******
* Added Filter.translate_gene_ids()
* Added CountFilter.differential_expression_deseq2()
* Added Filter.filter_by_kegg_annotations()
* Added Filter.filter_by_go_annotations()
* Added CountFilter.average_replicate_samples()
* Added fastq module that contains adapter-trimming functions utilizing CutAdapt, and mRNA-sequencing quantification using kallisto.

Changed
*******
* Added additional plotting parameters to visualization functions.
* Improved performance of some aspects of the graphical user interface.
* RNAlysis' basic features are now supported on Python versions 3.9 and 3.10.
* CountFilter.pca() now generates a plot for *every* pair of Principal Components requested by the user.
* CountFilter.split_clicom() now supports clustering each batch of replicates separately, using the 'replicates_grouping' parameter
* Biotype-based filtering and summary can now be done based on GTF annotation files instead of a Biotype Reference Table.
* Filter.biotypes() was refactored into Filter.biotypes_from_ref_table()
* Filter.filter_biotype() was refactored into Filter.filter_biotype_from_ref_table()

Fixed
******
* Users can now queue multiple computationally-intense enrichment/clustering tasks while another task is running.
* Fixed a bug where sometimes some function parameters would disappear from the graphical user interface.
* Fixed a bug where exceptions during computationally-intense tasks would cause *RNAlysis* to crash.
* Auxillary windows are now properly minimized when analysis starts, and restored when analysis ends or encounters an error.

3.1.0 (2022-10-16)
------------------
* This version introduces new count matrix normalization methods, as well as MA plots and minor bug fixes.

Added
******
* Added the visualization function ma_plot() for CountFilter
* Added functions for the normalization functions Relative Log Ratio (RLE), Trimmed Mean of M-values (TMM), Median of Ratios (MRN), Quantile normalization (quantile)

Changed
*******
* CountFilter.normalize_to_rpm() was renamed to CountFilter.normalize_to_rpm_htseqcount(), and was supplemented by the more general function for normalizing to Reads Per Million CountFilter.normalize_to_rpm()

Fixed
******
* Fixed a bug where some elements of the graphical user interface would not display correctly

3.0.1 (2022-10-12)
------------------
* This version fixes a bug with displaying the tutorial videos in the graphical user interface.


3.0.0 (2022-10-10)
------------------
* This version introduces a graphical user interface for RNAlysis, as well as new functions for KEGG Pathways enrichment analysis.


Added
******
* RNAlysis now includes a graphical user interface
* Pipelines can now be imported and exported
* Enrichment and single-set-enrichment for KEGG pathway data

Changed
*******
* Added function FeatureSet.user_defined_enrichment(), which will replace FeatureSet.enrich_hypergeometric() and FeatureSet.enrich_randomization()
* Updated signature of enrichment.venn_diagram
* enrichment.venn_diagram and enrichment.upset_plot can now be generated on a user-supplied FIgure
* Clustering functions now apply a power transform to count data prior to clustering by default
* Non-deprecated enrichment functions no longer filter the background set by biotype by default
* Changed signature of CountFilter.pca, CountFilter.box_plot, CountFilter.enhanced_box_plot, CountFilter.clustergram, and CountFilter.pairplot to ensure consistency among visualization functions.

Fixed
******
* enrichment.venn_diagram can now be plotted with outlines when the circles are unweighted
* Fixed bug in Pipeline.apply_to() where a Filter object would be returned even when the Pipeline was applied inplace


2.1.1 (2022-07-05)
------------------
* This version fixes issues with running GO enrichment that resulted from recent changes to UniProt's API.  Moreover, this version slightly improves the performance of some functions.

Changed
*******
* Fixed issues with running GO enrichment that resulted from changes to UniProt's API.
* Some functions that fetch annotations now cache their results, leading to improved runtimes.
* Updated the documentation of some functions to better reflect their usage and input parameters.

2.1.0 (2022-04-16)
------------------
* This version introduces multiple new features, as well as generally improved graphs and quality-of-life changes.

Added
******
* GO enrichment can now generate Ontology Graphs for the statistically significant GO terms.
* Added CountFilter.split_clicom(), an implementation of the CLICOM ensemble-based clustering method (Mimaroglu and Yagci 2012).
* Added Filter.transform(), a method that can transform your data tables with either predefined or user-defined transformations.

Changed
*******
* CountFilter.pairplot() now uses a logarithmic scale by default.
* Visually improved the graphs generated by many functions, including CountFilter.pairplot() and CountFilter.plot_expression().
* The clusters resulting from all clustering functions are now sorted by size instead of being sorted randomly.

Fixed
******
* Minor bug fixes.


2.0.1 (2022-04-02)
------------------
* This version introduces small bug fixes, as well as a new function in the Filtering module.

Added
******
* Added Filter.majority_vote_intersection(), which returns a set/string of the features that appear in at least (majority_threhold * 100)% of the given Filter objects/sets.

Changed
*******
* When mapping/inferring taxon IDs during GO enrichment analysis, organisms will now be prioritized based on their taxon ID values (numerically lower IDs will be considered to be more relevant).

Fixed
******
* Fixed bug that occured when mapping/inferring taxon IDs during GO enrichment analysis, where integer taxon IDs would be matched by name similarity before trying an exact ID match, leading to spurious matches.
* Fixed bug that occursed when plotting clustering results with style='all' on Python 3.8.

2.0.0 (2021-12-05)
------------------
* This version introduces new method to cluster your read count matrices, including K-Means/Medoids clustering, Hierarchical clustering, and HDBSCAN.
* This version introduces many new ways to perform enrichment analysis and to visualize your results, including highly customizable GO Enrichment, enrichment based on ranked lists of genes, and enrichment for non-categorical attributes.
* This version introduces Pipelines - a quicker and more convenient way to apply a particular analysis pipeline to multiple Filter objects.
* This version improves the performance of many functions in RNAlysis, and in particular the performance of randomization tests.
* This version includes changes to names and signatures of some functions in the module, as elaborated below.


Added
******
* Added class Pipeline to filtering module, which applies a series of filter functions to specified Filter objects.
* Added CountFilter.split_kmeans(), CountFilter.split_kmedoids(), CountFilter.split_hierarchical() and CountFilter.split_hdbscan(), which split your read count matrices into clusters with similar expression patterns.
* Added class RankedSet to enrichment module, which accepts a ranked list of genes/features, and can perform single-list enrichment analysis
* Added RankedSet.single_set_enrichment(), which can perfofm single-list enrichment analysis of user-defined attributes using XL-mHG test (see `Eden et al. (PLoS Comput Biol, 2007) <https://dx.doi.org/10.1371/journal.pcbi.0030039>`_  and `Wagner (PLoS One, 2015) <https://dx.doi.org/10.1371/journal.pone.0143196>`_ ).
* Added FeatureSet.go_enrichment() and RankedSet.single_set_go_enrichment(), which let you compute Gene Ontology enrichment for any organism of your choice, and filter the GO annotations used according to your preferences.
* Added FeatureSet.enrich_hypergeometric(), which can perform enrichment analysis using the Hypergeometric Test.
* Added more visualization functions, such CountFilter.enhanced_box_plot().
* Added FeatureSet.change_set_name(), to give a new 'set_name' to a FeatureSet object.


Changed
*******
* FeatureSet.enrich_randomization_parallel() was deprecated. Instead, you can compute your enrichment analysis with parallel computing by calling FeatureSet.enrich_randomization() with the argument 'parallel_processing=True'. Moreover, parallel session will now start automatically if one was not already active.
* Improved running time of enrich_randomization() about six-fold.
* Filter objects can be created from any delimiter-separated file format (.csv, .tsv, .txt, etc).
* CountFilter.pca() can now be plotted without labeled points.
* Filter.index_string is now sorted by the current order of indices in the Filter object, instead of by alphabetical order.
* CountFilter.violin_plot() now accepts a y_title argument.
* Added more optional arguments to visualization functions such as CountFilter.violin_plot() and CountFilter.clustergram().
* Automatic filenames for Filter objects should now reflect more clearly the operations that were performed.
* The DataFrame returned by enrich_randomization() and enrich_randomization_parallel() now contains the additional column 'data_scale', determined by the new optional argument 'data_scale'.
* The columns 'n obs' and 'n exp' in the DataFrame returned by enrich_randomization() and enrich_randomization_parallel() were renamed to 'obs' and 'exp' respectively.
* FeatureSets no longer support in-place set operations (intersection, union, difference, symmetric difference). Instead, these functions return a new FeatureSet.
* Filter.biotypes_from_ref_table() now accepts the boolean parameter 'long_format' instead of the str parameter 'format'.
* Filter.biotypes_from_ref_table() and FeatureSet.biotypes_from_ref_table() now count features which do not appear in the Biotype Reference Table as '_missing_from_biotype_reference' instead of 'not_in_biotype_reference'.

Fixed
******
* Updated type-hinting of specific functions.
* Filter.biotypes_from_ref_table() and FeatureSet.biotypes_from_ref_table() now support Biotype Reference Tables with different column names.
* Generally improved performance of RNAlysis.
* Fixed bug in Filter.filter_percentile() where the value at the exact percentile speficied (e.g. the median for percentile=0.5) would be removed from the Filter object.
* Fixed bug in enrichment.FeatureSet, where creating a FeatureSet from input string would result in an empty set.
* Various minor bug fixes.





1.3.5 (2020-05-27)
------------------
* This version introduces minor bug fixes and a few more visualization options.

Added
******
* Added Filter.filter_missing_values(), which can remove rows with NaN values in some (or all) columns.
* Added the visualization function CountFilter.box_plot().

Changed
*******
* Updated docstrings and printouts of several functions.
* Slightly improved speed and performance across the board.
* Filter.feature_string() is now sorted alphabetically.
* Enrichment randomization functions in the enrichment module now accept a 'random_seed' argument, to be able to generate consistent results over multiple sessions.
* Enrichment randomization functions can now return the matplotlib Figure object, in addition to the results table.


Fixed
******
* Fixed DepracationWarning on parsing functions from the general module.
* Fixed bug where saving csv files on Linux systems would save the files under the wrong directory.
* Fixed a bug where UTF-8-encoded Reference Tables won't be loaded correctly
* Fixed a bug where enrichment.upsetplot() and enrichment.venn_diagram() would sometimes modify the user dict input 'objs'.
* Fixed a bug in CountFilter.pairplot where log2 would be calculated without a pseudocount added, leading to division by 0.




1.3.4 (2020-04-07)
------------------
* This version fixed a bug that prevented installation of the package.


Changed
*******
* Updated docstrings and printouts of several functions


Fixed
******
* Fixed a bug with installation of the previous version






1.3.3 (2020-03-28)
------------------
* First stable release on PyPI.


Added
******
* Added Filter.sort(), and upgraded the functionality of Filter.filter_top_n().
* Added UpSet plots and Venn diagrams to enrichment module.
* User-defined biotype reference tables can now be used.
* Filter operations now print out the result of the operation.
* Enrichment randomization tests now also support non-WBGene indexing.
* Filter.biotypes_from_ref_table() and FeatureSet.biotypes_from_ref_table() now report genes that don't appear in the biotype reference table.
* Filter.biotypes_from_ref_table() can now give a long-form report with descriptive statistics of all columns, grouped by biotype.
* Added code examples to the user guide and to the docstrings of most functions.


Changed
*******
* Changed argument order and default values in filtering.CountFilter.from_folder_htseqcount().
* Changed default title in scatter_sample_vs_sample().
* Changed default filename in CountFilter.fold_change().
* Settings are now saved in a .yaml format. Reading and writing of settings have been modified.
* Changed argument name 'deseq_highlight' to 'highlight' in scatter_sample_vs_sample(). It can now accept any Filter object.
* Updated documentation and default 'mode' value for FeatureSet.go_enrichment().
* Updated the signature and function of general.load_csv() to be clearer and more predictable.
* Changed argument names in CountFilter.from_folder_htseqcount().
* Modified names and signatures of .csv test files functions to make them more comprehensible.
* Renamed 'Filter.filter_by_ref_table_attr()' to 'Filter.filter_by_attribute()'.
* Renamed 'Filter.split_by_ref_table_attr()' to 'Filter.split_by_attribute()'.
* Renamed 'Filter.norm_reads_with_size_factor()' to 'Filter.normalize_with_scaling_factors()'. It can now use any set of scaling factors to normalize libraries.
* Renamed 'Filter.norm_reads_to_rpm()' to 'Filter.normalize_to_rpm()'.
* Made some functions in the general module hidden.


Fixed
******
* Various bug fixes


Removed
********
* Removed the 'feature_name_to_wbgene' module from RNAlysis.






1.3.2 (2019-12-11)
------------------

* First beta release on PyPI.
