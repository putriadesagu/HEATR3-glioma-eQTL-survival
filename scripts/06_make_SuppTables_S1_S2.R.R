# Run in RStudio. Requires data.table; readxl only for the Excel fallback.
suppressPackageStartupMessages(library(data.table))

BASE <- "D:/RESIDENSI/Thesis-PPDS/R"
S1_RUN <- ""     # Optional: exact results/HEATR3_* directory.
S2_RUN <- ""     # Optional: exact results/HEATR3_* directory.
GTEX_PAIRS <- "" # Optional: specify exact paths if duplicate files differ.
GTEX_EGENES <- ""
SENTINEL <- "rs10852606" # same sentinel specified in HEATR3_final_project.R
SURVIVAL_GENES <- c("HEATR3", "LMF1", "EGFR", "RTEL1", "SLC16A8", "MIR4300HG")

build_supp_tables <- function() {
  if (!dir.exists(BASE)) stop("Project folder not found: ", BASE)
  OUT <- file.path(BASE, "results_final")
  dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
  log <- character()
  note <- function(...) {
    msg <- paste0(...)
    log <<- c(log, msg)
    message(msg)
  }
  on.exit(writeLines(c(paste("Run:", Sys.time()), "", log),
                    file.path(OUT, "SuppTables_found_files_and_checks.txt")), add = TRUE)
  save <- function(x, stem, caption, footer = "") {
    fwrite(x, file.path(OUT, paste0(stem, ".csv")), na = "")
    escape <- function(x) gsub("|", "&#124;", gsub("[\r\n]+", " ", x), fixed = TRUE)
    line <- function(x) paste0("| ", paste(escape(x), collapse = " | "), " |")
    cells <- as.data.frame(x)
    cells[] <- lapply(cells, function(x) ifelse(is.na(x), "", as.character(x)))
    rows <- if (nrow(cells)) apply(cells, 1L, line) else character()
    writeLines(c(caption, "", line(names(cells)), line(rep("---", ncol(cells))), rows,
                 "", footer), file.path(OUT, paste0(stem, ".md")))
  }

  find_project <- function(name, dir = BASE) {
    files <- list.files(dir, recursive = TRUE, full.names = TRUE)
    files[basename(files) == name]
  }
  recent_run <- function(name, explicit) {
    if (nzchar(explicit)) return(explicit)
    hits <- find_project(name, file.path(BASE, "results"))
    if (!length(hits)) return("")
    dirname(hits[which.max(file.info(hits)$mtime)])
  }
  resolve_gtex <- function(explicit, name) {
    if (nzchar(explicit)) {
      if (!file.exists(explicit)) stop("GTEx path not found: ", explicit)
      return(explicit)
    }
    hits <- sort(find_project(name))
    if (length(hits) < 2L) return(if (length(hits)) hits[[1L]] else "")
    hashes <- unname(tools::md5sum(hits))
    if (all(!is.na(hashes)) && length(unique(hashes)) == 1L) return(hits[[1L]])
    note("Different GTEx copies found: ", paste(hits, collapse = " | "))
    ""
  }
  rsids <- function(x) {
    x <- tolower(as.character(x))
    x[is.na(x)] <- ""
    unique(unlist(regmatches(x, gregexpr("rs[0-9]+", x)), use.names = FALSE))
  }
  ptext <- function(p) ifelse(is.na(p), NA_character_, formatC(p, format = "e", digits = 2))
  val <- function(d, options, default = NA_character_) {
    hit <- intersect(options, names(d))
    if (length(hit)) d[[hit[[1L]]]] else rep(default, nrow(d))
  }

  tryCatch({
    run <- recent_run("S1_GWAS_input.tsv", S1_RUN)
    s1_source <- if (nzchar(run)) file.path(run, "S1_GWAS_input.tsv") else ""
    if (nzchar(s1_source) && file.exists(s1_source)) {
      gw <- fread(s1_source, sep = "\t", check.names = FALSE)
    } else {
      s1_source <- file.path(BASE, "Glioblastoma_only_GWAS_Raw_Catalog_Export.xlsx")
      if (!file.exists(s1_source)) stop("Missing original GWAS Excel and S1_GWAS_input.tsv.")
      if (!requireNamespace("readxl", quietly = TRUE)) stop("Install readxl for the original GWAS Excel.")
      gw <- as.data.table(readxl::read_excel(s1_source, .name_repair = "minimal"))
    }
    note("S1 source: ", s1_source)
    required <- c("P-VALUE", "REPORTED GENE(S)", "STRONGEST SNP-RISK ALLELE", "SNPS")
    missing <- setdiff(required, names(gw))
    if (length(missing)) stop("GWAS columns missing: ", paste(missing, collapse = ", "))
    p <- suppressWarnings(as.numeric(gw[["P-VALUE"]]))
    keep <- which(is.finite(p) & p < 5e-8)
    sig <- gw[keep[order(p[keep])]]
    if (!nrow(sig)) stop("No associations meet P < 5e-8.")
    sig[, P_numeric := suppressWarnings(as.numeric(get("P-VALUE")))]
    export <- data.table(
      Risk_allele = sig[["STRONGEST SNP-RISK ALLELE"]], SNPs = sig$SNPS,
      Catalog_chromosome = val(sig, "CHR_ID"), Catalog_position = val(sig, "CHR_POS"),
      Reported_genes = sig[["REPORTED GENE(S)"]],
      Mapped_genes = val(sig, c("MAPPED_GENE", "MAPPED GENE(S)")),
      P = sig$P_numeric, OR_or_beta = val(sig, "OR or BETA"),
      Trait = val(sig, "DISEASE/TRAIT"), Author = val(sig, "FIRST AUTHOR"),
      Publication_date = val(sig, "DATE"), PMID = val(sig, "PUBMEDID"))
    # Retain the original attachment's filename for the expected 64-record run.
    stem <- if (nrow(sig) == 64L) "SuppTable_S1b_GWAS_associations_all64" else
      "SuppTable_S1b_GWAS_significant_associations"
    fwrite(export, file.path(OUT, paste0(stem, ".csv")), na = "")

    split_genes <- strsplit(as.character(export$Reported_genes), "[,;]")
    hits <- rbindlist(lapply(seq_along(split_genes), function(i) {
      genes <- unique(trimws(split_genes[[i]]))
      genes <- genes[!is.na(genes) & !toupper(genes) %in%
                       c("", "NR", "NA", "INTERGENIC", "NOT REPORTED", "-")]
      data.table(row = rep.int(i, length(genes)), Gene = genes)
    }))
    if (!nrow(hits)) stop("No usable reported-gene labels.")
    s1 <- hits[, {
      x <- export[row]
      .(Best_association_risk_allele = x$Risk_allele[1L], Best_P = ptext(x$P[1L]),
        OR_or_beta = x$OR_or_beta[1L], Significant_associations = .N,
        Distinct_rsIDs = length(rsids(x$SNPs)),
        Study_PMIDs = paste(unique(na.omit(x$PMID)), collapse = ", "))
    }, by = Gene]
    s1[, GTEx_cortex_pairs := NA_integer_]
    tryCatch({
      pairs_path <- resolve_gtex(GTEX_PAIRS, "Brain_Cortex.v8.signif_variant_gene_pairs.txt.gz")
      genes_path <- resolve_gtex(GTEX_EGENES, "Brain_Cortex.v8.egenes.txt.gz")
      if (!nzchar(pairs_path) || !nzchar(genes_path)) stop("No unambiguous pair of raw GTEx input files.")
      note("GTEx pairs: ", pairs_path, " | eGenes: ", genes_path)
      annot <- unique(fread(genes_path, select = c("gene_id", "gene_name")))
      pairs <- unique(fread(pairs_path, select = c("gene_id", "variant_id")))
      cnt <- merge(pairs, annot, by = "gene_id")[, .N, by = gene_name]
      s1[, GTEx_cortex_pairs := cnt$N[match(Gene, cnt$gene_name)]]
      s1[is.na(GTEx_cortex_pairs), GTEx_cortex_pairs := 0L]
    }, error = function(e) note("GTEx pair counts unavailable: ", conditionMessage(e)))
    s1[, Carried_to_survival_analysis := ifelse(Gene %in% SURVIVAL_GENES, "Yes", "No")]
    setorder(s1, -Carried_to_survival_analysis, Gene)
    check <- sprintf("Catalog: %d records; %d associations at P < 5e-8; %d distinct rsIDs; %d reported-gene labels.",
                     nrow(gw), nrow(sig), length(rsids(sig$SNPS)), nrow(s1))
    note(check)
    save(s1, "SuppTable_S1_GWAS_genes",
         "**Supplementary Table S1.** Reported genes for significant glioblastoma GWAS Catalog associations.",
         paste(check, "Distinct rsIDs do not imply independent loci. Blank GTEx counts indicate unavailable counts."))
  }, error = function(e) note("S1 pending: ", conditionMessage(e)))

  tryCatch({
    run <- recent_run("S2_GTEx_curated_top_eQTL.tsv", S2_RUN)
    audit <- if (nzchar(run)) file.path(run, "S2_GTEx_curated_top_eQTL.tsv") else ""
    identity_path <- file.path(OUT, "P2b_GTEx_variant_identity_verified.csv")
    if (nzchar(audit) && file.exists(audit)) {
      v <- fread(audit)[gene == "HEATR3"]
      s2 <- unique(v[, .(variant_id, sentinel_rsid, slope_ALT = slope,
                         p_eQTL = pval_nominal, S2_identity = identity)])
      note("S2 audit: ", audit)
    } else {
      if (!file.exists(identity_path)) stop("Neither the S2 audit nor P2b identity output was found.")
      s2 <- data.table(variant_id = unique(fread(identity_path)$variant_id))
    }
    if (file.exists(identity_path)) {
      v <- fread(identity_path)
      keep <- intersect(c("variant_id", "rsID", "GTEx_REF", "GTEx_ALT",
                          "r2_with_sentinel_EUR", "identity", "slope_ALT", "p_eQTL"), names(v))
      v <- unique(v[, ..keep])
      setnames(v, "identity", "Ensembl_identity", skip_absent = TRUE)
      # Keep alternative eQTL estimates separately and avoid overwriting audit values.
      overlap <- intersect(c("slope_ALT", "p_eQTL"), intersect(names(s2), names(v)))
      if (length(overlap)) setnames(v, overlap, paste0("P2b_", overlap))
      s2 <- merge(s2, v, by = "variant_id", all = TRUE)
      note("S2 Ensembl/LD source: ", identity_path)
    }
    harmonized <- if (nzchar(run)) file.path(run, "S2_eQTL_harmonized.tsv") else ""
    if (nzchar(harmonized) && file.exists(harmonized)) {
      h <- fread(harmonized)[gene == "HEATR3" & is.finite(risk_slope),
                             .(variant_id, harmonized_rsid = rsid, GWAS_risk_allele = risk_allele,
                               slope_GWAS_risk_allele = risk_slope)]
      s2 <- merge(s2, unique(h), by = "variant_id", all = TRUE, allow.cartesian = TRUE)
      note("S2 harmonized allele source: ", harmonized)
    }
    if (!nrow(s2)) stop("No HEATR3 variant rows were available.")
    note("S2 unique variants: ", uniqueN(s2$variant_id, na.rm = TRUE),
         "; rows with aligned risk-allele slope: ",
         sum(is.finite(val(s2, "slope_GWAS_risk_allele", NA_real_))))
    save(s2, "SuppTable_S2_GTEx_HEATR3_variant_identity",
         "**Supplementary Table S2.** GTEx v8 brain-cortex HEATR3 eQTL variant identity.",
         paste0("slope_ALT is oriented to the GTEx ALT allele. An rsID or LD r2 with ",
                SENTINEL, " does not establish GWAS risk-allele phase. Blank values were unavailable."))
  }, error = function(e) note("S2 pending: ", conditionMessage(e)))
  message("Output folder: ", OUT)
}

build_supp_tables()
