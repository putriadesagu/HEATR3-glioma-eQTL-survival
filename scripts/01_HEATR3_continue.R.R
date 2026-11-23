# RStudio / R 4.2.2. Save with UTF-8 encoding and run Source.
# Install once if needed: install.packages(c("survival", "readxl"))
# Use the original TCGA merged datasets for poster reproduction.
ROOT <- "D:/RESIDENSI/Thesis-PPDS/R"
# Input files from the project inventory.
LGG_FILE <- file.path(ROOT, "TCGA_LGG_merged_selected_genes_with_IDH1.csv")
GBM_FILE <- file.path(ROOT, "TCGA_GBM_merged.csv")
GWAS_FILE <- file.path(ROOT, "Glioblastoma_only_GWAS_Raw_Catalog_Export.xlsx")
GTEX_MATCHED_FILE <- file.path(
  ROOT, "FINAL_OUTPUT", "GTEx_Matching", "Curated_GWAS_GBM_with_GTEx.csv"
)
GTEX_PAIRS_FILE <- file.path(ROOT, "Brain_Cortex.v8.signif_variant_gene_pairs.txt.gz")
GTEX_EGENES_FILE <- file.path(ROOT, "data", "Brain_Cortex.v8.egenes.txt.gz")
# variant_id and rsid from the GTEx v8 variant lookup (rename columns as needed).
GTEX_VARIANT_LOOKUP <- file.path(ROOT, "data", "gtex_v8_variant_lookup.tsv.gz")
# Optional raw CGGA import. Use paths to decompressed .txt or .tsv files.
CGGA_RAW <- list(
  "693" = c(expression = "", clinical = ""),
  "325" = c(expression = "", clinical = "")
)
CGGA_COLUMNS <- c(patient_id = "CGGA_ID", os_days = "OS", event = "Censor",
                  age = "Age", sex = "Gender", type = "PRS_type", grade = "Grade")
CGGA_LOG2P1 <- FALSE  # Set to TRUE only for verified nonnegative, linear RSEM.
# Enter original column names here if detected names are ambiguous.
COLUMN_MAP <- list(LGG = c(), GBM = c())
# Example: COLUMN_MAP$GBM <- c(age = "age_at_diagnosis", event = "vital_status")
# Set units only after checking the original data dictionary or preprocessing code.
TCGA_AGE_UNIT <- c(LGG = "years", GBM = "days")  # years or days
TCGA_OS_UNIT <- c(LGG = "days", GBM = "days")    # days, months, or years
# Example only: TCGA_EVENT_MAP$GBM <- c("1" = 0, "2" = 1)
# Confirm which raw code means death before setting an event map.
TCGA_EVENT_MAP <- list(LGG = c(), GBM = c())
# Normalized supplementary inputs go in ROOT/data. Expression must already be
# on the specified analysis scale. OS is in days, event is 1=death / 0=censored.
# cns5.tsv: patient_id, eligible, h3_excluded, idh, codeletion, hist_gbm,
# tert_promoter, egfr_amp, chr7_gain, chr10_loss (flags 0/1 or NA).
# eligible=1 means adult diffuse glioma eligible for this classification;
# for IDH-wt, explicitly confirm diffuse astrocytic context and exclusions.
# molecular.tsv: patient_id and any of idh, grade, codeletion, MKI67, proliferation.
# proliferation is a prespecified continuous score excluding HEATR3.
# lgg_transcriptome.tsv: patient_id, then one numeric column per gene.
# heatr3_copy_number.tsv: patient_id, copy_number (one consistent assay/scale).
# gtex_matched.tsv: rsid, gene, variant_id, slope; optional pval_nominal.
# Build from GTEx v8 significant cortex pairs + verified GRCh38 rsID annotation
# and Ensembl-to-symbol mapping. variant_id is chr_pos_REF_ALT_b38.
# cgga_693.tsv and cgga_325.tsv: patient_id, os_days, event, age, sex, HEATR3,
# primary (1 initial primary tumor), hist_group (LGG for grades 2-3, GBM).
# These are joined, harmonized CGGA clinical + expression tables.
# Set TRUE only after checking patient overlap across the two CGGA releases.
CGGA_RELEASES_INDEPENDENT <- FALSE

run <- function() {
  if (!dir.exists(ROOT)) stop("Set ROOT to your R project folder.")
  if (!requireNamespace("survival", quietly = TRUE)) stop("Install survival.")
  library(survival)
  genes <- c("HEATR3", "LMF1", "RTEL1", "EGFR", "SLC16A8", "MIR4300HG")
  out <- file.path(ROOT, "results")
  dir.create(out, showWarnings = FALSE)
  # A fresh subdirectory prevents mixing outputs from different runs.
  out <- tempfile("HEATR3_", out)
  dir.create(out)
  writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
  notes <- character()
  note <- function(...) {
    line <- paste0(...)
    notes <<- c(notes, line)
    message(line)
  }
  save_table <- function(x, name) write.table(x, file.path(out, name),
    sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  on.exit(writeLines(notes, file.path(out, "MANUSCRIPT_NUMBERS.txt")))
  step <- function(name, work) tryCatch({
    work
    note(name, " completed")
  }, error = function(e) note(name, " pending: ", conditionMessage(e)))

  read_table <- function(path) {
    if (!length(path) || !nzchar(path) || !file.exists(path))
      stop("Missing input: ", paste(path, collapse = ", "))
    if (grepl("\\.xlsx?$", path, ignore.case = TRUE)) {
      if (!requireNamespace("readxl", quietly = TRUE)) stop("Install readxl.")
      return(as.data.frame(readxl::read_excel(path)))
    }
    if (grepl("\\.csv(\\.gz)?$", path, ignore.case = TRUE))
      return(read.csv(path, check.names = FALSE))
    read.delim(path, check.names = FALSE, quote = "", comment.char = "")
  }
  find_file <- function(explicit, stem) {
    if (nzchar(explicit)) return(explicit)
    hits <- list.files(ROOT, paste0("^", stem, "\\.(csv|tsv|txt|xlsx?)$"),
                       full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
    if (length(hits) != 1L) stop("See INPUT_CANDIDATES.txt. Expected one ", stem, " file; found ",
                                length(hits), ". Set the explicit path.")
    hits
  }
  extra <- function(name) read_table(file.path(ROOT, "data", name))
  need <- function(d, cols) {
    missing <- setdiff(cols, names(d))
    if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "))
  }
  num <- function(x) {
    y <- suppressWarnings(as.numeric(as.character(x)))
    if (any(is.na(y) & !is.na(x) & !trimws(x) %in% c("", "NA", "[Not Available]")))
      stop("Unexpected numeric values. Check input coding.")
    y
  }
  normalize_id <- function(x) {
    x <- trimws(as.character(x))
    tcga <- grepl("^TCGA[-.]", x, ignore.case = TRUE)
    x[tcga] <- substr(gsub("\\.", "-", toupper(x[tcga])), 1, 12)
    x
  }
  unique_id <- function(d) {
    need(d, "patient_id")
    d$patient_id <- normalize_id(d$patient_id)
    if (anyNA(d$patient_id) || any(!nzchar(d$patient_id)) || anyDuplicated(d$patient_id))
      stop("Expected one row per patient with a valid identifier.")
    d
  }
  attach <- function(d, e) {
    e <- unique_id(e)
    for (name in setdiff(names(e), names(d)))
      d[[name]] <- e[[name]][match(d$patient_id, e$patient_id)]
    d
  }
  normalize <- function(d, map = c(), cohort = NULL) {
    selected <- character()
    aliases <- list(patient_id = c("patient_id", "submitter_id", "bcr_patient_barcode", "patient", "sample"),
      os_days = c("os_days", "os_time", "OS.time", "OS_days", "survival_days"),
      event = c("event", "OS", "OS_event", "os_status", "vital_status", "status"),
      age = c("age", "age_at_diagnosis_years", "age_at_initial_pathologic_diagnosis", "Age_at_diagnosis"),
      sex = c("sex", "gender"), idh = c("idh", "IDH_status", "IDH1_status", "IDH1"))
    for (name in c(names(aliases), genes, "MKI67")) {
      candidates <- if (name %in% names(map)) map[[name]] else c(name, aliases[[name]])
      index <- match(tolower(candidates), tolower(names(d)))
      index <- index[!is.na(index)]
      if (name %in% names(map) && !length(index))
        stop("Mapped column for ", name, " was not found: ", map[[name]])
      if (length(index)) {
        selected[name] <- names(d)[index[1]]
        if (names(d)[index[1]] != name && name %in% names(d))
          d[[name]] <- d[[index[1]]]
        else names(d)[index[1]] <- name
      }
    }
    core <- c("patient_id", "os_days", "event", "age", "sex", "idh")
    if (!is.null(cohort)) {
      choices <- selected[intersect(core, names(selected))]
      note(cohort, " selected columns: ", paste(names(choices), choices, sep = " <- ", collapse = "; "))
      for (name in intersect(c("os_days", "age"), names(d))) {
        values <- suppressWarnings(as.numeric(as.character(d[[name]])))
        observed <- values[is.finite(values)]
        note(cohort, " ", name, " raw finite range=",
             if (length(observed)) paste(range(observed), collapse = " to ") else "none",
             "; negative=", sum(observed < 0), "; missing/non-numeric=", sum(is.na(values)))
      }
      if ("event" %in% names(d)) {
        counts <- sort(table(d$event, useNA = "ifany"), decreasing = TRUE)
        if (length(counts) <= 12) {
          note(cohort, " event raw counts: ",
               paste(ifelse(is.na(names(counts)), "<NA>", names(counts)), counts,
                     sep = "=", collapse = "; "))
        } else note(cohort, " event column has ", length(counts), " distinct values. Check mapping.")
      }
    }
    need(d, c("patient_id", "os_days", "event", "age", "sex"))
    d <- unique_id(d)
    e <- tolower(trimws(as.character(d$event)))
    if (!is.null(cohort) && length(TCGA_EVENT_MAP[[cohort]])) {
      lookup <- TCGA_EVENT_MAP[[cohort]]
      if (is.null(names(lookup)) || anyNA(lookup) || any(!lookup %in% 0:1) ||
          anyDuplicated(tolower(trimws(names(lookup)))))
        stop("TCGA_EVENT_MAP must be a named mapping to 0/1.")
      names(lookup) <- tolower(trimws(names(lookup)))
      observed <- unique(e[!is.na(e) & nzchar(e)])
      if (!all(observed %in% names(lookup)))
        stop("TCGA_EVENT_MAP does not cover all observed event codes.")
      e <- unname(lookup[e])
    }
    e[e %in% c("dead", "deceased", "1:deceased")] <- "1"
    e[e %in% c("alive", "living", "0:living")] <- "0"
    d$event <- num(e)
    d$os_days <- num(d$os_days)
    d$age <- num(d$age)
    if (!is.null(cohort)) {
      age_unit <- match.arg(TCGA_AGE_UNIT[[cohort]], c("years", "days"))
      os_unit <- match.arg(TCGA_OS_UNIT[[cohort]], c("days", "months", "years"))
      if (age_unit == "days") d$age <- d$age / 365.25
      d$os_days <- d$os_days * c(days = 1, months = 365.25 / 12, years = 365.25)[[os_unit]]
    }
    issues <- c(
      if (any(!is.na(d$event) & !d$event %in% 0:1)) "event has values outside 0/1",
      if (any(!is.na(d$os_days) & (!is.finite(d$os_days) | d$os_days < 0)))
        "OS days has negative or non-finite values",
      if (any(!is.na(d$age) & (!is.finite(d$age) | d$age > 120 | d$age < 0)))
        "age in years has values outside 0-120")
    if (length(issues)) stop(paste(issues, collapse = "; "))
    d$sex <- factor(tolower(trimws(d$sex)), levels = c("female", "male"))
    d
  }
  idh_factor <- function(x) {
    x <- tolower(trimws(as.character(x)))
    y <- ifelse(x %in% c("mut", "mutant", "mutated", "idh-mutant", "idh1-mutant", "idhmut"), "mut",
          ifelse(x %in% c("wt", "wildtype", "wild-type", "wild type", "idh-wildtype", "idh1-wildtype", "idhwt"), "wt", NA))
    factor(y, levels = c("mut", "wt"))
  }
  fit <- function(d, rhs) {
    f <- as.formula(paste("Surv(os_days, event) ~", rhs))
    m <- model.frame(f, d, na.action = na.omit)
    if (!nrow(m) || !sum(m[[1]][, 2])) stop("No complete cases with events.")
    d <- d[match(rownames(m), rownames(d)), , drop = FALSE]
    for (variable in all.vars(f[[3]])) {
      if (is.factor(d[[variable]])) d[[variable]] <- droplevels(d[[variable]])
      if (length(unique(d[[variable]])) < 2L) stop("No variation in ", variable)
    }
    model <- withCallingHandlers(coxph(f, data = d, x = TRUE, singular.ok = FALSE),
                                 warning = function(w) stop(conditionMessage(w)))
    if (any(!is.finite(coef(model)))) stop("Non-finite coefficient.")
    model
  }
  tidy <- function(model, label) {
    s <- summary(model)
    ph <- tryCatch(cox.zph(model, terms = FALSE)$table, error = function(e) NULL)
    p_ph <- if (is.null(ph)) rep(NA_real_, nrow(s$coefficients)) else
      ph[match(rownames(s$coefficients), rownames(ph)), "p"]
    data.frame(model = label, term = rownames(s$coefficients), n = model$n,
      events = model$nevent, beta = s$coefficients[, "coef"],
      se = s$coefficients[, "se(coef)"], HR = s$conf.int[, 1],
      lower = s$conf.int[, "lower .95"], upper = s$conf.int[, "upper .95"],
      p = s$coefficients[, "Pr(>|z|)"], c_index = unname(s$concordance[1]),
      PH_term = p_ph, PH_global = if (is.null(ph)) NA_real_ else ph["GLOBAL", "p"],
      row.names = NULL)
  }
  safe_fit <- function(d, rhs, label) tryCatch(tidy(fit(d, rhs), label),
    error = function(e) {
      note(label, " unavailable: ", conditionMessage(e))
      NULL
    })
  combine <- function(xs) {
    ans <- do.call(rbind, xs)
    if (is.null(ans) || !nrow(ans)) stop("No estimable results.")
    ans
  }
  headline <- function(x, label, term = "expression") {
    x <- x[x$term == term, , drop = FALSE]
    for (i in seq_len(nrow(x))) {
      z <- x[i, ]
      note(sprintf("%s | %s | n=%d events=%d HR=%.2f (%.2f-%.2f), p=%.4g, C=%.3f",
        label, z$model, z$n, z$events, z$HR, z$lower, z$upper, z$p, z$c_index))
    }
  }

  input_files <- list.files(ROOT, recursive = TRUE)
  candidates <- input_files[grepl("lgg|low.?grade|gwas|association", basename(input_files),
                                  ignore.case = TRUE) &
                              grepl("\\.(csv|tsv|txt|xlsx?)(\\.gz)?$", input_files,
                                    ignore.case = TRUE) &
                              !grepl("^results[/\\\\]", input_files)]
  writeLines(candidates, file.path(out, "INPUT_CANDIDATES.txt"))
  note("Input filenames for LGG/GWAS: see INPUT_CANDIDATES.txt.")
  gwas <- NULL
  step("S1 GWAS", {
    if (!nzchar(GWAS_FILE)) stop("Set GWAS_FILE to the glioblastoma Associations TSV.")
    input <- read_table(GWAS_FILE)
    save_table(data.frame(column = names(input)), "INPUT_HEADERS_GWAS.tsv")
    need(input, c("P-VALUE", "REPORTED GENE(S)"))
    gwas <- input
    save_table(gwas, "S1_GWAS_input.tsv")
    sig <- gwas[which(num(gwas[["P-VALUE"]]) < 5e-8), , drop = FALSE]
    symbols <- trimws(unlist(strsplit(as.character(sig[["REPORTED GENE(S)"]]), "[,;]")))
    symbols <- sort(unique(symbols[!is.na(symbols) & !toupper(symbols) %in%
      c("", "NR", "NA", "INTERGENIC", "NOT REPORTED", "-")]))
    save_table(sig, "S1_GWAS_significant_associations.tsv")
    save_table(data.frame(reported_gene = symbols), "S1_GWAS_reported_genes.tsv")
    note("GWAS input rows=", nrow(gwas), "; P<5e-8 rows=", nrow(sig),
         "; unique reported gene labels=", length(symbols))
    if ("SNPS" %in% names(sig)) {
      ids <- unique(unlist(regmatches(as.character(sig$SNPS),
                                    gregexpr("rs[0-9]+", as.character(sig$SNPS)))))
      note("Unique lead rsIDs=", length(ids), " (not LD-independent loci).")
    }
  })
  step("S2 GTEx", {
    if (is.null(gwas)) stop("Load the GWAS export in S1 first.")
    need(gwas, c("P-VALUE", "STRONGEST SNP-RISK ALLELE"))
    cached <- file.path(ROOT, "data", "gtex_matched.tsv")
    if (nzchar(GTEX_MATCHED_FILE) || file.exists(cached)) {
      e <- read_table(if (nzchar(GTEX_MATCHED_FILE)) GTEX_MATCHED_FILE else cached)
      save_table(data.frame(column = names(e)), "INPUT_HEADERS_GTEX.tsv")
      curated_cols <- c("Sentinel_rsID", "gene_clean", "top_variant_id", "top_slope")
      if (all(curated_cols %in% names(e))) {
        curated <- unique(data.frame(
          sentinel_rsid = tolower(trimws(e$Sentinel_rsID)),
          gene = trimws(e$gene_clean),
          variant_id = trimws(e$top_variant_id),
          slope = num(e$top_slope),
          pval_nominal = if ("top_pval" %in% names(e)) num(e$top_pval) else rep(NA_real_, nrow(e))
        ))
        curated <- curated[!is.na(curated$variant_id) & nzchar(curated$variant_id) &
                           !is.na(curated$gene) & nzchar(curated$gene), , drop = FALSE]
        curated$identity <- rep("unverified", nrow(curated))
        if (file.exists(GTEX_VARIANT_LOOKUP)) {
          lookup <- read_table(GTEX_VARIANT_LOOKUP)
          need(lookup, c("variant_id", "rsid"))
          lookup$rsid <- tolower(trimws(lookup$rsid))
          known <- curated$variant_id %in% lookup$variant_id
          matches <- paste(curated$variant_id, curated$sentinel_rsid) %in%
                     paste(lookup$variant_id, lookup$rsid)
          curated$identity[known] <- "different_rsid"
          curated$identity[matches & grepl("^rs[0-9]+$", curated$sentinel_rsid)] <- "verified_same_rsid"
        }
        save_table(curated, "S2_GTEx_curated_top_eQTL.tsv")
        e <- curated[curated$identity == "verified_same_rsid", , drop = FALSE]
        note("S2 curated variant/gene records=", nrow(curated),
             "; top eQTL verified as sentinel=", nrow(e),
             "; different rsID=", sum(curated$identity == "different_rsid"),
             "; unverified=", sum(curated$identity == "unverified"))
        if (!nrow(e))
          stop("Curated summary saved. Allele direction needs a GTEx v8 variant_id-to-rsid lookup or verified matched pairs.")
        e$rsid <- e$sentinel_rsid
      }
    } else {
      e <- read_table(GTEX_PAIRS_FILE)
      annotation <- read_table(GTEX_EGENES_FILE)
      lookup <- read_table(GTEX_VARIANT_LOOKUP)
      need(e, c("gene_id", "variant_id", "slope", "pval_nominal"))
      need(annotation, c("gene_id", "gene_name"))
      need(lookup, c("variant_id", "rsid"))
      annotation <- unique(annotation[, c("gene_id", "gene_name")])
      lookup <- unique(lookup[, c("variant_id", "rsid")])
      e <- merge(e[, c("gene_id", "variant_id", "slope", "pval_nominal")],
                 lookup, by = "variant_id")
      e <- merge(e, annotation, by = "gene_id")
      names(e)[names(e) == "gene_name"] <- "gene"
    }
    need(e, c("rsid", "gene", "variant_id", "slope"))
    e$rsid <- tolower(trimws(e$rsid))
    e$slope <- num(e$slope)
    sig <- gwas[which(num(gwas[["P-VALUE"]]) < 5e-8), , drop = FALSE]
    key <- toupper(trimws(sig[["STRONGEST SNP-RISK ALLELE"]]))
    valid <- which(grepl("^RS[0-9]+-[ACGT]$", key))
    risks <- data.frame(rsid = tolower(sub("-.*", "", key[valid])),
      risk_allele = sub(".*-", "", key[valid]), gwas_row = valid,
      gwas_p = sig[["P-VALUE"]][valid])
    save_table(sig[setdiff(seq_len(nrow(sig)), valid), , drop = FALSE],
               "S2_GWAS_allele_unresolved.tsv")
    h <- merge(risks, e, by = "rsid")
    h$ref <- sub("^[^_]+_[0-9]+_([^_]+)_.*$", "\\1", h$variant_id)
    h$alt <- sub("^[^_]+_[0-9]+_[^_]+_([^_]+)_.*$", "\\1", h$variant_id)
    simple <- grepl("^[^_]+_[0-9]+_[ACGT]_[ACGT]_b38$", h$variant_id)
    ambiguous <- paste0(h$ref, h$alt) %in% c("AT", "TA", "CG", "GC")
    h$risk_slope <- ifelse(simple & !ambiguous & h$risk_allele == h$alt, h$slope,
                    ifelse(simple & !ambiguous & h$risk_allele == h$ref, -h$slope, NA))
    h$alignment <- ifelse(is.na(h$risk_slope), "unresolved",
                          ifelse(h$risk_allele == h$alt, "ALT", "REF"))
    save_table(h, "S2_eQTL_harmonized.tsv")
    z <- unique(h[h$gene == "HEATR3" & !is.na(h$risk_slope),
                  c("rsid", "risk_allele", "risk_slope")])
    note("HEATR3 harmonized GWAS allele-eQTL records=", nrow(z),
         "; higher expression=", sum(z$risk_slope > 0),
         "; lower expression=", sum(z$risk_slope < 0))
  })

  cohorts <- list()
  for (cohort in c("LGG", "GBM")) {
    step(paste("Read TCGA", cohort), {
      explicit <- if (cohort == "LGG") LGG_FILE else GBM_FILE
      path <- find_file(explicit, paste0("TCGA_", cohort, "_merged"))
      raw <- read_table(path)
      save_table(data.frame(column = names(raw)), paste0("INPUT_HEADERS_", cohort, ".tsv"))
      d <- normalize(raw, COLUMN_MAP[[cohort]], cohort)
      for (gene in intersect(c(genes, "MKI67"), names(d))) d[[gene]] <- num(d[[gene]])
      cohorts[[cohort]] <- d
      expected <- if (cohort == "LGG") 512L else 268L
      note(cohort, " source=", basename(path), "; rows=", nrow(d),
           "; expected poster rows=", expected, "; events=", sum(d$event, na.rm = TRUE))
    })
  }

  step("S4 Cox Table 1", {
    rows <- list()
    for (cohort in names(cohorts)) for (gene in genes) {
      d <- cohorts[[cohort]]
      for (adjusted in c(FALSE, TRUE)) {
        label <- paste(cohort, gene, if (adjusted) "age+sex" else "unadjusted")
        if (gene %in% names(d)) {
          d$expression <- d[[gene]]
          ans <- safe_fit(d, if (adjusted) "expression + age + sex" else "expression", label)
        } else {
          note(label, " unavailable: expression column missing.")
          ans <- NULL
        }
        if (is.null(ans)) ans <- data.frame(model = label, term = "expression",
          n = NA_integer_, events = NA_integer_, beta = NA_real_, se = NA_real_,
          HR = NA_real_, lower = NA_real_, upper = NA_real_, p = NA_real_,
          c_index = NA_real_, PH_term = NA_real_, PH_global = NA_real_)
        ans$cohort <- cohort
        ans$gene <- gene
        ans$adjusted <- adjusted
        rows[[length(rows) + 1]] <- ans
      }
    }
    tab <- combine(rows)
    tab$FDR <- NA_real_
    for (cohort in names(cohorts)) for (adjusted in c(FALSE, TRUE)) {
      idx <- which(tab$cohort == cohort & tab$adjusted == adjusted & tab$term == "expression")
      tab$FDR[idx] <- p.adjust(tab$p[idx], "BH", n = length(genes))
    }
    save_table(tab, "S4_Cox_Table1.tsv")
    headline(tab, "S4")
    if (all(c("LGG", "GBM") %in% names(cohorts)) &&
        all(vapply(cohorts, function(d) "HEATR3" %in% names(d), logical(1)))) {
      pooled <- combine(lapply(names(cohorts), function(cohort) {
        d <- cohorts[[cohort]][, c("patient_id", "os_days", "event", "HEATR3", "age", "sex")]
        d$cohort <- cohort
        d$expression <- as.numeric(scale(d$HEATR3))
        d
      }))
      pooled <- unique_id(pooled)
      pooled$cohort <- factor(pooled$cohort, levels = c("LGG", "GBM"))
      ans <- safe_fit(pooled, "expression + expression:cohort + age + sex + strata(cohort)",
                      "HEATR3 per SD by legacy TCGA cohort")
      if (!is.null(ans)) {
        save_table(ans, "S4_HEATR3_cohort_interaction.tsv")
        headline(ans, "S4 cohort interaction", "expression:cohortGBM")
      }
    }
  })

  lgg <- NULL
  step("Prepare LGG models", {
    if (is.null(cohorts$LGG)) stop("LGG merged input unavailable.")
    d <- cohorts$LGG
    need(d, "HEATR3")
    path <- file.path(ROOT, "data", "molecular.tsv")
    if (file.exists(path)) {
      supplement <- unique_id(read_table(path))
      for (name in intersect(setdiff(names(supplement), "patient_id"), names(d))) {
        values <- supplement[[name]][match(d$patient_id, supplement$patient_id)]
        gap <- is.na(d[[name]]) | trimws(as.character(d[[name]])) == ""
        d[[name]][gap] <- values[gap]
      }
      d <- attach(d, supplement)
    }
    d$expression <- d$HEATR3
    if ("idh" %in% names(d)) {
      d$idh <- idh_factor(d$idh)
      note("LGG IDH recognized: mutant=", sum(d$idh == "mut", na.rm = TRUE),
           "; wild-type=", sum(d$idh == "wt", na.rm = TRUE),
           "; missing or unmapped=", sum(is.na(d$idh)))
    }
    lgg <- d
  })
  step("S5 IDH and proliferation", {
    need(lgg, "idh")
    note("Slide benchmark: IDH+age model n=512, events=126, HR=1.6, p=0.025, C=0.80.")
    rows <- list(
      safe_fit(lgg, "expression + idh + age", "Poster specification"),
      safe_fit(lgg, "expression + idh + age + sex", "IDH age sex"),
      safe_fit(lgg, "expression * idh + age + sex", "HEATR3 by IDH interaction")
    )
    for (status in c("mut", "wt")) {
      sub <- lgg[which(lgg$idh == status), , drop = FALSE]
      rows <- c(rows, list(safe_fit(sub, "expression + age + sex", paste("IDH", status))))
    }
    for (v in list("grade", "codeletion", "MKI67", "proliferation",
                   c("grade", "codeletion", "MKI67"),
                   c("grade", "codeletion", "proliferation"))) {
      label <- paste(v, collapse = " + ")
      if (!all(v %in% names(lgg))) {
        note("S5 ", label, " missing input: ", paste(setdiff(v, names(lgg)), collapse = ", "))
        next
      }
      cols <- c("expression", "age", "sex", "idh", v)
      sub <- lgg
      for (name in intersect(v, c("grade", "codeletion")))
        sub[[name]] <- factor(sub[[name]])
      for (name in intersect(v, c("MKI67", "proliferation")))
        sub[[name]] <- num(sub[[name]])
      sub <- sub[complete.cases(sub[, c("os_days", "event", cols)]), , drop = FALSE]
      rows <- c(rows, list(safe_fit(sub, "expression + idh + age + sex",
                                  paste("Matched baseline for", label)),
                           safe_fit(sub, paste("expression + idh + age + sex +", label),
                                  paste("Adjusted for", label))))
    }
    tab <- combine(rows)
    save_table(tab, "S5_HEATR3_models.tsv")
    headline(tab, "S5")
    headline(tab, "S5 interaction", "expression:idhwt")
  })

  km <- function(d, four, name, cut) {
    d <- d[complete.cases(d[, c("os_days", "event", "expression")]), ]
    d$group <- factor(ifelse(d$expression > cut, "High", "Low"), levels = c("Low", "High"))
    if (four) d$group <- interaction(d$idh, d$group, drop = TRUE, sep = " / ")
    d <- d[!is.na(d$group), , drop = FALSE]
    d$group <- droplevels(d$group)
    if (nlevels(d$group) < 2) stop("Fewer than two KM groups.")
    d$years <- d$os_days / 365.25
    f <- Surv(years, event) ~ group
    s <- survfit(f, data = d)
    logrank <- survdiff(f, data = d)
    p <- pchisq(logrank$chisq, length(logrank$n) - 1, lower.tail = FALSE)
    ticks <- seq(0, floor(max(d$years)), length.out = min(6, floor(max(d$years)) + 1))
    times <- sort(unique(round(ticks, 1)))
    risk <- sapply(times, function(t) vapply(levels(d$group),
      function(g) sum(d$group == g & d$years >= t), integer(1)))
    risk <- matrix(risk, nrow = nlevels(d$group), dimnames = list(levels(d$group), times))
    save_table(data.frame(group = rownames(risk), risk, check.names = FALSE),
               paste0(name, "_at_risk.tsv"))
    pdf(file.path(out, paste0(name, ".pdf")), width = 8, height = 6)
    on.exit(dev.off())
    layout(matrix(1:2, ncol = 1), heights = c(4, 1.5))
    par(mar = c(4, 5, 3, 2))
    palette <- c("#007C91", "#C74343", "#364C86", "#A66A22")
    plot(s, col = palette, lwd = 2, mark.time = TRUE, xlab = "Years", ylab = "Overall survival",
         main = paste(name, "| log-rank p =", format.pval(p, digits = 3)))
    legend("topright", levels(d$group), col = palette, lwd = 2, bty = "n")
    par(mar = c(1, 5, 2, 2))
    plot.new()
    plot.window(xlim = c(-1.8, length(times) + 0.5), ylim = c(0.3, nrow(risk) + 1.3))
    text(seq_along(times), nrow(risk) + 1, times)
    text(0, nrow(risk) + 1, "Years", adj = 1)
    for (i in seq_len(nrow(risk))) {
      y <- nrow(risk) + 1 - i
      text(0, y, rownames(risk)[i], adj = 1, col = palette[i])
      text(seq_along(times), y, risk[i, ])
    }
    title("Number at risk", line = 0.5)
    note(name, " n=", nrow(d), "; median cutoff=", signif(cut, 5), "; log-rank p=", signif(p, 5))
  }
  if (!is.null(lgg)) {
    cut <- median(lgg$expression[complete.cases(lgg[, c("os_days", "event", "expression")])])
    step("S6 LGG", km(lgg, FALSE, "S6_LGG", cut))
    for (status in c("wt", "mut")) {
      step(paste("S6 IDH", status), {
        need(lgg, "idh")
        km(lgg[which(lgg$idh == status), ], FALSE, paste0("S6_IDH_", status), cut)
      })
    }
    step("S6 four groups", {
      need(lgg, "idh")
      km(lgg, TRUE, "S6_four_groups", cut)
    })
  } else note("S6 pending: LGG HEATR3 input unavailable.")

  step("S7 Transcriptome Cox", {
    need(lgg, "HEATR3")
    wide <- unique_id(extra("lgg_transcriptome.tsv"))
    idx <- match(lgg$patient_id, wide$patient_id)
    if (!any(!is.na(idx))) stop("No matched LGG patients.")
    candidates <- setdiff(names(wide), "patient_id")
    if (!"HEATR3" %in% candidates) stop("Transcriptome needs a HEATR3 symbol column.")
    base <- lgg[, c("os_days", "event", "age", "sex")]
    rows <- lapply(candidates, function(gene) {
      base$expression <- num(wide[[gene]][idx])
      result <- tryCatch({
        m <- fit(base, "expression + age + sex")
        if (m$nevent < 10) stop("Fewer than 10 events")
        c(n = m$n, events = m$nevent, p = summary(m)$coefficients["expression", "Pr(>|z|)"])
      }, error = function(e) c(n = NA_real_, events = NA_real_, p = NA_real_))
      data.frame(gene = gene, t(result), row.names = NULL)
    })
    tab <- combine(rows)
    tab$FDR <- p.adjust(tab$p, "BH")
    tab$rank <- rank(tab$p, ties.method = "min", na.last = "keep")
    tested <- sum(!is.na(tab$p))
    tab$top_percent <- 100 * tab$rank / tested
    tab$percentile <- 100 * (tested - tab$rank + 1) / tested
    save_table(tab[order(tab$rank), ], "S7_Transcriptome_ranks.tsv")
    h <- tab[tab$gene == "HEATR3", ]
    note("S7 HEATR3 rank=", h$rank, "/", tested, "; top %=", signif(h$top_percent, 4),
         "; percentile=", signif(h$percentile, 4), "; FDR=", signif(h$FDR, 4))
  })

  step("S8 WHO CNS5 sensitivity analysis", {
    markers <- unique_id(extra("cns5.tsv"))
    flags <- c("eligible", "h3_excluded", "codeletion", "hist_gbm", "tert_promoter",
               "egfr_amp", "chr7_gain", "chr10_loss")
    need(markers, c(flags, "idh"))
    markers[flags] <- lapply(markers[flags], num)
    if (any(!unlist(markers[flags]) %in% c(0, 1, NA))) stop("CNS5 flags must be 0/1/NA.")
    markers$idh <- idh_factor(markers$idh)
    markers$entity <- "Unresolved"
    eligible <- markers$eligible == 1
    mut <- eligible & markers$idh == "mut"
    wt <- eligible & markers$idh == "wt" & markers$h3_excluded == 1
    molecular <- markers$tert_promoter == 1 | markers$egfr_amp == 1 |
      (markers$chr7_gain == 1 & markers$chr10_loss == 1)
    markers$entity[which(mut & markers$codeletion == 1)] <- "Oligodendroglioma_IDHmut"
    markers$entity[which(mut & markers$codeletion == 0)] <- "Astrocytoma_IDHmut"
    markers$entity[which(wt & (markers$hist_gbm == 1 | molecular))] <- "GBM_IDHwt"
    markers$entity[which(wt & markers$hist_gbm == 0 & !molecular)] <- "IDHwt_no_GBM_criteria"
    pooled <- combine(lapply(names(cohorts), function(cohort) {
      d <- cohorts[[cohort]]
      need(d, "HEATR3")
      d <- d[, c("patient_id", "os_days", "event", "age", "sex", "HEATR3")]
      d$legacy <- cohort
      d
    }))
    pooled <- unique_id(pooled)
    pooled <- attach(pooled, markers)
    pooled$entity[is.na(pooled$entity)] <- "Unresolved"
    save_table(as.data.frame(table(pooled$legacy, pooled$entity)), "S8_reclassification_counts.tsv")
    rows <- lapply(setdiff(unique(pooled$entity), "Unresolved"), function(entity) {
      d <- pooled[pooled$entity == entity, ]
      d$expression <- as.numeric(scale(d$HEATR3))
      safe_fit(d, "expression + age + sex", entity)
    })
    tab <- combine(rows)
    save_table(tab, "S8_HEATR3_by_entity.tsv")
    headline(tab, "S8 per SD")
  })

  step("S9 Copy number", {
    cn <- unique_id(extra("heatr3_copy_number.tsv"))
    need(cn, "copy_number")
    cn$copy_number <- num(cn$copy_number)
    rows <- lapply(names(cohorts), function(cohort) {
      d <- cohorts[[cohort]]
      if (!"HEATR3" %in% names(d)) return(NULL)
      d <- attach(d[, c("patient_id", "HEATR3")], cn)
      d <- d[complete.cases(d[, c("HEATR3", "copy_number")]), ]
      if (nrow(d) < 3L || length(unique(d$HEATR3)) < 2L ||
          length(unique(d$copy_number)) < 2L) return(NULL)
      test <- cor.test(d$HEATR3, d$copy_number, method = "spearman", exact = FALSE)
      note("S9 ", cohort, " n=", nrow(d), "; rho=", signif(test$estimate, 4),
           "; p=", signif(test$p.value, 4))
      data.frame(cohort = cohort, n = nrow(d), rho = unname(test$estimate), p = test$p.value)
    })
    save_table(combine(rows), "S9_copy_number_expression.tsv")
  })

  step("S10 External replication", {
    rows <- list()
    for (release in c("693", "325")) {
      step(paste("S10 CGGA", release), {
        paths <- CGGA_RAW[[release]]
        if (all(nzchar(paths))) {
          raw <- read_table(paths[["clinical"]])
          need(raw, unname(CGGA_COLUMNS))
          d <- raw[, unname(CGGA_COLUMNS)]
          names(d) <- names(CGGA_COLUMNS)
          d$primary <- as.integer(tolower(trimws(d$type)) == "primary")
          grade <- toupper(trimws(d$grade))
          d$hist_group <- ifelse(grade %in% c("II", "III", "2", "3", "WHO II", "WHO III"),
                                 "LGG", ifelse(grade %in% c("IV", "4", "WHO IV"), "GBM", NA))
          matrix <- read_table(paths[["expression"]])
          h <- which(trimws(matrix[[1]]) == "HEATR3")
          if (length(h) != 1L) stop("Expected one HEATR3 row in the CGGA expression matrix.")
          columns <- match(d$patient_id, names(matrix)[-1])
          values <- num(unlist(matrix[h, -1, drop = FALSE], use.names = FALSE))
          d$HEATR3 <- values[columns]
          if (CGGA_LOG2P1) {
            if (any(d$HEATR3 < 0, na.rm = TRUE)) stop("Negative RSEM values.")
            d$HEATR3 <- log2(d$HEATR3 + 1)
          }
          d <- normalize(d)
        } else {
          d <- normalize(extra(paste0("cgga_", release, ".tsv")))
        }
        need(d, c("primary", "hist_group", "HEATR3"))
        d <- d[which(num(d$primary) == 1 & toupper(d$hist_group) == "LGG"), ]
        d$HEATR3 <- num(d$HEATR3)
        d <- d[complete.cases(d[, c("os_days", "event", "age", "sex", "HEATR3")]), ]
        d$expression <- as.numeric(scale(d$HEATR3))
        rows[[release]] <- safe_fit(d, "expression + age + sex", paste0("CGGA_", release))
      })
    }
    if (!is.null(lgg)) {
      d <- lgg[complete.cases(lgg[, c("os_days", "event", "age", "sex", "HEATR3")]), ]
      d$expression <- as.numeric(scale(d$HEATR3))
      rows$TCGA <- safe_fit(d, "expression + age + sex", "TCGA_LGG")
    }
    tab <- combine(rows)
    tab <- tab[tab$term == "expression", ]
    save_table(tab, "S10_HEATR3_per_SD.tsv")
    headline(tab, "S10 per SD")
    if (!"TCGA_LGG" %in% tab$model || !any(grepl("^CGGA", tab$model)))
      stop("At least one estimable CGGA set and TCGA LGG are needed for meta-analysis.")
    if (sum(grepl("^CGGA", tab$model)) > 1L && !CGGA_RELEASES_INDEPENDENT)
      stop("Pooled estimate pending verification of CGGA release independence.")
    w <- 1 / tab$se^2
    fixed <- sum(w * tab$beta) / sum(w)
    q <- sum(w * (tab$beta - fixed)^2)
    df <- nrow(tab) - 1
    i2 <- if (q > 0) max(0, (q - df) / q) * 100 else 0
    tau2 <- max(0, (q - df) / (sum(w) - sum(w^2) / sum(w)))
    meta <- function(weights, label) {
      beta <- sum(weights * tab$beta) / sum(weights)
      se <- sqrt(1 / sum(weights))
      data.frame(method = label, k = nrow(tab), HR = exp(beta),
        lower = exp(beta - 1.96 * se), upper = exp(beta + 1.96 * se),
        p = 2 * pnorm(-abs(beta / se)), Q = q, I2 = i2, tau2_DL = tau2)
    }
    result <- rbind(meta(w, "fixed inverse variance"),
                    meta(1 / (tab$se^2 + tau2), "random DL sensitivity"))
    save_table(result, "S10_meta_analysis.tsv")
    note(sprintf("S10 fixed effect k=%d HR=%.2f (%.2f-%.2f), p=%.4g, I2=%.1f%%",
      result$k[1], result$HR[1], result$lower[1], result$upper[1], result$p[1], result$I2[1]))
  })
  note("S11 sessionInfo.txt and MANUSCRIPT_NUMBERS.txt recorded in ", out)
  invisible(out)
}

run()
