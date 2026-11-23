###############################################################################
#  HEATR3 / glioma — STRENGTHENING ANALYSES FOR THE MANUSCRIPT
#
#  A1  Proliferation: does HEATR3 survive adjustment for a proliferation
#      score and MKI67?  (TCGA-LGG + CGGA-693, then meta-analysis)
#  A2  Genome-wide null: where does HEATR3 rank among ALL genes in each
#      cohort, and among genes that replicate in both?
#  A3  Risk-allele harmonisation: does the GWAS risk allele at 16q12.1
#      INCREASE HEATR3 expression in GTEx brain cortex?
#  B7  Added prognostic value: does HEATR3 add to age + sex + IDH (+ grade,
#      1p/19q)?  LR test, delta C-index, time-dependent AUC at 1/3/5 years
#  C9  Pathways: which pathways co-vary with HEATR3 (IDH-adjusted)?  GSEA
#
# INPUT: native TCGA-LGG and CGGA-693 project files configured below.
# install.packages(c("survival", "data.table", "ggplot2", "readxl", "timeROC"))
# Optional GSEA: install.packages("msigdbr"); BiocManager::install("fgsea")
# A2 and B7 constant-HR estimates are exploratory given TCGA non-PH.
# A1 additionally fits time-varying HEATR3 and age, reporting HR at 1 year.

suppressPackageStartupMessages({
  library(survival); library(data.table); library(ggplot2)
})

## ============================== CONFIG =====================================
ROOT <- "D:/RESIDENSI/Thesis-PPDS/R"
if (!dir.exists(ROOT)) stop("Set ROOT to your project folder.")
DATA_DIR <- file.path(ROOT, "data")
OUT <- file.path(ROOT, "results_strengthen", format(Sys.time(), "%Y%m%d_%H%M%S"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
locate <- function(filename) {
  hits <- list.files(ROOT, recursive = TRUE, full.names = TRUE)
  hits <- hits[basename(hits) == filename]
  if (length(hits) == 1L) hits else file.path(DATA_DIR, filename)
}

COHORTS <- list(
  TCGA_LGG = list(clin = file.path(ROOT, "TCGA_LGG_merged_selected_genes_with_IDH1.csv"),
                  expr = locate("LGG_expression_xena.tsv.gz")),
  CGGA_693 = list(clin = locate("CGGA.mRNAseq_693_clinical.20200506.txt"),
                  expr = locate("CGGA.mRNAseq_693.RSEM-genes.20200506.txt"))
)
CGGA_LOG2P1_FOR_SCREEN <- TRUE  # CGGA STAR+RSEM FPKM; applies to this script only.
TCGA_EXPR_LOG2P1 <- FALSE       # Set only if verified Xena file contains linear expression.
# The TCGA selected-gene CSV supplies the exact HEATR3 used in the original models.

GENE         <- "HEATR3"
PROLIF_GENES <- c("MKI67", "TOP2A", "CCNB1", "CDK1", "BUB1", "PCNA", "MCM2",
                  "AURKA", "CCNA2", "TPX2")          # proliferation score
MIN_MEAN_EXPR <- 1      # genome-wide filter: mean log2 expression
MIN_SD_EXPR   <- 0.1
TIMES_DAYS    <- c(1, 3, 5) * 365.25

# A3 uses the existing S2 identity audit and verified harmonization output.
S2_RESULTS <- "" # Optionally set the exact TCGA results subdirectory.

set.seed(2026)

## ============================== HELPERS ====================================
summary_lines <- character()
note <- function(...) { s <- paste0(...); message(s); summary_lines <<- c(summary_lines, s) }
fmt_hr <- function(hr, lo, hi) sprintf("%.2f (%.2f–%.2f)", hr, lo, hi)
fmt_p  <- function(p) ifelse(is.na(p), "NA", ifelse(p < 0.001, formatC(p, format = "e", digits = 2), sprintf("%.3f", p)))
zs     <- function(x) as.numeric(scale(x))

canonical_id <- function(x) {
  x <- trimws(as.character(x))
  tcga <- grepl("^TCGA[-.]", x, ignore.case = TRUE)
  x[tcga] <- substr(gsub("\\.", "-", toupper(x[tcga])), 1L, 12L)
  x
}

read_expr <- function(f, nm) {
  if (grepl("\\.gz$", f, ignore.case = TRUE)) {
    con <- gzfile(f, "rt")
    tmp <- tempfile(fileext = ".tsv")
    on.exit({ close(con); unlink(tmp) })
    out <- file(tmp, "wt")
    repeat {
      lines <- readLines(con, n = 1000L, warn = FALSE)
      if (!length(lines)) break
      writeLines(lines, out)
    }
    close(out)
    f <- tmp
  }
  d <- fread(f, data.table = FALSE, check.names = FALSE)
  raw_ids <- names(d)[-1L]
  if (nm == "TCGA_LGG") {
    barcodes <- gsub("\\.", "-", toupper(raw_ids))
    primary <- grepl("^TCGA-.{2}-.{4}-01", barcodes)
    if (any(primary)) {
      keep <- primary | nchar(barcodes) == 12L
      d <- d[, c(TRUE, keep), drop = FALSE]
    }
  }
  m <- as.matrix(d[, -1L, drop = FALSE])
  storage.mode(m) <- "numeric"
  rownames(m) <- toupper(trimws(as.character(d[[1L]])))
  colnames(m) <- canonical_id(names(d)[-1L])
  if (anyDuplicated(colnames(m))) stop("Multiple expression samples per patient. Resolve barcodes before analysis.")
  genes <- rownames(m)
  unique_symbol <- nzchar(genes) & !is.na(genes) & !duplicated(genes) &
    !duplicated(genes, fromLast = TRUE)
  m <- m[unique_symbol, , drop = FALSE]
  if ((nm == "CGGA_693" && CGGA_LOG2P1_FOR_SCREEN) ||
      (nm == "TCGA_LGG" && TCGA_EXPR_LOG2P1)) {
    if (any(m < 0, na.rm = TRUE)) stop("log2p1 requested for negative expression.")
    m <- log2(m + 1)
  }
  m
}

read_clin <- function(f, nm) {
  d <- fread(f, data.table = FALSE, check.names = FALSE)
  cols <- if (nm == "TCGA_LGG")
    c(id = "patient_id", OS_time = "OS_time", OS_event = "OS_event", age = "age_at_diagnosis",
      sex = "gender", IDH = "IDH1_status", HEATR3_original = "HEATR3")
  else c(id = "CGGA_ID", OS_time = "OS", OS_event = "Censor (alive=0; dead=1)",
         age = "Age", sex = "Gender", IDH = "IDH_mutation_status",
         grade = "Grade", codel = "1p19q_codeletion_status", type = "PRS_type")
  miss <- setdiff(unname(cols), names(d))
  if (length(miss)) stop("Missing clinical columns: ", paste(miss, collapse = ", "))
  z <- d[, unname(cols), drop = FALSE]
  names(z) <- names(cols)
  if (nm == "CGGA_693") {
    grade <- toupper(trimws(z$grade))
    z <- z[which(tolower(trimws(z$type)) == "primary" & grade %in% c("II", "III", "2", "3")), ]
  }
  z$id <- canonical_id(z$id)
  if (anyNA(z$id) || any(!nzchar(z$id)) || anyDuplicated(z$id)) stop("Clinical IDs must be valid and unique.")
  k <- gsub("[^a-z0-9]", "", tolower(z$IDH))
  z$IDH <- factor(ifelse(k %in% c("mutant", "idh1mutant", "mut"), "mutant",
                    ifelse(k %in% c("wildtype", "idh1wildtype", "wt"), "wildtype", NA)),
                  levels = c("mutant", "wildtype"))
  sx <- tolower(trimws(z$sex))
  sx[sx %in% "f"] <- "female"
  sx[sx %in% "m"] <- "male"
  z$sex <- factor(sx, levels = c("female", "male"))
  if ("grade" %in% names(z)) z$grade <- factor(z$grade)
  if ("codel" %in% names(z)) {
    k <- gsub("[^a-z0-9]", "", tolower(trimws(z$codel)))
    z$codel <- factor(ifelse(k %in% c("codel", "codeletion", "codeleted"), "codel",
                       ifelse(k %in% c("noncodel", "noncodeletion", "noncodeleted"), "non-codel", NA)))
  }
  for (name in c("OS_time", "OS_event", "age"))
    z[[name]] <- suppressWarnings(as.numeric(z[[name]]))
  if (any(!is.na(z$OS_event) & !z$OS_event %in% 0:1) ||
      any(z$OS_time < 0, na.rm = TRUE) || any(z$age < 0 | z$age > 120, na.rm = TRUE))
    stop("Verify OS (days), age (years) and event coding.")
  z[complete.cases(z[, c("OS_time", "OS_event", "age", "sex")]), , drop = FALSE]
}

usable_rows <- function(ex) {
  apply(ex, 1L, function(x) all(is.finite(x)) &&
          mean(x) > MIN_MEAN_EXPR && sd(x) > MIN_SD_EXPR)
}

cox_row <- function(fit, term) {
  s <- summary(fit); i <- grep(paste0("^", term, "$"), rownames(s$coefficients))
  data.frame(HR = s$conf.int[i, 1], lo = s$conf.int[i, 3], hi = s$conf.int[i, 4],
             p = s$coefficients[i, "Pr(>|z|)"], logHR = s$coefficients[i, "coef"],
             se = s$coefficients[i, "se(coef)"], n = s$n, events = s$nevent)
}

meta_fe <- function(b, se) {
  w <- 1 / se^2; bm <- sum(w * b) / sum(w); sm <- sqrt(1 / sum(w))
  Q <- sum(w * (b - bm)^2); df <- length(b) - 1
  c(HR = exp(bm), lo = exp(bm - 1.96 * sm), hi = exp(bm + 1.96 * sm),
    p = 2 * pnorm(-abs(bm / sm)), Q_p = if (df > 0) 1 - pchisq(Q, df) else NA,
    I2 = if (df > 0 && Q > 0) max(0, (Q - df) / Q) * 100 else 0)
}

## ============================== LOAD =======================================
DAT <- list(); EXPR <- list()
for (nm in names(COHORTS)) {
  inputs <- unlist(COHORTS[[nm]], use.names = FALSE)
  missing <- inputs[!file.exists(inputs)]
  if (length(missing)) { note("Skipping ", nm, "; missing: ", paste(missing, collapse = " | ")); next }
  note("Inputs ", nm, ": ", paste(inputs, collapse = " | "))
  cl <- read_clin(COHORTS[[nm]]$clin, nm)
  ex <- read_expr(COHORTS[[nm]]$expr, nm)
  ids <- intersect(cl$id, colnames(ex))
  if (!length(ids)) { note("No matched expression IDs: ", nm); next }
  cl <- cl[match(ids, cl$id), ]; ex <- ex[, ids, drop = FALSE]
  if (nm == "TCGA_LGG") {
    ex <- rbind(ex[rownames(ex) != GENE, , drop = FALSE], HEATR3 = cl$HEATR3_original)
  }
  if (!GENE %in% rownames(ex)) { note("HEATR3 missing: ", nm); next }
  ok <- is.finite(ex[GENE, ])
  cl <- cl[ok, , drop = FALSE]; ex <- ex[, ok, drop = FALSE]
  if (nrow(cl) < 2L || sd(ex[GENE, ]) == 0) {
    note("Insufficient HEATR3 variation: ", nm)
    next
  }
  cl$H3 <- zs(ex[GENE, ])                                   # HEATR3 per SD
  pg <- intersect(PROLIF_GENES, rownames(ex))
  pg <- pg[vapply(pg, function(g) all(is.finite(ex[g, ])) && sd(ex[g, ]) > 0, logical(1))]
  cl$prolif <- if (length(pg) >= 3L) zs(colMeans(t(scale(t(ex[pg, , drop = FALSE]))))) else NA_real_
  cl$MKI67z <- if ("MKI67" %in% rownames(ex)) zs(ex["MKI67", ]) else NA
  DAT[[nm]] <- cl; EXPR[[nm]] <- ex
  note(sprintf("LOAD %s: n=%d, deaths=%d, IDH mut/wt/NA = %d/%d/%d, proliferation genes used: %s",
               nm, nrow(cl), sum(cl$OS_event), sum(cl$IDH %in% "mutant"), sum(cl$IDH %in% "wildtype"),
               sum(is.na(cl$IDH)), paste(pg, collapse = ",")))
}
if (!length(DAT)) stop("No cohort loaded — check CONFIG paths")

###############################################################################
## A1. PROLIFERATION ADJUSTMENT
###############################################################################
a1 <- list()
for (nm in names(DAT)) {
  d <- DAT[[nm]][!is.na(DAT[[nm]]$IDH), ]
  forms <- list(
    base           = Surv(OS_time, OS_event) ~ H3 + age + sex + strata(IDH),
    plus_prolif    = Surv(OS_time, OS_event) ~ H3 + prolif + age + sex + strata(IDH),
    plus_MKI67     = Surv(OS_time, OS_event) ~ H3 + MKI67z + age + sex + strata(IDH),
    prolif_only    = Surv(OS_time, OS_event) ~ prolif + age + sex + strata(IDH))
  if (all(c("grade", "codel") %in% names(d)))
    forms$plus_prolif_grade_codel <- Surv(OS_time, OS_event) ~ H3 + prolif + age + sex + grade + codel + strata(IDH)
  for (f in names(forms)) {
    vars <- setdiff(all.vars(forms[[f]]), c("OS_time", "OS_event"))
    dd <- d[complete.cases(d[, vars]), ]
    if (!nrow(dd) || !sum(dd$OS_event)) next
    fit <- tryCatch(coxph(forms[[f]], data = dd, singular.ok = FALSE),
                    error = function(e) { note(nm, " ", f, ": ", conditionMessage(e)); NULL })
    if (is.null(fit)) next
    if (f %in% c("plus_prolif", "plus_MKI67", "plus_prolif_grade_codel")) {
      matched <- coxph(forms$base, data = dd)
      b <- cox_row(matched, "H3")
      b$cohort <- nm; b$model <- paste0("matched_base_for_", f); b$term <- "H3"
      a1[[paste(nm, b$model)]] <- b
    }
    if (f != "prolif_only") {
      tv_fun <- function(x, t, ...) x * (log1p(t / 365.25) - log(2))
      tv_forms <- list(update(forms[[f]], . ~ . + tt(H3) + tt(age)))
      tv_names <- f
      if (f != "base") {
        tv_forms <- c(tv_forms, list(update(forms$base, . ~ . + tt(H3) + tt(age))))
        tv_names <- c(tv_names, paste0("matched_base_for_", f))
      }
      for (j in seq_along(tv_forms)) {
        tv <- tryCatch(coxph(tv_forms[[j]], data = dd, tt = tv_fun, singular.ok = FALSE),
                       error = function(e) { note(nm, " time-varying ", tv_names[j], ": ", conditionMessage(e)); NULL })
        if (!is.null(tv)) {
          tr <- cox_row(tv, "H3")
          tr$cohort <- nm; tr$model <- paste0(tv_names[j], "_time_varying_1y"); tr$term <- "H3"
          a1[[paste(nm, tr$model)]] <- tr
        }
      }
    }
    term <- if (f == "prolif_only") "prolif" else "H3"
    r <- cox_row(fit, term); r$cohort <- nm; r$model <- f; r$term <- term
    if (f == "plus_prolif") { rp <- cox_row(fit, "prolif"); r$prolif_HR <- rp$HR; r$prolif_p <- rp$p }
    a1[[paste(nm, f)]] <- r
  }
  if (sum(complete.cases(d[, c("H3", "prolif")])) > 2L) {
    ct <- cor.test(d$H3, d$prolif, method = "spearman", exact = FALSE)
    note(sprintf("A1 %s HEATR3 vs proliferation: rho = %.2f, p = %s",
                 nm, ct$estimate, fmt_p(ct$p.value)))
  }
}
a1 <- rbindlist(a1, fill = TRUE)
a1[, HR_fmt := fmt_hr(HR, lo, hi)]
fwrite(a1, file.path(OUT, "A1_HEATR3_proliferation_adjustment.csv"))
for (i in seq_len(nrow(a1)))
  note(sprintf("A1 %s %-24s %s HR/SD %s p=%s (n=%d, deaths=%d)", a1$cohort[i], a1$model[i], a1$term[i],
               a1$HR_fmt[i], fmt_p(a1$p[i]), a1$n[i], a1$events[i]))
for (m in c("base", "matched_base_for_plus_prolif", "plus_prolif")) {
  x <- a1[model == m & term == "H3" & is.finite(logHR) & is.finite(se) & se > 0]
  if (nrow(x) == 2L) {
    mm <- meta_fe(x$logHR, x$se)
    note(sprintf("A1 exploratory constant-HR meta %s: HR/SD %s p=%s", m,
                 fmt_hr(mm["HR"], mm["lo"], mm["hi"]), fmt_p(mm["p"])))
  }
}
note("A1: Compare matched rows on identical cases. Residual association does not rule out proliferation.")

###############################################################################
## A2. GENOME-WIDE NULL — how special is HEATR3?
###############################################################################
gw_all <- list()
for (nm in names(DAT)) {
  d  <- DAT[[nm]]; ok <- !is.na(d$IDH); d <- d[ok, ]
  ex <- EXPR[[nm]][, d$id, drop = FALSE]
  keep <- usable_rows(ex)
  ex <- ex[keep, , drop = FALSE]
  if (!nrow(ex)) next
  message("A2 ", nm, ": fitting ", nrow(ex), " Cox models (several minutes)...")
  S <- Surv(d$OS_time, d$OS_event)
  res <- t(vapply(seq_len(nrow(ex)), function(i) {
    z <- zs(ex[i, ])
    f <- tryCatch(coxph(S ~ z + d$age + d$sex + strata(d$IDH)), error = function(e) NULL)
    if (is.null(f)) return(c(NA_real_, NA_real_, NA_real_))
    co <- summary(f)$coefficients
    c(exp(co[1, "coef"]), co[1, "z"], co[1, "Pr(>|z|)"])
  }, numeric(3)))
  g <- data.table(gene = rownames(ex), HR = res[, 1], z = res[, 2], p = res[, 3])
  g[, FDR := p.adjust(p, "BH")]
  g[, rank := frank(-abs(z), ties.method = "min", na.last = "keep")]
  g[, cohort := nm]
  gw_all[[nm]] <- g
  h <- g[gene == GENE]
  note(sprintf("A2 %s: %d genes; FDR<0.05: %d (%.1f%%). HEATR3 z=%.2f, rank %d/%d (top %.1f%%), FDR=%s",
               nm, nrow(g), sum(g$FDR < 0.05, na.rm = TRUE), 100 * mean(g$FDR < 0.05, na.rm = TRUE),
               h$z, h$rank, sum(!is.na(g$z)), 100 * h$rank / sum(!is.na(g$z)), fmt_p(h$FDR)))
  if (!nrow(h) || !is.finite(h$z)) { note("A2 HEATR3 not ranked in ", nm); next }
  p <- ggplot(g, aes(z)) + geom_histogram(bins = 100, fill = "grey75") +
    geom_vline(xintercept = if (nrow(h)) h$z else NA_real_, colour = "#C0392B", linewidth = 0.9) +
    annotate("text", x = h$z, y = Inf, vjust = 1.8, hjust = -0.1, label = GENE, colour = "#C0392B") +
    labs(title = paste0(nm, ": HEATR3 against all expressed genes"),
         x = "Cox z per SD (age + sex, stratified by IDH)", y = "Genes") + theme_classic()
  ggsave(file.path(OUT, paste0("A2_genomewide_null_", nm, ".pdf")), p, width = 7, height = 4)
}
gw <- rbindlist(gw_all)
fwrite(gw, file.path(OUT, "A2_genomewide_cox_all_genes.csv"))

if (length(gw_all) >= 2) {
  w <- dcast(gw, gene ~ cohort, value.var = c("z", "p"))
  zc <- grep("^z_", names(w), value = TRUE); pc <- grep("^p_", names(w), value = TRUE)
  w <- w[complete.cases(w[, c(zc, pc), with = FALSE])]
  # same direction and nominal p<0.05 in BOTH cohorts
  w[, replicated := (sign(get(zc[1])) == sign(get(zc[2]))) & get(pc[1]) < 0.05 & get(pc[2]) < 0.05]
  # Stouffer combined z (same sign required to stay large)
  w[, z_comb := (get(zc[1]) + get(zc[2])) / sqrt(2)]
  w[, rank_comb := frank(-abs(z_comb), ties.method = "min")]
  fwrite(w[order(rank_comb)], file.path(OUT, "A2_two_cohort_gene_ranking.csv"))
  h <- w[gene == GENE]
  note(sprintf("A2 BOTH cohorts: %d genes tested in both; %d (%.1f%%) replicate (same sign, p<0.05 in both). ",
               nrow(w), sum(w$replicated), 100 * mean(w$replicated)),
       sprintf("HEATR3 replicated=%s, combined z=%.2f, combined rank %d/%d (top %.1f%%)",
               h$replicated, h$z_comb, h$rank_comb, nrow(w), 100 * h$rank_comb / nrow(w)))
  r <- if (nrow(w) > 1L) cor(w[[zc[1]]], w[[zc[2]]], method = "spearman") else NA_real_
  note(sprintf("A2 Genome-wide agreement of z between cohorts: Spearman rho = %.2f", r))
  note("A2 READ: a high proportion of prognostic genes in LGG is expected. HEATR3 is 'special' if it ranks ",
       "in the top few % within each cohort AND in the two-cohort ranking; otherwise frame it as one of many.")
}

###############################################################################
## A3. RISK-ALLELE HARMONISATION AT 16q12.1
###############################################################################
audit_paths <- if (nzchar(S2_RESULTS)) file.path(S2_RESULTS, "S2_GTEx_curated_top_eQTL.tsv") else
  list.files(file.path(ROOT, "results"), "^S2_GTEx_curated_top_eQTL[.]tsv$",
             recursive = TRUE, full.names = TRUE)
audit_paths <- audit_paths[file.exists(audit_paths)]
if (length(audit_paths)) {
  audit_path <- audit_paths[which.max(file.info(audit_paths)$mtime)]
  audit <- fread(audit_path)[gene == GENE]
  fwrite(audit, file.path(OUT, "A3_HEATR3_GTEx_identity_audit.csv"))
  verified <- file.path(dirname(audit_path), "S2_eQTL_harmonized.tsv")
  if (file.exists(verified)) {
    h <- fread(verified)[gene == GENE & !is.na(risk_slope)]
    fwrite(h, file.path(OUT, "A3_risk_allele_harmonisation.csv"))
    note("A3 HEATR3 aligned records: ", nrow(h))
  } else note("A3 allele direction pending: supply verified GTEx v8 variant-to-rsID mapping in the original S2 pipeline.")
} else note("A3 pending: no S2 identity audit located.")

###############################################################################
## B7. ADDED PROGNOSTIC VALUE
###############################################################################
b7 <- list()
for (nm in names(DAT)) {
  d <- DAT[[nm]]
  covs <- c("age", "sex")
  if ("grade" %in% names(d) && nlevels(d$grade) > 1L && sum(!is.na(d$grade)) > 0.8 * nrow(d)) covs <- c(covs, "grade")
  if ("codel" %in% names(d) && nlevels(d$codel) > 1L && sum(!is.na(d$codel)) > 0.8 * nrow(d)) covs <- c(covs, "codel")
  dd <- d[complete.cases(d[, c("OS_time", "OS_event", "H3", "IDH", covs)]), ]
  if (!nrow(dd) || !sum(dd$OS_event)) next
  f0 <- coxph(as.formula(paste("Surv(OS_time, OS_event) ~", paste(c(covs, "strata(IDH)"), collapse = " + "))), data = dd)
  f1 <- update(f0, . ~ . + H3)
  lr <- anova(f0, f1); lr_p <- lr$`Pr(>|Chi|)`[2]
  cc <- concordance(f0, f1)                       # paired C-index with covariance
  dC <- unname(diff(cc$concordance)); seC <- sqrt(sum(c(1, -1) %*% cc$var %*% c(1, -1)))
  r <- data.table(cohort = nm, n = nrow(dd), events = sum(dd$OS_event),
                  base_covariates = paste(c(covs, "strata(IDH)"), collapse = "+"),
                  C_base = cc$concordance[1], C_plusHEATR3 = cc$concordance[2],
                  delta_C = dC, delta_C_lo = dC - 1.96 * seC, delta_C_hi = dC + 1.96 * seC,
                  delta_C_p = 2 * pnorm(-abs(dC / seC)), LR_p = lr_p,
                  AIC_base = AIC(f0), AIC_plusHEATR3 = AIC(f1))
  # Nested time-varying test adds both HEATR3 terms after allowing age to vary.
  tv_fun <- function(x, t, ...) x * (log1p(t / 365.25) - log(2))
  r$LR_time_varying_p <- tryCatch({
    tv0 <- coxph(update(formula(f0), . ~ . + tt(age)), data = dd, tt = tv_fun)
    tv1 <- coxph(update(formula(f0), . ~ . + tt(age) + H3 + tt(H3)),
                 data = dd, tt = tv_fun)
    ll0 <- logLik(tv0)
    ll1 <- logLik(tv1)
    pchisq(2 * (as.numeric(ll1) - as.numeric(ll0)),
           attr(ll1, "df") - attr(ll0, "df"), lower.tail = FALSE)
  }, error = function(e) { note("B7 time-varying ", nm, ": ", conditionMessage(e)); NA_real_ })
  # Stratified baseline hazards differ. Compare LP-based AUC within IDH strata.
  if (requireNamespace("timeROC", quietly = TRUE)) {
    for (status in levels(dd$IDH)) {
      ix <- which(dd$IDH == status)
      if (length(ix) < 5L || !any(dd$OS_event[ix] == 1)) next
      tryCatch({
        roc <- lapply(list(f0, f1), function(fit)
          timeROC::timeROC(dd$OS_time[ix], dd$OS_event[ix],
                           predict(fit, type = "lp")[ix], cause = 1,
                           times = TIMES_DAYS, iid = FALSE))
        fwrite(data.table(cohort = nm, IDH = status, days = TIMES_DAYS,
                           AUC_base = roc[[1]]$AUC[match(TIMES_DAYS, roc[[1]]$times)],
                           AUC_plus = roc[[2]]$AUC[match(TIMES_DAYS, roc[[2]]$times)]),
               file.path(OUT, paste0("B7_in_sample_AUC_", nm, "_", status, ".csv")))
      }, error = function(e) note("B7 AUC ", nm, " ", status, ": ", conditionMessage(e)))
    }
  }
  b7[[nm]] <- r
  note(sprintf("B7 %s base(%s): C %.3f -> %.3f with HEATR3 (delta %.3f, 95%% CI %.3f to %.3f, p=%s); LR test p=%s; AIC %.1f -> %.1f",
               nm, r$base_covariates, r$C_base, r$C_plusHEATR3, dC, r$delta_C_lo, r$delta_C_hi,
               fmt_p(r$delta_C_p), fmt_p(lr_p), r$AIC_base, r$AIC_plusHEATR3))
}
b7 <- rbindlist(b7, fill = TRUE)
fwrite(b7, file.path(OUT, "B7_added_prognostic_value.csv"))
note("B7: apparent (in-sample) discrimination and exploratory constant-HR models.")

###############################################################################
## C9. PATHWAYS CO-VARYING WITH HEATR3 (IDH-adjusted)
###############################################################################
if (requireNamespace("fgsea", quietly = TRUE) && requireNamespace("msigdbr", quietly = TRUE)) {
  get_sets <- function(coll, sub = NULL) {
    m <- tryCatch(msigdbr::msigdbr(species = "Homo sapiens", collection = coll, subcollection = sub),
                  error = function(e) msigdbr::msigdbr(species = "Homo sapiens", category = coll, subcategory = sub))
    split(m$gene_symbol, m$gs_name)
  }
  sets <- c(get_sets("H"), get_sets("C2", "CP:REACTOME"))
  c9 <- list(); stats_all <- list()
  for (nm in names(DAT)) {
    d <- DAT[[nm]][!is.na(DAT[[nm]]$IDH), ]
    ex <- EXPR[[nm]][, d$id, drop = FALSE]
    ex <- ex[usable_rows(ex), , drop = FALSE]
    if (nrow(ex) < 15L) { note("C9 insufficient eligible genes: ", nm); next }
    # partial correlation with HEATR3, removing IDH (residuals of each on IDH)
    rh  <- resid(lm(d$H3 ~ d$IDH))
    rex <- t(apply(ex, 1, function(x) resid(lm(x ~ d$IDH))))
    rho <- apply(rex, 1, function(x) cor(x, rh, method = "spearman"))
    rho <- rho[is.finite(rho) & names(rho) != GENE]
    if (length(rho) < 15L) { note("C9 insufficient finite correlations: ", nm); next }
    rho <- rho + rnorm(length(rho), 0, 1e-8)        # break ties
    stats_all[[nm]] <- data.table(gene = names(rho), rho = rho, cohort = nm)
    fg <- fgsea::fgsea(pathways = sets, stats = sort(rho, decreasing = TRUE), minSize = 15, maxSize = 500)
    fg <- as.data.table(fg)[order(padj)]
    fg[, leadingEdge := vapply(leadingEdge, function(x) paste(head(x, 15), collapse = ","), "")]
    fg[, cohort := nm]
    c9[[nm]] <- fg
    top <- fg[padj < 0.05][order(-NES)]
    note(sprintf("C9 %s: %d pathways FDR<0.05. Top positive: %s", nm, nrow(top),
                 paste(head(top[NES > 0]$pathway, 5), collapse = "; ")))
    note(sprintf("C9 %s: top negative: %s", nm, paste(head(top[NES < 0][order(NES)]$pathway, 5), collapse = "; ")))
    focus <- fg[grepl("RIBOSOM|RRNA|MYC|P53|E2F|G2M|NUCLEOL|TRANSLATION", pathway)][order(padj)][seq_len(min(.N, 8L))]
    for (i in seq_len(nrow(focus)))
      note(sprintf("C9 %s focus %-55s NES %.2f FDR %s", nm, focus$pathway[i], focus$NES[i], fmt_p(focus$padj[i])))
  }
  if (length(c9)) {
  c9 <- rbindlist(c9)
  fwrite(c9, file.path(OUT, "C9_GSEA_HEATR3_IDH_adjusted.csv"))
  fwrite(rbindlist(stats_all), file.path(OUT, "C9_gene_partial_correlations_with_HEATR3.csv"))
  # pathways significant in BOTH cohorts, same direction
  if (length(unique(c9$cohort)) >= 2) {
    w <- dcast(c9, pathway ~ cohort, value.var = c("NES", "padj"))
    nc <- grep("^NES_", names(w), value = TRUE); pc <- grep("^padj_", names(w), value = TRUE)
    both <- w[get(pc[1]) < 0.05 & get(pc[2]) < 0.05 & sign(get(nc[1])) == sign(get(nc[2]))]
    fwrite(both[order(-abs(get(nc[1])))], file.path(OUT, "C9_GSEA_replicated_both_cohorts.csv"))
    note(sprintf("C9 pathways FDR<0.05 in both cohorts, same direction: %d", nrow(both)))
    # dot plot of top replicated pathways
    if (nrow(both)) {
      topb <- head(both[order(-abs(get(nc[1])))]$pathway, 20)
      pd <- c9[pathway %in% topb]
      pd[, pathway := factor(pathway, levels = rev(topb))]
      p <- ggplot(pd, aes(NES, pathway, colour = cohort, size = -log10(padj))) + geom_point() +
        geom_vline(xintercept = 0, colour = "grey60") + theme_classic() +
        labs(x = "Normalised enrichment score (genes ranked by IDH-adjusted correlation with HEATR3)",
             y = NULL, size = "-log10 FDR", title = "Pathways co-varying with HEATR3 in both cohorts")
      ggsave(file.path(OUT, "C9_GSEA_replicated_dotplot.pdf"), p, width = 10, height = 7)
    }
  }
  }
} else note("C9 skipped: install fgsea (Bioconductor) and msigdbr")

###############################################################################
## SUMMARY
###############################################################################
writeLines(c(paste("Run:", Sys.time()), paste("R", getRversion(), "| survival", packageVersion("survival")),
             "", summary_lines), file.path(OUT, "SUMMARY.txt"))
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
message("\nDone. Send results_strengthen/SUMMARY.txt (and the CSVs if convenient).")