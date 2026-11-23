###############################################################################
#  FIX: Table 1 + Fig. 3 — GBM expression on raw scale (HR 1.00 (1.00–1.00))
#
#  Reads the original merged TCGA cohorts and refits Table 1 and Figure 3.
#  HR per SD uses explicitly configured preprocessing, with original values
#  by default for LGG and log2(x + 1) for GBM. The observed range is logged for auditing, not used to infer scale.
#  Needs: survival, data.table, ggplot2; optional R.utils for Xena gzip fallback.
###############################################################################

suppressPackageStartupMessages({ library(survival); library(data.table); library(ggplot2) })

## ================================ CONFIG ===================================
BASE <- "D:/RESIDENSI/Thesis-PPDS/R"
OUT  <- file.path(BASE, "results_final"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

GENES <- c("HEATR3", "LMF1", "EGFR", "RTEL1", "SLC16A8", "MIR4300HG")

# Each cohort: a clinical CSV (may already contain the gene columns) and,
# optionally, a genes x samples expression file (first column = gene symbol).
CORE_COLUMNS <- c(id = "patient_id", OS_time = "OS_time", OS_event = "OS_event",
                  age = "age_at_diagnosis", sex = "gender")
COHORTS <- list(
  LGG = list(clin = file.path(BASE, "TCGA_LGG_merged_selected_genes_with_IDH1.csv"),
             expr = file.path(BASE, "LGG_expression_xena.tsv.gz"),
             colmap = CORE_COLUMNS, age_unit = "years", transform = "as_provided"),
  GBM = list(clin = file.path(BASE, "TCGA_GBM_merged.csv"),
             expr = NA_character_,
             colmap = CORE_COLUMNS, age_unit = "days", transform = "log2p1")
)
# Set transform = "log2p1" only for verified linear nonnegative expression,
# or "already_log2" for verified log2 inputs. Affects merged expression only.
# The optional Xena fallback keeps its source scale.
EXPECTED <- list(LGG = c(n = 512, deaths = 126), GBM = c(n = 257, deaths = 217))
DPI <- 600

## ================================ HELPERS ==================================
LOG <- character(); note <- function(...) { s <- paste0(...); message(s); LOG <<- c(LOG, s) }
fhr <- function(h, l, u) sprintf("%.2f (%.2f\u2013%.2f)", h, l, u)
fp  <- function(p) ifelse(is.na(p), "\u2014", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))

load_cohort <- function(cfg) {
  d <- as.data.frame(fread(cfg$clin, check.names = FALSE))
  miss <- setdiff(cfg$colmap, names(d))
  if (length(miss)) stop(basename(cfg$clin), ": columns not found: ", paste(miss, collapse = ", "),
                         "\n  available: ", paste(head(names(d), 60), collapse = ", "), "\n  -> edit colmap in CONFIG")
  cl <- setNames(d[, unname(cfg$colmap), drop = FALSE], names(cfg$colmap))
  genes_in_csv <- intersect(GENES, names(d))
  for (g in genes_in_csv) cl[[g]] <- suppressWarnings(as.numeric(d[[g]]))
  # genes not in the CSV -> take them from the expression matrix
  need <- setdiff(GENES, genes_in_csv)
  for (g in genes_in_csv) {
    x <- cl[[g]]
    if (cfg$transform == "log2p1") {
      if (any(x < 0, na.rm = TRUE)) stop("Negative expression in ", g)
      cl[[g]] <- log2(x + 1)
    }
  }
  canonical <- function(id) substr(gsub("[.]", "-", toupper(trimws(id))), 1L, 12L)
  cl$id <- canonical(cl$id)
  if (anyNA(cl$id) || any(!nzchar(cl$id)) || anyDuplicated(cl$id))
    stop("Invalid or duplicate clinical patient IDs in ", cfg$clin)
  if (length(need) && !is.na(cfg$expr) && file.exists(cfg$expr)) {
    tryCatch({
      e <- fread(cfg$expr, data.table = FALSE, check.names = FALSE)
      e <- e[toupper(trimws(e[[1L]])) %in% need, , drop = FALSE]
      if (nrow(e)) {
        symbols <- toupper(trimws(e[[1L]]))
        barcodes <- gsub("[.]", "-", toupper(trimws(names(e)[-1L])))
        keep <- grepl("^TCGA-.{2}-.{4}(-01.*)?$", barcodes)
        ids <- canonical(barcodes[keep])
        if (anyDuplicated(ids)) stop("Multiple primary expression samples per patient.")
        if (!any(cl$id %in% ids)) stop("No matching TCGA patient barcodes.")
        for (g in names(which(table(symbols) == 1L))) {
          values <- as.numeric(unlist(e[symbols == g, -1L, drop = FALSE], use.names = FALSE))
          cl[[g]] <- values[keep][match(cl$id, ids)]
          note("Expression fallback for ", g, ": ", basename(cfg$expr), " (as provided)")
        }
      } else note("None of the missing genes was found by symbol in ", basename(cfg$expr))
    }, error = function(e) note("Optional Xena input skipped: ", conditionMessage(e)))
  }
  cl$OS_time <- as.numeric(cl$OS_time)
  cl$OS_event <- as.numeric(cl$OS_event)
  cl$age <- as.numeric(cl$age)
  if (cfg$age_unit == "days") cl$age <- cl$age / 365.25
  cl$sex <- factor(tolower(trimws(cl$sex)), levels = c("female", "male"))
  if (any(!is.na(cl$OS_event) & !cl$OS_event %in% 0:1) ||
      any(cl$OS_time < 0, na.rm = TRUE) ||
      any(cl$age < 0 | cl$age > 120, na.rm = TRUE))
    stop("Check survival coding or age units in ", cfg$clin)
  cl <- cl[is.finite(cl$OS_time) & !is.na(cl$OS_event), , drop = FALSE]
  attr(cl, "verified_log2_genes") <- if (cfg$transform %in% c("log2p1", "already_log2")) genes_in_csv else character()
  cl
}

cox_gene <- function(d, g, adjust) {
  dd <- d[is.finite(d[[g]]), ]
  if (adjust) dd <- dd[!is.na(dd$age) & !is.na(dd$sex), ]
  reference <- d[[g]][is.finite(d[[g]])]
  sigma <- sd(reference)
  if (nrow(dd) < 2L || !sum(dd$OS_event) || !is.finite(sigma) || sigma <= 0)
    stop("Insufficient events or expression variation.")
  dd$x_sd <- (dd[[g]] - mean(reference)) / sigma
  rhs <- if (adjust) " + age + sex" else ""
  f_sd  <- coxph(as.formula(paste0("Surv(OS_time, OS_event) ~ x_sd", rhs)),  data = dd)
  s1 <- summary(f_sd)
  log_ci <- if (g %in% attr(d, "verified_log2_genes"))
    exp(s1$coefficients[1, "coef"] / sigma +
          c(0, -1, 1) * qnorm(.975) * s1$coefficients[1, "se(coef)"] / sigma) else rep(NA_real_, 3L)
  z <- tryCatch(cox.zph(f_sd)$table, error = function(e) NULL)
  data.table(HR_sd = s1$conf.int[1, 1], lo_sd = s1$conf.int[1, 3], hi_sd = s1$conf.int[1, 4],
             p = s1$coefficients[1, "Pr(>|z|)"],
             HR_log2 = log_ci[1], lo_log2 = log_ci[2], hi_log2 = log_ci[3], SD_reference = sigma,
             n = s1$n, deaths = s1$nevent,
             PH_gene_p = if (!is.null(z)) z[1, "p"] else NA_real_,
             PH_global_p = if (!is.null(z)) z["GLOBAL", "p"] else NA_real_)
}

## ================================ LOAD + SCALE CHECK =======================
DAT <- list()
for (coh in names(COHORTS)) {
  d <- load_cohort(COHORTS[[coh]])
  for (g in intersect(GENES, names(d))) {
    x <- d[[g]]
    x <- x[is.finite(x)]
    if (length(x)) note(sprintf("SCALE %s %s n=%d min=%g median=%g max=%g; merged transform=%s",
                                coh, g, length(x), min(x), median(x), max(x), COHORTS[[coh]]$transform))
  }
  DAT[[coh]] <- d
}

## ================================ TABLE 1 ==================================
rows <- list()
for (coh in names(DAT)) for (g in GENES) {
  d <- DAT[[coh]]
  if (!g %in% names(d) || sum(is.finite(d[[g]])) < 2) {
    rows[[paste(coh, g)]] <- data.table(cohort = coh, gene = g, note = "not available"); next }
  safe_fit <- function(adjust) tryCatch(cox_gene(d, g, adjust), error = function(e) {
    note(coh, " ", g, " adjusted=", adjust, ": ", conditionMessage(e))
    data.table(HR_sd = NA_real_, lo_sd = NA_real_, hi_sd = NA_real_, p = NA_real_,
               HR_log2 = NA_real_, lo_log2 = NA_real_, hi_log2 = NA_real_,
               SD_reference = NA_real_, n = NA_integer_, deaths = NA_integer_,
               PH_gene_p = NA_real_, PH_global_p = NA_real_)
  })
  u <- safe_fit(FALSE)
  a <- safe_fit(TRUE)
  rows[[paste(coh, g)]] <- data.table(cohort = coh, gene = g,
    uni_HR_sd = u$HR_sd, uni_lo_sd = u$lo_sd, uni_hi_sd = u$hi_sd, uni_p = u$p,
    adj_HR_sd = a$HR_sd, adj_lo_sd = a$lo_sd, adj_hi_sd = a$hi_sd, adj_p = a$p,
    adj_HR_log2 = a$HR_log2, adj_lo_log2 = a$lo_log2, adj_hi_log2 = a$hi_log2,
    SD_reference = a$SD_reference, uni_n = u$n, uni_deaths = u$deaths,
    uni_PH_gene_p = u$PH_gene_p, uni_PH_global_p = u$PH_global_p,
    uni_HR_log2 = u$HR_log2, uni_lo_log2 = u$lo_log2, uni_hi_log2 = u$hi_log2,
    n = a$n, deaths = a$deaths, PH_gene_p = a$PH_gene_p, PH_global_p = a$PH_global_p)
}
t1 <- rbindlist(rows, fill = TRUE)
if (!"adj_p" %in% names(t1)) stop("No candidate gene expression available in either cohort.")
t1[, adj_FDR := p.adjust(adj_p, "BH", n = length(GENES)), by = cohort]
t1[, uni_FDR := p.adjust(uni_p, "BH", n = length(GENES)), by = cohort]
t1[, `:=`(uni_SD = ifelse(is.na(uni_HR_sd), "Not available", fhr(uni_HR_sd, uni_lo_sd, uni_hi_sd)),
          adj_SD = ifelse(is.na(adj_HR_sd), "Not available", fhr(adj_HR_sd, adj_lo_sd, adj_hi_sd)),
          adj_log2 = ifelse(is.na(adj_HR_log2), "Not available", fhr(adj_HR_log2, adj_lo_log2, adj_hi_log2)))]
fwrite(t1, file.path(OUT, "Table1_per_SD_full.csv"))

for (coh in names(DAT)) {
  x <- t1[cohort == coh & !is.na(n)]
  if (nrow(x)) note(sprintf("N %s: n=%d, deaths=%d (expected %d / %d)", coh, max(x$n), max(x$deaths),
                            EXPECTED[[coh]]["n"], EXPECTED[[coh]]["deaths"]))
}
for (i in seq_len(nrow(t1))) note(sprintf("T1 %s %-9s uni %s | adj/SD %s p=%s FDR=%s | adj/log2 %s",
  t1$cohort[i], t1$gene[i], t1$uni_SD[i], t1$adj_SD[i], fp(t1$adj_p[i]), fp(t1$adj_FDR[i]), t1$adj_log2[i]))

# Manuscript-ready wide table (one row per gene)
wide <- dcast(t1, gene ~ cohort, value.var = c("uni_SD", "adj_SD", "adj_p", "adj_FDR"))
wide <- wide[match(GENES, gene)]
for (cc in grep("_p_|_FDR_", names(wide), value = TRUE)) set(wide, j = cc, value = fp(wide[[cc]]))
setcolorder(wide, c("gene", grep("LGG$", names(wide), value = TRUE), grep("GBM$", names(wide), value = TRUE)))
fwrite(wide, file.path(OUT, "Table1_per_SD_manuscript.csv"))

## ================================ RATIO OF HRs =============================
if (all(c("LGG", "GBM") %in% names(DAT)) && all(sapply(DAT, function(d) "HEATR3" %in% names(d)))) {
  P <- rbindlist(lapply(names(DAT), function(coh) {
    d <- DAT[[coh]]
    d$H3 <- as.numeric(scale(d$HEATR3))
    d <- d[is.finite(d$H3) & !is.na(d$age) & !is.na(d$sex), ]
    data.table(OS_time = d$OS_time, OS_event = d$OS_event, age = d$age, sex = d$sex,
               H3 = d$H3, cohort = coh)       # SD within each cohort
  }))
  P[, cohort := factor(cohort, levels = c("LGG", "GBM"))]
  f <- coxph(Surv(OS_time, OS_event) ~ H3 + H3:cohort + age + sex + strata(cohort), data = P)
  s <- summary(f); i <- match("H3:cohortGBM", rownames(s$coefficients))
  fwrite(data.table(ratio_HR = s$conf.int[i, 1], low = s$conf.int[i, 3],
                    high = s$conf.int[i, 4], p = s$coefficients[i, "Pr(>|z|)"],
                    n = s$n, deaths = s$nevent),
         file.path(OUT, "Table1_per_SD_HEATR3_interaction.csv"))
  note(sprintf("RATIO HEATR3 HR per SD, GBM vs LGG (age + sex, stratified by cohort): %s, p = %s (n = %d, deaths = %d)",
               fhr(s$conf.int[i, 1], s$conf.int[i, 3], s$conf.int[i, 4]),
               formatC(s$coefficients[i, "Pr(>|z|)"], format = "g", digits = 3), s$n, s$nevent))
}

## ================================ FIG. 3 ===================================
fd <- t1[is.finite(adj_HR_sd) & is.finite(adj_lo_sd) & is.finite(adj_hi_sd) & adj_lo_sd > 0]
if (!nrow(fd)) stop("No estimable adjusted HR available for Figure 3.")
lab_n <- fd[, .(n = max(n)), by = cohort]
fd[, panel := factor(paste0("TCGA-", cohort, " (n = ", lab_n$n[match(cohort, lab_n$cohort)], ")"),
                     levels = paste0("TCGA-", c("LGG", "GBM"), " (n = ", lab_n$n[match(c("LGG", "GBM"), lab_n$cohort)], ")"))]
fd[, gene := factor(gene, levels = rev(GENES))]
fd[, sig := !is.na(adj_FDR) & adj_FDR < 0.05]
fd[, lab := sprintf("%.2f", adj_HR_sd)]
# keep MIR4300HG row visible in LGG as "n/a"
na_rows <- t1[is.na(adj_HR_sd), .(cohort, gene)]
if (nrow(na_rows)) {
  na_rows[, panel := factor(paste0("TCGA-", cohort, " (n = ", lab_n$n[match(cohort, lab_n$cohort)], ")"), levels = levels(fd$panel))]
  na_rows[, gene := factor(gene, levels = rev(GENES))]
}
xr <- range(c(fd$adj_lo_sd, fd$adj_hi_sd), na.rm = TRUE)
p <- ggplot(fd, aes(adj_HR_sd, gene)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey55") +
  geom_errorbarh(aes(xmin = adj_lo_sd, xmax = adj_hi_sd), height = 0.22, colour = "grey25") +
  geom_point(aes(colour = sig), size = 2.8, show.legend = FALSE) +
  geom_text(aes(label = lab), vjust = -0.9, size = 3, colour = "grey20") +
  { if (nrow(na_rows)) geom_text(data = na_rows, aes(x = 1, y = gene, label = "not available"),
                                 inherit.aes = FALSE, size = 3, colour = "grey50", fontface = "italic") } +
  scale_colour_manual(values = c(`TRUE` = "#B03A2E", `FALSE` = "grey40")) +
  scale_x_log10(breaks = c(0.6, 0.8, 1, 1.25, 1.5, 2), limits = c(min(0.6, xr[1]), max(2, xr[2]))) +
  facet_wrap(~panel, nrow = 1) +
  labs(x = "Hazard ratio per 1 SD higher expression (age- and sex-adjusted, log scale)", y = NULL,
       title = "Risk-gene expression and overall survival: LGG vs GBM",
       subtitle = "Red: Benjamini\u2013Hochberg FDR < 0.05 within cohort") +
  theme_classic(base_size = 11) +
  theme(axis.text.y = element_text(face = "italic"), strip.background = element_blank(),
        strip.text = element_text(face = "bold"), plot.title = element_text(face = "bold"),
        panel.spacing = unit(1.2, "lines"))
ggsave(file.path(OUT, "Fig3_forest_LGG_GBM_perSD.pdf"), p, width = 8.5, height = 3.8, device = if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf)
ggsave(file.path(OUT, "Fig3_forest_LGG_GBM_perSD.tiff"), p, width = 8.5, height = 3.8, dpi = DPI, compression = "lzw")
note("FIG3 written: Fig3_forest_LGG_GBM_perSD.pdf / .tiff (", DPI, " dpi)")

writeLines(c(paste("Run:", Sys.time()), "", LOG), file.path(OUT, "SUMMARY_Table1_fix.txt"))
message("\nDone -> ", file.path(OUT, "SUMMARY_Table1_fix.txt"))