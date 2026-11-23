###############################################################################
#  HEATR3 / glioma — FINAL ANALYSES BEFORE SUBMISSION
#
#  PART 1  CGGA-693 rerun (fixed IDs): HEATR3 vs proliferation, added value,
#          time-varying effect; TCGA side-by-side; meta-analysis
#  PART 2  Administrative data
#          2a  Table 1 final: 6 genes x (TCGA-LGG, TCGA-GBM n=257), FDR,
#              HEATR3 ratio of HRs GBM vs LGG
#          2b  GTEx variant identity check via Ensembl REST (rsID, alleles,
#              LD with GWAS sentinel rs10852606) -> keep or drop the
#              "risk allele increases HEATR3" claim
#          2c  Publication figures, 600 dpi TIFF + vector PDF
#  PART 3  GSEA (Hallmark + Reactome) of genes co-varying with HEATR3
#
#  Needs: survival, data.table, ggplot2   (already installed)
#         fgsea (installed), msigdbr OR a downloaded MSigDB .gmt file (Part 3)
#         jsonlite (Part 2b; install.packages("jsonlite") if missing)
#  R.utils is required by fread for compressed GTEx inputs.
#  Internet is used for the Ensembl audit and may be used by msigdbr.
#
#  Run section by section in RStudio, or source() the whole file.
#  Every number for the manuscript is written to results_final/SUMMARY_FINAL.txt
###############################################################################

suppressPackageStartupMessages({ library(survival); library(data.table); library(ggplot2) })

## ================================ CONFIG ===================================
BASE <- "D:/RESIDENSI/Thesis-PPDS/R"           # <- your project folder
OUT <- file.path(BASE, "results_final")


# Use the filenames and column schema confirmed in the project runs.
if (!dir.exists(BASE)) stop("Project folder not found: ", BASE)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
find_input <- function(name, use_latest = FALSE) {
  hits <- list.files(BASE, recursive = TRUE, full.names = TRUE)
  hits <- hits[basename(hits) == name]
  if (length(hits) == 1L) return(hits)
  if (length(hits) && use_latest) return(hits[which.max(file.info(hits)$mtime)])
  if (length(hits)) stop("More than one copy of ", name, ". Set its path explicitly in CONFIG.")
  file.path(BASE, name)
}
TCGA_LGG_CLIN <- file.path(BASE, "TCGA_LGG_merged_selected_genes_with_IDH1.csv")
TCGA_LGG_EXPR <- find_input("LGG_expression_xena.tsv.gz")
TCGA_COLMAP <- c(id = "patient_id", OS_time = "OS_time", OS_event = "OS_event",
                 age = "age_at_diagnosis", sex = "gender", IDH = "IDH1_status",
                 HEATR3_original = "HEATR3")

CGGA_CLIN_RAW <- find_input("CGGA.mRNAseq_693_clinical.20200506.txt")
CGGA_EXPR_RAW <- find_input("CGGA.mRNAseq_693.RSEM-genes.20200506.txt")
CGGA_LOG2P1 <- FALSE # FALSE reproduces the original CGGA HEATR3 scale.

# Reuse S4 to preserve the original GBM dataset and age conversion.
TCGA_S4 <- find_input("S4_Cox_Table1.tsv", use_latest = TRUE)
TCGA_INTERACTION <- file.path(dirname(TCGA_S4), "S4_HEATR3_cohort_interaction.tsv")
# Optional full paths. Leave blank to detect the GTEx files automatically.
GTEX_PAIRS_PATH <- ""
GTEX_EGENES_PATH <- ""
find_gtex <- function(name, explicit = "") {
  if (nzchar(explicit)) {
    if (file.exists(explicit)) return(explicit)
    warning("GTEx input not found: ", explicit, call. = FALSE)
    return(NA_character_)
  }
  hits <- list.files(BASE, recursive = TRUE, full.names = TRUE)
  hits <- sort(hits[basename(hits) == name])
  if (!length(hits)) {
    message("GTEx input not found: ", name)
    return(NA_character_)
  }
  if (length(hits) == 1L) return(hits[[1L]])
  hashes <- unname(tools::md5sum(hits))
  if (all(!is.na(hashes)) && length(unique(hashes)) == 1L) {
    message("Identical GTEx copies found for ", name, ". Using ", hits[[1L]])
    return(hits[[1L]])
  }
  warning("Multiple GTEx copies differ or cannot be read: ", name, "\n",
          paste(hits, collapse = "\n"),
          "\nSet its explicit path in CONFIG to enable raw GTEx analyses.",
          call. = FALSE)
  NA_character_
}
GTEX_PAIRS <- find_gtex("Brain_Cortex.v8.signif_variant_gene_pairs.txt.gz", GTEX_PAIRS_PATH)
GTEX_EGENES <- find_gtex("Brain_Cortex.v8.egenes.txt.gz", GTEX_EGENES_PATH)
have_raw_gtex <- all(!is.na(c(GTEX_PAIRS, GTEX_EGENES))) &&
  requireNamespace("R.utils", quietly = TRUE)
GTEX_AUDIT <- find_input("A3_HEATR3_GTEx_identity_audit.csv", use_latest = TRUE)
if (!file.exists(GTEX_AUDIT))
  GTEX_AUDIT <- find_input("S2_GTEx_curated_top_eQTL.tsv", use_latest = TRUE)
GTEX_HARMONIZED <- file.path(dirname(TCGA_S4), "S2_eQTL_harmonized.tsv")
SENTINEL <- "rs10852606"
RUN_ENSEMBL <- TRUE

# ---- GSEA gene sets: msigdbr if installed, otherwise .gmt files from
#      https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp (symbols .gmt)
GMT_FILES <- c(file.path(BASE, "h.all.v2024.1.Hs.symbols.gmt"),
               file.path(BASE, "c2.cp.reactome.v2024.1.Hs.symbols.gmt"))

GENES  <- c("HEATR3", "LMF1", "EGFR", "RTEL1", "SLC16A8", "MIR4300HG")
GENE   <- "HEATR3"
PROLIF <- c("MKI67", "TOP2A", "CCNB1", "CDK1", "BUB1", "PCNA", "MCM2", "AURKA", "CCNA2", "TPX2")
YEARS  <- c(1, 3, 5)
DPI    <- 600
set.seed(2026)

## ================================ HELPERS ==================================
LOG <- character()
note <- function(...) { s <- paste0(...); message(s); LOG <<- c(LOG, s) }
fhr  <- function(h, l, u) sprintf("%.2f (%.2f\u2013%.2f)", h, l, u)
fp   <- function(p) ifelse(is.na(p), "NA", ifelse(p < 0.001, formatC(p, format = "e", digits = 2), sprintf("%.3f", p)))
zs   <- function(x) as.numeric(scale(x))

read_expr <- function(f) {
  if (!file.exists(f)) stop("Missing expression input: ", f)
  if (grepl("[.]gz$", f, ignore.case = TRUE)) {
    tmp <- tempfile(fileext = ".tsv")
    on.exit(unlink(tmp), add = TRUE)
    input <- gzfile(f, "rb")
    output <- file(tmp, "wb")
    tryCatch({
      repeat {
        chunk <- readBin(input, "raw", n = 1048576L)
        if (!length(chunk)) break
        writeBin(chunk, output)
      }
    }, finally = { close(input); close(output) })
    f <- tmp
  }
  d <- fread(f, data.table = FALSE, check.names = FALSE)
  g <- toupper(trimws(as.character(d[[1]]))); m <- as.matrix(d[, -1, drop = FALSE]); storage.mode(m) <- "numeric"; rownames(m) <- g
  m[!duplicated(g) & !duplicated(g, fromLast = TRUE) & !is.na(g) & g != "", , drop = FALSE]
}
std_clin <- function(f, colmap) {
  d <- as.data.frame(fread(f, check.names = FALSE))
  miss <- setdiff(colmap, names(d))
  if (length(miss)) stop(basename(f), ": columns not found: ", paste(miss, collapse = ", "),
                         "\n  available: ", paste(names(d), collapse = ", "),
                         "\n  -> edit the COLMAP in CONFIG")
  out <- setNames(d[, colmap, drop = FALSE], names(colmap))
  out$OS_time <- as.numeric(out$OS_time); out$OS_event <- as.numeric(out$OS_event)
  out$age <- as.numeric(out$age); out$sex <- factor(tolower(out$sex))
  if ("IDH" %in% names(out))
    out$IDH <- factor(ifelse(grepl("wild|^wt", out$IDH, ignore.case = TRUE), "wildtype",
                      ifelse(grepl("mut", out$IDH, ignore.case = TRUE), "mutant", NA)),
                      levels = c("mutant", "wildtype"))
  out[!is.na(out$OS_time) & out$OS_time > 0 & !is.na(out$OS_event), ]
}
attach_expr <- function(cl, ex, genes) {
  cl$id <- trimws(as.character(cl$id))
  colnames(ex) <- trimws(colnames(ex))
  if (any(grepl("^TCGA[-.]", cl$id, ignore.case = TRUE))) {
    patient <- function(x) substr(gsub("[.]", "-", toupper(x)), 1L, 12L)
    barcode <- gsub("[.]", "-", toupper(colnames(ex)))
    keep <- nchar(barcode) == 12L | grepl("^TCGA-.{2}-.{4}-01", barcode)
    ex <- ex[, keep, drop = FALSE]
    colnames(ex) <- patient(colnames(ex))
    cl$id <- patient(cl$id)
  }
  if (anyNA(cl$id) || any(!nzchar(cl$id)) || anyDuplicated(cl$id) || anyDuplicated(colnames(ex)))
    stop("Missing or duplicate patient/sample IDs after alignment.")
  ids <- intersect(cl$id, colnames(ex))
  if (!length(ids)) stop("No matching IDs between clinical and expression. Example clinical: ",
                         paste(head(cl$id, 3), collapse = ", "), " | expression: ",
                         paste(head(colnames(ex), 3), collapse = ", "))
  cl <- cl[match(ids, cl$id), ]; ex <- ex[, ids, drop = FALSE]
  for (g in intersect(genes, rownames(ex))) cl[[g]] <- as.numeric(ex[g, ])
  list(cl = cl, ex = ex)
}
add_scores <- function(d, ex) {
  if (!GENE %in% rownames(ex)) stop("HEATR3 expression row not found.")
  if (any(!is.na(d$OS_event) & !d$OS_event %in% 0:1)) stop("Unexpected event coding.")
  d$sex <- factor(trimws(tolower(as.character(d$sex))), levels = c("female", "male"))
  ok <- complete.cases(d[, c("OS_time", "OS_event", "age", "sex")]) & is.finite(ex[GENE, ])
  d$H3 <- rep(NA_real_, nrow(d))
  d$H3[ok] <- zs(ex[GENE, ok])
  pg <- intersect(PROLIF, rownames(ex))
  pg <- pg[vapply(pg, function(g) all(is.finite(ex[g, ])) && sd(ex[g, ]) > 0, logical(1))]
  d$prolif <- if (length(pg) >= 3L) zs(colMeans(t(scale(t(ex[pg, , drop = FALSE]))))) else NA_real_
  d$MKI67z <- if ("MKI67" %in% rownames(ex)) zs(ex["MKI67", ]) else NA
  attr(d, "prolif_genes") <- pg
  d
}
cox1 <- function(fit, term) {
  s <- summary(fit); i <- which(rownames(s$coefficients) == term)
  data.frame(HR = s$conf.int[i, 1], lo = s$conf.int[i, 3], hi = s$conf.int[i, 4],
             p = s$coefficients[i, "Pr(>|z|)"], b = s$coefficients[i, "coef"],
             se = s$coefficients[i, "se(coef)"], n = s$n, events = s$nevent)
}
meta_fe <- function(b, se) {
  w <- 1 / se^2; bm <- sum(w * b) / sum(w); sm <- sqrt(1 / sum(w)); Q <- sum(w * (b - bm)^2); k <- length(b) - 1
  c(HR = exp(bm), lo = exp(bm - 1.96 * sm), hi = exp(bm + 1.96 * sm), p = 2 * pnorm(-abs(bm / sm)),
    Qp = if (k > 0) 1 - pchisq(Q, k) else NA, I2 = if (k > 0 && Q > 0) max(0, (Q - k) / Q) * 100 else 0)
}
# time-varying HEATR3: time basis centered at one year, matching the TCGA and CGGA runs; age also time-varying
tv_fit <- function(d, extra = NULL) {
  rhs <- paste(c("H3", "tt(H3)", "age", "tt(age)", "sex", extra, "strata(IDH)"), collapse = " + ")
  coxph(as.formula(paste("Surv(OS_time, OS_event) ~", rhs)), data = d,
        tt = function(x, t, ...) x * (log1p(t / 365.25) - log(2)))
}
tv_hr <- function(fit, years = YEARS) {
  b <- coef(fit); V <- vcov(fit)
  i1 <- which(names(b) == "H3"); i2 <- which(names(b) == "tt(H3)")
  rbindlist(lapply(years, function(y) {
    L <- numeric(length(b)); L[i1] <- 1; L[i2] <- log1p(y) - log(2)
    est <- sum(L * b); se <- sqrt(drop(t(L) %*% V %*% L))
    data.table(year = y, HR = exp(est), lo = exp(est - 1.96 * se), hi = exp(est + 1.96 * se),
               p = 2 * pnorm(-abs(est / se)))
  }))
}
save_fig <- function(p, name, w, h) {
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h, device = if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf)
  ggsave(file.path(OUT, paste0(name, ".tiff")), p, width = w, height = h, dpi = DPI, compression = "lzw")
}
theme_pub <- theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12), legend.position = "top")

###############################################################################
## PART 1 — CGGA RERUN WITH FIXED IDs
###############################################################################
message("\n===== PART 1: CGGA =====")
cc <- as.data.frame(fread(CGGA_CLIN_RAW, check.names = FALSE))
names(cc) <- make.names(names(cc))
cens_col <- grep("^Censor", names(cc), value = TRUE)[1]
codel_col <- grep("1p19q", names(cc), value = TRUE)[1]
cc <- cc[which(tolower(trimws(cc$PRS_type)) == "primary" & toupper(trimws(cc$Grade)) %in% c("II", "III", "WHO II", "WHO III", "2", "3")), ]
cgga_cl <- data.frame(id = cc$CGGA_ID, OS_time = as.numeric(cc$OS), OS_event = as.numeric(cc[[cens_col]]),
                      age = as.numeric(cc$Age), sex = factor(tolower(cc$Gender)),
                      IDH = factor(ifelse(grepl("wild", cc$IDH_mutation_status, ignore.case = TRUE), "wildtype",
                                   ifelse(grepl("mut", cc$IDH_mutation_status, ignore.case = TRUE), "mutant", NA)),
                                   levels = c("mutant", "wildtype")),
                      grade = factor(cc$Grade),
                      codel = factor(ifelse(grepl("non", cc[[codel_col]], ignore.case = TRUE), "non-codel",
                                     ifelse(grepl("codel", cc[[codel_col]], ignore.case = TRUE), "codel", NA))))
cgga_cl <- cgga_cl[!is.na(cgga_cl$OS_time) & cgga_cl$OS_time > 0 & !is.na(cgga_cl$OS_event), ]
cgga_ex <- read_expr(CGGA_EXPR_RAW)
if (CGGA_LOG2P1) cgga_ex <- log2(cgga_ex + 1)
x <- attach_expr(cgga_cl, cgga_ex, c(GENES, PROLIF)); cgga <- add_scores(x$cl, x$ex); cgga_ex <- x$ex
note(sprintf("P1 CGGA-693 primary grade II-III with expression: n=%d, deaths=%d; IDH mut/wt/NA=%d/%d/%d",
             nrow(cgga), sum(cgga$OS_event), sum(cgga$IDH %in% "mutant"), sum(cgga$IDH %in% "wildtype"), sum(is.na(cgga$IDH))))

tl_cl <- std_clin(TCGA_LGG_CLIN, TCGA_COLMAP)
tl_ex <- read_expr(TCGA_LGG_EXPR)
x <- attach_expr(tl_cl, tl_ex, c(GENES, PROLIF))
x$ex <- rbind(x$ex[rownames(x$ex) != GENE, , drop = FALSE], HEATR3 = x$cl$HEATR3_original)
x$cl$HEATR3 <- x$cl$HEATR3_original
tcga <- add_scores(x$cl, x$ex); tcga_ex <- x$ex
note(sprintf("P1 TCGA-LGG: n=%d, deaths=%d", nrow(tcga), sum(tcga$OS_event)))

COH <- list(TCGA_LGG = tcga, CGGA_693 = cgga)
p1 <- list(); tv_all <- list(); b7 <- list()
for (nm in names(COH)) {
  d <- COH[[nm]]; d <- d[complete.cases(d[, c("IDH", "age", "sex", "H3")]), ]
  if (sum(complete.cases(d[, c("H3", "prolif")])) > 2L) {
  rho <- cor.test(d$H3, d$prolif, method = "spearman")
  note(sprintf("P1 %s: HEATR3 vs proliferation Spearman rho=%.2f p=%s (n=%d)", nm, rho$estimate, fp(rho$p.value), nrow(d)))
  }
  M <- list(base        = Surv(OS_time, OS_event) ~ H3 + age + sex + strata(IDH),
            plus_prolif = Surv(OS_time, OS_event) ~ H3 + prolif + age + sex + strata(IDH),
            plus_MKI67  = Surv(OS_time, OS_event) ~ H3 + MKI67z + age + sex + strata(IDH))
  for (m in names(M)) {
    dd <- d[complete.cases(d[, all.vars(M[[m]]), drop = FALSE]), , drop = FALSE]
    if (!nrow(dd) || !sum(dd$OS_event)) { note("P1 ", nm, " missing inputs: ", m); next }
    if (m != "base") {
      matched <- cox1(coxph(M$base, data = dd), "H3")
      matched$cohort <- nm; matched$model <- paste0("matched_base_for_", m)
      p1[[paste(nm, matched$model)]] <- matched
    }
    r <- cox1(coxph(M[[m]], data = dd), "H3"); r$cohort <- nm; r$model <- m; p1[[paste(nm, m)]] <- r
    note(sprintf("P1 %s %-12s HEATR3 HR/SD %s p=%s (n=%d, deaths=%d)", nm, m, fhr(r$HR, r$lo, r$hi), fp(r$p), r$n, r$events))
  }
  if (any(is.finite(d$prolif))) {
  fp_ <- coxph(M$plus_prolif, data = d); rp <- cox1(fp_, "prolif")
  note(sprintf("P1 %s proliferation score in same model: HR/SD %s p=%s", nm, fhr(rp$HR, rp$lo, rp$hi), fp(rp$p)))
  }
  # PH check for HEATR3 in base model
  z <- cox.zph(coxph(M$base, data = d)); note(sprintf("P1 %s Schoenfeld HEATR3 p=%s, global p=%s", nm,
                                                     fp(z$table["H3", "p"]), fp(z$table["GLOBAL", "p"])))
  # time-varying, without and with proliferation
  for (ex_ in list(NULL, "prolif")) {
    dd <- d[complete.cases(d[, c("OS_time", "OS_event", "H3", "age", "sex", "IDH", ex_), drop = FALSE]), ]
    if (!nrow(dd) || !sum(dd$OS_event)) next
    f <- tv_fit(dd, ex_); lab <- if (is.null(ex_)) "time-varying" else "time-varying + prolif"
    f0 <- coxph(as.formula(paste("Surv(OS_time, OS_event) ~ H3 + age + tt(age) + sex +",
                                 if (is.null(ex_)) "" else "prolif +", "strata(IDH)")),
                data = dd, tt = function(x, t, ...) x * (log1p(t / 365.25) - log(2)))
    lr <- 2 * (f$loglik[2] - f0$loglik[2]); lrp <- 1 - pchisq(lr, 1)
    h <- tv_hr(f); h[, cohort := nm]; h[, model := lab]; h[, LR_time_p := lrp]
    h[, at_risk := vapply(year, function(y) sum(dd$OS_time >= y * 365.25), integer(1))]
    h[year * 365.25 > max(dd$OS_time[dd$OS_event == 1]), c("HR", "lo", "hi", "p") := NA_real_]
    h[, `:=`(n = nrow(dd), events = sum(dd$OS_event))]
    tv_all[[paste(nm, lab)]] <- h
    note(sprintf("P1 %s %-22s HR/SD at 1/3/5 y: %s | HEATR3 x time LR p=%s", nm, lab,
                 paste(fhr(h$HR, h$lo, h$hi), collapse = " ; "), fp(lrp)))
  }
  # added prognostic value over age+sex+IDH+proliferation
  if (!nrow(d) || !sum(d$OS_event) || !any(is.finite(d$prolif))) next
  d <- d[complete.cases(d[, c("prolif", "H3")]), ]
  if (!nrow(d) || !sum(d$OS_event)) next
  b0 <- coxph(Surv(OS_time, OS_event) ~ age + sex + prolif + strata(IDH), data = d)
  b1 <- update(b0, . ~ . + H3)
  cc2 <- concordance(b0, b1); dC <- as.numeric(diff(cc2$concordance))
  seC <- as.numeric(sqrt(c(-1, 1) %*% cc2$var %*% c(-1, 1)))
  lrp <- anova(b0, b1)$`Pr(>|Chi|)`[2]
  b7[[nm]] <- data.table(cohort = nm, C_base_prolif = cc2$concordance[1], C_plus_HEATR3 = cc2$concordance[2],
                         delta_C = dC, lo = dC - 1.96 * seC, hi = dC + 1.96 * seC, LR_p = lrp)
  note(sprintf("P1 %s added value over age+sex+IDH+prolif: C %.3f -> %.3f (delta %.3f, %.3f to %.3f); LR p=%s",
               nm, cc2$concordance[1], cc2$concordance[2], dC, dC - 1.96 * seC, dC + 1.96 * seC, fp(lrp)))
}
p1 <- rbindlist(p1); tv_all <- rbindlist(tv_all); b7 <- rbindlist(b7)
for (m in unique(p1$model)) {
  x <- p1[model == m & is.finite(b) & is.finite(se) & se > 0]
  if (nrow(x) != 2L) next
  mm <- meta_fe(x$b, x$se)
  p1 <- rbind(p1, data.table(HR = mm["HR"], lo = mm["lo"], hi = mm["hi"], p = mm["p"], b = log(mm["HR"]), se = NA,
                             n = sum(x$n), events = sum(x$events), cohort = "Pooled (fixed effect)", model = m), fill = TRUE)
  note(sprintf("P1 META %-12s HEATR3 HR/SD %s p=%s I2=%.0f%% (Q p=%s)", m, fhr(mm["HR"], mm["lo"], mm["hi"]),
               fp(mm["p"]), mm["I2"], fp(mm["Qp"])))
}
fwrite(p1, file.path(OUT, "P1_HEATR3_proliferation_TCGA_CGGA_meta.csv"))
fwrite(tv_all, file.path(OUT, "P1_time_varying_HR.csv"))
fwrite(b7, file.path(OUT, "P1_added_value_over_proliferation.csv"))
note("P1 constant HR and pooled HR are supplementary because of TCGA non-PH. Compare matched baseline rows.")

###############################################################################
## PART 2a — TABLE 1 FINAL (LGG + GBM n = 257) AND RATIO OF HRs
###############################################################################
message("\n===== PART 2a: Table 1 =====")
t1 <- data.table()
if (file.exists(TCGA_S4)) {
  s4 <- fread(TCGA_S4)[term == "expression"]
  a <- s4[adjusted == TRUE, .(cohort, gene, adj_HR = HR, adj_lo = lower,
                             adj_hi = upper, adj_p = p, n, deaths = events, FDR, PH_term, PH_global)]
  u <- s4[adjusted == FALSE, .(cohort, gene, uni_HR = HR, uni_lo = lower, uni_hi = upper,
                              uni_p = p, uni_n = n, uni_deaths = events)]
  t1 <- merge(a, u, by = c("cohort", "gene"), all = TRUE)
  t1[, adj := ifelse(is.na(adj_HR), NA_character_, fhr(adj_HR, adj_lo, adj_hi))]
  t1[, uni := ifelse(is.na(uni_HR), NA_character_, fhr(uni_HR, uni_lo, uni_hi))]
  fwrite(t1, file.path(OUT, "Table1_final_LGG_GBM.csv"))
  note("P2a Table 1 copied from cohort-specific S4 models: ", TCGA_S4)
} else note("P2a unavailable: run HEATR3_continue.R to generate S4_Cox_Table1.tsv.")
if (file.exists(TCGA_INTERACTION)) {
  ri <- fread(TCGA_INTERACTION)[term == "expression:cohortGBM"]
  fwrite(ri, file.path(OUT, "Table1_HEATR3_cohort_interaction.csv"))
  if (nrow(ri)) note("P2a GBM/LGG ratio of HRs per within-cohort SD: ",
                     fhr(ri$HR, ri$lower, ri$upper), " P=", fp(ri$p))
}

###############################################################################
## PART 2b — GTEx VARIANT IDENTITY (Ensembl REST, GRCh38)
###############################################################################
message("\n===== PART 2b: GTEx identity =====")
if (file.exists(GTEX_HARMONIZED)) {
  h <- fread(GTEX_HARMONIZED)[gene == GENE & !is.na(risk_slope)]
  fwrite(h, file.path(OUT, "P2b_HEATR3_verified_allele_alignment.csv"))
  note("P2b verified HEATR3 risk-allele slopes in S2: ", nrow(h))
}
if (RUN_ENSEMBL && requireNamespace("jsonlite", quietly = TRUE)) {
  ens <- function(path) tryCatch(jsonlite::fromJSON(paste0("https://rest.ensembl.org/", path,
                                 ifelse(grepl("\\?", path), ";", "?"), "content-type=application/json")),
                                 error = function(e) { message("  Ensembl: ", conditionMessage(e)); NULL })
  vids <- character()
  if (file.exists(GTEX_AUDIT)) vids <- unique(fread(GTEX_AUDIT)$variant_id)
  top_pairs <- NULL
  if (have_raw_gtex) {
    eg <- fread(GTEX_EGENES); gid <- eg[gene_name == GENE, gene_id]
    top_pairs <- fread(GTEX_PAIRS)[gene_id %in% gid][order(pval_nominal)]
    vids <- unique(c(vids, head(top_pairs$variant_id, 5)))
  }
  vids <- vids[!is.na(vids) & grepl("^chr[^_]+_[0-9]+_[ACGT]+_[ACGT]+_b38$", vids)]
  # sentinel location + alleles
  sv <- ens(paste0("variation/human/", SENTINEL))
  s_pos <- if (!is.null(sv)) sv$mappings[sv$mappings$assembly_name == "GRCh38", ][1, ] else NULL
  if (!is.null(s_pos)) note(sprintf("P2b sentinel %s: chr%s:%s alleles %s", SENTINEL, s_pos$seq_region_name,
                                    s_pos$start, s_pos$allele_string))
  idt <- rbindlist(lapply(vids, function(v) {
    p <- strsplit(v, "_")[[1]]; chr <- sub("chr", "", p[1]); pos <- as.integer(p[2]); ref <- p[3]; alt <- p[4]
    Sys.sleep(0.15)
    ov <- ens(sprintf("overlap/region/human/%s:%d-%d?feature=variation", chr, pos, pos))
    rs <- NA_character_; al <- NA_character_
    if (is.data.frame(ov) && all(c("id", "alleles", "start", "end") %in% names(ov)) && nrow(ov)) {
      ov <- as.data.frame(ov); ov <- ov[which(grepl("^rs", ov$id) & ov$start == pos & ov$end == pos), , drop = FALSE]
      ok <- vapply(ov$alleles, function(a) all(c(ref, alt) %in% unlist(a)), logical(1))
      if (length(unique(ov$id[ok])) == 1L) { rs <- ov$id[ok][1]; al <- paste(unlist(ov$alleles[ok][1]), collapse = "/") }
    }
    r2 <- NA_real_
    if (!is.na(rs)) {
      if (rs == SENTINEL) r2 <- 1 else {
        Sys.sleep(0.15)
        ld <- ens(sprintf("ld/human/pairwise/%s/%s?population_name=1000GENOMES:phase_3:EUR", SENTINEL, rs))
        if (!is.null(ld) && NROW(ld)) r2 <- as.numeric(ld$r2[1])
      }
    }
    sl <- if (!is.null(top_pairs)) top_pairs[variant_id == v, slope][1] else NA
    pv <- if (!is.null(top_pairs)) top_pairs[variant_id == v, pval_nominal][1] else NA
    data.table(variant_id = v, rsID = rs, ensembl_alleles = al, GTEx_REF = ref, GTEx_ALT = alt,
               slope_ALT = sl, p_eQTL = pv, r2_with_sentinel_EUR = r2,
               identity = ifelse(is.na(rs), "unverified (no rsID with matching alleles)", "verified"))
  }))
  fwrite(idt, file.path(OUT, "P2b_GTEx_variant_identity_verified.csv"))
  for (i in seq_len(nrow(idt))) note(sprintf("P2b %s -> %s (%s), slope ALT=%s, r2 with %s=%s", idt$variant_id[i],
                                             idt$rsID[i], idt$identity[i], ifelse(is.na(idt$slope_ALT[i]), "NA", sprintf("%.2f", idt$slope_ALT[i])),
                                             SENTINEL, ifelse(is.na(idt$r2_with_sentinel_EUR[i]), "NA", sprintf("%.2f", idt$r2_with_sentinel_EUR[i]))))
  sentinel_verified <- if (nrow(idt))
    any(idt$rsID == SENTINEL & is.finite(idt$p_eQTL), na.rm = TRUE) else FALSE
  note("P2b sentinel with verified identity among audited HEATR3 eQTL pairs: ", sentinel_verified)
  note("P2b DECISION RULE: keep 'risk allele increases HEATR3' ONLY if the sentinel is itself an eQTL (then align its ",
       "GWAS risk allele to GTEx ALT) or a verified proxy has r2 >= 0.8 AND allele phase is confirmed (LDlink LDpair). ",
       "Otherwise write: 'the 16q12.1 locus harbours strong cis-eQTLs for HEATR3 in brain cortex' without a direction claim.")
} else note("P2b skipped: install.packages('jsonlite') (needs internet)")

###############################################################################
## PART 2c — PUBLICATION FIGURES (600 dpi TIFF + PDF)
###############################################################################
message("\n===== PART 2c: figures =====")
# Fig 2: eQTL pairs per gene
if (have_raw_gtex) {
  eg <- fread(GTEX_EGENES); pr <- fread(GTEX_PAIRS)
  pr <- merge(pr, eg[, .(gene_id, gene_name)], by = "gene_id")[gene_name %in% GENES]
  cnt <- pr[, .(pairs = .N, top = -log10(min(pval_nominal))), by = gene_name][order(pairs)]
  cnt[, gene_name := factor(gene_name, levels = gene_name)]
  p <- ggplot(cnt, aes(pairs, gene_name)) +
    geom_col(aes(fill = gene_name == GENE), width = 0.65, show.legend = FALSE) +
    geom_text(aes(label = pairs), hjust = -0.2, size = 3.5) +
    scale_fill_manual(values = c(`TRUE` = "#1F5A96", `FALSE` = "grey65")) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(x = "Significant cis-eQTL variant\u2013gene pairs (GTEx v8 brain cortex)", y = NULL,
         title = "Germline regulatory evidence for GWAS-nominated glioma risk genes") +
    theme_pub + theme(axis.text.y = element_text(face = "italic"))
  save_fig(p, "Fig2_eQTL_pairs", 7, 3.6)
}
# Fig 3: forest plot LGG vs GBM (adjusted)
if (nrow(t1) && any(is.finite(t1$adj_p))) {
  fd <- t1[!is.na(adj_p)]
  fd[, `:=`(HR = adj_HR, lo = adj_lo, hi = adj_hi)]
  fd[, cohort := factor(cohort, levels = c("LGG", "GBM"), labels = c("TCGA-LGG", "TCGA-GBM"))]
  fd[, gene := factor(gene, levels = rev(GENES))]
  p <- ggplot(fd, aes(HR, gene)) + geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2) +
    geom_point(aes(colour = !is.na(FDR) & FDR < 0.05), size = 2.6, show.legend = FALSE) +
    scale_colour_manual(values = c(`TRUE` = "#B03A2E", `FALSE` = "grey40")) +
    scale_x_log10() + facet_wrap(~cohort, nrow = 1) +
    labs(x = "Hazard ratio per original expression unit (age- and sex-adjusted, log scale)", y = NULL,
         title = "Risk-gene expression and overall survival: LGG vs GBM",
         subtitle = "Red indicates BH FDR < 0.05 within cohort") +
    theme_pub + theme(axis.text.y = element_text(face = "italic"), strip.background = element_blank(),
                      strip.text = element_text(face = "bold"))
  save_fig(p, "Fig3_forest_LGG_GBM", 8, 3.8)
}
# KM helper (no survminer)
km_gg <- function(d, title, cut) {
  d <- d[is.finite(d$HEATR3), , drop = FALSE]
  d$grp <- factor(ifelse(d$HEATR3 > cut, "High", "Low"), levels = c("Low", "High"))
  if (length(unique(d$grp)) < 2L) return(NULL)
  d$t <- d$OS_time / 365.25
  f <- survfit(Surv(t, OS_event) ~ grp, data = d)
  lr <- survdiff(Surv(t, OS_event) ~ grp, data = d)
  p <- pchisq(lr$chisq, df = 1, lower.tail = FALSE)
  groups <- sub("^grp=", "", names(f$strata))
  curve <- data.frame(
    time = f$time, surv = f$surv,
    grp = rep(groups, f$strata), n.censor = f$n.censor
  )
  origins <- data.frame(time = 0, surv = 1, grp = groups)
  df <- rbind(origins, curve[, c("time", "surv", "grp")])
  cen <- curve[curve$n.censor > 0, , drop = FALSE]
  br <- seq(0, max(2, 2 * ceiling(max(d$t) / 2)), by = 2)
  rk <- summary(f, times = br, extend = TRUE)
  risk <- data.frame(time = rk$time, n = rk$n.risk, grp = sub("grp=", "", rk$strata))
  cols <- c(Low = "#2C7FB8", High = "#C0392B")
  g1 <- ggplot(df, aes(time, surv, colour = grp)) + geom_step(linewidth = 0.8) +
    geom_point(data = cen, shape = 3, size = 1.4, show.legend = FALSE) +
    scale_colour_manual(values = cols, name = "HEATR3") +
    scale_y_continuous(limits = c(0, 1), labels = function(x) paste0(x * 100, "%")) +
    scale_x_continuous(breaks = br, limits = c(0, max(br))) +
    annotate("text", x = max(br) * 0.65, y = 0.92, size = 3.6,
             label = sprintf("Log-rank p = %s\nn: Low %d, High %d", formatC(p, format = "g", digits = 2),
                             sum(d$grp == "Low"), sum(d$grp == "High"))) +
    labs(title = title, x = NULL, y = "Overall survival") + theme_pub
  g2 <- ggplot(risk, aes(time, factor(grp, levels = c("High", "Low")), label = n, colour = grp)) +
    geom_text(size = 3.3, show.legend = FALSE) + scale_colour_manual(values = cols) +
    scale_x_continuous(breaks = br, limits = c(0, max(br))) +
    labs(x = "Years since diagnosis", y = "At risk") + theme_classic(base_size = 10)
  list(g1 = g1, g2 = g2, p = p)
}
save_km <- function(k, name) {
  if (is.null(k)) { note("KM unavailable: ", name); return(invisible(NULL)) }
  if (requireNamespace("cowplot", quietly = TRUE)) {
    p <- cowplot::plot_grid(k$g1, k$g2, ncol = 1, rel_heights = c(3.2, 1), align = "v", axis = "lr")
    save_fig(p, name, 6.5, 6)
  } else { save_fig(k$g1, name, 6.5, 4.8); save_fig(k$g2, paste0(name, "_risktable"), 6.5, 1.6) }
}
for (nm in names(COH)) {
  d <- COH[[nm]]; lab <- gsub("_", "-", nm)
  cut <- median(d$HEATR3, na.rm = TRUE)
  save_km(km_gg(d, paste0(lab, ": all"), cut), paste0("Fig4_KM_", nm, "_all"))
  save_km(km_gg(d[d$IDH %in% "wildtype", ], paste0(lab, ": IDH-wildtype"), cut), paste0("Fig4_KM_", nm, "_IDHwt"))
  save_km(km_gg(d[d$IDH %in% "mutant", ],   paste0(lab, ": IDH-mutant"), cut),   paste0("Fig4_KM_", nm, "_IDHmut"))
}
# Fig 5: time-varying HR, TCGA vs CGGA, with and without proliferation
tvp <- copy(tv_all); tvp[, cohort := gsub("_", "-", cohort)]
p <- ggplot(tvp, aes(factor(year), HR, colour = cohort, shape = cohort)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey50") +
  geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(0.45), size = 0.45) +
  scale_y_log10() + facet_wrap(~model) + scale_colour_manual(values = c("#1F3A5F", "#138D75")) +
  labs(x = "Years of follow-up", y = "HR per 1 SD higher HEATR3 (log scale)",
       title = "HEATR3 hazard ratio over follow-up, before and after adjustment for proliferation",
       colour = NULL, shape = NULL) +
  theme_pub + theme(strip.background = element_blank(), strip.text = element_text(face = "bold"))
save_fig(p, "Fig5C_time_varying_HR_TCGA_CGGA", 8, 4)
# Fig 5 forest: base vs +prolif, per cohort and pooled
fd <- copy(p1[model %in% c("base", "plus_MKI67", "plus_prolif")]); fd[, model := factor(model, levels = c("base", "plus_MKI67", "plus_prolif"),
                                      labels = c("Age + sex, IDH-stratified", "+ MKI67", "+ proliferation score"))]
fd[, cohort := factor(gsub("_", "-", cohort), levels = rev(c("TCGA-LGG", "CGGA-693", "Pooled (fixed effect)")))]
p <- ggplot(fd, aes(HR, cohort)) + geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2) +
  geom_point(aes(shape = cohort == "Pooled (fixed effect)"), size = 2.8, show.legend = FALSE) +
  scale_shape_manual(values = c(`TRUE` = 18, `FALSE` = 16)) + scale_x_log10() + facet_wrap(~model, nrow = 1) +
  geom_text(aes(label = sprintf("%.2f", HR)), vjust = -1, size = 3) +
  labs(x = "HR per 1 SD higher HEATR3 (log scale)", y = NULL,
       title = "HEATR3 and overall survival in LGG: exploratory constant-HR models") +
  theme_pub + theme(strip.background = element_blank(), strip.text = element_text(face = "bold"))
save_fig(p, "Fig5_forest_proliferation_TCGA_CGGA_meta", 9, 3.4)
if (!have_raw_gtex)
  note("P2c GTEx eQTL count figure skipped: check GTEx input paths and R.utils.")
note("P2c figures written to ", OUT, " (PDF vector + ", DPI, " dpi TIFF)")

###############################################################################
## PART 3 — GSEA
###############################################################################
message("\n===== PART 3: GSEA =====")
sets <- NULL
if (requireNamespace("msigdbr", quietly = TRUE)) {
  get_sets <- function(coll, sub = NULL) {
    m <- tryCatch(msigdbr::msigdbr(species = "Homo sapiens", collection = coll, subcollection = sub),
                  error = function(e) msigdbr::msigdbr(species = "Homo sapiens", category = coll, subcategory = sub))
    split(m$gene_symbol, m$gs_name)
  }
  sets <- tryCatch(c(get_sets("H"), get_sets("C2", "CP:REACTOME")), error = function(e) NULL)
}
if (is.null(sets) && requireNamespace("fgsea", quietly = TRUE) && any(file.exists(GMT_FILES)))
  sets <- do.call(c, lapply(GMT_FILES[file.exists(GMT_FILES)], fgsea::gmtPathways))
if (!is.null(sets) && requireNamespace("fgsea", quietly = TRUE)) {
  EX <- list(TCGA_LGG = tcga_ex, CGGA_693 = cgga_ex); gs <- list()
  for (nm in names(COH)) {
    d <- COH[[nm]]; ok <- !is.na(d$IDH); d <- d[ok, ]
    ex <- EX[[nm]][, d$id, drop = FALSE]
    ex <- ex[which(rowMeans(ex) > 1 & apply(ex, 1, sd) > 0.1), , drop = FALSE]
    X  <- model.matrix(~ IDH, data = d)
    R  <- ex - t(X %*% solve(crossprod(X), crossprod(X, t(ex))))        # residualise every gene on IDH (fast)
    rh <- resid(lm(d$H3 ~ d$IDH))
    rho <- apply(R, 1, function(x) cor(x, rh, method = "spearman"))
    rho <- rho[names(rho) != GENE & is.finite(rho)] + rnorm(sum(names(rho) != GENE & is.finite(rho)), 0, 1e-9)
    fg <- as.data.table(fgsea::fgsea(sets, sort(rho, decreasing = TRUE), minSize = 15, maxSize = 500))
    fg[, leadingEdge := vapply(leadingEdge, function(x) paste(head(x, 15), collapse = ","), "")]
    fg[, cohort := nm]; gs[[nm]] <- fg[order(padj)]
    up <- fg[padj < 0.05 & NES > 0][order(-NES)]
    note(sprintf("P3 %s: %d pathways FDR<0.05; top positive: %s", nm, nrow(fg[padj < 0.05]),
                 paste(head(up$pathway, 6), collapse = "; ")))
    foc <- fg[grepl("RIBOSOM|RRNA|MYC|P53|NUCLEOL|TRANSLATION|E2F|G2M|MITOTIC", pathway)][order(padj)][seq_len(min(.N, 10L))]
    for (i in seq_len(nrow(foc))) note(sprintf("P3 %s focus %-55s NES %5.2f FDR %s", nm, foc$pathway[i], foc$NES[i], fp(foc$padj[i])))
  }
  gs <- rbindlist(gs); fwrite(gs, file.path(OUT, "P3_GSEA_HEATR3_IDH_adjusted.csv"))
  w <- dcast(gs, pathway ~ cohort, value.var = c("NES", "padj"))
  both <- w[padj_TCGA_LGG < 0.05 & padj_CGGA_693 < 0.05 & sign(NES_TCGA_LGG) == sign(NES_CGGA_693)]
  fwrite(both[order(-NES_TCGA_LGG)], file.path(OUT, "P3_GSEA_replicated_both.csv"))
  note(sprintf("P3 pathways FDR<0.05 in both cohorts, same direction: %d", nrow(both)))
  if (nrow(both)) {
    top <- head(both[order(-abs(NES_TCGA_LGG))]$pathway, 20)
    pd <- gs[pathway %in% top]; pd[, pathway := factor(gsub("_", " ", pathway), levels = rev(gsub("_", " ", top)))]
    pd[, cohort := gsub("_", "-", cohort)]
    p <- ggplot(pd, aes(NES, pathway, colour = cohort, size = -log10(padj))) + geom_point(alpha = 0.85) +
      geom_vline(xintercept = 0, colour = "grey60") + scale_colour_manual(values = c("#138D75", "#1F3A5F")) +
      labs(x = "Normalised enrichment score (genes ranked by IDH-adjusted correlation with HEATR3)", y = NULL,
           size = expression(-log[10]~FDR), colour = NULL,
           title = "Pathways co-varying with HEATR3 in both LGG cohorts") +
      theme_pub + theme(axis.text.y = element_text(size = 7.5))
    save_fig(p, "FigS4_GSEA_HEATR3_both_cohorts", 10, 7)
  }
} else note("P3 skipped: install msigdbr (install.packages('msigdbr')) or download MSigDB .gmt files (see CONFIG)")

###############################################################################
writeLines(c(paste("Run:", Sys.time()), paste("R", getRversion(), "| survival", packageVersion("survival")),
             "", LOG), file.path(OUT, "SUMMARY_FINAL.txt"))
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
message("\nDone -> ", file.path(OUT, "SUMMARY_FINAL.txt"))