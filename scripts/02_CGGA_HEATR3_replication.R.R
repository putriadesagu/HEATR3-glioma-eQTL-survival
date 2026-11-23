# CGGA external replication of the TCGA-LGG HEATR3 Cox model.
# Open this file in RStudio and click Source.
# If needed, run install.packages("survival") once.
ROOT <- "D:/RESIDENSI/Thesis-PPDS/R"
RELEASE <- "693"
CLINICAL_FILE <- ""    # Blank opens a file chooser.
EXPRESSION_FILE <- ""  # Blank opens a file chooser.
EXPRESSION_TRANSFORM <- "as_provided"  # Or "log2p1" if appropriate for linear FPKM.

run_cgga <- function() {
  if (!requireNamespace("survival", quietly = TRUE)) stop("Install the survival package.")
  if (!dir.exists(ROOT)) stop("Project folder not found. Check ROOT.")
  if (!EXPRESSION_TRANSFORM %in% c("as_provided", "log2p1"))
    stop("EXPRESSION_TRANSFORM must be as_provided or log2p1.")
  library(survival)

  choose_input <- function(path, label) {
    if (!nzchar(path)) {
      message("Pilih file ", label, " CGGA ", RELEASE, " yang sudah diekstrak.")
      path <- file.choose()
    }
    if (!file.exists(path) || !grepl("\\.(txt|tsv|csv)(\\.gz)?$", path, ignore.case = TRUE))
      stop("Expected an extracted TXT, TSV or CSV file for ", label, ": ", path)
    path
  }

  read_input <- function(path) {
    delim <- if (grepl("\\.csv(\\.gz)?$", path, ignore.case = TRUE)) "," else "\t"
    read.table(path, header = TRUE, sep = delim, quote = "\"", comment.char = "",
               check.names = FALSE, stringsAsFactors = FALSE,
               na.strings = c("", "NA", "N/A", "--"))
  }

  number <- function(x, label) {
    x <- trimws(as.character(x))
    x[x %in% c("", "NA", "N/A", "--")] <- NA_character_
    values <- suppressWarnings(as.numeric(x))
    if (any(!is.na(x) & !is.finite(values))) stop("Invalid numeric values in ", label)
    values
  }

  clinical_path <- choose_input(CLINICAL_FILE, "KLINIS")
  clinical <- read_input(clinical_path)
  names(clinical) <- trimws(sub("^\ufeff", "", names(clinical)))
  columns <- c(id = "CGGA_ID", type = "PRS_type", grade = "Grade",
               os = "OS", event = "Censor (alive=0; dead=1)", age = "Age", sex = "Gender")
  missing <- setdiff(unname(columns), names(clinical))
  if (length(missing))
    stop("Kolom klinis tidak ditemukan: ", paste(missing, collapse = ", "),
         "\nHeader: ", paste(names(clinical), collapse = ", "))
  d <- clinical[, unname(columns), drop = FALSE]
  names(d) <- names(columns)
  d$raw_idh <- clinical[["IDH_mutation_status"]]
  d$id <- trimws(as.character(d$id))
  tumor_type <- tolower(trimws(d$type))
  grade <- toupper(trimws(d$grade))
  grade <- sub("^(CNS[[:space:]]+)?WHO[[:space:]]*", "", grade)
  grade <- sub("^GRADE[[:space:]]*", "", grade)
  print(table(tumor_type, grade, useNA = "ifany"))
  d <- d[which(tumor_type == "primary" & grade %in% c("II", "III", "2", "3")), , drop = FALSE]
  if (!nrow(d)) stop("No primary grade II–III samples. Check the type and grade table.")
  if (anyNA(d$id) || any(!nzchar(d$id)) || anyDuplicated(d$id))
    stop("Primary grade II–III clinical IDs must be present and unique.")

  # The clinical header explicitly defines 0 = alive/censored, 1 = dead.
  d$event <- number(d$event, columns[["event"]])
  d$os <- number(d$os, "OS")
  d$age <- number(d$age, "Age")
  if (any(!is.na(d$event) & !d$event %in% 0:1))
    stop("Event codes must be 0 (censored) or 1 (death).")
  if (any(d$os < 0, na.rm = TRUE) || any(d$age < 0 | d$age > 120, na.rm = TRUE))
    stop("Check OS and age values.")
  print(table(event = d$event, useNA = "ifany"))

  sex <- tolower(trimws(d$sex))
  sex[sex %in% "f"] <- "female"
  sex[sex %in% "m"] <- "male"
  d$sex <- factor(sex, levels = c("female", "male"))

  expression_path <- choose_input(EXPRESSION_FILE, "EKSPRESI RSEM")
  expression <- read_input(expression_path)
  sample_ids <- trimws(names(expression)[-1])
  if (!length(sample_ids) || anyNA(sample_ids) || any(!nzchar(sample_ids)) ||
      anyDuplicated(sample_ids)) stop("Check sample ID columns in the expression matrix.")
  heatr3_row <- which(toupper(trimws(expression[[1]])) == "HEATR3")
  if (length(heatr3_row) != 1L)
    stop("Expected one HEATR3 row with the gene symbol in the first column.")
  values <- number(unlist(expression[heatr3_row, -1, drop = FALSE], use.names = FALSE),
                   "HEATR3 expression")
  idx <- match(d$id, sample_ids)
  matched <- sum(!is.na(idx))
  if (!matched) stop("No matched IDs. Select clinical and expression files from the same release.")
  d$HEATR3 <- values[idx]
  if (EXPRESSION_TRANSFORM == "log2p1") {
    if (any(d$HEATR3 < 0, na.rm = TRUE)) stop("log2p1 requires nonnegative expression.")
    d$HEATR3 <- log2(d$HEATR3 + 1)
  }

  covariates <- c("os", "event", "age", "sex", "HEATR3")
  analysis <- droplevels(d[complete.cases(d[, covariates]), , drop = FALSE])
  if (nrow(analysis) < 2L || !sum(analysis$event) ||
      nlevels(analysis$sex) < 2L || length(unique(analysis$age)) < 2L ||
      sd(analysis$HEATR3) == 0) stop("Insufficient complete cases, events, or predictor variation.")
  analysis$HEATR3_SD <- as.numeric(scale(analysis$HEATR3))

  unadjusted <- coxph(Surv(os, event) ~ HEATR3_SD, data = analysis, x = TRUE)
  adjusted <- coxph(Surv(os, event) ~ HEATR3_SD + age + sex, data = analysis, x = TRUE)
  summarize <- function(fit, model) {
    s <- summary(fit)
    term <- "HEATR3_SD"
    data.frame(
      cohort = paste0("CGGA_", RELEASE), model = model,
      expression_transform = EXPRESSION_TRANSFORM, n = fit$n, events = fit$nevent,
      beta = s$coefficients[term, "coef"], SE = s$coefficients[term, "se(coef)"],
      HR_per_SD = s$conf.int[term, "exp(coef)"],
      CI95_low = s$conf.int[term, "lower .95"], CI95_high = s$conf.int[term, "upper .95"],
      p = s$coefficients[term, "Pr(>|z|)"], C_index = unname(s$concordance[1])
    )
  }
  results <- rbind(summarize(unadjusted, "Unadjusted, same complete cases"),
                   summarize(adjusted, "Adjusted for age and sex"))
  flow <- data.frame(
    clinical_file = basename(clinical_path), expression_file = basename(expression_path),
    clinical_rows = nrow(clinical), primary_grade_II_III = nrow(d),
    matched_to_expression = matched, complete_cases = nrow(analysis),
    deaths = sum(analysis$event)
  )
  ph <- cox.zph(adjusted)$table
  ph <- data.frame(term = rownames(ph), ph, row.names = NULL)
  out <- file.path(ROOT, "results", paste0("CGGA_", RELEASE, "_HEATR3"))
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  write.csv(results, file.path(out, "CGGA_HEATR3_Cox.csv"), row.names = FALSE)
  write.csv(flow, file.path(out, "CGGA_HEATR3_sample_flow.csv"), row.names = FALSE)
  write.csv(ph, file.path(out, "CGGA_HEATR3_PH_test.csv"), row.names = FALSE)
  print(flow, row.names = FALSE)
  print(results, row.names = FALSE, digits = 4)
  idh_column <- "IDH_mutation_status"
  if (!idh_column %in% names(clinical)) {
    warning("IDH analysis unavailable: missing ", idh_column, call. = FALSE)
  } else {
    raw_idh <- trimws(as.character(analysis$raw_idh))
    key <- gsub("[^a-z0-9]", "", tolower(raw_idh))
    idh <- rep(NA_character_, length(key))
    idh[key %in% c("mutant", "mutated", "mut", "idhmutant", "idh1mutant")] <- "Mutant"
    idh[key %in% c("wildtype", "wt", "idhwildtype", "idh1wildtype")] <- "Wild-type"
    analysis$idh <- factor(idh, levels = c("Mutant", "Wild-type"))
    counts <- as.data.frame(table(raw_IDH = raw_idh, mapped_IDH = analysis$idh, useNA = "ifany"))
    counts <- counts[counts$Freq > 0, , drop = FALSE]
    write.csv(counts, file.path(out, "CGGA_HEATR3_IDH_coding.csv"), row.names = FALSE)
    print(counts, row.names = FALSE)

    run_idh <- function() {
      subset <- droplevels(analysis[!is.na(analysis$idh), , drop = FALSE])
      if (nlevels(subset$idh) != 2L || !sum(subset$event) ||
          nlevels(subset$sex) != 2L || length(unique(subset$age)) < 2L ||
          length(unique(subset$HEATR3)) < 2L)
        stop("Insufficient cases, events or predictor variation after IDH filtering.")
      # Keep the SD from the original age/sex complete-case cohort in both models.
      base <- coxph(Surv(os, event) ~ HEATR3_SD + age + sex,
                    data = subset, x = TRUE, singular.ok = FALSE)
      full <- coxph(Surv(os, event) ~ HEATR3_SD + age + sex + idh,
                    data = subset, x = TRUE, singular.ok = FALSE)
      comparison <- rbind(summarize(base, "Age and sex, IDH complete cases"),
                          summarize(full, "Age, sex and IDH"))
      idh_flow <- data.frame(
        age_sex_complete = nrow(analysis), idh_complete = nrow(subset),
        missing_or_unmapped_idh = sum(is.na(analysis$idh)),
        mutant = sum(subset$idh == "Mutant"),
        wild_type = sum(subset$idh == "Wild-type"), deaths = sum(subset$event)
      )
      write.csv(comparison, file.path(out, "CGGA_HEATR3_IDH_Cox.csv"), row.names = FALSE)
      write.csv(idh_flow, file.path(out, "CGGA_HEATR3_IDH_sample_flow.csv"), row.names = FALSE)
      ph_idh <- cox.zph(full)$table
      write.csv(data.frame(term = rownames(ph_idh), ph_idh, row.names = NULL),
                file.path(out, "CGGA_HEATR3_IDH_PH_test.csv"), row.names = FALSE)
      print(idh_flow, row.names = FALSE)
      print(comparison[, c("model", "n", "events", "HR_per_SD", "CI95_low",
                            "CI95_high", "p")], row.names = FALSE, digits = 4)

      fit_stratified <- function(data) {
        coxph(Surv(os, event) ~ HEATR3_SD + age + sex + strata(idh),
              data = data, x = TRUE, singular.ok = FALSE)
      }
      save_ph <- function(fit, filename) {
        tab <- cox.zph(fit)$table
        write.csv(data.frame(term = rownames(tab), tab, row.names = NULL),
                  file.path(out, filename), row.names = FALSE)
      }
      model_label <- "Age + sex, stratified by IDH"
      fit_cgga <- fit_stratified(subset)
      cgga <- summarize(fit_cgga, model_label)
      write.csv(rbind(comparison, cgga),
                file.path(out, "CGGA_HEATR3_IDH_model_comparison.csv"), row.names = FALSE)
      save_ph(fit_cgga, "CGGA_HEATR3_IDH_stratified_PH.csv")
      print(cgga, row.names = FALSE)

      # Match the TCGA time-varying model. CGGA OS is interpreted in days.
      tryCatch({
        tv <- subset
        tv$age10 <- (tv$age - mean(tv$age)) / 10
        time_basis <- function(days) log1p(days / 365.25) - log(2)
        time_interaction <- function(x, t, ...) x * time_basis(t)
        age_tv <- coxph(
          Surv(os, event) ~ HEATR3_SD + age10 + sex + strata(idh) + tt(age10),
          data = tv, tt = time_interaction, singular.ok = FALSE
        )
        both_tv <- coxph(
          Surv(os, event) ~ HEATR3_SD + age10 + sex + strata(idh) +
            tt(age10) + tt(HEATR3_SD),
          data = tv, tt = time_interaction, singular.ok = FALSE
        )
        ll_age <- logLik(age_tv)
        ll_both <- logLik(both_tv)
        test_stat <- 2 * (as.numeric(ll_both) - as.numeric(ll_age))
        test_df <- attr(ll_both, "df") - attr(ll_age, "df")
        test <- data.frame(
          n = both_tv$n, events = both_tv$nevent, chisq = test_stat, df = test_df,
          p = pchisq(test_stat, df = test_df, lower.tail = FALSE)
        )

        years <- c(1, 3, 5)
        days <- years * 365.25
        at_risk <- vapply(days, function(t) sum(tv$os >= t), integer(1))
        terms <- c("HEATR3_SD", "tt(HEATR3_SD)")
        design <- cbind(1, time_basis(days))
        estimate <- drop(design %*% coef(both_tv)[terms])
        variance <- vcov(both_tv)[terms, terms, drop = FALSE]
        se <- sqrt(diag(design %*% variance %*% t(design)))
        # Suppress extrapolated HRs beyond the last observed event.
        supported <- at_risk > 0 & days <= max(tv$os[tv$event == 1])
        estimate[!supported] <- NA_real_
        result <- data.frame(
          year = years, patients_at_risk = at_risk,
          HR_per_SD = exp(estimate),
          CI95_low = exp(estimate - qnorm(.975) * se),
          CI95_high = exp(estimate + qnorm(.975) * se)
        )
        write.csv(result, file.path(out, "CGGA_693_HEATR3_time_varying_HR.csv"),
                  row.names = FALSE)
        write.csv(test, file.path(out, "CGGA_693_HEATR3_time_varying_LRT.csv"),
                  row.names = FALSE)
        message("CGGA exploratory HEATR3-by-time analysis")
        print(test, row.names = FALSE, digits = 4)
        print(result, row.names = FALSE, digits = 4)
      }, error = function(e) message("CGGA time-varying model: ", conditionMessage(e)))

      wt <- droplevels(subset[subset$idh == "Wild-type", , drop = FALSE])
      message("CGGA IDH-wildtype: n = ", nrow(wt), "; deaths = ", sum(wt$event))
      tryCatch({
        fit_wt <- coxph(Surv(os, event) ~ HEATR3_SD + age + sex,
                        data = wt, x = TRUE, singular.ok = FALSE)
        wt_result <- summarize(fit_wt, "IDH-wildtype: age + sex")
        write.csv(wt_result, file.path(out, "CGGA_HEATR3_IDH_wildtype_Cox.csv"),
                  row.names = FALSE)
        save_ph(fit_wt, "CGGA_HEATR3_IDH_wildtype_PH.csv")
        print(wt_result, row.names = FALSE)
      }, error = function(e) message("Wild-type model: ", conditionMessage(e)))

      tcga_path <- file.path(ROOT, "TCGA_LGG_merged_selected_genes_with_IDH1.csv")
      tcga_raw <- read_input(tcga_path)
      required <- c("patient_id", "OS_time", "OS_event", "age_at_diagnosis",
                    "gender", "IDH1_status", "HEATR3")
      missing <- setdiff(required, names(tcga_raw))
      if (length(missing)) stop("Missing TCGA columns: ", paste(missing, collapse = ", "))
      ids <- trimws(as.character(tcga_raw$patient_id))
      if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids))
        stop("TCGA patient IDs must be valid and unique.")

      tcga_sex <- tolower(trimws(tcga_raw$gender))
      tcga_sex[tcga_sex %in% "f"] <- "female"
      tcga_sex[tcga_sex %in% "m"] <- "male"
      tcga <- data.frame(
        os = number(tcga_raw$OS_time, "TCGA OS_time"),
        event = number(tcga_raw$OS_event, "TCGA OS_event"),
        age = number(tcga_raw$age_at_diagnosis, "TCGA age_at_diagnosis"),
        sex = factor(tcga_sex, levels = c("female", "male")),
        HEATR3 = number(tcga_raw$HEATR3, "TCGA HEATR3"),
        idh = factor(tolower(trimws(tcga_raw$IDH1_status)),
                     levels = c("idh1-mutant", "idh1-wildtype"),
                     labels = c("Mutant", "Wild-type"))
      )
      if (any(!is.na(tcga$event) & !tcga$event %in% 0:1) ||
          any(tcga$os < 0, na.rm = TRUE) ||
          any(tcga$age < 0 | tcga$age > 120, na.rm = TRUE))
        stop("Check TCGA event codes, survival time and age units.")
      core <- c("os", "event", "age", "sex", "HEATR3")
      tcga <- tcga[complete.cases(tcga[, core]), , drop = FALSE]
      tcga_sd <- sd(tcga$HEATR3)
      if (!is.finite(tcga_sd) || tcga_sd <= 0) stop("Invalid TCGA HEATR3 SD.")
      tcga$HEATR3_SD <- as.numeric(scale(tcga$HEATR3))
      tcga <- droplevels(tcga[!is.na(tcga$idh), , drop = FALSE])
      if (nlevels(tcga$idh) != 2L) stop("Expected both TCGA IDH strata.")

      fit_tcga <- fit_stratified(tcga)
      tcga_result <- summarize(fit_tcga, model_label)
      tcga_result$cohort <- "TCGA_LGG"
      tcga_result$expression_transform <- "TCGA merged file, as provided"
      save_ph(fit_tcga, "TCGA_LGG_HEATR3_IDH_stratified_PH.csv")

      meta_input <- rbind(tcga_result, cgga)
      if (any(!is.finite(meta_input$beta)) ||
          any(!is.finite(meta_input$SE) | meta_input$SE <= 0))
        stop("Invalid coefficient or SE for pooling.")
      weights <- 1 / meta_input$SE^2
      beta <- weighted.mean(meta_input$beta, weights)
      se <- sqrt(1 / sum(weights))
      Q <- sum(weights * (meta_input$beta - beta)^2)
      df <- nrow(meta_input) - 1L
      meta <- data.frame(
        model = "Fixed-effect inverse-variance; both cohorts stratified by IDH",
        cohorts = nrow(meta_input), HR_per_SD = exp(beta),
        CI95_low = exp(beta - qnorm(0.975) * se),
        CI95_high = exp(beta + qnorm(0.975) * se),
        p = 2 * pnorm(-abs(beta / se)), Q = Q,
        p_heterogeneity = pchisq(Q, df, lower.tail = FALSE),
        I2_percent = if (Q > 0) 100 * max(0, (Q - df) / Q) else 0
      )
      write.csv(meta_input, file.path(out, "TCGA_CGGA_HEATR3_IDH_meta_input.csv"),
                row.names = FALSE)
      write.csv(meta, file.path(out, "TCGA_CGGA_HEATR3_IDH_meta_fixed.csv"),
                row.names = FALSE)
      message("TCGA HEATR3 SD before IDH filtering: ", signif(tcga_sd, 5))
      print(meta_input, row.names = FALSE)
      print(meta, row.names = FALSE)

    }
    tryCatch(run_idh(), error = function(e)
      warning("IDH analysis: ", conditionMessage(e), call. = FALSE))
  }
  message("Hasil tersimpan di ", out)
}

run_cgga()
