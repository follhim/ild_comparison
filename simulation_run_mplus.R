# =============================================================================
# simulation_run_mplus.R
# Mplus DSEM simulation using MplusAutomation
# Mirrors simulation_run.R (Blimp) with identical data-generation functions
# and the same N × T conditions.
#
# USAGE:
#   source(here::here("simulation_run_mplus.R"))   # load functions
#   plan(multisession, workers = 4)                # set workers AFTER source
#   results_sim1 <- run_simulation_mplus_parallel(sim_num = 1, n_reps = 1000)
#
# REQUIREMENTS:
#   - Mplus installed and on system PATH (or set options(MplusAutomation.Mplus.exe = "..."))
#   - Data generation functions loaded (they live in analysis.qmd; source that
#     file's generate_* functions, or copy them here)
# =============================================================================

library(tidyverse)
library(here)
library(furrr)
library(future)
library(progressr)
library(MplusAutomation)

# ── Mplus model specifications ────────────────────────────────────────────────
# Each list element matches one simulation (keyed "1"–"4")
# VARIABLE and MODEL blocks are Mplus syntax strings.

.mplus_variable_block <- "
    CLUSTER = id;
    LAGGED = craving_smoking(1) stressed(1);
    TINTERVAL = time(1);"

.mplus_models <- list(

  # ── Sim 1: AR(1) DSEM ──────────────────────────────────────────────────────
  `1` = list(
    VARIABLE = .mplus_variable_block,
    MODEL = "
      %WITHIN%
        s_AR | craving_smoking ON craving_smoking&1;
        s_CL | craving_smoking ON stressed&1;

      %BETWEEN%
        craving_smoking stressed;
        [s_AR s_CL];
        s_AR s_CL craving_smoking stressed WITH
        s_AR s_CL craving_smoking stressed;"
  ),

  # ── Sim 2: VAR(1) DSEM ─────────────────────────────────────────────────────
  `2` = list(
    VARIABLE = .mplus_variable_block,
    MODEL = "
      %WITHIN%
        s_AR1 | stressed        ON stressed&1;
        s_CL1 | stressed        ON craving_smoking&1;
        s_AR2 | craving_smoking ON craving_smoking&1;
        s_CL2 | craving_smoking ON stressed&1;

      %BETWEEN%
        craving_smoking stressed;
        [s_AR1 s_CL1 s_AR2 s_CL2];
        s_AR1 s_CL1 s_AR2 s_CL2 craving_smoking stressed WITH
        s_AR1 s_CL1 s_AR2 s_CL2 craving_smoking stressed;"
  ),

  # ── Sim 3: AR(1) RDSEM ─────────────────────────────────────────────────────
  `3` = list(
    VARIABLE = .mplus_variable_block,
    MODEL = "
      %WITHIN%
        s_AR    | craving_smoking ON craving_smoking&1;
        s_CL    | craving_smoking ON stressed&1;
        s_trend | craving_smoking ON time;

      %BETWEEN%
        craving_smoking stressed;
        [s_AR s_CL s_trend];
        s_AR s_CL s_trend craving_smoking stressed WITH
        s_AR s_CL s_trend craving_smoking stressed;"
  ),

  # ── Sim 4: VAR(1) RDSEM ────────────────────────────────────────────────────
  `4` = list(
    VARIABLE = .mplus_variable_block,
    MODEL = "
      %WITHIN%
        s_AR1     | stressed        ON stressed&1;
        s_CL1     | stressed        ON craving_smoking&1;
        s_AR2     | craving_smoking ON craving_smoking&1;
        s_CL2     | craving_smoking ON stressed&1;
        s_trend_s | stressed        ON time;
        s_trend_c | craving_smoking ON time;

      %BETWEEN%
        craving_smoking stressed;
        [s_AR1 s_CL1 s_AR2 s_CL2 s_trend_s s_trend_c];
        s_AR1 s_CL1 s_AR2 s_CL2 s_trend_s s_trend_c
          craving_smoking stressed WITH
        s_AR1 s_CL1 s_AR2 s_CL2 s_trend_s s_trend_c
          craving_smoking stressed;"
  )
)

# ── Extraction helpers ────────────────────────────────────────────────────────

# Extract from a readModels() result object.
# Returns a one-row tibble with the same column names used by the Blimp
# extraction functions, so simulation_tables.R works on both.
.extract_mplus <- function(res, sim_num) {
  params <- res$parameters$unstandardized
  if (is.null(params) || nrow(params) == 0) return(NULL)

  get_mean_ci <- function(p) {
    row <- params %>%
      filter(paramHeader == "Means", param == toupper(p))
    if (nrow(row) == 0)
      return(list(est = NA_real_, lo = NA_real_, hi = NA_real_))
    list(est = row$est[1],
         lo  = row$lower_2.5ci[1],
         hi  = row$upper_2.5ci[1])
  }

  get_var_est <- function(p) {
    row <- params %>%
      filter(paramHeader == "Variances", param == toupper(p))
    if (nrow(row) == 0) return(NA_real_)
    row$est[1]
  }

  switch(
    as.character(sim_num),

    "1" = {
      ar <- get_mean_ci("s_ar"); cl <- get_mean_ci("s_cl")
      tibble(
        AR_mean    = ar$est, AR_mean_lo = ar$lo, AR_mean_hi = ar$hi,
        CL_mean    = cl$est, CL_mean_lo = cl$lo, CL_mean_hi = cl$hi,
        AR_var     = get_var_est("s_ar"),
        CL_var     = get_var_est("s_cl")
      )
    },

    "2" = {
      ar1 <- get_mean_ci("s_ar1"); cl1 <- get_mean_ci("s_cl1")
      ar2 <- get_mean_ci("s_ar2"); cl2 <- get_mean_ci("s_cl2")
      tibble(
        AR_1_mean = ar1$est, AR_1_mean_lo = ar1$lo, AR_1_mean_hi = ar1$hi,
        AR_2_mean = ar2$est, AR_2_mean_lo = ar2$lo, AR_2_mean_hi = ar2$hi,
        CL_1_mean = cl1$est, CL_1_mean_lo = cl1$lo, CL_1_mean_hi = cl1$hi,
        CL_2_mean = cl2$est, CL_2_mean_lo = cl2$lo, CL_2_mean_hi = cl2$hi,
        AR_1_var  = get_var_est("s_ar1"),
        AR_2_var  = get_var_est("s_ar2"),
        CL_1_var  = get_var_est("s_cl1"),
        CL_2_var  = get_var_est("s_cl2")
      )
    },

    "3" = {
      ar <- get_mean_ci("s_ar"); cl <- get_mean_ci("s_cl")
      tr <- get_mean_ci("s_trend")
      tibble(
        AR_mean      = ar$est, AR_mean_lo = ar$lo, AR_mean_hi = ar$hi,
        CL_mean      = cl$est, CL_mean_lo = cl$lo, CL_mean_hi = cl$hi,
        trend_mean   = tr$est, trend_mean_lo = tr$lo, trend_mean_hi = tr$hi,
        AR_var       = get_var_est("s_ar"),
        CL_var       = get_var_est("s_cl"),
        trend_var    = get_var_est("s_trend")
      )
    },

    "4" = {
      ar1 <- get_mean_ci("s_ar1"); cl1 <- get_mean_ci("s_cl1")
      ar2 <- get_mean_ci("s_ar2"); cl2 <- get_mean_ci("s_cl2")
      trc <- get_mean_ci("s_trend_c"); trs <- get_mean_ci("s_trend_s")
      tibble(
        AR_1_mean           = ar1$est, AR_1_mean_lo = ar1$lo, AR_1_mean_hi = ar1$hi,
        AR_2_mean           = ar2$est, AR_2_mean_lo = ar2$lo, AR_2_mean_hi = ar2$hi,
        CL_1_mean           = cl1$est, CL_1_mean_lo = cl1$lo, CL_1_mean_hi = cl1$hi,
        CL_2_mean           = cl2$est, CL_2_mean_lo = cl2$lo, CL_2_mean_hi = cl2$hi,
        craving_trend_mean  = trc$est, craving_trend_mean_lo = trc$lo, craving_trend_mean_hi = trc$hi,
        stressed_trend_mean = trs$est, stressed_trend_mean_lo = trs$lo, stressed_trend_mean_hi = trs$hi,
        AR_1_var            = get_var_est("s_ar1"),
        AR_2_var            = get_var_est("s_ar2"),
        CL_1_var            = get_var_est("s_cl1"),
        CL_2_var            = get_var_est("s_cl2"),
        craving_trend_var   = get_var_est("s_trend_c"),
        stressed_trend_var  = get_var_est("s_trend_s")
      )
    },

    stop("sim_num must be 1–4")
  )
}

# ── Main simulation function ──────────────────────────────────────────────────

run_simulation_mplus_parallel <- function(
    sim_num,
    N_vec    = c(50, 100, 200),
    T_vec    = c(30, 50),
    n_reps   = 1000,
    biter    = 50000,    # total MCMC iterations per chain
    warmup   = 10000,    # burn-in (BITERATIONS minimum)
    chains   = 2,        # MCMC chains per model (2 keeps each run fast)
    procs    = 2,        # Mplus PROCESSORS per model (chains × procs ≤ cores/workers)
    save_dir = here("simulation_results")
) {

  dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)

  # -- Select data-generation function
  local_gen_fn <- switch(
    as.character(sim_num),
    "1" = generate_ar1_dsem,
    "2" = generate_var1_dsem,
    "3" = generate_ar1_rdsem,
    "4" = generate_var1_rdsem,
    stop("sim_num must be 1–4")
  )

  local_sim      <- sim_num
  local_model    <- .mplus_models[[as.character(sim_num)]]
  local_analysis <- sprintf(
    "TYPE = TWOLEVEL RANDOM;\n    ESTIMATOR = BAYES;\n    BITERATIONS = %d (%d);\n    CHAINS = %d;\n    PROCESSORS = %d;",
    biter, warmup, chains, procs
  )

  all_results <- list()
  n_cond      <- length(N_vec) * length(T_vec)
  cond_idx    <- 0L

  handlers(handler_progress(
    format = "[:bar] :percent | :current/:total reps | ETA :eta | Elapsed :elapsed | :message",
    width  = 80,
    clear  = FALSE
  ))

  for (N in N_vec) {
    for (T in T_vec) {
      cond_idx <- cond_idx + 1L
      cat(sprintf(
        "\n=== Mplus Sim %d | N=%d T=%d (%d/%d) | %s ===\n",
        sim_num, N, T, cond_idx, n_cond, format(Sys.time(), "%H:%M:%S")
      ))
      start_t  <- proc.time()
      local_N  <- N
      local_T  <- T

      with_progress({
        p <- progressor(steps = n_reps)

        condition_results <- furrr::future_map_dfr(
          1:n_reps,
          function(rep) {

            # Each rep gets its own temp directory to prevent file collisions
            tmp_dir <- file.path(
              tempdir(),
              sprintf("mplus_s%d_N%d_T%d_r%04d_p%d",
                      local_sim, local_N, local_T, rep, Sys.getpid())
            )
            dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
            on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

            dat <- local_gen_fn(N = local_N, T = local_T,
                                seed = rep * 1000 + local_sim)

            result <- tryCatch({

              fit_time <- system.time({
                suppressMessages(
                  mplusModeler(
                    mplusObject(
                      TITLE    = sprintf("s%d_N%d_T%d_r%04d",
                                         local_sim, local_N, local_T, rep),
                      VARIABLE = local_model$VARIABLE,
                      ANALYSIS = local_analysis,
                      MODEL    = local_model$MODEL,
                      OUTPUT   = "CINTERVAL(HPD);",
                      usevariables = c("id", "time", "craving_smoking", "stressed"),
                      rdata    = dat
                    ),
                    dataout      = file.path(tmp_dir, "data.dat"),
                    modelout     = file.path(tmp_dir, "model.inp"),
                    run          = 1L,
                    writeData    = "always",
                    hashfilename = FALSE,
                    quiet        = TRUE
                  )
                )
              })

              out_file <- file.path(tmp_dir, "model.out")
              if (!file.exists(out_file)) stop("No .out file produced")

              res_read  <- readModels(out_file, quiet = TRUE)
              extracted <- .extract_mplus(res_read, local_sim)
              if (is.null(extracted)) stop("Parameter extraction returned NULL")

              extracted %>%
                mutate(converged = TRUE, rep = rep,
                       elapsed   = unname(fit_time["elapsed"]))

            }, error = function(e) {
              tibble(converged = FALSE, rep = rep,
                     elapsed   = NA_real_,
                     error_msg = conditionMessage(e))
            })

            p(sprintf("rep %d | %.0fs", rep,
                      if (isTRUE(result$converged[1])) result$elapsed[1] else 0))
            result
          },
          .options = furrr::furrr_options(seed = TRUE)
        )
      })

      elapsed_min <- unname(proc.time() - start_t)["elapsed"] / 60
      n_conv      <- sum(condition_results$converged, na.rm = TRUE)
      cat(sprintf("  Converged: %d/%d | Elapsed: %.1f min\n",
                  n_conv, n_reps, elapsed_min))

      # Save per-condition RDS
      fname <- sprintf("sim%d_N%d_T%d_mplus.Rds", sim_num, N, T)
      saveRDS(condition_results, file.path(save_dir, fname))

      # Rough ETA for remaining conditions
      remaining <- n_cond - cond_idx
      if (remaining > 0)
        cat(sprintf("  Est. remaining: ~%.1f min\n", elapsed_min * remaining))

      all_results[[sprintf("N%d_T%d", N, T)]] <- condition_results
    }
  }

  full <- dplyr::bind_rows(all_results) %>%
    mutate(N = rep(rep(N_vec, each = length(T_vec) * n_reps),
                   times = 1),
           .before = 1)

  # Recalculate N and T columns properly
  full <- dplyr::bind_rows(
    Map(function(key, df) {
      parts <- strsplit(key, "_")[[1]]
      df %>% mutate(
        N = as.integer(sub("N", "", parts[1])),
        T = as.integer(sub("T", "", parts[2]))
      )
    },
    names(all_results),
    all_results)
  )

  saveRDS(full,
          file.path(save_dir, sprintf("sim%d_full_mplus.Rds", sim_num)))
  cat(sprintf(
    "\nAll conditions done — saved: sim%d_full_mplus.Rds\n", sim_num
  ))
  invisible(full)
}

# =============================================================================
# Example calls (run interactively after setting plan):
#
#   source(here::here("simulation_run_mplus.R"))
#   plan(multisession, workers = 4)
#   # workers × (chains × procs) should not exceed total cores;
#   # e.g., 4 workers × 2 chains × 2 procs = 16 cores
#
#   results_sim1 <- run_simulation_mplus_parallel(
#     sim_num = 1, N_vec = c(50, 100, 200), T_vec = c(30, 50),
#     n_reps = 1000, biter = 50000, warmup = 10000, chains = 2, procs = 2
#   )
#   results_sim2 <- run_simulation_mplus_parallel(sim_num = 2, n_reps = 1000)
#   results_sim3 <- run_simulation_mplus_parallel(sim_num = 3, n_reps = 1000)
#   results_sim4 <- run_simulation_mplus_parallel(sim_num = 4, n_reps = 1000)
# =============================================================================
