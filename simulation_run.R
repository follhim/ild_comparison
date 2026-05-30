# =============================================================================
# simulation_run.R
# Blimp DSEM simulation — with progressr ETA bars
# =============================================================================

library(tidyverse)
library(rblimp)
library(furrr)
library(future)
library(here)
library(progressr)   # progress with ETA
library(progress)    # needed by handler_progress

plan(multisession, workers = 3)   # 3 workers × 4 Blimp chains = 12 cores

# ── Progress handler: shows bar + % + ETA + current message ──────────────────
handlers(handler_progress(
  format = "[:bar] :percent | :current/:total reps | ETA :eta | Elapsed :elapsed | :message",
  width  = 80,
  clear  = FALSE   # keep finished bars visible so you can read timing
))

# ── Core simulation function ──────────────────────────────────────────────────
run_simulation_blimp_parallel <- function(
    sim_num,
    N_vec    = c(50, 100, 200),
    T_vec    = c(30, 50),
    n_reps   = 500,
    burn     = 2000,
    iter     = 5000,
    save_dir = here("simulation_results")
) {

  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

  gen_fn <- switch(sim_num,
    `1` = generate_ar1_dsem,
    `2` = generate_var1_dsem,
    `3` = generate_ar1_rdsem,
    `4` = generate_var1_rdsem
  )

  model_syntax <- switch(sim_num,
    `1` = list(
      latent = 'id = stressed_mean craving_mean AR CL',
      model  = '
        level2.models:
          1 -> stressed_mean craving_mean AR CL;
        level1.models:
          stressed ~ 1@stressed_mean;
          craving_smoking ~
            1@craving_mean
            (craving_smoking.lag - craving_mean)@AR
            (stressed - stressed_mean)@CL;'
    ),
    `2` = list(
      latent = 'id = stressed_mean craving_mean AR_1 CL_1 AR_2 CL_2',
      model  = '
        level2.models:
          1 -> stressed_mean craving_mean AR_1 CL_1 AR_2 CL_2;
        level1.models:
          stressed ~
            1@stressed_mean
            (stressed.lag - stressed_mean)@AR_1
            (craving_smoking.lag - craving_mean)@CL_1;
          craving_smoking ~
            1@craving_mean
            (craving_smoking.lag - craving_mean)@AR_2
            (stressed.lag - stressed_mean)@CL_2;'
    ),
    `3` = list(
      latent = 'id = stressed_mean craving_mean craving_trend AR CL',
      model  = '
        lag_hat = craving_mean + (time - 1)*craving_trend;
        lag_res = ifelse(time <= 1, 0, craving_smoking.lag - lag_hat);
        level2.models:
          1 -> stressed_mean craving_mean AR CL craving_trend;
        level1.models:
          stressed ~ 1@stressed_mean;
          craving_smoking ~
            1@craving_mean
            lag_res@AR
            (stressed - stressed_mean)@CL
            time@craving_trend;'
    ),
    `4` = list(
      latent = 'id = craving_mean stressed_mean
                CL_1 AR_1 CL_2 AR_2
                craving_trend stressed_trend',
      model  = '
        craving_lag_hat = craving_mean + (time - 1)*craving_trend;
        craving_lag_res = ifelse(time <= 1, 0,
                            craving_smoking.lag - craving_lag_hat);
        stressed_lag_hat = stressed_mean + (time - 1)*stressed_trend;
        stressed_lag_res = ifelse(time <= 1, 0,
                            stressed.lag - stressed_lag_hat);
        level2.models:
          1 -> craving_mean stressed_mean
               CL_1 AR_1 CL_2 AR_2
               craving_trend stressed_trend;
        level1.models:
          stressed ~
            1@stressed_mean
            craving_lag_res@CL_1
            stressed_lag_res@AR_1
            time@stressed_trend;
          craving_smoking ~
            1@craving_mean
            craving_lag_res@AR_2
            stressed_lag_res@CL_2
            time@craving_trend;'
    )
  )

  extract_inline <- function(fit, sim_num) {
    est <- fit@estimates
    get_param <- function(pattern) {
      rows <- est[grep(pattern, rownames(est), ignore.case = TRUE), , drop = FALSE]
      if (nrow(rows) == 0) return(c(est = NA, lo = NA, hi = NA))
      c(est = rows[1, "Estimate"], lo = rows[1, "2.5%"], hi = rows[1, "97.5%"])
    }
    if (sim_num == 1) {
      tibble::tibble(
        AR_mean    = get_param("^AR ~ Intercept")[["est"]],
        AR_mean_lo = get_param("^AR ~ Intercept")[["lo"]],
        AR_mean_hi = get_param("^AR ~ Intercept")[["hi"]],
        CL_mean    = get_param("^CL ~ Intercept")[["est"]],
        CL_mean_lo = get_param("^CL ~ Intercept")[["lo"]],
        CL_mean_hi = get_param("^CL ~ Intercept")[["hi"]],
        AR_var     = get_param("^AR residual variance")[["est"]],
        CL_var     = get_param("^CL residual variance")[["est"]]
      )
    } else if (sim_num == 2) {
      tibble::tibble(
        AR_1_mean    = get_param("^AR_1 ~ Intercept")[["est"]],
        AR_1_mean_lo = get_param("^AR_1 ~ Intercept")[["lo"]],
        AR_1_mean_hi = get_param("^AR_1 ~ Intercept")[["hi"]],
        AR_2_mean    = get_param("^AR_2 ~ Intercept")[["est"]],
        AR_2_mean_lo = get_param("^AR_2 ~ Intercept")[["lo"]],
        AR_2_mean_hi = get_param("^AR_2 ~ Intercept")[["hi"]],
        CL_1_mean    = get_param("^CL_1 ~ Intercept")[["est"]],
        CL_1_mean_lo = get_param("^CL_1 ~ Intercept")[["lo"]],
        CL_1_mean_hi = get_param("^CL_1 ~ Intercept")[["hi"]],
        CL_2_mean    = get_param("^CL_2 ~ Intercept")[["est"]],
        CL_2_mean_lo = get_param("^CL_2 ~ Intercept")[["lo"]],
        CL_2_mean_hi = get_param("^CL_2 ~ Intercept")[["hi"]],
        AR_1_var  = get_param("^AR_1 residual variance")[["est"]],
        AR_2_var  = get_param("^AR_2 residual variance")[["est"]],
        CL_1_var  = get_param("^CL_1 residual variance")[["est"]],
        CL_2_var  = get_param("^CL_2 residual variance")[["est"]]
      )
    } else if (sim_num == 3) {
      tibble::tibble(
        AR_mean       = get_param("^AR ~ Intercept")[["est"]],
        AR_mean_lo    = get_param("^AR ~ Intercept")[["lo"]],
        AR_mean_hi    = get_param("^AR ~ Intercept")[["hi"]],
        CL_mean       = get_param("^CL ~ Intercept")[["est"]],
        CL_mean_lo    = get_param("^CL ~ Intercept")[["lo"]],
        CL_mean_hi    = get_param("^CL ~ Intercept")[["hi"]],
        trend_mean    = get_param("^craving_trend ~ Intercept")[["est"]],
        trend_mean_lo = get_param("^craving_trend ~ Intercept")[["lo"]],
        trend_mean_hi = get_param("^craving_trend ~ Intercept")[["hi"]],
        AR_var    = get_param("^AR residual variance")[["est"]],
        CL_var    = get_param("^CL residual variance")[["est"]],
        trend_var = get_param("^craving_trend residual variance")[["est"]]
      )
    } else if (sim_num == 4) {
      tibble::tibble(
        AR_1_mean              = get_param("^AR_1 ~ Intercept")[["est"]],
        AR_1_mean_lo           = get_param("^AR_1 ~ Intercept")[["lo"]],
        AR_1_mean_hi           = get_param("^AR_1 ~ Intercept")[["hi"]],
        AR_2_mean              = get_param("^AR_2 ~ Intercept")[["est"]],
        AR_2_mean_lo           = get_param("^AR_2 ~ Intercept")[["lo"]],
        AR_2_mean_hi           = get_param("^AR_2 ~ Intercept")[["hi"]],
        CL_1_mean              = get_param("^CL_1 ~ Intercept")[["est"]],
        CL_1_mean_lo           = get_param("^CL_1 ~ Intercept")[["lo"]],
        CL_1_mean_hi           = get_param("^CL_1 ~ Intercept")[["hi"]],
        CL_2_mean              = get_param("^CL_2 ~ Intercept")[["est"]],
        CL_2_mean_lo           = get_param("^CL_2 ~ Intercept")[["lo"]],
        CL_2_mean_hi           = get_param("^CL_2 ~ Intercept")[["hi"]],
        craving_trend_mean     = get_param("^craving_trend ~ Intercept")[["est"]],
        craving_trend_mean_lo  = get_param("^craving_trend ~ Intercept")[["lo"]],
        craving_trend_mean_hi  = get_param("^craving_trend ~ Intercept")[["hi"]],
        stressed_trend_mean    = get_param("^stressed_trend ~ Intercept")[["est"]],
        stressed_trend_mean_lo = get_param("^stressed_trend ~ Intercept")[["lo"]],
        stressed_trend_mean_hi = get_param("^stressed_trend ~ Intercept")[["hi"]],
        AR_1_var           = get_param("^AR_1 residual variance")[["est"]],
        AR_2_var           = get_param("^AR_2 residual variance")[["est"]],
        CL_1_var           = get_param("^CL_1 residual variance")[["est"]],
        CL_2_var           = get_param("^CL_2 residual variance")[["est"]],
        craving_trend_var  = get_param("^craving_trend residual variance")[["est"]],
        stressed_trend_var = get_param("^stressed_trend residual variance")[["est"]]
      )
    }
  }

  # Localize everything for furrr workers
  local_gen_fn  <- gen_fn
  local_extract <- extract_inline
  local_syntax  <- model_syntax
  local_burn    <- burn
  local_iter    <- iter
  local_sim     <- sim_num

  all_results  <- list()
  sim_start    <- proc.time()

  n_conditions <- length(N_vec) * length(T_vec)
  condition_i  <- 0L

  for (N in N_vec) {
    for (T in T_vec) {

      condition_i <- condition_i + 1L
      cond_label  <- sprintf("Sim %d | N=%d T=%d [%d/%d]",
                             sim_num, N, T, condition_i, n_conditions)

      cat("\n========================================\n")
      cat(cond_label, "\n")
      cat(sprintf("Running %d replications...\n", n_reps))
      cat("========================================\n")

      # ── progressr: one bar per condition, step per completed rep ────────────
      with_progress({
        p <- progressor(steps = n_reps)

        condition_results <- furrr::future_map_dfr(
          1:n_reps,
          function(rep) {

            dat <- local_gen_fn(N = N, T = T, seed = rep * 1000 + local_sim)

            fit_time <- system.time({
              fit <- tryCatch({
                rblimp::rblimp(
                  data         = dat,
                  clusterid    = 'id',
                  timeid       = 'time',
                  latent       = local_syntax$latent,
                  model        = local_syntax$model,
                  burn         = local_burn,
                  iter         = local_iter,
                  seed         = rep,
                  print_output = FALSE
                )
              }, error = function(e) {
                message(sprintf("Rep %d fit failed: %s", rep, e$message))
                return(NULL)
              })
            })

            # Signal progress — this travels back from worker to main process
            p(sprintf("rep %d | %.0fs", rep, fit_time["elapsed"]))

            if (is.null(fit)) {
              return(tibble::tibble(rep = rep, converged = FALSE, time_sec = NA_real_))
            }

            estimates <- tryCatch(
              local_extract(fit, local_sim),
              error = function(e) {
                message(sprintf("Rep %d extraction failed: %s", rep, e$message))
                return(NULL)
              }
            )

            if (is.null(estimates)) {
              return(tibble::tibble(
                rep = rep, converged = FALSE,
                time_sec = as.numeric(fit_time["elapsed"])
              ))
            }

            estimates %>%
              dplyr::mutate(
                rep       = rep,
                converged = TRUE,
                time_sec  = as.numeric(fit_time["elapsed"])
              )
          },
          .options = furrr::furrr_options(seed = TRUE)
          # Note: NO .progress = TRUE here — progressr handles it
        ) %>%
          dplyr::mutate(sim = sim_num, N = N, T = T)
      })  # end with_progress

      all_results[[paste0("N", N, "_T", T)]] <- condition_results

      saveRDS(
        condition_results,
        here::here(save_dir,
             sprintf("sim%d_N%d_T%d_blimp.Rds", sim_num, N, T))
      )

      n_conv     <- sum(condition_results$converged, na.rm = TRUE)
      mean_t     <- mean(condition_results$time_sec, na.rm = TRUE)
      min_t      <- min(condition_results$time_sec,  na.rm = TRUE)
      max_t      <- max(condition_results$time_sec,  na.rm = TRUE)

      cat(sprintf("\n  Saved: sim%d_N%d_T%d_blimp.Rds\n", sim_num, N, T))
      cat(sprintf("  Converged: %d / %d\n", n_conv, n_reps))
      cat(sprintf("  Rep timing — Mean: %.1f sec | Min: %.1f | Max: %.1f\n",
                  mean_t, min_t, max_t))

      # Rough projection for remaining conditions
      elapsed_so_far <- (proc.time() - sim_start)["elapsed"]
      secs_per_cond  <- elapsed_so_far / condition_i
      remaining_cond <- n_conditions - condition_i
      eta_min        <- (secs_per_cond * remaining_cond) / 60
      if (remaining_cond > 0) {
        cat(sprintf("  ETA for remaining %d condition(s): ~%.0f min\n",
                    remaining_cond, eta_min))
      }
    }
  }

  full_results <- dplyr::bind_rows(all_results)

  saveRDS(
    full_results,
    here::here(save_dir, sprintf("sim%d_full_blimp.Rds", sim_num))
  )

  total_elapsed <- (proc.time() - sim_start)["elapsed"]
  n_total_conv  <- sum(full_results$converged, na.rm = TRUE)
  n_total_reps  <- nrow(full_results)

  cat(sprintf("\n========================================\n"))
  cat(sprintf("Simulation %d complete in %.1f min.\n", sim_num, total_elapsed / 60))
  cat(sprintf("Total converged: %d / %d\n", n_total_conv, n_total_reps))
  cat(sprintf("Overall mean time per rep: %.1f sec\n",
              mean(full_results$time_sec, na.rm = TRUE)))
  cat(sprintf("========================================\n"))

  return(full_results)
}


# ── Example calls ─────────────────────────────────────────────────────────────
# Run sims individually — each saves its own .Rds files independently.

# results_sim1 <- run_simulation_blimp_parallel(sim_num = 1, N_vec = c(50, 100, 200), T_vec = c(30, 50), n_reps = 500, burn = 2000, iter = 5000)
# results_sim2 <- run_simulation_blimp_parallel(sim_num = 2, N_vec = c(50, 100, 200), T_vec = c(30, 50), n_reps = 500, burn = 2000, iter = 5000)
# results_sim3 <- run_simulation_blimp_parallel(sim_num = 3, N_vec = c(50, 100, 200), T_vec = c(30, 50), n_reps = 500, burn = 2000, iter = 5000)
# results_sim4 <- run_simulation_blimp_parallel(sim_num = 4, N_vec = c(50, 100, 200), T_vec = c(30, 50), n_reps = 500, burn = 2000, iter = 5000)

# ── Quick test (3 reps, low iter) ─────────────────────────────────────────────
# results_sim1 <- run_simulation_blimp_parallel(sim_num = 1, N_vec = c(50), T_vec = c(30), n_reps = 3, burn = 200,
#   iter     = 100,
#   save_dir = here("simulation_results/test"),
#   email    = FALSE
# )
