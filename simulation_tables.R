# =============================================================================
# simulation_tables.R
# Summarize Blimp DSEM simulation results and produce gt tables
# =============================================================================

library(tidyverse)
library(here)
library(gt)

# ── Truth values ──────────────────────────────────────────────────────────────
truth <- list(
  `1` = list(AR_mean = 0.30, CL_mean = 0.10, AR_var = 0.01, CL_var = 0.02),
  `2` = list(AR_1_mean = 0.30, AR_2_mean = 0.10,
             CL_1_mean = 0.05, CL_2_mean = 0.10,
             AR_1_var  = 0.02, AR_2_var  = 0.01,
             CL_1_var  = 0.002, CL_2_var = 0.001),
  `3` = list(AR_mean = 0.30, CL_mean = 0.10, trend_mean = 0.10,
             AR_var  = 0.01, CL_var  = 0.02,  trend_var  = 0.01),
  `4` = list(AR_1_mean = 0.30, AR_2_mean = 0.10,
             CL_1_mean = 0.05, CL_2_mean = 0.10,
             craving_trend_mean = 0.10, stressed_trend_mean = 0.10,
             AR_1_var  = 0.02, AR_2_var  = 0.01,
             CL_1_var  = 0.002, CL_2_var = 0.001,
             craving_trend_var = 0.01, stressed_trend_var = 0.01)
)

# ── Parameter sets with CI columns ───────────────────────────────────────────
params_with_ci <- list(
  `1` = c("AR_mean", "CL_mean"),
  `2` = c("AR_1_mean", "AR_2_mean", "CL_1_mean", "CL_2_mean"),
  `3` = c("AR_mean", "CL_mean", "trend_mean"),
  `4` = c("AR_1_mean", "AR_2_mean", "CL_1_mean", "CL_2_mean",
          "craving_trend_mean", "stressed_trend_mean")
)

params_no_ci <- list(
  `1` = c("AR_var", "CL_var"),
  `2` = c("AR_1_var", "AR_2_var", "CL_1_var", "CL_2_var"),
  `3` = c("AR_var", "CL_var", "trend_var"),
  `4` = c("AR_1_var", "AR_2_var", "CL_1_var", "CL_2_var",
          "craving_trend_var", "stressed_trend_var")
)

# ── Core summarization function ───────────────────────────────────────────────
summarize_sim <- function(results, sim_num) {
  tr   <- truth[[as.character(sim_num)]]
  pci  <- params_with_ci[[as.character(sim_num)]]
  pnci <- params_no_ci[[as.character(sim_num)]]

  results %>%
    filter(converged) %>%
    group_by(N, T) %>%
    summarise(
      n_converged = n(),
      # Parameters with CIs: bias, rel_bias, rmse, coverage
      across(
        all_of(pci),
        list(
          bias     = ~ mean(.x, na.rm = TRUE) - tr[[cur_column()]],
          rel_bias = ~ (mean(.x, na.rm = TRUE) - tr[[cur_column()]]) /
                       abs(tr[[cur_column()]]) * 100,
          rmse     = ~ sqrt(mean((.x - tr[[cur_column()]])^2, na.rm = TRUE)),
          coverage = ~ {
            lo <- cur_data()[[paste0(cur_column(), "_lo")]]
            hi <- cur_data()[[paste0(cur_column(), "_hi")]]
            mean(lo <= tr[[cur_column()]] & hi >= tr[[cur_column()]], na.rm = TRUE) * 100
          }
        ),
        .names = "{.col}__{.fn}"
      ),
      # Parameters without CIs: bias, rel_bias, rmse only
      across(
        all_of(pnci),
        list(
          bias     = ~ mean(.x, na.rm = TRUE) - tr[[cur_column()]],
          rel_bias = ~ (mean(.x, na.rm = TRUE) - tr[[cur_column()]]) /
                       abs(tr[[cur_column()]]) * 100,
          rmse     = ~ sqrt(mean((.x - tr[[cur_column()]])^2, na.rm = TRUE))
        ),
        .names = "{.col}__{.fn}"
      ),
      .groups = "drop"
    ) %>%
    arrange(N, T)
}

# ── Load full results ─────────────────────────────────────────────────────────
sim1 <- readRDS(here("simulation_results/sim1_full_blimp.Rds"))
sim2 <- readRDS(here("simulation_results/sim2_full_blimp.Rds"))
sim3 <- readRDS(here("simulation_results/sim3_full_blimp.Rds"))
sim4 <- readRDS(here("simulation_results/sim4_full_blimp.Rds"))

summary1 <- summarize_sim(sim1, 1)
summary2 <- summarize_sim(sim2, 2)
summary3 <- summarize_sim(sim3, 3)
summary4 <- summarize_sim(sim4, 4)

# ── gt table builder ──────────────────────────────────────────────────────────
make_gt_table <- function(summary, sim_num) {

  sim_labels <- c(
    "1" = "Sim 1: AR(1) DSEM",
    "2" = "Sim 2: VAR(1) DSEM",
    "3" = "Sim 3: AR(1) RDSEM",
    "4" = "Sim 4: VAR(1) RDSEM"
  )

  # Pivot to long: one row per condition × parameter × metric
  long <- summary %>%
    pivot_longer(
      cols      = -c(N, T, n_converged),
      names_to  = c("parameter", "metric"),
      names_sep = "__"
    ) %>%
    pivot_wider(names_from = metric, values_from = value) %>%
    mutate(
      condition = sprintf("N=%d, T=%d", N, T),
      rel_bias  = round(rel_bias, 1),
      bias      = round(bias, 4),
      rmse      = round(rmse, 4),
      coverage  = if ("coverage" %in% names(.)) round(coverage, 1) else NA_real_
    ) %>%
    select(condition, parameter, bias, rel_bias, rmse, any_of("coverage"))

  # Build gt
  tbl <- long %>%
    gt(groupname_col = "condition") %>%
    tab_header(
      title    = sim_labels[as.character(sim_num)],
      subtitle = "Bias, Relative Bias (%), RMSE, and Coverage (%) across conditions"
    ) %>%
    cols_label(
      parameter = "Parameter",
      bias      = "Bias",
      rel_bias  = "Rel. Bias (%)",
      rmse      = "RMSE",
      coverage  = "Coverage (%)"
    ) %>%
    fmt_number(columns = c(bias, rmse), decimals = 4) %>%
    fmt_number(columns = c(rel_bias), decimals = 1) %>%
    sub_missing(columns = coverage, missing_text = "—") %>%
    # Highlight coverage outside 93–97%
    tab_style(
      style     = cell_fill(color = "#fff3cd"),
      locations = cells_body(
        columns = coverage,
        rows    = !is.na(coverage) & (coverage < 93 | coverage > 97)
      )
    ) %>%
    # Highlight large relative bias
    tab_style(
      style     = cell_fill(color = "#f8d7da"),
      locations = cells_body(
        columns = rel_bias,
        rows    = abs(rel_bias) > 10
      )
    ) %>%
    tab_footnote("Yellow = coverage outside 93–97%. Red = |relative bias| > 10%.") %>%
    opt_stylize(style = 1) %>%
    tab_options(
      row_group.font.weight  = "bold",
      heading.title.font.size = px(16)
    )

  tbl
}

# ── Print tables ──────────────────────────────────────────────────────────────
gt1 <- make_gt_table(summary1, 1)
gt2 <- make_gt_table(summary2, 2)
gt3 <- make_gt_table(summary3, 3)
gt4 <- make_gt_table(summary4, 4)

# View any table:
gt1
gt2
gt3
gt4
