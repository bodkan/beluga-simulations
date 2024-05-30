# beluga diversity ----------------------------------------------------------------------------

library(slendr)
init_env(quiet = TRUE)

library(ggplot2)
library(ggpubr)
library(dplyr)
library(broom)
library(tidyr)
library(parallel)
library(cowplot)

NE_START <- c(40e3, 30e3, 20e3, 10e3)
NE_HUNTED <- c(100, 250, 500, 750, 1000, 1500, 2000, 3000)
CENSUS_RATIO <- c(1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1)

GENERATION_TIME <- 32
SEQUENCE_LENGTH <- 100e6
RECOMBINATION_RATE <- 1e-8
MUTATION_RATE <- 1e-8

T_START <- 1150
T_END <- 2200
T_STEP <- 32

generate_model <- function(Ne_start, Ne_hunted, census_ratio) {
  # create the initial population of belugas
  t <- T_START
  N <- Ne_start * census_ratio
  demography <- population("Beluga", time = t, N = N)

  # every t_step years, decrease the population by N_hunted
  while (TRUE) {
    if (t > T_END) break

    t <- t + T_STEP
    N <- N - Ne_hunted * census_ratio
    # cat(t, "\t", N, "\n")
    if (N <= 0)
      return(NULL)

    demography <- resize(demography, time = t, N = N, how = "step")
  }

  model <- compile_model(demography, generation_time = GENERATION_TIME, simulation_length = 1000, time_units = "years C.E.")

  model
}

simulate_ts <- function(model, engine_fun, Ne_start) {
  # schedule sampling at regular time points
  population <- model$populations[[1]]
  samples <- schedule_sampling(model, times = c(seq(T_START, T_END, by = 50)), list(population, 100))

  # simulate tree sequence
  ts <-
    engine_fun(model, sequence_length = SEQUENCE_LENGTH, recombination_rate = RECOMBINATION_RATE, samples = samples, random_seed = 42)

  if (identical(engine_fun, slim)) {
    ts <- ts_recapitate(ts, Ne = Ne_start, recombination_rate = RECOMBINATION_RATE, random_seed = 42) %>%
      ts_simplify()
  }

  ts <- ts_mutate(ts, mutation_rate = MUTATION_RATE, random_seed = 42)

  ts
}

compute_results <- function(ts, Ne_start, Ne_hunted, census_ratio) {
  result <- tibble(Ne_start, Ne_hunted, census_ratio)

  if (is.null(ts)) {
    df <- NULL
    df_lm <- NULL
    model <- NULL
  } else {
    model <- attr(ts, "model")

    df <- ts_samples(ts)
    df$pi <- ts_diversity(ts, df$name)$diversity

    res <- lm(pi ~ time, data = df)
    summary(res)

    df_lm <- cbind(
      tidy(res) %>% filter(term == "time") %>% select(slope = estimate),
      glance(res) %>% select(r.squared, p.value)
    )
  }

  result <- result %>% mutate(
    model = list(model),
    pi = list(df),
    lm = list(df_lm)
  )

  result
}

# test run ------------------------------------------------------------------------------------

# Ne_start <- NE_START[1]
# Ne_hunted <- NE_HUNTED[1]
# census_ratio <- CENSUS_RATIO[1]
#
# model <- generate_model(Ne_start, Ne_hunted, census_ratio)
#
# plot_model(model)
#
# ts <- simulate_ts(model)
#
# results <- compute_results(ts, Ne_start = Ne_start, Ne_hunted = Ne_hunted, census_ratio = census_ratio)
#
# plot_diversity(results, Ne_start = Ne_start, Ne_hunted = Ne_hunted, census_ratio = census_ratio)
# plot_slopes(results)

# grid run ------------------------------------------------------------------------------------


grid <- expand_grid(NE_START, NE_HUNTED, CENSUS_RATIO)

t_start <- Sys.time()

results <- mclapply((1:nrow(grid)), function(i) {
  Ne_start <- grid[i, ]$NE_START
  Ne_hunted <- grid[i, ]$NE_HUNTED
  census_ratio <- grid[i, ]$CENSUS_RATIO

  model <- generate_model(Ne_start, Ne_hunted, census_ratio)

  if (is.null(model)) {
    ts_slim <- NULL
    ts_msprime <- NULL
  } else {
    ts_slim <- simulate_ts(model, engine_fun = slim, Ne_start = Ne_start)
    ts_msprime <- simulate_ts(model, engine_fun = msprime, Ne_start = Ne_start)
  }

  results_slim <- compute_results(ts_slim, Ne_start = Ne_start, Ne_hunted = Ne_hunted, census_ratio = census_ratio)
  results_msprime <- compute_results(ts_msprime, Ne_start = Ne_start, Ne_hunted = Ne_hunted, census_ratio = census_ratio)
  results_slim$engine <- "SLiM"
  results_msprime$engine <- "msprime"

  results <- rbind(results_slim, results_msprime)
  results
}, mc.cores = detectCores() / 2)

t_end <- Sys.time()

t_end - t_start
# Time difference of 1.650988 hours

results_df <- bind_rows(results)

saveRDS(results_df, sprintf("pi_beluga_%sMb.rds", SEQUENCE_LENGTH / 1e6))
