devtools::load_all("~/Projects/slendr")
init_env(quiet = TRUE)

library(ggplot2)
library(ggpubr)
library(dplyr)
library(broom)
library(tidyr)
library(parallel)

NE_START <- c(10e3, 20e3, 30e3, 40e3)
NE_HUNTED <- c(100, 250, 500, 750, 1000, 1500, 2000, 3000)
CENSUS_RATIO <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)

GENERATION_TIME <- 32
SEQUENCE_LENGTH <- 10e6
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

simulate_ts <- function(model) {
  # schedule sampling at regular time points
  population <- model$populations[[1]]
  samples <- schedule_sampling(model, times = c(seq(T_START, T_END, by = 50)), list(population, 100))

  # simulate tree sequence
  ts <-
    msprime(model, sequence_length = SEQUENCE_LENGTH, recombination_rate = RECOMBINATION_RATE, samples = samples) %>%
    ts_mutate(mutation_rate = MUTATION_RATE)

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

plot_diversity <- function(results, Ne_start, Ne_hunted, census_ratio) {
  df <- filter(results, Ne_start == !!Ne_start, Ne_hunted == !!Ne_hunted, census_ratio == !!census_ratio)
  if (nrow(df) == 0) stop("No simulation corresponding to the given parameter combination", call. = FALSE)
  if (is.null(df$model[[1]])) stop("This parameter combination lead to an invalid model", call. = FALSE)

  ggplot(df$pi[[1]], aes(time, pi, group = time)) +
    geom_jitter(size = 0.2, width = 5) +
    geom_violin(trim = FALSE, alpha = 0.5) +
    geom_smooth(aes(group = 1), method = lm, formula = y ~ x) +
    # stat_cor(aes(group = 1, label = paste(after_stat(rr.label), after_stat(p.label), sep = "~~~~"))) +
    labs(
      x = "time [years C.E.]",
      y = "nucleotide diversity",
      title = sprintf("starting Ne = %s, Ne = %s hunted per generation, Ne / census ratio = %s",
                      Ne_start, Ne_hunted, census_ratio),
      subtitle = paste0("linear regression p-value = ", format(df$lm[[1]][, "p.value"], digits = 2),
                       ", R-squared = ", format(df$lm[[1]][, "r.squared"], digits = 2),
                       ", slope = ", format(df$lm[[1]][, "slope"], digits = 2))
    )
}

plot_demography <- function(results, Ne_start, Ne_hunted, census_ratio) {
  df <- filter(results, Ne_start == !!Ne_start, Ne_hunted == !!Ne_hunted, census_ratio == !!census_ratio)
  if (nrow(df) == 0) stop("No simulation corresponding to the given parameter combination", call. = FALSE)
  if (is.null(df$model[[1]])) stop("This parameter combination lead to an invalid model", call. = FALSE)

  plot_model(df$model[[1]])
}

plot_slopes <- function(results) {
  results_df %>%
    unnest(lm) %>%
    mutate(Ne_start = paste("Ne start =", Ne_start)) %>%
    ggplot(aes(Ne_hunted, slope, color = census_ratio)) +
    geom_point() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    facet_grid(. ~ Ne_start) +
    guides(color = guide_legend("Ne / census\nsize ratio")) +
    labs(x = "Ne hunted", y = "slope of the nucleotide diversity trajectory")
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

results <- mclapply((1:nrow(grid)), function(i) {
  # cat(sprintf("%d / %d", i, nrow(grid)))
  Ne_start <- grid[i, ]$NE_START
  Ne_hunted <- grid[i, ]$NE_HUNTED
  census_ratio <- grid[i, ]$CENSUS_RATIO

  model <- generate_model(Ne_start, Ne_hunted, census_ratio)

  if (is.null(model))
    ts <- NULL
  else
    ts <- simulate_ts(model)

  results <- compute_results(ts, Ne_start = Ne_start, Ne_hunted = Ne_hunted, census_ratio = census_ratio)
  cat("\r")

  results
}, mc.cores = detectCores())

results_df <- bind_rows(results)

unique(results_df$Ne_start)
unique(results_df$Ne_hunted)
unique(results_df$census_ratio)

# > unique(results_df$Ne_start)
# [1] 10000 20000 30000 40000
# > unique(results_df$Ne_hunted)
# [1]  100  250  500  750 1000 1500 2000 3000
# > unique(results_df$census_ratio)
# [1] 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0

plot_slopes(results_df)

plot_demography(results_df, Ne_start = 40000, Ne_hunted = 1000, census_ratio = 0.1)
plot_diversity(results_df, Ne_start = 40000, Ne_hunted = 1000, census_ratio = 0.1)

ggplot(df, aes(time, pi, group = time)) +
  geom_jitter(size = 0.2, width = 5) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_smooth(aes(group = 1), method = "lm") +
  stat_cor(aes(group = 1, label = paste(after_stat(rr.label), after_stat(p.label), sep = "~~~~"))) +
  labs(
    x = "time [years C.E.]",
    y = "nucleotide diversity",
    title = sprintf("starting Ne = %s, %s hunted per generation)",
                    Ne_start, N_hunted),
    subtitle = paste0("linear regression p-value = ", format(df_lm[, "p.value"], digits = 2),
                     ", R-squared = ", format(df_lm[, "r.squared"], digits = 2),
                     ", slope = ", format(df_lm[, "slope"], digits = 2))
  )
