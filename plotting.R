plot_diversity <- function(results, N_start, N_hunted, census_ratio, engine) {
  df <- filter(results, N_start == !!N_start, N_hunted == !!N_hunted, census_ratio == !!census_ratio, engine == !!engine)
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
      title = sprintf("start Ne = %s, Ne = %s hunted, Ne / census ratio = %s",
                      N_start, N_hunted, census_ratio),
      subtitle = paste0("regression p-value = ", format(df$lm[[1]][, "p.value"], digits = 2),
                        ", R-squared = ", format(df$lm[[1]][, "r.squared"], digits = 2),
                        ", slope = ", format(df$lm[[1]][, "slope"], digits = 2))
    ) +
    expand_limits(y = 0)
}

plot_demography <- function(results, N_start, N_hunted, census_ratio, engine) {
  df <- filter(results, N_start == !!N_start, N_hunted == !!N_hunted, census_ratio == !!census_ratio, engine == !!engine)
  model <- df$model[[1]]

  if (nrow(df) == 0) stop("No simulation corresponding to the given parameter combination", call. = FALSE)
  if (is.null(model[[1]])) stop("This parameter combination lead to an invalid model", call. = FALSE)

  start_N <- model %>% extract_parameters() %>% .$splits %>% .$N
  end_N <- model %>% extract_parameters() %>% .$resizes %>% .[nrow(.), ] %>% .$N

  plot_model(model) + ggtitle(paste0("start Ne = ", start_N, "\nend Ne = ", end_N))
}

# plot_slopes <- function(df, engine) {
#   df %>%
#     filter(engine == !!engine) %>% unnest(pi) %>%
#     group_by(N_start, N_hunted, census_ratio, time) %>%
#     summarise(pi = mean(pi)) %>%
#     do({
#       model <- lm(pi ~ time, data = .)
#       tidy_result <- tidy(model)
#       glance_result <- glance(model)
#       tidy_result %>%
#         filter(term == "time") %>%
#         mutate(r.squared = glance_result$r.squared,
#                p.value = glance_result$p.value)
#     }) %>%
#     mutate(slope = estimate, N_start = paste("Ne start =", N_start)) %>%
#     ggplot(aes(N_hunted, slope, color = census_ratio)) +
#     geom_point() +
#     geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
#     facet_grid(. ~ N_start) +
#     guides(color = guide_legend("Ne / census\nsize ratio")) +
#     labs(x = "Ne hunted", y = "slope of linear model",
#          title = paste("Slope of the nucleotide diversity trajectory", paste0("(", engine, " engine)")))
# }

plot_slopes <- function(df, engine) {
  df %>%
    filter(engine == !!engine) %>%
    unnest(lm) %>%
    mutate(N_start = paste("Ne start =", N_start)) %>%
    ggplot(aes(N_hunted, slope, color = census_ratio)) +
    geom_point() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    facet_grid(. ~ N_start) +
    guides(color = guide_legend("Ne / census\nsize ratio")) +
    labs(x = "Ne hunted", y = "slope of linear model",
         title = paste("Slope of the nucleotide diversity trajectory", paste0("(", engine, " engine)")))
}

plot_panels <- function(df, N_start, N_hunted, census_ratio, engine) {
  plot_grid(
    plot_slopes(pi_beluga, engine = engine) + coord_cartesian(ylim = c(-1e-8, 1e-8)),
    plot_grid(
      plot_demography(pi_beluga, N_start = N_start, N_hunted = N_hunted, census_ratio = census_ratio, engine),
      plot_diversity(pi_beluga, N_start = N_start, N_hunted = N_hunted, census_ratio = census_ratio, engine),
      nrow = 1, rel_widths = c(0.4, 1)
    ),
    nrow = 2, rel_heights = c(0.5, 1)
  )
}