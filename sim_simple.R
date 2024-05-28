# simple Ne vs diversity analysis -------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(parallel)

library(slendr)
init_env(quiet = TRUE)

pi_df <- mclapply(seq(1000, 50000, by = 1000), function(Ne) {
  pop <- population("pop", time = 10000, N = Ne) %>%
    resize(time = 5000, N = 1000, how = "step")
  model <- compile_model(pop, generation_time = 1, serialize = FALSE)
  samples <- schedule_sampling(model, times = c(5000, 4000, 3000, 2000, 1000, 0), list(pop, 50))
  ts <- msprime(model, sequence_length = 10e6, recombination_rate = 1e-8, samples = samples) %>%
    ts_mutate(1e-8)
  result <- ts_samples(ts) %>% mutate(Ne = Ne)
  result$diversity <- ts_diversity(ts, sample_sets = ts_names(ts))$diversity
  result
}, mc.cores = detectCores()) %>%
  do.call(rbind, .)

pi_df <-
  pi_df %>%
  mutate(bottleneck_time = 5000 - time,
         model = ifelse(bottleneck_time == 0,
                        "constant Ne",
                        paste0("bottleneck at ", bottleneck_time, " gens ago")),
         model = factor(model, levels = unique(sort(model, decreasing = TRUE))))

saveRDS(pi_df, "pi_simple.rds")

pi_df %>% filter(model == "constant Ne") %>%
  ggplot(aes(factor(Ne), diversity, color = model)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +
  # geom_jitter(position = position_dodge(width = 0.8), alpha = 0.2, size = 0.75) +
  geom_smooth(aes(group = model), method = "lm", linewidth = 0.5, color = "black", linetype = 2) +
  scale_color_discrete(drop = FALSE) +
  theme(axis.text.x = element_text(hjust = 1, angle = 45), legend.position = "bottom") +
  labs(x = "Ne", y = "nucleotide diversity",
       title = "Expected nucleotide diversity as a function of Ne")

pi_df %>% filter(model == "constant Ne" | grepl("1000 gens", model)) %>%
  ggplot(aes(factor(Ne), diversity, color = model)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +
  # geom_jitter(position = position_dodge(width = 0.8), alpha = 0.2, size = 0.75) +
  geom_smooth(aes(group = model), method = "lm", linewidth = 0.5, color = "black", linetype = 2) +
  scale_color_discrete(drop = FALSE) +
  theme(axis.text.x = element_text(hjust = 1, angle = 45), legend.position = "bottom") +
  labs(x = "Ne", y = "nucleotide diversity",
       title = "Expected nucleotide diversity as a function of Ne")

pi_df %>% filter(model == "constant Ne" | grepl("1000 gens", model) | grepl("2000 gens", model)) %>%
  ggplot(aes(factor(Ne), diversity, color = model)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +
  # geom_jitter(position = position_dodge(width = 0.8), alpha = 0.2, size = 0.75) +
  geom_smooth(aes(group = model), method = "lm", linewidth = 0.5, color = "black", linetype = 2) +
  scale_color_discrete(drop = FALSE) +
  theme(axis.text.x = element_text(hjust = 1, angle = 45), legend.position = "bottom") +
  labs(x = "Ne", y = "nucleotide diversity",
       title = "Expected nucleotide diversity as a function of Ne")

ggplot(pi_df, aes(factor(Ne), diversity, color = model)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +
  # geom_jitter(position = position_dodge(width = 0.8), alpha = 0.2, size = 0.75) +
  geom_smooth(aes(group = model), method = "lm", linewidth = 0.5, color = "black", linetype = 2) +
  scale_color_discrete(drop = FALSE) +
  theme(axis.text.x = element_text(hjust = 1, angle = 45), legend.position = "bottom") +
  labs(x = "Ne", y = "nucleotide diversity",
       title = "Expected nucleotide diversity as a function of Ne")
