Beluga nucleotide diversity
================
2024-05-30

## Note

Most of the analyses in this document are exploratory. **The code behind
the final results for the paper can be found below under [this
heading](#results-for-the-paper). The rendered result files themselves
are [here](results/).**

## Setup

First clone the entire project repository:

    git clone https://github.com/bodkan/beluga-simulations
    cd beluga-simulations

Then, running this in your R console will install all project
dependencies:

``` r
renv::restore()
```

## Toy simulations of $Ne_e$ vs $\pi$

``` r
library(dplyr)
library(ggplot2)
```

![](figures/unnamed-chunk-3-1.png)<!-- -->

### Running the simulations

When run, this script will generate the file `pi_simple.rds`:

    Rscript pi_simple.R

### Results

``` r
pi_simple <- readRDS("pi_simple.rds")
```

``` r
pi_simple %>% filter(model == "constant Ne") %>%
  ggplot(aes(factor(Ne), diversity, color = model)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +
  # geom_jitter(position = position_dodge(width = 0.8), alpha = 0.2, size = 0.75) +
  geom_smooth(aes(group = model), method = "lm", linewidth = 0.5, color = "black", linetype = 2) +
  theme(axis.text.x = element_text(hjust = 1, angle = 45), legend.position = "bottom") +
  labs(x = "Ne", y = "nucleotide diversity",
       title = "Expected nucleotide diversity as a function of Ne")
#> `geom_smooth()` using formula = 'y ~ x'
```

![](figures/unnamed-chunk-5-1.png)<!-- -->

``` r
pi_simple %>% filter(model == "constant Ne" | grepl("100 gens", model)) %>%
  ggplot(aes(factor(Ne), diversity, color = model)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +
  # geom_jitter(position = position_dodge(width = 0.8), alpha = 0.2, size = 0.75) +
  geom_smooth(aes(group = model), method = "lm", linewidth = 0.5, color = "black", linetype = 2) +
  scale_color_discrete(drop = FALSE) +
  theme(axis.text.x = element_text(hjust = 1, angle = 45), legend.position = "bottom") +
  labs(x = "Ne", y = "nucleotide diversity",
       title = "Expected nucleotide diversity as a function of Ne")
#> `geom_smooth()` using formula = 'y ~ x'
```

![](figures/unnamed-chunk-6-1.png)<!-- -->

``` r
pi_simple %>% filter(model == "constant Ne" | grepl("1000 gens", model)) %>%
  ggplot(aes(factor(Ne), diversity, color = model)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +
  # geom_jitter(position = position_dodge(width = 0.8), alpha = 0.2, size = 0.75) +
  geom_smooth(aes(group = model), method = "lm", linewidth = 0.5, color = "black", linetype = 2) +
  scale_color_discrete(drop = FALSE) +
  theme(axis.text.x = element_text(hjust = 1, angle = 45), legend.position = "bottom") +
  labs(x = "Ne", y = "nucleotide diversity",
       title = "Expected nucleotide diversity as a function of Ne")
#> `geom_smooth()` using formula = 'y ~ x'
```

![](figures/unnamed-chunk-7-1.png)<!-- -->

``` r
pi_simple %>% filter(model == "constant Ne" | grepl("1000 gens", model) | grepl("2000 gens", model)) %>%
  ggplot(aes(factor(Ne), diversity, color = model)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +
  # geom_jitter(position = position_dodge(width = 0.8), alpha = 0.2, size = 0.75) +
  geom_smooth(aes(group = model), method = "lm", linewidth = 0.5, color = "black", linetype = 2) +
  scale_color_discrete(drop = FALSE) +
  theme(axis.text.x = element_text(hjust = 1, angle = 45), legend.position = "bottom") +
  labs(x = "Ne", y = "nucleotide diversity",
       title = "Expected nucleotide diversity as a function of Ne")
#> `geom_smooth()` using formula = 'y ~ x'
```

![](figures/unnamed-chunk-8-1.png)<!-- -->

``` r
ggplot(pi_simple, aes(factor(Ne), diversity, color = model)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +
  # geom_jitter(position = position_dodge(width = 0.8), alpha = 0.2, size = 0.75) +
  geom_smooth(aes(group = model), method = "lm", linewidth = 0.5, color = "black", linetype = 2) +
  scale_color_discrete(drop = FALSE) +
  theme(axis.text.x = element_text(hjust = 1, angle = 45), legend.position = "bottom") +
  labs(x = "Ne", y = "nucleotide diversity",
       title = "Expected nucleotide diversity as a function of Ne")
#> `geom_smooth()` using formula = 'y ~ x'
```

![](figures/unnamed-chunk-9-1.png)<!-- -->

## Simulations of Beluga diversity given a level of hunting

![](figures/unnamed-chunk-10-1.png)<!-- -->

``` r
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(cowplot)
library(slendr)
```

### Running the simulations

When run, this script will generate the file `pi_simple.rds`:

    Rscript pi_beluga.R

### Results

``` r
pi_beluga <- readRDS("pi_beluga_50Mb.rds")# %>% filter(engine == "msprime")

source("plotting.R")
```

``` r
unique(pi_beluga$N_start)
#> [1] 40000 30000 20000 10000
```

``` r
unique(pi_beluga$N_hunted)
#> [1]  100  250  500  750 1000 1500 2000 3000
```

``` r
unique(pi_beluga$census_ratio)
#>  [1] 1.0 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1
```

    > unique(pi_beluga$N_start)
    10000 20000 30000 40000
    # use: 40000
    > unique(pi_beluga$N_hunted)
    100  250  500  750 1000 1500 2000 3000
    # use: 250  500  1000 2500
    > unique(pi_beluga$census_ratio)
    0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0
    # use: 0.3 0.5 0.7

``` r
plot_grid(
  plot_slopes(pi_beluga, "msprime") + coord_cartesian(ylim = c(-1e-8, 1e-8)),
  plot_slopes(pi_beluga, "SLiM") + coord_cartesian(ylim = c(-1e-8, 1e-8)),
  nrow = 2
)
```

![](figures/unnamed-chunk-14-1.png)<!-- -->

**As a sanity check, all models have been run with both `msprime()` and
`slim()` engines of *slendr*. All results are expected to be the same
regardless of whether they are run in a coalescent or forward-time
setting.**

#### Starting $N_e$ = 40000, $N_e$ hunted = 100, census ratio = 1.0

``` r
plot_panels(pi_beluga, N_start = 40000, N_hunted = 100, census_ratio = 1.0, engine = "msprime")
```

![](figures/unnamed-chunk-15-1.png)<!-- -->

``` r
plot_panels(pi_beluga, N_start = 40000, N_hunted = 100, census_ratio = 1.0, engine = "SLiM")
```

![](figures/unnamed-chunk-16-1.png)<!-- -->

#### Starting $N_e$ = 40000, $N_e$ hunted = 1000, census ratio = 1.0

``` r
plot_panels(pi_beluga, N_start = 40000, N_hunted = 1000, census_ratio = 1.0, engine = "msprime")
```

![](figures/unnamed-chunk-17-1.png)<!-- -->

``` r
plot_panels(pi_beluga, N_start = 40000, N_hunted = 1000, census_ratio = 1.0, engine = "SLiM")
```

![](figures/unnamed-chunk-18-1.png)<!-- -->

#### Starting $N_e$ = 10000, $N_e$ hunted = 100, census ratio = 1.0

``` r
plot_panels(pi_beluga, N_start = 10000, N_hunted = 100, census_ratio = 1.0, engine = "msprime")
```

![](figures/unnamed-chunk-19-1.png)<!-- -->

``` r
plot_panels(pi_beluga, N_start = 10000, N_hunted = 100, census_ratio = 1.0, engine = "SLiM")
```

![](figures/unnamed-chunk-20-1.png)<!-- -->

#### Starting $N_e$ = 10000, $N_e$ hunted = 300, census ratio = 1.0

``` r
plot_panels(pi_beluga, N_start = 10000, N_hunted = 250, census_ratio = 1.0, engine = "msprime")
```

![](figures/unnamed-chunk-21-1.png)<!-- -->

``` r
plot_panels(pi_beluga, N_start = 10000, N_hunted = 250, census_ratio = 1.0, engine = "SLiM")
```

![](figures/unnamed-chunk-22-1.png)<!-- -->

## Linear fit metrics across all models

``` r
markers <- tibble(metric = c("slope", "p.value"), value = c(0, 0.05))

pi_beluga %>%
  unnest(lm) %>% 
  pivot_longer(cols = c("r.squared", "slope", "p.value"), names_to = "metric") %>%
  ggplot(aes(value, color = engine)) +
  geom_density() +
  geom_vline(data = markers, aes(xintercept = value), linetype = "dashed") +
  facet_wrap(metric ~ N_start, scales = "free") +
  ggtitle("p-value, R-squared, and slope vs starting N and hunted N")
```

![](figures/unnamed-chunk-23-1.png)<!-- -->

In all cases, $R^2$ values as well as slope estimates of the linear
models are miniscule, indicating basically no effect of time on
nucleotide diversity. If there’s any effect, it pops up when the
starting population size is towards the lower end of the range (10.000)
and even then there’s practically no measurable effect size.

In fact, we can zoom in on the $p < 0.05$ models and inspect that the
slopes and $R^2$ are really mostly just noise:

``` r
p_p_vs_rsquared <-
  pi_beluga %>%
  unnest(lm) %>%
  ggplot(aes(p.value, r.squared * 100)) +
  geom_point(alpha = 0.5) +
  labs(x = "p-value of the linear regression: nucleotide diversity ~ time",
       y = expression(R^2~~"[% of variance explained]")) +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "red") +
  coord_cartesian(ylim = c(0, 5)) +
  ggtitle("p-value of linear regression vs goodness of fit") +
  facet_wrap(~ engine)

p_p_vs_slope <-
  pi_beluga %>%
  unnest(lm) %>%
  ggplot(aes(p.value, slope, aes(color = N_start))) +
  geom_point(alpha = 0.5) +
  labs(x = "p-value of the linear regression: nucleotide diversity ~ time",
       y = "slope of the linear regression") +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_cartesian(ylim = c(-2e-8, 2e-8)) +
  ggtitle("p-value of linear regression vs effect size (slope)") +
  facet_wrap(~ engine)

cowplot::plot_grid(p_p_vs_rsquared, p_p_vs_slope, nrow = 2)
```

![](figures/unnamed-chunk-24-1.png)<!-- -->

*msprime* vs SLiM sanity check – if slendr correctly interprets a
compiled model, both simulation engines should give the same result in
terms of simulated nucleotide diversities:

``` r
pi_beluga %>%
  unnest(pi) %>%
  select(-model, -lm) %>%
  group_by(N_start, N_hunted, census_ratio, engine) %>%
  summarise(pi = mean(pi)) %>%
  ungroup %>%
  pivot_wider(names_from = "engine", values_from = "pi") %>%
  mutate(N_start = paste("N start =", N_start)) %>%
  # sample_frac(0.1) %>%
  ggplot(aes(SLiM, msprime, color = factor(N_hunted))) +
  geom_point(alpha = 0.8) +
  geom_abline(slope = 1, linetype = "dashed", color = "black") +
  facet_wrap(~ N_start, scales = "free", nrow = 1) +
  theme(legend.position = "bottom") +
  labs(x = expression(paste(pi, " as simulated by SLiM")),
       y = expression(paste(pi, " as simulated by msprime")))
#> `summarise()` has grouped output by 'N_start', 'N_hunted', 'census_ratio'. You
#> can override using the `.groups` argument.
```

![](figures/unnamed-chunk-25-1.png)<!-- -->

## Empirical comparison

``` r
library(dplyr)
library(readr)

pi_empirical <- read_tsv("pi_empirical.txt") %>%
  pivot_longer(
    cols = everything(),
    names_to = "snapshot", names_prefix = "pi_",
    values_to = "pi"
  ) %>% 
  mutate(snapshot = ifelse(snapshot == "Present", "Contemporary", snapshot),
         snapshot = factor(snapshot,
                           levels = c("Time1", "Time2", "Time3", "Contemporary", "Anadyr", "Bristol", "Cook")))
#> Rows: 1275 Columns: 7
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: "\t"
#> dbl (7): pi_Time1, pi_Time2, pi_Time3, pi_Present, pi_Anadyr, pi_Cook, pi_Br...
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

### Raw results from Figure 3a

``` r
pi_empirical%>%
ggplot(aes(snapshot, pi, fill = snapshot)) +
  geom_violin()
```

![](figures/unnamed-chunk-27-1.png)<!-- -->

Subset to empirical diversities only in Mackenzie belugas:

``` r
pi_empirical <- filter(pi_empirical, snapshot %in% c("Time1", "Time2", "Time3", "Contemporary"))
```

Normalize diversities relative to the first time snapshot:

``` r
pi_empirical$pi_relative <- pi_empirical$pi / mean(pi_empirical[pi_empirical$snapshot == "Time1", ]$pi)
```

``` r
pi_empirical%>%
ggplot(aes(snapshot, pi_relative, fill = snapshot)) +
  geom_violin()
```

![](figures/unnamed-chunk-30-1.png)<!-- -->

### Empirical vs simulated diversities

Subset simulation results to those roughly matching the empirical data:

``` r
pi_simulated <- pi_beluga %>%
  filter(engine == "msprime") %>%
  unnest(pi) %>%
  select(-model, -engine) %>%
  mutate(snapshot = case_when(
    time >= 1290 & time <= 1440 ~ "Snapshot #1",
    time >= 1450 & time <= 1650 ~ "Snapshot #2",
    time >= 1800 & time <= 1870 ~ "Snapshot #3",
    time >= 1950 & time <= 2025 ~ "Contemporary",
    TRUE ~ "other"
  )) %>%
  filter(snapshot != "other") %>%
  mutate(snapshot = factor(snapshot,
                         levels = c("Snapshot #1", "Snapshot #2", "Snapshot #3", "Contemporary"))) %>%
  select(-name, -pop, -time) %>%
  filter( # filter based on Mikkel's ideas for the parameter grid
    N_start == 40000,
    N_hunted %in% c(250, 500, 1000, 2500),
    census_ratio %in% c(0.3, 0.5, 0.7)
  )
```

Normalize nucleotide diversities in each simulation vs the mean in the
first respective snapshot:

``` r
pi_simulated <-
  pi_simulated %>%
  group_by(N_start, N_hunted, census_ratio) %>%
  mutate(mean_pi = mean(pi[snapshot == "Snapshot #1"])) %>%
  ungroup() %>%
  mutate(pi_relative = pi / mean_pi) %>%
  select(-mean_pi)
```

## Results for the paper

### Supplementary Figure \#x

**Rendered version is [here](results/sim_violins.pdf)**

**Caption:** *Distributions of nucleotide diversity across four temporal
snapshots of a simulated beluga population corresponding to time points
in Figure 3a. Each panel displays results from each individual
simulation scenario. The starting size of the population was always
40.000. The number of individuals removed from the population each
generation due to hunting is indicated in the title of each panel, as is
the ratio used to convert these two quantities to the $N_e$ simulated in
each scenario.*

``` r
pi_simulated_labels <-
  pi_simulated %>%
  mutate(
    census_ratio = factor(paste("Ne : census ratio =", census_ratio)),
    N_hunted = paste("N hunted =", N_hunted) %>% factor(., levels = gtools::mixedsort(unique(.))),
  )

p_sim_violins <-
  pi_simulated_labels %>%
  ggplot(aes(snapshot, pi_relative)) +
  geom_violin(aes(fill = snapshot), color = NA, alpha = 0.8) +
  geom_smooth(method = "lm", aes(group = 1), fullrange = TRUE,
              color = "black", se = FALSE, linetype = "dashed") +
  facet_wrap(N_hunted ~ census_ratio, nrow = 3) +
  theme(legend.position = "right", axis.text.x = element_text(hjust = 1, angle = 45)) +
  guides(fill = guide_legend("Time point")) +
  labs(
    x = "",
    y = expression(paste(pi, " relative to time point #1")),
    # title = paste("Simulated nucleotide diversity across time assuming\nstarting size of 40.000 individuals")
  ) +
  coord_cartesian(ylim = c(0.8, 1.15)); p_sim_violins
#> `geom_smooth()` using formula = 'y ~ x'
```

![](figures/unnamed-chunk-33-1.png)<!-- -->

``` r

ggsave("results/sim_violins.pdf", p_sim_violins, width = 7, height = 8)
#> `geom_smooth()` using formula = 'y ~ x'
```

Potentially a Figure \#2 but the table below is the rare case where a
plain text table works better:

``` r
markers <- tibble(metric = c("slope", "p-value", "R-squared"), value = c(0, 0.05, 0))

lm_metrics <- pi_simulated_labels %>%
  select(-pi, -pi_relative, -snapshot) %>%
  unnest(lm) %>% 
  pivot_longer(cols = c("r.squared", "slope", "p.value"), names_to = "metric") %>%
  mutate(metric = case_when(
    metric == "r.squared" ~ "R-squared",
    metric == "p.value" ~ "p-value",
    metric == "slope" ~ "slope"
  )) %>%
  distinct()

p_metrics <- lm_metrics %>%
  ggplot(aes(value, color = metric)) +
  geom_density() +
  geom_vline(data = markers, aes(xintercept = value), linetype = "dashed") +
  facet_wrap(metric ~ ., scales = "free", nrow = 1) +
  theme(legend.position = "none") +
  ggtitle("Distribution of metrics of the linear fit 'time ~ nucleotide diversity'")
```

### Supplementary Table \#y

**Rendered version is [here](results/table.png)**

**Data is [here](results/table.tsv)**

**Caption:** *Metrics of the linear regression fit of simulated
nucleotide diversity as a function of time for each simulation scenario
shown in panels in Figure Sx.*

``` r
lm_table <- pi_simulated_labels %>%
  select(-pi, -pi_relative, -snapshot) %>%
  unnest(lm) %>% 
  mutate(N_hunted = as.numeric(gsub("N hunted = ", "", N_hunted)),
         census_ratio = as.numeric(gsub("Ne : census ratio = ", "", census_ratio))) %>%
  distinct() %>%
  select(N_start, N_hunted, census_ratio, p.value, r.squared, slope) %>%
  mutate(
    p.value = format(p.value, digits = 2, nsmall = 2),
    r.squared = format(r.squared, digits = 1),
    slope = format(slope, digits = 1)
  ) %>%
  arrange(N_start, N_hunted, census_ratio)

lm_table
#> # A tibble: 9 × 6
#>   N_start N_hunted census_ratio p.value r.squared slope   
#>     <dbl>    <dbl>        <dbl> <chr>   <chr>     <chr>   
#> 1   40000      250          0.3 0.42    3e-04     "-1e-09"
#> 2   40000      250          0.5 0.92    5e-06     "-3e-10"
#> 3   40000      250          0.7 0.25    7e-04     "-3e-09"
#> 4   40000      500          0.3 0.66    1e-04     "-8e-10"
#> 5   40000      500          0.5 0.91    7e-06     "-2e-10"
#> 6   40000      500          0.7 0.32    5e-04     "-2e-09"
#> 7   40000     1000          0.3 0.63    1e-04     " 1e-09"
#> 8   40000     1000          0.5 1.00    8e-09     "-9e-12"
#> 9   40000     1000          0.7 0.81    3e-05     " 6e-10"
```

``` r
library(gt)
library(webshot2)
#> 
#> Attaching package: 'webshot2'
#> The following object is masked from 'package:slendr':
#> 
#>     resize
```

``` r

gt_table <- lm_table %>%
  gt(auto_align = FALSE) %>%
  cols_label(
    N_start = md("**N start**"),
    N_hunted = md("**N hunted**"),
    census_ratio = md("**Ne : census ratio**"),
    p.value = md("**p-value**"),
    r.squared = md("**R-squared**"),
    slope = md("**slope**")
  )

write_tsv(lm_table, "results/table.tsv")
gtsave(gt_table, "results/table.png")
```

### Supplementary Figure \#z

**Rendered version is [here](results/model.pdf)**

**Caption:** *A schematic representation of a *slendr* model used to
simulate the effect of hunting on the beluga demographic history. The
starting size of the population was always 40.000. The number of
individuals removed from the population each generation due to hunting
is indicated by each “step” decrease in population size.*

``` r
p_model <- pi_beluga %>%
  filter(N_start == 40000, N_hunted == 250, census_ratio == 1.0, engine == "msprime") %>%
  { .$model[[1]] } %>%
  plot_model(); p_model
```

![](figures/unnamed-chunk-37-1.png)<!-- -->

``` r

ggsave("results/model.pdf", p_model, width = 6, height = 7)
```

### Total count of simulation scenarios

``` r
pi_simulated %>%
  group_by(N_start, N_hunted, census_ratio) %>%
  slice(1) %>%
  select(-c(pi, lm, snapshot, pi_relative))
#> # A tibble: 9 × 3
#> # Groups:   N_start, N_hunted, census_ratio [9]
#>   N_start N_hunted census_ratio
#>     <dbl>    <dbl>        <dbl>
#> 1   40000      250          0.3
#> 2   40000      250          0.5
#> 3   40000      250          0.7
#> 4   40000      500          0.3
#> 5   40000      500          0.5
#> 6   40000      500          0.7
#> 7   40000     1000          0.3
#> 8   40000     1000          0.5
#> 9   40000     1000          0.7
```
