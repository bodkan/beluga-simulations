Beluga nucleotide diversity
================
2024-05-30

## Setup

Running this will make sure all dependencies of the project will be
setup before running anything else:

``` r
install.packages("renv")

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
