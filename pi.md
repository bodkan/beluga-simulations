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

![](pi_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

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
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](pi_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

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
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](pi_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

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
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](pi_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
ggplot(pi_simple, aes(factor(Ne), diversity, color = model)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +
  # geom_jitter(position = position_dodge(width = 0.8), alpha = 0.2, size = 0.75) +
  geom_smooth(aes(group = model), method = "lm", linewidth = 0.5, color = "black", linetype = 2) +
  scale_color_discrete(drop = FALSE) +
  theme(axis.text.x = element_text(hjust = 1, angle = 45), legend.position = "bottom") +
  labs(x = "Ne", y = "nucleotide diversity",
       title = "Expected nucleotide diversity as a function of Ne")
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](pi_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

## Simulations of Beluga diversity given a level of hunting

![](pi_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

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
unique(pi_beluga$Ne_start)
```

    ## [1] 10000 20000 30000 40000

``` r
unique(pi_beluga$Ne_hunted)
```

    ## [1]  100  250  500  750 1000 1500 2000 3000

``` r
unique(pi_beluga$census_ratio)
```

    ##  [1] 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0

    > unique(pi_beluga$Ne_start)
    10000 20000 30000 40000
    > unique(pi_beluga$Ne_hunted)
    100  250  500  750 1000 1500 2000 3000
    > unique(pi_beluga$census_ratio)
    0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0

``` r
plot_grid(
  plot_slopes(pi_beluga, "msprime") + labs(subtitle = "msprime engine") + coord_cartesian(ylim = c(-1e-8, 1e-8)),
  plot_slopes(pi_beluga, "SLiM") + labs(subsubtitle ="SLiM engine") + coord_cartesian(ylim = c(-1e-8, 1e-8)),
  nrow = 2
)
```

![](pi_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

#### Starting $N_e$ = 40000, $N_e$ hunted = 100, census ratio = 1.0

``` r
plot_panels(pi_beluga, Ne_start = 40000, Ne_hunted = 100, census_ratio = 1.0, engine = "msprime")
```

![](pi_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

<!-- ```{r} -->
<!-- plot_panels(pi_beluga, Ne_start = 40000, Ne_hunted = 100, census_ratio = 1.0, engine = "SLiM") -->
<!-- ``` -->

#### Starting $N_e$ = 40000, $N_e$ hunted = 1000, census ratio = 1.0

``` r
plot_panels(pi_beluga, Ne_start = 40000, Ne_hunted = 1000, census_ratio = 1.0, engine = "msprime")
```

![](pi_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

<!-- ```{r} -->
<!-- plot_panels(pi_beluga, Ne_start = 40000, Ne_hunted = 1000, census_ratio = 1.0, engine = "SLiM") -->
<!-- ``` -->

#### Starting $N_e$ = 10000, $N_e$ hunted = 100, census ratio = 1.0

``` r
plot_panels(pi_beluga, Ne_start = 10000, Ne_hunted = 100, census_ratio = 1.0, engine = "msprime")
```

![](pi_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

<!-- ```{r} -->
<!-- plot_panels(pi_beluga, Ne_start = 10000, Ne_hunted = 100, census_ratio = 1.0, engine = "SLiM") -->
<!-- ``` -->

#### Starting $N_e$ = 10000, $N_e$ hunted = 300, census ratio = 1.0

``` r
plot_panels(pi_beluga, Ne_start = 10000, Ne_hunted = 250, census_ratio = 1.0, engine = "msprime")
```

![](pi_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

<!-- ```{r} -->
<!-- plot_panels(pi_beluga, Ne_start = 10000, Ne_hunted = 250, census_ratio = 1.0, engine = "SLiM") -->
<!-- ``` -->

``` r
markers <- tibble(statistic = c("slope", "p.value"), value = c(0, 0.05))

bind_rows(pi_beluga$lm) %>%
  as_tibble %>%
  pivot_longer(cols = c("r.squared", "slope", "p.value"), names_to = "statistic") %>%
  ggplot(aes(value, fill = statistic)) +
  geom_density() +
  geom_vline(data = markers, aes(xintercept = value), linetype = "dashed") +
  facet_wrap(~ statistic, scales = "free")
```

![](pi_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

``` r
p_p_vs_rsquared <-
  filter(pi_beluga, engine == "msprime") %>%
  unnest(lm) %>%
  ggplot(aes(p.value, r.squared * 100)) +
  geom_point(alpha = 0.5) +
  labs(x = "p-value of the linear regression: nucleotide diversity ~ time",
       y = expression(R^2~~"[% of variance explained]")) +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "red") +
  coord_cartesian(ylim = c(0, 5)) +
  ggtitle("p-value of linear regression vs goodness of fit")

p_p_vs_slope <-
  filter(pi_beluga, engine == "msprime") %>%
  unnest(lm) %>%
  ggplot(aes(p.value, slope, aes(color = Ne_start))) +
  geom_point(alpha = 0.5) +
  labs(x = "p-value of the linear regression: nucleotide diversity ~ time",
       y = "slope of the linear regression") +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_cartesian(ylim = c(-2e-8, 2e-8)) +
  ggtitle("p-value of linear regression vs effect size (slope)")

cowplot::plot_grid(p_p_vs_rsquared, p_p_vs_slope, nrow = 2)
```

![](pi_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->
