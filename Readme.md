
# copulagraphicr

This package implements the following research paper: MING ZHENG , JOHN
P. KLEIN, Estimates of marginal survival for dependent competing risks
based on an assumed copula, Biometrika, Volume 82, Issue 1, March 1995,
Pages 127â€“138, <https://doi.org/10.1093/biomet/82.1.127> The package
allows the user to analyze a dataset consisting of observations of two
competing risks. The package contains a function to load the competing
risk data, a function to calculate characteristic quantities of the
competing data and a function to calculate the corresponding
Copula-Graphic estimator for a given Copula. Further the package
contains predefined plots of the sample data and the Copula-Graphic
estimator.

| Version    | Maintainer                                 |
|:-----------|:-------------------------------------------|
| 0.0.0.9000 | Claudio Longo <ClaudioLongo2212@gmail.com> |

# 1 Installation

## 1.1 Released Version

All released versions can be found at
<https://github.com/CLongoTUM/copulagraphicr>.

From this location the released version can be installed using package
drat:

<!-- note that the package must be pushed there, using -->
<!-- drat::insertPackage(file = devtools::build(),
repodir = paste0("M:/", desc::desc()$get_field("Organisationseinheit"))) -->

``` r
## install latest released package version from git repository
remotes::install_github("https://github.com/CLongoTUM/copulagraphicr.git")
```

## 1.2 Initialize object of `copula_graphic` class

``` r
library(copulagraphicr)
copula_graphic <- new("copula_graphic")
```

## 1.3 Initialize sample data and estimates for k, p_1 and p_2.

``` r
## load package sample data (competing risk data, package data)
copula_graphic@sample_data <- copulagraphicr::load_data()

## visualize competing risk data
copulagraphicr::visualize_data(
  sample_data = copula_graphic@sample_data
)
```

![](README_files/figure-gfm/plot%20sample%20data-1.png)<!-- -->

``` r
## estimate empirical quantities from sample data
est <- copulagraphicr::estimate_quantities(
  sample_data = copula_graphic@sample_data
)

copula_graphic@empirical_quantities$t_grid <- est[[1]]
copula_graphic@empirical_quantities$k <- est[[2]]
copula_graphic@empirical_quantities$p_1 <- est[[3]]
copula_graphic@empirical_quantities$p_2 <- est[[4]]
```

![](README_files/figure-gfm/plot%20initialize%20quantities-1.png)<!-- -->

## 1.4 The Copula-Graphic estimator for different assumed copulas

``` r
## initialize error terms and tau for Copula-Graphic estimator
error_A <- 1e-3
error_B <- 1e-3
tau <- 0.5

## calculate Copula-Graphic estimator for different assumed copulas
copula_available <- c("Frank", "Gumbel", "Joe", "Independence")
est_fit <- lapply(copula_available, function(copula) {
  return(
    copulagraphicr::CG_fit(
      k = copula_graphic@empirical_quantities$k,
      p_1 = copula_graphic@empirical_quantities$p_1,
      t_grid = copula_graphic@empirical_quantities$t_grid,
      copula = copula,
      error_A = error_A,
      error_B = error_B,
      tau = tau
    )
  )
})

est_Frank <- est_fit[[1]]
est_Gumbel <- est_fit[[2]]
est_Joe <- est_fit[[3]]
est_Ind <- est_fit[[4]]

## extract Copula-Graphic estimate of survival function
S_Frank <- 1 - est_Frank[[1]]
S_Gumbel <- 1 - est_Gumbel[[1]]
S_Joe <- 1 - est_Joe[[1]]
S_Ind <- 1 - est_Ind[[1]]
```

## 1.5 Kaplan-Meier estimator as benchmark

``` r
## standard estimator: survival fit, Kaplan-Meier for censored data
d <- rep(0, nrow(copula_graphic@sample_data))
x <- copula_graphic@sample_data$T
e <- copula_graphic@sample_data$delta
fit <- survival::survfit(survival::Surv(d, x, e) ~ 1)
```

## 1.6 Plot of Copula-Graphic estimates

``` r
## Plot
y_lim_lower <- round(min(S_Frank, S_Gumbel, S_Frank, S_Ind) - 0.15, 1)

coord <- graphics::par("usr")
graphics::plot(
  copula_graphic@empirical_quantities$t_grid[1:length(S_Ind)],
  S_Ind,
  type = "s",
  ylim = c(y_lim_lower, 1),
  xlab = "time",
  ylab = "survival function",
  col = "white"
)
graphics::lines(
  fit,
  col = "red",
  conf.int = FALSE,
  type = "s"
)
graphics::lines(
  copula_graphic@empirical_quantities$t_grid[1:length(S_Frank)],
  S_Frank,
  col = "orange",
  type = "s"
)
graphics::lines(
  copula_graphic@empirical_quantities$t_grid[1:length(S_Gumbel)],
  S_Gumbel,
  col = "purple",
  type = "s"
)
graphics::lines(
  copula_graphic@empirical_quantities$t_grid[1:length(S_Joe)],
  S_Joe,
  col = "green",
  type = "s"
)
graphics::legend(
  0,
  y_lim_lower+0.25,
  legend = c("Independence", "Frank's", "Gumbel", "Joe"),
  c("red", "orange", "purple", "green"),
  cex = 1,
  box.lty = 0,
  lty = 1,
  fill = 0,
  border = 0,
  x.intersp = 0.2,
  y.intersp = 0.8
)
graphics::title(
  "Copula-Graphic estimator for different assumed copulas"
)
```

![](README_files/figure-gfm/plot%20CG%20estimate-1.png)<!-- -->
