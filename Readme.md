
# copulagraphicr

The package allows the user to analyze a dataset consisting of
observations of two competing risks. The package contains a function to
load the competing risk data, a function to calculate characteristic
quantities of the data and a function to calculate the corresponding
Copula-Graphic estimator for a given Copula.

| Version    | Maintainer                                 |
|:-----------|:-------------------------------------------|
| 0.0.0.9000 | Claudio Longo <ClaudioLongo2212@gmail.com> |

# 1 Installation

## 1.1 Released Version

All released versions can be found at xxx.

From this location the released version can be installed using package
drat:

<!-- note that the package must be pushed there, using -->
<!-- drat::insertPackage(file = devtools::build(),
repodir = paste0("M:/", desc::desc()$get_field("Organisationseinheit"))) -->

``` r
# install latest released package version from central repo
install.packages("copulagraphicr,
                 type = "source",
                 repos = paste0("file:M:/", params$pkg_org))
```
