# knobi

Authors: Cousido-Rocha, M., Paz, A., Pennino, M.G., Cerviño, S.

Application of a KBPM (Known Biomass Production Model):
(1) the fitting of KBPM to each stock;
(2) the retrospective analysis to identify regime shifts;  and 
(3) the estimation of the effects of environmental variability.

First step is to install `devtools` package using `install.packages("devtools")`.

There are two options for installing our package:

## 1. Version with vignettes

```
devtools::install_github("MERVEX-group/knobi",build_vignettes = TRUE)
```

This option needs to install previously the following packages:

```
install.packages(c("corrplot", "ggplot2",  "gridExtra", "grDevices",  "optimr", "plot3D"))
```

## 2. Version without vignettes

```
devtools::install_github("MERVEX-group/knobi")
```

This is a faster option.

It is recommended to restart R after the installation of the package.
