# regcell

- This package provides the functions to compute the CR-Lasso (cellwise regularized Lasso) proposed by Peng Su, Samuel Muller, Garth Tarr and Suojin Wang. The manuscript could be found soon on Arxiv.

- We added a demonstration (demo) in vignettes.

- We also created an online R repository with some example scripts.  https://posit.cloud/content/6051440

To get started, you can install the package using:

```r
remotes::install_github("PengSU517/regcell", build = FALSE)
```

For macOS users, if there are some problems with `gfortran`, you can try install the GNU Fortran compiler from this page: https://mac.r-project.org/tools/.


If there are still some errors, you could extract functions from `R` and `src` folders.

