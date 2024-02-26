# CR-Lasso: Robust cellwise regularized sparse regression

The `regcell` package provides the functions to compute the CR-Lasso (cellwise regularized Lasso) proposed by Peng Su, Samuel Muller, Garth Tarr and Suojin Wang. The manuscript can be found [here](https://arxiv.org/abs/2307.05234).

We have included a demonstration (demo), a simulation demonstration (simu) and a real data demonstration (realdata, Bone mineral density data) in vignettes.

We also created an [online R repository](https://posit.cloud/content/7571075) with some example scripts.

### Installation

You can install the package using:

```r
remotes::install_github("PengSU517/regcell", build = FALSE)
```

For macOS users, you may see this error related to `gfortran` when installing:

```
ld: warning: search path '/opt/gfortran/lib/gcc/aarch64-apple-darwin20.0/12.2.0' not found
ld: warning: search path '/opt/gfortran/lib' not found
ld: library 'gfortran' not found
clang: error: linker command failed with exit code 1 (use -v to see invocation)
make: *** [regcell.so] Error 1
ERROR: compilation failed for package ‘regcell’
```

If that happens, please install the GNU Fortran compiler from this page: https://mac.r-project.org/tools/ (for example, this direct link [gfortran-12.2-universal.pkg](https://mac.r-project.org/tools/gfortran-12.2-universal.pkg)) and then try

```r
remotes::install_github("PengSU517/regcell", build = TRUE)
```

