# regcell

- This package provides the functions to compute the CR-Lasso (cellwise regularized Lasso) proposed by Peng Su, Samuel Muller, Garth Tarr and Suojin Wang. The manuscript could be found here soon on Arxiv.

- We added a demonstration (demo) in vignettes

To get started, you can install the package using:

```r
remotes::install_github("PengSU517/regcell", build = FALSE)
```

#### For Mac users

If there are some problems with `gfortran`, you can try the following steps:

Install `gfortran`: Ensure that the gfortran compiler and associated libraries are installed on your system. You can typically install it through a package manager specific to your operating system. For example, on macOS with Homebrew, you can run the following command to install gfortran:

```
brew install gcc
```

Note that gfortran is usually included as part of the GCC (GNU Compiler Collection) package.

Set environment variable: If the gfortran library is installed in a non-standard location, you may need to set the `LD_LIBRARY_PATH` environment variable to include the directory where the library is located. This will help the linker locate the library during the package installation process.

On Unix-like systems (macOS, Linux, etc.), use the following command to set the `LD_LIBRARY_PATH` variable in the current session:

```
export LD_LIBRARY_PATH=/path/to/gfortran/lib:$LD_LIBRARY_PATH
```

Replace `/path/to/gfortran/lib` with the actual path where the gfortran library is installed on your system.


On Windows, you can set the environment variable through the system settings or by using the setx command:

```
setx LD_LIBRARY_PATH "C:\path\to\gfortran\lib;%LD_LIBRARY_PATH%"
```

Replace `C:\path\to\gfortran\lib` with the actual path where the gfortran library is installed on your system.

Retry package installation: After installing gfortran and setting the LD_LIBRARY_PATH environment variable, attempt to install the 'regcell' package again using the following command:

```r
remotes::install_github("PengSU517/regcell", build = FALSE)
```

With the gfortran library properly installed and the environment variable set, the linker should be able to find and link against the library, resolving the previous linker error.



If there are still some errors, you could extract functions from `R` and `src` folders.

