## RAutoDiff: an R package for first and second order Automatic Differentiation of R code.

Easy to use first and second order Automatic Differentiation (AD) based on forward mode AD and operator overloading. The package is mainly intended for algorithm prototyping and smaller problems, and not really for production code.
 

### Install
To install the package, run the line
```
devtools::install_github("https://github.com/torekleppe/RAutoDiff",
,build_vignettes=TRUE)
```
in the R console.

### Documentation
The documentation is mainly contained in two vignettes: A "getting started guide":
```
vignette("introduction",package="RAutoDiff")
```
and an overview of the functions currently implemented in the package
```
vignette("functions",package="RAutoDiff")
```
