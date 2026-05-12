## CRAN Package Check Results

On the r-devel Windows win-builder there is a pre-existing compile error in Rcpp/Function.h                  
(R_NamespaceRegistry not declared) that affects every Rcpp-based package on that builder and is unrelated to   
this package. It builds cleanly on r-release and r-oldrel Windows, and on macOS/Linux.

## R CMD check results

0 errors | 0 warnings | 1 note

## revdepcheck results

revdepcheck not available for R version 4.5.3
