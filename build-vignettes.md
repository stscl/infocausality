Due to the time-consuming computations involved in the vignettes of the *infocausality* package, 
it is necessary to pre-build the vignettes prior to package submission.

``` r
.prebuild_vignettes = \(name){
  out = paste0("vignettes/",name,".Rmd")
  inp = paste0(out,".orig")
  knitr::knit(inp,out)
}
```

-   Compile all vignettes at once

``` r
for (v in c("surd")) {
  .prebuild_vignettes(v)
}
```

-   Build vignettes separately

``` r
.prebuild_vignettes("surd")
```
