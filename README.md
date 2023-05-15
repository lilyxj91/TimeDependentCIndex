
# TimeDependentCIndex

<!-- badges: start -->
<!-- badges: end -->

This contains R-code and example data for the statistical method for
evaluating a time-dependent C-index for recurent event data, described
in:“Wang J, Jiang X, Ning J. Evaluating dynamic and predictive
discrimination for recurrent event models: use of a time-dependent
C-index. Under Revision.”

## Installation

You can install the development version of TimeDependentCIndex from
[GitHub](https://github.com/) with:

``` r
library(devtools)
#> Loading required package: usethis
devtools::install_github("lilyxj91/TimeDependentCIndex",force = TRUE)
#> Downloading GitHub repo lilyxj91/TimeDependentCIndex@HEAD
#> 
#> ── R CMD build ─────────────────────────────────────────────────────────────────
#>          checking for file 'C:\Users\Xinyang\AppData\Local\Temp\RtmpE3V507\remotes734448d82289\lilyxj91-TimeDependentCIndex-b71d157/DESCRIPTION' ...  ✔  checking for file 'C:\Users\Xinyang\AppData\Local\Temp\RtmpE3V507\remotes734448d82289\lilyxj91-TimeDependentCIndex-b71d157/DESCRIPTION'
#>       ─  preparing 'TimeDependentCIndex':
#>    checking DESCRIPTION meta-information ...     checking DESCRIPTION meta-information ...   ✔  checking DESCRIPTION meta-information
#>       ─  checking for LF line-endings in source and make files and shell scripts
#>   ─  checking for empty or unneeded directories
#>   ─  building 'TimeDependentCIndex_0.1.0.tar.gz'
#>      
#> 
#> Installing package into 'C:/Users/Xinyang/Documents/R/win-library/packages'
#> (as 'lib' is unspecified)
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(TimeDependentCIndex)

### Define the parameters
n.var<-3  ## number of predictors
B<-10  ## number of perturbations. B=10 for illustration purpose, in practice, we use N=200 in general.
train.data = TrainData
validate.data = ValidateData
```

``` r
### Calculate alpha values and perturbation alpha values

results<-TimeDependentCIndex(train.data,validate.data,n.var,B)
```

<img src="man/figures/README-example-1.png" width="100%" />

``` r
alpha<-results[[1]] ## alpha estimation
alpha.pb<-results[[2]] ## alpha perturbation results
alpha
#> [1]  0.70346819  0.11383066 -0.04545526
alpha.pb
#>           X1            X2          X3
#> 1  0.7723826  0.1076185935 -0.17003306
#> 2  0.7554317 -0.0007518938 -0.03216336
#> 3  0.6865039  0.1453919472 -0.07349922
#> 4  0.6255758  0.0850488637 -0.09529885
#> 5  0.7336176  0.1447695312 -0.15584926
#> 6  0.5146179  0.1313476591  0.04232732
#> 7  0.8155319  0.0490439693 -0.13285944
#> 8  0.7305976  0.1041085593 -0.10626370
#> 9  0.7141166  0.2037132461 -0.07484825
#> 10 0.6780435  0.2539866460 -0.09277914
```

``` r

head(train.data)
#>   id         t t.obs indi Num X1         X2        X3
#> 1  1  4.410887    15    0   6  1 -0.3020925 0.4300591
#> 2  1  5.158452    15    0   6  1 -0.3020925 0.4300591
#> 3  1  6.274329    15    0   6  1 -0.3020925 0.4300591
#> 4  1  6.740334    15    0   6  1 -0.3020925 0.4300591
#> 5  1  8.379347    15    0   6  1 -0.3020925 0.4300591
#> 6  1 13.156354    15    0   6  1 -0.3020925 0.4300591
# id is the subject ID; 
# t is the time to recurrent events;
# t.obs is the censoring time;
# indi is the censoring indicator (1-observed; 0-censored);
# Num is the number of recurrent events;
# X1, X2 and X3 are three predictors.  
# 
# Please arrange your data into this format with the same column names.  
```
