### to also run skip_on_cran() tests, uncomment:
Sys.setenv(NOT_CRAN="true")

library(testthat)
library(bumblebee)

test_check("bumblebee")
