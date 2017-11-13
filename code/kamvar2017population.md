Introduction
============

Here we are loading the data and packages necessary for the analyses.

    library("poppr")
    library("tidyverse")

    data_cols <- cols(
      .default = col_integer(),
      Severity = col_double(),
      Region = col_character(),
      Source = col_character(),
      Host = col_character()
    )
    dat <- readr::read_csv(here::here("data/kamvar2017population.csv"), col_types = data_cols)
    dat

    ## # A tibble: 366 x 23
    ##    Severity   MCG Region Source  Year  Host Isolate `5-2(F)` `5-3(F)`
    ##       <dbl> <int>  <chr>  <chr> <int> <chr>   <int>    <int>    <int>
    ##  1      3.9     4     NE    unk  2003    GH     152      320      328
    ##  2      5.4    45     NE    unk  2003    GH     274      320      328
    ##  3      6.3     5     NY    unk  2003    GH     443      324      308
    ##  4      4.4     4     MN    wmn  2003  G122     444      320      328
    ##  5      4.7     4     MN    wmn  2003 Beryl     445      320      328
    ##  6      6.1     3     MI    wmn  2003 Beryl     446      322      339
    ##  7      5.5     5     MI    wmn  2003 Beryl     447      322      308
    ##  8      5.0     3     MI    wmn  2003 Beryl     448      324      339
    ##  9      5.2     3     MI    wmn  2003 Bunsi     449      322      339
    ## 10      5.3     5     MI    wmn  2003 Bunsi     450      322      308
    ## # ... with 356 more rows, and 14 more variables: `6-2(F)` <int>,
    ## #   `7-2(F)` <int>, `8-3(H)` <int>, `9-2(F)` <int>, `12-2(H)` <int>,
    ## #   `17-3(H)` <int>, `20-3(F)` <int>, `36-4(F)` <int>, `50-4(F)` <int>,
    ## #   `55-4(F)` <int>, `92-4(F)` <int>, `106-4(H)` <int>, `110-4(H)` <int>,
    ## #   `114-4(H)` <int>

    options(width = 100)
    devtools::session_info()

    ## Session info --------------------------------------------------------------------------------------

    ##  setting  value                       
    ##  version  R version 3.4.2 (2017-09-28)
    ##  system   x86_64, darwin15.6.0        
    ##  ui       X11                         
    ##  language (EN)                        
    ##  collate  en_US.UTF-8                 
    ##  tz       America/Chicago             
    ##  date     2017-11-13

    ## Packages ------------------------------------------------------------------------------------------

    ##  package    * version     date       source                             
    ##  ade4       * 1.7-8       2017-11-02 Github (sdray/ade4@4348204)        
    ##  adegenet   * 2.1.1       2017-10-30 local                              
    ##  ape          5.0         2017-10-30 CRAN (R 3.4.2)                     
    ##  assertthat   0.2.0       2017-04-11 CRAN (R 3.4.0)                     
    ##  backports    1.1.1       2017-09-25 CRAN (R 3.4.2)                     
    ##  base       * 3.4.2       2017-10-04 local                              
    ##  bindr        0.1         2016-11-13 CRAN (R 3.4.0)                     
    ##  bindrcpp     0.2         2017-06-17 CRAN (R 3.4.0)                     
    ##  boot         1.3-20      2017-07-30 CRAN (R 3.4.1)                     
    ##  broom        0.4.2       2017-02-13 CRAN (R 3.4.0)                     
    ##  cellranger   1.1.0       2016-07-27 CRAN (R 3.4.0)                     
    ##  cli          1.0.0       2017-11-05 CRAN (R 3.4.2)                     
    ##  cluster      2.0.6       2017-03-16 CRAN (R 3.4.0)                     
    ##  coda         0.19-1      2016-12-08 CRAN (R 3.4.0)                     
    ##  colorspace   1.3-3       2017-08-16 R-Forge (R 3.4.1)                  
    ##  compiler     3.4.2       2017-10-04 local                              
    ##  crayon       1.3.4       2017-09-23 Github (gaborcsardi/crayon@b5221ab)
    ##  datasets   * 3.4.2       2017-10-04 local                              
    ##  deldir       0.1-14      2017-04-22 CRAN (R 3.4.0)                     
    ##  devtools     1.13.3      2017-08-02 CRAN (R 3.4.1)                     
    ##  digest       0.6.12      2017-01-27 CRAN (R 3.4.0)                     
    ##  dplyr      * 0.7.4       2017-09-28 CRAN (R 3.4.1)                     
    ##  evaluate     0.10.1      2017-06-24 CRAN (R 3.4.1)                     
    ##  expm         0.999-2     2017-03-29 CRAN (R 3.4.0)                     
    ##  fastmatch    1.1-0       2017-01-28 CRAN (R 3.4.0)                     
    ##  forcats    * 0.2.0       2017-01-23 CRAN (R 3.4.0)                     
    ##  foreign      0.8-69      2017-06-21 CRAN (R 3.4.0)                     
    ##  gdata        2.18.0      2017-06-06 CRAN (R 3.4.0)                     
    ##  ggplot2    * 2.2.1       2016-12-30 CRAN (R 3.4.0)                     
    ##  glue         1.2.0       2017-10-29 CRAN (R 3.4.2)                     
    ##  gmodels      2.16.2      2015-07-22 CRAN (R 3.4.0)                     
    ##  graphics   * 3.4.2       2017-10-04 local                              
    ##  grDevices  * 3.4.2       2017-10-04 local                              
    ##  grid         3.4.2       2017-10-04 local                              
    ##  gtable       0.2.0       2016-02-26 CRAN (R 3.4.0)                     
    ##  gtools       3.5.0       2015-05-29 CRAN (R 3.4.0)                     
    ##  haven        1.1.0       2017-07-09 CRAN (R 3.4.1)                     
    ##  here         0.1         2017-05-28 CRAN (R 3.4.0)                     
    ##  hms          0.3         2016-11-22 CRAN (R 3.4.0)                     
    ##  htmltools    0.3.6       2017-04-28 CRAN (R 3.4.0)                     
    ##  httpuv       1.3.5       2017-07-04 CRAN (R 3.4.1)                     
    ##  httr         1.3.1       2017-08-20 cran (@1.3.1)                      
    ##  igraph       1.1.2       2017-07-21 cran (@1.1.2)                      
    ##  jsonlite     1.5         2017-06-01 CRAN (R 3.4.0)                     
    ##  knitr        1.17        2017-08-10 cran (@1.17)                       
    ##  lattice      0.20-35     2017-03-25 CRAN (R 3.4.0)                     
    ##  lazyeval     0.2.1       2017-10-29 CRAN (R 3.4.2)                     
    ##  LearnBayes   2.15        2014-05-29 CRAN (R 3.4.0)                     
    ##  lubridate    1.7.1       2017-11-03 CRAN (R 3.4.2)                     
    ##  magrittr     1.5         2014-11-22 CRAN (R 3.4.0)                     
    ##  MASS         7.3-47      2017-04-21 CRAN (R 3.4.0)                     
    ##  Matrix       1.2-11      2017-08-16 CRAN (R 3.4.1)                     
    ##  memoise      1.1.0       2017-04-21 CRAN (R 3.4.0)                     
    ##  methods    * 3.4.2       2017-10-04 local                              
    ##  mgcv         1.8-22      2017-09-19 CRAN (R 3.4.2)                     
    ##  mime         0.5         2016-07-07 CRAN (R 3.4.0)                     
    ##  mnormt       1.5-5       2016-10-15 CRAN (R 3.4.0)                     
    ##  modelr       0.1.1       2017-07-24 CRAN (R 3.4.1)                     
    ##  munsell      0.4.3       2016-02-13 CRAN (R 3.4.0)                     
    ##  nlme         3.1-131     2017-02-06 CRAN (R 3.4.0)                     
    ##  parallel     3.4.2       2017-10-04 local                              
    ##  pegas        0.10        2017-05-03 CRAN (R 3.4.0)                     
    ##  permute      0.9-4       2016-09-09 CRAN (R 3.4.0)                     
    ##  phangorn     2.3.1       2017-11-01 CRAN (R 3.4.2)                     
    ##  pkgconfig    2.0.1       2017-03-21 CRAN (R 3.4.0)                     
    ##  plyr         1.8.4       2016-06-08 CRAN (R 3.4.0)                     
    ##  poppr      * 2.5.0.99-12 2017-11-09 local                              
    ##  psych        1.7.8       2017-09-09 CRAN (R 3.4.1)                     
    ##  purrr      * 0.2.4       2017-10-18 cran (@0.2.4)                      
    ##  quadprog     1.5-5       2013-04-17 CRAN (R 3.4.0)                     
    ##  R6           2.2.2       2017-06-17 cran (@2.2.2)                      
    ##  Rcpp         0.12.13.1   2017-10-10 Github (RcppCore/Rcpp@136d50f)     
    ##  readr      * 1.1.1       2017-05-16 CRAN (R 3.4.0)                     
    ##  readxl       1.0.0       2017-04-18 CRAN (R 3.4.0)                     
    ##  reshape2     1.4.2       2016-10-22 CRAN (R 3.4.0)                     
    ##  rlang        0.1.4       2017-11-05 CRAN (R 3.4.2)                     
    ##  rmarkdown    1.7         2017-11-10 cran (@1.7)                        
    ##  rprojroot    1.2         2017-01-16 CRAN (R 3.4.0)                     
    ##  rstudioapi   0.7         2017-09-07 CRAN (R 3.4.1)                     
    ##  rvest        0.3.2       2016-06-17 CRAN (R 3.4.0)                     
    ##  scales       0.5.0.9000  2017-08-28 Github (hadley/scales@d767915)     
    ##  seqinr       3.4-5       2017-08-01 CRAN (R 3.4.1)                     
    ##  shiny        1.0.5       2017-08-23 cran (@1.0.5)                      
    ##  sp           1.2-5       2017-06-29 CRAN (R 3.4.1)                     
    ##  spdep        0.6-15      2017-09-01 CRAN (R 3.4.1)                     
    ##  splines      3.4.2       2017-10-04 local                              
    ##  stats      * 3.4.2       2017-10-04 local                              
    ##  stringi      1.1.5       2017-04-07 CRAN (R 3.4.0)                     
    ##  stringr    * 1.2.0       2017-02-18 CRAN (R 3.4.0)                     
    ##  tibble     * 1.3.4       2017-08-22 cran (@1.3.4)                      
    ##  tidyr      * 0.7.2       2017-10-16 CRAN (R 3.4.2)                     
    ##  tidyverse  * 1.2.0       2017-11-08 local                              
    ##  tools        3.4.2       2017-10-04 local                              
    ##  utils      * 3.4.2       2017-10-04 local                              
    ##  vegan        2.4-4       2017-08-24 cran (@2.4-4)                      
    ##  withr        2.1.0       2017-11-01 CRAN (R 3.4.2)                     
    ##  xml2         1.1.1       2017-01-24 CRAN (R 3.4.0)                     
    ##  xtable       1.8-2       2016-02-05 CRAN (R 3.4.0)                     
    ##  yaml         2.1.14      2016-11-12 CRAN (R 3.4.0)
