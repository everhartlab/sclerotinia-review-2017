---
title: "MCG assessment of Kamvar et al 2017"
author: "Zhian N. Kamvar"
date: "2017-11-21"
output: github_document
bibliography: bibliography.bib
editor_options: 
  chunk_output_type: inline
---



Introduction
=============

This document will test the hypothesis that mycelial compatibility groups will
deliniate sexually recombining populations of *Sclerotinia sclerotiorum*. In
this document, we are using data from @kamvar2017data, which consists of 366
isolates of *S. sclerotiorum* sampled over 11 states in the United States of
America, as well as Australia, France, and Mexico.

If this hypothesis is true, we expect to find no signatures of linkage within
the data before or after clone correction. We will be measuring clone-correction
via the standardized index of association, $\bar{r}_d$ as implemented in the
*poppr* package [@agapow2001indices; @kamvar2014poppr].

Packages and Data
-----------------

Here we are loading the data and packages necessary for the analyses. 


```r
library("poppr")
library("tidyverse")
```

Now that the packages are loaded, we can load the data:


```r
data_cols <- cols(
  .default = col_integer(),
  Severity = col_double(),
  Region = col_character(),
  Source = col_character(),
  Host = col_character()
)
dat <- readr::read_csv(here::here("data/kamvar2017population.csv"), col_types = data_cols)
dat
```

```
## # A tibble: 366 x 23
##    Severity   MCG Region Source  Year  Host Isolate `5-2(F)` `5-3(F)` `6-2(F)` `7-2(F)`
##       <dbl> <int>  <chr>  <chr> <int> <chr>   <int>    <int>    <int>    <int>    <int>
##  1      3.9     4     NE    unk  2003    GH     152      320      328      489      172
##  2      5.4    45     NE    unk  2003    GH     274      320      328      489      172
##  3      6.3     5     NY    unk  2003    GH     443      324      308      483      172
##  4      4.4     4     MN    wmn  2003  G122     444      320      328      489      172
##  5      4.7     4     MN    wmn  2003 Beryl     445      320      328      489      172
##  6      6.1     3     MI    wmn  2003 Beryl     446      322      339      483      172
##  7      5.5     5     MI    wmn  2003 Beryl     447      322      308      483      172
##  8      5.0     3     MI    wmn  2003 Beryl     448      324      339      483      172
##  9      5.2     3     MI    wmn  2003 Bunsi     449      322      339      483      172
## 10      5.3     5     MI    wmn  2003 Bunsi     450      322      308      483      172
## # ... with 356 more rows, and 12 more variables: `8-3(H)` <int>, `9-2(F)` <int>,
## #   `12-2(H)` <int>, `17-3(H)` <int>, `20-3(F)` <int>, `36-4(F)` <int>, `50-4(F)` <int>,
## #   `55-4(F)` <int>, `92-4(F)` <int>, `106-4(H)` <int>, `110-4(H)` <int>,
## #   `114-4(H)` <int>
```

Since these are the data for the 16 loci, we only want to keep the 11 that were
used for the study:


```r
replen <- c(
"5-2(F)" = 2,
"5-3(F)" = 4,
"6-2(F)" = 5.99999,
"7-2(F)" = 2,
"8-3(H)" = 2,
"9-2(F)" = 2,
"12-2(H)" = 2,
"17-3(H)" = 3,
"20-3(F)" = 2,
"36-4(F)" = 4,
"50-4(F)" = 4,
"55-4(F)" = 4,
"92-4(F)" = 2,
"106-4(H)" = 4,
"110-4(H)" = 3.99999,
"114-4(H)" = 4
)
loci_to_keep <- c("5-2(F)", "6-2(F)", "7-2(F)", "8-3(H)", "9-2(F)", "12-2(H)", 
"17-3(H)", "20-3(F)", "55-4(F)", "110-4(H)", "114-4(H)")
```

Now we can use these to subset our data:


```r
dat11 <- dat %>% 
  select(loci_to_keep) %>%
  df2genind(strata = select(dat, Region, Source, Host, MCG, Year), 
            ind.names = dat$Isolate,
            ploidy = 1) %>%
  as.genclone()
stopifnot(nLoc(dat11) == 11L)
stopifnot(nmll(dat11) == 165L)
other(dat11)$REPLEN <- replen
other(dat11)$meta <- dat %>% select(Severity, Isolate)
setPop(dat11) <- ~Region
dat11
```

```
## 
## This is a genclone object
## -------------------------
## Genotype information:
## 
##    165 original multilocus genotypes 
##    366 haploid individuals
##     11 codominant loci
## 
## Population information:
## 
##      5 strata - Region, Source, Host, MCG, Year
##     14 populations defined - NE, NY, MN, ..., France, Mexico, ND
```

Index of Association
====================

By MCG
--------

When we assess this by MCG, we need to first ensure that we are not drastically
reducing our sample size because there is evidence that small sample sizes
reduces the power of the index of association and can make clonal populations
appear to be sexual.


```r
dat11 %>% 
  setPop(~MCG) %>%                   # set population to MCGs
  selPopSize(n = 10) %>%             # constrain to 10 samples per population
  poppr(sample = 999, total = FALSE) # test index of association for each pop
```

![plot of chunk by-mcg](./figures/kamvar2017population//by-mcg-1.png)

```
##   Pop  N MLG eMLG    SE    H     G lambda   E.5  Hexp   Ia  p.Ia rbarD  p.rD File
## 1   4 14   9 6.93 0.849 1.97  5.44  0.816 0.724 0.431 2.96 0.001 0.302 0.001    .
## 2  45 16   7 5.25 0.881 1.56  3.37  0.703 0.630 0.192 2.52 0.001 0.286 0.001    .
## 3   5 73  37 7.72 1.266 3.05 10.68  0.906 0.481 0.408 2.74 0.001 0.292 0.001    .
## 4   2 10   9 9.00 0.000 2.16  8.33  0.880 0.952 0.590 1.43 0.001 0.144 0.001    .
## 5  44 36  19 7.37 1.201 2.57  8.53  0.883 0.625 0.338 2.93 0.001 0.330 0.001    .
## 6   1 15  10 7.47 0.876 2.15  7.26  0.862 0.822 0.451 1.94 0.001 0.243 0.001    .
## 7  49 11   6 5.73 0.445 1.64  4.48  0.777 0.836 0.274 3.87 0.001 0.487 0.001    .
## 8   9 15   8 5.67 0.943 1.60  3.17  0.684 0.549 0.390 5.23 0.001 0.592 0.001    .
```

In terms of our hypothesis, that was underwhelming... What happens if we clone
correct our data? We'll use the scheme in [@kamvar2017data], adding MCG as the
highest level:


```r
dat11 %>% 
  setPop(~MCG) %>%                   # set population to MCGs
  selPopSize(n = 10) %>%             # constrain to 10 samples per population
  clonecorrect(~MCG/Region/Source/Host/Year, keep = 1) %>%
  poppr(sample = 999, total = FALSE) # test index of association for each pop
```

![plot of chunk by-mcg-cc](./figures/kamvar2017population//by-mcg-cc-1.png)

```
##   Pop  N MLG eMLG    SE    H     G lambda   E.5  Hexp   Ia  p.Ia rbarD  p.rD File
## 1   4 14   9 6.93 0.849 1.97  5.44  0.816 0.724 0.431 2.96 0.001 0.302 0.001    .
## 2  45 12   7 6.15 0.657 1.70  4.24  0.764 0.724 0.240 2.34 0.002 0.264 0.002    .
## 3   5 60  37 8.69 1.017 3.32 19.57  0.949 0.696 0.451 2.27 0.001 0.239 0.001    .
## 4   2  9   9 9.00 0.000 2.20  9.00  0.889 1.000 0.604 1.41 0.001 0.142 0.001    .
## 5  44 27  19 8.26 1.065 2.72 11.22  0.911 0.717 0.404 2.45 0.001 0.274 0.001    .
## 6   1 12  10 8.50 0.584 2.21  8.00  0.875 0.862 0.471 1.50 0.001 0.188 0.001    .
## 7  49  9   6 6.00 0.000 1.68  4.76  0.790 0.866 0.225 3.30 0.001 0.413 0.001    .
## 8   9 13   8 6.38 0.788 1.74  3.93  0.746 0.625 0.434 4.91 0.001 0.553 0.001    .
```

Okay, maybe that heirarchy is a bit... detailed. What happens if we go the
opposite way? What if we simply clone-corrected just on MCG?


```r
dat11 %>% 
  setPop(~MCG) %>%                   # set population to MCGs
  selPopSize(n = 10) %>%             # constrain to 10 samples per population
  clonecorrect(~MCG) %>% 
  poppr(sample = 999, total = FALSE) # test index of association for each pop
```

![plot of chunk by-mcg-cc-bomb](./figures/kamvar2017population//by-mcg-cc-bomb-1.png)

```
##   Pop  N MLG eMLG       SE    H  G lambda E.5  Hexp   Ia  p.Ia  rbarD  p.rD File
## 1   4  9   9    9 0.00e+00 2.20  9  0.889   1 0.551 1.26 0.002 0.1279 0.002    .
## 2  45  7   7    7 0.00e+00 1.95  7  0.857   1 0.372 1.51 0.009 0.1680 0.009    .
## 3   5 37  37   10 0.00e+00 3.61 37  0.973   1 0.509 1.61 0.001 0.1663 0.001    .
## 4   2  9   9    9 0.00e+00 2.20  9  0.889   1 0.604 1.41 0.001 0.1416 0.001    .
## 5  44 19  19   10 2.51e-07 2.94 19  0.947   1 0.489 1.28 0.001 0.1425 0.001    .
## 6   1 10  10   10 0.00e+00 2.30 10  0.900   1 0.489 0.78 0.003 0.0977 0.003    .
## 7  49  6   6    6 0.00e+00 1.79  6  0.833   1 0.315 3.12 0.003 0.3909 0.003    .
## 8   9  8   8    8 0.00e+00 2.08  8  0.875   1 0.575 3.04 0.001 0.3431 0.001    .
```


Okay. So far, we have no evidence for this hypothesis. But one of the issues
that we saw in Kamvar et al. 2017 was that these data reflected a clonal
population structure. One thing that I'm curious about is the results that were
obtained by @prugnolle2010apparent showing that individuals sampled from 
differentiated populations will have a lower index of association. I can test
this by simulating data.


Data Simulations
----------------

I've written a simulator in the accompanying R script `pop_generator.R` and am
loading it here.


```r
source("code/pop_generator.R")
```

```
## pop_generator is loaded.
## 
## USAGE:
```

```
## pop_generator(n = 100, ploidy = 1, freq = NULL, nall = sample(4:12, 11, replace = TRUE), clone_gen = 0, mu = 0.05, verbose = TRUE, genclone = TRUE)
```

First step is to create some populations. I'm choosing to go through 7 rounds of
clonal reproduction to create 10 populations.


```r
set.seed(2017-11-21)
test <- replicate(10, pop_generator(clone_gen = 7, mu = 0.5)) %>%
  repool() %>%
  as.genclone()
```

```
## I recorded 3 mutation events
## I recorded 2 mutation events
## I recorded 5 mutation events
## I recorded 3 mutation events
## I recorded 0 mutation events
## I recorded 1 mutation events
## I recorded 3 mutation events
## I recorded 2 mutation events
## I recorded 3 mutation events
## I recorded 1 mutation events
```

```r
strata(test) <- data.frame(pop = pop(test))
test
```

```
## 
## This is a genclone object
## -------------------------
## Genotype information:
## 
##     221 original multilocus genotypes 
##    1000 haploid individuals
##      11 codominant loci
## 
## Population information:
## 
##       1 stratum - pop
##      10 populations defined - 
## unknown_1, unknown_2, unknown_3, ..., unknown_8, unknown_9, unknown_10
```

```r
nAll(test)
```

```
##  locus 1  locus 2  locus 3  locus 4  locus 5  locus 6  locus 7  locus 8  locus 9 locus 10 
##       11       12       12        9       12       12       11        9       12       10 
## locus 11 
##       12
```

We can see from this that the clonal reproduction reduced the number of unique
individuals quite a bit. If we test the index of association for these
populations, we can see that they are indeed clonal:


```r
poppr(test, total = FALSE, sample = 999)
```

![plot of chunk ia-sims](./figures/kamvar2017population//ia-sims-1.png)

```
##           Pop   N MLG eMLG SE    H    G lambda   E.5  Hexp   Ia  p.Ia rbarD  p.rD File
## 1   unknown_1 100  24   24  0 3.02 18.2  0.945 0.887 0.763 1.45 0.001 0.147 0.001 test
## 2   unknown_2 100  23   23  0 2.79 12.5  0.920 0.754 0.749 2.41 0.001 0.243 0.001 test
## 3   unknown_3 100  22   22  0 2.73 12.7  0.921 0.813 0.771 2.66 0.001 0.268 0.001 test
## 4   unknown_4 100  28   28  0 3.11 18.6  0.946 0.824 0.799 1.77 0.001 0.179 0.001 test
## 5   unknown_5 100  20   20  0 2.71 12.9  0.923 0.848 0.717 2.13 0.001 0.218 0.001 test
## 6   unknown_6 100  17   17  0 2.50 10.0  0.900 0.803 0.723 2.49 0.001 0.252 0.001 test
## 7   unknown_7 100  17   17  0 2.63 12.4  0.919 0.879 0.748 2.37 0.001 0.242 0.001 test
## 8   unknown_8 100  25   25  0 2.97 16.0  0.938 0.810 0.767 2.01 0.001 0.202 0.001 test
## 9   unknown_9 100  20   20  0 2.59 11.1  0.910 0.812 0.709 2.29 0.001 0.232 0.001 test
## 10 unknown_10 100  25   25  0 2.98 16.5  0.939 0.833 0.763 1.77 0.001 0.180 0.001 test
```

```r
poppr(test, clonecorrect = TRUE, total = FALSE, strata = ~pop, sample = 999)
```

![plot of chunk ia-sims](./figures/kamvar2017population//ia-sims-2.png)

```
##           Pop  N MLG eMLG       SE    H  G lambda E.5  Hexp      Ia  p.Ia    rbarD  p.rD
## 1   unknown_1 24  24   17 0.00e+00 3.18 24  0.958   1 0.796 -0.0937 0.882 -0.00956 0.882
## 2   unknown_2 23  23   17 0.00e+00 3.14 23  0.957   1 0.797  0.0718 0.192  0.00732 0.192
## 3   unknown_3 22  22   17 3.40e-07 3.09 22  0.955   1 0.824  0.3369 0.001  0.03424 0.001
## 4   unknown_4 28  28   17 0.00e+00 3.33 28  0.964   1 0.842  0.1791 0.070  0.01825 0.070
## 5   unknown_5 20  20   17 1.07e-07 3.00 20  0.950   1 0.777  0.0658 0.244  0.00693 0.244
## 6   unknown_6 17  17   17 0.00e+00 2.83 17  0.941   1 0.793 -0.0844 0.768 -0.00884 0.768
## 7   unknown_7 17  17   17 0.00e+00 2.83 17  0.941   1 0.810  0.0829 0.220  0.00866 0.220
## 8   unknown_8 25  25   17 4.73e-07 3.22 25  0.960   1 0.797  0.1176 0.086  0.01194 0.086
## 9   unknown_9 20  20   17 1.07e-07 3.00 20  0.950   1 0.769  0.2008 0.031  0.02075 0.031
## 10 unknown_10 25  25   17 4.73e-07 3.22 25  0.960   1 0.806 -0.0764 0.851 -0.00791 0.851
##    File
## 1  test
## 2  test
## 3  test
## 4  test
## 5  test
## 6  test
## 7  test
## 8  test
## 9  test
## 10 test
```

We can see how they all are related (or not so) to each other


```r
aboot(test, ~pop, sample = 999, dist = "nei.dist")
```

```
## Running bootstraps:       100 / 999Running bootstraps:       200 / 999Running bootstraps:       300 / 999Running bootstraps:       400 / 999Running bootstraps:       500 / 999Running bootstraps:       600 / 999Running bootstraps:       700 / 999Running bootstraps:       800 / 999Running bootstraps:       900 / 999
## Calculating bootstrap values... done.
```

![plot of chunk ia-aboot](./figures/kamvar2017population//ia-aboot-1.png)

```
## 
## Phylogenetic tree with 10 tips and 9 internal nodes.
## 
## Tip labels:
## 	unknown_1, unknown_2, unknown_3, unknown_4, unknown_5, unknown_6, ...
## Node labels:
## 	100, 45.7457457457458, 20.7207207207207, 8.20820820820821, 4.8048048048048, 35.4354354354354, ...
## 
## Rooted; includes branch lengths.
```

Now that we've set up the popualtions, we can randomly sample two individuals
from each and pool them together.


```r
sample_two <- quote(flatten_dbl(map(popNames(test), ~{ sample(which(pop(test) == .), 2) })))
subsamples <- replicate(100, test[eval(sample_two), ])
```


```r
p <- dplyr::progress_estimated(length(subsamples))
res <- map(subsamples, ~{
  p$tick()$print()
  ia(., sample = 999, valuereturn = TRUE, plot = FALSE, quiet = TRUE)
})
p$stop()

resdf <- map(res, 1) %>% 
  map(as.list) %>% 
  bind_rows()
(resdf <- bind_cols(resdf, data_frame(sims = map(res, 2))))
```

```
## # A tibble: 100 x 5
##             Ia  p.Ia        rbarD  p.rD                   sims
##          <dbl> <dbl>        <dbl> <dbl>                 <list>
##  1  0.21039302 0.032  0.021247657 0.032 <data.frame [999 x 2]>
##  2  0.42522718 0.001  0.042789229 0.001 <data.frame [999 x 2]>
##  3  0.45001788 0.001  0.045507524 0.001 <data.frame [999 x 2]>
##  4  0.19283401 0.033  0.019530412 0.033 <data.frame [999 x 2]>
##  5  0.08062331 0.209  0.008183942 0.209 <data.frame [999 x 2]>
##  6  0.15492112 0.374  0.015778779 0.374 <data.frame [999 x 2]>
##  7 -0.01118965 0.529 -0.001126808 0.529 <data.frame [999 x 2]>
##  8  0.10354491 0.170  0.010551664 0.170 <data.frame [999 x 2]>
##  9  0.67430751 0.001  0.068084058 0.001 <data.frame [999 x 2]>
## 10  0.22415784 0.021  0.022824349 0.021 <data.frame [999 x 2]>
## # ... with 90 more rows
```

Now we can see what fraction of the simulations resulted in a significant value
of $\bar{r}_d$


```r
sum(resdf$p.rD <= 0.05)/nrow(resdf)
```

```
## [1] 0.56
```

What do the data look like:


```r
ccia <- poppr(test, quiet = TRUE, total = FALSE)
resdf %>% 
  ggplot(aes(x = rbarD, fill = p.rD >= 0.05)) + 
  geom_histogram(binwidth = 0.01, position = "stack", color = "black") +
  geom_rug(aes(color = p.rD)) +
  viridis::scale_color_viridis() +
  scale_fill_manual(values = c("grey30", "grey80")) +
  geom_vline(xintercept = ccia$rbarD, lty = 3) +
  theme_bw(base_size = 16, base_family = "Helvetica") +
  labs(list(
    fill = expression(paste("p >= 0.05")),
    color = "p-value",
    x = expression(paste(bar(r)[d])),
    title = "The index of association of mixed population samples",
    subtitle = "100 replicates; 10 populations; 2 individuals from each population",
    caption = expression(paste("(dashed lines: observed ", bar(r)[d], " values)"))
  ))
```

![plot of chunk unnamed-chunk-6](./figures/kamvar2017population//unnamed-chunk-6-1.png)


<details>
<summary>Session Information</summary>


```r
devtools::session_info()
```

```
## Session info ----------------------------------------------------------------------------
```

```
##  setting  value                       
##  version  R version 3.4.2 (2017-09-28)
##  system   x86_64, darwin15.6.0        
##  ui       X11                         
##  language (EN)                        
##  collate  en_US.UTF-8                 
##  tz       America/Chicago             
##  date     2017-11-21
```

```
## Packages --------------------------------------------------------------------------------
```

```
##  package     * version     date       source                             
##  ade4        * 1.7-8       2017-11-19 Github (sdray/ade4@2ee45cb)        
##  adegenet    * 2.1.1       2017-10-30 local                              
##  ape           5.0         2017-10-30 CRAN (R 3.4.2)                     
##  assertthat    0.2.0       2017-04-11 CRAN (R 3.4.0)                     
##  backports     1.1.1       2017-09-25 CRAN (R 3.4.2)                     
##  base        * 3.4.2       2017-10-04 local                              
##  bindr         0.1         2016-11-13 CRAN (R 3.4.0)                     
##  bindrcpp      0.2         2017-06-17 CRAN (R 3.4.0)                     
##  boot          1.3-20      2017-07-30 CRAN (R 3.4.1)                     
##  broom         0.4.2       2017-02-13 CRAN (R 3.4.0)                     
##  cellranger    1.1.0       2016-07-27 CRAN (R 3.4.0)                     
##  cli           1.0.0       2017-11-05 CRAN (R 3.4.2)                     
##  cluster       2.0.6       2017-03-16 CRAN (R 3.4.0)                     
##  coda          0.19-1      2016-12-08 CRAN (R 3.4.0)                     
##  codetools     0.2-15      2016-10-05 CRAN (R 3.4.0)                     
##  colorspace    1.3-3       2017-08-16 R-Forge (R 3.4.1)                  
##  compiler      3.4.2       2017-10-04 local                              
##  crayon        1.3.4       2017-09-23 Github (gaborcsardi/crayon@b5221ab)
##  datasets    * 3.4.2       2017-10-04 local                              
##  deldir        0.1-14      2017-04-22 CRAN (R 3.4.0)                     
##  devtools      1.13.3      2017-08-02 CRAN (R 3.4.1)                     
##  digest        0.6.12      2017-01-27 CRAN (R 3.4.0)                     
##  dplyr       * 0.7.4       2017-09-28 CRAN (R 3.4.1)                     
##  evaluate      0.10.1      2017-06-24 CRAN (R 3.4.1)                     
##  expm          0.999-2     2017-03-29 CRAN (R 3.4.0)                     
##  ezknitr       0.6         2016-09-16 CRAN (R 3.4.0)                     
##  fastmatch     1.1-0       2017-01-28 CRAN (R 3.4.0)                     
##  forcats     * 0.2.0       2017-01-23 CRAN (R 3.4.0)                     
##  foreign       0.8-69      2017-06-21 CRAN (R 3.4.0)                     
##  gdata         2.18.0      2017-06-06 CRAN (R 3.4.0)                     
##  ggplot2     * 2.2.1       2016-12-30 CRAN (R 3.4.0)                     
##  glue          1.2.0       2017-10-29 CRAN (R 3.4.2)                     
##  gmodels       2.16.2      2015-07-22 CRAN (R 3.4.0)                     
##  graphics    * 3.4.2       2017-10-04 local                              
##  grDevices   * 3.4.2       2017-10-04 local                              
##  grid          3.4.2       2017-10-04 local                              
##  gridExtra     2.3         2017-09-09 CRAN (R 3.4.1)                     
##  gtable        0.2.0       2016-02-26 CRAN (R 3.4.0)                     
##  gtools        3.5.0       2015-05-29 CRAN (R 3.4.0)                     
##  haven         1.1.0       2017-07-09 CRAN (R 3.4.1)                     
##  here          0.1         2017-05-28 CRAN (R 3.4.0)                     
##  highr         0.6         2016-05-09 CRAN (R 3.4.0)                     
##  hms           0.3         2016-11-22 CRAN (R 3.4.0)                     
##  htmltools     0.3.6       2017-04-28 CRAN (R 3.4.0)                     
##  httpuv        1.3.5       2017-07-04 CRAN (R 3.4.1)                     
##  httr          1.3.1       2017-08-20 cran (@1.3.1)                      
##  igraph        1.1.2       2017-07-21 cran (@1.1.2)                      
##  jsonlite      1.5         2017-06-01 CRAN (R 3.4.0)                     
##  knitr         1.17        2017-08-10 cran (@1.17)                       
##  labeling      0.3         2014-08-23 CRAN (R 3.4.0)                     
##  lattice       0.20-35     2017-03-25 CRAN (R 3.4.0)                     
##  lazyeval      0.2.1       2017-10-29 CRAN (R 3.4.2)                     
##  LearnBayes    2.15        2014-05-29 CRAN (R 3.4.0)                     
##  lubridate     1.7.1       2017-11-03 CRAN (R 3.4.2)                     
##  magrittr      1.5         2014-11-22 CRAN (R 3.4.0)                     
##  MASS          7.3-47      2017-04-21 CRAN (R 3.4.0)                     
##  Matrix        1.2-11      2017-08-16 CRAN (R 3.4.1)                     
##  memoise       1.1.0       2017-04-21 CRAN (R 3.4.0)                     
##  methods     * 3.4.2       2017-10-04 local                              
##  mgcv          1.8-22      2017-09-19 CRAN (R 3.4.2)                     
##  mime          0.5         2016-07-07 CRAN (R 3.4.0)                     
##  mnormt        1.5-5       2016-10-15 CRAN (R 3.4.0)                     
##  modelr        0.1.1       2017-07-24 CRAN (R 3.4.1)                     
##  munsell       0.4.3       2016-02-13 CRAN (R 3.4.0)                     
##  nlme          3.1-131     2017-02-06 CRAN (R 3.4.0)                     
##  parallel      3.4.2       2017-10-04 local                              
##  pegas         0.10        2017-05-03 CRAN (R 3.4.0)                     
##  permute       0.9-4       2016-09-09 CRAN (R 3.4.0)                     
##  phangorn      2.3.1       2017-11-01 CRAN (R 3.4.2)                     
##  pkgconfig     2.0.1       2017-03-21 CRAN (R 3.4.0)                     
##  plyr          1.8.4       2016-06-08 CRAN (R 3.4.0)                     
##  poppr       * 2.5.0.99-29 2017-11-19 Github (grunwaldlab/poppr@21440e6) 
##  psych         1.7.8       2017-09-09 CRAN (R 3.4.1)                     
##  purrr       * 0.2.4       2017-10-18 cran (@0.2.4)                      
##  quadprog      1.5-5       2013-04-17 CRAN (R 3.4.0)                     
##  R.methodsS3   1.7.1       2016-02-16 CRAN (R 3.4.0)                     
##  R.oo          1.21.0      2016-11-01 CRAN (R 3.4.0)                     
##  R.utils       2.6.0       2017-11-05 CRAN (R 3.4.2)                     
##  R6            2.2.2       2017-06-17 cran (@2.2.2)                      
##  Rcpp          0.12.13     2017-09-28 CRAN (R 3.4.2)                     
##  readr       * 1.1.1       2017-05-16 CRAN (R 3.4.0)                     
##  readxl        1.0.0       2017-04-18 CRAN (R 3.4.0)                     
##  reshape2      1.4.2       2016-10-22 CRAN (R 3.4.0)                     
##  rlang         0.1.4       2017-11-05 CRAN (R 3.4.2)                     
##  rprojroot     1.2         2017-01-16 CRAN (R 3.4.0)                     
##  rstudioapi    0.7         2017-09-07 CRAN (R 3.4.1)                     
##  rvest         0.3.2       2016-06-17 CRAN (R 3.4.0)                     
##  scales        0.5.0.9000  2017-08-28 Github (hadley/scales@d767915)     
##  seqinr        3.4-5       2017-08-01 CRAN (R 3.4.1)                     
##  shiny         1.0.5       2017-08-23 cran (@1.0.5)                      
##  sp            1.2-5       2017-06-29 CRAN (R 3.4.1)                     
##  spdep         0.6-15      2017-09-01 CRAN (R 3.4.1)                     
##  splines       3.4.2       2017-10-04 local                              
##  stats       * 3.4.2       2017-10-04 local                              
##  stringi       1.1.5       2017-04-07 CRAN (R 3.4.0)                     
##  stringr     * 1.2.0       2017-02-18 CRAN (R 3.4.0)                     
##  tibble      * 1.3.4       2017-08-22 cran (@1.3.4)                      
##  tidyr       * 0.7.2       2017-10-16 CRAN (R 3.4.2)                     
##  tidyverse   * 1.2.0       2017-11-08 local                              
##  tools         3.4.2       2017-10-04 local                              
##  utils       * 3.4.2       2017-10-04 local                              
##  vegan         2.4-4       2017-08-24 cran (@2.4-4)                      
##  viridis       0.4.0       2017-03-27 CRAN (R 3.4.0)                     
##  viridisLite   0.2.0       2017-03-24 CRAN (R 3.4.0)                     
##  withr         2.1.0       2017-11-01 CRAN (R 3.4.2)                     
##  xml2          1.1.1       2017-01-24 CRAN (R 3.4.0)                     
##  xtable        1.8-2       2016-02-05 CRAN (R 3.4.0)
```

</details>

References
==========
