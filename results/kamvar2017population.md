---
title: "MCG assessment of Kamvar et al 2017"
author: "Zhian N. Kamvar"
date: "2017-12-05"
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
library("ggridges")
library("ape")
library("cowplot")
library("gridGraphics")
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
## 4   2 10   9 9.00 0.000 2.16  8.33  0.880 0.952 0.590 1.43 0.002 0.144 0.002    .
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
## 2  45 12   7 6.15 0.657 1.70  4.24  0.764 0.724 0.240 2.34 0.001 0.264 0.001    .
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
## 2  45  7   7    7 0.00e+00 1.95  7  0.857   1 0.372 1.51 0.006 0.1680 0.006    .
## 3   5 37  37   10 0.00e+00 3.61 37  0.973   1 0.509 1.61 0.001 0.1663 0.001    .
## 4   2  9   9    9 0.00e+00 2.20  9  0.889   1 0.604 1.41 0.001 0.1416 0.001    .
## 5  44 19  19   10 2.51e-07 2.94 19  0.947   1 0.489 1.28 0.001 0.1425 0.001    .
## 6   1 10  10   10 0.00e+00 2.30 10  0.900   1 0.489 0.78 0.002 0.0977 0.002    .
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

I've written a simulator in an R package called "kop" and am loading it here:


```r
library("doParallel")
```

```
## Loading required package: foreach
```

```
## 
## Attaching package: 'foreach'
```

```
## The following objects are masked from 'package:purrr':
## 
##     accumulate, when
```

```
## Loading required package: iterators
```

```
## Loading required package: parallel
```

```r
library("kop")
```

First step is to create some populations. Because I want to get a representative
from each population, I'm going to create 20 populations. These populations are
initially simulated from a multinomial distribution, but this often results in
extremely long branches for a tree:


```r
kop::pop_generator(mate_gen = 0, mu = 0.5) %>% 
  aboot(sample = 1, tree = "nj", dist = "dist", quiet = TRUE) %>%
  invisible()
```

```
## Beginning mating
```

```
## 
## I recorded 0 mutation events
```

![plot of chunk example_sim](./figures/kamvar2017population//example_sim-1.png)


To avoid this, I'm parameterizing the simulations thusly:

| Parameter | Value |
| --------- | ----- |
| Census Size | 1000 |
| Sample Size | 100 |
| Generations of random mating | 400 |
| Generations of clonal reproduction (after random mating) | 100 |
| Mutations/Generation (over all samples/loci) | 0.5 |

The random mating serves to shorten those long terminal branches. This is
important for ensuring that the populations we DO simulate are sufficiently 
different from each other.


```r
cl <- makeCluster(4)
registerDoParallel(cl, cores = 4)
set.seed(2017-11-21)
test <- foreach(seq_len(20), .combine = c, .packages = c("kop", "poppr", "dplyr", "purrr", "tibble")) %dopar% 
  pop_generator(n = 1000, mate_gen = 400, clone_gen = 100, mu = 0.5, verbose = FALSE)
stopCluster(cl)
test <- test %>%
  repool() %>%
  as.genclone()
strata(test) <- data.frame(pop = pop(test))
```


```r
test
```

```
## 
## This is a genclone object
## -------------------------
## Genotype information:
## 
##     341 original multilocus genotypes 
##    2000 haploid individuals
##      11 codominant loci
## 
## Population information:
## 
##       1 stratum - pop
##      20 populations defined - 
## unknown_1, unknown_2, unknown_3, ..., unknown_18, unknown_19, unknown_20
```

```r
nAll(test)
```

```
##  locus 1  locus 2  locus 3  locus 4  locus 5  locus 6  locus 7  locus 8  locus 9 locus 10 
##       12       11       12       10       11       11       10       10       10       12 
## locus 11 
##       10
```

```r
plot(ape::nj(diss.dist(test, percent = TRUE)), lab4ut = "axial", type = "unrooted", tip.col = adegenet::funky(nPop(test))[pop(test)])
```

![plot of chunk poptree](./figures/kamvar2017population//poptree-1.png)


We can see from this that the clonal reproduction reduced the number of unique
individuals quite a bit. If we test the index of association for these
populations, we can see that they are indeed clonal:


```r
poppr(test, total = FALSE, sample = 999) 
```

![plot of chunk ia-sims](./figures/kamvar2017population//ia-sims-1.png)

```
##           Pop   N MLG eMLG SE    H     G lambda   E.5  Hexp    Ia  p.Ia  rbarD  p.rD File
## 1   unknown_1 100  19   19  0 2.56 10.37  0.904 0.787 0.507 1.208 0.001 0.1345 0.001 test
## 2   unknown_2 100  15   15  0 2.40  9.31  0.893 0.829 0.541 1.276 0.001 0.1419 0.001 test
## 3   unknown_3 100  16   16  0 2.50 10.31  0.903 0.830 0.511 1.342 0.001 0.1365 0.001 test
## 4   unknown_4 100  17   17  0 2.63 12.47  0.920 0.887 0.487 1.078 0.001 0.1207 0.001 test
## 5   unknown_5 100  20   20  0 2.71 12.53  0.920 0.823 0.469 0.453 0.001 0.0455 0.001 test
## 6   unknown_6 100  15   15  0 2.48 10.42  0.904 0.859 0.467 1.003 0.001 0.1116 0.001 test
## 7   unknown_7 100  16   16  0 2.33  7.40  0.865 0.688 0.400 1.068 0.001 0.1232 0.001 test
## 8   unknown_8 100  16   16  0 2.54 10.99  0.909 0.853 0.512 1.192 0.001 0.1235 0.001 test
## 9   unknown_9 100  14   14  0 2.09  5.68  0.824 0.661 0.426 1.415 0.001 0.1425 0.001 test
## 10 unknown_10 100  16   16  0 2.44 10.10  0.901 0.865 0.525 1.452 0.001 0.1615 0.001 test
## 11 unknown_11 100  18   18  0 2.42  7.69  0.870 0.655 0.435 1.324 0.001 0.1332 0.001 test
## 12 unknown_12 100  17   17  0 2.48  9.11  0.890 0.745 0.484 1.079 0.001 0.1201 0.001 test
## 13 unknown_13 100  17   17  0 2.55 10.64  0.906 0.813 0.600 1.717 0.001 0.1722 0.001 test
## 14 unknown_14 100  16   16  0 2.46 10.14  0.901 0.850 0.429 0.767 0.001 0.0773 0.001 test
## 15 unknown_15 100  23   23  0 2.85 14.29  0.930 0.815 0.560 0.901 0.001 0.0906 0.001 test
## 16 unknown_16 100  17   17  0 2.42  8.18  0.878 0.699 0.482 1.097 0.001 0.1102 0.001 test
## 17 unknown_17 100  18   18  0 2.49  8.77  0.886 0.705 0.459 1.293 0.001 0.1298 0.001 test
## 18 unknown_18 100  18   18  0 2.55 10.20  0.902 0.783 0.597 1.691 0.001 0.1891 0.001 test
## 19 unknown_19 100  18   18  0 2.58 10.82  0.908 0.805 0.442 0.907 0.001 0.0967 0.001 test
## 20 unknown_20 100  15   15  0 2.37  8.14  0.877 0.739 0.416 1.265 0.001 0.1481 0.001 test
```

```r
poppr(test, clonecorrect = TRUE, total = FALSE, strata = ~pop, sample = 999)
```

![plot of chunk ia-sims](./figures/kamvar2017population//ia-sims-2.png)

```
##           Pop  N MLG eMLG       SE    H  G lambda E.5  Hexp       Ia  p.Ia    rbarD  p.rD
## 1   unknown_1 19  19   14 3.79e-07 2.94 19  0.947   1 0.537 -0.01019 0.537 -0.00114 0.537
## 2   unknown_2 15  15   14 0.00e+00 2.71 15  0.933   1 0.587 -0.08677 0.761 -0.00967 0.761
## 3   unknown_3 16  16   14 0.00e+00 2.77 16  0.938   1 0.590  0.14437 0.278  0.01452 0.278
## 4   unknown_4 17  17   14 1.60e-07 2.83 17  0.941   1 0.524  0.30389 0.033  0.03390 0.033
## 5   unknown_5 20  20   14 0.00e+00 3.00 20  0.950   1 0.512 -0.24573 0.994 -0.02461 0.994
## 6   unknown_6 15  15   14 0.00e+00 2.71 15  0.933   1 0.532  0.13975 0.189  0.01557 0.189
## 7   unknown_7 16  16   14 0.00e+00 2.77 16  0.938   1 0.433 -0.04173 0.578 -0.00470 0.578
## 8   unknown_8 16  16   14 0.00e+00 2.77 16  0.938   1 0.580  0.06271 0.296  0.00635 0.296
## 9   unknown_9 14  14   14 0.00e+00 2.64 14  0.929   1 0.517 -0.00459 0.582 -0.00046 0.582
## 10 unknown_10 16  16   14 0.00e+00 2.77 16  0.938   1 0.575  0.19450 0.086  0.02166 0.086
## 11 unknown_11 18  18   14 0.00e+00 2.89 18  0.944   1 0.465  0.29236 0.046  0.02933 0.046
## 12 unknown_12 17  17   14 1.60e-07 2.83 17  0.941   1 0.541 -0.11771 0.876 -0.01315 0.876
## 13 unknown_13 17  17   14 1.60e-07 2.83 17  0.941   1 0.632  0.45957 0.002  0.04632 0.002
## 14 unknown_14 16  16   14 0.00e+00 2.77 16  0.938   1 0.485  0.34701 0.032  0.03509 0.032
## 15 unknown_15 23  23   14 0.00e+00 3.14 23  0.957   1 0.587  0.05140 0.360  0.00516 0.360
## 16 unknown_16 17  17   14 1.60e-07 2.83 17  0.941   1 0.581 -0.17684 0.928 -0.01774 0.928
## 17 unknown_17 18  18   14 0.00e+00 2.89 18  0.944   1 0.572  0.07258 0.293  0.00728 0.293
## 18 unknown_18 18  18   14 0.00e+00 2.89 18  0.944   1 0.633  0.11256 0.176  0.01267 0.176
## 19 unknown_19 18  18   14 0.00e+00 2.89 18  0.944   1 0.527 -0.07293 0.698 -0.00740 0.698
## 20 unknown_20 15  15   14 0.00e+00 2.71 15  0.933   1 0.483  0.25481 0.068  0.02861 0.068
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
## 11 test
## 12 test
## 13 test
## 14 test
## 15 test
## 16 test
## 17 test
## 18 test
## 19 test
## 20 test
```

We can see how they all are related (or not so) to each other


```r
aboot(test, ~pop, sample = 1000, dist = "nei.dist", tree = "nj")
```

```
## Warning in infinite_vals_replacement(D, warning): Infinite values detected.

## Warning in infinite_vals_replacement(D, warning): Infinite values detected.

## Warning in infinite_vals_replacement(D, warning): Infinite values detected.

## Warning in infinite_vals_replacement(D, warning): Infinite values detected.

## Warning in infinite_vals_replacement(D, warning): Infinite values detected.

## Warning in infinite_vals_replacement(D, warning): Infinite values detected.

## Warning in infinite_vals_replacement(D, warning): Infinite values detected.

## Warning in infinite_vals_replacement(D, warning): Infinite values detected.
```

```
## Running bootstraps:       100 / 1000
```

```
## Warning in infinite_vals_replacement(D, warning): Infinite values detected.

## Warning in infinite_vals_replacement(D, warning): Infinite values detected.

## Warning in infinite_vals_replacement(D, warning): Infinite values detected.

## Warning in infinite_vals_replacement(D, warning): Infinite values detected.

## Warning in infinite_vals_replacement(D, warning): Infinite values detected.
```

```
## Running bootstraps:       200 / 1000
```

```
## Warning in infinite_vals_replacement(D, warning): Infinite values detected.

## Warning in infinite_vals_replacement(D, warning): Infinite values detected.

## Warning in infinite_vals_replacement(D, warning): Infinite values detected.

## Warning in infinite_vals_replacement(D, warning): Infinite values detected.

## Warning in infinite_vals_replacement(D, warning): Infinite values detected.

## Warning in infinite_vals_replacement(D, warning): Infinite values detected.

## Warning in infinite_vals_replacement(D, warning): Infinite values detected.

## Warning in infinite_vals_replacement(D, warning): Infinite values detected.
```

```
## Running bootstraps:       300 / 1000
```

```
## Warning in infinite_vals_replacement(D, warning): Infinite values detected.

## Warning in infinite_vals_replacement(D, warning): Infinite values detected.
```

```
## Running bootstraps:       400 / 1000
```

```
## Warning in infinite_vals_replacement(D, warning): Infinite values detected.

## Warning in infinite_vals_replacement(D, warning): Infinite values detected.

## Warning in infinite_vals_replacement(D, warning): Infinite values detected.

## Warning in infinite_vals_replacement(D, warning): Infinite values detected.

## Warning in infinite_vals_replacement(D, warning): Infinite values detected.

## Warning in infinite_vals_replacement(D, warning): Infinite values detected.

## Warning in infinite_vals_replacement(D, warning): Infinite values detected.
```

```
## Running bootstraps:       500 / 1000
```

```
## Warning in infinite_vals_replacement(D, warning): Infinite values detected.

## Warning in infinite_vals_replacement(D, warning): Infinite values detected.

## Warning in infinite_vals_replacement(D, warning): Infinite values detected.

## Warning in infinite_vals_replacement(D, warning): Infinite values detected.

## Warning in infinite_vals_replacement(D, warning): Infinite values detected.

## Warning in infinite_vals_replacement(D, warning): Infinite values detected.
```

```
## Running bootstraps:       600 / 1000
```

```
## Warning in infinite_vals_replacement(D, warning): Infinite values detected.

## Warning in infinite_vals_replacement(D, warning): Infinite values detected.

## Warning in infinite_vals_replacement(D, warning): Infinite values detected.

## Warning in infinite_vals_replacement(D, warning): Infinite values detected.

## Warning in infinite_vals_replacement(D, warning): Infinite values detected.

## Warning in infinite_vals_replacement(D, warning): Infinite values detected.

## Warning in infinite_vals_replacement(D, warning): Infinite values detected.
```

```
## Running bootstraps:       700 / 1000
```

```
## Warning in infinite_vals_replacement(D, warning): Infinite values detected.

## Warning in infinite_vals_replacement(D, warning): Infinite values detected.

## Warning in infinite_vals_replacement(D, warning): Infinite values detected.

## Warning in infinite_vals_replacement(D, warning): Infinite values detected.
```

```
## Running bootstraps:       800 / 1000
```

```
## Warning in infinite_vals_replacement(D, warning): Infinite values detected.

## Warning in infinite_vals_replacement(D, warning): Infinite values detected.

## Warning in infinite_vals_replacement(D, warning): Infinite values detected.

## Warning in infinite_vals_replacement(D, warning): Infinite values detected.
```

```
## Running bootstraps:       900 / 1000
```

```
## Warning in infinite_vals_replacement(D, warning): Infinite values detected.

## Warning in infinite_vals_replacement(D, warning): Infinite values detected.

## Warning in infinite_vals_replacement(D, warning): Infinite values detected.

## Warning in infinite_vals_replacement(D, warning): Infinite values detected.

## Warning in infinite_vals_replacement(D, warning): Infinite values detected.

## Warning in infinite_vals_replacement(D, warning): Infinite values detected.
```

```
## Running bootstraps:       1000 / 1000
## Calculating bootstrap values... done.
```

![plot of chunk ia-aboot](./figures/kamvar2017population//ia-aboot-1.png)

```
## 
## Phylogenetic tree with 20 tips and 18 internal nodes.
## 
## Tip labels:
## 	unknown_1, unknown_2, unknown_3, unknown_4, unknown_5, unknown_6, ...
## Node labels:
## 	100, 0.2, 3.3, 7.1, 8.2, 3.1, ...
## 
## Unrooted; includes branch lengths.
```

Now that we've set up the popualtions, we can randomly sample two individuals
from each and pool them together.


```r
# Sample one individual per population
sample_one <- quote(flatten_dbl(map(popNames(test), ~{ sample(which(pop(test) == .), 1) })))
subsamples <- replicate(100, test[eval(sample_one), ])
```


```r
cl <- makeCluster(4)
registerDoParallel(cl, cores = 4)
set.seed(2017-11-21)
res <- foreach(i = seq(subsamples), .combine = c, .packages = c("poppr")) %dopar% 
  list(ia(subsamples[[i]], sample = 999, valuereturn = TRUE, plot = FALSE, quiet = TRUE))
stopCluster(cl)

# reshaping everything into a data frame
resdf <- map(res, 1) %>% 
  map(as.list) %>% 
  bind_rows()
(resdf <- bind_cols(resdf, data_frame(sims = map(res, 2))))
```

```
## # A tibble: 100 x 5
##              Ia  p.Ia         rbarD  p.rD                   sims
##           <dbl> <dbl>         <dbl> <dbl>                 <list>
##  1 -0.120020796 0.899 -0.0122073685 0.899 <data.frame [999 x 2]>
##  2  0.069219049 0.250  0.0070207221 0.250 <data.frame [999 x 2]>
##  3  0.016096664 0.819  0.0016269359 0.819 <data.frame [999 x 2]>
##  4  0.111602639 0.127  0.0113122685 0.127 <data.frame [999 x 2]>
##  5 -0.005363237 0.508 -0.0005461089 0.508 <data.frame [999 x 2]>
##  6  0.021070984 0.415  0.0021449863 0.415 <data.frame [999 x 2]>
##  7  0.088731109 0.574  0.0090058169 0.574 <data.frame [999 x 2]>
##  8  0.011790492 0.435  0.0011895302 0.435 <data.frame [999 x 2]>
##  9 -0.061584697 0.739 -0.0062729773 0.739 <data.frame [999 x 2]>
## 10  0.024161602 0.408  0.0024630431 0.408 <data.frame [999 x 2]>
## # ... with 90 more rows
```


> This simulation take a long time to run. When you come back to it, use the 
> function `lazierLoad("cache")` to reload all of the saved results:
> `source("https://gist.githubusercontent.com/zkamvar/655be26d1d5bec1ac948/raw/97db57030d016a0e7a8156aa84ad02a7c8bc644c/lazierLoad.R")`

Now we can see what fraction of the simulations resulted in a significant value
of $\bar{r}_d$


```r
sum(resdf$p.rD <= 0.05)/nrow(resdf)
```

```
## [1] 0.05
```

That's not very many.

What do the data look like:


```r
ccia <- poppr(test, quiet = TRUE, total = FALSE)
the_plot <- resdf %>% 
  ggplot(aes(x = rbarD, fill = p.rD >= 0.05)) + 
  geom_histogram(binwidth = 0.01, position = "stack", color = "black") +
  geom_rug(aes(color = p.rD)) +
  viridis::scale_color_viridis() +
  scale_fill_manual(values = c("grey30", "grey80")) +
  theme_bw(base_size = 16, base_family = "Helvetica") +
  labs(list(
    fill = expression(paste("p >= 0.05")),
    color = "p-value",
    x = expression(paste(bar(r)[d])),
    title = "The index of association of mixed population samples",
    subtitle = "100 replicates; 20 populations; 1 individual from each population",
    caption = expression(paste("(dashed lines: observed ", bar(r)[d], " values)"))
  ))
the_plot + geom_vline(xintercept = ccia$rbarD, lty = 3)
```

![plot of chunk unnamed-chunk-6](./figures/kamvar2017population//unnamed-chunk-6-1.png)


What happens when we randomly sample 20 individuals from each population?


```r
set.seed(2017-12-04)
resampled_rbarD <- seppop(test) %>%
  map_df(resample.ia, n = 20, .id = "Population") %>%
  as_tibble()
resampled_rbarD_res <- resampled_rbarD %>% 
  group_by(Population) %>% 
  summarize(rbarD = mean(rbarD)) %>% 
  pull(rbarD)
the_plot + 
  geom_vline(xintercept = resampled_rbarD_res, lty = 3) +
  labs(caption = expression(paste("(dashed lines: resampled ", bar(r)[d], " values)")))
```

![plot of chunk randsamp](./figures/kamvar2017population//randsamp-1.png)


# Attempting to make pretty plots


Because the ulitmate goal here is to make graphics that may be understandable,
I'm goint to attempt to create a ridgeline density plot of $\bar{r}_d$ values
and combine that with a tree of the 20 random samples overlayed on random
samples from all populations.


```r
ridgelines <- resdf %>% 
  select(Ia, rbarD) %>% 
  add_column(Population = "pooled") %>% 
  bind_rows(resampled_rbarD, .id = "Data") %>% 
  mutate(Data = case_when(Population == "pooled" ~ "pooled", TRUE ~ "single")) %>% 
  # mutate(Population = case_when(
  #          grepl("unknown", Population) ~ sprintf("pop %2d", as.integer(gsub("unknown_", "", Population))),
  #          TRUE ~ Population
  #          )) %>% 
  group_by(Population) %>% 
  mutate(m = mean(rbarD)) %>% 
  arrange(m) %>% 
  select(-m) %>% 
  ungroup() %>% 
  mutate(Population = fct_inorder(Population))

## Function for desaturating colors by specified proportion
desat <- function(cols, sat=0.5) {
    X <- diag(c(1, sat, 1)) %*% rgb2hsv(col2rgb(cols))
    hsv(X[1,], X[2,], X[3,])
}

popcolors <- setNames(adegenet::funky(nPop(test)), adegenet::popNames(test))

p_ridge <- ggplot(ridgelines, aes(x = rbarD, y = Population, fill = Population, height = ..density..)) + 
  geom_density_ridges(scale = 5) +
  theme_ridges(font_size = 16, font_family = "Helvetica", grid = TRUE, center_axis_labels = TRUE) +
  scale_x_continuous(limits = c(NA, 0.4), breaks = c(0, 0.2, 0.4)) +
  scale_y_discrete(expand = c(0.01, 0)) +
  theme(aspect.ratio = 1.25) +
  theme(legend.position = "top") +
  theme(legend.justification = "center") +
  theme(axis.text.y = element_blank()) +
  theme(axis.ticks.y = element_blank()) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(panel.grid.major.y = element_blank()) +
  theme(panel.grid.major.x = element_line(linetype = 3, color = "grey50")) +
  scale_fill_manual(values = c(pooled = "black", popcolors)) +
  labs(list(
    x = expression(paste(italic(bar(r)[d])))
  ))
p_ridge
```

```
## Picking joint bandwidth of 0.0088
```

```
## Warning: Removed 2 rows containing non-finite values (stat_density_ridges).
```

![plot of chunk unnamed-chunk-7](./figures/kamvar2017population//unnamed-chunk-7-1.png)

```r
ggsave(p_ridge, filename = here::here("results/figures/p-ridge.pdf"))
```

```
## Saving 3.5 x 5 in image
## Picking joint bandwidth of 0.0088
```

```
## Warning: Removed 2 rows containing non-finite values (stat_density_ridges).
```


Here we need to create a tree that is digestable, which means creating one that
we can easily label. I am first taking one of the data sets that had random
subsamples and am resampling 19 individuals from each of the 20 populations and
creating a tree from the resulting 400 individuals. From this tree, I will
highlight the edges that are found in the subsample and compare those to the
trees created from within-population subsample.


```r
s1 <- subsamples[[1]]
set.seed(2017-12-04)
samps <- data_frame(pop = pop(test), ind = indNames(test)) %>% 
  filter(!ind %in% indNames(s1)) %>%
  group_by(pop) %>%
  sample_n(19) %>%
  bind_rows(unsampled = ., sampled = data_frame(pop = pop(s1), ind = indNames(s1)), .id = "tree")

ftree  <- diss.dist(test[samps$ind], percent = TRUE) %>% nj()
cols    <- setNames(adegenet::funky(nPop(test)), popNames(test))
fullcol  <- samps %>% mutate(cols = case_when(tree == "sampled" ~ cols[pop], TRUE ~ transp(cols[pop], 0.5))) %>% pull(cols)
fullsize <- samps %>% mutate(cols = case_when(tree == "sampled" ~ 2, TRUE ~ 0.5)) %>% pull(cols)
newedge <- seq(ftree$edge) %in% which.edge(ftree, tail(samps$ind, 20))
ew      <- ifelse(newedge, 2, 1)
ec      <- ifelse(newedge, "black", transp("grey50", 0.5))
maintree <- function(){
  plot.phylo(ftree, 
             type = "unrooted", 
             show.tip.label = FALSE, 
             no.margin = TRUE, 
             edge.width = ew, 
             edge.color = ec)
  tiplabels(pch = 21, bg = fullcol, col = NA, cex = fullsize)
  add.scale.bar(length = 0.1, lwd = 2)
}
subtrees <- function(){
  for (i in popNames(test)){
    inds <- samps$ind[samps$pop == i]
    test[inds] %>% 
      diss.dist(percent = TRUE) %>%
      nj() %>%
      plot.phylo(type = "unrooted", 
                 show.tip.label = FALSE, 
                 no.margin = TRUE, 
                 edge.width = 2,
                 edge.color = cols[i])
    add.scale.bar(length = 0.1, lwd = 2)
  }
}
treeplot <- function(){
  m <- matrix(c(rep(1, 20), 2:21), nrow = 4, ncol = 10)
  l <- layout(m)
  maintree()
  subtrees()
  layout(matrix(1, 1, 1))
}
treeplot()
```

![plot of chunk treecomparison](./figures/kamvar2017population//treecomparison-1.png)

```r
pdf(here::here("results/figures/treeplot.pdf"), width = 6, height = 3.5)
treeplot()
dev.off()
```

```
## quartz_off_screen 
##                 2
```

# Example populations

I can use base graphics to draw example populations:


```r
rando_points <- function(N = 20, x, y, xm, ym, spread = 25){
  rX <- range(x)
  rX <- rX + ( rX * c(-spread, spread) * 1/xm )
  rY <- range(y)
  rY <- rY + ( rY * c(-spread, spread) * 1/ym)
  X <- runif(N, rX[1], rX[2])
  Y <- runif(N, rY[1], rY[2])
  points(X, Y, pch = 21, col = c(adegenet::transp(rep("grey30", N - 1), 0.5), "black"), 
         bg = c(adegenet::transp(rep("white", N - 1), 0.5), "grey30"))
}

plot_fields <- function(N = 100, P = 20, S = 20, space = 500, spread = 25){
  res <- matrix(nrow = P * N, ncol = 4)
  rows <- floor(sqrt(P))
  columns <- ceiling(P/rows)
  counter <- 1
  for (i in space * (1:columns)){
    for (j in space * (1:rows)){
      stop <- counter + N - 1
      res[counter:stop, 1] <- rnorm(N, mean = i, sd = 10)
      res[counter:stop, 2] <- rnorm(N, mean = j, sd = 10)
      res[counter:stop, 3] <- i
      res[counter:stop, 4] <- j
      counter <- stop + 1
    }
  }
  
  pops <- res %>% 
    as_data_frame() %>%
    set_names(c("x", "y", "xm", "ym")) %>%
    add_column(group = rep(seq.int(P), each = N))
  
  islands <- pops %>%
    group_by(group) %>%
    summarize(X = mean(x), Y = mean(y), rad = sd(x))
  
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))
  par(mar = rep(0.1, 4))
  
  # creating background
  symbols(islands$X, islands$Y, 
          circles = islands$rad * 10, 
          inches = FALSE, 
          bg = adegenet::funky(P),
          asp = 1,
          bty = "n",
          xaxt = "n",
          xlab = NA,
          yaxt = "n",
          ylab = NA)
  rect(
    par("usr")[1],
    par("usr")[3],
    par("usr")[2],
    par("usr")[4],
    col = "white",
    border = NA
    )
  
  # Drawing "islands"
  symbols(islands$X, islands$Y, 
          circles = islands$rad * 10, 
          inches = FALSE, 
          bg = adegenet::funky(P),
          asp = 1,
          bty = "n",
          xaxt = "n",
          xlab = NA,
          yaxt = "n",
          ylab = NA,
          add = TRUE)

  # Drawing samples
  for (i in unique(pops$group)){
    res <- filter(pops, group == i)
    rando_points(S, res$x, res$y, res$xm, res$ym, spread = spread)
  }
}
set.seed(2017-12-05)
plot_fields(space = 500)
```

![plot of chunk fieldplot](./figures/kamvar2017population//fieldplot-1.png)


I can use plot_grid from cowplot to combine the main tree and the ridgeline density plot.


```r
# p <- plot_grid(function() plot_fields(space = 1000), maintree, p_ridge + theme(legend.position = "none"), 
#                ncol = 2,
#                rel_heights = c(0.5, 0.5, 0.5),
#                rel_widths = c(1, 0.4, 0.6),
#                labels = "AUTO")
# p
set.seed(2017-12-05)
p <-  cowplot::ggdraw(xlim = c(0, 1), ylim = c(-0.25, 1)) +
  cowplot::draw_plot(function() plot_fields(N = 1000, space = 500, spread = 70), 
                     x = 0.125, 
                     y = 0.33, 
                     height = 0.66, 
                     width = 0.75) +
  cowplot::draw_plot(p_ridge + theme(legend.position = "none") + theme(asp = 1), 
                     x = 0.5, 
                     y = -0.25, 
                     height = 0.33 + 0.25, 
                     width = 0.6) +
  cowplot::draw_plot(maintree, 
                     x = 0.025, 
                     y = -0.225, 
                     width = 0.5, 
                     height = 0.33 + 0.25) +
  cowplot::draw_plot_label(c("A", "B", "C"), 
                           x = c(0.025, 0.025, 0.55), 
                           y = c(0.975, 0.35, 0.35), size = 25)
```

```
## Warning in par(inlinePars): "asp" is not a graphical parameter
```

```
## Warning in par(opar): "asp" is not a graphical parameter
```

```
## Warning in par(inlinePars): "asp" is not a graphical parameter
```

```
## Warning in par(opar): "asp" is not a graphical parameter
```

```
## Picking joint bandwidth of 0.0088
```

```
## Warning: Removed 2 rows containing non-finite values (stat_density_ridges).
```

```r
p
```

![plot of chunk unnamed-chunk-8](./figures/kamvar2017population//unnamed-chunk-8-1.png)

```r
ggsave(filename = here::here("results/figures/iatree.pdf"), plot = p, width = 6, height = 6)
ggsave(filename = here::here("results/figures/iatree.png"), plot = p, width = 6, height = 6, dpi = 600)
```



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
##  date     2017-12-05
```

```
## Packages --------------------------------------------------------------------------------
```

```
##  package      * version    date       source                             
##  ade4         * 1.7-8      2017-12-05 Github (sdray/ade4@bc33825)        
##  adegenet     * 2.1.0      2017-10-12 CRAN (R 3.4.2)                     
##  ape          * 5.0        2017-10-30 CRAN (R 3.4.2)                     
##  assertthat     0.2.0      2017-04-11 CRAN (R 3.4.0)                     
##  backports      1.1.1      2017-09-25 CRAN (R 3.4.2)                     
##  base         * 3.4.2      2017-10-04 local                              
##  bindr          0.1        2016-11-13 CRAN (R 3.4.0)                     
##  bindrcpp     * 0.2        2017-06-17 CRAN (R 3.4.0)                     
##  boot           1.3-20     2017-07-30 CRAN (R 3.4.1)                     
##  broom          0.4.2      2017-02-13 CRAN (R 3.4.0)                     
##  cellranger     1.1.0      2016-07-27 CRAN (R 3.4.0)                     
##  cli            1.0.0      2017-11-22 Github (r-lib/cli@ab1c3aa)         
##  cluster        2.0.6      2017-03-16 CRAN (R 3.4.0)                     
##  coda           0.19-1     2016-12-08 CRAN (R 3.4.0)                     
##  codetools      0.2-15     2016-10-05 CRAN (R 3.4.0)                     
##  colorspace     1.4-0      2017-11-23 R-Forge (R 3.4.2)                  
##  compiler       3.4.2      2017-10-04 local                              
##  cowplot      * 0.9.1      2017-12-04 Github (wilkelab/cowplot@89b795c)  
##  crayon         1.3.4      2017-09-23 Github (gaborcsardi/crayon@b5221ab)
##  datasets     * 3.4.2      2017-10-04 local                              
##  DBI            0.7        2017-06-18 CRAN (R 3.4.0)                     
##  deldir         0.1-14     2017-04-22 CRAN (R 3.4.0)                     
##  devtools       1.13.3     2017-08-02 CRAN (R 3.4.1)                     
##  digest         0.6.12     2017-01-27 CRAN (R 3.4.0)                     
##  doParallel   * 1.0.11     2017-09-28 CRAN (R 3.4.1)                     
##  dplyr        * 0.7.4      2017-09-28 CRAN (R 3.4.1)                     
##  evaluate       0.10.1     2017-06-24 CRAN (R 3.4.1)                     
##  expm           0.999-2    2017-03-29 CRAN (R 3.4.0)                     
##  ezknitr        0.6        2016-09-16 CRAN (R 3.4.0)                     
##  fastmatch      1.1-0      2017-01-28 CRAN (R 3.4.0)                     
##  forcats      * 0.2.0      2017-01-23 CRAN (R 3.4.0)                     
##  foreach      * 1.4.3      2015-10-13 CRAN (R 3.4.0)                     
##  foreign        0.8-69     2017-06-21 CRAN (R 3.4.0)                     
##  gdata          2.18.0     2017-06-06 CRAN (R 3.4.0)                     
##  ggplot2      * 2.2.1      2016-12-30 CRAN (R 3.4.0)                     
##  ggridges     * 0.4.1      2017-09-15 CRAN (R 3.4.1)                     
##  glue           1.2.0      2017-10-29 CRAN (R 3.4.2)                     
##  gmodels        2.16.2     2015-07-22 CRAN (R 3.4.0)                     
##  graphics     * 3.4.2      2017-10-04 local                              
##  grDevices    * 3.4.2      2017-10-04 local                              
##  grid         * 3.4.2      2017-10-04 local                              
##  gridExtra      2.3        2017-09-09 CRAN (R 3.4.1)                     
##  gridGraphics * 0.2        2017-06-06 CRAN (R 3.4.2)                     
##  gtable         0.2.0      2016-02-26 CRAN (R 3.4.0)                     
##  gtools         3.5.0      2015-05-29 CRAN (R 3.4.0)                     
##  haven          1.1.0      2017-07-09 CRAN (R 3.4.1)                     
##  here           0.1        2017-05-28 CRAN (R 3.4.0)                     
##  highr          0.6        2016-05-09 CRAN (R 3.4.0)                     
##  hms            0.3        2016-11-22 CRAN (R 3.4.0)                     
##  htmltools      0.3.6      2017-04-28 CRAN (R 3.4.0)                     
##  httpuv         1.3.5      2017-07-04 CRAN (R 3.4.1)                     
##  httr           1.3.1      2017-08-20 cran (@1.3.1)                      
##  igraph         1.1.2      2017-07-21 cran (@1.1.2)                      
##  iterators    * 1.0.8      2015-10-13 CRAN (R 3.4.0)                     
##  jsonlite       1.5        2017-06-01 CRAN (R 3.4.0)                     
##  knitr          1.17       2017-08-10 cran (@1.17)                       
##  kop          * 0.0.0.9000 2017-12-04 local (@0.0.0.9)                   
##  labeling       0.3        2014-08-23 CRAN (R 3.4.0)                     
##  lattice        0.20-35    2017-03-25 CRAN (R 3.4.0)                     
##  lazyeval       0.2.1      2017-10-29 CRAN (R 3.4.2)                     
##  LearnBayes     2.15       2014-05-29 CRAN (R 3.4.0)                     
##  lubridate      1.7.1      2017-11-03 CRAN (R 3.4.2)                     
##  magrittr       1.5        2014-11-22 CRAN (R 3.4.0)                     
##  MASS           7.3-47     2017-04-21 CRAN (R 3.4.0)                     
##  Matrix         1.2-11     2017-08-16 CRAN (R 3.4.1)                     
##  memoise        1.1.0      2017-04-21 CRAN (R 3.4.0)                     
##  methods      * 3.4.2      2017-10-04 local                              
##  mgcv           1.8-22     2017-09-19 CRAN (R 3.4.2)                     
##  mime           0.5        2016-07-07 CRAN (R 3.4.0)                     
##  mnormt         1.5-5      2016-10-15 CRAN (R 3.4.0)                     
##  modelr         0.1.1      2017-07-24 CRAN (R 3.4.1)                     
##  munsell        0.4.3      2016-02-13 CRAN (R 3.4.0)                     
##  nlme           3.1-131    2017-02-06 CRAN (R 3.4.0)                     
##  parallel     * 3.4.2      2017-10-04 local                              
##  pegas          0.10       2017-05-03 CRAN (R 3.4.0)                     
##  permute        0.9-4      2016-09-09 CRAN (R 3.4.0)                     
##  phangorn       2.3.1      2017-11-01 CRAN (R 3.4.2)                     
##  pkgconfig      2.0.1      2017-03-21 CRAN (R 3.4.0)                     
##  plyr           1.8.4      2016-06-08 CRAN (R 3.4.0)                     
##  poppr        * 2.5.0      2017-09-11 CRAN (R 3.4.1)                     
##  psych          1.7.8      2017-09-09 CRAN (R 3.4.1)                     
##  purrr        * 0.2.4      2017-10-18 cran (@0.2.4)                      
##  quadprog       1.5-5      2013-04-17 CRAN (R 3.4.0)                     
##  R.methodsS3    1.7.1      2016-02-16 CRAN (R 3.4.0)                     
##  R.oo           1.21.0     2016-11-01 CRAN (R 3.4.0)                     
##  R.utils        2.6.0      2017-11-05 CRAN (R 3.4.2)                     
##  R6             2.2.2      2017-06-17 cran (@2.2.2)                      
##  Rcpp           0.12.13    2017-09-28 CRAN (R 3.4.2)                     
##  readr        * 1.1.1      2017-05-16 CRAN (R 3.4.0)                     
##  readxl         1.0.0      2017-04-18 CRAN (R 3.4.0)                     
##  reshape2       1.4.2      2016-10-22 CRAN (R 3.4.0)                     
##  rlang          0.1.4      2017-11-05 CRAN (R 3.4.2)                     
##  rprojroot      1.2        2017-01-16 CRAN (R 3.4.0)                     
##  rstudioapi     0.7        2017-09-07 CRAN (R 3.4.1)                     
##  rvest          0.3.2      2016-06-17 CRAN (R 3.4.0)                     
##  scales         0.5.0.9000 2017-08-28 Github (hadley/scales@d767915)     
##  seqinr         3.4-5      2017-08-01 CRAN (R 3.4.1)                     
##  shiny          1.0.5      2017-08-23 cran (@1.0.5)                      
##  sp             1.2-5      2017-06-29 CRAN (R 3.4.1)                     
##  spdep          0.6-15     2017-09-01 CRAN (R 3.4.1)                     
##  splines        3.4.2      2017-10-04 local                              
##  stats        * 3.4.2      2017-10-04 local                              
##  stringi        1.1.5      2017-04-07 CRAN (R 3.4.0)                     
##  stringr      * 1.2.0      2017-02-18 CRAN (R 3.4.0)                     
##  tibble       * 1.3.4      2017-08-22 cran (@1.3.4)                      
##  tidyr        * 0.7.2      2017-10-16 CRAN (R 3.4.2)                     
##  tidyverse    * 1.2.1      2017-11-14 cran (@1.2.1)                      
##  tools          3.4.2      2017-10-04 local                              
##  utils        * 3.4.2      2017-10-04 local                              
##  vegan          2.4-4      2017-08-24 cran (@2.4-4)                      
##  viridis        0.4.0      2017-03-27 CRAN (R 3.4.0)                     
##  viridisLite    0.2.0      2017-03-24 CRAN (R 3.4.0)                     
##  withr          2.1.0      2017-11-01 CRAN (R 3.4.2)                     
##  xml2           1.1.1      2017-01-24 CRAN (R 3.4.0)                     
##  xtable         1.8-2      2016-02-05 CRAN (R 3.4.0)
```

</details>

References
==========
