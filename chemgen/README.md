-   [Introduction](#introduction)
-   [Look at raw data](#look-at-raw-data)
-   [Data analysis](#data-analysis)
-   [Step 1 | Adjust outer rows/cols](#step-1-adjust-outer-rowscols)
-   [Step 2 | Calculate fitness](#step-2-calculate-fitness)
-   [Step 3 | Change per condition](#step-3-change-per-condition)
    -   [Opacity](#opacity)
    -   [Size](#size)

# Introduction

Analysis for nitroxoline chemical genetics.

# Look at raw data

![](output/fig/unnamed-chunk-1-1.png)

The effect of the conditions (compared to LB) was small.

# Data analysis

Working in parallel with both opacity (which we prefer) and colony size:

1.  Bring the median of two outermost columns/rows multiplicatively to
    plate median
2.  Calculate the mutant fitness as a fraction of the plate median
3.  Estimate the multiplicative change in mutant fitness per condition
    using empirical Bayes’ moderated T-test

# Step 1 | Adjust outer rows/cols

![](output/fig/unnamed-chunk-3-1.png)

As the adjustment was multiplicative, larger deviations from plate
median appear corrected “more strongly” in linear scale, i.e. were
y-axis in log scale, dark and light grey points would be parallel.

# Step 2 | Calculate fitness

![](output/fig/unnamed-chunk-5-1.png)

There’s some interesting subset with lower opacity and bigger size than
on average.

# Step 3 | Change per condition

In the following, names are given to genes of mutants with more than
2-fold and significantly changed fitness.

One could speculate that the somewhat stronger effect of a drug, as for
nitro2 vs nitro1 was beneficial to elucidate the phenotype.

## Opacity

Find the results for each condition as a csv file [here](output/csv).

![](output/fig/unnamed-chunk-6-1.png)

## Size

![](output/fig/unnamed-chunk-7-1.png)
