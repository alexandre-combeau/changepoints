# Change-points detection in a time series

A central goal in time-series analysis is to identify the sequence of structural changes underlying the observed data.
This problem has a long-standing history, with major developments since the mid-twentieth century.

This problem presents two distinct challenges :
  - **SINGLE** changepoint detection, which involves identifying a change within an online data stream
  - **MULTIPLE** changepoints detection, which seeks to partition a fixed-length time series into homogeneous segments.

For clarity, we focus on the canonical case of detecting changes in the **MEAN** of the data distribution.

Dynamic programming plays a central role in multiple change-point detection, as it provides the optimal segmentation by exactly minimizing a global cost function – typically formulated as a penalized likelihood
  - $R_t=\{R_s + C,\ 0\le s<t\}$

## Installation via the `devtools` package :

```r
install.packages("devtools")
devtools::install_github("alexandre-combeau/changepoints")
```

| Runtime of OP | Runtime of PELT |
|--------------|------------------|
| ![OP_time_execution2](https://raw.githubusercontent.com/alexandre-combeau/changepoints/refs/heads/main/simulations/OP_time_execution2.svg) | ![PELT_time_execution](https://raw.githubusercontent.com/alexandre-combeau/changepoints/refs/heads/main/simulations/PELT_time_execution.svg) |

## Runtime SVP / PELT with a fixed changepoint density (0.6%)
![SVP_PELT](https://raw.githubusercontent.com/alexandre-combeau/changepoints/refs/heads/main/simulations/SVP_PELT.svg)

## Runtime SVP / PELT when data has no changepoint
![SVP_PELT_nochange](https://raw.githubusercontent.com/alexandre-combeau/changepoints/refs/heads/main/simulations/SVP_PELT_nochange.svg)

