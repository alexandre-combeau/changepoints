# Changepoints detection in a time series

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

## Runtime of the OP algorithm without PELT pruning
<!-- ![Execution_Time_OP](https://github.com/user-attachments/assets/2a1f5e1a-ee14-4781-beba-1a59ea38fd7f) -->
<!-- ![OP_time_execution](https://github.com/user-attachments/assets/99c0535a-0245-4262-adcc-1f66ef8d95c1) -->
![OP_time_execution](https://raw.githubusercontent.com/alexandre-combeau/changepoints/refs/heads/main/simulations/OP_time_execution2.svg)

## Runtime of the OP algorithm with PELT pruning
<!-- ![Execution_Time_PELT](https://github.com/user-attachments/assets/b435a7f7-f7b9-47ce-9268-04037abaac08) -->
![PELT_time_execution](https://raw.githubusercontent.com/alexandre-combeau/changepoints/refs/heads/main/simulations/PELT_time_execution.svg)

## Runtime SVP / PELT
![SVP_PELT](https://raw.githubusercontent.com/alexandre-combeau/changepoints/refs/heads/main/simulations/SVP_PELT.svg)
