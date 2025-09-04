# Changepoints detection in a time series

A central goal in time-series analysis is to identify the sequence of structural changes underlying the observed data.
This problem has a long-standing history, with major developments since the mid-twentieth century.

This problem presents two distinct challenges :
  - **SINGLE** changepoint detection, which involves identifying a change within an online data stream
  - **MULTIPLE** changepoints detection, which seeks to partition a fixed-length time series into homogeneous segments.

For clarity, we focus on the canonical case of detecting changes in the **MEAN** of the data distribution.

Dynamic programming plays a central role in multiple change-point detection, as it provides the optimal segmentation by exactly minimizing a global cost function – typically formulated as a penalized likelihood

## Installation via the `devtools` package :

```r
install.packages("devtools")
devtools::install_github("alexandre-combeau/changepoints")
```

## Temps d'exécution pour OP sans Pruning
![Execution_Time_OP](https://github.com/user-attachments/assets/2a1f5e1a-ee14-4781-beba-1a59ea38fd7f)

## Temps d'exécution pour OP avec Pruning
![Execution_Time_PELT](https://github.com/user-attachments/assets/b435a7f7-f7b9-47ce-9268-04037abaac08)

## Temps
![test](https://github.com/user-attachments/assets/c190bcfb-d4a2-4658-bca6-b6af1eb8cd66)
