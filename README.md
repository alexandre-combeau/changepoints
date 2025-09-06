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
<!-- ![Execution_Time_OP](https://github.com/user-attachments/assets/2a1f5e1a-ee14-4781-beba-1a59ea38fd7f) -->
![OP_time_execution](https://github.com/user-attachments/assets/99c0535a-0245-4262-adcc-1f66ef8d95c1)
![test](https://github.com/user-attachments/assets/6361b77f-394d-457a-997c-efd318abedab)

## Temps d'exécution pour OP avec Pruning
<!-- ![Execution_Time_PELT](https://github.com/user-attachments/assets/b435a7f7-f7b9-47ce-9268-04037abaac08) -->
![PELT_time_execution](https://github.com/user-attachments/assets/fb5aa17a-12d9-4741-9b76-f1c6a48b8e99)

## Temps d'exécution SVP / PELT
![svp_pelt](https://github.com/user-attachments/assets/d8fa062a-10f1-473b-a36e-44131cc50da7)
