# changepoints

## Comparaison des algorithmes PELT et Optimal Partitioning en R et Rcpp

A ce jour, ce package compare deux algorithmes de segmentation de séries temporelles :
  - **PELT** (Pruned Exact Linear Time)
  - **Optimal Partitioning**

implémentés à la fois en **R** et en **Rcpp**. L'objectif est d'évaluer leurs performances en termes de **vitesse d'exécution**, de **précision des ruptures détectées**, et de **complexité algorithmique**.

Equation : $R_n$

## Installation

### 1. Avec `devtools` :
```r
install.packages("devtools")
devtools::install_github("alexandre-combeau/changepoints")
```



![Execution_Time_OP](https://github.com/user-attachments/assets/4ceb774c-ae42-43b0-be17-b45d26d1642f)
