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


![Execution_Time_OP](https://github.com/user-attachments/assets/63042279-f82c-4b8e-ab38-b752cf738d66)
