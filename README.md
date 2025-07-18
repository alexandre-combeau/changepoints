# Changepoints detection in a time series

A central goal in time-series analysis is to identify the sequence of structural changes underlying the observed data.

This problem has a long-standing history, with major developments since the mid-twentieth century...

## Comparaison des algorithmes PELT et Optimal Partitioning en R et Rcpp

A ce jour, ce package compare deux algorithmes de segmentation de séries temporelles :
  - **PELT** (Pruned Exact Linear Time)
  - **Optimal Partitioning**

implémentés à la fois en **R** et en **Rcpp**. L'objectif est d'évaluer leurs performances en termes de **vitesse d'exécution**, de **précision des ruptures détectées**, et de **complexité algorithmique**.

Equation :

$Q_t=\min_{0\le s<t}[Q_s+C(y_{s..t})+\beta]$

## Installation avec `devtools` :

```r
install.packages("devtools")
devtools::install_github("alexandre-combeau/changepoints")
```

## Temps d'exécution pour OP sans Pruning
![Execution_Time_OP](https://github.com/user-attachments/assets/2a1f5e1a-ee14-4781-beba-1a59ea38fd7f)

## Temps d'exécution pour OP avec Pruning
![Execution_Time_PELT](https://github.com/user-attachments/assets/b435a7f7-f7b9-47ce-9268-04037abaac08)
