# Détection de ruptures dans une série temporelle

## Comparaison des algorithmes PELT et Optimal Partitioning en R et Rcpp

A ce jour, ce package compare deux algorithmes de segmentation de séries temporelles :
  - **PELT** (Pruned Exact Linear Time)
  - **Optimal Partitioning**

implémentés à la fois en **R** et en **Rcpp**. L'objectif est d'évaluer leurs performances en termes de **vitesse d'exécution**, de **précision des ruptures détectées**, et de **complexité algorithmique**.

Equation :

$Q_t=\min_1\{Q_s+C(y_{s..t})+\beta\}$

## Installation avec `devtools` :
```r
install.packages("devtools")
devtools::install_github("alexandre-combeau/changepoints")
```
## Temps d'exécution entre R et Rcpp
![Execution_Time_OP](https://github.com/user-attachments/assets/2a1f5e1a-ee14-4781-beba-1a59ea38fd7f)
