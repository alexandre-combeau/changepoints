# changepoints

## Comparaison des algorithmes PELT et Optimal Partitioning en R et Rcpp

A ce jour, ce package compare deux algorithmes de segmentation de séries temporelles :
  - **PELT** (Pruned Exact Linear Time)
  - **Optimal Partitioning**

implémentés à la fois en **R** et en **Rcpp**. L'objectif est d'évaluer leurs performances en termes de **vitesse d'exécution**, de **précision des ruptures détectées**, et de **complexité algorithmique**.

Equation : $R_n=\min_1\min_2\Bigg\{\Big(\sum_{k=0}^{K-1}\C(y_{\tau_{k}..\tau_{k+1}}),K\Big)\,,\, f(y_{\tau_{k}..\tau_{k+1}})\le\gamma,\ k=0,\dots,K-1\Bigg\}$

## Installation

### 1. Avec `devtools` :
```r
install.packages("devtools")
devtools::install_github("alexandre-combeau/changepoints")
```

![expected1](https://github.com/user-attachments/assets/5d756dd6-a3ba-40f0-b825-ad0b2f26d2b9)
