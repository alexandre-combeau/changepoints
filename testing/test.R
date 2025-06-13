#install.packages("devtools")
#devtools::install_github("alexandre-combeau/changepoints")
#devtools::install_github("vrunge/dust")
library(changepoints)
library(dust)
gap <- 10
#gap <- 1
data <- dataGenerator_1D(chpts = c(50,100,200,300),
                         parameters = c(0,gap,0,gap),
                         sdNoise = 0)
plot(data, type = 'b')
n <- 300
resOP <- optimal_partitioning(data, beta = 2*log(n))
resDUST <- dust.1D(data, penalty = log(n))

resOP$changepoints ### OK
resDUST$changepoints ### OK
(resOP$costQ - 2*(resDUST$costQ + cumsum(data^2)/2)) ### OK


# Example of SVP in OP mode
resOP <- optimal_partitioning(data, beta = 2*log(n))
resSVP_OP <- smallest_valid_partitioning_op_mode(data, gamma = 2*log(n))
resOP$changepoints # OK
resSVP_OP$changepoints # OK
resOP$costQ # OK
resSVP_OP$costQ # OK

##########################################################################

resSVP <- smallest_valid_partitioning(data, gamma = 2*log(n), valid_test = validity_function1)
resSVP_OP$costQ
resOP$costQ
resSVP$costQ #### WHAT ??? WITH gap <- 10

resOP$costQ - resSVP_OP$costQ
resOP$costQ[n]
resSVP$costQ[n]

#####################

# FAIS UN PETIT EXEMPLE AVEC DES PRINTS !!!

data <- c(rep(5,5), rep(15,5))
plot(data, type = 'b')
n <- length(data)

resSVP <- smallest_valid_partitioning(data, gamma = 2*log(n), valid_test = validity_function_range) # does not work
resSVP <- smallest_valid_partitioning(data, gamma = 2*log(n), valid_test = validity_function2) # cost problem
resSVP$changepoints
resSVP$costQ


#####################