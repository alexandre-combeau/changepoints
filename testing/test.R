library(dust)
library(changepoints)

n <- 100
x <- dataGenerator_1D(chpts = c(2, 4, 6) * n,
                      parameters = c(5,10, 2),
                      sdNoise = 1,
                      type = "gauss")
best_beta <- 2 * sdDiff(x, "HALL")^2 * log(length(x))
OP_CPP <- rcpp_optimal_partitioning(x, best_beta)

# Visualization
cpts <- OP_CPP$changepoints
plot(x, type = "o", main = "Optimal Partitioning avec R",
     ylab = "Values", xlab = "Time", col = "blue")
abline(v = cpts, col = "red", lwd = 2, lty = 2)
cpts2 <- c(0,cpts)
for (i in 1:(length(cpts2) - 1))
{
  segment_mean <- mean(x[(cpts2[i] + 1):cpts2[i + 1]])
  segments(x0 = cpts2[i] + 1, y0 = segment_mean,
           x1 = cpts2[i+1], y1 = segment_mean,
           col = "green", lwd = 2)
}
