##  GPL-3 License
## Copyright (c) 2025 Alexandre Combeau

#' @title Plot Segmentation Results
#' @description
#' Visualizes the result of a segmentation algorithm by plotting the time series,
#' vertical lines at detected changepoints, and horizontal segment means.
#'
#' @param data A numeric vector representing the time series data.
#' @param changepoints An integer vector of detected changepoint indices.
#' @param title A character string for the plot title. Defaults to "Segmentation Plot".
#'
#' @return This function returns a base R plot as a side effect, with:
#' \itemize{
#'   \item the original data in blue,
#'   \item vertical dashed red lines at each changepoint,
#'   \item green horizontal lines representing the segment-wise mean.
#' }
#' It does not return a value.
#'
#' @examples
#' n <- 100
#' data <- rep(c(0, 5, 2.5, 7), each = n) + rnorm(4 * n)
#' changepoints <- c(n, 2 * n, 3 * n, length(data))
#' plot_segmentation(data, changepoints, title = "Example Segmentation")
#'
#' @export
plot_segmentation <- function(data, changepoints, title = "Segmentation Plot")
{
  plot(data, type = "o", main = title, xlab = "Time", ylab = "Values", col = "blue")
  
  abline(v = changepoints, col = "red", lwd = 2, lty = 2)
  
  cpts <- c(0, changepoints)
  for (i in seq_len(length(cpts) - 1))
  {
    segment_start <- cpts[i] + 1
    segment_end <- cpts[i + 1]
    seg_mean <- mean(data[segment_start:segment_end])
    segments(x0 = segment_start, y0 = seg_mean,
             x1 = segment_end, y1 = seg_mean,
             col = "green", lwd = 2)
  }
}
