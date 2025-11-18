#' Plot method for quasar objects
#'
#' Produces a plot of a \code{quasar} object,
#' typically returned by the \code{\link{closedTesting}} function.
#' It shows the estimated coefficients by quantile level, highlighting
#' statistically significant coefficients based on adjusted p-values.
#'
#' @param x An object of class \code{quasar}.
#' @param alpha Significance level.
#' @param legend.position Position of the legend.
#' @param main Main plot title.
#' @param xlab Label for the x-axis.
#' @param ylab Label for the y-axis.
#' @param col.line Color of the connecting line.
#' @param col.sig Color for significant points.
#' @param col.nonsig Color for non-significant points.
#' @param pch.sig Point character for significant points.
#' @param pch.nonsig Point character for non-significant points.
#' @param show.legend Logical; whether to display a legend.
#' @param ... Additional graphical parameters passed to \code{plot()}.
#' @return A base R plot.
#' @seealso \code{\link{closedTesting}}
#' @importFrom graphics points
#' @importFrom graphics abline
#' @importFrom graphics legend
#' @author Anna Vesely
#' @export
plot.quasar <- function(
    x,
    alpha = 0.05,
    legend.position = "topright",
    main = NULL,
    xlab = "Quantile level",
    ylab = "Coefficient",
    col.line = "darkgrey",
    col.sig = "darkred",
    col.nonsig = "darkgrey",
    pch.sig = 19,
    pch.nonsig = 17,
    show.legend = TRUE,
    ...
) {
  colors <- ifelse(x$p.value.adjusted <= alpha, col.sig, col.nonsig)
  shapes <- ifelse(x$p.value.adjusted <= alpha, pch.sig, pch.nonsig)

  plot(
    x$Quantile, x$Coefficient,
    type = "b",
    pch = NA,
    col = col.line,
    xlab = xlab,
    ylab = ylab,
    main = main,
    ...
  )
  points(x$Quantile, x$Coefficient, pch = shapes, col = colors)
  abline(h = 0, lty = 2, col = "darkgrey")
  if (show.legend) {
    legend(
      legend.position,
      legend = c(paste0("p.adj <= ", alpha), paste0("p.adj > ", alpha)),
      col = c(col.sig, col.nonsig),
      pch = c(pch.sig, pch.nonsig),
      bty = "n"
    )
  }
}
