
#' Compute Intersection Points Between an Ellipse and a Line
#'
#' This function calculates the intersection points between an ellipse and a straight line
#' that passes through the center of the ellipse with a given slope. It handles both regular
#' slopes and vertical lines (slope = Inf).
#'
#' @param center A numeric vector of length 2 specifying the center of the ellipse \code{(h, k)}.
#' @param radius A numeric vector of length 2 specifying the semi-axes of the ellipse \code{(a, b)}.
#' @param lineslop A numeric value specifying the slope \code{m} of the line. Use \code{Inf} for vertical lines.
#'
#' @return A 2x2 matrix with row names \code{"p1"} and \code{"p2"}, and column names \code{"x"} and \code{"y"},
#' representing the coordinates of the two intersection points.
#'
#' @examples
#' \dontrun{
#'   # Regular slope
#'   EllipseLine_intersection(center = c(0, 0), radius = c(3, 2), lineslop = 1)
#'
#'   # Vertical line
#'   EllipseLine_intersection(center = c(0, 0), radius = c(3, 2), lineslop = Inf)
#' }
#'
#' @export

EllipseLine_intersection = function(center,radius,lineslop){

  h = center[1] #x
  k = center[2] #y
  a = radius[1] #x
  b = radius[2] #y
  m = lineslop

  if(m == Inf){

    y_offset <- b

    y1 <- k - y_offset
    y2 <- k + y_offset
    x1 <- h
    x2 <- h

  }else{

    # calculate x_offset
    denominator <- 1 / a^2 + m^2 / b^2
    x_offset <- sqrt(1 / denominator)

    # calculate intersection
    x1 <- h - x_offset
    x2 <- h + x_offset

    y1 <- m * (x1 - h) + k
    y2 <- m * (x2 - h) + k

  }

  mat = matrix(c(x1,x2,y1,y2),nrow = 2,ncol = 2,dimnames = list(c("p1","p2"),c("x","y")))

  return(mat)

}
