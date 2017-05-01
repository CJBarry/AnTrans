# Domenico 1987

#' @rdname Domenico
#'
#' @title Domenico solution
#'
#' @description
#' The Domenico 1987 analytical solution to contaminant transport in
#'  uniform groundwater flow, including longitudinal and transverse
#'  dispersion and first-order degradation.
#'
#' @param C0 source concentration
#' @param x,y,z evaluation locations
#' @param t evaluation time
#' @param vx velocity in x direction (it is assumed that velocity only has
#'  an x component); this can be divided by the retardation factor to model
#'  sorption
#' @param ax,ay,az dispersivities (units of length) in x, y and z
#'  directions
#' @param lambda first-order decay constant
#' @param sY,sZ source dimensions in y and z directions
#'
#' @return
#' Numeric result, of same form as x, y or z.  x, y and z may be recycled.
#'
#' @details
#' \code{const_Dom} returns the solution from a constant source which
#'  starts at time 0.  \code{square_Dom} returns the solution from a finite
#'  duration constant source that starts at time 0 and stops after
#'  \code{dur}.
#'
#' @note
#' Some special cases may be achieved with certain parameter sets.
#'
#' If lambda = 0 (no degradation), the Domenico-Robbins (1985) solution
#'  results.\cr
#' If ax = 0 (no longitudinal dispersion), the Domenico-Palciauskas (1982)
#'  solution results.\cr
#' If ay = az = 0 (no transverse dispersion), the Bear (1979) solution
#'  results.\cr
#' If ay = az = 0 and lambda = 0, the Ogata-Banks (1961) solution
#'  results.\cr
#'
#' @references
#' Domenico, P. A. (1987). An analytical model for multidimensional transport of a decaying contaminant species. Journal of Hydrology, 91(1), 49-58. https://doi.org/10.1016/0022-1694(87)90127-2
#'
#' @importFrom pracma erf
#' @importFrom pracma erfc
#' @export
#'
const_Dom <- function(C0 = 1, x, y, z, t, vx, ax, ay, az, lambda, sY, sZ){
  # recycle ax so that the comparisons are fully vectorised
  N <- max(length(x), length(vx), length(t))
  tmp <- double(N)
  tmp[] <- ax; ax <- tmp; rm(tmp)

  (C0/8)*
    ifelse(ax == 0, ifelse(x < vx*t, 2, 0), {
      exp((x/(2*ax))*(1 - (1 + (4*lambda*ax)/vx)^.5))*
        erfc((x - vx*t*(1 + (4*lambda*ax)/vx)^.5)/(2*(ax*vx*t)^.5))
    })*
    (erf((y + sY/2)/(2*(ay*x)^.5)) - erf((y - sY/2)/(2*(ay*x)^.5)))*
    (erf((z + sZ/2)/(2*(az*x)^.5)) - erf((z - sZ/2)/(2*(az*x)^.5)))
}

#' @rdname Domenico
#'
#' @param dur duration of pulse
#'
#' @export
#'
square_Dom <- function(C0 = 1, x, y, z, t, vx, ax, ay, az, lambda, sY, sZ,
                       dur){
  # recycle t so that the comparisons are fully vectorised
  tmp <- double(length(x))
  tmp[] <- t; t <- tmp; rm(tmp)

  ifelse(t < dur,
         const_Dom(C0, x, y, z, t, vx, ax, ay, az, lambda, sY, sZ),
         const_Dom(C0, x, y, z, t, vx, ax, ay, az, lambda, sY, sZ) -
           const_Dom(C0, x, y, z, t - dur, vx, ax, ay, az, lambda, sY, sZ))
}
