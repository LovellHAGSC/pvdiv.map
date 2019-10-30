#' @title Move points that are too close
#'
#' @description
#' \code{buffer_ovlpts} Buffer overlapping points along a
#' grid.
#'
#' @param spacing numeric of length 1, specifying minimum distance
#' (in units of bounds/pt.coords) that points should be apart.
#' @param pt.buffer numeric of length 1, specifying minimum distance
#' (in units of bounds/pt.coords) that points should be away from true
#' locations
#' @param bounds either object of class extent, or a numeric vector of
#' length 4, specifying the ymin, ymax, xmin and xmax bounds of the
#' points.
#' @param pt.coords Coordinates of points (x,y)
#' @param ... Not currently in use
#'
#' @details ...
#'
#' @return A four column data.table specifying the
#' original and buffered point locations
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @importFrom dbscan frNN
#' @importFrom raster ymin ymax xmin xmax pointDistance
#' @export
buffer_ovlpts <- function(spacing,
                          pt.buffer = NULL,
                          bounds,
                          pt.coords){
  spacey <- spacing
  spacex <- spacing/2

  ############################################
  # 1. Parse the plot extent
  if (class(bounds) == "Extent") {
    ymin = ymin(bounds)
    ymax = ymax(bounds)
    xmin = xmin(bounds)
    xmax = xmax(bounds)
  }else{
    if (length(bounds) != 4)
      stop("bounds must either be an object of class extent, or a numeric vector of length 4 (ymin, ymax, xmin, xmax)")
    ymin = bounds[1]
    ymax = bounds[2]
    xmin = bounds[3]
    xmax = bounds[4]
  }

  ############################################
  # 2. Build the grid coordinates
  ygrd1 <- seq(
    from = ymin,
    to = ymax,
    by = spacey)
  ygrd2 <- seq(
    from = ymin+spacex,
    to = ymax,
    by = spacey)
  xgrd1 <- seq(
    from = xmin,
    to = xmax,
    by = spacex)
  xgrd2 <- seq(
    from = xmin + spacex/2,
    to = xmax,
    by = spacex)

  nx <- min(c(
    length(xgrd1),
    length(xgrd2)))
  ny <- min(c(
    length(ygrd1),
    length(ygrd2)))

  ############################################
  # 3. Build the grid
  grd1 <- rbindlist(lapply(1:nx, function(i){
    rbind(
      data.table(x = xgrd1[i],
                 y = ygrd1[1:ny]),
      data.table(x = xgrd2[i],
                 y = ygrd2[1:ny]))
  }))

  ############################################
  # 4. Cull grid near points, if desired
  pts.xy <- data.table(pt.coords)
  setnames(pts.xy, c("x","y"))

  if (!is.null(pt.buffer)) {
    which.is.pts <- 1:nrow(pts.xy)
    grd.pts <- rbind(
      pts.xy,
      grd1)

    nn <- frNN(
      x = grd.pts,
      eps = pt.buffer)

    near <- unique(unlist(nn$id[which.is.pts]))
    near <- near[order(near)]
    near <- near[!near %in% which.is.pts]
    near <- near - nrow(pts.xy)

    if (length(near) == 0) {
      grd2 <- grd1
    }else{
      grd2 <- grd1[-near]
    }
  }else{
    grd2 <- grd1
  }

  ############################################
  # 5. Calculate point-grid distances
  ds <- apply(pts.xy, 1, function(x)
    pointDistance(
      x,
      grd2,
      lonlat = F))

  colnames(ds) <- 1:nrow(pts.xy)
  rownames(ds) <- 1:nrow(grd2)
  dt <- cbind(
    ds,
    rep(max(ds) + 1, nrow(ds)))

  ############################################
  # 6. Pull out the closest distances for each
  i <- ncol(dt)
  bests <- NULL
  while(i >= 3){
    i <- ncol(dt)
    wh <- which(dt == min(dt), arr.ind = T)
    if(length(wh) > 2){
      nr <- floor(mean(1:nrow(wh)))
      wh <- wh[nr,]
    }
    who <-  data.table(
      which.grd = rownames(dt)[wh[1]],
      which.pt = colnames(dt)[wh[2]])
    if(is.null(bests)){
      bests = who
    }else{
      bests <- rbind(bests,who)
    }
    dt <- dt[-wh[1],-wh[2]]
  }


  bests[,which.grd := as.numeric(which.grd)]
  bests[,which.pt := as.numeric(which.pt)]
  setkey(bests, which.pt)

  out <- data.table(
    true.x = pts.xy$x,
    true.y = pts.xy$y,
    plot.x = grd2$x[bests$which.grd],
    plot.y = grd2$y[bests$which.grd]
  )
  return(out)
}
