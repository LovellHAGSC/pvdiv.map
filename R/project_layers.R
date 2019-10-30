#' @title Project geodata layers
#'
#' @description
#' \code{project_layers} Project layers for plotting
#'
#' @param proj Projection to use. Otherwise defaults to WGS84,
#' centered on the middle of the extent
#' @param geodata named list containing geographic data layers. For a
#' layer to be projected and returned, it must have a class of either
#' Raster or Spatial.
#' @param extent a numeric vector of length 4, specifying
#' the ymin, ymax, xmin and xmax bounds of the map
#' @param crop.raster logical, should rasters be cropped? Speeds
#' things up a lot.
#' @param layers vector of geodata names to use. Defaults to all.
#' @param ... Not currently in use
#'
#' @details ...
#'
#' @return A four column data.table specifying the
#' original and buffered point locations
#'
#' @examples
#' \dontrun{
#' data("geodata", package = "pvdiv.map")
#' proj.layers <- project_layers(geodata = geodata)
#' projection <- proj.layers$projection
#' proj.layers <- proj.layers$proj.layers
#' plot(proj.layers[[1]])
#' }
#' @import data.table
#' @importFrom sp CRS spTransform
#' @importFrom raster projectRaster crop
#' @export
project_layers <- function(
  proj = NULL,
  extent = c(-108.53000, -64.51618, 20.76168, 52.53620),
  geodata,
  crop.raster = TRUE,
  layers = "all"){

  if (is.null(proj)) {
    if (is.null(extent))
      stop("if proj is not supplied, must give extent\n")
    proj.parse <- paste0(
      "+proj=poly +lat_0=",
      mean(extent[3:4]),
      " +lon_0=",
      mean(extent[1:2]),
      " +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
    proj <- CRS(proj.parse)
  }

  if (layers == "all")
    layers <- names(geodata)

  geodata <- geodata[layers]
  geonames <- sapply(geodata, class)
  raster.data <- geodata[grep("Raster", geonames)]
  if (length(raster.data) > 0) {
    if (crop.raster)
      raster.data <- lapply(raster.data, crop, extent)
    proj.rasterdata <- lapply(raster.data, projectRaster, crs = proj)
  }else{
    proj.rasterdata <- NULL
  }

  sp.data <- geodata[grep("^Spatial", geonames)]
  if (length(sp.data) > 0) {
    proj.spdata <- lapply(sp.data, spTransform, proj)
  }else{
    proj.spdata <- NULL
  }

  proj.out <- c(proj.spdata, proj)
  return(list(projection = proj,
              proj.layers = proj.out))
}



