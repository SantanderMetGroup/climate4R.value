#     valueIndex.R Bias correction methods
#
#     Copyright (C) 2017 Santander Meteorology Group (http://www.meteo.unican.es)
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.

#' @title VALUE measure calculation for climate4R grids
#' @description VALUE measure calculation for climate4R grids
#' @param grid Grid (also station data) of observations
#' @param index.code Character of the index code to be computed (use VALUE::show.indices). 
#' @param return.NApercentage Logical to also return or not a grid containing NA percentage information.
#' @param ... Additional arguments pased to FUNCTION.CALL (use VALUE::show.indices to check the available arguments).
#' @param condition Inequality operator to be applied to the given \code{"threshold"}. Only the days that satisfy the condition will be used for validation the model.
#' \code{"GT"} = greater than the value of \code{threshold}, \code{"GE"} = greater or equal,
#' \code{"LT"} = lower than, \code{"LE"} = lower or equal than.
#' @param threshold Numeric value. Threshold used as reference for the condition. Default is NULL. If a threshold value is supplied with no specificaction of the argument \code{condition}. Then condition is set to \code{"GE"}.
#' @param which.wetdays A string, default to NULL. Infer the measure/index taking into account only the wet days of the temporal serie. 
#' In this case, set which.wetdays = "Independent".
#' @template templateParallelParams 
#' @return A grid of the index or a list containing the grid of the index and the 
#' grid of NA percenatage
#' @importFrom abind adrop abind
#' @import transformeR
#' @importFrom VALUE valueIndex1D show.indices
#' @author M. Iturbide
#' @author M. N. Legasa
#' @export
#' @examples 
#' library(transformeR)
#' require(climate4R.datasets)
#' data(EOBS_Iberia_tas)
#' y <- EOBS_Iberia_tas
#' m <- valueIndex(y, "mean")
#' str(m$Index)
#' str(m$NApercentage)

valueIndex <- function(grid = NULL, index.code = NULL,
                       return.NApercentage = FALSE,
                       ...,
                       parallel = FALSE,
                       max.ncores = 16,
                       condition = NULL, 
                       threshold = NULL,
                       which.wetdays = NULL,
                       ncores = NULL){
  index.arg.list <- list(...)
  if (is.null(index.code) || length(index.code) > 1) stop("Please provide a single index.code.")
  if (!is.null(threshold) & is.null(condition)) condition = "GE"
  if (!is.null(threshold) & !is.null(condition) & is.null(which.wetdays)) stop("Please select the wet days subset with the which.wetdays parameter.")
  if (!is.null(which.wetdays)) { 
    if (which.wetdays != "Independent") stop("The only valid which.wetdays value is 'Independent' ")
  }
  if (!is.null(condition)) {
    if (is.null(threshold)) stop("Please specify the threshold value with parameter 'threshold'")
    ineq <- switch(condition,
                   "GT" = ">",
                   "GE" = ">=",
                   "LT" = "<",
                   "LE" = "<=")
  }
  station <- FALSE
  if ("loc" %in% getDim(grid)) station <- TRUE
  xy <- grid$xyCoords
  dimNames <- attr(grid$Data, "dimensions")
  grid <- redim(grid, drop = TRUE)
  grid <- redim(grid, member = TRUE, runtime = TRUE)
  out <- grid
  n.run <- getShape(grid)["runtime"]
  n.mem <- getShape(grid)["member"]
  runarr <- lapply(1:n.run, function(l) {
    memarr <- lapply(1:n.mem, function(m) {
      message("[", Sys.time(), "] Computing member ", m, " out of ", n.mem)
      p = adrop(grid$Data[l, m, , , , drop = FALSE], drop = c(T, T, F, F, F))
      if (!station) {
          attr(p, "dimensions") <- c("time", "lat", "lon")
          p <- abind(array3Dto2Dmat(p), along = 3)
      }
      p <- lapply(seq_len(ncol(p)), function(i) {
        yy_p <- p[,i,1]
        if (is.null(which.wetdays)) {
          yy_p
        } else if (which.wetdays >= "Independent") {
          ind_p = eval(parse(text = paste("yy_p", ineq, "threshold")))
          yy_p[ind_p]
        }
      })
      nona.p <- lapply(p, function(x) which(!is.na(x)))
      sea <- unlist(lapply(1:length(nona.p), function(x) {
        if (length(nona.p[[x]]) == 0)  x
      }))
      nonaind <- lapply(1:length(nona.p), function(x) {
        if (length(nona.p[[x]]) != 0)  {
          nona.p[[x]]
        } else {
          1:length(grid$Dates$start)
        }
      })
      dates <- list()
      
      # NAs are handled internally for spatial indices
      # if (!show.indices(index.code = index.code)[1,"SPATIAL"]){
      #   for (i in 1:length(nonaind)) { 
      #     dates[[i]] <- grid$Dates$start[nonaind[[i]]]
      #     p[[i]] <- p[[i]][nonaind[[i]]]
      #   }
      # } else {
        dates <- grid$Dates$start
      # }
      
      na.mat <- do.call("abind", 
                        list(lapply(nonaind,
                                    function(x) 
                                      100 - (length(x)/length(grid$Dates$start)*100))
                             , along = 2)
                        )
      na.mat[1, sea] <- NA
      args <- c(list(ts = p, dates = dates, 
                   index.code = index.code, 
                   parallel = parallel,
                   max.ncores = max.ncores, 
                   ncores = ncores), index.arg.list)
      idx <- do.call("valueIndexXD", args)
      if (is.matrix(idx)){
        if (is.null(attributes(idx)$dimnames.names)) attr(idx, "dimnames.names") <- c("loc", "loc")
        tbr <- list(idx = idx, na.mat = na.mat)
      }
      else{
        idx <- abind(idx, along = 0)
        if (!station) idx <- mat2Dto3Darray(idx, xy$x, xy$y)
        if (!station) na.mat <- mat2Dto3Darray(na.mat, xy$x, xy$y)
        tbr <- list(idx = idx, na.mat = na.mat)
      }
      return(tbr)
    })
    tbr <-
      list(idx = unname(do.call("abind", list(lapply(memarr, "[[", 1), along = 0))),
           na = unname(do.call("abind", list(lapply(memarr, "[[", 2), along = 0))))
    if (!is.null(attributes(memarr[[1]]$idx)$dimnames.names)){
      attr(tbr, "dimnames.names") <- 
        attributes(memarr[[1]]$idx)$dimnames.names
    }
    return(tbr)
  })

  if (return.NApercentage) out.na <- out
  
  # if (show.indices(index.code = index.code)[1,"SPATIAL"]){
  #   dimnames.names <- attributes(runarr[[1]])$dimnames.names
  #   out$Data <- unname(do.call("abind", list(lapply(runarr, "[[", 1), along = 0)))
  #   
  #   attr(out$Data, "dimensions") <- 
  #     c(c("runtime", "member"), dimnames.names)
  #   out$Dates$start <- grid$Dates$start[1]
  #   out$Dates$end <- range(grid$Dates$end)[2]
  #   out$Variable$varName <- index.code
  #   
  #   if (!station) {
  #     xyCoords.attr <- attributes(out$xyCoords)
  #     out$xyCoords <- transformeR::get2DmatCoordinates(out)
  #     attr(out$xyCoords, "resX.grid") <- xyCoords.attr$resX
  #     attr(out$xyCoords, "resX") <- 0
  #     attr(out$xyCoords, "resY.grid") <- xyCoords.attr$resY
  #     attr(out$xyCoords, "resY") <- 0
  #     dimnames(out$Data)[[length(dim(out$Data))]] <-
  #       mapply(FUN = paste, out$xyCoords$x, out$xyCoords$y,
  #              MoreArgs = list(sep = ",")
  #     )
  #     dimnames(out$Data)[[length(dim(out$Data))-1]] <- 
  #       mapply(FUN = paste, out$xyCoords$x, out$xyCoords$y,
  #              MoreArgs = list(sep = ",")
  #     )
  #   } else {
  #     dimnames(out$Data)[[length(dim(out$Data))]] <- 
  #       grid$Metadata$station_id
  #     dimnames(out$Data)[[length(dim(out$Data))-1]] <- 
  #       grid$Metadata$station_id
  #   }
  #       
  #   out <- redim(out, drop = TRUE)
  # } else {
    out$Data <- unname(do.call("abind", list(lapply(runarr, "[[", 1), along = 0)))
    attr(out$Data, "dimensions") <- unique(c("runtime", "member", dimNames))
    out$Dates$start <- grid$Dates$start[1]
    out$Dates$end <- range(grid$Dates$end)[2]

    if (station) out <- redim(out, loc = TRUE)
    out$Variable$varName <- index.code
    out <- redim(out, drop = TRUE)
  # }
  message("[", Sys.time(), "] Done.")
  
  if (return.NApercentage) {
    out.na$Data <- unname(do.call("abind", list(lapply(runarr, "[[", 2), along = 0)))
    attr(out.na$Data, "dimensions") <- unique(c("runtime", "member", dimNames))
    out.na$Dates$start <- grid$Dates$start[1]
    out.na$Dates$end <- range(grid$Dates$end)[2]
    if (station) out.na <- redim(out.na, loc = TRUE)
    out.na$Variable$varName <- "NApercentage"
    out.na <- redim(out.na, drop = TRUE)
    
    return(list("Index" = out, "NApercentage" = out.na))
  } else {
    return(out)
  }
}
#end


#' @title VALUE atomic and spatial index calculation
#' @description VALUE index calculation helper function
#' @param ts list containing time series for each gridpoint/station
#' @param dates A character (or \code{POSIXct}) vector following the format \dQuote{YYYY-MM-DD}
#'  (i.e., \code{format = "\%Y-\%m-\%d"} as in \code{\link{strptime}}). Note that the dates are not required 
#'  by all indices, so the default is \code{NULL}. The function will yield an error if dates are required but not supplied.
#' @param index.code Characher of the index code to be computed (use VALUE::show.indices). 
#' @param ... Additional arguments pased to FUNCTION.CALL (use VALUE::show.indices to check the available arguments).
#' @param parallel Use parallel computing?
#' @param max.ncores Maximum number of cores to use
#' @param ncores Number of cores to use for parallel computing. Default to ncores-1.
#' @return A grid of the index or a list containing the grid of the index and the 
#' grid of NA percenatage
#' @importFrom transformeR parallelCheck
#' @importFrom VALUE valueIndex1D show.indices
#' @author M. Iturbide, M.N. Legasa

valueIndexXD <- function(ts = NULL,
                         dates = NULL, 
                         index.code = NULL, 
                         ...,
                         parallel = FALSE,
                         max.ncores = 16,
                         ncores = NULL){
  index.arg.list <- list(...)
  # if (show.indices(index.code = index.code)[1,"SPATIAL"]){
  #   # parallelization is carried out internally
  #   valueIndex2D(tsl = ts, dates = dates, index.code = index.code,
  #                parallel = parallel, max.ncores = max.ncores, ncores = ncores)
  # } else {
    parallel.pars <- parallelCheck(parallel, max.ncores, ncores)
    mapply_fun <- selectPar.pplyFun(parallel.pars, .pplyFUN = "mapply")
    if (parallel.pars$hasparallel) on.exit(parallel::stopCluster(parallel.pars$cl))
    mapply_fun(valueIndex1D, ts, dates, MoreArgs = c(list(index.codes = index.code), index.arg.list))
  # }
}
#end
