#     valueMeasure.R Bias correction methods
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
#' @param y Grid (also station data) of observations
#' @param x Grid (also station data) of the grid that is being validated
#' @param measure.code characher of the measure code to be computed (use VALUE::show.measures)
#' @param index.code Default is NULL. characher of the index code to be computed (use VALUE::show.indices). 
#' @param return.NApercentage Logical to also return or not a grid containing NA percentage information.
#' @template templateParallelParams 
#' @details Some measures are computed directly from the original time series (e.g. temporal correlation),
#' whereas others are computed upon previouly computed indices (e.g. mean bias). 
#' Thus, argument \code{index.code} must be provided for the latter case.
#' @return A grid of the index or a list containing the grid of the index and the 
#' grid of NA percenatage
#' @importFrom abind abind adrop
#' @import transformeR
#' @importFrom VALUE valueMeasure1D valueIndex1D
#' @author M. Iturbide
#' @export
#' @examples 
# library(transformeR)
# y <- EOBS_Iberia_tas
# x <- CFS_Iberia_tas
# bias <- valueMeasure(y, x, measure.code = "bias", index.code = "mean")
# str(bias$Measure)
# str(bias$NAmeanPercentage)

valueMeasure <- function(y, x, 
                         measure.code, 
                         index.code = NULL,
                         return.NApercentage = TRUE,
                         parallel = FALSE,
                         max.ncores = 16,
                         ncores = NULL){
  if (any(c("biasCirc", "bias", "biasRel", "ratio") %in% measure.code)) {
    if (is.null(index.code)) stop("This measure requires previous index calculation. Please provide an index.code (Use VALUE::show.indices() to choose an index).")
    index.y <- valueIndex(grid = y, index.code = index.code, parallel = parallel, max.ncores = max.ncores, ncores = ncores)
    index.x <- valueIndex(grid = x, index.code = index.code, parallel = parallel, max.ncores = max.ncores, ncores = ncores)
    y <- index.y[["Index"]]
    x <- index.x[["Index"]]
    na.percent <- index.y[["NApercentage"]]
    rm(index.y, index.x)
  }
  station <- FALSE
  if ("loc" %in% getDim(y)) station <- TRUE
  xy <- y$xyCoords
  dimNames <- attr(redim(y, member = FALSE, loc = station)$Data, "dimensions")
  suppressWarnings(suppressMessages(ix <- interpGrid(x, getGrid(y))))
  y <- redim(y, drop = TRUE)
  y <- redim(y, member = FALSE, runtime = FALSE)
  ix <- redim(ix, member = TRUE, runtime = TRUE)
  out <- ix
  n.run <- getShape(ix)["runtime"]
  n.mem <- getShape(ix)["member"]
  runarr <- lapply(1:n.run, function(l) {
    memarr <- lapply(1:n.mem, function(m) {
      message("[", Sys.time(), "] Computing member ", m, " out of ", n.mem)
      o = y$Data[, , , drop = FALSE]
      p = adrop(ix$Data[l, m, , , , drop = FALSE], drop = c(T, T, F, F, F))
      data <- list(o, p)
      if (!station) {
        data <- lapply(1:length(data), function(x) {
          attr(data[[x]], "dimensions") <- c("time", "lat", "lon")
          abind(array3Dto2Dmat(data[[x]]), along = 3)
        }) 
      }
      o <- lapply(seq_len(ncol(data[[1]])), function(i) data[[1]][,i,1])
      p <- lapply(seq_len(ncol(data[[2]])), function(i) data[[2]][,i,1])
      nona.o <- lapply(o, function(x) which(!is.na(x)))
      nona.p <- lapply(p, function(x) which(!is.na(x)))
      sea <- unlist(lapply(1:length(nona.o), function(x) {
        if (length(nona.o[[x]]) == 0 | length(nona.p[[x]]) == 0)  x
      }))
      nonaind <- lapply(1:length(nona.o), function(x) {
        if (length(nona.o[[x]]) != 0 && length(nona.p[[x]]) != 0)  {
          intersect(nona.o[[x]] , nona.p[[x]])
        } else {
          1:length(y$Dates$start)
        }
      })
      dates <- list()
      for (i in 1:length(nonaind)) {
        dates[[i]] <- y$Dates$start[nonaind[[i]]]
        o[[i]] <- o[[i]][nonaind[[i]]]
        p[[i]] <- p[[i]][nonaind[[i]]]
      }
      mat <- abind(valueMeasureXD(o = o, p = p, io = o, ip = p, dates = dates, measure.code = measure.code, 
                                  parallel = parallel, max.ncores = max.ncores, ncores = ncores), 
                   along = 0)
      na.mat <- do.call("abind", 
                        list(lapply(nonaind, function(x) 100 - (length(x)/length(y$Dates$start)*100))
                             , along = 2))
      na.mat[1, sea] <- NA
      if (!station) mat <- mat2Dto3Darray(mat, xy$x, xy$y)
      if (!station) na.mat <- mat2Dto3Darray(na.mat, xy$x, xy$y)
      list(mat, na.mat)
    })
    list(unname(do.call("abind", list(lapply(memarr, "[[", 1), along = 0))),
         unname(do.call("abind", list(lapply(memarr, "[[", 2), along = 0))))
  })
  out.na <- out
  out$Data <- unname(do.call("abind", list(lapply(runarr, "[[", 1), along = 0)))
  out.na$Data <- unname(do.call("abind", list(lapply(runarr, "[[", 2), along = 0)))
  attr(out$Data, "dimensions") <- unique(c("runtime", "member", dimNames))
  attr(out.na$Data, "dimensions") <- unique(c("runtime", "member", dimNames))
  out$Dates$start <- y$Dates$start[1]
  out$Dates$end <- range(y$Dates$end)[2]
  out.na$Dates$start <- y$Dates$start[1]
  out.na$Dates$end <- range(y$Dates$end)[2]
  if (station) out <- redim(out, loc = TRUE)
  if (station) out.na <- redim(out.na, loc = TRUE)
  out.na$Variable$varName <- "NApercentage"
  out <- redim(out, drop = TRUE)
  out.na <- redim(out.na, drop = TRUE)
  message("[", Sys.time(), "] Done.")
  if (return.NApercentage) {
    if (any(c("biasCirc", "bias", "biasRel", "ratio") %in% measure.code)) { 
      return(list("Measure" = out, "NAmeanPercentage" = na.percent))
    } else {
      return(list("Measure" = out, "NApercentage" = out.na))
    }
  } else {
    return(out)
  }
}
#end

valueMeasureXD <- function(o = NULL, p = NULL, 
                           io = NULL, ip = NULL,
                           dates = NULL, measure.code = NULL, 
                           parallel = FALSE,
                           max.ncores = 16,
                           ncores = NULL){
  parallel.pars <- parallelCheck(parallel, max.ncores, ncores)
  mapply_fun <- selectPar.pplyFun(parallel.pars, .pplyFUN = "mapply")
  if (parallel.pars$hasparallel) on.exit(parallel::stopCluster(parallel.pars$cl))
  mapply_fun(valueMeasure1D, obs = o, prd = p, indexObs = io, indexPrd = ip, dates, MoreArgs = list(measure.codes = measure.code))
}
#end
