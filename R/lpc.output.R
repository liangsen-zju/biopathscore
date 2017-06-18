plot.lpc <- function(x, type, unscale=TRUE, lwd=1, datcol="grey60",
                     datpch=21, masscol=NULL, masspch=15, curvecol=1, splinecol=3,
                     projectcol=4, startcol=NULL, startpch=NULL,...){
  object <- x
  if (class(object)=="lpc.spline"){
    splineobject <- object
    lpcobject    <- object$lpcobject
  } else {
    lpcobject <- object
    splineobject <-NULL
  }

  if (missing(type)){
    if (class(object)=="lpc"){type="curve"}
    else if (class(object)=="lpc.spline" && object$project==FALSE){type="spline"}
    else if (class(object)=="lpc.spline" && object$project){type=c("spline","project")}
  }


  if ("project" %in% type  && class(object)=="lpc" || "project" %in% type &&  class(object)=="lpc.spline" && object$closest.branch=="none"  ){
    splineobject <-    lpc.spline(lpcobject, project=TRUE)
  } else if ("spline" %in% type  && class(object)=="lpc"){
    splineobject <-    lpc.spline(lpcobject)
  }

  n        <- dim(lpcobject$data)[1]
  d        <- dim(lpcobject$data)[2]
  branch   <- as.factor(lpcobject$P[,"branch"])
  nbranch  <- nlevels(branch)
  lc       <- length(curvecol)
  sc       <- length(splinecol)
  pc       <- length(projectcol)
  stnames  <- dimnames(lpcobject$start)[[1]]
  stdepth  <- as.numeric(stnames)
  maxdepth <- max(stdepth)


  if (is.null(masscol)){
    masscol<-as.numeric(c(2,4,1)[1:maxdepth])
  }
  if (is.null(startcol)){
    startcol<-as.numeric(c(6,5,"grey20")[1:maxdepth])
  }

  mc       <- length(masscol)
  stc      <- length(startcol)

  if (!is.null(startpch)){
    if (length(stnames)!=length(startpch)){
      stnames<-rep(startpch[1], length(stnames))
    } else {
      stnames<-startpch
    }
  }

  if (lpcobject$scaled && unscale){
    u       <- if (is.null(splineobject)) {unscale(lpcobject)} else  {unscale(splineobject)}
    Xi      <- u$data
    fit     <- u$LPC
    start   <- u$starting.points
    if ("spline" %in% type) knots.coords <- u$knots.coords
    if ("project" %in% type) closest.coords <- u$closest.coords
  } else {
    Xi    <-  lpcobject$data
    fit   <-  lpcobject$LPC
    start <-  lpcobject$starting.points
    if ("spline" %in% type) knots.coords <- splineobject$knots.coords
    if ("project" %in% type) closest.coords <- splineobject$closest.coords
  }


  if (d==2){
    # eqscplot(Xi, col=datcol, pch=datpch,...)
    plot(Xi, col=datcol, pch=datpch,...)
    if ("curve" %in% type){
      if (maxdepth==1){
        for (j in 0:(nbranch-1)){
          lines(fit[branch==j,], lwd=lwd, col= if (lc>=nbranch) curvecol[j+1] else if (lc>1) rep(curvecol,ceiling(nbranch/lc))[j+1] else curvecol)
        }
      } else {
        for (j in 0:(nbranch-1)){
          lines(fit[branch==j,], lwd=lwd, col= if (lc>=maxdepth) curvecol[stdepth[j+1]] else if (lc>1) rep(curvecol,ceiling(nbranch/lc))[j+1] else curvecol)
        }
      }
    }


    if ("mass" %in% type){
      if (mc== dim(fit)[1]){
        points(fit, col=masscol,  pch=masspch )
      } else {
        for (j in 0:(nbranch-1)){
          points(fit[branch==j,], pch=masspch, col=  if (mc==maxdepth)  masscol[as.numeric(dimnames(start)[[1]])[j+1] ]  else  if  (mc>=nbranch) masscol[j+1] else if
                 (mc>1) rep(masscol,ceiling(nbranch/mc))[j+1] else masscol)
        }
      }
    }
    if ("start" %in% type){
      points(start, col = if (stc == maxdepth) startcol[stdepth] else startcol, pch=stnames)
    }

    if ("spline" %in% type){
      for (j in 0:(nbranch-1)){
        lines(knots.coords[[j+1]][1,],knots.coords[[j+1]][2,], lwd=lwd, col= if (sc>=nbranch) splinecol[j+1] else if (sc>1) rep(splinecol,ceiling(nbranch/sc))[j+1] else splinecol)
      }
    }
    if ("project" %in% type){
      for (i in 1:n){
        x1 <-closest.coords[i,1]
        y1<- closest.coords[i,2]
        x2 <- Xi[i,1]
        y2 <- Xi[i,2]
        segments(x1,y1,x2,y2,  col= if (pc==n) projectcol[i] else if (pc>=nbranch) projectcol[splineobject$closest.branch[i]+1] else if (pc>1)  rep(projectcol,ceiling(nbranch/pc))[splineobject$closest.branch[i]+1] else projectcol  ) # 09/10/09
      }
    }



  } else if (d==3){

    #require(scatterplot3d)
    plotlpc3 <- scatterplot3d::scatterplot3d(Xi, color=datcol, pch=datpch,...)
    if ("curve" %in% type){
      for (j in 0:(nbranch-1)){
        plotlpc3$points3d(fit[branch==j,], lwd=lwd, col=if (lc>=nbranch) curvecol[j+1] else if (lc>1) rep(curvecol,ceiling(nbranch/lc))[j+1] else curvecol, type="l")
      }
    }



    if ("mass" %in% type){
      if (mc== dim(fit)[1]){
        plotlpc3$points3d(fit, col=masscol,  pch=masspch )
      } else {
        for (j in 0:(nbranch-1)){
          plotlpc3$points3d(fit[branch==j,], pch=masspch, col= if (mc==maxdepth)  masscol[as.numeric(dimnames(start)[[1]])[j+1] ]  else  if (mc>=nbranch) masscol[j+1] else if (mc>1) rep(masscol,ceiling(nbranch/mc))[j+1] else masscol, type="p")
        }
      }
    }

    if ("start" %in% type){
      plotlpc3$points3d(start, col = if (stc == maxdepth) startcol[stdepth] else startcol, pch=stnames)
    }



    if ("spline" %in% type){
      for (j in 0:(nbranch-1)){
        plotlpc3$points3d(t(knots.coords[[j+1]]), lwd=lwd, col= if (sc>=nbranch) splinecol[j+1] else if (sc>1) rep(splinecol,ceiling(nbranch/sc))[j+1] else splinecol, type="l")
      }
    }
    if ("project" %in% type){

      for (i in 1:n){
        x1 <- closest.coords[i,1]
        y1<-  closest.coords[i,2]
        z1 <- closest.coords[i,3]
        x2 <- Xi[i,1]
        y2 <- Xi[i,2]
        z2 <- Xi[i,3]
        plotlpc3$points3d(c(x1,x2),c(y1,y2),c(z1,z2), type="l",
                          col= if (pc==n) projectcol[i]  else if (pc>=nbranch) projectcol[splineobject$closest.branch[i]+1] else if (pc>1)  rep(projectcol,ceiling(nbranch/pc))[splineobject$closest.branch[i]+1] else projectcol)  # 09/10/2009
      }
    }



  }  else if (d<=16){
    pairs(fit,
          panel= if ("curve" %in%  type) "lines" else "points",
          labels= dimnames(lpcobject$data)[[2]],
          lwd=lwd,
          col=curvecol,
          ...
    )

    #X11()
    #require(lattice)
    #splom(lpcobject$LPC)
    ##splom(lpcobject$LPC)#, type=if ("curve" %in%  type) "l" else "p",labels= dimnames(lpcobject$data)[[2]] )

  }  else {
    cat("Data set has too many dimensions to be plotted in a pairs plot.\n\n")
  }


  invisible()

}



plot.lpc.spline <- plot.lpc


print.lpc.spline <- function(x, digits=max(3,getOption('digits')-3), ...){

  sx <- as.character(substitute(x))
  if (sum(nchar(sx))>200){ sx<-"name.of.this.object"}

  cat("\n")
  cat("Type plot(", sx, ") to see a graphical display of the fitted object. \n\n")
  cat("Type names(", sx, ") to see an overview of items available. \n\n")

  if(class(x)=="lpc"){
    if (x$scaled){
      cat("The data have been scaled by dividing through \n")
      cat(x$Misc$scaled.by)
      cat("\n")
    } else {
      cat("The data have not been scaled. \n")
    }
  }
  else if(class(x)=="lpc.spline"){
    cat("A cubic spline with ", dim(x$knots.coords[[1]])[2], " knots and total arc length ", diff(range(x$knots.pi[[1]])), " has been laid through the  local centers of mass representing the local principal curve. \n")
  }

}

#' Plot samples LPC line and scatter
#'
#'
#' @param res.lpc data after run get_biopath_score
#' @param normals a vector of the samples, if samples is Normal then TRUE, else Tumor then FALSE
#' @param outputHTML if need output this plot as a html file is TRUE, if not then FALSE
#' @param outputDir if outputHTML=T, this is availiable, the output file directory
#' @param fileName if outputHTML=T, this is availiable, the output file name.
#'
#' @return not have a return
#' @export
#' @importFrom rgl plot3d lines3d points3d writeWebGL
#'
#' @examples
#' library(biopathscore)
#' data("dat_brca")
#' data('v_brca_normals')
#' data('l_kegg_gs_min')
#' test_pathway = l_kegg_gs_min$KEGG_DNA_REPLICATION
#' pds = get_biopath_score(data = dat_brca, biopath_genes = test_pathway, normals = v_brca_normals, maximize_stability = F)
#' pds2 = get_biopath_score(data = dat_brca, biopath_genes = test_pathway, normals = v_brca_normals, maximize_stability = T)
#'
#' plot3D_lpc(pds, v_brca_normals,drawLine = F, outputHTML = F )
#' plot3D_lpc(pds2, v_brca_normals,drawLine = F, outputHTML = F )
#'
plot3D_lpc = function(res.lpc, normals, drawLine = F, outputHTML = F, outputDir = "./plot", fileName = "plot") {

  # res.lpc = pds
  p.pca = res.lpc$z
  p.curve = res.lpc$curves$closest.coords
  p.curve.order = res.lpc$curves_order
  p.curve.o = res.lpc$or.curve$LPC

  plot3d(p.pca[,c(1:3)], col = ifelse(normals, "green", "red"))
  lines3d(p.curve[p.curve.order, 1:3], color="blue")
  lines3d(p.curve.o[,c(1:3)], color="black")

  v.normal.mean = colMeans(p.pca[normals,])
  v.tumor.mean = colMeans(p.pca[!normals,])

  points3d(v.normal.mean[1], v.normal.mean[2],v.normal.mean[3], color="green", size=15.0)
  points3d(v.tumor.mean[1], v.tumor.mean[2],v.tumor.mean[3], color="red", size=15.0)

  if(drawLine == TRUE) {
    for(d in 1:nrow(p.pca)){
      lines3d(rbind(p.pca[d, c(1:3)], p.curve[d, 1:3]), color=ifelse(normals[d], "green", "red"))
    }
  }


  if(outputHTML == TRUE){
    if(!dir.exists(outputDir)) dir.create(outputDir, showWarnings = T, recursive = T)
    writeWebGL(filename=paste0(outputDir,"/",fileName, ".html"), snapshot=F)
  }

}
