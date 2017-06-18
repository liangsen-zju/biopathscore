lpc.splinefun <- function(lpcobject){
  nbranch <- max(lpcobject$Parametrization[,2])
  result <- list()
  for (r in 0:nbranch) {
    #print(range(lpcobject$Parametrization[lpcobject$Parametrization[,2]==r,1]))
    element <- list(branch=r, range=range(lpcobject$Parametrization[lpcobject$Parametrization[,2]==r,1]), splinefun=list())
    for (j in 1:ncol(lpcobject$LPC))
      element$splinefun[[j]] <- splinefun(lpcobject$Parametrization[lpcobject$Parametrization[,2]==r,1], lpcobject$LPC[lpcobject$Parametrization[,2]==r,j])
    result[[r+1]] <- element
  }
  result
}


# Evaluates the splines of the indicated branch at the original projection index or.pi
lpc.spline.eval <- function(lpcsl,or.pi,branch=0){

  if (length(branch)==1)
    branch <- rep(branch,length(branch))
  result <- matrix(nrow=length(or.pi), ncol=length(lpcsl[[1]]$splinefun))
  branches <- unique(branch)
  for (cur.branch in branches) {
    cur.or.pi <- or.pi[branch==cur.branch]
    for (j in 1:length(lpcsl[[cur.branch+1]]$splinefun))
      result[branch==cur.branch,j] <- lpcsl[[cur.branch+1]]$splinefun[[j]](cur.or.pi)  # here cur.op.pi instead of or.pi
  }
  result
}


lpc.spline <- function(lpcobject, optimize=TRUE, compute.Rc=FALSE,  project=FALSE, ...){
  if (class(lpcobject)=="lpc.spline"){
     lpcobject <-lpcobject$lpcobject
  }
  lpcsl    <- suppressWarnings(lpc.splinefun(lpcobject))
  spline   <- lpc.fit.spline(lpcsl,...)
  if (project || compute.Rc){
       proj           <- lpc.project.spline(lpcsl, newdata=lpcobject$data, optimize=optimize,...)
       closest.coords <- proj$closest.coords
  } else {
       closest.coords<- "none"
     }
  if (project){
       closest.pi     <- round(proj$closest.pi, digits=6)
       closest.dist   <- proj$closest.dist
       closest.branch <- proj$closest.branch
  } else {
      closest.pi <- closest.dist <-  closest.branch <-"none"
  }
  if (compute.Rc){
       R <- Rc(lpcobject$data, closest.coords)
  } else {
       R<-"none"
  }

  fit <- list(knots.pi=spline$knots.pi, knots.coords= spline$knots.coords, closest.pi=closest.pi, closest.coords=closest.coords, closest.dist=closest.dist,closest.branch=closest.branch, Rc=R, project=project, lpcobject=lpcobject, splinefun=lpcsl)
  class(fit)<-"lpc.spline"
  return(fit)
}


lpc.fit.spline <- function(lpcsl,  num.knots=100){
  nbranch <- length(lpcsl)-1
  # First of all evaluate the spline at each knot
  knots.pi <- list()
  knots.coords <- list()
  for (r in 0:nbranch) {
    knots.pi[[r+1]] <- round(lpcsl[[r+1]]$range[1] + 0:(num.knots-1)/(num.knots-1) * (lpcsl[[r+1]]$range[2]-lpcsl[[r+1]]$range[1]), digits=6)
    coords <- matrix(ncol=num.knots,nrow=length(lpcsl[[r+1]]$splinefun))
    for (j in 1:length(lpcsl[[r+1]]$splinefun))
      coords[j,] <- lpcsl[[r+1]]$splinefun[[j]](knots.pi[[r+1]])
    knots.coords[[r+1]] <- coords
  }
  list(knots.coords=knots.coords, knots.pi=knots.pi)
}



