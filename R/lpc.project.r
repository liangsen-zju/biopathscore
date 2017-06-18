
lpc.project <- function(object, newdata, ...){
  if (class(object)=="lpc"){
       lpcobject  <- object
       lpcsl  <- lpc.splinefun(lpcobject)
   } else {
       lpcobject  <- object$lpcobject
       lpcsl  <- object$splinefun
   }

   if (lpcobject$scaled){
      data <-sweep(as.matrix(newdata),2, lpcobject$Misc$scaled.by, "/")
   }   else {
      data <-newdata
    }


   result <- lpc.project.spline(lpcsl, data,...)[c(2,3,4,5,1)]
  result
}




# This function computes the projections of new data points onto the curve, VERY SLOW!
lpc.project.spline <- function(lpcsl, newdata, num.knots=100, optimize=TRUE) {

  squared.distance <- function(t,coords,sps) {
    dist <- 0
    for (i in 1:length(sps))
      dist <- dist + (coords[i]-sps[[i]](t))^2
    dist
  }

  nbranch <- length(lpcsl)-1

  fit.spline <- lpc.fit.spline(lpcsl, num.knots=num.knots)
  knots.coords  <- fit.spline[[1]]
  knots.pi      <- fit.spline[[2]]
  newdata <- as.matrix(newdata)

  closest.branch <- numeric(nrow(newdata))
  closest.or.pi <- numeric(nrow(newdata))
  closest.idx <- numeric(nrow(newdata))
  closest.dist <- numeric(nrow(newdata))

  for (i in 1:nrow(newdata)) {
    sqdist <- Inf
    for (r in 0:nbranch) {

      cur.sqdist <- apply((knots.coords[[r+1]]-newdata[i,])^2,2,sum)
      if (min(cur.sqdist)<sqdist) {
        closest.branch[i] <- r
        closest.idx[i] <- which.min(cur.sqdist)
        closest.or.pi[i] <- knots.pi[[r+1]][closest.idx[i]]
        sqdist <- min(cur.sqdist)
      }
    }
    if (optimize) {
      from <- knots.pi[[closest.branch[i]+1]][max(1,closest.idx[i]-1)]
      to <- knots.pi[[closest.branch[i]+1]][min(num.knots,closest.idx[i]+1)]
      if (from==to){
        closest.dist[i] <- sqrt(squared.distance(closest.or.pi[i],coords=newdata[i,],sps=lpcsl[[closest.branch[i]+1]]$splinefun))
      }  else {
        opt <- optimize(squared.distance,lower=from,upper=to,coords=newdata[i,],sps=lpcsl[[closest.branch[i]+1]]$splinefun)
        closest.or.pi[i] <- opt$minimum
        closest.dist[i] <- sqrt(opt$objective)
      }
    } else {
      closest.dist[i] <- sqrt(squared.distance(closest.or.pi[i],coords=newdata[i,],sps=lpcsl[[closest.branch[i]+1]]$splinefun))
    }
  }
  closest.pi <- lpc.curve.length(lpcsl,closest.or.pi,closest.branch)
  closest.coords <- lpc.spline.eval(lpcsl, closest.or.pi, closest.branch)

  result <- list(closest.branch=closest.branch,closest.pi=closest.pi,closest.or.pi=closest.or.pi,closest.coords=closest.coords,closest.dist=closest.dist)
  result
}


lpc.curve.length <- function(lpcsl,or.pi,branch=0,total.subdivisions=1e6,min.subdivisions=100){
  gradnorm <- function(tt,sps) {
    sum <- numeric(length(tt))
    for (j in 1:length(sps))
      sum <- sum + sps[[j]](tt,deriv=1)^2
    sqrt(sum)
  }
  if (length(branch)==1)
    branch <- rep(branch,length(or.pi))
  if (length(or.pi)!=length(branch))
    stop("Length mismatch between or.pi and branch")
  result <- numeric(length(or.pi))
  branches <- unique(branch)
  for (cur.branch in branches) {
    cur.or.pi <- or.pi[branch==cur.branch]
    length.between.points <- numeric(length(cur.or.pi))
    cur.or.pi.order <- order(cur.or.pi)
    cur.or.pi <- cur.or.pi[cur.or.pi.order]
    from <- c(0,cur.or.pi)
    for (i in 1:length(cur.or.pi.order)) {
      subdivisions <- max(min.subdivisions,ceiling(total.subdivisions*(cur.or.pi[i]-from[i])/(lpcsl[[cur.branch+1]]$range[2]-lpcsl[[cur.branch+1]]$range[1])))
      if(from[i]==cur.or.pi[i]){
        length.between.points[i] <-0
      } else {
        # base::save(gradnorm, i, from, cur.or.pi, lpcsl, cur.branch,  file="snap.RData",envir = environment())
        length.between.points[i] <- integrate(gradnorm,lower=from[i],upper=cur.or.pi[i],subdivisions=subdivisions, rel.tol=5e-03,sps=lpcsl[[cur.branch+1]]$splinefun)$value
      }
    }
    total.length <- cumsum(length.between.points)
    result[branch==cur.branch][cur.or.pi.order] <- total.length
  }
  result
}
