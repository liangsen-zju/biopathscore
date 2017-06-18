# ********************************************
# * Name: sen@uga.edu
# * Date: Sun Jun 18 12:14:34 2017
# * Detail: null
# *
# ********************************************

#' @export
get_biopath_score = function(data, biopath_genes,normals=NULL, ranks=NULL,attempts = 50, maximize_stability=F, use_min = T, min_exp = 2, min_std = 0.1)
{

  logfile = format(Sys.time(), "%Y%m%d-%H%M%S.log")
  cat(file=logfile,append=FALSE,'robust_score_bydist. min_exp=',min_exp,', min_std=',min_std,'\n')

  v_genes = rownames(data)
  v_normals = normals
  v_ranks = ranks
  v_biopath_genes = biopath_genes

  dat_min_exp = round(min(data),2)
  dat_min_std = round(min(apply(data, 1, sd)),2)

  # ***********
  # Min_exp &&& min_std
  # cut-off by min_std
  min_exp = max(dat_min_exp, min_exp)
  data[data<min_exp]=min_exp		# the minime expression data

  min_std = max(dat_min_std, min_std)
  v.data.sd = apply(data, 1, sd)
  data = data[which(v.data.sd >= min_std), ]

  n.samples=ncol(data)	   # the number of Samples

  ## check rank && normals
  ## ??? rank?? by pca, ranks??, if not have normal, than all are normal???
  if (is.null(normals)) { normals=rep(TRUE,n.samples);start="by pca"; } else {start="by ranks"}

  if (is.null(ranks)){ ranks=!normals;	}
  ranks=rank(ranks)


  if ((length(normals)!=n.samples)||(length(ranks)!=n.samples)) { stop("invalid dimentions")}

  n.normal=sum(normals)                                     # the number of Normal Samples
  n.sampling=floor(0.8*(n.samples-n.normal))+n.normal	      # sampling: 0.8 Tumor Samples + All Normal Samples

  ## sample n.sampling samples by `attempts` times-> sampling matrix
  ##if (is.null(samplings)) {
  samplings=matrix(0,attempts,n.sampling)
  w=which(!normals)
  for(a in 1:attempts) {
    samplings[a,]=sort(c(w[sample(n.samples-n.normal,n.sampling-n.normal)],which(normals)))
  }
  ##}

  pathwayindata = get_biopath_data(data, v_biopath_genes)


  k1=sum(pathwayindata$isin)
  if (k1<3) {
    # control the vaild pathway gene
    si=NULL
    cat(file=logfile,append=TRUE,'skipping pathway ',i,' k1=', k1,'\n')
    return(si)
  }

  d.pathway=pathwayindata$x
  v.pathway.id=v_biopath_genes[pathwayindata$isin]		# update pathway

  ## normal samples' gene's mean and sd
  ## normal's mean && sd is not same, Tumor's mean && sd are average mean && sd
  xm=apply(d.pathway[normals,],2,mean)
  xs=apply(d.pathway[normals,],2,sd)

  ## +LS:2017-06-10 add a cut-off by local pathway data
  v.gene.del = which(xs < min_std)
  v.pathway.id = v.pathway.id[-v.gene.del]
  d.pathway = d.pathway[,-v.gene.del]
  k1=length(v.pathway.id)

  xm=apply(d.pathway[normals,],2,mean)
  xs=apply(d.pathway[normals,],2,sd)

  xm_i=t(sapply(1:length(which(normals)), function(j){colMeans(d.pathway[which(normals)[-j],])}))
  xm=matrix(rep(xm,each=n.samples), nrow=n.samples)
  xm[normals,]=xm_i


  xs_i=t(sapply(1:length(which(normals)), function(j){apply(d.pathway[which(normals)[-j],],2,sd)}))
  xs=matrix(rep(xs,each=n.samples), nrow=n.samples)
  xs[normals,]=xs_i

  xs[xs<min_std]=min_std


  if (0 %in% xs) {
    si=NULL
    cat(file=logfile,append=TRUE,'skipping pathway ',i,' (0 in xs)\n')
    return(si)
  }

  ## normalization expression data with (d.pathway-u)/sd
  z = (d.pathway - xm ) / xs
  t=prcomp(z, center=T, scale=T)

  k0=which(summary(t)$importance[3,]>0.85)[1]
  k2=max(sum(t$sdev>1.1),4)		  # t$sdev > 1.1
  k2=min(k1,k2,k0,0.75*dim(d.pathway)[1],sum(t$sdev>0.25))

  if (k2<3) {
    si=NULL
    cat(file=logfile,append=TRUE,'skipping pathway ',i,' k2=', k2,'\n')
    return(si)
  }

  pca=t$x[,1:k2]                  #oringal
  # pca=sweep(pca, 2, sd(pca), "/")
  res=.score_all_pathways_helper(pca, ranks, samplings, attempts, maximize_stability, logfile, start=start)

  if (is.null(res)) {
    si=NULL
    cat(file=logfile,append=TRUE,'skipping pathway ',i,'\n')
    return(si)
  }

  si=list(scores = res$score,genesinpathway = v.pathway.id,newmeanstd = res$sig,origmeanstd = res$origsig,pathwaysize = res$k,
          curves = res$thecurve,curves_order = res$tag,or.curve = res$thecurve.o ,z = res$z,compin = res$isin,
          xm = xm, xs = xs,k2 = k2, pv = sum(summary(t)$importance[2, res$isin]),
          pca = t,samplings=samplings,logfile=logfile) #1 modify

  return(si)
}


#' @importFrom R.oo trim
get_biopath_data = function(data, v_biopath_genes)
{

  l=length(v_biopath_genes)
  x=NULL
  isin=rep(FALSE,l);
  allgenes = rownames(data)
  for(i in 1:l) {
    ind=unique(grep(paste('\\b',trim(v_biopath_genes[i]),'\\b',sep=""),allgenes))
    n=length(ind)
    if (n>0) {
      if (n==1)
        t=data[ind,]
      else
        t=colMeans(data[ind,])
      if (var(t)>0) {
        isin[i]=TRUE;
        x=c(x,t)
      }
    }
  }
  if (is.null(x)) {
    list(x=NULL,isin=isin)
  } else {
    list(x=matrix(x,nrow=ncol(data)),isin=isin)
  }
}
