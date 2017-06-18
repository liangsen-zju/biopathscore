".getmeasuredgenesinpathway"=
  function(syms,allgenes)
  {
    l=length(syms)
    pathways=vector('list',l)
    for (i in 1:l) {
      n=length(syms[[i]])
      isin=matrix(FALSE,n)
      for(j in 1:n) {
        isin[j]=(length(grep(paste('\\b',trim(syms[[i]][j]),'\\b',sep=""),allgenes))>0)
      }
      pathways[[i]]=unique(syms[[i]][isin])
    }
    pathways
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

".score_pathway"=
  function(pca,m,ranks,calcerr=FALSE,thresh=0.0005,maxit=200,start,logfile="")
  {
    #browser()
    # cutoff the gene whose sd lower than 0.001
    # pca=pca[,apply(pca,2,sd)>0.001] # not need

    #print(pca)
    k=dim(pca)[2]
    if (k<3) {
      c=NULL
      cat(file=logfile,append=TRUE,'scoring failed (k=',k,').\n')
    } else {
      d=matrix(0,1,m)

      if (start=="by pca") {
        start=NULL
      } else if (start=="by ranks") {
        # start=aggregate(pca, by=list(ranks), FUN=mean) #### Test 06/22/2016, stable version
        # start=as.matrix(start[,-1])
        start=apply(pca, 2, function(pca){weighted.mean(pca, ranks)}) # select centroid as start point avoiding mapping all tumor samples to one end of curve
      }


      # c=principal.curve(pca, ranks, start=start, thresh=thresh, maxit=maxit, plot.true=FALSE)
      c=lpc(pca, scaled=F, weights=ranks, x0=start, pen=2, control=lpc.control(iter=100, mult=1,cross=F, ms.sub=30))
      #invisible(c=lpc(pca, scaled=F, weights=ranks, x0=start, pen=2, control=lpc.control(iter=100, mult=1,cross=F, ms.sub=30)))
      # invisible(c=lpc(pca, scaled=F, weight=ranks, x0=start[2, ], control=lpc.control(iter=100, mult=1)))

    }
    ########################################################################################
    if (!is.null(c)) {
      p=lpc.project(c, pca)
      # tag=sort(p$closest.pi, index.return=TRUE, decreasing=FALSE)$ix

      tag=order(sapply(c(1:length(p$closest.pi)), function(i){ifelse(p$closest.branch[i]==1, p$closest.pi[i]+max(p$closest.pi[which(p$closest.branch==0)]), p$closest.pi[i])}), decreasing=FALSE)
      d[tag[1]]=0

      for (j in 2:m) {
        d[tag[j]]=d[tag[j-1]]+dist(p$closest.coords[tag[(j-1):j],])
      }
      d=d/d[tag[m]]

      if (calcerr) {
        e=matrix(0,1,k)
        for (i in 1:k) {
          e[i]=mean((p$closest.coords[,i]-pca[,i])^2)
        }
      } else {
        e=FALSE;
      }
      list(score=d, error=e, thecurve=p, thecurve.o=c, tag=tag)

    } else {
      cat(file=logfile,append=TRUE,'scoring failed.\n')
      NULL
    }
    ########################################################################################
    # if (!is.null(c)) {
    # 	d[c$tag[1]]=0
    # 	for (j in 2:m) {
    # 		d[c$tag[j]]=d[c$tag[j-1]]+dist(c$s[c$tag[(j-1):j],])
    # 	}
    # 	d=d/d[c$tag[m]]
    # 	if (calcerr) {
    # 		e=matrix(0,1,k)
    # 		for (i in 1:k) {
    # 			e[i]=mean((c$s[,i]-pca[,i])^2)
    # 		}
    # 	} else {
    # 		e=FALSE;
    # 	}
    # 	list(score=d,error=e,thecurve=c)
    # } else {
    # 	cat(file=logfile,append=TRUE,'scoring failed.\n')
    # 	NULL
    # }
  }

".samplings_stdev"=
  function(m,n,attempts,z,ranks,samplings,start,logfile="")
  {
    dall=array(dim=c(attempts,n))
    skip=0

    for(a in 1:attempts) {
      #print(a)
      res=.score_pathway(z[samplings[a,],],m,ranks[samplings[a,]],start=start,logfile=logfile)

      if (!is.null(res)) {
        dall[a,samplings[a,]]=res$score
      } else {
        skip=skip+1
      }
    }
    if (skip < attempts/2) {
      # mean of every column's sd
      mean(apply(dall,2,sd,'na.rm'=TRUE), 'na.rm'=TRUE)
    } else {
      Inf
    }

  }

".score_all_pathways_helper"=
  function(pca, ranks, samplings, attempts, maximize_stability, logfile="",start)
  {
    n=dim(pca)[1]	       ## samples' number
    k=dim(pca)[2]	       ## PCA dims
    m=dim(samplings)[2]   ## sampling the #samples

    mincheck=5
    kmin=max(floor(0.8*k),mincheck+1)
    mindelta=min(0.009,max(0.002,1.5/k))
    sig=matrix(0,1,k)

    ## compute the score of the pathway
    res=.score_pathway(pca,n,ranks,calcerr=TRUE,start=start,logfile=logfile)

    if (is.null(res)) {
      cat(file=logfile,append=TRUE,'pathway > scoring failed 1.\n')
    } else {

      ## using sampling: sig=average sd score by attempts times sampling
      sig=.samplings_stdev(m,n,attempts,pca,ranks,samplings,start=start)

      if (sig>10000) {
        ## sig=Inf represent skips(skip < attempts/2) too much
        cat(file=logfile,append=TRUE,'pathway > scoring failed 2 (sig:', sig, ').\n')
        res=NULL
      } else {

        ## vaild sampling
        origsig=sig
        cat(file=logfile,append=TRUE,'pathway > sig:', sig, '\n')

        isin=1:k
        if (maximize_stability) {
          testsig=max(mincheck,floor(0.1*k))
          newsig=rep(0,testsig)

          while ((k>=kmin)&(sig>0.05)) {
            se=sort(res$error,index.return=TRUE,decreasing=TRUE)

            for (j in 1:testsig) {
              newsig[j]=.samplings_stdev(m,n,attempts,pca[,-se$ix[j]],ranks,samplings,start=start)
            }

            ## why get the minime
            wj=which.min(newsig)
            # cat(file=logfile,append=TRUE,'pathway  k=', k, '(', ncol(res$thecurve$s), ') wj=', wj, '>new sig:', newsig[wj])

            cat(file=logfile,append=TRUE,'pathway k=', k, '(', ncol(res$thecurve$closest.coords), ') wj=', wj, '>new sig:', newsig[wj])

            if (sig-newsig[wj]<mindelta) {
              cat(file=logfile,append=TRUE,' x rejected\n')
              break
            }
            cat(file=logfile,append=TRUE,' | accepted!\n')

            ## update
            sig=newsig[wj]
            isin=isin[-se$ix[wj]];
            pca=pca[,-se$ix[wj]]
            k=k-1

            ## do compute score again.
            res=.score_pathway(pca,n,ranks,calcerr=TRUE,start=start,logfile=logfile)
            if (is.null(res)) {
              cat(file=logfile,append=TRUE,'pathway > scoring failed 3.\n')
              break;
            }
          }
        }
      }
    }

    ## return
    if (is.null(res)) {
      NULL
    } else {
      list(score=res$score,thecurve=res$thecurve,thecurve.o=res$thecurve.o,tag=res$tag,z=pca,isin=isin,sig=sig,origsig=origsig,k=k)
    }
  }

"quantify_pathways_deregulation"=
  function(data, allgenes, geneset, normals=NULL, ranks=NULL, attempts=100, maximize_stability=TRUE, logfile="", plot=T, samplings=NULL, min_exp=4, min_std=0.4)
  {
    cat(file=logfile,append=FALSE,'robust_score_bydist. min_exp=',min_exp,', min_std=',min_std,'\n')


    # data[data<min_exp]=min_exp		# the minime expression data

    n.samples=ncol(data)	   # the number of Samples

    ## check rank && normals
    ## ??? rank?? by pca, ranks??, if not have normal, than all are normal???
    if (is.null(normals)) { normals=rep(TRUE,n.samples);start="by pca"; } else {start="by ranks"}

    if (is.null(ranks)){ ranks=!normals;	}
    ranks=rank(ranks)


    if ((length(normals)!=n.samples)||(length(ranks)!=n.samples)) { stop("invalid dimentions")}

    n.geneset=length(geneset)	                                # the number of pathways
    n.normal=sum(normals)                                     # the number of Normal Samples
    n.sampling=floor(0.8*(n.samples-n.normal))+n.normal	      # sampling: 0.8 Tumor Samples + All Normal Samples

    pathwaynames=names(geneset)

    ## sample n.sampling samples by `attempts` times-> sampling matrix
    if (is.null(samplings)) {
      samplings=matrix(0,attempts,n.sampling)
      w=which(!normals)
      for(a in 1:attempts) {
        samplings[a,]=sort(c(w[sample(n.samples-n.normal,n.sampling-n.normal)],which(normals)))
      }
    }

    s=NULL
    ind=NULL

    s=foreach (i=1:n.geneset, .options.multicore=list(preschedule=FALSE), .combine=rbind, .inorder=FALSE, .verbose=FALSE, .errorhandling='stop', .multicombine=TRUE, .export=as.vector(lsf.str(all.names=T, envir=.GlobalEnv)), .packages='R.oo') %dopar% {
      # for (i in 1:n.geneset) {
      # if using pathifier for large number of pathways, you might want to use the doMC library to parallelize your code, in that case replace the above for with: s=foreach (i=1:n.geneset, .options.multicore=list(preschedule=FALSE), .combine=rbind, .inorder=FALSE, .verbose=FALSE, .errorhandling='stop', .multicombine=TRUE) %dopar% {
      # s=foreach (i=1:n.geneset, .options.multicore=list(preschedule=FALSE), .combine=rbind, .inorder=FALSE, .verbose=FALSE, .errorhandling='stop', .multicombine=TRUE, .export=as.vector(lsf.str(all.names=T, envir=.GlobalEnv)), .packages='R.oo') %dopar% {

      #i=1
      print(i)
      v.pathway.id=geneset[[i]]
      pathwayindata=.getpathway(v.pathway.id,allgenes,data)

      k1=sum(pathwayindata$isin)
      if (k1<3) {
        # control the vaild pathway gene
        si=NULL
        cat(file=logfile,append=TRUE,'skipping pathway ',i,' k1=', k1,'\n')
        return(si)
      }

      d.pathway=pathwayindata$x
      v.pathway.id=v.pathway.id[pathwayindata$isin]		# update pathway

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


      #	pca=scale(z, center=T, scale=T)
      #	pca=pca[, order(apply(pca,2,sd), decreasing=T)]
      #	k2=max(sum(apply(pca,2,sd)>0.5),4)
      #	k2=min(k2,k1,0.75*dim(d.pathway)[1],sum(apply(pca,2,sd)>0.25))

      k0=which(summary(t)$importance[3,]>0.85)[1]
      k2=max(sum(t$sdev>1.1),4)		  # t$sdev > 1.1
      k2=min(k1,k2,k0,0.75*dim(d.pathway)[1],sum(t$sdev>0.25))

      ############# 06/29/16 #
      # z=t(apply(z, 1, function(zi){(zi-min(zi))/(max(zi)-min(zi))}))
      # t=prcomp(z, center=T, scale=T) # oringal with pca
      # k2=max(sum(t$sdev>quantile(t$sdev, probs=0.6),4)) # oringal with pca
      # k2=min(k2,k1,0.75*dim(d.pathway)[1],sum(t$sdev>quantile(t$sdev, probs=0.4))) # oringal with pca

      # t=prcomp(z, center=T, scale=T) #1 modify with pca
      # k2=max(which(summary(t)$importance[3,]>0.85)[1], 4) #1 modify with pca
      # k2=min(k2,k1,0.75*dim(d.pathway)[1],sum(t$sdev>0.05))
      # print(summary(t)$importance[2:3,])

      if (k2<3) {
        si=NULL
        cat(file=logfile,append=TRUE,'skipping pathway ',i,' k2=', k2,'\n')
        return(si)
      }

      pca=t$x[,1:k2]                  #oringal
      # pca=sweep(pca, 2, sd(pca), "/")
      res=.score_all_pathways_helper(pca, ranks, samplings, i, attempts, maximize_stability, logfile, start=start)

      if (is.null(res)) {
        si=NULL
        cat(file=logfile,append=TRUE,'skipping pathway ',i,'\n')
        return(si)
      }

      ind=c(ind,i)
      # si=list(res$score,pathway,res$sig,res$origsig,res$k,res$thecurve$s,res$thecurve$tag,res$z,res$isin,xm,xs,t$center,t$rotation,k2) #oringal

      si=list(res$score,v.pathway.id,res$sig,res$origsig,res$k,
              res$thecurve$closest.coords,res$tag,res$z,res$isin,xm,
              xs,t$center,t$rotation,k2,sum(summary(t)$importance[2, res$isin]),
              res$thecurve.o,t$sdev, t$x) #1 modify

      # si=list(res$score,v.pathway.id,res$sig,res$origsig,res$k,res$thecurve$closest.coords,res$tag,res$z,res$isin,xm,xs,"t$center","t$rotation",k2,"sum(summary(t)$importance[2, res$isin])")
      # s=rbind(s,si) # if using doMC foreach above replace this line with simply: si

      return(si)
    }


    cat(file=logfile,append=TRUE,length(ind),'pathways processed with start=',start,'\n')
    rownames(s)=pathwaynames[ind]
    list(scores=s[,1], genesinpathway=s[,2], newmeanstd=s[,3], origmeanstd=s[,4], pathwaysize=s[,5],
         curves=s[,6], curves_order=s[,7], z=s[,8],compin=s[,9],xm=s[,10],
         xs=s[,11],center=s[,12],rot=s[,13],pctaken=s[,14],pv=s[,15],
         or.curve=s[,16],sdev=s[,17],pca = s[,18],samplings=samplings,sucess=ind,logfile=logfile)
  }



#' @export
get_biopath_score = function(data, biopath_genes,normals=NULL, ranks=NULL,attempts = 100, maximize_stability=F, use_min = T, min_exp = 2, min_std = 0.1)
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
  if (is.null(samplings)) {
    samplings=matrix(0,attempts,n.sampling)
    w=which(!normals)
    for(a in 1:attempts) {
      samplings[a,]=sort(c(w[sample(n.samples-n.normal,n.sampling-n.normal)],which(normals)))
    }
  }

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



