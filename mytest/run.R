library("R.oo")
#library("foreach")
#library("doParallel")
library(rgl)
library(scatterplot3d)

#library("doMC")
library(devtools)
library(roxygen2)
library(testthat)
library(knitr)
library(stringr)

root.path = "E:/Programming/CODE/biopathscore/biopathscore/mytest/"
data.path = "E:/Programming/DATA/TCGA-RData/"
data.path2 = paste0(root.path, "_data/")

origin_wd = getwd()
setwd(root.path)


load_all()

#load("./_data/tcga_data.RData",  verbose= T)
#load("./_data/l.hallmark.gs.RData",  verbose= T)
#load("./_data/l.kegg.gs.RData", verbose = T)
#load("./_data/l.go.bp.gs.RData", verbose = T)
#load("./_data/l.kegg.min.gs.RData", verbose = T)


#tcga = 'LIHC'
#print(tcga)

# registerDoMC(25)
# cl <- makeCluster(25)
# registerDoParallel(cl)
#
# geneset = l_kegg_gs_min$KEGG_P53_SIGNALING_PATHWAY
# attempts = 100
# maximize_stability=FALSE
# samplings=NULL
# data = dat_brca
# normals = v_normal
# ranks = v_normal
# logfile= paste0("./_result/tcga.brca.kegg.20170616.log")
#
#
#
#
#
# PDS <- quantify_pathways_deregulation( data = data, allgenes = allgenes,geneset = geneset,normals = normals, ranks=ranks,
#                                        attempts = attempts,logfile=logfile, min_exp=min_exp, min_std=min_std,
#                                        maximize_stability = maximize_stability,plot=plot,samplings=samplings)
#
#
# save(PDS, file=paste0("./_result/PDS.",tolower(tcga),".kegg.mini.RData"))
# setwd(origin_wd)
#
#


data = dat_brca
biopath_genes = l_kegg_gs_min$KEGG_DNA_REPLICATION
attempts = 100
maximize_stability=FALSE
samplings=NULL
normals = v_normal
ranks = v_normal

#attempts=100, maximize_stability=TRUE, logfile="", samplings=NULL,
pds = get_biopath_score(data, biopath_genes,normals)
pds2 = get_biopath_score(data, biopath_genes,normals,maximize_stability = T)

save(pds, file="./_result/KEGG_P53_SIGNALING_PATHWAY.pds.RData")

plot3D.lpc2(pds, normals, fileName = "kegg-dna-50")
plot3D.lpc2(pds2, normals, fileName = "kegg-dna-50-stb")


