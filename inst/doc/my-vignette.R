## ---- warning=FALSE------------------------------------------------------
citation("biopathscore")

## ---- eval=FALSE---------------------------------------------------------
#  install.packages("devtools")
#  devtools::install_github("liangsen35/biopathscore")
#  library("biopathscore")
#  

## ---- eval=FALSE---------------------------------------------------------
#  # install.packages("biopathscore")

## ---- results='hide'-----------------------------------------------------
library(biopathscore)
data("dat_brca")
data('v_brca_normals')
data('l_kegg_gs_min')

dat_brca[1:10, 1:10]

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(dat_brca[1:10, 1:10])

## ------------------------------------------------------------------------
dim(dat_brca)
length(v_brca_normals)
head(v_brca_normals)
length(l_kegg_gs_min)
head(l_kegg_gs_min$KEGG_P53_SIGNALING_PATHWAY)


## ---- eval = F-----------------------------------------------------------
#  library(biopathscore)
#  data("dat_brca")
#  data('v_brca_normals')
#  data('l_kegg_gs_min')
#  
#  
#  test_pathway = l_kegg_gs_min$KEGG_DNA_REPLICATION
#  
#  
#  ### Examples 1: not use maximize_stability = F, (if run slow, set attempts lower number)
#  pds = get_biopath_score(data = dat_brca, biopath_genes = test_pathway, normals = v_brca_normals,
#                          attempts = 100, maximize_stability = F)
#  
#  
#  plot3D_lpc(pds, v_brca_normals,drawLine = F, outputHTML = F )
#  plot3D_lpc(pds, v_brca_normals,drawLine = T, outputHTML = F )
#  
#  
#  ### Examples 1: use maximize_stability = T, very slow
#  pds2 = get_biopath_score(data = dat_brca, biopath_genes = test_pathway, normals = v_brca_normals,
#                          attempts = 100, maximize_stability = T)
#  
#  
#  plot3D_lpc(pds2, v_brca_normals,drawLine = F, outputHTML = F )
#  plot3D_lpc(pds2, v_brca_normals,drawLine = T, outputHTML = F )
#  
#  

