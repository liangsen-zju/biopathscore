# ********************************************
# * Name: sen@uga.edu
# * Date: Thu Jun 15 15:48:34 2017
# * Detail: Tide TCGA data
# *
# ********************************************

root.path = "E:/Programming/CODE/biopathscore/biopathscore/_src/"
data.path = "E:/Programming/DATA/TCGA-RData/"
data.path2 = paste0(root.path, "_data/")

origin_wd = getwd()
setwd(data.path)
v.tcga.project = unique(sapply(list.files("."), substr, 1, 4))

tcga_data <- lapply(v.tcga.project, function(c){

  #c= "BRCA"
  cc = tolower(c)

  ndata = get(load(paste0("./", c,"/",cc, ".fpkm.normal.RData", sep="")))
  tdata = get(load(paste0("./", c,"/",cc, ".fpkm.tumor.RData", sep="")))

  colnames(ndata) = paste0(colnames(ndata),"-N")
  colnames(tdata) = paste0(colnames(tdata),"-T")

  data<-cbind(ndata, tdata)
  normals<-c(rep(T, ncol(ndata)), rep(F, ncol(tdata)))

  data = as.matrix(round(data,2))
  data[which(data<4)]<-4                         # if data < 4 then data = 4
  data<-log2(data)                               # data = log(data)

  v.data.gene.sd = apply(data, 1, sd)
  std_cutoff<-quantile(v.data.gene.sd, probs=c(0.25, 0.95))
  min_std_cutoff = max(std_cutoff[1], 0.1)

  data<-data[which(v.data.gene.sd<=std_cutoff[2] & v.data.gene.sd>=min_std_cutoff), ]

  # data<-t(apply(data, 1, function(x){(x-min(x))/(max(x)-min(x))}))
  return(list("data"=data, "samples"=colnames(data), "normals"=normals, "min_exp"=2, "min_std"=min_std_cutoff, "allgenes"=rownames(data)))
})

names(tcga_data) <- v.tcga.project
save(tcga_data, file = paste0(data.path2, "tcga_data.RData"))


#
setwd(origin_wd)
load("./_src/_data/tcga_data.RData")
load("./_data/l.kegg.gs.RData", verbose = T)
load("./_data/l.kegg.min.gs.RData", verbose = T)

list_brca = tcga_data$BRCA

# ********************************************
# * Save for /data
# ********************************************
dat_brca = list_brca$data
dat_brca = round(dat_brca,2)
v_brca_normal = list_brca$normals
l_kegg_gs = l.kegg.gs
l_kegg_gs_min = l.kegg.min.gs


v_brca_normals = normals
names(v_brca_normals) = colnames(dat_brca)

devtools::use_data(v_brca_normals, overwrite = T)
devtools::use_data(dat_brca, v_brca_normals, overwrite = T)
devtools::use_data(l_kegg_gs, overwrite = T)
devtools::use_data(l_kegg_gs_min, overwrite = T)
