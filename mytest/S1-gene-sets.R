# ********************************************
# * Name: sen@uga.edu
# * Date: Tue Jun 06 21:30:13 2017
# * Detail: null
# *         
# ********************************************
library(stringr)
root.path = "E:/Programming/CODE/LPC/"
data.path = "E:/Programming/DATA/TCGA-RData/"
data.path2 = paste0(root.path, "_data/")

################ MsigDB
### gene set from http://software.broadinstitute.org/gsea/downloads.jsp

gmt2list = function(file.path, ...)
{
  gmt = str_split(readLines(con = file(file.path)), "\t")
  v.names = sapply(gmt, function(line){line[1]}, simplify = T)
  gmt = lapply(gmt, function(line){line = line[-c(1:2)]; line = str_trim(line)})
  names(gmt) = v.names
  
  return(gmt)
}



l.hallmark.gs = gmt2list("E:/Programming/DATA/msigdb/h.all.v6.0.entrez.gmt")
save(l.hallmark.gs, file = paste0(data.path2, "l.hallmark.gs.RData"))

l.kegg.gs = gmt2list("E:/Programming/DATA/msigdb/c2.cp.kegg.v6.0.entrez.gmt")
save(l.kegg.gs, file = paste0(data.path2, "l.kegg.gs.RData"))

l.go.bp.gs = gmt2list("E:/Programming/DATA/msigdb/c5.bp.v6.0.entrez.gmt")
save(l.go.bp.gs, file = paste0(data.path2, "l.go.bp.gs.RData"))


l.kegg.min.gs = vector('list',5)
names(l.kegg.min.gs) = names(l.kegg.gs)[c(95,114,93,76,162)]
l.kegg.min.gs[[1]] = l.kegg.gs[[95]]
l.kegg.min.gs[[2]] = l.kegg.gs[[114]]
l.kegg.min.gs[[3]] = l.kegg.gs[[93]]
l.kegg.min.gs[[4]] = l.kegg.gs[[76]]
l.kegg.min.gs[[5]] = l.kegg.gs[[162]]
save(l.kegg.min.gs, file = paste0(data.path2, "l.kegg.min.gs.RData"))

