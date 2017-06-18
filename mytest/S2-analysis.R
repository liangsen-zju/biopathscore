# ********************************************
# * Name: sen@uga.edu
# * Date: Fri Jun 09 14:28:09 2017
# * Detail: null
# *         
# ********************************************

library(rgl)
library(scatterplot3d)

load('./_result/PDS.lihc.kegg.mini.RData',verbose = T)

v.pathway.name = c("KEGG_P53_SIGNALING_PATHWAY","KEGG_ECM_RECEPTOR_INTERACTION", "KEGG_CELL_CYCLE",
                   "KEGG_DNA_REPLICATION","KEGG_PATHWAYS_IN_CANCER")

pathway.name = "KEGG_ECM_RECEPTOR_INTERACTION";

p.pca = PDS$z[[which(names(PDS$z) == pathway.name)]]
p.curve = PDS$curves[[which(names(PDS$curves) == pathway.name)]]
p.curve.order = PDS$curves_order[[which(names(PDS$curves_order) == pathway.name)]]

p.curve.o = PDS$or.curve[[which(names(PDS$or.curve) == pathway.name)]]

plot3d(p.pca[,c(1:3)], col = ifelse(normals, "green", "red"))
lines3d(p.curve[p.curve.order, 1:3], color="blue")
lines3d(p.curve.o$LPC, color="black")


points3d(colMeans(p.pca[normals,c(1:3)]), color="black", size=26.0)

points3d(colMeans(pca[!normals,c(1:3)]), color="grey", size=16.0)


# scatterplotData=as.data.frame(pca)
# scatterplotMatrix(~.|factor(normals), data=pca, main=paste("pathway_", i, sep=""))

#plot3d(pca[,c(1:3)], col=clab)# sapply(normals, function(x)ifelse(x, "green", "red")))
plot3d(pca[,c(1:3)])# sapply(normals, function(x)ifelse(x, "green", "red")))
lines3d(res$thecurve$closest.coords[res$tag, 1:3], color="blue")
lines3d(res$thecurve.o$LPC, color="black")
#lines3d(res$thecurve.o$LPC, color="red")

for(d in 1:nrow(pca)){
  lines3d(rbind(pca[d, c(1:3)], res$thecurve$closest.coords[d, 1:3]), color=ifelse(normals[d], "green", "red"))
}
writeWebGL(filename=paste("webGL", as.character(paste("pathway_", i, ".html", sep="")), sep="/"), snapshot=F)