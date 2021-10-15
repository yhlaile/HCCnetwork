setwd("/hccnetwork");options(stringsAsFactors = F);library(WGCNA);enableWGCNAThreads(12)
data=read.table("limor.fpkm.log2.txt",header=T,row.names = 1)[,-c(1:2)]
#filter FPKM <10
aa=apply(data,1,mean);data=subset(data,aa>3.323)
data=as.matrix(apply(data,2,rank, ties.method= "max"))
datSummary=rownames(data)
no.samples = dim(data)[[1]];
library(preprocessCore)
datExpr=t(data)
GeneName= datSummary
ArrayName= colnames(data)
powers=c(seq(1,10,by=1),seq(12,18,by=2));
sft=pickSoftThreshold(datExpr, powerVector=powers,networkType ="signed",corFnc =cor, corOptions =list(use = 'p'))
RpowerTable=sft[[2]]
sizeGrWindow(9, 5);
pdf('choosing power.pdf');
par(mfrow = c(1,2));cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],     xlab="Soft Threshold (power)",     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red");
dev.off()
# Mean connectivity as a function of the soft-thresholding power
sizeGrWindow(9, 5);
pdf('mean connectivity.pdf');
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",     ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red");
dev.off()
softPower =12
Connectivity=softConnectivity(datExpr,corFnc = "cor", corOptions ="use = 'p'",power=softPower,type="signed")
pdf("scale-free.pdf");
scaleFreePlot(Connectivity,nBreaks = 10,truncated = FALSE,removeFirst = FALSE, main = "");
dev.off()
adjacency = adjacency(datExpr,corFnc = "cor", corOptions ="use = 'p'",type = "signed", power = softPower)
TOM = TOMsimilarity(adjacency,TOMType="signed")
dissTOM = 1-TOM
geneTree = hclust(as.dist(dissTOM), method = "average")
minModuleSize =30;
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,deepSplit = 0,  pamRespectsDendro = FALSE,minClusterSize = minModuleSize,  cutHeight=0.99);
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

MEList = moduleEigengenes(datExpr, colors = dynamicMods)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average");
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",xlab = "", sub = "")
MEDissThres = 0.2
abline(h=MEDissThres, col = "red")
merge = mergeCloseModules(datExpr, dynamicMods, cutHeight = MEDissThres, verbose = 3);
mergedColors = merge$colors;
mergedMEs = merge$newMEs;
sizeGrWindow(12, 9)
pdf("DendroAndColors.pdf")
plotDendroAndColors(geneTree, cbind(dynamicMods, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"),dendroLabels = FALSE,   hang = 0.03,addGuide = TRUE, guideHang = 0.05)
dev.off()
moduleColors = mergedColors
colorOrder = c("grey", standardColors(unique(moduleColors)));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average")
pdf("METree.pdf")
plot(METree, main = "Clustering of module eigengenes",xlab = "", sub = "")
dev.off()
MEList = moduleEigengenes(datExpr, colors = dynamicMods)
nSamples=nrow(datExpr)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = cbind.data.frame(datSummary,corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
write.table(data.frame(ArrayName,MEs),"MEs.csv",row.name=F)
kMEdat=data.frame(geneModuleMembership,MMPvalue)
write.table(data.frame(datSummary,kMEdat),"kME-MMPvalue.csv",row.names=FALSE)
k.in=intramodularConnectivity(adjacency(datExpr,corFnc = "cor", corOptions = "use ='p'",type = "signed", power = softPower),moduleColors,scaleByMax = FALSE)
datout=data.frame(datSummary, colorNEW=moduleColors, k.in)
write.table(datout, file="OutputCancerNetwork.csv", sep=",", row.names=F)
hubs    = chooseTopHubInEachModule(datExpr, moduleColors)
write.csv(data.frame(module=names(hubs),moduleColor=c("grey",standardColors(length(hubs)-1)),hub=hubs), "num2color.csv",row.names=F)

gene=read.csv("OutputCancerNetwork.csv",header=T)
library(gProfileR)
for (i in unique(gene$colorNEW)[-4]){
  genes=subset(gene$datSummary,gene$colorNEW==i)
  go=gprofiler(as.vector(genes), organism = "hsapiens",numeric_ns="ENTREZGENE_ACC")[,-14]
  write.table(go,"module_enrichment.csv",append =T,row.names=rep(i,nrow(go)),sep=",")}

moduleColors=gene$colorNEW 
modules=unique(moduleColors)
n=length(modules)
pb <- txtProgressBar(min = 0, max = n, style = 3)
for (p in (1:n)[-4]){
  inModule = is.finite(match(moduleColors, modules[p]));
  dat2=data.frame(t(datExpr))[inModule,] 
  resamples=lapply(1:1000,function(i) a=t(sample(dat2[,1:nSamples],round(nSamples/2),replace=F)))
  K1=sapply(resamples,softConnectivity,power= softPower,type="signed")
  K=softConnectivity(t(dat2[,1:nSamples]),power= softPower,type="signed")
  #outfile=paste(modules[p],"-edit.txt",sep="")
  write.table(data.frame(mean(cor(K,K1)),apply(cor(K,K1),1,sd)), file = "module-stability-5.csv", row.names = modules[p], append = TRUE, col.names = FALSE, sep = ", ")
  setTxtProgressBar(pb, p)}
close(pb) 

#automatic finish the Cytoscape mods exporting
probes = dat0[,1]
n=length(unique(moduleColors))
pb <- txtProgressBar(min = 0, max = n, style = 3)
for (p in 1:n){ modules=unique(moduleColors)[p]
inModule = is.finite(match(moduleColors,modules));modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,threshold = quantile(abs(modTOM),probs=0.8),nodeNames = modProbes ,nodeAttr = moduleColors[inModule]);
setTxtProgressBar(pb, p)}
close(pb)

##############module preservation SELF
setLabels = c("Female", "Male");
datSummaryFemale=rownames(data)
datSummaryMale=rownames(data)
datExprFemale= datExpr
no.samplesFemale <- dim(datExprFemale)[[1]]
dim(datExprFemale)
datExprMale= t(read.table("limor.fpkm.log2.txt",header=T,row.names = 1)[,-c(1:2)]) 
colorsFemale = gene$colorNEW
colnames(datExprMale)=colnames(datExpr)
colnames(datExprFemale)=colnames(datExpr)
nSets = 2
ref = 1
test = 2
multiExpr = list(Female = list(data = datExprFemale), Male = list(data = datExprMale));
multiColor = list(Female = labels2colors(colorsFemale));
mp = modulePreservation(multiExpr, multiColor,referenceNetworks = 1,networkType="signed",nPermutations = 100,randomSeed = 1,parallelCalculation=T,quickCor = 0,verbose = 3)
save(mp, file = "modulePreservation.RData");#load("modulePreservation.RData")
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1]);
# Compare preservation to quality:
print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )
# Module labels and module sizes are also contained in the results
modColors = rownames(mp$preservation$observed[[ref]][[test]])
moduleSizes = mp$preservation$Z[[ref]][[test]][, 1];

plotMods = !(modColors %in% c("grey", "gold"));
# Text labels for points
labs = match(modColors[plotMods], standardColors(unique(modColors)-2))
# Auxiliary convenience variable
plotData = cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2])
# Main titles for the plot
mains = c("Preservation Median rank", "Preservation Zsummary");
# Start the plot
sizeGrWindow(10, 5);
pdf(file="FemaleOnly-modulePreservation-Zsummary-medianRank.pdf", wi=10, h=5,onefile=TRUE)
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))
for (p in 1:2)
{
  min = min(plotData[, p], na.rm = TRUE);
  max = max(plotData[, p], na.rm = TRUE);
  # Adjust ploting ranges appropriately
  if (p==2)
  {
    if (min > -max/10) min = -max/10
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  } else
    ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
  plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
       main = mains[p],
       cex = 1.3,
       ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
  labs = match(modColors[plotMods], standardColors(length(unique(modColors))))
  write.table(data.frame(mod))
  #replace text to labs as number labeling: labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1, offs = 0.08); 
  labelPoints(moduleSizes[plotMods], plotData[plotMods, p], labs, cex = 1, offs = 0.08)
  # For Zsummary, add threshold lines
  if (p==2)
  {
    abline(h=0)
    abline(h=2, col = "blue", lty = 2)
    abline(h=10, col = "darkgreen", lty = 2)
  }
}
# If plotting into a file, close it
dev.off();


# Re-initialize module color labels and sizes
modColors = rownames(statsZ)
moduleSizes = mp$quality$Z[[ref]][[test]][, 1];
# Exclude improper modules
plotMods = !(modColors %in% c("grey", "gold"));
# Create numeric labels for each module
labs = match(modColors[plotMods], standardColors(length(unique(modColors))));  #50 should larger than module number
# Start the plot: open a suitably sized graphical window and set sectioning and margins. Alternatively,
# plot into a pdf file.
sizeGrWindow(10, 9);
pdf(file="PreservationZStatistics.pdf", w=10, h=9)
par(mfrow = c(4,4))
par(mar = c(3,3,2,1))
par(mgp = c(1.6, 0.4, 0));
for (s in 1:ncol(statsZ))
{
  min = min(statsZ[plotMods, s], na.rm = TRUE);
  max = max(statsZ[plotMods, s], na.rm = TRUE);
  if (min > -max/12) min = -max/12
  plot(moduleSizes[plotMods], statsZ[plotMods, s], col = 1, bg = modColors[plotMods], pch = 21,
       main = colnames(statsZ)[s],
       cex = 1.1,
       ylab = colnames(statsZ)[s], xlab = "Module size", log = "x",
       ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min)),
       xlim = c(30, 1200),
       cex.lab = 1.2, cex.axis = 1.2)
  labelPoints(moduleSizes[plotMods], statsZ[plotMods, s], labs, cex = 1, offs = 0.06);
  #text(moduleSizes[-1], statsZ[-c(1:2), s], labels = letter[-c(1:2)], col = "black"); #modColors[-2]);
  abline(h=0)
  abline(h=2, col = "blue", lty = 2)
  abline(h=10, col = "darkgreen", lty = 2)
}
# If plotting into a file, close it, otherwise it is unreadable.
dev.off();


##############module preservation lccl
lccl=read.table("LCCL_RNAseqFPKM.txt",header=T,row.names = 1)
a=(duplicated(lccl[,1]))  
b=(lccl[a==FALSE,])
lccl <- data.frame(b[,-1], row.names = b[,1])
commongenes=intersect(rownames(data),rownames(lccl)) #rownames(data) as reference
lccl=lccl[which(rownames(lccl) %in% commongenes),]
lccl=as.matrix(apply(lccl,2,rank, ties.method= "max"))
data2=data[which(rownames(data) %in% commongenes),]
data2=as.matrix(apply(data2,2,rank, ties.method= "max"))
moduleColors
setLabels = c("Female", "Male");
datSummaryFemale=rownames(data2)
datSummaryMale=rownames(lccl)
datExprFemale= t(data2)
no.samplesFemale <- dim(datExprFemale)[[1]]
dim(datExprFemale)
datExprMale= t(lccl) 
colorsFemale = moduleColors[which(rownames(data) %in% commongenes)]
nSets = 2
ref = 1
test = 2
multiExpr = list(Female = list(data = datExprFemale), Male = list(data = datExprMale));
multiColor = list(Female = labels2colors(colorsFemale));
mp = modulePreservation(multiExpr, multiColor,referenceNetworks = 1,networkType="signed",nPermutations = 100,randomSeed = 1,parallelCalculation=T,quickCor = 0,verbose = 3)
save(mp, file = "modulePreservation.RData");
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1]);
# Compare preservation to quality:
print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )
# Module labels and module sizes are also contained in the results
modColors = rownames(mp$preservation$observed[[ref]][[test]])
moduleSizes = mp$preservation$Z[[ref]][[test]][, 1];

plotMods = !(modColors %in% c("grey", "gold"));
# Text labels for points
labs = match(modColors[plotMods], standardColors(unique(modColors)-2))
# Auxiliary convenience variable
plotData = cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2])
# Main titles for the plot
mains = c("Preservation Median rank", "Preservation Zsummary");
# Start the plot
sizeGrWindow(10, 5);
pdf(file="FemaleOnly-modulePreservation-Zsummary-medianRank.pdf", wi=10, h=5,onefile=TRUE)
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))
for (p in 1:2)
{
  min = min(plotData[, p], na.rm = TRUE);
  max = max(plotData[, p], na.rm = TRUE);
  # Adjust ploting ranges appropriately
  if (p==2)
  {
    if (min > -max/10) min = -max/10
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  } else
    ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
  plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
       main = mains[p],
       cex = 1.4,
       ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
  labs = match(modColors[plotMods], standardColors(length(unique(modColors))))
  write.table(data.frame(mod))
  #replace text to labs as number labeling: labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1, offs = 0.08); 
  labelPoints(moduleSizes[plotMods], plotData[plotMods, p], labs, cex = 1, offs = 0.08)
  # For Zsummary, add threshold lines
  if (p==2)
  {
    abline(h=0)
    abline(h=2, col = "blue", lty = 2)
    abline(h=10, col = "darkgreen", lty = 2)
  }
}
# If plotting into a file, close it
dev.off();


# Re-initialize module color labels and sizes
modColors = rownames(statsZ)
moduleSizes = mp$quality$Z[[ref]][[test]][, 1];
# Exclude improper modules
plotMods = !(modColors %in% c("grey", "gold"));
# Create numeric labels for each module
labs = match(modColors[plotMods], standardColors(length(unique(modColors))));  #50 should larger than module number
# Start the plot: open a suitably sized graphical window and set sectioning and margins. Alternatively,
# plot into a pdf file.
sizeGrWindow(10, 9);
pdf(file="PreservationZStatistics.pdf", w=10, h=9)
par(mfrow = c(4,4))
par(mar = c(3,3,2,1))
par(mgp = c(1.6, 0.4, 0));
for (s in 1:ncol(statsZ))
{
  min = min(statsZ[plotMods, s], na.rm = TRUE);
  max = max(statsZ[plotMods, s], na.rm = TRUE);
  if (min > -max/12) min = -max/12
  plot(moduleSizes[plotMods], statsZ[plotMods, s], col = 1, bg = modColors[plotMods], pch = 21,
       main = colnames(statsZ)[s],
       cex = 1.2,
       ylab = colnames(statsZ)[s], xlab = "Module size", log = "x",
       ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min)),
       xlim = c(30, 1200),
       cex.lab = 1.2, cex.axis = 1.2)
  labelPoints(moduleSizes[plotMods], statsZ[plotMods, s], labs, cex = 1, offs = 0.06);
  #text(moduleSizes[-1], statsZ[-c(1:2), s], labels = letter[-c(1:2)], col = "black"); #modColors[-2]);
  abline(h=0)
  abline(h=2, col = "blue", lty = 2)
  abline(h=10, col = "darkgreen", lty = 2)
}
# If plotting into a file, close it, otherwise it is unreadable.
dev.off();

##############module preservation ccle
setLabels = c("Female", "Male");
datSummaryFemale=rownames(data)
no.samplesFemale <- dim(datExprFemale)[[1]]
dim(datExprFemale)
ccle=read.table("ccle.txt",header=T)[-1,]
a=(duplicated(ccle[,1]))  
b=(ccle[a==FALSE,])
ccle<- data.frame(b[,-1], row.names = b[,1])
commongenes=intersect(rownames(data),rownames(ccle)) #rownames(data) as reference
ccle=ccle[which(rownames(ccle) %in% commongenes),]
ccle=as.matrix(apply(ccle,2,rank, ties.method= "max"))
data2=data[which(rownames(data) %in% commongenes),]
data2=as.matrix(apply(data2,2,rank, ties.method= "max"))
moduleColors
setLabels = c("Female", "Male");
datSummaryFemale=rownames(data2)
datSummaryMale=rownames(ccle)
datExprFemale= t(data2)
no.samplesFemale <- dim(datExprFemale)[[1]]
dim(datExprFemale)
datExprMale= t(ccle) 
colorsFemale = moduleColors[which(rownames(data) %in% commongenes)]
nSets = 2
ref = 1
test = 2
multiExpr = list(Female = list(data = datExprFemale), Male = list(data = datExprMale));
multiColor = list(Female = labels2colors(colorsFemale));
mp = modulePreservation(multiExpr, multiColor,referenceNetworks = 1,networkType="signed",nPermutations = 100,randomSeed = 1,parallelCalculation=T,quickCor = 0,verbose = 3)
save(mp, file = "modulePreservation.RData");
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1]);
# Compare preservation to quality:
print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )
# Module labels and module sizes are also contained in the results
modColors = rownames(mp$preservation$observed[[ref]][[test]])
moduleSizes = mp$preservation$Z[[ref]][[test]][, 1];

plotMods = !(modColors %in% c("grey", "gold"));
# Text labels for points
labs = match(modColors[plotMods], standardColors(unique(modColors)-2))
# Auxiliary convenience variable
plotData = cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2])
# Main titles for the plot
mains = c("Preservation Median rank", "Preservation Zsummary");
# Start the plot
sizeGrWindow(10, 5);
pdf(file="FemaleOnly-modulePreservation-Zsummary-medianRank.pdf", wi=10, h=5,onefile=TRUE)
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))
for (p in 1:2)
{
  min = min(plotData[, p], na.rm = TRUE);
  max = max(plotData[, p], na.rm = TRUE);
  # Adjust ploting ranges appropriately
  if (p==2)
  {
    if (min > -max/10) min = -max/10
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  } else
    ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
  plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
       main = mains[p],
       cex = 1.4,
       ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
  labs = match(modColors[plotMods], standardColors(length(unique(modColors))))
  write.table(data.frame(mod))
  #replace text to labs as number labeling: labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1, offs = 0.08); 
  labelPoints(moduleSizes[plotMods], plotData[plotMods, p], labs, cex = 1, offs = 0.08)
  # For Zsummary, add threshold lines
  if (p==2)
  {
    abline(h=0)
    abline(h=2, col = "blue", lty = 2)
    abline(h=10, col = "darkgreen", lty = 2)
  }
}
# If plotting into a file, close it
dev.off();


# Re-initialize module color labels and sizes
modColors = rownames(statsZ)
moduleSizes = mp$quality$Z[[ref]][[test]][, 1];
# Exclude improper modules
plotMods = !(modColors %in% c("grey", "gold"));
# Create numeric labels for each module
labs = match(modColors[plotMods], standardColors(length(unique(modColors))));  #50 should larger than module number
# Start the plot: open a suitably sized graphical window and set sectioning and margins. Alternatively,
# plot into a pdf file.
sizeGrWindow(10, 9);
pdf(file="PreservationZStatistics.pdf", w=10, h=9)
par(mfrow = c(4,4))
par(mar = c(3,3,2,1))
par(mgp = c(1.6, 0.4, 0));
for (s in 1:ncol(statsZ))
{
  min = min(statsZ[plotMods, s], na.rm = TRUE);
  max = max(statsZ[plotMods, s], na.rm = TRUE);
  if (min > -max/12) min = -max/12
  plot(moduleSizes[plotMods], statsZ[plotMods, s], col = 1, bg = modColors[plotMods], pch = 21,
       main = colnames(statsZ)[s],
       cex = 1.2,
       ylab = colnames(statsZ)[s], xlab = "Module size", log = "x",
       ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min)),
       xlim = c(30, 1200),
       cex.lab = 1.2, cex.axis = 1.2)
  labelPoints(moduleSizes[plotMods], statsZ[plotMods, s], labs, cex = 1, offs = 0.06);
  #text(moduleSizes[-1], statsZ[-c(1:2), s], labels = letter[-c(1:2)], col = "black"); #modColors[-2]);
  abline(h=0)
  abline(h=2, col = "blue", lty = 2)
  abline(h=10, col = "darkgreen", lty = 2)
}
# If plotting into a file, close it, otherwise it is unreadable.
dev.off();

#####other datase projection,note match genes and moduleColors
gse10140=read.table("GSE10140.txt",sep=",",header=T)
a=(duplicated(gse10140[,2]))  
b=(gse10140[a==FALSE,])
gse10140 <- data.frame(b[,-c(1:3)], row.names = b[,2])
commongenes=intersect(rownames(data),rownames(gse10140)) #rownames(data) as reference
gse10140=as.matrix(apply(gse10140[commongenes,],2,rank, ties.method= "max"))
ME_gse10140 = moduleEigengenes(t(gse10140[commongenes,]), colors = moduleColors[which(rownames(data) %in% commongenes)])
ME_gse10140=ME_gse10140$eigengenes
ArrayName_gse10140=colnames(gse10140)
write.table(data.frame(ArrayName_gse10140,ME_gse10140),"ME_gse10140.csv",row.name=F)

gse10141=read.table("GSE10141.txt",sep=",",header=T)
a=(duplicated(gse10141[,2]))  
b=(gse10141[a==FALSE,])
gse10141 <- data.frame(b[,-c(1:3)], row.names = b[,2])
commongenes=intersect(rownames(data),rownames(gse10141)) #rownames(data) as reference
gse10141=as.matrix(apply(gse10141[commongenes,],2,rank, ties.method= "max"))
ME_gse10141 = moduleEigengenes(t(gse10141[commongenes,]), colors = moduleColors[which(rownames(data) %in% commongenes)])
ME_gse10141=ME_gse10141$eigengenes
ArrayName_gse10141=colnames(gse10141)
write.table(data.frame(ArrayName_gse10141,ME_gse10141),"ME_gse10141.csv",row.name=F)

gse107170=read.table("GSE107170.txt",sep=",",header=T)
a=(duplicated(gse107170[,2]))  
b=(gse107170[a==FALSE,])
c=which(is.na(b[,2]))
gse107170 <- data.frame(b[-c,-c(1:3)], row.names = b[,2][-c])
commongenes=intersect(rownames(data),rownames(gse107170)) #rownames(data) as reference
gse107170=as.matrix(apply(gse107170[commongenes,],2,rank, ties.method= "max"))
ME_gse107170 = moduleEigengenes(t(gse107170[commongenes,]), colors = moduleColors[which(rownames(data) %in% commongenes)])
ME_gse107170=ME_gse107170$eigengenes
ArrayName_gse107170=colnames(gse107170)
write.table(data.frame(ArrayName_gse107170,ME_gse107170),"ME_gse107170.csv",row.name=F)

gse102079=read.table("GSE102079.txt",sep=",",header=T)
a=(duplicated(gse102079[,2]))  
b=(gse102079[a==FALSE,])
c=which(is.na(b[,2]))
gse102079 <- data.frame(b[-c,-c(1:3)], row.names = b[,2][-c])
commongenes=intersect(rownames(data),rownames(gse102079)) #rownames(data) as reference
gse102079=as.matrix(apply(gse102079[commongenes,],2,rank, ties.method= "max"))
ME_gse102079 = moduleEigengenes(t(gse102079[commongenes,]), colors = moduleColors[which(rownames(data) %in% commongenes)])
ME_gse102079=ME_gse102079$eigengenes
ArrayName_gse102079=colnames(gse102079)
write.table(data.frame(ArrayName_gse102079,ME_gse102079),"ME_gse102079.csv",row.name=F)

gse14323=read.table("GSE14323.txt",sep=",",header=T)
a=(duplicated(gse14323[,2]))  
b=(gse14323[a==FALSE,])
c=which(is.na(b[,2]))
gse14323 <- data.frame(b[-c,-c(1:3)], row.names = b[,2][-c])
commongenes=intersect(rownames(data),rownames(gse14323)) #rownames(data) as reference
gse14323=as.matrix(apply(gse14323[commongenes,],2,rank, ties.method= "max"))
ME_gse14323 = moduleEigengenes(t(gse14323[commongenes,]), colors = moduleColors[which(rownames(data) %in% commongenes)])
ME_gse14323=ME_gse14323$eigengenes
ArrayName_gse14323=colnames(gse14323)
write.table(data.frame(ArrayName_gse14323,ME_gse14323),"ME_gse14323.csv",row.name=F)

gse69844=read.table("GSE69844.txt",sep="\t",header=T)
tes=gProfileR::gconvert(gse69844[,1],organism = "hsapiens", target = "HGNC",mthreshold = 1,filter_na = F)
gse69844[,1]=tes[,4]
a=(duplicated(gse69844[,1]))  
b=(gse69844[a==FALSE,])
c=which(is.na(b[,1]))
gse69844 <- data.frame(b[,-1], row.names = b[,1])
commongenes=intersect(rownames(data),rownames(gse69844)) #rownames(data) as reference
gse69844=as.matrix(apply(gse69844[commongenes,],2,rank, ties.method= "max"))
ME_gse69844 = moduleEigengenes(t(gse69844[commongenes,]), colors = moduleColors[which(rownames(data) %in% commongenes)])
ME_gse69844=ME_gse69844$eigengenes
ArrayName_gse69844=colnames(gse69844)
write.table(data.frame(ArrayName_gse69844,ME_gse69844),"ME_gse69844.csv",row.name=F)

lccl=read.table("LCCL_RNAseqFPKM.txt",sep="\t",header=T)
a=(duplicated(lccl[,2]))  
b=(lccl[a==FALSE,])
c=which(is.na(b[,1]))
lccl <- data.frame(b[,-(1:2)], row.names = b[,2])
commongenes=intersect(rownames(data),rownames(lccl)) #rownames(data) as reference
lccl=as.matrix(apply(lccl[commongenes,],2,rank, ties.method= "max"))
ME_lccl = moduleEigengenes(t(lccl[commongenes,]), colors = moduleColors[which(rownames(data) %in% commongenes)])
ME_lccl=ME_lccl$eigengenes
ArrayName_lccl=colnames(lccl)
write.table(data.frame(ArrayName_lccl,ME_lccl),"ME_lccl.csv",row.name=F,sep=",")

ccle=read.table("ccle.txt",header=T)[-1,]
a=(duplicated(ccle[,1]))  
b=(ccle[a==FALSE,])
ccle<- data.frame(b[,-1], row.names = b[,1])
commongenes=intersect(rownames(data),rownames(ccle)) #rownames(data) as reference
ccle=as.matrix(apply(ccle[commongenes,],2,rank, ties.method= "max"))
ME_ccle = moduleEigengenes(t(ccle[commongenes,]), colors = moduleColors[which(rownames(data) %in% commongenes)])
ME_ccle=ME_ccle$eigengenes
ArrayName_ccle=colnames(ccle)
write.table(data.frame(ArrayName_ccle,ME_ccle),"ME_ccle.csv",row.name=F,sep=",")
