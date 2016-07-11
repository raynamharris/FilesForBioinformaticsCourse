### WGCNA ###
# the documentation provided by the developers if VERY useful and user friendly!
# https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/
# this will run pretty well on a regular laptop using 4000 - 5000 genes (probably 8000 max). if you really want to do more than that need to run it on a cluster. you can select the top DE genes by modifying the DESeq2 code we used yesterday. i also recommend using DESeq2 to variance transform the counts prior to running WGCNA (vsd transformation is an easy step in DESeq2) 

# first you will need to install it first. you will only need to do this once. 
source("http://bioconductor.org/biocLite.R")
biocLite("impute")
install.packages("WGCNA")

# load the package
library(WGCNA)
options(stringsAsFactors=FALSE)


# read in file with counts. CAUTION: counts must be NUMERIC type and vsd transformed.
counts <- read.csv(file = "vsd-counts.csv", row.names=1, stringsAsFactors = FALSE) 
head(counts)

# transpose file with counts so that each column is a gene and each row is a sample
data = as.data.frame(t(counts[]))


### CLEAN UP THE DATA ###

# check for samples with too many missing values 
gsg = goodSamplesGenes(data, verbose=3)
gsg$allOK # this should return TRUE. if this does NOT come back as true there is code in the documentation to return the names of the 'bad' samples and remove them from the data

# cluster the samples to see if there are outlies.
sampleTree = hclust(dist(data), method = "average")
par(cex=0.6)
par(mar=c(0,4,2,0))
plot(sampleTree, main = "sample clustering to detect outlies", sub="", xlab="", cex.lab=1.5, cex.axis=1.5, cex.main=2)

# if there are bad outliers they can be removed using the code below. 
abline(h=30, col="red") # plot a line adjusting the height to get rid of the offending sample(s)
clust = cutreeStatic(sampleTree, cutHeight=30, minSize=5)
table(clust)
keepSamples = (clust==1) # cluster 1 has samples we want to keep after above code
data = data[keepSamples,] # might want to change the name of new data set here if you want to keep the unfiltered one, too
nGenes = ncol(data)
nSamples = nrow(data)


### LOAD TRAIT INFO ###

traits <- read.csv(file = "traits.csv", stringsAsFactors = TRUE)
samples = rownames(data)
rowName <- match(samples, traits$fish)
names(traits) # if you need to you can use this info to remove any samples you don't need
datTraits = traits[rowName,c(3:11)] # this is a data frame of JUST the categorical variables
rownames(datTraits)=traits[rowName,1]


# we can now visualize the relationships between the counts and samples
sampleTree = hclust(dist(data), method = "average") # re-cluster in case you trimmed earlier
traitColors = numbers2colors(datTraits, signed=TRUE)
plotDendroAndColors(sampleTree, traitColors, groupLabels=names(datTraits), main="sample dendrogram and trait heatmap")


save(data, datTraits, file="WGCNA-BO-dataInput.RData") #this way you can come back to the loaded data more easily without have to reconstruct the dataset



### NETWORK CONSTRUCTION ###
# there are multiple ways to do this step. we are doing the step-by-step so we can see how things work, but the tutorials will show you an easier single step version as well. 

# choosing soft-thresholding power: analysis of network topology
powers = c(c(1:10), seq(from=12, to=40, by=2)) #choose a set of soft-thresholding powers
sft = pickSoftThreshold(data, powerVector = powers, verbose = 5, networkType="signed") #call the network topology analysis function


# plot the results
par(mfrow = c(1,2))
cex1 = 0.9
# scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="soft threshold (power)", ylab="scale free topology model fit, signed R^2", type="n", main=paste("scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, cex=cex1, col="red")
abline(h=0.90, col="red") #this line corresponds to using an R^2 cut-off of h
# mean connectivity as a funciton of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="soft threshold (power)", ylab="mean connectivity", type="n", main = paste("mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")


# based on the graphs above chose a soft thresholding power and fill it in here. 
softPower = 11
adjacency = adjacency(data, power = softPower, type="signed")

# turns adjacency into topological overlap. this minimizes effects of noise and spurious associations
TOM = TOMsimilarity(adjacency, TOMType="signed")
dissTOM = 1 - TOM

# cluster using TOM and plot clustering tree
geneTree = hclust(as.dist(dissTOM), method = "average")
plot(geneTree, xlab="", ylab="", sub="", main="gene clustering on TOM-based dissimilarity", labels=FALSE, hang=0.04)

# module identification using dynamic tree cut
minModuleSize = 30 # module size can be modified, developers prefer larger modules so they suggest 30
dynamicMods = cutreeDynamic(dendro = geneTree, distM=dissTOM, deepSplit=2, pamRespectsDendro=FALSE, minClusterSize=minModuleSize)

#output modules and numbers of genes in them. label "0" is for genes not assigned to any modules
table(dynamicMods)

#convert numeric labels into colors
dynamicColors = labels2colors(dynamicMods) 

# output a table of module colors and number of genes in them (grey always represents the unassigned module)
table(dynamicColors)

# plot dendrogram with colors for modules underneath
plotDendroAndColors(geneTree, dynamicColors, "dynamic tree cut", dendroLabels=FALSE, hang=0.03, addGuide=TRUE, guideHang=0.05, main="gene dendrogram and module colors")


# merging of modules whose expression profiles are very similar. 
MEList = moduleEigengenes(data, colors=dynamicColors) #calculate eigengenes
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs) # calculate dissimilarity of module eigengenes
METree = hclust(as.dist(MEDiss), method="average") #cluster module eigengenes
sizeGrWindow(7,6)
plot(METree, main="clustering of module eigengenes", xlab="", ylab="", sub="") #plot the results

#from the plot above you can choose a threshold to cut your eigengene cluster tree at
MEDissThres = 0.25 
abline(h=MEDissThres, col="red") #this plots the cutline chosen above

# call the automatic merging function
merge = mergeCloseModules(data, dynamicMods, cutHeight=MEDissThres, verbose=3) 

# relabel things so they line up post cutting, i.e. module merging
moduleLabels = merge$colors # numeric labels
mergedColors = labels2colors(moduleLabels) # the new merged colors
mergedMEs = merge$newMEs # the new merged module

# replot the dendrogram with the old module and the new merged modules
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("dynamic tree cut", "merged dynamic"), dendroLabels=FALSE, hang=0.03, addGuide=TRUE, guideHang = 0.05)

# rename and save things for subsequent analyses (using merged modules)
moduleColors = mergedColors # rename module colors to merged colors
colorOrder = c("grey", standardColors(50)) #construct numerical labels corresponding to the colors
moduleLabels = match(moduleColors, colorOrder) - 1
MEs = mergedMEs


save(MEs, moduleLabels, moduleColors, geneTree, file="WGCNA-BO-networkConstruction.RData")



### RELATE MODULES TO EXTERNAL TRAITS ###

# load data from both previous parts
lnames=load(file="WGCNA-BO-dataInput.RData") 
lnames 
lnames = load(file="WGCNA-BO-networkConstruction.RData")
lnames

# define numbers of genes and samples
nGenes = ncol(data)
nSamples = ncol(data)

# recalculate MEs with color labels
MEsO = moduleEigengenes(data, moduleColors)$eigengenes
MEs = orderMEs(MEsO)
moduleTraitCor = cor(MEs, datTraits, use="p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

# table with each association color coded by the correlation value.
sizeGrWindow(10,6)
# display correlations and their pvalues
textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep="")
dim(textMatrix) = dim(moduleTraitCor)
par(mar=c(6,8.5,3,3))
#display correlation values within a heatmap plot
labeledHeatmap(Matrix=moduleTraitCor, xLabels=names(datTraits), yLabels=names(MEs), ySymbols=names(MEs), colorLabels=FALSE, colors=blueWhiteRed(50), textMatrix=textMatrix, setStdMargins=FALSE, cex.text=0.5, zlim=c(-1,1), main=paste("module-trait relationships"))


# gene relationship to trait of important modules. example with weight as the trait of interest. 

trait = as.data.frame(datTraits$weight) #define variable trait containing the desired trait of interest (here weight). wrote the code this way so only have to substitue the trait name in this line of code and in the axis for the graph
names(trait) = "weight" #change name here to get the right annotation
modNames = substring(names(MEs), 3) #names (colors) of the modules

geneModuleMembership = as.data.frame(cor(data, MEs, use="p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

geneTraitSignificance = as.data.frame(cor(data, trait, use="p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) = paste("GS.", names(trait), sep="")
names(GSPvalue) = paste("p.GS", names(trait), sep="")


# intramodular analysis: identifying genes with high gene significance (GS) and module membership (MM)

module = "brown" # to get the names of the modules call 'modNames'
column = match(module, modNames)
moduleGenes = moduleColors==module

sizeGrWindow = (7,7)
par(mfrwo=c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]), abs(geneTraitSignificance[moduleGenes, 1]), xlab=paste("module membership in", module, "module"), ylab="gene significance for body weight", main=paste("module membership vs gene significance\n"), cex.main=1.2, cex.lab=1.2, cex.axis=1.2, col=module)


