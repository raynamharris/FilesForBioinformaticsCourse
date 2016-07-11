### Plotting with EdgeR ###

# install the packages you will need. you will only need to do this once. 
source("http://bioconductor.org/biocLite.R")
biocLite("edgeR")
biocLite("cluster")
biocLite("Biobase")
biocLite("qvalue")

# now that the packages are installed we can load them
library(edgeR)
library(cluster)
library(Biobase)
library(qvalue)

# read your matrix table in & define your groups. be sure to have your working directory set to the location of your files:
samples <- read.delim(file = "samples.txt", stringsAsFactors = TRUE)
rnaseqMatrix = read.delim(file = "counts.txt", row.names=1, stringsAsFactors = FALSE)
rnaseqMatrix = rnaseqMatrix[rowSums(rnaseqMatrix)>=2,]
conditions = samples$pop


# next, run the edgeR statistical program, which will run an ANOVA between each group.
exp_study = DGEList(counts=rnaseqMatrix, group=conditions)
exp_study = calcNormFactors(exp_study)
exp_study = estimateCommonDisp(exp_study)
exp_study = estimateTagwiseDisp(exp_study)
et = exactTest(exp_study)
tTags = topTags(et,n=NULL)

# write your results to a file:
write.table(tTags, file='edgeR.DE_results', sep=' ', quote=F, row.names=T)

# now we are going to make a Volcano Plot to visualize our data:
source("rnaseq_plot_funcs.R")
pdf("edgeR.DE_results.MA_n_Volcano.pdf")
result_table = tTags$table
plot_MA_and_Volcano(result_table$logCPM, result_table$logFC, result_table$FDR)
dev.off()


# great. we have a volcano plot, but it would be great to generate a heat map. lets do that now:

# load some scripts and data we will need:
source("heatmap.3.R")
source("misc_rnaseq_funcs.R")
source("pairs3.R")

primary_data = read.csv("DEgenes_p0.05.csv", header=T, row.names=1)
primary_data = as.matrix(primary_data)

# now we are going to create a correlation matrix:
data = primary_data
sample_types = colnames(data)
nsamples = length(sample_types)
sample_colors = rainbow(nsamples)
sample_type_list = list()
for (i in 1:nsamples) {
    sample_type_list[[sample_types[i]]] = sample_types[i]
}
sample_factoring = colnames(data)
for (i in 1:nsamples) {
    sample_type = sample_types[i]
    replicates_want = sample_type_list[[sample_type]]
    sample_factoring[ colnames(data) %in% replicates_want ] = sample_type
}
data = log2(data+1)
sample_factoring = colnames(data)
for (i in 1:nsamples) {
    sample_type = sample_types[i]
    replicates_want = sample_type_list[[sample_type]]
    sample_factoring[ colnames(data) %in% replicates_want ] = sample_type
}
sampleAnnotations = matrix(ncol=ncol(data),nrow=nsamples)
for (i in 1:nsamples) {
  sampleAnnotations[i,] = colnames(data) %in% sample_type_list[[sample_types[i]]]
}
sampleAnnotations = apply(sampleAnnotations, 1:2, function(x) as.logical(x))
sampleAnnotations = sample_matrix_to_color_assignments(sampleAnnotations, col=sample_colors)
rownames(sampleAnnotations) = as.vector(sample_types)
colnames(sampleAnnotations) = colnames(data)
data = as.matrix(data) # convert to matrix
write.table(data, file="diffExpr.p0.05.log2.dat", quote=F, sep='	');
sample_cor = cor(data, method='pearson', use='pairwise.complete.obs')
sample_dist = dist(t(data), method='euclidean')
hc_samples = hclust(sample_dist, method='complete')
pdf("diffExpr.p0.05.log2.sample_cor_matrix.pdf")
heatmap.3(sample_cor, dendrogram='both', Rowv=as.dendrogram(hc_samples), Colv=as.dendrogram(hc_samples), col = greenred(75), scale='none', symm=TRUE, key=TRUE,density.info='none', trace='none', symkey=FALSE, margins=c(10,10), cexCol=1, cexRow=1, cex.main=0.75, main=paste("sample correlation matrix
", "diffExpr.P0.05_C2.matrix.log2") )
dev.off()


# and next we are going to create a heatmap:
gene_cor = NULL
gene_dist = dist(data, method='euclidean')
if (nrow(data) <= 1) { message('Too few genes to generate heatmap'); quit(status=0); }
hc_genes = hclust(gene_dist, method='complete')
myheatcol = greenred(75)
data = t(scale(t(data), scale=F)) # center rows, mean substracted
write.table(data, file="diffExpr.p0.05.log2.centered.dat", quote=F, sep='	');
heatmap_data = data
pdf("diffExpr.p0.05.log2.centered.genes_vs_samples_heatmap.pdf")
heatmap.3(heatmap_data, dendrogram='both', Rowv=as.dendrogram(hc_genes), Colv=as.dendrogram(hc_samples), col=myheatcol, scale="none", density.info="none", trace="none", key=TRUE, keysize=1.2, cexCol=1, margins=c(10,10), cex.main=0.75, main=paste("samples vs. features
", "diffExpr.P0.05_C2.matrix.log2.centered" ) )
dev.off()


# finally, lets save all our data:
save(list=ls(all=TRUE), file="diffExpr.P0.05_C2.matrix.RData")