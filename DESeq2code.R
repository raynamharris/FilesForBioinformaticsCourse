### ANALYSIS OF DIFFERENTIAL GENE EXPRESSION WITH DESEQ2 ###


# before running this you will need to install the DESeq2 package if you have not already. you will only need to do this once. 
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")

# you will need to access the module each time you start a new R session
library(DESeq2)


# read in file with sample information
samples <- read.delim(file = "samples.txt", stringsAsFactors = TRUE)

#set reference levels (ie one that other will be compared to), default is alphabetical
samples$pop <- relevel(samples$pop, ref="HP") 
head(samples)

# read in file with counts. default is tab delimiter, have to specify if you have a different delimiter. make sure the column with the contig names DOES NOT have a header and that all values are integers.
counts <- read.delim(file = "counts.txt", row.names=1, stringsAsFactors = FALSE)
head(counts)



# construct DESeq database. we are using a simple design here, but much more complex designs are also possible (see documentation)
d <- DESeqDataSetFromMatrix(countData = counts, colData = samples, design = ~pop )
d


# run DE analysis!
d <- DESeq(d)


# we can have a look at our results and also get some basic summary information
res <- results(d)
res
summary(res)

# now we can zoom in on specific results and save some information.
#to get names for the effects. the order and names help you get the specific results you're interested in below.
resultsNames(d)

# look at desired effect. this gets more complex with a more complex design. here we are just going to pull out differences between populations and look at the change in the low-predation (LP) as compared to the high-predation (HP population)
# put the results into a new variable
resPop <- results(d, name="popLP") 
# sort the results by adjusted p-value (ie to see the most differentially expressed genes)
resPopOrdered <- resPop[order(resPop$padj),] 
head(resPopOrdered)
# note you can change the p-value here to your desired cutoff
table(resPop$padj<0.1)

# make single file with LFC and pvals by pulling those variables and then binding them together (if you have other effects of interested you could employ a similar approach)
pop.LFC = resPop$log2FoldChange
pop.pval = resPop$pvalue
pop.padj = resPop$padj

result = data.frame(cbind(pop.LFC,pop.pval,pop.padj))
row.names(result) = row.names(resPop)
head(result)


# write your results into a single output file (for later, for other analyses, etc.)
write.csv(as.data.frame(result), file="DESeq2Summary.csv") 


# let's also have a look at the data in a few different ways ...

# the MA plot shows the log2 fold changes attributable to a given variable over the mean of normalized counts. points with an adjusted p-value < 0.1 are colored red. 
plotMA(res, main="DESeq2", ylim=c(-2,2))


# we may also want to look at expression differences in a specific gene (i.e. one gene at a time). we can call genes by name, but here we will just plot the gene with the greatest expression change (i.e. the minimum p-value)
plotCounts(d, gene=which.min(res$padj), intgroup="pop")

# PCA
rld <- rlog(d, blind=FALSE)
plotPCA(rld, intgroup=c("pop"))
