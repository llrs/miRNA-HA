# Read the files

library("WGCNA")
enableWGCNAThreads()
pheno <- read.csv("miRNA_phenotype.csv")
exprs <- read.csv("miRNA_normQ.csv", row.names = 1)

pheno <- cbind(pheno, "samplename" = substring(pheno$sampleid, 8),
               "group" = as.factor(substring(pheno$GROUP, 1, 2)))

data.wgcna <- t(exprs[, pheno$group  == "AH"])
gsg <- goodSamplesGenes(data.wgcna, verbose = 3)
if (!gsg$allOK)
{ # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", 
                     paste(names(data.wgcna)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples) > 0)
    printFlush(paste("Removing samples:", paste(rownames(data.wgcna)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  data.wgcna <- data.wgcna[gsg$goodSamples, gsg$goodGenes]
}


pdf("samples_dendrogram.pdf")
sampleTree <- hclust(dist(data.wgcna), method = "average")
pars <- par(mar = c(0, 4, 2, 0), cex = 0.6)
plot(sampleTree, main = "Sample clustering to detect outliers", 
     sub = "", xlab = "", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()

save(pheno, exprs, file = "Input.RData")
