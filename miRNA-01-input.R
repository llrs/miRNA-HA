# Read the files

library("WGCNA")
enableWGCNAThreads()
pheno <- read.csv("miRNA_phenotype.csv")
exprs <- read.csv("miRNA_normQ.csv", row.names = 1)
vclin <- read.csv("DB_MIR-Alcoholic_Hepatitis.csv")

AH <- paste0("HA-", c(23, 24, 30, 36, 39, 51, 54, 57, 60, 77, 91, 84, 10, 90,
                      18, 21, 25, 34, 40, 50, 2, 15, 26, 28, 33, 37, 52, 86,
                      87))
normal <- paste0("NQ-", c(2:5))
ids <- c(AH, normal)

pheno <- cbind(pheno, "samplename" = substring(pheno$sampleid, 8),
               "group" = as.factor(substring(pheno$GROUP, 1, 2)),
               "patientid" = ids)

vclin <- merge(pheno, vclin, by.x = "patientid", by.y = "id")

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
