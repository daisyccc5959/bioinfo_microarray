#WGCNA #17025
if (!requireNamespace("BiocManager", quietly = TRUE))
  BiocManager::install("WGCNA")
# if wgcna cannot install then install these: 
#BiocManager::install(c("GO.db", "preprocessCore", "impute"))
library(WGCNA)

#data acquire: the DE genes expression of each samples(case, control)
library(readr)
GSE17025_exprSet <- read_csv("~/GSE17025_exprSet.csv") #loadin data
GSE17025_exprSet #row: all genes  #col: gene list + patient expression
all_gene<- GSE17025_exprSet$...1
row.names(GSE17025_exprSet)<- all_gene
all_gene<- as.data.frame(all_gene)
all_gene$gene_id= row.names(all_gene)
colnames(all_gene)<- c("gene", "gene_id")

deg_17025 <- read_csv("~/deg_17025.csv") #loadin data
colnames(deg_17025)<- c("gene_id", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "change"   )
deg_17025<- merge(deg_17025, all_gene, by= "gene_id")
write.csv(deg_17025, "C:/Users/admin/Desktop/deg_17025_wname.csv")
deg_17025_wgene

#merge
GSE17025_exprSet$gene_id<- row.names(GSE17025_exprSet)
deg_GSE17025_sample_expr<- merge(GSE17025_exprSet, deg_17025, by= "gene_id")
deg_GSE17025_sample_expr<- deg_GSE17025_sample_expr[,-c(106:113)]

#corrplot have to add script at here!!!####



#network construction and module detection by automatic (2.a)
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
dat_expr<- deg_GSE17025_sample_expr[ ,-c(1:2)]

sft = pickSoftThreshold(dat_expr, powerVector = powers, verbose = 5)
#pickSoftThreshold: will use block size 6
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=1,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=1,col="red")
#network construction for each module
net = blockwiseModules(t(dat_expr), power = 6,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "GSE17025TOM",
                       verbose = 3)
table(net$colors) #there are 3 color of modules from automatic construction
mergedColors = labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)

#3. Relating modules to external clinical traits and identifying important genes
# Define numbers of genes and samples
nGenes = ncol(dat_expr)
nSamples = nrow(dat_expr)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(t(dat_expr), moduleColors)$eigengenes
MEs = orderMEs(MEs0)

#case-control(0101contrast matrix) comparison
patient<- as.data.frame(colnames(GSE17025_exprSet))
patient<- patient[-1,]
patient<- as.data.frame(patient)
patient$group<- rep("case")
patient[92:103,2]<- "control"
#normal
patient[1:103,3]= 0
patient[92:103,3]= 1
#cancer
patient[1:103,4]= 0
patient[1:91,4]= 1
colnames(patient)<- c("patient_id", "group", "normal", "cancer")
contrast_matrix<- patient[ ,3:4]
#character to numeric (class)
normal<- as.numeric(contrast_matrix[,1])
cancer<- as.numeric(contrast_matrix[,2])
contrast_matrix<- cbind(normal, cancer)
row.names(contrast_matrix)<- patient[,1]

moduleTraitCor = cor(MEs, contrast_matrix, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
# Will display correlations and their p-values
sizeGrWindow(10,6)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(contrast_matrix),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

#to find the importance of modules and eigengene
normal = as.data.frame(contrast_matrix[,1]);
names(normal) = "normal"
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(t(dat_expr), MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames=5, sep="");
names(MMPvalue) = paste("p.MM", modNames=5, sep="");

geneTraitSignificance = as.data.frame(cor(dat_expr, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(normal), sep="")
names(GSPvalue) = paste("p.GS.", names(normal), sep="")
#
module = "brown" #change module color
column = match(module, modNames)
moduleGenes = moduleColors==module

sizeGrWindow(7, 7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for case_control",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, pch=20, col = module)
#abs() from module membership may effect the correlation of traits and egiengenes
table(moduleGenes)

#module module comparsion
head(MEs)
# Plot the dendrogram and correlation heatmap
plotEigengeneNetworks(MEs, "Module-Module relationships",
                      excludeGrey = FALSE,  #optional to draw w/o grey
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = TRUE,
                      xLabelsAngle = 90)

#output
#combine with gene name  deg_17025_wgene
gene_module <- data.frame(ID=deg_17025_wgene$gene, module=moduleColors)
gene_module = gene_module[order(gene_module$module),]

#export
#all
write.csv(gene_module,file="~/deg17025_module.csv") # ~ means ur output directory

#others (separate)
M_color<- unique(gene_module$module)
for (i in 1:length(M_color)){
  write.csv(subset(gene_module, module==M_color[i]),
            paste("~/deg17025_module_", # ~ means ur output directory
                  M_color[i],".csv")
            )
}

#2-TOM
#load matrix
dissTOM = 1-TOMsimilarityFromExpr(t(dat_expr), power = 6)
plotTOM = dissTOM^7
diag(plotTOM) = NA
TOMplot(plotTOM, net$dendrograms, moduleColors, 
        main = "Network heatmap plot, all genes")
