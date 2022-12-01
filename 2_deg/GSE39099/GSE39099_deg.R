#install packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GEOquery")
BiocManager::install("hgu133plus2.db")
BiocManager::install("edgeR")

#BiocManager::install("hgu133plus2.db")
#utils::install.packages('hgu133plus2.db', repos=NULL)
#-->Error in utils::install.packages("hgu133plus2.db", repos = NULL) : 
#     type == "both" cannot be used with 'repos = NULL'
#-->install for .tar in document or install package in R console(bypass R studio)


#data input
setwd("~") #set ur directory
GSE39099_matrix <- read.csv("~GSE39099_series_matrix.txt") #loadin data

#delete series and characteristics
GSE39099_matrix <- GSE39099_matrix[-(1:61),]
a<-unlist(strsplit(as.character(GSE39099_matrix), "\t"))
a<- t(as.data.frame(a))
colnames(a)<-a[1,]
rownames(a)<-a[,1]
a=a[-1,]
GSE39099_expr_probes=a[,-1]
rm(a)
rm(GSE39099_matrix)       #GSE39099_expr_probes: row:probes, col:sample

#probes referred to genes

#array platform: 
#GPL570	[HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array
library('hgu133plus2.db')
library(GEOquery)

#info of platform id to symbol(genes)
hgu133plus2_dbschema()
hgu133plus2_dbInfo()
xx <- as.list(hgu133plus2ALIAS2PROBE)

#totable: probes_id to symbol matrix
#symbol: mappings between manufacturer identifiers and gene abbreviations.
ID = toTable(hgu133plus2SYMBOL)
length(unique(ID$symbol)) #total gene 20824

# take unique gene list and reads frequency
ID_symbol<-sort(table(ID$symbol)) 
table(sort(table(ID$symbol)))
count_density<- as.data.frame(sort(table(ID$symbol))) 
library(ggplot2)
ggplot(count_density, aes(x = Freq)) +
geom_histogram(fill="black", colour="black", alpha = 0.25, binwidth=0.5) + 
geom_density(colour="dark red", adjust=4) 

ID_symbol=as.data.frame(ID_symbol,decreasing = TRUE)

#join GSE39099_expr_probes`probes_id` and ID`probes_id`
exprSet = GSE39099_expr_probes
table(rownames(exprSet) %in% ID$probe_id)
#FALSE  TRUE 
#11575 43101
exprSet=exprSet[rownames(exprSet) %in% ID$probe_id,]
#check_1: nrow(exprSet) [1] 43101
ID=ID[match(rownames(exprSet),ID$probe_id),]
head(ID)

#change into gene name
exprSet<-aggregate(exprSet, by=list(ID$symbol),FUN=max)
#having maximum based on : Microarray background correction: maximum likelihood estimation for the normal–exponential convolution
dim(exprSet)
#table: 20824(gene)*4(case)(the first row is gene name)
rownames(exprSet) = exprSet[,1]
exprSet = exprSet[,-1]
#all genes expression level in GSE39099
write.csv(exprSet, file="GSE39099_exprSet.csv"
          , row.name= TRUE, quote= FALSE)
#rm(ID) rm(ID_symbol) rm(GSE39099_expr_probes) rm(xx)

#change data structure
#is character
genes <- rownames(exprSet)
exprSet = as.data.frame(exprSet ,stringsAsFactors=FALSE)
exprSet = as.data.frame(sapply(exprSet, as.numeric))
sapply(exprSet,summary)
exprSet <- cbind(genes,exprSet)

#classification n clustering
expr_h_data<-t(exprSet[,-1])
E.dist <- dist(expr_h_data, method="euclidean")
M.dist <- dist(expr_h_data, method="manhattan")

# E
h.E.cluster <- hclust(E.dist) #default of hclust is complete
plot(h.E.cluster, xlab="euclidean")

# E + ward
h.E.cluster1 <- hclust(E.dist,method = "ward.D2")
plot(h.E.cluster1, xlab="euclidean")

# M
h.M.cluster <- hclust(M.dist) 
plot(h.M.cluster, xlab="manhattan")

# M + ward
h.M.cluster1 <- hclust(M.dist,method = "ward.D2")
plot(h.M.cluster1, xlab="manhattan")

#expression pattern b/w cases
#untransform data just for an overview
boxplot(exprSet[,2:5])

#RMA background correction: background correction, quantile normalization, and log2 transformation
#library(affy)
#exprSet <- rma2(Dilution)
#for it doesn't work, just with log2 transformation
lo2_exprSet<-log2(exprSet[,-1])

#exprssion distribution
boxplot(lo2_exprSet[,1:4])
library(reshape2)
melt_log<- melt(cbind(genes,lo2_exprSet))
colnames(melt_log)=c('genes','sample','value')
#plot
ggplot(melt_log,aes(x=sample,y=value,fill=sample))+ geom_boxplot() + scale_fill_brewer(palette="Pastel1")

ggplot(melt_log,aes(x=sample,y=value,fill=sample))+geom_violin() + scale_fill_brewer(palette="Pastel1")

ggplot(melt_log,aes(value,fill=sample))+ geom_histogram(bins = 200) + facet_wrap(~sample, nrow = 4) + scale_fill_brewer(palette="Pastel1")

ggplot(melt_log,aes(value))+ geom_histogram(bins = 100) 
plotDensities(melt_log$value)

rm(melt_log)




#for this sample w/o replicates
#followed with edgeR 2.12 to estimate DEGs (DEG List)

#method 1: extactTest
group<- 1:4
y<- DGEList(counts=exprSet[,3:6], group=(1:4))
y<- calcNormFactors(y)
#human recommend bcv(gene dispersion) = 0.4 (in document)
bcv<- 0.4
m1_exacttest<- exactTest(y, dispersion=bcv^2)
deg_test_1<- decideTestsDGE(m1_exacttest, p.value=0.05, lfc=0)
summary(deg_test_1)
#         2-1  toomany...
#Down    2700
#NotSig 15817
#Up      2307

#method 2: glmfit
#estimating dispersion value using control (GSM956133)
y1<- y #copy
y1$samples$group <- 1
#assume control is not DEGs--> to estimate dispersion value
y0 <- estimateDisp(y1[2,], trend="none", tagwise=FALSE)
y$common.dispersion <- y0$common.dispersion
design<- model.matrix(~group)
fit <- glmFit(y, design)
lrt <- glmLRT(fit)
deg_test_2<- decideTestsDGE(lrt, p.value=0.05, lfc=0)
summary(deg_test_2)
#       group
#Down     231
#NotSig 20384
#Up       209

#compare: venn plot
head(deg_test_1@.Data) #data:0, 1, -1
#up
up_gene1= names(deg_test_1@.Data[deg_test_1@.Data==1,])
up_gene2= names(deg_test_2@.Data[deg_test_2@.Data==1,])
a<-list(extactTest=up_gene1, glmfit=up_gene2)

#down
d_gene1= names(deg_test_1@.Data[deg_test_1@.Data==-1,])
d_gene2= names(deg_test_2@.Data[deg_test_2@.Data==-1,])
a1<- list(extactTest=d_gene1, glmfit=d_gene2)

library(gplots)
venn(a)
venn(a1)

#mapping to gene list and volcano plot (just use m2)
deg<- lrt$table
deg_test_2<- as.data.frame(deg_test_2)
deg_test_2$filter<- row.names(deg_test_2)
deg= cbind(deg_test_2, genes, deg)
#the output of m2
deg_filter<- subset(deg, group!=0)
#deg
write.csv(deg_filter, "deg39099_filter.csv")
ggplot(data= deg, aes(x= logFC, y= -log10(PValue)))+
  geom_point(alpha=0.4)+
  xlab("log2 fold change") + ylab("-log10 p-value")+
  theme_bw()

deg_probengenes<-cbind(genes,probe_id,deg)

logFC_cutoff <- with(deg_probengenes, mean(abs(logFC)) + 2*sd(abs(logFC)) )
deg_probengenes$change = as.factor(ifelse(deg_probengenes$PValue < 0.05 &
                              abs(deg_probengenes$logFC) > logFC_cutoff,
                                        ifelse(deg_probengenes$logFC > logFC_cutoff ,'UP'
                                               ,'DOWN')
                              ,'NOT'))

ggplot(data=deg_probengenes, aes(x= logFC, y= -log10(PValue),color=change))+
  geom_point(alpha=0.4, size=1.2) +
  theme_set(theme_set(theme_bw(base_size=20)))+
  labs(title = "vocalno plot", x= "log2 fold change", y= "-log10 p-value") +
  theme(plot.title = element_text(size=15,hjust = 0.5))+
  scale_colour_manual(values = c('blue','black','red'))+
  geom_hline(yintercept = -log10(0.05), lty= 4, lwd= 0.6, alpha= 0.8)+
  geom_vline(xintercept = c(1,-1), lty= 4, lwd= 0.6, alpha= 0.8)+
  theme_bw()+theme(panel.border = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   axis.line= element_line(color="black"))
#all gene
write.csv(deg_probengenes,'deg_39099_all.csv')
