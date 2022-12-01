# install package
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GEOquery")
BiocManager::install("hgu133plus2.db")
BiocManager::install("limma")

#BiocManager::install("hgu133plus2.db")
#utils::install.packages('hgu133plus2.db', repos=NULL)
#-->Error in utils::install.packages("hgu133plus2.db", repos = NULL) : 
#     type == "both" cannot be used with 'repos = NULL'
#-->install for .tar in document or install package in R console(bypass R studio)
#data input

setwd("~") #set ur directory
GSE17025_matrix <- read.csv("~GSE17025_series_matrix.txt")  #loadin data

#delete series and characteristics
GSE17025_matrix <- GSE17025_matrix[-(1:83),]
a<-strsplit(GSE17025_matrix, "\t")
a<- t(as.data.frame(a))
rownames(a)<-a[,1]
a=a[,-1]
GSE17025_expr_probes<-a #GSE17025_expr_probes: row:probes, col:sample
rm(a)
rm(GSE17025_matrix)       

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

#join GSE17025_expr_probes`probes_id` and ID`probes_id`
exprSet_17025 = GSE17025_expr_probes
table(rownames(exprSet_17025) %in% ID$probe_id)
#FALSE  TRUE 
#11575 43101
exprSet_17025=exprSet_17025[rownames(exprSet_17025) %in% ID$probe_id,]
#check_1: nrow(exprSet_17025) [1] 43101
ID=ID[match(rownames(exprSet_17025),ID$probe_id),]
head(ID)

#change into gene name
exprSet_17025<-aggregate(exprSet_17025, by=list(ID$symbol),FUN=max)
#having maximum based on : Microarray background correction: maximum likelihood estimation for the normal–exponential convolution
dim(exprSet_17025)
#table: 20824(gene)*4(case)(the first row is gene name)
rownames(exprSet_17025) = exprSet_17025[,1]
exprSet_17025 = exprSet_17025[,-1]
#all genes expression level in GSE17025
write.csv(exprSet_17025, file="GSE17025_exprSet.csv"
          , row.name= TRUE, quote= FALSE)
#rm(ID) rm(ID_symbol) rm(GSE17025_expr_probes) rm(xx)

#change data structure
#is character
patient_id<- colnames(exprSet_17025)
genes <- rownames(exprSet_17025)
exprSet_17025 = as.data.frame(exprSet_17025 ,stringsAsFactors=FALSE)
exprSet_17025 = as.data.frame(sapply(exprSet_17025, as.numeric))
sapply(exprSet_17025,summary)
exprSet_17025 <- cbind(genes,exprSet_17025)

#RMA background correction: background correction, quantile normalization, and log2 transformation
#library(affy)
#exprSet <- rma2(Dilution)
#for it doesn't work, just with log2 transformation
lo2_exprSet_17025<-log2(exprSet_17025[,-1])
#exprssion distribution
boxplot(lo2_exprSet_17025[,1:ncol(lo2_exprSet_17025)])

#add group
group<-as.data.frame(patient_id)
group$group<-rep("case")
group[,92:103]="control"

#DEG by limma ##single factor design(case n control)
library(limma)
##create case-control matrix and contrast matrix with function
#model matrix
matrix_1 <- model.matrix(~0+factor(group$group))
colnames(matrix_1)=levels(factor(group$group))
rownames(matrix_1)=colnames(lo2_exprSet_17025)
matrix_1

#contrast matrix (case-control)
contrast_matrix<- matrix(c(1,-1), ncol=1)
dimnames(contrast_matrix)<- list(c("case","control"),"diff")
#contrast_matrix
contrast_matrix<- makeContrasts(Diff= case - control,levels= matrix_1)

#fit limma
fit<-lmFit(lo2_exprSet_17025,matrix_1)
fit2<-contrasts.fit(fit,contrast_matrix)
fit3<-eBayes(fit2)

#all data just for volcano plot
deg_17025_all<- topTable(fit3,genelist=fit$genes,number = nrow(lo2_exprSet_17025))
#the column name is the id name of input, for gene name not unique
#in (https://rdrr.io/bioc/limma/man/toptable.html)
logFC_cutoff <- with(deg_17025_all,mean(abs( logFC)) + 2*sd(abs( logFC)) )
deg_17025_all$change = as.factor(ifelse(deg_17025_all$P.Value < 0.05 & abs(deg_17025_all$logFC) > logFC_cutoff,
                              ifelse(deg_17025_all$logFC > logFC_cutoff ,'UP','DOWN'),'NOT'))
ggplot(data=deg_17025_all, aes(x= logFC, y= -log10(P.Value),color=change))+
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

#data output deg 
#set: p-value=0.05, fdr, log2 fold change=2
deg_17025<-topTable(fit3, coef="Diff",p.value=0.05, adjust.method="fdr", lfc=2,number = nrow(lo2_exprSet_17025))
deg_17025$change = as.factor(ifelse(deg_17025$P.Value < 0.05 & abs(deg_17025$logFC) > logFC_cutoff,
                                        ifelse(deg_17025$logFC > logFC_cutoff ,'UP','DOWN'),'NOT'))
write.csv(deg_17025,"deg_17025.csv")
