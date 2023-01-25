##############
# DE_analyse #
##############

library(edgeR)
library(limma)
library(biomaRt)
library(plotrix)
library(BSDA)

##############
# DE_analyse #
##############


setwd("C:/Users/Louis Coussement/OneDrive/Louis/AAP/Hoofdstuk_Wim/")

counts_GEO <- read.table(file="countsGEO.txt",header=T,row.names = 1,stringsAsFactors = F)

counts <- read.table(file="counts.txt",header=F,row.names = 1,stringsAsFactors = F)
counts <- counts[-((nrow(counts)-4):nrow(counts)),]

counts_exon <- read.table(file="counts_exon.txt",header=F,row.names = 1,stringsAsFactors = F)
counts_exon <- counts_exon[-((nrow(counts_exon)-4):nrow(counts_exon)),]

counts_trim <- read.table(file="counts_trim_bis.txt",header=F,row.names = 1,stringsAsFactors = F)
counts_trim <- counts_trim[-((nrow(counts_trim)-4):nrow(counts_trim)),]

counts_exon_trim <- read.table(file="counts_exon_trim_bis.txt",header=F,row.names = 1,stringsAsFactors = F)
counts_exon_trim <- counts_exon_trim[-((nrow(counts_exon_trim)-4):nrow(counts_exon_trim)),]

samp1_bis <- read.table("test.txt",header=F,row.names = 1)
samp1_bis <- samp1_bis[-((nrow(samp1_bis)-4):nrow(samp1_bis)),]

not_equal <- cbind(counts_trim[counts_trim[,1]!=samp1_bis,1],samp1_bis[counts_trim[,1]!=samp1_bis])
rownames(not_equal) <- rownames(counts_trim)[counts_trim[,1]!=samp1_bis]
rRNA <- read.table(file="rrna.txt",header=T,stringsAsFactors = F,sep="\t")

colnames(counts) <- colnames(counts_GEO)
colnames(counts_exon) <- colnames(counts_GEO)
colnames(counts_trim) <- colnames(counts_GEO)
colnames(counts_exon_trim) <- colnames(counts_GEO)


dim(counts_GEO)
dim(counts)
dim(counts_exon)
head(counts)
head(counts_exon)

dim(counts_trim)
dim(counts_exon_trim)
head(counts_trim)
head(counts_exon_trim)

colSums(counts_GEO)
colSums(counts)
colSums(counts_exon)

colSums(counts)/colSums(counts_GEO)
colSums(counts_exon)/colSums(counts_GEO)

# counts[sort(rowMeans(counts),index.return=T,decreasing = T)$ix[1:15],]
# counts_GEO[sort(rowMeans(counts_GEO),index.return=T,decreasing = T)$ix[1:10],]
# which(names(head(sort(rowMeans(counts),decreasing = T)))%in%rRNA$Gene.stable.ID)

dim(counts_trim)
counts_trim <- counts_trim[!(rownames(counts_trim)%in%rRNA$Gene.stable.ID),]
dim(counts_trim)


# object_olivier <- DGEList(counts_trim)
# object_olivier <- calcNormFactors(object_olivier)
# counts_dge <- object_olivier
# save(counts_dge,file="Data_DGElist.Rda")
# object_olivier <- cpm(object_olivier,log=T)
# counts_norm <- object_olivier
# save(counts_norm,list=c("counts_norm"),file="Data_normCounts.Rda")

## Get annotation
if(!file.exists("Genes.Rda")){
  ensembl=useMart("ensembl")
  ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
  genes <- getBM(attributes = c('ensembl_gene_id','hgnc_symbol','chromosome_name','start_position','end_position','transcription_start_site'),mart = ensembl)
  print(dim(genes))
  save(genes,file="Genes.Rda")
} else {
  load("Genes.Rda")
}


## Make annotation for design later on
individuals <- factor(unlist(strsplit(colnames(counts),"\\_"))[1:44*3-1])
cvst <- factor(unlist(strsplit(colnames(counts),"\\_"))[1:44*3])

cutoff_filter <- 0.2


####################
####################
####            ####
#### Gene-level ####
####            ####
####################
####################

###############
##           ##
## Norm: TMM ##
##           ##
###############

##################
## Non filtered ##
##################

## Make DGElist
counts_TMM_nf <- counts_trim
y_TMM_nf <- DGEList(counts_TMM_nf)
y_TMM_nf <- calcNormFactors(y_TMM_nf)

## Density-plots
jpeg("density_nf.jpg",width=1000,height=1000)
for (i in 1:ncol(counts_trim)){
  if (cvst[i]=="normal"){
    color="blue"
  } else {
    color="orange"
  }
  if (i==1){
    plot(density(log(y_TMM_nf$counts[,i])),ylim=c(0,max(density(log(y_TMM_nf$counts[,i]))$y)*2.1),xlim=c(-1,12),col=color,
         xlab="Expression (log)",main="Density plot")
  } else {
    lines(density(log(y_TMM_nf$counts[,i])),col=color)
  }
}
dev.off()

## MA-plots
jpeg("MAplot_nf.jpg")
plot(log2(rowMeans(y_TMM_nf$counts)),log2(rowMeans(y_TMM_nf$counts[,1:22])/rowMeans(y_TMM_nf$counts[,23:44])),pch=".",main="MA plot",xlab="A",ylab="M")
dev.off()

## Density-plots normalized
countspm_TMM_nf <- cpm(y_TMM_nf,log=T)
jpeg("density_TMM_nf.jpg",width=1000,height=1000)
for (i in 1:ncol(counts_trim)){
  if (i==23){
    next
  }
  if (i==1){
    plot(density(countspm_TMM_nf[,i]),ylim=c(0,max(density(log(counts_trim[,i]))$y)*8),
         xlab="Expression (log)",main="Density plot")
  } else {
    lines(density(countspm_TMM_nf[,i]))
  }
}
dev.off()

## MA-plots
jpeg("MAplot_TMM_nf.jpg")
plot(rowMeans(countspm_TMM_nf),rowMeans(countspm_TMM_nf[,23:44])-rowMeans(countspm_TMM_nf[,1:22]),pch=".",main="MA plot",xlab="A",ylab="M")
dev.off()



## Unpaired analysis
##############

## Make Design
design_TMM_nf <- model.matrix(~cvst)
rownames(design_TMM_nf) <- colnames(y_TMM_nf)

## Build Model
jpeg("limmatrend_TMM_nf.jpg")
v_TMM_nf <- voomWithQualityWeights(y_TMM_nf, design_TMM_nf,plot=T)
dev.off()
fit_TMM_nf <- lmFit(v_TMM_nf, design_TMM_nf)
fit_TMM_nf <- eBayes(fit_TMM_nf)
res_TMM_nf <- topTable(fit_TMM_nf, coef="cvsttumor", n=nrow(v_TMM_nf))
# coef="cvsttumor" => tumor vs. controle, see "design" object 

## MAplot
jpeg("MAplot_TMM_res_nf.jpg")
# limma::plotMA(fit_TMM_nf)
# MAplot: all data points
plot(rowSums(v_TMM_nf$E)[rownames(res_TMM_nf)],
       res_TMM_nf$logFC,pch=16,cex=0.6)
# MA-plot: significant loci
points(rowSums(v_TMM_nf$E)[rownames(res_TMM_nf)[res_TMM_nf$adj.P.Val<0.05]],
       res_TMM_nf$logFC[res_TMM_nf$adj.P.Val<0.05],pch=16,col="red",cex=0.6)
# X-axis
abline(0,0) 
dev.off()

## Volcano plot
jpeg("Volcano_TMM_res_nf.jpg")
plot(res_TMM_nf$logFC,-log10(res_TMM_nf$adj.P.Val),pch=16,cex=0.2,xlab="logFC",ylab="-log10(FDR)") 
points(res_TMM_nf$logFC[abs(res_TMM_nf$logFC)>1 & -log10(res_TMM_nf$adj.P.Val)>4],
       -log10(res_TMM_nf$adj.P.Val)[abs(res_TMM_nf$logFC)>1 & -log10(res_TMM_nf$adj.P.Val)>4],pch=16,col="red",cex=0.6)
title("Volcano plot")
dev.off()

## P-value distribution
jpeg("Pdistr_TMM_res_nf.jpg")
hist(res_TMM_nf$P.Value,main="Unfiltered Data",xlab="P-values")
dev.off()


## Paired analysis
##############

## Make Design
designP_TMM_nf <- model.matrix(~individuals+cvst) #individuals as well
rownames(designP_TMM_nf) <- colnames(y_TMM_nf)

## Build Model
jpeg("limmatrendP_TMM_nf.jpg")
vP_TMM_nf <- voomWithQualityWeights(y_TMM_nf, designP_TMM_nf, plot=T)
dev.off()
fitP_TMM_nf <- lmFit(vP_TMM_nf, designP_TMM_nf)
fitP_TMM_nf <- eBayes(fitP_TMM_nf)
resP_TMM_nf <- topTable(fitP_TMM_nf, coef="cvsttumor", n=nrow(vP_TMM_nf))
# coef="cvsttumor" => tumor vs. control 

## MAplot
jpeg("MAplotP_TMM_res_nf.jpg")
# limma::plotMA(fit_TMM_nf)
# MAplot: all data points
plot(rowSums(vP_TMM_nf$E)[rownames(resP_TMM_nf)],
     resP_TMM_nf$logFC,pch=16,cex=0.6)
# MA-plot: significant loci
points(rowSums(vP_TMM_nf$E)[rownames(resP_TMM_nf)[resP_TMM_nf$adj.P.Val<0.05]],
       resP_TMM_nf$logFC[resP_TMM_nf$adj.P.Val<0.05],pch=16,col="red",cex=0.6)
# X-axis
abline(0,0) 
dev.off()

## Volcano plot
jpeg("VolcanoP_TMM_res_nf.jpg")
plot(resP_TMM_nf$logFC,-log10(resP_TMM_nf$adj.P.Val),pch=16,cex=0.2,xlab="logFC",ylab="-log10(FDR)") 
points(resP_TMM_nf$logFC[abs(resP_TMM_nf$logFC)>1 & -log10(resP_TMM_nf$adj.P.Val)>4],
       -log10(resP_TMM_nf$adj.P.Val)[abs(resP_TMM_nf$logFC)>1 & -log10(resP_TMM_nf$adj.P.Val)>4],pch=16,col="red",cex=0.6)
title("Volcano plot")
dev.off()

## P-value distribution
jpeg("PdistrP_TMM_res_nf.jpg")
hist(resP_TMM_nf$P.Value,main="Unfiltered Data",xlab="P-values")
dev.off()


##############
## Filtered ##
##############

## Make DGElist
y_TMM_f <- DGEList(counts_TMM_nf)
y_TMM_f <- calcNormFactors(y_TMM_f)

## Filtering
countspm_TMM_f <- cpm(y_TMM_f,log=T)
counts_TMM_f <- counts_trim[rowMeans(countspm_TMM_f)>=cutoff_filter,]
y_TMM_f <- DGEList(counts_TMM_f)
y_TMM_f <- calcNormFactors(y_TMM_f)

## Density-plots
jpeg("density_f.jpg",width=1000,height=1000)
for (i in 1:ncol(counts_trim)){
  if (cvst[i]=="normal"){
    color="blue"
  } else {
    color="orange"
  }
  if (i==1){
    plot(density(log(y_TMM_f$counts[,i])),ylim=c(0,max(density(log(y_TMM_f$counts[,i]))$y)*2.1),xlim=c(-1,12),col=color,
         xlab="Expression (log)",main="Density plot")
  } else {
    lines(density(log(y_TMM_f$counts[,i])),col=color)
  }
}
dev.off()

## MA-plots
jpeg("MAplot_f.jpg")
plot(log2(rowMeans(y_TMM_f$counts)),log2(rowMeans(y_TMM_f$counts[,1:22])/rowMeans(y_TMM_f$counts[,23:44])),pch=".",main="MA plot",xlab="A",ylab="M")
dev.off()

## Density-plots normalized
countspm_TMM_f <- cpm(y_TMM_f,log=T)
jpeg("density_TMM_f.jpg",width=1000,height=1000)
for (i in 1:ncol(counts_trim)){
  if (i==23){
    next
  }
  if (i==1){
    plot(density(countspm_TMM_f[,i]),ylim=c(0,max(density(log(counts_trim[,i]))$y)*1.5),
         xlab="Expression (log)",main="Density plot")
  } else {
    lines(density(countspm_TMM_f[,i]))
  }
}
dev.off()

## MA-plots
jpeg("MAplot_TMM_f.jpg")
plot(rowMeans(countspm_TMM_f),rowMeans(countspm_TMM_f[,23:44])-rowMeans(countspm_TMM_f[,1:22]),pch=".",main="MA plot",xlab="A",ylab="M")
dev.off()


## Unpaired analysis
##############

## Make Design
design_TMM_f <- model.matrix(~cvst)
rownames(design_TMM_f) <- colnames(y_TMM_f)

## Build Model
jpeg("limmatrend_TMM_f.jpg")
v_TMM_f <- voomWithQualityWeights(y_TMM_f, design_TMM_f,plot=T)
dev.off()
fit_TMM_f <- lmFit(v_TMM_f, design_TMM_f)
fit_TMM_f <- eBayes(fit_TMM_f)
res_TMM_f <- topTable(fit_TMM_f, coef="cvsttumor", n=nrow(v_TMM_f))
# coef="cvsttumor" => tumor vs. controle, see "design" object 

## MAplot
jpeg("MAplot_TMM_res_f.jpg")
# limma::plotMA(fit_TMM_nf)
# MAplot: all data points
plot(rowSums(v_TMM_f$E)[rownames(res_TMM_f)],
     res_TMM_f$logFC,pch=16,cex=0.6)
# MA-plot: significant loci
points(rowSums(v_TMM_f$E)[rownames(res_TMM_f)[res_TMM_f$adj.P.Val<0.05]],
       res_TMM_f$logFC[res_TMM_f$adj.P.Val<0.05],pch=16,col="red",cex=0.6)
# X-axis
abline(0,0) 
dev.off()

## Volcano plot
jpeg("Volcano_TMM_res_f.jpg")
plot(res_TMM_f$logFC,-log10(res_TMM_f$adj.P.Val),pch=16,cex=0.2,xlab="logFC",ylab="-log10(FDR)") 
points(res_TMM_f$logFC[abs(res_TMM_f$logFC)>1 & -log10(res_TMM_f$adj.P.Val)>4],
       -log10(res_TMM_f$adj.P.Val)[abs(res_TMM_f$logFC)>1 & -log10(res_TMM_f$adj.P.Val)>4],pch=16,col="red",cex=0.6)
title("Volcano plot")
dev.off()

## P-value distribution
jpeg("Pdistr_TMM_res_f.jpg")
hist(res_TMM_f$P.Value,main="Filtered Data",xlab="P-values")
dev.off()


## Paired analysis
##############

# Make Design
designP_TMM_f <- model.matrix(~individuals+cvst) #individuals as well
rownames(designP_TMM_f) <- colnames(y_TMM_f)

# Build Model
jpeg("limmatrendP_TMM_f.jpg")
vP_TMM_f <- voomWithQualityWeights(y_TMM_f, designP_TMM_f, plot=T)
dev.off()
fitP_TMM_f <- lmFit(vP_TMM_f, designP_TMM_f)
fitP_TMM_f <- eBayes(fitP_TMM_f)
resP_TMM_f <- topTable(fitP_TMM_f, coef="cvsttumor", n=nrow(vP_TMM_f))
# coef="cvsttumor" => tumor vs. control 


## MAplot
jpeg("MAplotP_TMM_res_f.jpg")
# limma::plotMA(fit_TMM_nf)
# MAplot: all data points
plot(rowSums(vP_TMM_f$E)[rownames(resP_TMM_f)],
     resP_TMM_f$logFC,pch=16,cex=0.6)
# MA-plot: significant loci
points(rowSums(vP_TMM_f$E)[rownames(resP_TMM_f)[resP_TMM_f$adj.P.Val<0.05]],
       resP_TMM_f$logFC[resP_TMM_f$adj.P.Val<0.05],pch=16,col="red",cex=0.6)
# X-axis
abline(0,0) 
dev.off()

## Volcano plot
jpeg("VolcanoP_TMM_res_f.jpg")
plot(resP_TMM_f$logFC,-log10(resP_TMM_f$adj.P.Val),pch=16,cex=0.2,xlab="logFC",ylab="-log10(FDR)") 
points(resP_TMM_f$logFC[abs(resP_TMM_f$logFC)>1 & -log10(resP_TMM_f$adj.P.Val)>4],
       -log10(resP_TMM_f$adj.P.Val)[abs(resP_TMM_f$logFC)>1 & -log10(resP_TMM_f$adj.P.Val)>4],pch=16,col="red",cex=0.6)
title("Volcano plot")
dev.off()

## P-value distribution
jpeg("PdistrP_TMM_res_f.jpg")
hist(resP_TMM_f$P.Value,main="Filtered Data",xlab="P-values")
dev.off()



####################
##                ##
## Norm: Quantile ##
##                ##
####################


##################
## Non filtered ##
##################

## Make DGElist
counts_quant_nf <- counts_trim
y_quant_nf <- DGEList(counts_quant_nf)
y_quant_nf <- calcNormFactors(y_quant_nf)

## Make Design
design_quant_nf <- model.matrix(~cvst)
rownames(design_quant_nf) <- colnames(y_quant_nf)

## Density-plots normalized: quantile
countspm_quant_nf <- voomWithQualityWeights(y_quant_nf, design_quant_nf,normalize.method = "quantile")$E
jpeg("density_quant_nf.jpg",width=1000,height=1000)
for (i in 1:ncol(counts_trim)){
  if (colSums(counts_trim)[i]<10000000){
    color <- "red"
  } else {
    color <- "black"
  }
  if (i==1){
    plot(density(countspm_quant_nf[,i]),col=color,xlab="Expression (log)",main="Density plot")
  } else {
    lines(density(countspm_quant_nf[,i]),col=color)
  }
}
dev.off()

## MA-plots
jpeg("MAplot_quant_nf.jpg")
plot(rowMeans(countspm_quant_nf),rowMeans(countspm_quant_nf[,23:44])-rowMeans(countspm_quant_nf[,1:22]),pch=".",main="MA plot",xlab="A",ylab="M")
dev.off()

# ## Object Olivier
# counts_nf <- counts_quant_nf
# save(counts_nf,file="Data_DGElist_nf.Rda")
# object_olivier <-voomWithQualityWeights(y_quant_nf, design_quant_nf,normalize.method = "quantile")$E
# save(object_olivier,list=c("counts_quant_nf"),file="Data_normCounts_nf.Rda")

## Unpaired analysis
##############

## Build Model
jpeg("limmatrend_quant_nf.jpg")
v_quant_nf <- voomWithQualityWeights(y_quant_nf, design_quant_nf,plot=T,normalize.method = "quantile")
dev.off()
fit_quant_nf <- lmFit(v_quant_nf, design_quant_nf)
fit_quant_nf <- eBayes(fit_quant_nf)
res_quant_nf <- topTable(fit_quant_nf, coef="cvsttumor", n=nrow(v_quant_nf))
# coef="cvsttumor" => tumor vs. controle, see "design" object 

## MAplot
jpeg("MAplot_quant_res_nf.jpg")
# limma::plotMA(fit_quant_nf)
# MAplot: all data points
plot(rowSums(v_quant_nf$E)[rownames(res_quant_nf)],
     res_quant_nf$logFC,pch=16,cex=0.6)
# MA-plot: significant loci
points(rowSums(v_quant_nf$E)[rownames(res_quant_nf)[res_quant_nf$adj.P.Val<0.05]],
       res_quant_nf$logFC[res_quant_nf$adj.P.Val<0.05],pch=16,col="red",cex=0.6)
# X-axis
abline(0,0) 
dev.off()

## Volcano plot
jpeg("Volcano_quant_res_nf.jpg")
plot(res_quant_nf$logFC,-log10(res_quant_nf$adj.P.Val),pch=16,cex=0.2,xlab="logFC",ylab="-log10(FDR)") 
points(res_quant_nf$logFC[abs(res_quant_nf$logFC)>1 & -log10(res_quant_nf$adj.P.Val)>4],
       -log10(res_quant_nf$adj.P.Val)[abs(res_quant_nf$logFC)>1 & -log10(res_quant_nf$adj.P.Val)>4],pch=16,col="red",cex=0.6)
title("Volcano plot")
dev.off()

## P-value distribution
jpeg("Pdistr_quant_res_nf.jpg")
hist(res_quant_nf$P.Value,main="Unfiltered Data",xlab="P-values")
dev.off()


## Paired analysis
##############

## Make Design
designP_quant_nf <- model.matrix(~individuals+cvst) #individuals as well
rownames(designP_quant_nf) <- colnames(y_quant_nf)

## Build Model
jpeg("limmatrendP_quant_nf.jpg")
vP_quant_nf <- voomWithQualityWeights(y_quant_nf, designP_quant_nf, plot=T,normalize.method = "quantile")
dev.off()
fitP_quant_nf <- lmFit(vP_quant_nf, designP_quant_nf)
fitP_quant_nf <- eBayes(fitP_quant_nf)
resP_quant_nf <- topTable(fitP_quant_nf, coef="cvsttumor", n=nrow(vP_quant_nf))
# coef="cvsttumor" => tumor vs. control 

## MAplot
jpeg("MAplotP_quant_res_nf.jpg")
# limma::plotMA(fit_quant_nf)
# MAplot: all data points
plot(rowSums(vP_quant_nf$E)[rownames(resP_quant_nf)],
     resP_quant_nf$logFC,pch=16,cex=0.6)
# MA-plot: significant loci
points(rowSums(vP_quant_nf$E)[rownames(resP_quant_nf)[resP_quant_nf$adj.P.Val<0.05]],
       resP_quant_nf$logFC[resP_quant_nf$adj.P.Val<0.05],pch=16,col="red",cex=0.6)
# X-axis
abline(0,0) 
dev.off()

## Volcano plot
jpeg("VolcanoP_quant_res_nf.jpg")
plot(resP_quant_nf$logFC,-log10(resP_quant_nf$adj.P.Val),pch=16,cex=0.2,xlab="logFC",ylab="-log10(FDR)") 
points(resP_quant_nf$logFC[abs(resP_quant_nf$logFC)>1 & -log10(resP_quant_nf$adj.P.Val)>4],
       -log10(resP_quant_nf$adj.P.Val)[abs(resP_quant_nf$logFC)>1 & -log10(resP_quant_nf$adj.P.Val)>4],pch=16,col="red",cex=0.6)
title("Volcano plot")
dev.off()

## P-value distribution
jpeg("PdistrP_quant_res_nf.jpg")
hist(resP_quant_nf$P.Value,main="Unfiltered Data",xlab="P-values")
dev.off()



##############
## Filtered ##
##############

## Make DGElist
y_quant_f <- DGEList(counts_quant_nf)
y_quant_f <- calcNormFactors(y_quant_f)

## Filtering
countspm_quant_f <- cpm(y_quant_f,log=T)
counts_quant_f <- counts_trim[rowMeans(countspm_quant_f)>=cutoff_filter,]
y_quant_f <- DGEList(counts_quant_f)
y_quant_f <- calcNormFactors(y_quant_f)

## Make Design
design_quant_f <- model.matrix(~cvst)
rownames(design_quant_f) <- colnames(y_quant_f)

## Density-plots normalized: quantile
countspm_quant_f <- voomWithQualityWeights(y_quant_f, design_quant_f,normalize.method = "quantile")$E
jpeg("density_quant_f.jpg",width=1000,height=1000)
for (i in 1:ncol(counts_trim)){
  if (colSums(counts_trim)[i]<10000000){
    color <- "red"
  } else {
    color <- "black"
  }
  if (i==1){
    plot(density(countspm_quant_f[,i]),col=color,ylim=c(0,max(density(countspm_quant_f[,i])$y)*2),
         xlab="Expression (log)",main="Density plot")
  } else {
    lines(density(countspm_quant_f[,i]),col=color)
  }
}
dev.off()

## MA-plots
jpeg("MAplot_quant_f.jpg")
plot(rowMeans(countspm_quant_f),rowMeans(countspm_quant_f[,23:44])-rowMeans(countspm_quant_f[,1:22]),pch=".",main="MA plot",xlab="A",ylab="M")
dev.off()

# ## Object Olivier
# counts_nf <- counts_quant_f
# save(counts_nf,file="Data_DGElist_f.Rda")
# object_olivier <-voomWithQualityWeights(y_quant_f, design_quant_f,normalize.method = "quantile")$E
# save(object_olivier,list=c("counts_quant_f"),file="Data_normCounts_f.Rda")


## Unpaired analysis
##############

## Build Model
jpeg("limmatrend_quant_f.jpg")
v_quant_f <- voomWithQualityWeights(y_quant_f, design_quant_f,plot=T,normalize.method = "quantile")
dev.off()
fit_quant_f <- lmFit(v_quant_f, design_quant_f)
fit_quant_f <- eBayes(fit_quant_f)
res_quant_f <- topTable(fit_quant_f, coef="cvsttumor", n=nrow(v_quant_f))
# coef="cvsttumor" => tumor vs. controle, see "design" object 

## MAplot
jpeg("MAplot_quant_res_f.jpg")
# limma::plotMA(fit_quant_f)
# MAplot: all data points
plot(rowSums(v_quant_f$E)[rownames(res_quant_f)],
     res_quant_f$logFC,pch=16,cex=0.6)
# MA-plot: significant loci
points(rowSums(v_quant_f$E)[rownames(res_quant_f)[res_quant_f$adj.P.Val<0.05]],
       res_quant_f$logFC[res_quant_f$adj.P.Val<0.05],pch=16,col="red",cex=0.6)
# X-axis
abline(0,0) 
dev.off()

## Volcano plot
jpeg("Volcano_quant_res_f.jpg")
plot(res_quant_f$logFC,-log10(res_quant_f$adj.P.Val),pch=16,cex=0.2,xlab="logFC",ylab="-log10(FDR)") 
points(res_quant_f$logFC[abs(res_quant_f$logFC)>1 & -log10(res_quant_f$adj.P.Val)>4],
       -log10(res_quant_f$adj.P.Val)[abs(res_quant_f$logFC)>1 & -log10(res_quant_f$adj.P.Val)>4],pch=16,col="red",cex=0.6)
title("Volcano plot")
dev.off()

## P-value distribution
jpeg("Pdistr_quant_res_f.jpg")
hist(res_quant_f$P.Value,main="Filtered Data",xlab="P-values")
dev.off()


## Paired analysis
##############

# Make Design
designP_quant_f <- model.matrix(~individuals+cvst) #individuals as well
rownames(designP_quant_f) <- colnames(y_quant_f)

# Build Model
jpeg("limmatrendP_quant_f.jpg")
vP_quant_f <- voomWithQualityWeights(y_quant_f, designP_quant_f, plot=T,normalize.method = "quantile")
dev.off()
fitP_quant_f <- lmFit(vP_quant_f, designP_quant_f)
fitP_quant_f <- eBayes(fitP_quant_f)
resP_quant_f <- topTable(fitP_quant_f, coef="cvsttumor", n=nrow(vP_quant_f))
# coef="cvsttumor" => tumor vs. control 


## MAplot
jpeg("MAplotP_quant_res_f.jpg")
# limma::plotMA(fit_quant_f)
# MAplot: all data points
plot(rowSums(vP_quant_f$E)[rownames(resP_quant_f)],
     resP_quant_f$logFC,pch=16,cex=0.6)
# MA-plot: significant loci
points(rowSums(vP_quant_f$E)[rownames(resP_quant_f)[resP_quant_f$adj.P.Val<0.05]],
       resP_quant_f$logFC[resP_quant_f$adj.P.Val<0.05],pch=16,col="red",cex=0.6)
# X-axis
abline(0,0) 
dev.off()

## Volcano plot
jpeg("VolcanoP_quant_res_f.jpg")
plot(resP_quant_f$logFC,-log10(resP_quant_f$adj.P.Val),pch=16,cex=0.2,xlab="logFC",ylab="-log10(FDR)") 
points(resP_quant_f$logFC[abs(resP_quant_f$logFC)>1 & -log10(resP_quant_f$adj.P.Val)>4],
       -log10(resP_quant_f$adj.P.Val)[abs(resP_quant_f$logFC)>1 & -log10(resP_quant_f$adj.P.Val)>4],pch=16,col="red",cex=0.6)
title("Volcano plot")
dev.off()

## P-value distribution
jpeg("PdistrP_quant_res_f.jpg")
hist(resP_quant_f$P.Value,main="Filtered Data",xlab="P-values")
dev.off()



##################
##              ##
## Library Size ##
##              ##
##################

v_libSize_f <- voomWithQualityWeights(equalizeLibSizes(y_TMM_f$counts)$pseudo.counts)
countspm_libSize_f <- v_libSize_f$E

jpeg("density_libSize_f.jpg",width=1000,height=1000)
for (i in 1:ncol(counts_trim)){
  if (colSums(counts_trim)[i]<10000000){
    color <- "red"
  } else {
    color <- "black"
  }
  if (i==1){
    plot(density(countspm_libSize_f[,i]),col=color,ylim=c(0,max(density(countspm_libSize_f[,i])$y)*1.5),
         xlab="Expression (log)",main="Density plot")
  } else {
    lines(density(countspm_libSize_f[,i]),col=color)
  }
}
dev.off()

## MA-plots
jpeg("MAplot_libSize_f.jpg")
plot(rowMeans(countspm_libSize_f),rowMeans(countspm_libSize_f[,23:44])-rowMeans(countspm_libSize_f[,1:22]),pch=".",main="MA plot",xlab="A",ylab="M")
dev.off()


## Paired analysis
##############

# Make Design
designP_libSize_f <- model.matrix(~individuals+cvst) #individuals as well
rownames(designP_libSize_f) <- colnames(y_TMM_f)

# Build Model
jpeg("limmatrendP_libSize_f.jpg")
vP_libSize_f <- voomWithQualityWeights(equalizeLibSizes(y_TMM_f$counts)$pseudo.counts,plot=T)
dev.off()
fitP_libSize_f <- lmFit(vP_libSize_f, designP_libSize_f)
fitP_libSize_f <- eBayes(fitP_libSize_f)
resP_libSize_f <- topTable(fitP_libSize_f, coef="cvsttumor", n=nrow(vP_libSize_f))
# coef="cvsttumor" => tumor vs. control 

## MAplot
jpeg("MAplotP_libSize_res_f.jpg")
# limma::plotMA(fit_libSize_f)
# MAplot: all data points
plot(rowSums(vP_libSize_f$E)[rownames(resP_libSize_f)],
     resP_libSize_f$logFC,pch=16,cex=0.6)
# MA-plot: significant loci
points(rowSums(vP_libSize_f$E)[rownames(resP_libSize_f)[resP_libSize_f$adj.P.Val<0.05]],
       resP_libSize_f$logFC[resP_libSize_f$adj.P.Val<0.05],pch=16,col="red",cex=0.6)
# X-axis
abline(0,0) 
dev.off()

## Volcano plot
jpeg("VolcanoP_libSize_res_f.jpg")
plot(resP_libSize_f$logFC,-log10(resP_libSize_f$adj.P.Val),pch=16,cex=0.2,xlab="logFC",ylab="-log10(FDR)") 
points(resP_libSize_f$logFC[abs(resP_libSize_f$logFC)>1 & -log10(resP_libSize_f$adj.P.Val)>4],
       -log10(resP_libSize_f$adj.P.Val)[abs(resP_libSize_f$logFC)>1 & -log10(resP_libSize_f$adj.P.Val)>4],pch=16,col="red",cex=0.6)
title("Volcano plot")
dev.off()

## P-value distribution
jpeg("PdistrP_libSize_res_f.jpg")
hist(resP_libSize_f$P.Value,main="Unfiltered Data",xlab="P-values")
dev.off()

colscale <- color.scale(colSums(counts_TMM_f),extremes=c("orange","blue"))
jpeg("MDS_libSize_f.jpeg",width = 720,height = 720)
plotMDS(v_libSize_f$E,col=colscale,pch=c(1:19,35:38,1:19,35:38))
dev.off()
jpeg("MDS_libSize_f_cvst.jpeg")
plotMDS(v_libSize_f$E,labels=substr(cvst,1,1))
dev.off()


## Annotation
###################

head(genes)

genes <- genes[genes$ensembl_gene_id%in%rownames(counts_TMM_f),]
genes <- genes[!duplicated(genes$ensembl_gene_id),]

dim(genes)
dim(counts_trim)

annotate_results <- function(res,genes){
  genes <- genes[sort(genes$ensembl_gene_id,index.return=T)$ix,]
  
  res_sorted <-  res[sort(rownames(res),index.return=T)$ix,]
  res_sorted$geneSymbol <- ""
  
  if (sum(rownames(res_sorted[rownames(res_sorted)%in%genes$ensembl_gene_id,])==genes$ensembl_gene_id)!=length(res_sorted$geneSymbol[rownames(res_sorted)%in%genes$ensembl_gene_id])){
    stop("Mismatch in annotation information of results and gene annotation object, this will result in faulty annotation")
  }
  
  res_sorted$geneSymbol[rownames(res_sorted)%in%genes$ensembl_gene_id] <- genes$hgnc_symbol
  res_out <- res_sorted[sort(res_sorted$P.Value,index.return=T)$ix,]
  
  return(res_out)
}

# Unpaired analyses
res_TMM_f <- annotate_results(res_TMM_f,genes)
write.table(res_TMM_f,file="res_TMM_f.txt",sep="\t",col.names = T,row.names = T,quote = F)

res_TMM_nf <- annotate_results(res_TMM_nf,genes)

res_quant_f <- annotate_results(res_quant_f,genes)
write.table(res_quant_f,file="res_TMM_f.txt",sep="\t",col.names = T,row.names = T,quote = F)

res_quant_nf <- annotate_results(res_quant_nf,genes)

# Paired analyses
resP_TMM_f <- annotate_results(resP_TMM_f,genes)
write.table(resP_TMM_f,file="resP_TMM_f.txt",sep="\t",col.names = T,row.names = T,quote = F)

resP_TMM_nf <- annotate_results(resP_TMM_nf,genes)

resP_quant_f <- annotate_results(resP_quant_f,genes)
write.table(resP_quant_f,file="resP_TMM_f.txt",sep="\t",col.names = T,row.names = T,quote = F)

resP_quant_nf <- annotate_results(resP_quant_nf,genes)




####################
####################
####            ####
#### Exon-level ####
####            ####
####################
####################



###############
##           ##
## Norm: TMM ##
##           ##
###############


##################
## Non filtered ##
##################

## Make DGElist
counts_TMM_exon_nf <- counts_exon_trim
y_TMM_exon_nf <- DGEList(counts_TMM_exon_nf)
y_TMM_exon_nf <- calcNormFactors(y_TMM_exon_nf)

## Density-plots
jpeg("density_exon_nf.jpg",width=1000,height=1000)
for (i in 1:ncol(counts_exon_trim)){
  if (cvst[i]=="normal"){
    color="blue"
  } else {
    color="orange"
  }
  if (i==1){
    plot(density(log(y_TMM_exon_nf$counts[,i])),ylim=c(0,max(density(log(y_TMM_exon_nf$counts[,i]))$y)*2.1),xlim=c(-1,12),col=color)
  } else {
    lines(density(log(y_TMM_exon_nf$counts[,i])),col=color)
  }
}
dev.off()

## MA-plots
jpeg("MAplot_exon_nf.jpg")
plot(log2(rowMeans(y_TMM_exon_nf$counts)),log2(rowMeans(y_TMM_exon_nf$counts[,1:22])/rowMeans(y_TMM_exon_nf$counts[,23:44])),pch=".",main="MA plot",xlab="A",ylab="M")
dev.off()

## Density-plots normalized
countspm_exon_TMM_nf <- cpm(y_TMM_exon_nf,log=T)
jpeg("density_TMM_exon_nf.jpg",width=1000,height=1000)
for (i in 1:ncol(counts_exon_trim)){
  if (i==23){
    next
  }
  if (i==1){
    plot(density(countspm_exon_TMM_nf[,i]),ylim=c(0,max(density(log(counts_trim[,i]))$y)*25))
  } else {
    lines(density(countspm_exon_TMM_nf[,i]))
  }
}
dev.off()

## MA-plots
jpeg("MAplot_TMM_exon_nf.jpg")
plot(rowMeans(countspm_exon_TMM_nf),rowMeans(countspm_exon_TMM_nf[,23:44])-rowMeans(countspm_exon_TMM_nf[,1:22]),pch=".",main="MA plot",xlab="A",ylab="M")
dev.off()


## Unpaired analysis
##############

## Make Design
design_TMM_exon_nf <- model.matrix(~cvst)
rownames(design_TMM_exon_nf) <- colnames(y_TMM_exon_nf)

## Build Model
jpeg("limmatrend_TMM_exon_nf.jpg")
v_TMM_exon_nf <- voomWithQualityWeights(y_TMM_exon_nf, design_TMM_exon_nf,plot=T)
dev.off()
fit_TMM_exon_nf <- lmFit(v_TMM_exon_nf, design_TMM_exon_nf)
fit_TMM_exon_nf <- eBayes(fit_TMM_exon_nf)
res_TMM_exon_nf <- topTable(fit_TMM_exon_nf, coef="cvsttumor", n=nrow(v_TMM_exon_nf))
# coef="cvsttumor" => tumor vs. controle, see "design" object 

## MAplot
jpeg("MAplot_TMM_res_exon_nf.jpg")
# limma::plotMA(fit_TMM_exon_nf)
# MAplot: all data points
plot(rowSums(v_TMM_exon_nf$E)[rownames(res_TMM_exon_nf)],
     res_TMM_exon_nf$logFC,pch=16,cex=0.6)
# MA-plot: significant loci
points(rowSums(v_TMM_exon_nf$E)[rownames(res_TMM_exon_nf)[res_TMM_exon_nf$adj.P.Val<0.05]],
       res_TMM_exon_nf$logFC[res_TMM_exon_nf$adj.P.Val<0.05],pch=16,col="red",cex=0.6)
# X-axis
abline(0,0) 
dev.off()

## Volcano plot
jpeg("Volcano_TMM_res_exon_nf.jpg")
plot(res_TMM_exon_nf$logFC,-log10(res_TMM_exon_nf$adj.P.Val),pch=16,cex=0.2,xlab="logFC",ylab="-log10(FDR)") 
points(res_TMM_exon_nf$logFC[abs(res_TMM_exon_nf$logFC)>1 & -log10(res_TMM_exon_nf$adj.P.Val)>4],
       -log10(res_TMM_exon_nf$adj.P.Val)[abs(res_TMM_exon_nf$logFC)>1 & -log10(res_TMM_exon_nf$adj.P.Val)>4],pch=16,col="red",cex=0.6)
title("Volcano plot")
dev.off()

## P-value distribution
jpeg("Pdistr_TMM_res_exon_nf.jpg")
hist(res_TMM_exon_nf$P.Value,main="Unfiltered Data",xlab="P-values")
dev.off()


## Paired analysis
##############

## Make Design
designP_TMM_exon_nf <- model.matrix(~individuals+cvst) #individuals as well
rownames(designP_TMM_exon_nf) <- colnames(y_TMM_exon_nf)

## Build Model
jpeg("limmatrendP_TMM_exon_nf.jpg")
vP_TMM_exon_nf <- voomWithQualityWeights(y_TMM_exon_nf, designP_TMM_exon_nf, plot=T)
dev.off()
fitP_TMM_exon_nf <- lmFit(vP_TMM_exon_nf, designP_TMM_exon_nf)
fitP_TMM_exon_nf <- eBayes(fitP_TMM_exon_nf)
resP_TMM_exon_nf <- topTable(fitP_TMM_exon_nf, coef="cvsttumor", n=nrow(vP_TMM_exon_nf))
# coef="cvsttumor" => tumor vs. control 

## MAplot
jpeg("MAplotP_TMM_res_exon_nf.jpg")
# limma::plotMA(fit_TMM_exon_nf)
# MAplot: all data points
plot(rowSums(vP_TMM_exon_nf$E)[rownames(resP_TMM_exon_nf)],
     resP_TMM_exon_nf$logFC,pch=16,cex=0.6)
# MA-plot: significant loci
points(rowSums(vP_TMM_exon_nf$E)[rownames(resP_TMM_exon_nf)[resP_TMM_exon_nf$adj.P.Val<0.05]],
       resP_TMM_exon_nf$logFC[resP_TMM_exon_nf$adj.P.Val<0.05],pch=16,col="red",cex=0.6)
# X-axis
abline(0,0) 
dev.off()

## Volcano plot
jpeg("VolcanoP_TMM_res_exon_nf.jpg")
plot(resP_TMM_exon_nf$logFC,-log10(resP_TMM_exon_nf$adj.P.Val),pch=16,cex=0.2,xlab="logFC",ylab="-log10(FDR)") 
points(resP_TMM_exon_nf$logFC[abs(resP_TMM_exon_nf$logFC)>1 & -log10(resP_TMM_exon_nf$adj.P.Val)>4],
       -log10(resP_TMM_exon_nf$adj.P.Val)[abs(resP_TMM_exon_nf$logFC)>1 & -log10(resP_TMM_exon_nf$adj.P.Val)>4],pch=16,col="red",cex=0.6)
title("Volcano plot")
dev.off()

## P-value distribution
jpeg("PdistrP_TMM_res_exon_nf.jpg")
hist(resP_TMM_exon_nf$P.Value,main="Unfiltered Data",xlab="P-values")
dev.off()


##############
## Filtered ##
##############

## Make DGElist
y_TMM_exon_f <- DGEList(counts_TMM_exon_nf)
y_TMM_exon_f <- calcNormFactors(y_TMM_exon_f)

## Filtering
countspm_TMM_exon_f <- cpm(y_TMM_exon_f,log=T)
counts_TMM_exon_f <- counts_exon_trim[rowMeans(countspm_TMM_exon_f)>=cutoff_filter,]
y_TMM_exon_f <- DGEList(counts_TMM_exon_f)
y_TMM_exon_f <- calcNormFactors(y_TMM_exon_f)

## Density-plots
jpeg("density_exon_f.jpg",width=1000,height=1000)
for (i in 1:ncol(counts_trim)){
  if (cvst[i]=="normal"){
    color="blue"
  } else {
    color="orange"
  }
  if (i==1){
    plot(density(log(y_TMM_exon_f$counts[,i])),ylim=c(0,max(density(log(y_TMM_exon_f$counts[,i]))$y)*2.1),xlim=c(-1,12),col=color)
  } else {
    lines(density(log(y_TMM_exon_f$counts[,i])),col=color)
  }
}
dev.off()

## MA-plots
jpeg("MAplot_exon_f.jpg")
plot(log2(rowMeans(y_TMM_exon_f$counts)),log2(rowMeans(y_TMM_exon_f$counts[,1:22])/rowMeans(y_TMM_exon_f$counts[,23:44])),pch=".",main="MA plot",xlab="A",ylab="M")
dev.off()

## Density-plots normalized
countspm_exon_TMM_f <- cpm(y_TMM_exon_f,log=T)
jpeg("density_TMM_exon_f.jpg",width=1000,height=1000)
for (i in 1:ncol(counts_trim)){
  if (i==23){
    next
  }
  if (i==1){
    plot(density(countspm_exon_TMM_f[,i]),ylim=c(0,max(density(log(counts_trim[,i]))$y)*5))
  } else {
    lines(density(countspm_exon_TMM_f[,i]))
  }
}
dev.off()

## MA-plots
jpeg("MAplot_TMM_exon_f.jpg")
plot(rowMeans(countspm_exon_TMM_f),rowMeans(countspm_exon_TMM_f[,23:44])-rowMeans(countspm_exon_TMM_f[,1:22]),pch=".",main="MA plot",xlab="A",ylab="M")
dev.off()


## Unpaired analysis
##############

## Make Design
design_TMM_exon_f <- model.matrix(~cvst)
rownames(design_TMM_exon_f) <- colnames(y_TMM_exon_f)

## Build Model
jpeg("limmatrend_TMM_exon_f.jpg")
v_TMM_exon_f <- voomWithQualityWeights(y_TMM_exon_f, design_TMM_exon_f,plot=T)
dev.off()
fit_TMM_exon_f <- lmFit(v_TMM_exon_f, design_TMM_exon_f)
fit_TMM_exon_f <- eBayes(fit_TMM_exon_f)
res_TMM_exon_f <- topTable(fit_TMM_exon_f, coef="cvsttumor", n=nrow(v_TMM_exon_f))
# coef="cvsttumor" => tumor vs. controle, see "design" object 

## MAplot
jpeg("MAplot_TMM_res_exon_f.jpg")
# limma::plotMA(fit_TMM_nf)
# MAplot: all data points
plot(rowSums(v_TMM_exon_f$E)[rownames(res_TMM_exon_f)],
     res_TMM_exon_f$logFC,pch=16,cex=0.6)
# MA-plot: significant loci
points(rowSums(v_TMM_exon_f$E)[rownames(res_TMM_exon_f)[res_TMM_exon_f$adj.P.Val<0.05]],
       res_TMM_exon_f$logFC[res_TMM_exon_f$adj.P.Val<0.05],pch=16,col="red",cex=0.6)
# X-axis
abline(0,0) 
dev.off()

## Volcano plot
jpeg("Volcano_TMM_res_exon_f.jpg")
plot(res_TMM_exon_f$logFC,-log10(res_TMM_exon_f$adj.P.Val),pch=16,cex=0.2,xlab="logFC",ylab="-log10(FDR)") 
points(res_TMM_exon_f$logFC[abs(res_TMM_exon_f$logFC)>1 & -log10(res_TMM_exon_f$adj.P.Val)>4],
       -log10(res_TMM_exon_f$adj.P.Val)[abs(res_TMM_exon_f$logFC)>1 & -log10(res_TMM_exon_f$adj.P.Val)>4],pch=16,col="red",cex=0.6)
title("Volcano plot")
dev.off()

## P-value distribution
jpeg("Pdistr_TMM_res_exon_f.jpg")
hist(res_TMM_exon_f$P.Value,main="Filtered Data",xlab="P-values")
dev.off()


## Paired analysis
##############

# Make Design
designP_TMM_exon_f <- model.matrix(~individuals+cvst) #individuals as well
rownames(designP_TMM_exon_f) <- colnames(y_TMM_exon_f)

# Build Model
jpeg("limmatrendP_TMM_exon_f.jpg")
vP_TMM_exon_f <- voomWithQualityWeights(y_TMM_exon_f, designP_TMM_exon_f, plot=T)
dev.off()
fitP_TMM_exon_f <- lmFit(vP_TMM_exon_f, designP_TMM_exon_f)
fitP_TMM_exon_f <- eBayes(fitP_TMM_exon_f)
resP_TMM_exon_f <- topTable(fitP_TMM_exon_f, coef="cvsttumor", n=nrow(vP_TMM_exon_f))
# coef="cvsttumor" => tumor vs. control 


## MAplot
jpeg("MAplotP_TMM_res_exon_f.jpg")
# limma::plotMA(fit_TMM_nf)
# MAplot: all data points
plot(rowSums(vP_TMM_exon_f$E)[rownames(resP_TMM_exon_f)],
     resP_TMM_exon_f$logFC,pch=16,cex=0.6)
# MA-plot: significant loci
points(rowSums(vP_TMM_exon_f$E)[rownames(resP_TMM_exon_f)[resP_TMM_exon_f$adj.P.Val<0.05]],
       resP_TMM_exon_f$logFC[resP_TMM_exon_f$adj.P.Val<0.05],pch=16,col="red",cex=0.6)
# X-axis
abline(0,0) 
dev.off()

## Volcano plot
jpeg("VolcanoP_TMM_res_exon_f.jpg")
plot(resP_TMM_exon_f$logFC,-log10(resP_TMM_exon_f$adj.P.Val),pch=16,cex=0.2,xlab="logFC",ylab="-log10(FDR)") 
points(resP_TMM_exon_f$logFC[abs(resP_TMM_exon_f$logFC)>1 & -log10(resP_TMM_exon_f$adj.P.Val)>4],
       -log10(resP_TMM_exon_f$adj.P.Val)[abs(resP_TMM_exon_f$logFC)>1 & -log10(resP_TMM_exon_f$adj.P.Val)>4],pch=16,col="red",cex=0.6)
title("Volcano plot")
dev.off()

## P-value distribution
jpeg("PdistrP_TMM_res_exon_f.jpg")
hist(resP_TMM_exon_f$P.Value,main="Filtered Data",xlab="P-values")
dev.off()



####################
##                ##
## Norm: Quantile ##
##                ##
####################


##################
## Non filtered ##
##################

## Make DGElist
counts_quant_exon_nf <- counts_exon_trim
y_quant_exon_nf <- DGEList(counts_quant_exon_nf)
y_quant_exon_nf <- calcNormFactors(y_quant_exon_nf)

## Make Design
design_quant_exon_nf <- model.matrix(~cvst)
rownames(design_quant_exon_nf) <- colnames(y_quant_exon_nf)

## Density-plots normalized: quantile
countspm_quant_exon_nf <- voomWithQualityWeights(y_quant_exon_nf, design_quant_exon_nf,normalize.method = "quantile")$E
jpeg("density_quant_exon_nf.jpg",width=1000,height=1000)
for (i in 1:ncol(counts_trim)){
  if (colSums(counts_trim)[i]<10000000){
    color <- "red"
  } else {
    color <- "black"
  }
  if (i==1){
    plot(density(countspm_quant_exon_nf[,i]),col=color)
  } else {
    lines(density(countspm_quant_exon_nf[,i]),col=color)
  }
}
dev.off()

## MA-plots
jpeg("MAplot_quant_exon_nf.jpg")
plot(rowMeans(countspm_quant_exon_nf),rowMeans(countspm_quant_exon_nf[,23:44])-rowMeans(countspm_quant_exon_nf[,1:22]),pch=".",main="MA plot",xlab="A",ylab="M")
dev.off()


## Unpaired analysis
##############

## Build Model
jpeg("limmatrend_quant_exon_nf.jpg")
v_quant_exon_nf <- voomWithQualityWeights(y_quant_exon_nf, design_quant_exon_nf,plot=T,normalize.method = "quantile")
dev.off()
fit_quant_exon_nf <- lmFit(v_quant_exon_nf, design_quant_exon_nf)
fit_quant_exon_nf <- eBayes(fit_quant_exon_nf)
res_quant_exon_nf <- topTable(fit_quant_exon_nf, coef="cvsttumor", n=nrow(v_quant_exon_nf))
# coef="cvsttumor" => tumor vs. controle, see "design" object 

## MAplot
jpeg("MAplot_quant_res_exon_nf.jpg")
# limma::plotMA(fit_quant_exon_nf)
# MAplot: all data points
plot(rowSums(v_quant_exon_nf$E)[rownames(res_quant_exon_nf)],
     res_quant_exon_nf$logFC,pch=16,cex=0.6)
# MA-plot: significant loci
points(rowSums(v_quant_exon_nf$E)[rownames(res_quant_exon_nf)[res_quant_exon_nf$adj.P.Val<0.05]],
       res_quant_exon_nf$logFC[res_quant_exon_nf$adj.P.Val<0.05],pch=16,col="red",cex=0.6)
# X-axis
abline(0,0) 
dev.off()

## Volcano plot
jpeg("Volcano_quant_res_exon_nf.jpg")
plot(res_quant_exon_nf$logFC,-log10(res_quant_exon_nf$adj.P.Val),pch=16,cex=0.2,xlab="logFC",ylab="-log10(FDR)") 
points(res_quant_exon_nf$logFC[abs(res_quant_exon_nf$logFC)>1 & -log10(res_quant_exon_nf$adj.P.Val)>4],
       -log10(res_quant_exon_nf$adj.P.Val)[abs(res_quant_exon_nf$logFC)>1 & -log10(res_quant_exon_nf$adj.P.Val)>4],pch=16,col="red",cex=0.6)
title("Volcano plot")
dev.off()

## P-value distribution
jpeg("Pdistr_quant_res_exon_nf.jpg")
hist(res_quant_exon_nf$P.Value,main="Unfiltered Data",xlab="P-values")
dev.off()


## Paired analysis
##############

## Make Design
designP_quant_exon_nf <- model.matrix(~individuals+cvst) #individuals as well
rownames(designP_quant_exon_nf) <- colnames(y_quant_exon_nf)

## Build Model
jpeg("limmatrendP_quant_exon_nf.jpg")
vP_quant_exon_nf <- voomWithQualityWeights(y_quant_exon_nf, designP_quant_exon_nf, plot=T,normalize.method = "quantile")
dev.off()
fitP_quant_exon_nf <- lmFit(vP_quant_exon_nf, designP_quant_exon_nf)
fitP_quant_exon_nf <- eBayes(fitP_quant_exon_nf)
resP_quant_exon_nf <- topTable(fitP_quant_exon_nf, coef="cvsttumor", n=nrow(vP_quant_exon_nf))
# coef="cvsttumor" => tumor vs. control 

## MAplot
jpeg("MAplotP_quant_res_exon_nf.jpg")
# limma::plotMA(fit_quant_exon_nf)
# MAplot: all data points
plot(rowSums(vP_quant_exon_nf$E)[rownames(resP_quant_exon_nf)],
     resP_quant_exon_nf$logFC,pch=16,cex=0.6)
# MA-plot: significant loci
points(rowSums(vP_quant_exon_nf$E)[rownames(resP_quant_exon_nf)[resP_quant_exon_nf$adj.P.Val<0.05]],
       resP_quant_exon_nf$logFC[resP_quant_exon_nf$adj.P.Val<0.05],pch=16,col="red",cex=0.6)
# X-axis
abline(0,0) 
dev.off()

## Volcano plot
jpeg("VolcanoP_quant_res_exon_nf.jpg")
plot(resP_quant_exon_nf$logFC,-log10(resP_quant_exon_nf$adj.P.Val),pch=16,cex=0.2,xlab="logFC",ylab="-log10(FDR)") 
points(resP_quant_exon_nf$logFC[abs(resP_quant_exon_nf$logFC)>1 & -log10(resP_quant_exon_nf$adj.P.Val)>4],
       -log10(resP_quant_exon_nf$adj.P.Val)[abs(resP_quant_exon_nf$logFC)>1 & -log10(resP_quant_exon_nf$adj.P.Val)>4],pch=16,col="red",cex=0.6)
title("Volcano plot")
dev.off()

## P-value distribution
jpeg("PdistrP_quant_res_exon_nf.jpg")
hist(resP_quant_exon_nf$P.Value,main="Unfiltered Data",xlab="P-values")
dev.off()


##############
## Filtered ##
##############

## Make DGElist
counts_quant_exon_nf <- counts_exon_trim
y_quant_exon_f <- DGEList(counts_quant_exon_nf)
y_quant_exon_f <- calcNormFactors(y_quant_exon_f)

## Filtering
countspm_quant_exon_f <- cpm(y_quant_exon_f,log=T)
counts_quant_exon_f <- counts_exon_trim[rowMeans(countspm_quant_exon_f)>=cutoff_filter,]
y_quant_exon_f <- DGEList(counts_quant_exon_f)
y_quant_exon_f <- calcNormFactors(y_quant_exon_f)

## Make Design
design_quant_exon_f <- model.matrix(~cvst)
rownames(design_quant_exon_f) <- colnames(y_quant_exon_f)

## Density-plots normalized: quantile
countspm_quant_exon_f <- voomWithQualityWeights(y_quant_exon_f, design_quant_exon_f,normalize.method = "quantile")$E
jpeg("density_quant_exon_f.jpg",width=1000,height=1000)
for (i in 1:ncol(counts_trim)){
  if (colSums(counts_trim)[i]<10000000){
    color <- "red"
  } else {
    color <- "black"
  }
  if (i==1){
    plot(density(countspm_quant_exon_f[,i]),col=color)
  } else {
    lines(density(countspm_quant_exon_f[,i]),col=color)
  }
}
dev.off()

## MA-plots
jpeg("MAplot_quant_exon_f.jpg")
plot(rowMeans(countspm_quant_exon_f),rowMeans(countspm_quant_exon_f[,23:44])-rowMeans(countspm_quant_exon_f[,1:22]),pch=".",main="MA plot",xlab="A",ylab="M")
dev.off()


## Unpaired analysis
##############

## Build Model
jpeg("limmatrend_quant_exon_f.jpg")
v_quant_exon_f <- voomWithQualityWeights(y_quant_exon_f, design_quant_exon_f,plot=T,normalize.method = "quantile")
dev.off()
fit_quant_exon_f <- lmFit(v_quant_exon_f, design_quant_exon_f)
fit_quant_exon_f <- eBayes(fit_quant_exon_f)
res_quant_exon_f <- topTable(fit_quant_exon_f, coef="cvsttumor", n=nrow(v_quant_exon_f))
# coef="cvsttumor" => tumor vs. controle, see "design" object 

## MAplot
jpeg("MAplot_quant_res_exon_f.jpg")
# limma::plotMA(fit_quant_exon_f)
# MAplot: all data points
plot(rowSums(v_quant_exon_f$E)[rownames(res_quant_exon_f)],
     res_quant_exon_f$logFC,pch=16,cex=0.6)
# MA-plot: significant loci
points(rowSums(v_quant_exon_f$E)[rownames(res_quant_exon_f)[res_quant_exon_f$adj.P.Val<0.05]],
       res_quant_exon_f$logFC[res_quant_exon_f$adj.P.Val<0.05],pch=16,col="red",cex=0.6)
# X-axis
abline(0,0) 
dev.off()

## Volcano plot
jpeg("Volcano_quant_res_exon_f.jpg")
plot(res_quant_exon_f$logFC,-log10(res_quant_exon_f$adj.P.Val),pch=16,cex=0.2,xlab="logFC",ylab="-log10(FDR)") 
points(res_quant_exon_f$logFC[abs(res_quant_exon_f$logFC)>1 & -log10(res_quant_exon_f$adj.P.Val)>4],
       -log10(res_quant_exon_f$adj.P.Val)[abs(res_quant_exon_f$logFC)>1 & -log10(res_quant_exon_f$adj.P.Val)>4],pch=16,col="red",cex=0.6)
title("Volcano plot")
dev.off()

## P-value distribution
jpeg("Pdistr_quant_res_exon_f.jpg")
hist(res_quant_exon_f$P.Value,main="Filtered Data",xlab="P-values")
dev.off()


## Paired analysis
##############

# Make Design
designP_quant_exon_f <- model.matrix(~individuals+cvst) #individuals as well
rownames(designP_quant_exon_f) <- colnames(y_quant_exon_f)

# Build Model
jpeg("limmatrendP_quant_exon_nf.jpg")
vP_quant_exon_f <- voomWithQualityWeights(y_quant_exon_f, designP_quant_exon_f, plot=T,normalize.method = "quantile")
dev.off()
fitP_quant_exon_f <- lmFit(vP_quant_exon_f, designP_quant_exon_f)
fitP_quant_exon_f <- eBayes(fitP_quant_exon_f)
resP_quant_exon_f <- topTable(fitP_quant_exon_f, coef="cvsttumor", n=nrow(vP_quant_exon_f))
# coef="cvsttumor" => tumor vs. control 

## MAplot
jpeg("MAplotP_quant_res_exon_f.jpg")
# limma::plotMA(fit_quant_exon_f)
# MAplot: all data points
plot(rowSums(vP_quant_exon_f$E)[rownames(resP_quant_exon_f)],
     resP_quant_exon_f$logFC,pch=16,cex=0.6)
# MA-plot: significant loci
points(rowSums(vP_quant_exon_f$E)[rownames(resP_quant_exon_f)[resP_quant_exon_f$adj.P.Val<0.05]],
       resP_quant_exon_f$logFC[resP_quant_exon_f$adj.P.Val<0.05],pch=16,col="red",cex=0.6)
# X-axis
abline(0,0) 
dev.off()

## Volcano plot
jpeg("VolcanoP_quant_res_exon_f.jpg")
plot(resP_quant_exon_f$logFC,-log10(resP_quant_exon_f$adj.P.Val),pch=16,cex=0.2,xlab="logFC",ylab="-log10(FDR)") 
points(resP_quant_exon_f$logFC[abs(resP_quant_exon_f$logFC)>1 & -log10(resP_quant_exon_f$adj.P.Val)>4],
       -log10(resP_quant_exon_f$adj.P.Val)[abs(resP_quant_exon_f$logFC)>1 & -log10(resP_quant_exon_f$adj.P.Val)>4],pch=16,col="red",cex=0.6)
title("Volcano plot")
dev.off()

## P-value distribution
jpeg("PdistrP_quant_res_exon_f.jpg")
hist(resP_quant_exon_f$P.Value,main="Unfiltered Data",xlab="P-values")
dev.off()


#########################
# Analysis without voom #
#########################

# logCPM_TMM_f_nov <- cpm(y_TMM_f,log=T)
y_quant_f_nov <- vP_quant_f$E
y_quant_f_nov <- log2(normalizeBetweenArrays(y_quant_f$counts,method="quantile")+0.5)
fitP_quant_f_nov <- lmFit(y_quant_f_nov, designP_TMM_f)
fitP_quant_f_nov <- eBayes(fitP_quant_f_nov)
resP_quant_f_nov <- topTable(fitP_quant_f_nov, coef="cvsttumor", n=nrow(y_quant_f_nov))

## MAplot
jpeg("MAplotP_quant_res_f_nov.jpg")
# limma::plotMA(fit_quant_f)
# MAplot: all data points
plot(rowSums(y_quant_f_nov)[rownames(resP_quant_f_nov)],
     resP_quant_f_nov$logFC,pch=16,cex=0.6)
# MA-plot: significant loci
points(rowSums(y_quant_f_nov)[rownames(resP_quant_f_nov)[resP_quant_f_nov$adj.P.Val<0.05]],
       resP_quant_f_nov$logFC[resP_quant_f_nov$adj.P.Val<0.05],pch=16,col="red",cex=0.6)
# X-axis
abline(0,0) 
dev.off()

## Volcano plot
jpeg("VolcanoP_quant_res_f_nov.jpg")
plot(resP_quant_f_nov$logFC,-log10(resP_quant_f_nov$adj.P.Val),pch=16,cex=0.2,xlab="logFC",ylab="-log10(FDR)") 
points(resP_quant_f_nov$logFC[abs(resP_quant_f_nov$logFC)>1 & -log10(resP_quant_f_nov$adj.P.Val)>4],
       -log10(resP_quant_f_nov$adj.P.Val)[abs(resP_quant_f_nov$logFC)>1 & -log10(resP_quant_f_nov$adj.P.Val)>4],pch=16,col="red",cex=0.6)
title("Volcano plot")
dev.off()

## P-value distribution
jpeg("PdistrP_quant_res_f_nov.jpg")
hist(resP_quant_f_nov$P.Value,main="Filtered Data",xlab="P-values")
dev.off()

## Annotation
resP_quant_f_nov_sorted <- resP_quant_f_nov[sort(rownames(resP_quant_f_nov),index.return=T)$ix,]
resP_quant_f_nov_sorted$geneSymbol <- ""
resP_quant_f_nov_sorted$geneSymbol[rownames(resP_quant_f_nov_sorted)%in%genes$ensembl_gene_id] <- genes$hgnc_symbol[genes$ensembl_gene_id%in%rownames(resP_quant_f_nov_sorted[rownames(resP_quant_f_nov_sorted)%in%genes$ensembl_gene_id,])]
resP_quant_f_nov <- resP_quant_f_nov_sorted[sort(resP_quant_f_nov_sorted$P.Value,index.return=T)$ix,]

# Comparison Tophits
list_quant <- rownames(head(resP_quant_f))
list_quant_nov <- rownames(head(resP_quant_f_nov))
jpeg("comp_quant_quantnov_tophits.jpg",width=600,height=800)
par(mfrow=c(4,3))
for (i in 1:length(list_quant)){
  summary(as.numeric(counts_quant_f[list_quant[i],]))
  
  plot(density(as.numeric(vP_quant_f$E[list_quant[i],])[1:22]-as.numeric(vP_quant_f$E[list_quant[i],])[23:44]),
       col="blue",main="Density",xlab="logFC")
  legend("topright",legend=c(resP_quant_f$geneSymbol[i],
                             paste("mean",summary(as.numeric(vP_quant_f$E[list_quant[i],])[1:22]-as.numeric(vP_quant_f$E[list_quant[i],])[23:44])[4]),
                             paste("3rd Q",summary(as.numeric(vP_quant_f$E[list_quant[i],])[1:22]-as.numeric(vP_quant_f$E[list_quant[i],])[23:44])[5]),
                             paste("max",summary(as.numeric(vP_quant_f$E[list_quant[i],])[1:22]-as.numeric(vP_quant_f$E[list_quant[i],])[23:44])[6])))
  
}
for (i in 1:length(list_quant_nov)){
  summary(as.numeric(counts_quant_f[list_quant_nov[i],]))
  
  plot(density(as.numeric(y_quant_f_nov[list_quant_nov[i],])[1:22]-as.numeric(y_quant_f_nov[list_quant_nov[i],])[23:44]),
       col="blue",main="Density")
  legend("topright",legend=c(resP_quant_f_nov$geneSymbol[i],
                             paste("mean",summary(as.numeric(y_quant_f_nov[list_quant_nov[i],])[1:22]-as.numeric(y_quant_f_nov[list_quant_nov[i],])[23:44])[4]),
                             paste("3rd Q",summary(as.numeric(y_quant_f_nov[list_quant_nov[i],])[1:22]-as.numeric(y_quant_f_nov[list_quant_nov[i],])[23:44])[5]),
                             paste("max",summary(as.numeric(y_quant_f_nov[list_quant_nov[i],])[1:22]-as.numeric(y_quant_f_nov[list_quant_nov[i],])[23:44])[6])))
  
}
dev.off()

## Analysis without sample weighting
###################################

vP_quant_f_nosw <- voom(y_quant_f, designP_quant_f, plot=T,normalize.method = "quantile")
y_quant_f_nosw <- vP_quant_f_nosw$E
fitP_quant_f_nosw <- lmFit(vP_quant_f_nosw, designP_TMM_f)
fitP_quant_f_nosw <- eBayes(fitP_quant_f_nosw)
resP_quant_f_nosw <- topTable(fitP_quant_f_nosw, coef="cvsttumor", n=nrow(y_quant_f))

## MAplot
jpeg("MAplotP_quant_res_f_nosw.jpg")
# limma::plotMA(fit_quant_f)
# MAplot: all data points
plot(rowSums(y_quant_f_nosw)[rownames(resP_quant_f_nosw)],
     resP_quant_f_nosw$logFC,pch=16,cex=0.6)
# MA-plot: significant loci
points(rowSums(y_quant_f_nosw)[rownames(resP_quant_f_nosw)[resP_quant_f_nosw$adj.P.Val<0.05]],
       resP_quant_f_nosw$logFC[resP_quant_f_nosw$adj.P.Val<0.05],pch=16,col="red",cex=0.6)
# X-axis
abline(0,0) 
dev.off()

## Volcano plot
jpeg("VolcanoP_quant_res_f_nosw.jpg")
plot(resP_quant_f_nosw$logFC,-log10(resP_quant_f_nosw$adj.P.Val),pch=16,cex=0.2,xlab="logFC",ylab="-log10(FDR)") 
points(resP_quant_f_nosw$logFC[abs(resP_quant_f_nosw$logFC)>1 & -log10(resP_quant_f_nosw$adj.P.Val)>4],
       -log10(resP_quant_f_nosw$adj.P.Val)[abs(resP_quant_f_nosw$logFC)>1 & -log10(resP_quant_f_nosw$adj.P.Val)>4],pch=16,col="red",cex=0.6)
title("Volcano plot")
dev.off()

## P-value distribution
jpeg("PdistrP_quant_res_f_nosw.jpg")
hist(resP_quant_f_nosw$P.Value,main="Filtered Data",xlab="P-values")
dev.off()

###########
# Numbers #
###########

dim(resP_quant_nf)
sum(resP_quant_nf$adj.P.Val<0.05)
sum(resP_quant_nf$adj.P.Val<0.05)/nrow(resP_quant_nf)

dim(resP_quant_exon_nf)
sum(resP_quant_exon_nf$adj.P.Val<0.05)
sum(resP_quant_exon_nf$adj.P.Val<0.05)/nrow(resP_quant_exon_nf)


dim(resP_quant_f)
sum(resP_quant_f$adj.P.Val<0.05)
sum(resP_quant_f$adj.P.Val<0.05)/nrow(resP_quant_f)

dim(resP_quant_exon_f)
sum(resP_quant_exon_f$adj.P.Val<0.05)
sum(resP_quant_exon_f$adj.P.Val<0.05)/nrow(resP_quant_exon_f)

dim(res_TMM_nf)
sum(res_TMM_nf$adj.P.Val<0.05)
dim(res_TMM_f)
sum(res_TMM_f$adj.P.Val<0.05)

dim(res_TMM_exon_nf)
sum(res_TMM_exon_nf$adj.P.Val<0.05)
dim(res_TMM_exon_f)
sum(res_TMM_exon_f$adj.P.Val<0.05)


#################
# Plots Chapter #
#################

letters <- as.character(c("a":"z"))
i <- 1
jpeg("ExonvsGene.jpeg",width=1200,height=1200)
par(mfrow=c(2,2),mar=c(6,6,6,6))
hist(resP_quant_exon_nf$P.Value,cex.lab=2.5, cex.axis=2.5, cex.main=2.5, cex.sub=2.5,main="Exon level",xlab="P-values")
mtext(toupper(letters[i]), side = 3, adj = 0.05, line = 2.4, cex=3, font=4)
i <- i + 1
hist(resP_quant_nf$P.Value,cex.lab=2.5, cex.axis=2.5, cex.main=2.5, cex.sub=2.5,main="Gene level",xlab="P-values")
mtext(toupper(letters[i]), side = 3, adj = 0.05, line = 2.4,cex=3,font=4)
i <- i + 1
hist(resP_quant_exon_f$P.Value,cex.lab=2.5, cex.axis=2.5, cex.main=2.5, cex.sub=2.5,main="Exon level (filtered)",xlab="P-values")
mtext(toupper(letters[i]), side = 3, adj = 0.05, line = 2.4,cex=3,font=4)
i <- i + 1
hist(resP_quant_f$P.Value,cex.lab=2.5, cex.axis=2.5, cex.main=2.5, cex.sub=2.5,main="Gene level (filtered)",xlab="P-values")
mtext(toupper(letters[i]), side = 3, adj = 0.05, line = 2.4,cex=3,font=4)
i <- i + 1
dev.off()

i <- 1
jpeg("Normalisation.jpg",width=2400,height=3200)
par(mfrow=c(4,3),mar=c(6,6,6,6))
colscale <- color.scale(colSums(counts_TMM_f),extremes=c("orange","blue")) 
# Raw
for (j in 1:ncol(counts_trim)){
  if (cvst[j]=="normal"){
    color="blue"
  } else {
    color="orange"
  }
  if (j==1){
    plot(density(log(y_TMM_f$counts[,j])),ylim=c(0,max(density(log(y_TMM_f$counts[,i]))$y)*2.1),xlim=c(-1,12),col=color,cex.lab=2.5, cex.axis=2.5, cex.main=2.5, cex.sub=2.5,main="Density plot")
  } else {
    lines(density(log(y_TMM_f$counts[,j])),col=color)
  }
}
mtext(toupper(letters[i]), side = 3, adj = 0.05, line = 2.4,cex=3,font=4)
i <- i + 1
plot(log2(rowMeans(y_TMM_f$counts)),log2(rowMeans(y_TMM_f$counts[,1:22])/rowMeans(y_TMM_f$counts[,23:44])),pch=".",main="MA plot",xlab="A",ylab="M",cex.lab=2.5, cex.axis=2.5, cex.main=2.5, cex.sub=2.5, cex=3.5)
abline(h=0,col="red")
mtext(toupper(letters[i]), side = 3, adj = 0.05, line = 2.4,cex=3,font=4)
i <- i + 1
plotMDS(y_TMM_f$counts,col=colscale,pch=c(1:19,35:38,1:19,35:38),main="MDS plot",cex.lab=2.5, cex.axis=2.5, cex.main=2.5, cex.sub=2.5, cex=3.5)
mtext(toupper(letters[i]), side = 3, adj = 0.05, line = 2.4,cex=3,font=4)
i <- i + 1

# LibSize
for (j in 1:ncol(counts_trim)){
  if (cvst[j]=="normal"){
    color="blue"
  } else {
    color="orange"
  }
  if (j==1){
    plot(density(countspm_libSize_f[,j]),col=color,ylim=c(0,max(density(countspm_libSize_f[,i])$y)*1.5),cex.lab=2.5, cex.axis=2.5, cex.main=2.5, cex.sub=2.5, cex=3.5,main="Density plot")
  } else {
    lines(density(countspm_libSize_f[,j]),col=color)
  }
}
mtext(toupper(letters[i]), side = 3, adj = 0.05, line = 2.4,cex=3,font=4)
i <- i + 1
plot(rowMeans(countspm_libSize_f),rowMeans(countspm_libSize_f[,23:44])-rowMeans(countspm_libSize_f[,1:22]),pch=".",main="MA plot",xlab="A",ylab="M",cex.lab=2.5, cex.axis=2.5, cex.main=2.5, cex.sub=2.5, cex=3.5)
abline(h=0,col="red")
mtext(toupper(letters[i]), side = 3, adj = 0.05, line = 2.4,cex=3,font=4)
i <- i + 1
plotMDS(v_libSize_f$E,col=colscale,pch=c(1:19,35:38,1:19,35:38),main="MDS plot",cex.lab=2.5, cex.axis=2.5, cex.main=2.5, cex.sub=2.5, cex=3.5)
mtext(toupper(letters[i]), side = 3, adj = 0.05, line = 2.4,cex=3,font=4)
i <- i + 1

# TMM
countspm_TMM_f <- cpm(y_TMM_f,log=T)
for (j in 1:ncol(counts_trim)){
  if (cvst[j]=="normal"){
    color="blue"
  } else {
    color="orange"
  }
  if (j==1){
    plot(density(countspm_TMM_f[,j]),col=color,ylim=c(0,max(density(log(counts_trim[,i]))$y)*2),cex.lab=2.5, cex.axis=2.5, cex.main=2.5, cex.sub=2.5, cex=3.5,main="Density plot")
  } else {
    lines(density(countspm_TMM_f[,j]),col=color)
  }
}
mtext(toupper(letters[i]), side = 3, adj = 0.05, line = 2.4,cex=3,font=4)
i <- i + 1
plot(rowMeans(countspm_TMM_f),rowMeans(countspm_TMM_f[,23:44])-rowMeans(countspm_TMM_f[,1:22]),pch=".",main="MA plot",xlab="A",ylab="M",cex.lab=2.5, cex.axis=2.5, cex.main=2.5, cex.sub=2.5, cex=3.5)
abline(h=0,col="red")
mtext(toupper(letters[i]), side = 3, adj = 0.05, line = 2.4,cex=3,font=4)
i <- i + 1
plotMDS(cpm(y_TMM_f,log=T),col=colscale,pch=c(1:19,35:38,1:19,35:38),main="MDS plot",cex.lab=2.5, cex.axis=2.5, cex.main=2.5, cex.sub=2.5, cex=3.5)
mtext(toupper(letters[i]), side = 3, adj = 0.05, line = 2.4,cex=3,font=4)
i <- i + 1

# Quantile
countspm_quant_f <- voomWithQualityWeights(y_quant_f, design_quant_f,normalize.method = "quantile")$E
for (j in 1:ncol(counts_trim)){
  if (cvst[j]=="normal"){
    color="blue"
  } else {
    color="orange"
  }
  if (j==1){
    plot(density(countspm_quant_f[,j]),col=color,ylim=c(0,max(density(countspm_quant_f[,i])$y)*2),cex.lab=2.5, cex.axis=2.5, cex.main=2.5, cex.sub=2.5, cex=3.5,main="Density plot")
  } else {
    lines(density(countspm_quant_f[,j]),col=color)
  }
}
mtext(toupper(letters[i]), side = 3, adj = 0.05, line = 2.4,cex=3,font=4)
i <- i + 1
plot(rowMeans(countspm_quant_f),rowMeans(countspm_quant_f[,23:44])-rowMeans(countspm_quant_f[,1:22]),pch=".",main="MA plot",xlab="A",ylab="M",cex.lab=2.5, cex.axis=2.5, cex.main=2.5, cex.sub=2.5, cex=3.5)
abline(h=0,col="red")
mtext(toupper(letters[i]), side = 3, adj = 0.05, line = 2.4,cex=3,font=4)
i <- i + 1
plotMDS(v_quant_f$E,col=colscale,pch=c(1:19,35:38,1:19,35:38),main="MDS plot",cex.lab=2.5, cex.axis=2.5, cex.main=2.5, cex.sub=2.5, cex=3.5)
mtext(toupper(letters[i]), side = 3, adj = 0.05, line = 2.4,cex=3,font=4)
i <- i + 1
dev.off()

a<-vennCounts(cbind(resP_quant_f_nov$adj.P.Val[sort(rownames(resP_quant_f_nov),index.return=T)$ix]<0.05,
                    resP_quant_f_nosw$adj.P.Val[sort(rownames(resP_quant_f_nosw),index.return=T)$ix]<0.05,
                    resP_quant_f$adj.P.Val[sort(rownames(resP_quant_f),index.return=T)$ix]<0.05))
colnames(a) <- c("no VT, no SW","VT, no SW","VT + SW","Counts")
vP_quant_f <- voomWithQualityWeights(y_quant_f, designP_quant_f, plot=T,normalize.method = "quantile",save.plot=T)
sx <- vP_quant_f$voom.xy$x
sy <- vP_quant_f$voom.xy$y
l <- vP_quant_f$voom.line
aw <- vP_quant_f$targets$sample.weights

i <- 1
jpeg("VoomvsNoVoom.jpg",width=2400,height=2400)
par(mfrow=c(2,2),mar=c(6,6,6,6))
vennDiagram(a,cex.lab=2.5, cex.axis=2.5, cex.main=2.5, cex.sub=2.5, cex=3.5)
title("Venn Diagram",cex.main=2.5)
mtext(toupper(letters[i]), side = 3, adj = 0.05, line = 2.4,cex=3,font=4)
i <- i + 1
plot(sx, sy, xlab = "log2( count size + 0.5 )", ylab = "Sqrt( standard deviation )", main="voom: Mean-variance trend",
     pch = 16,cex.lab=2.5, cex.axis=2.5, cex.main=2.5, cex.sub=2.5, cex=2)
lines(l, col = "red", cex=2)
mtext(toupper(letters[i]), side = 3, adj = 0.05, line = 2.4,cex=3,font=4)
i <- i + 1
barplot(aw, names = 1:length(aw), main = "Sample-specific weights", 
        ylab = "Weight", xlab = "Sample",cex.lab=2.5, cex.axis=2.5, cex.main=2.5, cex.sub=2.5, cex=3.5)
abline(h = 1, col = 2, lty = 2)
mtext(toupper(letters[i]), side = 3, adj = 0.05, line = 2.4,cex=3,font=4)
i <- i + 1
plot(resP_quant_f$logFC,-log10(resP_quant_f$adj.P.Val),pch=16,cex.lab=2.5, 
     cex.axis=2.5, cex.main=2.5, cex.sub=2.5, cex=2,xlab="log2FC",ylab="-log10(FDR)",main="Volcano plot") 
points(resP_quant_f$logFC[abs(resP_quant_f$logFC)>2 & -log10(resP_quant_f$adj.P.Val)>6],
       -log10(resP_quant_f$adj.P.Val)[abs(resP_quant_f$logFC)>2 & -log10(resP_quant_f$adj.P.Val)>6],pch=16,col="red",
       cex.lab=2.5, cex.axis=2.5, cex.sub=2.5, cex=2)
abline(h=1.30103,lty=2)
abline(v=-1,lty=2)
abline(v=1,lty=2)
mtext(toupper(letters[i]), side = 3, adj = 0.05, line = 2.4,cex=3,font=4)
i <- i + 1

dev.off()

list_quant <- rownames(head(resP_quant_f[abs(resP_quant_f$logFC)>2,],12))
jpeg("comp_quant_tophits.jpg",width=2400,height=3200)
par(mfrow=c(4,3),mar=c(6,8,6,4))
for (i in 1:length(list_quant)){
  plot(density(as.numeric(vP_quant_f$E[list_quant[i],])[23:44]-as.numeric(vP_quant_f$E[list_quant[i],])[1:22]),
       col="blue",lwd=2.5,main="",xlab="log2FC",ylab="Density",cex.main=2.5, cex.lab=2.5, cex.axis=2.5, cex.sub=2.5, cex=4)
  legend("topright",legend=c(resP_quant_f$geneSymbol[rownames(resP_quant_f)==list_quant[i]],
                             paste("mean:",round(summary(as.numeric(vP_quant_f$E[list_quant[i],])[23:44]-as.numeric(vP_quant_f$E[list_quant[i],])[1:22])[4],digits=3)),
                             paste("3rd Q:",round(summary(as.numeric(vP_quant_f$E[list_quant[i],])[23:44]-as.numeric(vP_quant_f$E[list_quant[i],])[1:22])[5],digits=3)),
                             paste("max:",round(summary(as.numeric(vP_quant_f$E[list_quant[i],])[23:44]-as.numeric(vP_quant_f$E[list_quant[i],])[1:22])[6],digits=3))), cex=3)
  
}
dev.off()



##################
# Additional Tim #
##################

y_TMM_f <- estimateDisp(y_TMM_f,design_TMM_f)
plotBCV(y_TMM_f)
test <- DGEList(v_quant_f$E)
test <- estimateDisp(test,design_quant_f)
plotBCV(test)

# Orange is low
colscale <- color.scale(colSums(counts_TMM_f),extremes=c("orange","blue"))
jpeg("MDS_TMM_f.jpeg",width = 720,height = 720)
plotMDS(cpm(y_TMM_f,log=T),col=colscale,pch=c(1:19,35:38,1:19,35:38))
dev.off()
jpeg("MDS_quant_f.jpeg",width = 720,height = 720)
plotMDS(v_quant_f$E,col=colscale,pch=c(1:19,35:38,1:19,35:38))
dev.off()
jpeg("MDS_TMM_nf.jpeg",width = 720,height = 720)
plotMDS(cpm(y_TMM_nf,log=T),col=colscale,pch=c(1:19,35:38,1:19,35:38))
dev.off()
jpeg("MDS_quant_nf.jpeg",width = 720,height = 720)
plotMDS(v_quant_nf$E,col=colscale,pch=c(1:19,35:38,1:19,35:38))
dev.off()

colscale <- color.scale(colSums(counts_TMM_f),extremes=c("orange","blue"))
jpeg("MDS_TMM_f_10K.jpeg",width = 720,height = 720)
plotMDS(cpm(y_TMM_f,log=T)[apply(cpm(y_TMM_f,log=T),1,var)>sort(apply(cpm(y_TMM_f,log=T),1,var),decreasing=T)[10000],],col=colscale,pch=c(1:19,35:38,1:19,35:38))
dev.off()
jpeg("MDS_quant_f_10K.jpeg",width = 720,height = 720)
plotMDS(v_quant_f$E[apply(v_quant_f$E,1,var)>sort(apply(v_quant_f$E,1,var),decreasing=T)[10000],],col=colscale,pch=c(1:19,35:38,1:19,35:38))
dev.off()
jpeg("MDS_TMM_nf_10K.jpeg",width = 720,height = 720)
plotMDS(cpm(y_TMM_nf[apply(cpm(y_TMM_nf,log=T),1,var)>sort(apply(cpm(y_TMM_nf,log=T),1,var),decreasing=T)[10000],],log=T),col=colscale,pch=c(1:19,35:38,1:19,35:38))
dev.off()
jpeg("MDS_quant_nf_10K.jpeg",width = 720,height = 720)
plotMDS(v_quant_nf$E[apply(v_quant_nf$E,1,var)>sort(apply(v_quant_nf$E,1,var),decreasing=T)[10000],],col=colscale,pch=c(1:19,35:38,1:19,35:38))
dev.off()


jpeg("MDS_TMM_f_cvst.jpeg")
plotMDS(cpm(y_TMM_f,log=T),labels=substr(cvst,1,1))
dev.off()
jpeg("MDS_quant_f_cvst.jpeg")
plotMDS(v_quant_f$E,labels=substr(cvst,1,1))
dev.off()
jpeg("MDS_TMM_nf_cvst.jpeg")
plotMDS(cpm(y_TMM_nf,log=T),labels=substr(cvst,1,1))
dev.off()
jpeg("MDS_quant_nf_cvst.jpeg")
plotMDS(v_quant_nf$E,labels=substr(cvst,1,1))
dev.off()



jpeg("MDS_TMM_f_sampleweights.jpeg")
plotMDS(cpm(y_TMM_f,log=T),labels=round(v_TMM_f$sample.weights,digits=2))
dev.off()
jpeg("MDS_quant_f_sampleweights.jpeg")
plotMDS(v_quant_f$E,labels=round(v_quant_f$sample.weights,digits=2))
dev.off()
jpeg("MDS_TMM_nf_sampleweights.jpeg")
plotMDS(cpm(y_TMM_nf,log=T),labels=round(v_TMM_nf$sample.weights,digits=2))
dev.off()
jpeg("MDS_quant_nf_sampleweights.jpeg")
plotMDS(v_quant_nf$E,labels=round(v_quant_nf$sample.weights,digits=2))
dev.off()



# jpeg("comp_quant_quantnov_outliers.jpg",width=600,height=800)
# par(mfrow=c(4,3))
# for (i in 1:length(head(list_quantv))){
#   summary(as.numeric(counts_quant_f[head(list_quantv)[i],]))
#   
#   plot(density(as.numeric(vP_quant_f$E[head(list_quantv)[i],])[1:22]-as.numeric(vP_quant_f$E[head(list_quantv)[i],])[23:44]),
#        col="blue",main="Density")
#   legend("topright",legend=c(resP_quant_f$geneSymbol[head(list_quantnov)[i]],
#                              paste("mean",summary(as.numeric(vP_quant_f$E[head(list_quantv)[i],])[1:22]-as.numeric(vP_quant_f$E[head(list_quantv)[i],])[23:44])[4]),
#                              paste("3rd Q",summary(as.numeric(vP_quant_f$E[head(list_quantv)[i],])[1:22]-as.numeric(vP_quant_f$E[head(list_quantv)[i],])[23:44])[5]),
#                              paste("max",summary(as.numeric(vP_quant_f$E[head(list_quantv)[i],])[1:22]-as.numeric(vP_quant_f$E[head(list_quantv)[i],])[23:44])[6])))
#   
# }
# for (i in 1:length(head(list_quantnov))){
#   summary(as.numeric(counts_quant_f[head(list_quantnov)[i],]))
#   
#   plot(density(as.numeric(y_quant_f_nov[head(list_quantnov)[i],])[1:22]-as.numeric(y_quant_f_nov[head(list_quantnov)[i],])[23:44]),
#        col="blue",main="Density")
#   legend("topright",legend=c(resP_quant_f_nov$geneSymbol[head(list_quantnov)[i]],
#                              paste("mean",summary(as.numeric(y_quant_f_nov[head(list_quantnov)[i],])[1:22]-as.numeric(y_quant_f_nov[head(list_quantnov)[i],])[23:44])[4]),
#                              paste("3rd Q",summary(as.numeric(y_quant_f_nov[head(list_quantnov)[i],])[1:22]-as.numeric(y_quant_f_nov[head(list_quantnov)[i],])[23:44])[5]),
#                              paste("max",summary(as.numeric(y_quant_f_nov[head(list_quantnov)[i],])[1:22]-as.numeric(y_quant_f_nov[head(list_quantnov)[i],])[23:44])[6])))
#   
# }
# dev.off()

##
length(rownames(resP_quant_f_nosw[resP_quant_f_nosw$adj.P.Val<0.05,])[!(rownames(resP_quant_f_nosw[resP_quant_f_nosw$adj.P.Val<0.05,])%in%rownames(resP_quant_f_nov[resP_quant_f_nov$adj.P.Val<0.05,]))])
length(rownames(resP_quant_f_nov[resP_quant_f_nov$adj.P.Val<0.05,])[!(rownames(resP_quant_f_nov[resP_quant_f_nov$adj.P.Val<0.05,])%in%rownames(resP_quant_f_nosw[resP_quant_f_nosw$adj.P.Val<0.05,]))])
length(rownames(resP_quant_f_nov[resP_quant_f_nov$adj.P.Val<0.05,])[(rownames(resP_quant_f_nov[resP_quant_f_nov$adj.P.Val<0.05,])%in%rownames(resP_quant_f_nosw[resP_quant_f_nosw$adj.P.Val<0.05,]))])

list_quantv <- rownames(resP_quant_f_nosw[resP_quant_f_nosw$adj.P.Val<0.05,])[!(rownames(resP_quant_f_nosw[resP_quant_f_nosw$adj.P.Val<0.05,])%in%rownames(resP_quant_f_nov[resP_quant_f_nov$adj.P.Val<0.05,]))]
list_quantnov <- rownames(resP_quant_f_nov[resP_quant_f_nov$adj.P.Val<0.05,])[!(rownames(resP_quant_f_nov[resP_quant_f_nov$adj.P.Val<0.05,])%in%rownames(resP_quant_f_nosw[resP_quant_f_nosw$adj.P.Val<0.05,]))]
list_quantboth <- rownames(resP_quant_f_nov[resP_quant_f_nov$adj.P.Val<0.05,])[(rownames(resP_quant_f_nov[resP_quant_f_nov$adj.P.Val<0.05,])%in%rownames(resP_quant_f_nosw[resP_quant_f_nosw$adj.P.Val<0.05,]))]

logFC <- NULL
for (i in 1:nrow(vP_quant_f$E)){
  logFC <- rbind(logFC,as.numeric(vP_quant_f$E[i,1:22])-as.numeric(vP_quant_f$E[i,23:44]))
}
rownames(logFC) <- rownames(vP_quant_f$E)
dim(logFC)
head(logFC)

p_val_quantv <- NULL
for (i in 1:length(list_quantv)){
  p_val_quantv <- c(p_val_quantv,SIGN.test(logFC[list_quantv[i],],md=0)$p.value)
}

p_val_quantnov <- NULL
for (i in 1:length(list_quantnov)){
  p_val_quantnov <- c(p_val_quantnov,SIGN.test(logFC[list_quantnov[i],],md=0)$p.value)
}

p_val_quantboth <- NULL
for (i in 1:length(list_quantboth)){
  p_val_quantboth <- c(p_val_quantboth,SIGN.test(logFC[list_quantboth[i],],md=0)$p.value)
}

sni
# # Comparison
# length(rownames(resP_quant_f[resP_quant_f$adj.P.Val<0.05,])[!(rownames(resP_quant_f[resP_quant_f$adj.P.Val<0.05,])%in%rownames(resP_quant_f_nov[resP_quant_f_nov$adj.P.Val<0.05,]))])
# list_quant_vs_quantnov <- head(rownames(resP_quant_f[resP_quant_f$adj.P.Val<0.05,])[!(rownames(resP_quant_f[resP_quant_f$adj.P.Val<0.05,])%in%rownames(resP_quant_f_nov[resP_quant_f_nov$adj.P.Val<0.05,]))])
# 
# length(rownames(resP_quant_f_nov[resP_quant_f_nov$adj.P.Val<0.05,])[!(rownames(resP_quant_f_nov[resP_quant_f_nov$adj.P.Val<0.05,])%in%rownames(resP_quant_f[resP_quant_f$adj.P.Val<0.05,]))])
# list_quantnov_vs_quant <- head(rownames(resP_quant_f_nov[resP_quant_f_nov$adj.P.Val<0.05,])[!(rownames(resP_quant_f_nov[resP_quant_f_nov$adj.P.Val<0.05,])%in%rownames(resP_quant_f[resP_quant_f$adj.P.Val<0.05,]))])
# 
# jpeg("comp_quant_quantnov_outliers.jpg",width=600,height=800)
# par(mfrow=c(4,3))
# for (i in 1:length(list_quant_vs_quantnov)){
#   summary(as.numeric(counts_quant_f[list_quant_vs_quantnov[i],]))
#   
#   plot(density(as.numeric(counts_quant_f[list_quant_vs_quantnov[i],])[1:22]),col="orange",
#        ylim=c(0,max(c(density(as.numeric(counts_quant_f[list_quant_vs_quantnov[i],])[1:22])$y,
#                       density(as.numeric(counts_quant_f[list_quant_vs_quantnov[i],])[23:44])$y))),
#        xlim=c(0,max(c(density(as.numeric(counts_quant_f[list_quant_vs_quantnov[i],])[1:22])$x,
#                       density(as.numeric(counts_quant_f[list_quant_vs_quantnov[i],])[23:44])$x))),
#        main="Density")
#   lines(density(as.numeric(counts_quant_f[list_quant_vs_quantnov[i],])[23:44]),col="blue")
#   legend("topright",legend=c(paste("rank",which(rownames(resP_quant_f)==list_quant_vs_quantnov[i])),
#                              paste("mean",summary(as.numeric(counts_quant_f[list_quant_vs_quantnov[i],]))[4]),
#                              paste("3rd Q",summary(as.numeric(counts_quant_f[list_quant_vs_quantnov[i],]))[5]),
#                              paste("max",summary(as.numeric(counts_quant_f[list_quant_vs_quantnov[i],]))[6])))
#   
# }
# for (i in 1:length(list_quantnov_vs_quant)){
#   summary(as.numeric(counts_quant_f[list_quantnov_vs_quant[i],]))
#   
#   plot(density(as.numeric(counts_quant_f[list_quantnov_vs_quant[i],])[1:22]),col="orange",
#        ylim=c(0,max(c(density(as.numeric(counts_quant_f[list_quantnov_vs_quant[i],])[1:22])$y,
#                       density(as.numeric(counts_quant_f[list_quantnov_vs_quant[i],])[23:44])$y))),
#        xlim=c(0,max(c(density(as.numeric(counts_quant_f[list_quantnov_vs_quant[i],])[1:22])$x,
#                       density(as.numeric(counts_quant_f[list_quantnov_vs_quant[i],])[23:44])$x))),
#        main="Density")
#   lines(density(as.numeric(counts_quant_f[list_quantnov_vs_quant[i],])[23:44]),col="blue",main="Density")
#   legend("topright",legend=c(paste("rank",which(rownames(resP_quant_f_nov)==list_quantnov_vs_quant[i])),
#                              paste("mean",summary(as.numeric(counts_quant_f[list_quantnov_vs_quant[i],]))[4]),
#                              paste("3rd Q",summary(as.numeric(counts_quant_f[list_quantnov_vs_quant[i],]))[5]),
#                              paste("max",summary(as.numeric(counts_quant_f[list_quantnov_vs_quant[i],]))[6])))
#   
# }
# dev.off()

# # Comparison
# length(rownames(resP_TMM_f[resP_TMM_f$adj.P.Val<0.05,])[!(rownames(resP_TMM_f[resP_TMM_f$adj.P.Val<0.05,])%in%rownames(resP_TMM_f_nov[resP_TMM_f_nov$adj.P.Val<0.05,]))])
# list_TMM_vs_TMMnov <- head(rownames(resP_TMM_f[resP_TMM_f$adj.P.Val<0.05,])[!(rownames(resP_TMM_f[resP_TMM_f$adj.P.Val<0.05,])%in%rownames(resP_TMM_f_nov[resP_TMM_f_nov$adj.P.Val<0.05,]))])
# 
# length(rownames(resP_TMM_f_nov[resP_TMM_f_nov$adj.P.Val<0.05,])[!(rownames(resP_TMM_f_nov[resP_TMM_f_nov$adj.P.Val<0.05,])%in%rownames(resP_TMM_f[resP_TMM_f$adj.P.Val<0.05,]))])
# list_TMMnov_vs_TMM <- head(rownames(resP_TMM_f_nov[resP_TMM_f_nov$adj.P.Val<0.05,])[!(rownames(resP_TMM_f_nov[resP_TMM_f_nov$adj.P.Val<0.05,])%in%rownames(resP_TMM_f[resP_TMM_f$adj.P.Val<0.05,]))])
# 
# jpeg("comp_TMM_TMMnov_outliers.jpg",width=600,height=800)
# par(mfrow=c(4,3))
# for (i in 1:length(list_TMM_vs_TMMnov)){
#   summary(as.numeric(counts_TMM_f[list_TMM_vs_TMMnov[i],]))
#   
#   plot(density(as.numeric(counts_TMM_f[list_TMM_vs_TMMnov[i],])[1:22]),col="orange",
#        ylim=c(0,max(c(density(as.numeric(counts_TMM_f[list_TMM_vs_TMMnov[i],])[1:22])$y,
#                       density(as.numeric(counts_TMM_f[list_TMM_vs_TMMnov[i],])[23:44])$y))),
#        xlim=c(0,max(c(density(as.numeric(counts_TMM_f[list_TMM_vs_TMMnov[i],])[1:22])$x,
#                       density(as.numeric(counts_TMM_f[list_TMM_vs_TMMnov[i],])[23:44])$x))),
#        main="Density")
#   lines(density(as.numeric(counts_TMM_f[list_TMM_vs_TMMnov[i],])[23:44]),col="blue")
#   legend("topright",legend=c(paste("rank",which(rownames(resP_TMM_f)==list_TMM_vs_TMMnov[i])),
#                              paste("mean",summary(as.numeric(counts_TMM_f[list_TMM_vs_TMMnov[i],]))[4]),
#                              paste("3rd Q",summary(as.numeric(counts_TMM_f[list_TMM_vs_TMMnov[i],]))[5]),
#                              paste("max",summary(as.numeric(counts_TMM_f[list_TMM_vs_TMMnov[i],]))[6])))
#   
# }
# for (i in 1:length(list_TMMnov_vs_TMM)){
#   summary(as.numeric(counts_TMM_f[list_TMMnov_vs_TMM[i],]))
#   
#   plot(density(as.numeric(counts_TMM_f[list_TMMnov_vs_TMM[i],])[1:22]),col="orange",
#        ylim=c(0,max(c(density(as.numeric(counts_TMM_f[list_TMMnov_vs_TMM[i],])[1:22])$y,
#                       density(as.numeric(counts_TMM_f[list_TMMnov_vs_TMM[i],])[23:44])$y))),
#        xlim=c(0,max(c(density(as.numeric(counts_TMM_f[list_TMMnov_vs_TMM[i],])[1:22])$x,
#                       density(as.numeric(counts_TMM_f[list_TMMnov_vs_TMM[i],])[23:44])$x))),
#        main="Density")
#   lines(density(as.numeric(counts_TMM_f[list_TMMnov_vs_TMM[i],])[23:44]),col="blue")
#   legend("topright",legend=c(paste("rank",which(rownames(resP_TMM_f_nov)==list_TMMnov_vs_TMM[i])),
#                              paste("mean",summary(as.numeric(counts_TMM_f[list_TMMnov_vs_TMM[i],]))[4]),
#                              paste("3rd Q",summary(as.numeric(counts_TMM_f[list_TMMnov_vs_TMM[i],]))[5]),
#                              paste("max",summary(as.numeric(counts_TMM_f[list_TMMnov_vs_TMM[i],]))[6])))
#   
# }
# dev.off()

# 
# ## Comp TMM-quant outliers
# ##########################
# 
# length(rownames(resP_TMM_f[resP_TMM_f$adj.P.Val<0.05,])[!(rownames(resP_TMM_f[resP_TMM_f$adj.P.Val<0.05,])%in%rownames(resP_quant_f_nov[resP_quant_f$adj.P.Val<0.05,]))])
# list_TMM <- head(rownames(resP_TMM_f[resP_TMM_f$adj.P.Val<0.05,])[!(rownames(resP_TMM_f[resP_TMM_f$adj.P.Val<0.05,])%in%rownames(resP_quant_f_nov[resP_quant_f$adj.P.Val<0.05,]))])
# 
# length(rownames(resP_quant_f[resP_quant_f$adj.P.Val<0.05,])[!(rownames(resP_quant_f[resP_quant_f$adj.P.Val<0.05,])%in%rownames(resP_TMM_f[resP_TMM_f$adj.P.Val<0.05,]))])
# list_quant <- head(rownames(resP_quant_f[resP_quant_f$adj.P.Val<0.05,])[!(rownames(resP_quant_f[resP_quant_f$adj.P.Val<0.05,])%in%rownames(resP_TMM_f[resP_TMM_f$adj.P.Val<0.05,]))])
# 
# jpeg("comp_TMM_quant_outliers.jpg",width=600,height=800)
# par(mfrow=c(4,3))
# for (i in 1:length(list_TMM)){
#   summary(as.numeric(counts_TMM_f[list_TMM[i],]))
#   
#   plot(density(as.numeric(counts_TMM_f[list_TMM[i],])[1:22]),col="orange",
#        ylim=c(0,max(c(density(as.numeric(counts_TMM_f[list_TMM[i],])[1:22])$y,
#                       density(as.numeric(counts_TMM_f[list_TMM[i],])[23:44])$y))),
#        xlim=c(0,max(c(density(as.numeric(counts_TMM_f[list_TMM[i],])[1:22])$x,
#                       density(as.numeric(counts_TMM_f[list_TMM[i],])[23:44])$x))),
#        main="Density")
#   lines(density(as.numeric(counts_TMM_f[list_TMM[i],])[23:44]),col="blue")
#   legend("topright",legend=c(paste("rank",which(rownames(resP_TMM_f)==list_TMM[i])),
#                              paste("mean",summary(as.numeric(counts_TMM_f[list_TMM[i],]))[4]),
#                              paste("3rd Q",summary(as.numeric(counts_TMM_f[list_TMM[i],]))[5]),
#                              paste("max",summary(as.numeric(counts_TMM_f[list_TMM[i],]))[6])))
#          
# }
# for (i in 1:length(list_quant)){
#   summary(as.numeric(counts_quant_f[list_quant[i],]))
#   
#   plot(density(as.numeric(counts_quant_f[list_quant[i],])[1:22]),col="orange",
#        ylim=c(0,max(c(density(as.numeric(counts_quant_f[list_quant[i],])[1:22])$y,
#                       density(as.numeric(counts_quant_f[list_quant[i],])[23:44])$y))),
#        xlim=c(0,max(c(density(as.numeric(counts_quant_f[list_quant[i],])[1:22])$x,
#                       density(as.numeric(counts_quant_f[list_quant[i],])[23:44])$x))),
#        main="Density")
#   lines(density(as.numeric(counts_quant_f[list_quant[i],])[23:44]),col="blue")
#   legend("topright",legend=c(paste("rank",which(rownames(resP_quant_f)==list_quant[i])),
#                              paste("mean",summary(as.numeric(counts_quant_f[list_quant[i],]))[4]),
#                              paste("3rd Q",summary(as.numeric(counts_quant_f[list_quant[i],]))[5]),
#                              paste("max",summary(as.numeric(counts_quant_f[list_quant[i],]))[6])))
#   
# }
# dev.off()




## Venn Diagram
a<-vennCounts(cbind(resP_quant_f_nov$adj.P.Val[sort(rownames(resP_quant_f_nov),index.return=T)$ix]<0.05,
                    resP_quant_f_nosw$adj.P.Val[sort(rownames(resP_quant_f_nosw),index.return=T)$ix]<0.05,
                    resP_quant_f$adj.P.Val[sort(rownames(resP_quant_f),index.return=T)$ix]<0.05))
colnames(a) <- c("no SW nor VT","SW, no VT","SW + VT","Counts")
vennDiagram(a)

## Coverage vs Sample weights
round(rbind(colSums(counts_quant_f)/1000000,vP_quant_f$sample.weights),digits=2)
round(rbind(colSums(counts_TMM_f)/1000000,vP_TMM_f$sample.weights),digits=2)
round(rbind(colSums(counts_quant_f)/1000000,v_quant_f$sample.weights),digits=2)



# # Results
# # SIRPA
# # MME
# # EGFR
# # CHST3
# # DMD
# # TINAGL1
# # AIF1L
# # SPTBN1
# # ACACB
# # CRYAB
# 
# cvst_test <- as.character(cvst)
# cvst_test[cvst_test=="tumor"] <- "T"
# cvst_test[cvst_test=="normal"] <- "N"
# test <- counts_trim[rownames(counts_trim)%in%rownames(head(resP_f)),]
# as.character(cvst_test)[sort(unlist(test[1,]),index.return=T)$ix]
# 


#############
# Local FDR #
#############


localFDR <- fdrtool(resP_quant_f$P.Value,statistic = "pvalue",plot=T)
sum(fdrtool(resP_quant_f$P.Value,statistic = "pvalue",plot=T)$lfdr<0.05)

plot(1-localFDR$pval,localFDR$lfdr,type="l")
lines(1-resP_quant_f$P.Value,resP_quant_f$adj.P.Val,col="red")



