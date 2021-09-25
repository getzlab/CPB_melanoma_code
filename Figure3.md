Figure3.Rmd
================
Sam Freeman

``` {.r}
library(plyr)
library(ggplot2)
library(survival)
library(survminer)
library(rms)
library(psych)
library(grid)
library(gridExtra)
library(ggmosaic)
library(reshape2)
library(gplots)
library(RColorBrewer)
library(ggrepel)
```

Load R objects

``` {.r}
skcm_clin <- readRDS("figure3_inputs/skcm_clin.3.rds")
skcm_clin$surv <- Surv(skcm_clin$overall_survival,skcm_clin$os_status)

clin <- readRDS("figure3_inputs/rna_clin.3.rds")
clin$surv <- Surv(clin$overall_survival,clin$DOD=="DOD")

cbpalette_reordered <- c("#E69F00","#999999","#56B4E9","#009E73","#F0E442")
```

Functions

``` {.r}
formatSF <- function(n, digits, format) {
         if(missing(format)) {
                        return(formatC(signif(n,digits=digits), digits=digits,format="fg"))
         } else {
                        return(formatC(signif(n,digits=digits), digits=digits-1,format=format, flag="#"))
                        #return(formatC(signif(n,digits=digits), digits=digits,format=format))
         }
}

make_marker_plot <- function(df,gene,pdf_name) {
kt <- kruskal.test(df[,gene],df$plot_cluster)
kt_p <- formatSF(kt$p.value,digits=2,format="e")
g_skcm_nmf <- ggplot(df,aes_string(x="plot_cluster",y=gene),environment=environment()) + geom_boxplot(aes(fill=plot_cluster),outlier.shape = NA) + geom_jitter() + theme_classic() + ylab(gene) +
scale_fill_manual(values=cbpalette_reordered) + ggtitle(paste0("Kruskal-Wallis p =",kt_p)) + xlab("") + theme(legend.position="none")
ggsave(pdf_name,g_skcm_nmf)
}

makeKMplot_2group <- function(df,var1,pdf_name) {
#note: need to pass in clin[!is.na(clin$surv),] or class counts will be wrong
df$class <- df[,var1]
class_table <- table(df$class)
s_coxph <- summary(coxph(surv ~ df[,var1],data=df))
if (s_coxph$coef[1,2]<1) {
   df$class <- factor(df$class,levels=c(names(class_table)[2],names(class_table)[1]))
   class_table <- table(df$class)
   s_coxph <- summary(coxph(surv ~ class,data=df))
}
hr <- s_coxph$coef[1,2]
hr_lo <- s_coxph$conf.int[1,3]
hr_hi <- s_coxph$conf.int[1,4]
hr_string <- paste0("HR=",formatSF(hr,digits=3)," (95% CI, ",formatSF(hr_lo,digits=3),"-",formatSF(hr_hi,digits=3),")","\n",
names(class_table)[1],"=",class_table[[1]]," ",names(class_table)[2],"=",class_table[[2]])
pdf(pdf_name,onefile=F)
#fit<- do.call(survfit,list(formula = surv ~ class, data = df))
fit<- surv_fit(surv ~ class, data = df)
print(ggsurvplot(fit, data = df,pval=TRUE,pval.method=TRUE,legend=c(0.75,0.75),censor=TRUE,palette=c("blue","red"),title=hr_string))
dev.off()
}

makeKMplot_5group <- function(df,var1,pdf_name) {
#note: need to pass in clin[!is.na(clin$surv),] or class counts will be wrong
df$class <- df[,var1]
class_table <- table(df$class)
count_string <- paste0(names(class_table)[1],"=",class_table[[1]]," ",names(class_table)[2],"=",class_table[[2]],"\n",
names(class_table)[3],"=",class_table[[3]]," ",names(class_table)[4],"=",class_table[[4]]," ",names(class_table)[5],"=",class_table[[5]])
pdf(pdf_name,onefile=F)
#fit<- do.call(survfit,list(formula = surv ~ class, data = df))
fit<- surv_fit(surv ~ class, data = df)
print(ggsurvplot(fit, data = df,pval=TRUE,pval.method=TRUE,legend=c(0.75,0.75),censor=TRUE,palette=cbpalette_reordered,title=count_string))
dev.off()
}

makeKMplot_5group_withp <- function(df,var1,pdf_name) {
#note: need to pass in clin[!is.na(clin$surv),] or class counts will be wrong
df$class <- df[,var1]
class_table <- table(df$class)
s_coxph <- summary(coxph(surv ~ df[,var1],data=df))
count_string <- paste0(names(class_table)[1],"=",class_table[[1]]," ",names(class_table)[2],"=",class_table[[2]]," ",names(class_table)[3],"=",class_table[[3]],"\n",
names(class_table)[4],"=",class_table[[4]]," ",names(class_table)[5],"=",class_table[[5]],
" Log-rank p=",as.character(formatSF(s_coxph$sctest["pvalue"],digits=3,format="e")))
pdf(pdf_name,onefile=F)
#fit<- do.call(survfit,list(formula = surv ~ class, data = df))
fit<- surv_fit(surv ~ class, data = df)
print(ggsurvplot(fit, data = df,pval=TRUE,pval.method=TRUE,legend=c(0.75,0.75),censor=TRUE,palette=cbpalette_reordered,title=count_string))
dev.off()
}

plot.expr.fold.heatmap0 <- function(mat,scale0,cut.fold,g.Bayes) {
#Jaegil Kim
          scale <- 1
        color.axis <- 'black'
        .theme_ss <- theme_bw(base_size=12) +
                theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size=12*scale, family="mono",face='bold',color=color.axis),
                axis.text.y = element_text(hjust = 0.5,size=8*scale0, family="mono",face='bold',color=color.axis),
                axis.text = element_text(size = 12*scale, family = "mono",color=color.axis),
                axis.title=element_text(face="bold",size=12*scale,color="black"),
                plot.title=element_text(face="bold",size=12*scale))
        #hc <- hclust(dist(mat,method="euclidean"),method="ward.D")
        #gene.ordering <- hc$labels[hc$order]
        sample.ordering <- colnames(mat)
        #x <- mat[gene.ordering,]
        x <- mat
        x[x > +cut.fold] <- +cut.fold
        x[x < -cut.fold] <- -cut.fold
        y <- data.frame(rownames(x),x)
        colnames(y) <- c("gene",colnames(x))
        df <- reshape2::melt(y,id="gene")
        colnames(df) <- c("gene","sample","activity")
        df$gene <- factor(as.character(df$gene),levels=as.character(rev(rownames(mat))))
        df$sample <- factor(df$sample,levels=sample.ordering)
        p = ggplot(df,aes(x=sample,y=gene,fill=activity))+geom_tile() #geom_tile(colour="yellow")
        #p = p + facet_grid(. ~ signature, scale = "free_y")
        p = p + scale_fill_gradient2(low="blue",mid="black",high ="red",name=paste("Log2(Fold",'\n',"Changes)",sep=""))
        p = p + .theme_ss
        #p = p + ggtitle("Log2(Fold Changes)")
        p = p + xlab("Sample") + ylab("") #+ ylab("Differentially Expressed Genes")
        p = p + theme(axis.title.x = element_text(face="bold",colour="black",size=12*scale))
        p = p + theme(axis.title.y = element_text(face="bold",colour="black",size=12*scale))
        p = p + theme(axis.text.x=element_text(angle = 90,vjust=0.5,hjust=1,size=12*0))
        #p = p + theme(axis.text.y=element_text(angle = 00,vjust=0.5,hjust=0,size=12*scale))
        p = p + theme(legend.position="right")
        n.subtype <- table(g.Bayes)
        tmp <- n.subtype[1]
    n.subtype[1] <- n.subtype[2]
    n.subtype[2] <- tmp
    for (i in 1:(length(n.subtype)-1)) {
                p = p + geom_vline(xintercept=sum(n.subtype[1:i])+0.5,col='yellow')
        }
        #p = p + geom_vline(xintercept=sum(g.Bayes==1)+0.5,col='yellow')
    p = p + geom_vline(xintercept=sum(g.Bayes==2)+0.5,col='yellow')
        p = p + geom_vline(xintercept=sum(g.Bayes%in%c(1,2))+0.5,col='yellow')
        p = p + geom_vline(xintercept=sum(g.Bayes%in%c(1,2,3))+0.5,col='yellow')
        p = p + geom_vline(xintercept=sum(g.Bayes%in%c(1,2,3,4))+0.5,col='yellow')
        p = p + geom_vline(xintercept=sum(g.Bayes%in%c(1,2,3,4,5))+0.5,col='yellow')
        return(p)
}

get.sample.association.heatmap <- function(H,g.Bayes,scale) {
#Jaegil Kim
          color.axis <- 'black'

            g.ordering <- rev(c("G2","G1","G3","G4","G5"))
        #g.ordering <- c("G5","G4","G3","G2","G1")

        H.norm <- apply(H,2,function(x) x/sum(x))

    ph.order1 <- names(H.norm[1,g.Bayes==1][order(H.norm[1,g.Bayes==1],decreasing=T)])
    ph.order2 <- names(H.norm[2,g.Bayes==2][order(H.norm[2,g.Bayes==2],decreasing=T)])
    ph.order3 <- names(H.norm[3,g.Bayes==3][order(H.norm[3,g.Bayes==3],decreasing=T)])
    ph.order4 <- names(H.norm[4,g.Bayes==4][order(H.norm[4,g.Bayes==4],decreasing=T)])
    ph.order5 <- names(H.norm[5,g.Bayes==5][order(H.norm[5,g.Bayes==5],decreasing=T)])

    sample.ordering <- c(ph.order2,ph.order1,ph.order3,ph.order4,ph.order5)
        #sample.ordering <- colnames(H)[order(g.Bayes,decreasing=F)]
        df <- t(H)
        df <- reshape2::melt(df)
        colnames(df) <- c("sample","cluster","value")
        df$sample <- factor(df$sample,sample.ordering)
        df$cluster <- factor(df$cluster,g.ordering)
        df1 <- df
        df1[,"type"] <- "H matrix"

        df <- t(H.norm)
        df <- reshape2::melt(df) 
        colnames(df) <- c("sample","cluster","value")
        df$sample <- factor(df$sample,sample.ordering)
        df$cluster <- factor(df$cluster,g.ordering)
        df2 <- df
        df2[,"type"] <- "Normalized H"

        p = ggplot(df1,aes(x=sample,y=cluster,fill=value))+geom_tile() #geom_tile(colour="yellow")
        #p = p + facet_grid(. ~ type, scale = "free_y")
        p = p + scale_fill_gradient2(low="blue",mid="white",high ="red",name=paste("Activity",sep=""))
#        p = p + scale_fill_gradient2(low="blue",mid="white",high ="darkblue",name=paste("Activity",sep=""))
#        p = p + .theme_ss
        p = p + ggtitle("H matrix")
        p = p + xlab("Sample") + ylab("mRNA Clusters")
        p = p + theme(axis.title.x = element_text(face="bold",colour="black",size=14*scale))
        p = p + theme(axis.title.y = element_text(face="bold",colour="black",size=14*scale))
#         p = p + theme(axis.title.x = element_text(face="bold",colour="black",size=14))
#         p = p + theme(axis.title.y = element_text(face="bold",colour="black",size=14))
        p = p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size=5*scale, family="mono",face='bold',color=color.axis))
        p1 = p + theme(legend.position="top")

        p = ggplot(df2,aes(x=sample,y=cluster,fill=value))+geom_tile() #geom_tile(colour="yellow")
        #p = p + facet_grid(. ~ type, scale = "free_y")
        p = p + scale_fill_gradient2(low="blue",mid="black",high ="red",name=paste("Activity",sep=""))
#        p = p + scale_fill_gradient2(low="blue",mid="white",high ="darkblue",name=paste("Activity",sep=""))
#        p = p + .theme_ss
        p = p + ggtitle("Normalized H matrix")
        p = p + xlab("Sample") + ylab("mRNA Clusters")
        p = p + theme(axis.title.x = element_text(face="bold",colour="black",size=14*scale))
        p = p + theme(axis.title.y = element_text(face="bold",colour="black",size=14*scale))
#        p = p + theme(axis.title.x = element_text(face="bold",colour="black",size=14))
#        p = p + theme(axis.title.y = element_text(face="bold",colour="black",size=14))
        p = p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size=5*scale, family="mono",face='bold',color=color.axis))
        p2 = p + theme(legend.position="top")

        df <- data.frame(t(H))
        colnames(df) <- paste("H",seq(1:nrow(H)),sep="")
        df[,"mRNA"] <- g.Bayes
        df$mRNA <- paste("G",df$mRNA,sep="")
        df <- reshape2::melt(df,id="mRNA")
        colnames(df) <- c("Subtype","cluster","Association")
        p = ggplot(df,aes(x=Association)) 
        p = p + geom_histogram(color="black",fill="gray75",binwidth=0.05) 
        p = p + facet_grid(Subtype ~ cluster,scale='free_y')
#        p = p + .theme_ss
        p = p + ggtitle("Sample Associations")
        p = p + xlab("Association") + ylab("Sample Counts")
        p = p + theme(axis.title.x = element_text(face="bold",colour="black",size=14*scale))
        p = p + theme(axis.title.y = element_text(face="bold",colour="black",size=14*scale))
        p = p + theme(axis.text.x = element_text(angle=0,vjust=0.0, size=12*scale, family="mono",face='bold',color=color.axis))
        p3 = p + theme(legend.position="top")
        return(list(p1,p2,p3))
}
```

Figure 3

``` {.r}
dir.create("figure3_outputs", showWarnings = FALSE)

makeKMplot_5group_withp(skcm_clin[!is.na(skcm_clin$surv)&skcm_clin$stage34==TRUE,],"plot_cluster","figure3_outputs/3a.skcm_34_subtype_survival.pdf")
```

    ## pdf 
    ##   2

``` {.r}
knitr::include_graphics("figure3_outputs/3a.skcm_34_subtype_survival.pdf")
```

<embed src="figure3_outputs/3a.skcm_34_subtype_survival.pdf" width="500px" height="500px" type="application/pdf" />

``` {.r}
load(file=paste(paste("figure3_inputs/res.L1EU.Bayes.2.RData",sep="."),sep="")) #### BayesNMF ouput used for the de-novo expression subtyping

res <- res.Bayes
W <- res[[1]]
H <- res[[2]]
W <- W[,colSums(W)!=0]
H <- H[rowSums(H)!=0,]
rownames(H) <- paste("G",seq(1:nrow(H)),sep="")
H.norm <- apply(H,2,function(x) x/sum(x))
g.Bayes <- apply(H.norm,2,function(x) which.max(x))
g.tmp <- g.Bayes  ## from the classification
#H.tmp <- H.Bayes  ## from the classification
H.tmp <- H
H.tmp.norm <- H.norm
g1 <- names(g.tmp)[g.tmp==1]
g2 <- names(g.tmp)[g.tmp==2]
g3 <- names(g.tmp)[g.tmp==3]
g4 <- names(g.tmp)[g.tmp==4]
g5 <- names(g.tmp)[g.tmp==5]
order1 <- g1[order(H.tmp.norm[1,g.tmp==1],decreasing=T)]
order2 <- g2[order(H.tmp.norm[2,g.tmp==2],decreasing=T)]
order3 <- g3[order(H.tmp.norm[3,g.tmp==3],decreasing=T)]
order4 <- g4[order(H.tmp.norm[4,g.tmp==4],decreasing=T)]
order5 <- g5[order(H.tmp.norm[5,g.tmp==5],decreasing=T)]
sample0 <- c(order2,order1,order3,order4,order5)
ordering <- match(sample0,names(g.tmp),nomatch=0)

gene1.marker <- read.table("figure3_inputs/SKCM.1.markers.txt",as.is=T)$x
gene2.marker <- read.table("figure3_inputs/SKCM.2.markers.txt",as.is=T)$x
gene3.marker <- read.table("figure3_inputs/SKCM.3.markers.txt",as.is=T)$x
gene4.marker <- read.table("figure3_inputs/SKCM.4.markers.txt",as.is=T)$x
gene5.marker <- read.table("figure3_inputs/SKCM.5.markers.txt",as.is=T)$x

cut.NA <- 0.1
OUTPUT <- "figure3_outputs/"
cohort <- "SKCM"

load(file=paste("figure3_inputs/tcga.gencode_v19.coding.tpm.txt.RData",sep=""))
#loads log2.RSEM.SKCM
mRNA.comm <- log2.RSEM.SKCM
index.gene <- rowSums(is.na(mRNA.comm))<round(cut.NA*ncol(mRNA.comm))
mRNA.comm <- mRNA.comm[index.gene,]
mRNA.fold <- t(apply(mRNA.comm,1,function(x) x-median(x,na.rm=T)))
expr.fold <- mRNA.fold

x <- c(gene2.marker,gene1.marker,gene3.marker,gene4.marker,gene5.marker) ## marker genes 
x <- x[x%in%rownames(expr.fold)]  ### expr.fold=fold-change based expression matrix
x0 <- expr.fold[match(x,rownames(expr.fold),nomatch=0),ordering]

gencode <- read.delim(paste("figure3_inputs/gencode.v19.genes.v7.patched_contigs.bed",sep=""),header=T,sep='\t',as.is=T,comment="#")
rownames(x0) <- gencode$Gene[match(rownames(x0),gencode$Gene.id,nomatch=0)]

pdf(file=paste(OUTPUT,paste("S16a",cohort,"marker0.expr.fold.Bayes",0,"ordered.full.pdf",sep="."),sep=""),width=12,height=12)
        p <- plot.expr.fold.heatmap0(x0,0.5,3.0,g.Bayes[ordering])
        plot(p)
dev.off()
```

    ## pdf 
    ##   2

``` {.r}
load("figure3_inputs/fcmeta.H.Bayes.RData")
#loads H.Bayes.new
H.norm.new <- apply(H.Bayes.new,2,function(x) x/sum(x))
g.Bayes.new <- apply(H.norm.new,2,function(x) which.max(x))

g1.new <- names(g.Bayes.new)[g.Bayes.new==1]
g2.new <- names(g.Bayes.new)[g.Bayes.new==2]
g3.new <- names(g.Bayes.new)[g.Bayes.new==3]
g4.new <- names(g.Bayes.new)[g.Bayes.new==4]
g5.new <- names(g.Bayes.new)[g.Bayes.new==5]
order1.new <- g1.new[order(H.norm.new[1,g.Bayes.new==1],decreasing=T)]
order2.new <- g2.new[order(H.norm.new[2,g.Bayes.new==2],decreasing=T)]
order3.new <- g3.new[order(H.norm.new[3,g.Bayes.new==3],decreasing=T)]
order4.new <- g4.new[order(H.norm.new[4,g.Bayes.new==4],decreasing=T)]
order5.new <- g5.new[order(H.norm.new[5,g.Bayes.new==5],decreasing=T)]
sample0.new <- c(order2.new,order1.new,order3.new,order4.new,order5.new)
ordering.new <- match(sample0.new,names(g.Bayes.new),nomatch=0)

cohort <- "CPB"

expr <- readRDS("figure3_inputs/cpb_primary_154_combat_corrected_log2TPM.rds")
expr[expr==0] <- NA

cut.NA <- 0.1
index.NA <- rowSums(is.na(expr))<= round(cut.NA*ncol(expr))
expr <- expr[index.NA,]
expr.fold <- t(apply(expr,1,function(x) x-median(x,na.rm=T)))

x <- c(gene2.marker,gene1.marker,gene3.marker,gene4.marker,gene5.marker) ## marker genes
x <- x[x%in%rownames(expr.fold)]  ### expr.fold=fold-change based expression matrix
x0.new <- expr.fold[match(x,rownames(expr.fold),nomatch=0),ordering.new]
rownames(x0.new) <- gencode$Gene[match(rownames(x0.new),gencode$Gene.id,nomatch=0)]

pdf(file=paste(OUTPUT,paste("3b",cohort,"marker0.expr.fold.Bayes",0.5,"ordered.full.pdf",sep="."),sep=""),width=12,height=12)
        p <- plot.expr.fold.heatmap0(x0.new,0.5,3,g.Bayes.new[ordering.new])
        plot(p)
dev.off()
```

    ## pdf 
    ##   2

``` {.r}
knitr::include_graphics("figure3_outputs/3b.CPB.marker0.expr.fold.Bayes.0.5.ordered.full.pdf")
```

<embed src="figure3_outputs/3b.CPB.marker0.expr.fold.Bayes.0.5.ordered.full.pdf" width="500px" height="500px" type="application/pdf" />

``` {.r}
makeKMplot_5group(clin[!is.na(clin$surv),],"plot_cluster","figure3_outputs/3c.cpb_subtype_survival.pdf")
```

    ## pdf 
    ##   2

``` {.r}
knitr::include_graphics("figure3_outputs/3c.cpb_subtype_survival.pdf")
```

<embed src="figure3_outputs/3c.cpb_subtype_survival.pdf" width="500px" height="500px" type="application/pdf" />

``` {.r}
clin$plot_Immune <- ifelse(clin$plot_cluster=="Immune","Immune","Others")
makeKMplot_2group(clin[!is.na(clin$surv),],"plot_Immune","figure3_outputs/3d.cpb_immune_survival.pdf")
```

    ## pdf 
    ##   2

``` {.r}
knitr::include_graphics("figure3_outputs/3d.cpb_immune_survival.pdf")
```

<embed src="figure3_outputs/3d.cpb_immune_survival.pdf" width="500px" height="500px" type="application/pdf" />

``` {.r}
library("tximport")
library("readr")
library("tximportData")
library("DESeq2")

load("figure3_inputs/txi.rsem.coding.RData")
#loads txi.rsem.coding

gencode_coding <- read.delim("figure4_inputs/gencode_coding_withentrez.txt",header=T,as.is=T)

clin$OS_status <- ifelse(clin$overall_survival>=365,"high OS","low OS")
clin_nona <- clin[!is.na(clin$OS_status),]

col_count_data_all_os <- clin_nona[,c("rna_id","cohort","OS_status")]
txi.rsem.coding.all.os <- txi.rsem.coding
txi.rsem.coding.all.os$abundance <- txi.rsem.coding.all.os$abundance[,rownames(col_count_data_all_os)]
txi.rsem.coding.all.os$counts <- txi.rsem.coding.all.os$counts[,rownames(col_count_data_all_os)]
txi.rsem.coding.all.os$length <- txi.rsem.coding.all.os$length[,rownames(col_count_data_all_os)]
txi.rsem.coding.all.os$length[txi.rsem.coding.all.os$length == 0] <- 1

dds_all_os <- DESeqDataSetFromTximport(txi= txi.rsem.coding.all.os, colData = col_count_data_all_os, design = ~ cohort + OS_status)
```

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

``` {.r}
dds_all_os <- DESeq(dds_all_os)
res_all_os <- results(dds_all_os, coef="OS_status_low.OS_vs_high.OS")
df.res.all.os <- data.frame(res_all_os,stringsAsFactors=F)
df.res.all.os$gene <- gencode_coding[rownames(df.res.all.os),]$Gene

tab <- readRDS("figure3_inputs/cpb_primary_154_combat_corrected_log2TPM.rds")
#tab has all samples

use_mediantpm1_all_os <- names(apply(tab,1,median)[apply(tab,1,median)>1])
tmpdf_deseq_mediantpm1filtered_all_os <- data.frame(res_all_os[use_mediantpm1_all_os,],gene=gencode_coding[rownames(res_all_os[use_mediantpm1_all_os,]),]$Gene,stringsAsFactors=F)
tmpdf_deseq_mediantpm1filtered_all_os$new_padj <- p.adjust(tmpdf_deseq_mediantpm1filtered_all_os$pvalue,method="fdr")

all_deseq_logpval_threshold_os <- 0.001
all_deseq_pval_threshold_os <- 3
all_deseq_lfc_threshold_os <- 0.5

g_all_volcano_os <- ggplot(tmpdf_deseq_mediantpm1filtered_all_os,aes(x=-1*log2FoldChange,y=-log10(pvalue))) + geom_point(color="#999999") + theme_bw() + theme(legend.position = "none") + geom_text_repel(data=tmpdf_deseq_mediantpm1filtered_all_os[tmpdf_deseq_mediantpm1filtered_all_os$pvalue<all_deseq_logpval_threshold_os&abs(tmpdf_deseq_mediantpm1filtered_all_os$log2FoldChange)> all_deseq_lfc_threshold_os,],aes(label=gene,x=-1*log2FoldChange,y=-1*log10(pvalue),color=log2FoldChange>0),max.overlaps=100) + geom_point(data=tmpdf_deseq_mediantpm1filtered_all_os[tmpdf_deseq_mediantpm1filtered_all_os$pvalue<all_deseq_logpval_threshold_os&abs(tmpdf_deseq_mediantpm1filtered_all_os$log2FoldChange)> all_deseq_lfc_threshold_os,],aes(x=-1*log2FoldChange,y=-1*log10(pvalue),color=log2FoldChange>0)) + scale_color_manual(values=c("blue","red")) + geom_vline(xintercept=c(-1* all_deseq_lfc_threshold_os, all_deseq_lfc_threshold_os), linetype="dotted") + geom_hline(yintercept=all_deseq_pval_threshold_os, linetype="dotted") + ylab("-log10(p value)") + xlab("log2(Fold Change)")
ggsave("figure3_outputs/3e.deseq2_all_volcano_os_withlabels.pdf",g_all_volcano_os,height=7,width=7)
knitr::include_graphics("figure3_outputs/3e.deseq2_all_volcano_os_withlabels.pdf")
```

<embed src="figure3_outputs/3e.deseq2_all_volcano_os_withlabels.pdf" width="500px" height="500px" type="application/pdf" />

``` {.r}
print(dim(tmpdf_deseq_mediantpm1filtered_all_os[tmpdf_deseq_mediantpm1filtered_all_os$new_padj<0.05,]))
```

    ## [1] 83  8

``` {.r}
print(table(tmpdf_deseq_mediantpm1filtered_all_os[tmpdf_deseq_mediantpm1filtered_all_os$new_padj<0.05,]$log2FoldChange>0))
```

    ## 
    ## FALSE  TRUE 
    ##    55    28

Figure S15

``` {.r}
plist_SKCM_original <- get.sample.association.heatmap(H,g.Bayes,scale=1)
ggsave("figure3_outputs/S15a.SKCM_H_norm_newcolors.pdf",plist_SKCM_original[[2]],height=7,width=7)
knitr::include_graphics("figure3_outputs/S15a.SKCM_H_norm_newcolors.pdf")
```

<embed src="figure3_outputs/S15a.SKCM_H_norm_newcolors.pdf" width="500px" height="500px" type="application/pdf" />

``` {.r}
pdf("figure3_outputs/S15b.skcm_cluster_vs_oldTCGA_cluster.pdf")
grid.table(table(skcm_clin$plot_TCGA_cluster,skcm_clin$plot_cluster))
dev.off()
```

    ## pdf 
    ##   2

``` {.r}
knitr::include_graphics("figure3_outputs/S15b.skcm_cluster_vs_oldTCGA_cluster.pdf")
```

<embed src="figure3_outputs/S15b.skcm_cluster_vs_oldTCGA_cluster.pdf" width="500px" height="500px" type="application/pdf" />

``` {.r}
pdf("figure3_outputs/S15c.skcm_cluster_site.pdf",height=7,width=10)
grid.table(table(skcm_clin$Tissue.Site,skcm_clin$plot_cluster))
dev.off()
```

    ## pdf 
    ##   2

``` {.r}
knitr::include_graphics("figure3_outputs/S15c.skcm_cluster_site.pdf")
```

<embed src="figure3_outputs/S15c.skcm_cluster_site.pdf" width="500px" height="500px" type="application/pdf" />

``` {.r}
pdf("figure3_outputs/S15dskcm_cluster_vs_tsoi_cluster.pdf")
grid.table(table(skcm_clin$tsoi_subtype,skcm_clin$plot_cluster))
dev.off()
```

    ## pdf 
    ##   2

``` {.r}
knitr::include_graphics("figure3_outputs/S15dskcm_cluster_vs_tsoi_cluster.pdf")
```

<embed src="figure3_outputs/S15dskcm_cluster_vs_tsoi_cluster.pdf" width="500px" height="500px" type="application/pdf" />

Figure S16

``` {.r}
knitr::include_graphics("figure3_outputs/S16a.SKCM.marker0.expr.fold.Bayes.0.ordered.full.pdf")
```

<embed src="figure3_outputs/S16a.SKCM.marker0.expr.fold.Bayes.0.ordered.full.pdf" width="500px" height="500px" type="application/pdf" />

``` {.r}
#MITF high
make_marker_plot(skcm_clin,"PMEL","figure3_outputs/S16b.1_1_skcm_nmf_cluster_PMEL.pdf")
make_marker_plot(skcm_clin,"TRPM1","figure3_outputs/S16b.1_2_skcm_nmf_cluster_TRPM1.pdf")
make_marker_plot(skcm_clin,"MITF","figure3_outputs/S16b.1_3_skcm_nmf_cluster_MITF.pdf")
make_marker_plot(skcm_clin,"MLANA","figure3_outputs/S16b.1_4_skcm_nmf_cluster_MLANA.pdf")
knitr::include_graphics("figure3_outputs/S16b.1_1_skcm_nmf_cluster_PMEL.pdf")
```

<embed src="figure3_outputs/S16b.1_1_skcm_nmf_cluster_PMEL.pdf" width="500px" height="500px" type="application/pdf" />

``` {.r}
knitr::include_graphics("figure3_outputs/S16b.1_2_skcm_nmf_cluster_TRPM1.pdf")
```

<embed src="figure3_outputs/S16b.1_2_skcm_nmf_cluster_TRPM1.pdf" width="500px" height="500px" type="application/pdf" />

``` {.r}
knitr::include_graphics("figure3_outputs/S16b.1_3_skcm_nmf_cluster_MITF.pdf")
```

<embed src="figure3_outputs/S16b.1_3_skcm_nmf_cluster_MITF.pdf" width="500px" height="500px" type="application/pdf" />

``` {.r}
knitr::include_graphics("figure3_outputs/S16b.1_4_skcm_nmf_cluster_MLANA.pdf")
```

<embed src="figure3_outputs/S16b.1_4_skcm_nmf_cluster_MLANA.pdf" width="500px" height="500px" type="application/pdf" />

``` {.r}
#Intermediate
make_marker_plot(skcm_clin,"EDNRB","figure3_outputs/S16b.2_1_skcm_nmf_cluster_EDNRB.pdf")
make_marker_plot(skcm_clin,"SOX6","figure3_outputs/S16b.2_2_skcm_nmf_cluster_SOX6.pdf")
make_marker_plot(skcm_clin,"PAX3","figure3_outputs/S16b.2_3_skcm_nmf_cluster_PAX3.pdf")
make_marker_plot(skcm_clin,"DCT","figure3_outputs/S16b.2_4_skcm_nmf_cluster_DCT.pdf")
knitr::include_graphics("figure3_outputs/S16b.2_1_skcm_nmf_cluster_EDNRB.pdf")
```

<embed src="figure3_outputs/S16b.2_1_skcm_nmf_cluster_EDNRB.pdf" width="500px" height="500px" type="application/pdf" />

``` {.r}
knitr::include_graphics("figure3_outputs/S16b.2_2_skcm_nmf_cluster_SOX6.pdf")
```

<embed src="figure3_outputs/S16b.2_2_skcm_nmf_cluster_SOX6.pdf" width="500px" height="500px" type="application/pdf" />

``` {.r}
knitr::include_graphics("figure3_outputs/S16b.2_3_skcm_nmf_cluster_PAX3.pdf")
```

<embed src="figure3_outputs/S16b.2_3_skcm_nmf_cluster_PAX3.pdf" width="500px" height="500px" type="application/pdf" />

``` {.r}
knitr::include_graphics("figure3_outputs/S16b.2_4_skcm_nmf_cluster_DCT.pdf")
```

<embed src="figure3_outputs/S16b.2_4_skcm_nmf_cluster_DCT.pdf" width="500px" height="500px" type="application/pdf" />

``` {.r}
#MITF low
make_marker_plot(skcm_clin,"AXL","figure3_outputs/S16b.3_1_skcm_nmf_cluster_AXL.pdf")
make_marker_plot(skcm_clin,"TGFBI","figure3_outputs/S16b.3_2_skcm_nmf_cluster_TGFBI.pdf")
make_marker_plot(skcm_clin,"NGFR","figure3_outputs/S16b.3_3_skcm_nmf_cluster_NGFR.pdf")
make_marker_plot(skcm_clin,"TGFA","figure3_outputs/S16b.3_4_skcm_nmf_cluster_TGFA.pdf")
knitr::include_graphics("figure3_outputs/S16b.3_1_skcm_nmf_cluster_AXL.pdf")
```

<embed src="figure3_outputs/S16b.3_1_skcm_nmf_cluster_AXL.pdf" width="500px" height="500px" type="application/pdf" />

``` {.r}
knitr::include_graphics("figure3_outputs/S16b.3_2_skcm_nmf_cluster_TGFBI.pdf")
```

<embed src="figure3_outputs/S16b.3_2_skcm_nmf_cluster_TGFBI.pdf" width="500px" height="500px" type="application/pdf" />

``` {.r}
knitr::include_graphics("figure3_outputs/S16b.3_3_skcm_nmf_cluster_NGFR.pdf")
```

<embed src="figure3_outputs/S16b.3_3_skcm_nmf_cluster_NGFR.pdf" width="500px" height="500px" type="application/pdf" />

``` {.r}
knitr::include_graphics("figure3_outputs/S16b.3_4_skcm_nmf_cluster_TGFA.pdf")
```

<embed src="figure3_outputs/S16b.3_4_skcm_nmf_cluster_TGFA.pdf" width="500px" height="500px" type="application/pdf" />

``` {.r}
#Immune
make_marker_plot(skcm_clin,"CD2","figure3_outputs/S16b.4_1_skcm_nmf_cluster_CD2.pdf")
make_marker_plot(skcm_clin,"MS4A1","figure3_outputs/S16b.4_2_skcm_nmf_cluster_MS4A1.pdf")
make_marker_plot(skcm_clin,"PTPRC","figure3_outputs/S16b.4_3_skcm_nmf_cluster_PTPRC.pdf")
make_marker_plot(skcm_clin,"PDCD1","figure3_outputs/S16b.4_4_skcm_nmf_cluster_PDCD1.pdf")
knitr::include_graphics("figure3_outputs/S16b.4_1_skcm_nmf_cluster_CD2.pdf")
```

<embed src="figure3_outputs/S16b.4_1_skcm_nmf_cluster_CD2.pdf" width="500px" height="500px" type="application/pdf" />

``` {.r}
knitr::include_graphics("figure3_outputs/S16b.4_2_skcm_nmf_cluster_MS4A1.pdf")
```

<embed src="figure3_outputs/S16b.4_2_skcm_nmf_cluster_MS4A1.pdf" width="500px" height="500px" type="application/pdf" />

``` {.r}
knitr::include_graphics("figure3_outputs/S16b.4_3_skcm_nmf_cluster_PTPRC.pdf")
```

<embed src="figure3_outputs/S16b.4_3_skcm_nmf_cluster_PTPRC.pdf" width="500px" height="500px" type="application/pdf" />

``` {.r}
knitr::include_graphics("figure3_outputs/S16b.4_4_skcm_nmf_cluster_PDCD1.pdf")
```

<embed src="figure3_outputs/S16b.4_4_skcm_nmf_cluster_PDCD1.pdf" width="500px" height="500px" type="application/pdf" />

``` {.r}
#keratin high
make_marker_plot(skcm_clin,"KRT10","figure3_outputs/S16b.5_1_skcm_nmf_cluster_KRT10.pdf")
make_marker_plot(skcm_clin,"KRT5","figure3_outputs/S16b.5_2_skcm_nmf_cluster_KRT5.pdf")
make_marker_plot(skcm_clin,"DMKN","figure3_outputs/S16b.5_3_skcm_nmf_cluster_DMKN.pdf")
make_marker_plot(skcm_clin,"DEGS2","figure3_outputs/S16b.5_4_skcm_nmf_cluster_DEGS2.pdf")
knitr::include_graphics("figure3_outputs/S16b.5_1_skcm_nmf_cluster_KRT10.pdf")
```

<embed src="figure3_outputs/S16b.5_1_skcm_nmf_cluster_KRT10.pdf" width="500px" height="500px" type="application/pdf" />

``` {.r}
knitr::include_graphics("figure3_outputs/S16b.5_2_skcm_nmf_cluster_KRT5.pdf")
```

<embed src="figure3_outputs/S16b.5_2_skcm_nmf_cluster_KRT5.pdf" width="500px" height="500px" type="application/pdf" />

``` {.r}
knitr::include_graphics("figure3_outputs/S16b.5_3_skcm_nmf_cluster_DMKN.pdf")
```

<embed src="figure3_outputs/S16b.5_3_skcm_nmf_cluster_DMKN.pdf" width="500px" height="500px" type="application/pdf" />

``` {.r}
knitr::include_graphics("figure3_outputs/S16b.5_4_skcm_nmf_cluster_DEGS2.pdf")
```

<embed src="figure3_outputs/S16b.5_4_skcm_nmf_cluster_DEGS2.pdf" width="500px" height="500px" type="application/pdf" />

``` {.r}
kt <- kruskal.test(log10(skcm_clin$nonsilent_snpdnpindel+1),skcm_clin$plot_cluster)
kt_p <- formatSF(kt$p.value,format="e",digits=2)
g_skcm_nmf_tmb <- ggplot(skcm_clin,aes(x=plot_cluster,y=nonsilent_snpdnpindel+1)) + geom_boxplot(aes(fill=plot_cluster),outlier.shape = NA) + geom_jitter() + scale_y_log10() + ylab("TMB") +
theme_classic() + scale_fill_manual(values= cbpalette_reordered) + xlab("") + annotate(geom="text",x=2.5,y=5000,label=paste0("Kruskal-Wallis p =",kt_p))
ggsave("figure3_outputs/S16c.skcm_nmf_cluster_tmb.pdf",g_skcm_nmf_tmb)
knitr::include_graphics("figure3_outputs/S16c.skcm_nmf_cluster_tmb.pdf")
```

<embed src="figure3_outputs/S16c.skcm_nmf_cluster_tmb.pdf" width="500px" height="500px" type="application/pdf" />

``` {.r}
kt <- kruskal.test(skcm_clin[!is.na(skcm_clin$purity),]$purity,skcm_clin[!is.na(skcm_clin$purity),]$plot_cluster)
kt_p <- formatSF(kt$p.value,format="e",digits=2)
g_skcm_nmf_purity <- ggplot(skcm_clin[!is.na(skcm_clin$purity),],aes(x=plot_cluster,y=purity)) + geom_boxplot(aes(fill=plot_cluster),outlier.shape = NA) + geom_jitter() + ylab("Tumor purity") +
theme_classic() + scale_fill_manual(values= cbpalette_reordered) + xlab("") + annotate(geom="text",x=2.5,y=0.25,label=paste0("Kruskal-Wallis p =",kt_p))
ggsave("figure3_outputs/S16d.skcm_nmf_cluster_purity.pdf",g_skcm_nmf_purity)
knitr::include_graphics("figure3_outputs/S16d.skcm_nmf_cluster_purity.pdf")
```

<embed src="figure3_outputs/S16d.skcm_nmf_cluster_purity.pdf" width="500px" height="500px" type="application/pdf" />

``` {.r}
kt <- kruskal.test(skcm_clin$rna_tcb,skcm_clin$plot_cluster)
kt_p <- formatSF(kt$p.value,format="e",digits=2)
g_skcm_nmf_tcb <- ggplot(skcm_clin,aes(x=plot_cluster,y=rna_tcb)) + geom_boxplot(aes(fill=plot_cluster),outlier.shape = NA) + geom_jitter() + ylab(expression(TCB[RNA])) + theme_classic() +
scale_fill_manual(values= cbpalette_reordered) + xlab("") + scale_y_log10() + annotate(geom="text",x=2.5,y=10,label=paste0("Kruskal-Wallis p =",kt_p))
ggsave("figure3_outputs/S16e.skcm_nmf_cluster_tcb.pdf",g_skcm_nmf_tcb)
knitr::include_graphics("figure3_outputs/S16e.skcm_nmf_cluster_tcb.pdf")
```

<embed src="figure3_outputs/S16e.skcm_nmf_cluster_tcb.pdf" width="500px" height="500px" type="application/pdf" />

``` {.r}
kt <- kruskal.test(skcm_clin$rna_bcb,skcm_clin$plot_cluster)
kt_p <- formatSF(kt$p.value,format="e",digits=2)
g_skcm_nmf_bcb <- ggplot(skcm_clin,aes(x=plot_cluster,y=rna_bcb)) + geom_boxplot(aes(fill=plot_cluster),outlier.shape = NA) + geom_jitter() + ylab(expression(BCB[DNA])) + theme_classic() +
scale_fill_manual(values= cbpalette_reordered) + xlab("") + scale_y_log10() + annotate(geom="text",x=2.5,y=3000,label=paste0("Kruskal-Wallis p =",kt_p))
ggsave("figure3_outputs/S16f.skcm_nmf_cluster_bcb.pdf",g_skcm_nmf_bcb)
knitr::include_graphics("figure3_outputs/S16f.skcm_nmf_cluster_bcb.pdf")
```

<embed src="figure3_outputs/S16f.skcm_nmf_cluster_bcb.pdf" width="500px" height="500px" type="application/pdf" />

``` {.r}
makeKMplot_5group_withp(skcm_clin[!is.na(skcm_clin$surv),],"plot_cluster","figure3_outputs/S16g.skcm_subtype_survival.pdf")
```

    ## pdf 
    ##   2

``` {.r}
knitr::include_graphics("figure3_outputs/S16g.skcm_subtype_survival.pdf")
```

<embed src="figure3_outputs/S16g.skcm_subtype_survival.pdf" width="500px" height="500px" type="application/pdf" />

Figure S17

``` {.r}
plist_CPB_proj <- get.sample.association.heatmap(H.Bayes.new,g.Bayes.new,scale=1)
ggsave("figure3_outputs/S17b.CPB_H_norm_newcolors.pdf",plist_CPB_proj[[2]],height=7,width=7)
knitr::include_graphics("figure3_outputs/S17b.CPB_H_norm_newcolors.pdf")
```

<embed src="figure3_outputs/S17b.CPB_H_norm_newcolors.pdf" width="500px" height="500px" type="application/pdf" />

``` {.r}
mdf_cohort_cluster <- reshape2::melt(table(clin$cohort,clin$plot_cluster))
g_cohort_cluster <- ggplot(mdf_cohort_cluster,aes(x=Var1,y=value,fill=Var2)) + geom_col() + scale_fill_manual(values=cbpalette_reordered) + theme_classic() + xlab("Cohort") + ylab("Frequency")
ggsave("figure3_outputs/S17c.barplot_rna_cluster_by_cohort.pdf",g_cohort_cluster)
knitr::include_graphics("figure3_outputs/S17c.barplot_rna_cluster_by_cohort.pdf")
```

<embed src="figure3_outputs/S17c.barplot_rna_cluster_by_cohort.pdf" width="500px" height="500px" type="application/pdf" />

``` {.r}
kt <- kruskal.test(log10(clin[!is.na(clin$purity),]$nonsilent_snpdnpindel),clin[!is.na(clin$purity),]$plot_cluster)
kt_p <- formatSF(kt$p.value,format="e",digits=2)
g_cpb_nmf_tmb <- ggplot(clin[!is.na(clin$purity),],aes(x=plot_cluster,y=nonsilent_snpdnpindel)) + geom_boxplot(aes(fill=plot_cluster),outlier.shape = NA) + geom_jitter() +
scale_y_log10() + ylab("TMB") + theme_classic() + scale_fill_manual(values= cbpalette_reordered) + xlab("") + annotate(geom="text",x=3,y=2000,label=paste0("Kruskal-Wallis p =",kt_p))
ggsave("figure3_outputs/S17d.cpb_nmf_cluster_tmb.pdf",g_cpb_nmf_tmb)
knitr::include_graphics("figure3_outputs/S17d.cpb_nmf_cluster_tmb.pdf")
```

<embed src="figure3_outputs/S17d.cpb_nmf_cluster_tmb.pdf" width="500px" height="500px" type="application/pdf" />

``` {.r}
kt <- kruskal.test(clin[!is.na(clin$purity),]$purity,clin[!is.na(clin$purity),]$plot_cluster)
kt_p <- formatSF(kt$p.value,format="e",digits=2)
g_cpb_nmf_purity <- ggplot(clin[!is.na(clin$purity),],aes(x=plot_cluster,y=purity)) + geom_boxplot(aes(fill=plot_cluster),outlier.shape = NA) + geom_jitter() +
ylab("Tumor purity") + theme_classic() + scale_fill_manual(values= cbpalette_reordered) + xlab("") + annotate(geom="text",x=2.5,y=0.15,label=paste0("Kruskal-Wallis p =",kt_p))
ggsave("figure3_outputs/S17e.cpb_nmf_cluster_purity.pdf",g_cpb_nmf_purity)
knitr::include_graphics("figure3_outputs/S17e.cpb_nmf_cluster_purity.pdf")
```

<embed src="figure3_outputs/S17e.cpb_nmf_cluster_purity.pdf" width="500px" height="500px" type="application/pdf" />

``` {.r}
kt <- kruskal.test(log10(clin$rna_tcb),clin$plot_cluster)
kt_p <- formatSF(kt$p.value,format="e",digits=2)
g_cpb_nmf_tcb <- ggplot(clin,aes(x=plot_cluster,y=rna_tcb)) + geom_boxplot(aes(fill=plot_cluster),outlier.shape = NA) + geom_jitter() +
scale_y_log10() + ylab(expression(TCB[RNA])) + theme_classic() + scale_fill_manual(values= cbpalette_reordered) + xlab("") + annotate(geom="text",x=3,y=10,label=paste0("Kruskal-Wallis p =",kt_p))
ggsave("figure3_outputs/S17f.cpb_nmf_cluster_tcb.pdf",g_cpb_nmf_tcb)
knitr::include_graphics("figure3_outputs/S17f.cpb_nmf_cluster_tcb.pdf")
```

<embed src="figure3_outputs/S17f.cpb_nmf_cluster_tcb.pdf" width="500px" height="500px" type="application/pdf" />

``` {.r}
kt <- kruskal.test(log10(clin$rna_bcb),clin$plot_cluster)
kt_p <- formatSF(kt$p.value,format="e",digits=2)
g_cpb_nmf_bcb <- ggplot(clin,aes(x=plot_cluster,y=rna_bcb)) + geom_boxplot(aes(fill=plot_cluster),outlier.shape = NA) + geom_jitter() +
scale_y_log10() + ylab(expression(BCB[RNA])) + theme_classic() + scale_fill_manual(values= cbpalette_reordered) + xlab("") + annotate(geom="text",x=3,y=3000,label=paste0("Kruskal-Wallis p =",kt_p))
ggsave("figure3_outputs/S17g.cpb_nmf_cluster_bcb.pdf",g_cpb_nmf_bcb)
knitr::include_graphics("figure3_outputs/S17g.cpb_nmf_cluster_bcb.pdf")
```

<embed src="figure3_outputs/S17g.cpb_nmf_cluster_bcb.pdf" width="500px" height="500px" type="application/pdf" />

``` {.r}
tmp_ft <- chisq.test(table(clin$Response,clin$plot_cluster))
```

    ## Warning in chisq.test(table(clin$Response, clin$plot_cluster)): Chi-squared
    ## approximation may be incorrect

``` {.r}
tmp_ft_pstring <- paste0("Chi-Squared p=",formatSF(tmp_ft$p.value,digits=3))
tdf <- melt(table(clin$Response,clin$plot_cluster))
g_rnr_temp <- ggplot(tdf) + geom_mosaic(data=tdf,aes(weight=value,x=ggmosaic::product(Var2),fill=Var1),offset=0.01) + scale_fill_manual(values=c("red","blue")) + theme_bw() +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.title=element_blank(),axis.text.x=element_text(size=12,angle=90),axis.text.y=element_text(size=12)) +
xlab("") + ylab("") + annotate(geom="text",x=0.5,y=1.02,label=tmp_ft_pstring) + scale_y_continuous(breaks=seq(0,1,by=0.25)) 
ggsave("figure3_outputs/S17h.cpb_cluster_response.pdf",g_rnr_temp + geom_text(data = ggplot_build(g_rnr_temp)$data[[1]], aes(x = (xmin+xmax)/2, y = (ymin+ymax)/2, label=.wt)),height=7,width=7)
knitr::include_graphics("figure3_outputs/S17h.cpb_cluster_response.pdf")
```

<embed src="figure3_outputs/S17h.cpb_cluster_response.pdf" width="500px" height="500px" type="application/pdf" />

Figure S18

``` {.r}
col_count_data_all <- clin[,c("rna_id","cohort","Response")]

dds_all_fix <- DESeqDataSetFromTximport(txi= txi.rsem.coding, colData = col_count_data_all, design = ~ cohort + Response)
```

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

``` {.r}
dds_all_fix <- DESeq(dds_all_fix)
res_all_fix <- results(dds_all_fix, coef="Response_responder_vs_nonresponder")
df.res.all.fix <- data.frame(res_all_fix,stringsAsFactors=F)
df.res.all.fix$gene <- gencode_coding[rownames(df.res.all.fix),]$Gene

use_mediantpm1 <- names(apply(tab,1,median)[apply(tab,1,median)>1])

tmpdf_deseq_mediantpm1filtered_all_fix <- data.frame(res_all_fix[use_mediantpm1,],gene=gencode_coding[rownames(res_all_fix[use_mediantpm1,]),]$Gene,stringsAsFactors=F)
tmpdf_deseq_mediantpm1filtered_all_fix$new_padj <- p.adjust(tmpdf_deseq_mediantpm1filtered_all_fix$pvalue,method="fdr")


all_deseq_logpval_threshold <- 0.01
all_deseq_pval_threshold <- 2
all_deseq_lfc_threshold <- 0.75

g_all_volcano_base_fix <- ggplot(tmpdf_deseq_mediantpm1filtered_all_fix,aes(x=log2FoldChange,y=-log10(pvalue))) + geom_point(color="#999999") + theme_bw() + theme(legend.position = "none") + geom_vline(xintercept=c(-1* all_deseq_lfc_threshold, all_deseq_lfc_threshold), linetype="dotted") + geom_hline(yintercept=all_deseq_pval_threshold, linetype="dotted") + ylab("-log10(p value)") + xlab("log2(Fold Change)") + geom_point(data=tmpdf_deseq_mediantpm1filtered_all_fix[tmpdf_deseq_mediantpm1filtered_all_fix$pvalue<all_deseq_logpval_threshold&abs(tmpdf_deseq_mediantpm1filtered_all_fix$log2FoldChange)> all_deseq_lfc_threshold,],aes(x=log2FoldChange,y=-1*log10(pvalue),color=log2FoldChange>0)) + scale_color_manual(values=c("red","blue"))

g_all_volcano_fix <- g_all_volcano_base_fix + geom_text_repel(data=tmpdf_deseq_mediantpm1filtered_all_fix[tmpdf_deseq_mediantpm1filtered_all_fix$pvalue<all_deseq_logpval_threshold&abs(tmpdf_deseq_mediantpm1filtered_all_fix$log2FoldChange)> all_deseq_lfc_threshold,],aes(label=gene,x=log2FoldChange,y=-1*log10(pvalue),color=log2FoldChange>0),max.overlaps=100)
ggsave("figure3_outputs/S18a.deseq2_RNR_volcano_withlabels.pdf",g_all_volcano_fix,height=7,width=7)
knitr::include_graphics("figure3_outputs/S18a.deseq2_RNR_volcano_withlabels.pdf")
```

<embed src="figure3_outputs/S18a.deseq2_RNR_volcano_withlabels.pdf" width="500px" height="500px" type="application/pdf" />

``` {.r}
print(dim(tmpdf_deseq_mediantpm1filtered_all_fix[tmpdf_deseq_mediantpm1filtered_all_fix$new_padj<0.05,]))
```

    ## [1] 101   8

``` {.r}
print(table(tmpdf_deseq_mediantpm1filtered_all_fix[tmpdf_deseq_mediantpm1filtered_all_fix$new_padj<0.05,]$log2FoldChange<0))
```

    ## 
    ## FALSE  TRUE 
    ##    75    26

``` {.r}
sessionInfo()
```

    ## R version 3.4.0 (2017-04-21)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Red Hat Enterprise Linux Server release 6.9 (Santiago)
    ## 
    ## Matrix products: default
    ## BLAS: /broad/software/free/Linux/redhat_6_x86_64/pkgs/r_3.4.0/lib64/R/lib/libRblas.so
    ## LAPACK: /broad/software/free/Linux/redhat_6_x86_64/pkgs/r_3.4.0/lib64/R/lib/libRlapack.so
    ## 
    ## locale:
    ## [1] C
    ## 
    ## attached base packages:
    ##  [1] parallel  stats4    grid      methods   stats     graphics  grDevices
    ##  [8] utils     datasets  base     
    ## 
    ## other attached packages:
    ##  [1] DESeq2_1.16.1              SummarizedExperiment_1.6.1
    ##  [3] DelayedArray_0.2.7         matrixStats_0.53.1        
    ##  [5] Biobase_2.36.2             GenomicRanges_1.28.6      
    ##  [7] GenomeInfoDb_1.12.3        IRanges_2.10.5            
    ##  [9] S4Vectors_0.14.7           BiocGenerics_0.22.1       
    ## [11] tximportData_1.4.0         readr_1.3.1               
    ## [13] tximport_1.4.0             ggrepel_0.9.1             
    ## [15] RColorBrewer_1.1-2         gplots_3.0.1              
    ## [17] reshape2_1.4.3             ggmosaic_0.2.0            
    ## [19] gridExtra_2.3              psych_1.7.5               
    ## [21] rms_5.1-1                  SparseM_1.77              
    ## [23] Hmisc_4.0-3                Formula_1.2-3             
    ## [25] lattice_0.20-35            survminer_0.4.3           
    ## [27] ggpubr_0.2                 magrittr_2.0.1            
    ## [29] survival_2.42-6            ggplot2_3.2.1             
    ## [31] plyr_1.8.6                
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] TH.data_1.0-9           colorspace_1.4-1       
    ##  [3] ellipsis_0.3.2          rprojroot_1.3-2        
    ##  [5] htmlTable_1.9           XVector_0.16.0         
    ##  [7] base64enc_0.1-3         MatrixModels_0.4-1     
    ##  [9] bit64_0.9-7             AnnotationDbi_1.38.2   
    ## [11] fansi_0.5.0             mvtnorm_1.0-8          
    ## [13] codetools_0.2-15        splines_3.4.0          
    ## [15] mnormt_1.5-5            geneplotter_1.54.0     
    ## [17] knitr_1.16              jsonlite_1.7.2         
    ## [19] broom_0.5.0             km.ci_0.5-2            
    ## [21] annotate_1.54.0         cluster_2.0.7-1        
    ## [23] compiler_3.4.0          httr_1.4.2             
    ## [25] productplots_0.1.1      backports_1.1.2        
    ## [27] assertthat_0.2.0        Matrix_1.2-14          
    ## [29] fastmap_1.1.0           lazyeval_0.2.2         
    ## [31] acepack_1.4.1           htmltools_0.5.2        
    ## [33] quantreg_5.36           tools_3.4.0            
    ## [35] gtable_0.3.0            glue_1.4.2             
    ## [37] GenomeInfoDbData_0.99.0 dplyr_0.8.4            
    ## [39] Rcpp_1.0.7              vctrs_0.3.8            
    ## [41] gdata_2.18.0            nlme_3.1-137           
    ## [43] stringr_1.3.1           lifecycle_1.0.0        
    ## [45] gtools_3.8.1            XML_3.98-1.12          
    ## [47] polspline_1.1.13        MASS_7.3-50            
    ## [49] zlibbioc_1.22.0         zoo_1.8-3              
    ## [51] scales_0.5.0            hms_1.1.0              
    ## [53] sandwich_2.3-4          yaml_2.2.1             
    ## [55] memoise_1.1.0           KMsurv_0.1-5           
    ## [57] rpart_4.1-13            RSQLite_2.1.1          
    ## [59] latticeExtra_0.6-28     stringi_1.2.3          
    ## [61] genefilter_1.58.1       checkmate_1.8.5        
    ## [63] caTools_1.17.1          BiocParallel_1.10.1    
    ## [65] rlang_0.4.11            pkgconfig_2.0.3        
    ## [67] bitops_1.0-6            evaluate_0.11          
    ## [69] purrr_0.3.3             htmlwidgets_1.5.4      
    ## [71] labeling_0.3            cmprsk_2.2-7           
    ## [73] bit_1.1-14              tidyselect_1.0.0       
    ## [75] R6_2.5.1                multcomp_1.4-6         
    ## [77] DBI_1.0.0               pillar_1.6.2           
    ## [79] foreign_0.8-70          withr_2.4.2            
    ## [81] RCurl_1.95-4.11         nnet_7.3-12            
    ## [83] tibble_2.1.3            crayon_1.4.1           
    ## [85] survMisc_0.5.4          KernSmooth_2.23-15     
    ## [87] utf8_1.2.2              plotly_4.9.4.1         
    ## [89] rmarkdown_1.5           locfit_1.5-9.1         
    ## [91] data.table_1.14.0       blob_1.1.1             
    ## [93] digest_0.6.27           xtable_1.8-2           
    ## [95] tidyr_1.0.2             munsell_0.5.0          
    ## [97] viridisLite_0.3.0
