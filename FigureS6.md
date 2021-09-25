FigureS6
================
Sam Freeman

``` {.r}
library(ggplot2)
library(plyr)
library(reshape2)
library(data.table)
library(xtable)
library(survival)
```

Load R objects

``` {.r}
dna_maf <- read.delim("figure12_inputs/meta2_cpb_singlesample_189.validated_mc3filtered_dnpfix.aggregated.maf.lite.gz",header=T,as.is=T)
print("read mafs")
```

    ## [1] "read mafs"

``` {.r}
coding_categ <- c(("De_novo_Start_InFrame"),("De_novo_Start_OutOfFrame"),("Frame_Shift_Del"),("Frame_Shift_Ins"),("In_Frame_Del"),("In_Frame_Ins"),("Missense_Mutation"),
("Nonsense_Mutation"),("Nonstop_Mutation"),("Silent"),("Splice_Site"),("Start_Codon_Del"),("Start_Codon_Ins"),("Start_Codon_SNP"),("Stop_Codon_Del"),("Stop_Codon_Ins"))
cf_maf_coding <- dna_maf[dna_maf$Variant_Classification %in% coding_categ,]
cf_maf_nonsilent <- cf_maf_coding[cf_maf_coding$Variant_Classification!="Silent",]
cf_maf_nonsilent$is_lof <- ifelse(cf_maf_nonsilent$Variant_Classification %in% c("Missense_Mutation","In_Frame_Del","In_Frame_Ins"),FALSE,TRUE)

cf_maf_coding$key <- paste0(cf_maf_coding$Chromosome,cf_maf_coding$Start_Position,cf_maf_coding$Tumor_Seq_Allele2,"_",cf_maf_coding$Tumor_Sample_Barcode)
tmp_cf_maf_coding <- data.table(cf_maf_coding)
print("making key")
```

    ## [1] "making key"

``` {.r}
setkey(tmp_cf_maf_coding,"key")

clin <- readRDS("figure12_inputs/dna_clin.rds")
rownames(clin) <- clin$pair_id
clin$surv <- Surv(clin$overall_survival,clin$DOD=="DOD")
```

Functions

``` {.r}
# qq_pval
# By: Marcin Imielinski
# plots qq plot for observed pval vs uniform
# https://github.com/mskilab/skitools
qq_pval = function(obs, highlight = c(), hexbins = NULL, samp = NULL, lwd = 1, bestfit=T, color='black', input.pch=18, input.cex=1, conf.lines=T, input.MAX=NULL, qvalues=NULL, genes=NULL, ...)
{   
    obs = -log10(obs[!is.na(obs)])
    obs = obs[!is.infinite(obs)]
    
    if (!is.null(samp))
        if (samp<length(obs))
            obs = sample(obs, samp)
    
    N <- length(obs)
    ## create the null distribution 
    ## (-log10 of the uniform)
    exp <- -log(1:N/N,10)
    
    if (is.null(input.MAX))
        MAX <- max(obs,exp) + 0.5
    else
        MAX <- input.MAX
    
    c95 <- rep(0,N)
    c05 <- rep(0,N)
    
    for(i in 1:N){
        c95[i] <- qbeta(0.975,i,N-i+1)
        c05[i] <- qbeta(0.025,i,N-i+1)
    }
    
    if (conf.lines){
        ## plot the two confidence lines
        plot(exp, -log(c95,10), ylim=c(0,MAX), xlim=c(0,MAX), type="l", axes=FALSE, xlab="", ylab="")
        par(new=T)
        plot(exp, -log(c05,10), ylim=c(0,MAX), xlim=c(0,MAX), type="l", axes=FALSE, xlab="", ylab="")
        par(new=T)
        p1 <- c(exp[1], exp[1])
        p2 <- c(-log(c95,10)[1], -log(c05,10)[1])
        lines(x=p1, y=p2)
        x.coords <- c(exp,rev(exp))
        y.coords <- c(-log(c95,10),rev(-log(c05,10)))
        polygon(x.coords, y.coords, col='light gray', border=NA)
        par(new=T)
    }
    
    ord = order(obs)
    
    colors = vector(mode = "character", length = length(obs)); colors[] = "black"; colors[highlight] = "red";
    colors = factor(colors[ord]);
    
        dat = data.frame(x = sort(exp), y = obs[ord]);
        plot(dat$x, dat$y, xlab = expression(Expected -log[10](italic(P))), ylab = expression(Observed -log[10](italic(P))), xlim = c(0, MAX), ylim = c(0, MAX), pch=input.pch, cex=input.cex, bg=color, ...);
        if (!is.null(qvalues)){
          genes.to.label <- which(qvalues[ord] <= 0.2);
          genes <- genes[ord];
          if (length(genes.to.label > 0)){
            text(dat$x[genes.to.label], dat$y[genes.to.label], labels=genes[genes.to.label], pos=3);
          }
        }
        lines(x=c(0, MAX), y = c(0, MAX), col = "black", lwd = lwd);
        lambda = lm(y ~ x-1, dat)$coefficients;
        if (bestfit)
          {
            lines(x=c(0, MAX), y = c(0, lambda*MAX), col = "red", lty = 2, lwd = lwd);
            legend('bottomright',sprintf('lambda = %.2f', lambda), text.col='red', bty='n')
        }        
}
```

Figure S6

``` {.r}
snp <- dna_maf[dna_maf$Variant_Type %in% c("SNP","DNP","TNP","ONP"),]
indel <- dna_maf[dna_maf$Variant_Type %in% c("INS","DEL"),]
snp$Tumor_Sample_Barcode <- factor(snp$Tumor_Sample_Barcode,levels=clin$tumor_id)
indel$Tumor_Sample_Barcode <- factor(indel$Tumor_Sample_Barcode,levels=clin$tumor_id)
snp_coding <- snp[snp$Variant_Classification %in% coding_categ & snp$Variant_Classification!="Silent",]
indel_coding <- indel[indel$Variant_Classification %in% coding_categ,]
df_variant_tables <- data.frame(ns_s_SNV=table(snp$Tumor_Sample_Barcode),ns_SNV=table(snp_coding$Tumor_Sample_Barcode),ns_s_SSNV=table(snp[snp$Variant_Type=="SNP",]$Tumor_Sample_Barcode),ns_SSNV=table(snp_coding[snp_coding$Variant_Type=="SNP",]$Tumor_Sample_Barcode),ns_s_DNV=table(snp[snp$Variant_Type=="DNP",]$Tumor_Sample_Barcode),ns_s_DNV=table(snp_coding[snp_coding$Variant_Type=="DNP",]$Tumor_Sample_Barcode),ns_s_INDEL=table(indel$Tumor_Sample_Barcode),ns_INDEL=table(indel_coding$Tumor_Sample_Barcode),stringsAsFactors=F)
df_variant_tables <- df_variant_tables[,c(1,seq(2,16,by=2))]
colnames(df_variant_tables) <- unlist(lapply(strsplit(colnames(df_variant_tables),"\\."),"[",1))
colnames(df_variant_tables)[1] <- "tumor_id"
df_variant_tables$ns_s_SNV_INDEL <- df_variant_tables$ns_s_SNV + df_variant_tables$ns_s_INDEL
df_variant_tables$ns_SNV_INDEL <- df_variant_tables$ns_SNV + df_variant_tables$ns_INDEL

clin <- merge(clin,df_variant_tables,by="tumor_id")
rownames(clin) <- clin$tumor_id
clin$mutrate <- log10(clin$nonsilent_snpdnpindel)

print("finished mutrate")
```

    ## [1] "finished mutrate"

``` {.r}
clin$Response_binary <- as.numeric(factor(clin$Response,levels=c("responder","nonresponder")))*1-1

print("LoF fisher")
```

    ## [1] "LoF fisher"

``` {.r}
cf_maf_lof <- cf_maf_nonsilent[cf_maf_nonsilent$is_lof==TRUE,]
cf_maf_lof$key3 <- paste(cf_maf_lof$Hugo_Symbol,cf_maf_lof$Tumor_Sample_Barcode,sep="$")
u_cf_maf_lof <- unique(cf_maf_lof$key3)
df_u_lof <- data.frame(cbind(unlist(lapply(strsplit(u_cf_maf_lof,"\\$"),"[",1)),unlist(lapply(strsplit(u_cf_maf_lof,"\\$"),"[",2))),stringsAsFactors=F)
colnames(df_u_lof) <- c("gene","tumor_id")
over3_lof_genes <- names(table(df_u_lof$gene)[table(df_u_lof$gene)>=3])
df_u_lof <- df_u_lof[df_u_lof$gene %in% over3_lof_genes,]
df_u_lof$RNR <- clin[df_u_lof$tumor_id,]$Response
df_lof_res <- data.frame(cbind(gene=over3_lof_genes,n_NR=-1,n_R=-1,fisher_p=-1,fisher_or=-1),stringsAsFactors=F)
tot_n_NR <- as.numeric(table(clin$Response)["nonresponder"])
tot_n_R <- as.numeric(table(clin$Response)["responder"])

for (i in 1:nrow(df_lof_res)) {
    tmp_g <- df_lof_res$gene[i]
    n_R_mut <- nrow(df_u_lof[df_u_lof$gene==tmp_g & df_u_lof$RNR=="responder",])
    n_NR_mut <- nrow(df_u_lof[df_u_lof$gene==tmp_g & df_u_lof$RNR=="nonresponder",])
    tmp_ft <- fisher.test(matrix(c(n_R_mut,n_NR_mut,tot_n_R-n_R_mut,tot_n_NR-n_NR_mut),nrow=2))
    df_lof_res[i,]$n_NR <- n_NR_mut
    df_lof_res[i,]$n_R <- n_R_mut
    df_lof_res[i,]$fisher_p <- tmp_ft$p.value
    df_lof_res[i,]$fisher_or <- tmp_ft$estimate
}
df_lof_res$fisher_q <- p.adjust(df_lof_res$fisher_p,method="BH")
df_lof_res$frac_NR <- as.numeric(df_lof_res$n_NR)/tot_n_NR
df_lof_res$frac_R <- as.numeric(df_lof_res$n_R)/tot_n_R

cf_maf_nonsilent$key3 <- paste(cf_maf_nonsilent$Hugo_Symbol,cf_maf_nonsilent$Tumor_Sample_Barcode,sep="$")
u_cf_maf_nonsilent <- unique(cf_maf_nonsilent$key3)
df_u_ns <- data.frame(cbind(unlist(lapply(strsplit(u_cf_maf_nonsilent,"\\$"),"[",1)),unlist(lapply(strsplit(u_cf_maf_nonsilent,"\\$"),"[",2))),stringsAsFactors=F)
colnames(df_u_ns) <- c("gene","tumor_id")
over3_ns_genes <- names(table(df_u_ns$gene)[table(df_u_ns$gene)>=3])
df_u_ns <- df_u_ns[df_u_ns$gene %in% over3_ns_genes,]
df_u_ns$RNR <- clin[df_u_ns$tumor_id,]$Response
df_ns_res <- data.frame(cbind(gene=over3_ns_genes,n_NR=-1,n_R=-1,fisher_p=-1,fisher_or=-1),stringsAsFactors=F)

print("nsSNV fisher")
```

    ## [1] "nsSNV fisher"

``` {.r}
for (i in 1:nrow(df_ns_res)) {
    if (i %% 500 == 0) {print(i)}
    tmp_g <- df_ns_res$gene[i]
    n_R_mut <- nrow(df_u_ns[df_u_ns$gene==tmp_g & df_u_ns$RNR=="responder",])
    n_NR_mut <- nrow(df_u_ns[df_u_ns$gene==tmp_g & df_u_ns$RNR=="nonresponder",])
    tmp_ft <- fisher.test(matrix(c(n_R_mut,n_NR_mut,tot_n_R-n_R_mut,tot_n_NR-n_NR_mut),nrow=2))
    df_ns_res[i,]$n_NR <- n_NR_mut
    df_ns_res[i,]$n_R <- n_R_mut
    df_ns_res[i,]$fisher_p <- tmp_ft$p.value
    df_ns_res[i,]$fisher_or <- tmp_ft$estimate
}
```

    ## [1] 500
    ## [1] 1000
    ## [1] 1500
    ## [1] 2000
    ## [1] 2500
    ## [1] 3000
    ## [1] 3500
    ## [1] 4000
    ## [1] 4500
    ## [1] 5000
    ## [1] 5500
    ## [1] 6000
    ## [1] 6500
    ## [1] 7000
    ## [1] 7500
    ## [1] 8000
    ## [1] 8500
    ## [1] 9000
    ## [1] 9500
    ## [1] 10000
    ## [1] 10500

``` {.r}
df_ns_res$fisher_q <- p.adjust(df_ns_res$fisher_p,method="BH")
df_ns_res$frac_NR <- as.numeric(df_ns_res$n_NR)/tot_n_NR
df_ns_res$frac_R <- as.numeric(df_ns_res$n_R)/tot_n_R

pdf("figure12_outputs/S6b.RNR_lof_qq.pdf")
qq_pval(as.numeric(df_lof_res$fisher_p))
dev.off()
```

    ## pdf 
    ##   2

``` {.r}
knitr::include_graphics("figure12_outputs/S6b.RNR_lof_qq.pdf")
```

<embed src="figure12_outputs/S6b.RNR_lof_qq.pdf" width="500px" height="500px" type="application/pdf" />

``` {.r}
pdf("figure12_outputs/S6d.RNR_ns_qq.pdf")
qq_pval(as.numeric(df_ns_res$fisher_p))
dev.off()
```

    ## pdf 
    ##   2

``` {.r}
knitr::include_graphics("figure12_outputs/S6d.RNR_ns_qq.pdf")
```

<embed src="figure12_outputs/S6d.RNR_ns_qq.pdf" width="500px" height="500px" type="application/pdf" />

``` {.r}
df_lof_res_cox <- data.frame(cbind(gene=over3_lof_genes,n_NR=-1,n_R=-1,cox_p=-1,cox_hr=-1,notmb_cox_p=-1,notmb_cox_hr=-1,logistic_p=-1,logistic_or=-1,notmb_logistic_p=-1,notmb_logistic_or=-1),stringsAsFactors=F)

print("lof cox")
```

    ## [1] "lof cox"

``` {.r}
clin$mutstatus <- 0
for (i in 1:nrow(df_lof_res_cox)) {
    if (i %% 100 == 0) {print(i)}
    clin$mutstatus <- 0
    tmp_g <- df_lof_res_cox$gene[i]
    n_R_mut <- nrow(df_u_lof[df_u_lof$gene==tmp_g & df_u_lof$RNR=="responder",])
    n_NR_mut <- nrow(df_u_lof[df_u_lof$gene==tmp_g & df_u_lof$RNR=="nonresponder",])
    df_lof_res_cox[i,]$n_NR <- n_NR_mut
    df_lof_res_cox[i,]$n_R <- n_R_mut
    clin[df_u_lof[df_u_lof$gene==tmp_g,]$tumor_id,]$mutstatus <- 1
    tmp_cox <- summary(coxph(clin$surv ~ clin$mutrate + clin$mutstatus))
    df_lof_res_cox[i,]$cox_p <- tmp_cox$coef[2,5]
    df_lof_res_cox[i,]$cox_hr <- tmp_cox$coef[2,2]
    tmp_logistic <- summary(glm(clin$Response_binary ~ clin$mutrate + clin$mutstatus,family=binomial(link="logit")))
    df_lof_res_cox[i,]$logistic_p <- tmp_logistic$coef[3,4]
    df_lof_res_cox[i,]$logistic_or <- exp(tmp_logistic$coef[3,1])
}
```

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## [1] 100

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## [1] 200

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## [1] 300

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## [1] 400

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## [1] 500

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## [1] 600

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## [1] 700

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## [1] 800

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## [1] 900

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

``` {.r}
df_lof_res_cox$cox_p <- as.numeric(df_lof_res_cox$cox_p)
df_lof_res_cox$logistic_p <- as.numeric(df_lof_res_cox$logistic_p)

df_lof_res_cox$cox_q <- p.adjust(df_lof_res_cox$cox_p,method="BH")
df_lof_res_cox$logistic_q <- p.adjust(df_lof_res_cox$logistic_p,method="BH")

pdf("figure12_outputs/S6a.OS_lof_qq.pdf")
qq_pval(as.numeric(df_lof_res_cox$cox_p))
dev.off()
```

    ## pdf 
    ##   2

``` {.r}
knitr::include_graphics("figure12_outputs/S6a.OS_lof_qq.pdf")
```

<embed src="figure12_outputs/S6a.OS_lof_qq.pdf" width="500px" height="500px" type="application/pdf" />

``` {.r}
write.table(df_lof_res_cox[,c("gene","cox_hr","cox_p","cox_q","logistic_or","logistic_p","logistic_q")],gzfile("figure12_outputs/lof_table.txt.gz"),quote=F,sep="\t")

print("nsSNV cox")
```

    ## [1] "nsSNV cox"

``` {.r}
df_ns_res_cox <- data.frame(cbind(gene=over3_ns_genes,n_NR=-1,n_R=-1,cox_p=-1,cox_hr=-1,notmb_cox_p=-1,notmb_cox_hr=-1,logistic_p=-1,logistic_or=-1,notmb_logistic_p=-1,notmb_logistic_or=-1),stringsAsFactors=F)

for (i in 1:nrow(df_ns_res)) {
    if (i %% 1000 == 0) {print(i)}
    clin$mutstatus <- 0
    tmp_g <- df_ns_res_cox$gene[i]
    n_R_mut <- nrow(df_u_ns[df_u_ns$gene==tmp_g & df_u_ns$RNR=="responder",])
    n_NR_mut <- nrow(df_u_ns[df_u_ns$gene==tmp_g & df_u_ns$RNR=="nonresponder",])
    df_ns_res_cox[i,]$n_NR <- n_NR_mut
    df_ns_res_cox[i,]$n_R <- n_R_mut
    clin[df_u_ns[df_u_ns$gene==tmp_g,]$tumor_id,]$mutstatus <- 1
    tmp_cox <- summary(coxph(clin$surv ~ clin$mutrate + clin$mutstatus))
    df_ns_res_cox[i,]$cox_p <- tmp_cox$coef[2,5]
    df_ns_res_cox[i,]$cox_hr <- tmp_cox$coef[2,2]
    tmp_logistic <- summary(glm(clin$Response_binary ~ clin$mutrate + clin$mutstatus,family=binomial(link="logit")))
    df_ns_res_cox[i,]$logistic_p <- tmp_logistic$coef[3,4]
    df_ns_res_cox[i,]$logistic_or <- exp(tmp_logistic$coef[3,1])
}
```

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## [1] 1000

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## [1] 2000

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## [1] 3000

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## [1] 4000

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## [1] 5000

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## [1] 6000

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## [1] 7000

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## [1] 8000

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## [1] 9000

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## [1] 10000

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 2 ; beta may be infinite.

``` {.r}
df_ns_res_cox$cox_p <- as.numeric(df_ns_res_cox$cox_p)
df_ns_res_cox$logistic_p <- as.numeric(df_ns_res_cox$logistic_p)

df_ns_res_cox$cox_q <- p.adjust(df_ns_res_cox$cox_p,method="BH")
df_ns_res_cox$logistic_q <- p.adjust(df_ns_res_cox$logistic_p,method="BH")


pdf("figure12_outputs/S6c.OS_ns_qq.pdf")
qq_pval(as.numeric(df_ns_res_cox$cox_p))
dev.off()
```

    ## pdf 
    ##   2

``` {.r}
knitr::include_graphics("figure12_outputs/S6c.OS_ns_qq.pdf")
```

<embed src="figure12_outputs/S6c.OS_ns_qq.pdf" width="500px" height="500px" type="application/pdf" />

``` {.r}
write.table(df_ns_res_cox[,c("gene","cox_hr","cox_p","cox_q","logistic_or","logistic_p","logistic_q")],gzfile("figure12_outputs/ns_table.txt.gz"),quote=F,sep="\t")
```

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
    ## [1] stats     graphics  grDevices utils     datasets  base     
    ## 
    ## other attached packages:
    ## [1] survival_2.42-6   xtable_1.8-2      data.table_1.14.0 reshape2_1.4.3   
    ## [5] plyr_1.8.6        ggplot2_3.2.1    
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.7       pillar_1.6.2     compiler_3.4.0   methods_3.4.0   
    ##  [5] tools_3.4.0      digest_0.6.27    lattice_0.20-35  evaluate_0.11   
    ##  [9] lifecycle_1.0.0  tibble_2.1.3     gtable_0.3.0     pkgconfig_2.0.3 
    ## [13] rlang_0.4.11     Matrix_1.2-14    yaml_2.2.1       fastmap_1.1.0   
    ## [17] withr_2.4.2      stringr_1.3.1    dplyr_0.8.4      knitr_1.16      
    ## [21] vctrs_0.3.8      rprojroot_1.3-2  grid_3.4.0       tidyselect_1.0.0
    ## [25] glue_1.4.2       R6_2.5.1         fansi_0.5.0      rmarkdown_1.5   
    ## [29] purrr_0.3.3      magrittr_2.0.1   splines_3.4.0    backports_1.1.2 
    ## [33] scales_0.5.0     ellipsis_0.3.2   htmltools_0.5.2  assertthat_0.2.0
    ## [37] colorspace_1.4-1 utf8_1.2.2       stringi_1.2.3    lazyeval_0.2.2  
    ## [41] munsell_0.5.0    crayon_1.4.1
