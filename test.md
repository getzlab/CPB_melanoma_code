test
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
library(ggpubr)
library(knitr)
library(ggsignif)
```

Load R objects

``` {.r}
clin <- readRDS("figure12_inputs/rna_clin.rds")
clin$surv <- Surv(clin$overall_survival,clin$DOD=="DOD")

dna_clin <- readRDS("figure12_inputs/dna_clin.rds")
dna_clin$surv <- Surv(dna_clin$overall_survival,dna_clin$DOD=="DOD")

skcm_clin <- readRDS("figure12_inputs/tcga_clin.rds")
skcm_clin$surv <- Surv(skcm_clin$overall_survival,skcm_clin$os_status)
skcm_clin_34 <- skcm_clin[skcm_clin$stage34==TRUE,]
```

Functions

``` {.r}
formatSF <- function(n, digits, format) {
         if(missing(format)) {
                        return(formatC(signif(n,digits=digits), digits=digits,format="fg", flag="#"))
         } else {
                        #return(formatC(signif(n,digits=digits), digits=digits,format=format, flag="#"))
                        return(formatC(signif(n,digits=digits), digits=digits-1,format=format, flag="#"))
         }
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

dna_clin$tmb_binary <- ifelse(dna_clin$nonsilent_snpdnpindel>median(dna_clin$nonsilent_snpdnpindel),"high_tmb","low_tmb")
dna_clin$plot_tmb_binary <- ifelse(dna_clin$nonsilent_snpdnpindel>median(dna_clin$nonsilent_snpdnpindel),"High TMB","Low TMB")

dir.create("tmp", showWarnings = FALSE)
#makeKMplot_2group(dna_clin,"plot_tmb_binary","figure12_outputs/1b.1.dna_tmb_survival.pdf")
knitr::include_graphics("figure12_outputs/1b.1.dna_tmb_survival.pdf")
```

<embed src="figure12_outputs/1b.1.dna_tmb_survival.pdf" width="500px" height="500px" type="application/pdf" />

![](figure12_outputs/1b.1.dna_tmb_survival.pdf)

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
    ## [1] grid      methods   stats     graphics  grDevices utils     datasets 
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] ggsignif_0.2.0  knitr_1.16      reshape2_1.4.3  ggmosaic_0.2.0 
    ##  [5] gridExtra_2.3   psych_1.7.5     rms_5.1-1       SparseM_1.77   
    ##  [9] Hmisc_4.0-3     Formula_1.2-3   lattice_0.20-35 survminer_0.4.3
    ## [13] ggpubr_0.2      magrittr_2.0.1  survival_2.42-6 ggplot2_3.2.1  
    ## [17] plyr_1.8.6     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] httr_1.4.2          tidyr_1.0.2         viridisLite_0.3.0  
    ##  [4] jsonlite_1.7.2      splines_3.4.0       assertthat_0.2.0   
    ##  [7] latticeExtra_0.6-28 yaml_2.2.1          pillar_1.6.2       
    ## [10] backports_1.1.2     quantreg_5.36       glue_1.4.2         
    ## [13] digest_0.6.27       RColorBrewer_1.1-2  checkmate_1.8.5    
    ## [16] colorspace_1.4-1    sandwich_2.3-4      productplots_0.1.1 
    ## [19] cmprsk_2.2-7        htmltools_0.5.2     Matrix_1.2-14      
    ## [22] pkgconfig_2.0.3     broom_0.5.0         purrr_0.3.3        
    ## [25] xtable_1.8-2        mvtnorm_1.0-8       scales_0.5.0       
    ## [28] km.ci_0.5-2         KMsurv_0.1-5        MatrixModels_0.4-1 
    ## [31] htmlTable_1.9       tibble_2.1.3        ellipsis_0.3.2     
    ## [34] TH.data_1.0-9       withr_2.4.2         nnet_7.3-12        
    ## [37] lazyeval_0.2.2      mnormt_1.5-5        crayon_1.4.1       
    ## [40] polspline_1.1.13    evaluate_0.11       fansi_0.5.0        
    ## [43] nlme_3.1-137        MASS_7.3-50         foreign_0.8-70     
    ## [46] tools_3.4.0         data.table_1.14.0   lifecycle_1.0.0    
    ## [49] multcomp_1.4-6      stringr_1.3.1       plotly_4.9.4.1     
    ## [52] munsell_0.5.0       cluster_2.0.7-1     compiler_3.4.0     
    ## [55] rlang_0.4.11        htmlwidgets_1.5.4   base64enc_0.1-3    
    ## [58] rmarkdown_1.5       gtable_0.3.0        codetools_0.2-15   
    ## [61] R6_2.5.1            zoo_1.8-3           dplyr_0.8.4        
    ## [64] fastmap_1.1.0       survMisc_0.5.4      utf8_1.2.2         
    ## [67] rprojroot_1.3-2     stringi_1.2.3       parallel_3.4.0     
    ## [70] Rcpp_1.0.7          vctrs_0.3.8         rpart_4.1-13       
    ## [73] acepack_1.4.1       tidyselect_1.0.0
