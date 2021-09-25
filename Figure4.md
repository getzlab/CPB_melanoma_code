Figure4
================
Sam Freeman

``` {.r}
library(ggplot2)
library(survival)
library(survminer)
library(ggmosaic)
#note: tricky dependency issues with ggmosaic/ggplot2 versions
#see "Version compatibility issues with ggplot2" in https://cran.r-project.org/web/packages/ggmosaic/vignettes/ggmosaic.html
library(reshape2)
library(corrplot)
library(ggsignif)
library(meta)
library(plyr)
library(ggrepel)
library(pheatmap)
library(verification)
library(pROC)
library(RColorBrewer)
library(psych)
library(data.table)
library(grid)
library(gridExtra)
library(RAM)
```

Load R objects

``` {.r}
fcmeta_tab <- readRDS("figure4_inputs/small_fcmeta_tab.rds")

clin <- readRDS("figure4_inputs/primary_clin.rds")
clin$surv <- Surv(clin$overall_survival,clin$DOD=="DOD")

val_df_res_sig <- readRDS("figure4_inputs/val_df_res_sig.rds")
val_df_res_os_sig <- readRDS("figure4_inputs/val_df_res_os_sig.rds")

clin[,val_df_res_sig$gene] <- t(fcmeta_tab[rownames(val_df_res_sig),])
clin[,val_df_res_os_sig$gene] <- t(fcmeta_tab[rownames(val_df_res_os_sig),])

pcaexl_ccle_melanoma_tab <- readRDS("figure4_inputs/small_ccle_tab.rds")
hpa_tab <- readRDS("figure4_inputs/small_hpa_tab.rds")

val_df_res_os_sig$CCLE_median_TPM <- as.numeric(apply(pcaexl_ccle_melanoma_tab[rownames(val_df_res_os_sig),],1,median))
val_df_res_os_sig$human_protein_atlas_median <- log2(as.numeric(apply(hpa_tab[unlist(lapply(strsplit(rownames(val_df_res_os_sig),"\\."),"[",1)),],1,median))+1)
```

Functions

``` {.r}
StatSignif <- ggplot2::ggproto("StatSignif", ggplot2::Stat,
                  required_aes = c("x", "y", "group"),
                  setup_params = function(data, params) {
                    # if(any(data$group == -1)|| any(data$group != data$x)){
                    if(any(data$group == -1)){
                      stop("Can only handle data with groups that are plotted on the x-axis")
                    }

                    if (is.character(params$test)) params$test <- match.fun(params$test)
                    params$complete_data <- data
                    if(is.null(params$xmin) != is.null(params$xmax) || length(params$xmin) != length(params$xmax))
                      stop("If xmin or xmax is set, the other one also needs to be set and they need to contain the same number of values")
                    if(!is.null(params$xmin) && !is.null(params$comparisons))
                      stop("Set either the xmin, xmax values or the comparisons")
                    if(!is.null(params$xmin) && is.null(params$y_position))
                      stop("If xmin, xmax are defined also define y_position")
                    if(! is.null(params$y_position) && length(params$y_position) == 1)
                      params$y_position <- rep(params$y_position, max(length(params$comparisons), length(params$xmin), 1))
                    if(length(params$margin_top) == 1) params$margin_top <- rep(params$margin_top, max(length(params$comparisons),length(params$xmin), 1))
                    if(length(params$step_increase) == 1) params$step_increase <- rep(params$step_increase, max(length(params$comparisons),length(params$xmin), 1))
                    if(length(params$tip_length) == 1) params$tip_length <- rep(params$tip_length, max(length(params$comparisons),length(params$xmin), 1) * 2)
                    if(length(params$tip_length) == length(params$comparisons)) params$tip_length <- rep(params$tip_length, each=2)
                    if(length(params$tip_length) == length(params$xmin)) params$tip_length <- rep(params$tip_length, each=2)
                    if(! is.null(params$annotations) && length(params$annotations) == 1)
                      params$annotations <- rep(params$annotations, max(length(params$comparisons),length(params$xmin), 1))
                    if(! is.null(params$annotations) && length(params$annotations) != max(length(params$comparisons),length(params$xmin), 1))
                      stop(paste0("annotations contains a different number of elements (", length(params$annotations),
                                  ") than comparisons or xmin (", max(length(params$comparisons),length(params$xmin), 1), ")."))

                    if(all(is.logical(params$map_signif_level)) && all(params$map_signif_level == TRUE)){
                      params$map_signif_level <- c("***"=0.001, "**"=0.01, "*"=0.05)
                    }else if(is.numeric(params$map_signif_level)){
                      if(is.null(names(params$map_signif_level)) ){
                        if(length(params$map_signif_level) <= 3){
                          names(params$map_signif_level) <- tail(c("***", "**", "*"), n=length(params$map_signif_level))
                        }else{
                          stop('Cannot handle un-named map for significance values, please provide in the following format: c("***"=0.001, "**"=0.01, "*"=0.05)')
                        }
                      }
                    }
                    return(params)
                  },
                  compute_group = function(data, scales, comparisons, test, test.args, complete_data,
                                           annotations, map_signif_level, y_position, xmax, xmin,
                                           margin_top, step_increase, tip_length, manual) {

                    if("annotations" %in% colnames(data)){
                      annotations <- data[["annotations"]]
                    }
                    if("y_position" %in% colnames(data)){
                      y_position <- data[["y_position"]]
                    }
                    if("xmax" %in% colnames(data)){
                      xmax <- data[["xmax"]]
                    }
                    if("xmin" %in% colnames(data)){
                      xmin <- data[["xmin"]]
                    }
                    if("map_signif_level" %in% colnames(data)){
                      map_signif_level <- data[["map_signif_level"]]
                    }
                    if("tip_length" %in% colnames(data)){
                      tip_length <-  rep(data[["tip_length"]], each=2)
                    }

                    if(! is.null(comparisons)){
                      i <- 0
                      result <- lapply(comparisons, function(comp){
                        i <<- i + 1
                        # All entries in group should be the same
                        if(scales$x$map(comp[1]) == data$group[1] | manual){
                          test_result <- if(is.null(annotations)){
                            group_1 <- complete_data$y[complete_data$x == scales$x$map(comp[1]) & complete_data$PANEL == data$PANEL[1]]
                            group_2 <- complete_data$y[complete_data$x == scales$x$map(comp[2]) & complete_data$PANEL == data$PANEL[1]]
                            p_value <- do.call(test, c(list(group_1, group_2), test.args))$p.value
                            if(is.numeric(map_signif_level)){
                              temp_value <- names(which.min(map_signif_level[which(map_signif_level > p_value)]))
                              if(is.null(temp_value)){
                                "NS."
                              }else{
                                temp_value
                              }
                            }else if(is.function(map_signif_level)){
                              map_signif_level(p_value)
                            }else{
                              if(is.numeric(p_value)){
                                  as.character(sprintf("%.2g", p_value))
                              }else{
                                as.character(p_value)
                              }

                            }
                          }else{
                            annotations[i]
                          }
                          y_scale_range <- (scales$y$range$range[2] - scales$y$range$range[1])
                          if(is.null(y_position)){
                            y_pos <- scales$y$range$range[2] + y_scale_range * margin_top[i] + y_scale_range * step_increase[i] * (i-1)
                          }else{
                            y_pos <- y_position[i] + y_scale_range * margin_top[i] + y_scale_range * step_increase[i] * (i-1)
                          }
                          data.frame(x=c(min(comp[1],comp[2]),min(comp[1],comp[2]),max(comp[1],comp[2])),
                                     xend=c(min(comp[1],comp[2]),max(comp[1],comp[2]),max(comp[1],comp[2])),
                                     y=c(y_pos - y_scale_range*tip_length[(i-1)*2+1], y_pos, y_pos),
                                     yend=c(y_pos, y_pos, y_pos-y_scale_range*tip_length[(i-1)*2+2]),
                                     annotation=test_result, group=paste(c(comp, i), collapse = "-"))
                        }
                      })
                      do.call(rbind, result)
                    }else{
                      if((data$x[1] == min(complete_data$x) & data$group[1] == min(complete_data$group)) | manual) {
                        y_scale_range <- (scales$y$range$range[2] - scales$y$range$range[1])
                        if(is.character(xmin)){
                          xmin <- scales$x$map(xmin)
                        }
                        if(is.character(xmax)){
                          xmax <- scales$x$map(xmax)
                        }
                        if("expression" %in% class(annotations)){
                          stop("annotations must be a character vector. To use plotmath set parse=TRUE.")
                        }
                        data.frame(x=c(xmin, xmin, xmax),
                                   xend=c(xmin, xmax, xmax),
                                   y=c(y_position - y_scale_range*tip_length[seq_len(length(tip_length))%% 2 == 1], y_position, y_position),
                                   yend=c(y_position, y_position, y_position-y_scale_range*tip_length[seq_len(length(tip_length))%% 2 == 0]),
                                   annotation=rep(annotations, times=3), group=if(manual){rep(data$group, times=3)}else{rep(seq_along(xmin), times=3)})
                      }
                    }
                  }
)


#' Create significance layer
#'
#' @param comparisons A list of length-2 vectors.
#'   The entries in the vector are either the names of 2 values on the x-axis
#'   or the 2 integers that correspond to the index of the columns of interest
#' @param test the name of the statistical test that is applied to the values of the 2 columns (e.g. `t.test`, `wilcox.test` etc.).
#'   If you implement a custom test make sure that it returns a list that has an entry called `p.value`.
#' @param test.args additional arguments for the test method
#' @param annotations character vector with alternative annotations, if not null test is ignored
#' @param map_signif_level boolean value, if the p-value are directly written as annotation or asterisks are used instead.
#'   Alternatively one can provide a named numeric vector to create custom mappings from p-values to annotation:
#'   For example: c("***"=0.001, "**"=0.01, "*"=0.05)
#'   Alternatively, one can provide a function that takes a numeric argument (the p-value) and returns a string
#' @param xmin numeric vector with the positions of the left sides of the brackets
#' @param xmax numeric vector with the positions of the right sides of the brackets
#' @param y_position numeric vector with the y positions of the brackets
#' @param size change the width of the lines of the bracket
#' @param textsize change the size of the text
#' @param family change the font used for the text
#' @param vjust move the text up or down relative to the bracket
#' @param margin_top numeric vector how much higher that the maximum value that bars start as fraction of total height
#' @param step_increase numeric vector with the increase in fraction of total height for every additional comparison to
#'   minimize overlap.
#' @param tip_length numeric vector with the fraction of total height that the bar goes down to indicate the precise column
#' @param parse If `TRUE`, the labels will be parsed into expressions and
#'   displayed as described in `?plotmath`.
#' @param manual boolean flag that indicates that the parameters are provided with a data.frame. This option is necessary if
#'   one wants to plot different annotations per facet.
#' @param na.rm If \code{FALSE} (the default), removes missing values with
#'    a warning.  If \code{TRUE} silently removes missing values.
#' @param ... other arguments passed on to \code{\link{layer}}. These are
#'   often aesthetics, used to set an aesthetic to a fixed value, like
#'   \code{color = "red"} or \code{size = 3}. They may also be parameters
#'   to the paired geom/stat.
#' @inheritParams ggplot2::layer
#' @examples
#' \dontrun{
#' library(ggplot2)
#' library(ggsignif)
#' ggplot(mpg, aes(class, hwy)) +
#'  geom_boxplot() +
#'  geom_signif(comparisons = list(c("compact", "pickup"),
#'                                 c("subcompact", "suv")))
#' ggplot(mpg, aes(class, hwy)) +
#'  geom_boxplot() +
#'  geom_signif(comparisons = list(c("compact", "pickup"),
#'                                 c("subcompact", "suv")),
#'              map_signif_level=function(p)sprintf("p = %.2g", p))
#'
#' ggplot(mpg, aes(class, hwy)) +
#'   geom_boxplot() +
#'   geom_signif(annotations = c("First", "Second"),
#'               y_position = c(30, 40), xmin=c(4,1), xmax=c(5,3))
#' }
#'
#' @export
stat_signif <- function(mapping = NULL, data = NULL,
                    position = "identity", na.rm = FALSE, show.legend = NA,
                    inherit.aes = TRUE, comparisons=NULL, test="wilcox.test", test.args=NULL,
                    annotations=NULL, map_signif_level=FALSE,y_position=NULL,xmin=NULL, xmax=NULL,
                    margin_top=0.05, step_increase=0, tip_length=0.03,
                    size=0.5, textsize = 3.88, family="", vjust = 0,
                    parse = FALSE, manual=FALSE,
                    ...) {
  if(manual){
    if(! is.null(data) & ! is.null(mapping)){
      if(! "x" %in% names(data)) mapping$x <- 1
      if(! "y" %in% names(data)) mapping$y <- 1
    }else{
      stop("If manual mode is selected you need to provide the data and mapping parameters")
    }
  }
  ggplot2::layer(
    stat = StatSignif, data = data, mapping = mapping, geom = "signif",
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(comparisons=comparisons, test=test, test.args=test.args,
                  annotations=annotations, map_signif_level=map_signif_level,
                  y_position=y_position,xmin=xmin, xmax=xmax,
                  margin_top=margin_top, step_increase=step_increase,
                  tip_length=tip_length, size=size, textsize=textsize,
                  family=family, vjust=vjust, parse=parse, manual=manual, na.rm = na.rm, ...)
  )
}


GeomSignif <- ggplot2::ggproto("GeomSignif", ggplot2::Geom,
                           required_aes = c("x", "xend", "y", "yend", "annotation"),
                           default_aes = ggplot2::aes(shape = 19, colour = "black", textsize = 3.88, angle = 0, hjust = 0.5,
                                             vjust = 0, alpha = NA, family = "", fontface = 1, lineheight = 1.2, linetype=1, size=0.5),
                           draw_key = function(...){grid::nullGrob()},

                           draw_group = function(data, panel_params, coord, parse=FALSE) {
                             lab <- as.character(data$annotation)
                             if (parse) {
                               lab <- parse_safe(as.character(lab))
                             }
                             coords <- coord$transform(data, panel_params)
                             grid::gList(
                               grid::textGrob(
                                 label=lab,
                                 x=mean(c(coords$x[1], tail(coords$xend, n=1))),
                                 y=max(c(coords$y, coords$yend))+0.01,
                                 default.units = "native",
                                 hjust = coords$hjust, vjust = coords$vjust,
                                 rot = coords$angle,
                                 gp = grid::gpar(
                                   col = scales::alpha(coords$colour, coords$alpha),
                                   fontsize = coords$textsize * ggplot2::.pt,
                                   fontfamily = coords$family,
                                   fontface = coords$fontface,
                                   lineheight = coords$lineheight
                                 )
                               ),
                               grid::segmentsGrob(
                                 coords$x, coords$y,
                                 default.units = "native",
                                 coords$xend, coords$yend,
                                 gp = grid::gpar(
                                   col = scales::alpha(coords$colour, coords$alpha),
                                   lty = coords$linetype,
                                   lwd = coords$size * ggplot2::.pt
                                 )
                               )
                             )
                           }
)

#' @rdname stat_signif
#' @export
geom_signif <- function(mapping = NULL, data = NULL, stat = "signif",
                        position = "identity", na.rm = FALSE, show.legend = NA,
                        inherit.aes = TRUE, comparisons=NULL, test="wilcox.test", test.args=NULL,
                        annotations=NULL, map_signif_level=FALSE,y_position=NULL,xmin=NULL, xmax=NULL,
                        margin_top=0.05, step_increase=0, tip_length=0.03,
                        size=0.5, textsize = 3.88, family="", vjust = 0,
                        parse = FALSE,
                        manual=FALSE,
                        ...) {
  params <- list(na.rm = na.rm, ...)
  if (identical(stat, "signif")) {
    if(! is.null(data) & ! is.null(mapping) & ! manual) warning("You have set data and mapping, are you sure that manual = FALSE is correct?")
    if(manual){
      if(is.null(mapping$annotations)) stop("Manual mode only works if with 'annotations' is provided in mapping")
      if(! is.null(data) & ! is.null(mapping)){
        if(! "x" %in% names(mapping)){
          if("xmin" %in% names(mapping)){
            mapping$x <- mapping$xmin
          }else{
            mapping$x <- xmin
          }
        }
        if(! "y" %in% names(mapping)){
          if("y_position" %in% names(mapping)){
            mapping$y <- mapping$y_position
          }else{
            mapping$y <- y_position
          }
        }
      }else{
        stop("If manual mode is selected you need to provide the data and mapping parameters")
      }
    }
    params <- c(params, list(comparisons=comparisons, test=test, test.args=test.args,
                   annotations=annotations, map_signif_level=map_signif_level,
                   y_position=y_position,xmin=xmin, xmax=xmax,
                   margin_top=margin_top, step_increase=step_increase,
                   tip_length=tip_length, size=size, textsize=textsize,
                   family=family, vjust=vjust, parse=parse, manual=manual))
  }
  ggplot2::layer(
    stat = stat, geom = GeomSignif, mapping = mapping,  data = data,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = params
  )
}
```

Figure 3

``` {.r}
dir.create("figure3_outputs", showWarnings = FALSE)

g_val_os_CCLE_HPA <- ggplot(val_df_res_os_sig[!is.na(val_df_res_os_sig$human_protein_atlas_median),],aes(x=human_protein_atlas_median,y=CCLE_median_TPM)) +
geom_point(aes(color=log2FoldChange<0)) + geom_vline(xintercept=1) +  geom_hline(yintercept=1) + theme_classic(base_size=16) + scale_color_manual(values=c("#D55E00","#0072B2")) +
theme(legend.position="bottom") + geom_text_repel(aes(label=gene),size=3) + geom_abline(slope=1,intercept=0) + xlim(0,10) + ylim(0,7) +
xlab("HPA Blood median log2(TPM+1)") + ylab("CCLE Melanoma median log2(TPM+1)")
ggsave("figure3_outputs/3f.CCLE_HPA_expression_OSgenes.pdf",g_val_os_CCLE_HPA,height=7,width=7)
```

    ## Warning: ggrepel: 10 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

``` {.r}
knitr::include_graphics("figure3_outputs/3f.CCLE_HPA_expression_OSgenes.pdf")
```

<embed src="figure3_outputs/3f.CCLE_HPA_expression_OSgenes.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
tmp_cor_res_os <- cor(t(fcmeta_tab[rownames(val_df_res_os_sig),]))
rownames(tmp_cor_res_os) <- val_df_res_os_sig$gene
colnames(tmp_cor_res_os) <- val_df_res_os_sig$gene
m_tmp_cor_res_os <- reshape2::melt(tmp_cor_res_os)
m_tmp_cor_res_os$gene_pair <- paste0(m_tmp_cor_res_os$Var1,"_",m_tmp_cor_res_os$Var2)

pdf("figure3_outputs/S19a.val_gene_cor_os.pdf")
corrplot(tmp_cor_res_os, tl.cex=0.50,type = "upper", order = "hclust",tl.col = "black", tl.srt = 90,col=colorRampPalette(c("blue","white","red"))(200))
dev.off()
```

    ## pdf 
    ##   2

``` {.r}
knitr::include_graphics("figure3_outputs/S19a.val_gene_cor_os.pdf")
```

<embed src="figure3_outputs/S19a.val_gene_cor_os.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
dt_hpa_os <- data.table(t(log2(hpa_tab[unlist(lapply(strsplit(rownames(val_df_res_os_sig[val_df_res_os_sig$log2FoldChange<0,]),"\\."),"[",1)),]+1)))
colnames(dt_hpa_os) <- val_df_res_os_sig[val_df_res_os_sig$log2FoldChange<0,]$gene
dt_hpa_os$cell_type <- unlist(lapply(strsplit(colnames(hpa_tab),"\\.\\.\\."),"[",1))
tmp_hpa_os_rank <- data.frame(dt_hpa_os[ , lapply(.SD, mean) , by=cell_type],stringsAsFactors=F)
tmp_hpa_os_rank <- data.frame(apply(tmp_hpa_os_rank[,2:ncol(tmp_hpa_os_rank)],2,rank),row.names=tmp_hpa_os_rank[,1],stringsAsFactors=F)
m_hpa_os_rank <- reshape2::melt(as.matrix(tmp_hpa_os_rank))
m_hpa_os_rank$Var1 <- gsub("\\."," ",m_hpa_os_rank$Var1)
m_hpa_os_rank$Var1 <- factor(as.character(m_hpa_os_rank$Var1),levels=rev(as.character(arrange(data.table(m_hpa_os_rank)[,median(value),by=Var1],V1)$Var1)))

g_osgene_hpa_ranks <- ggplot(m_hpa_os_rank,aes(x=Var1,y=value)) + geom_boxplot(aes(fill=Var1),outlier.shape=NA) + geom_jitter(height=0) + theme_classic(base_size=16) +
theme(legend.position="none",axis.text.x=element_text(angle=90,hjust=1)) + xlab("") + ylab("HPA Blood cell type expression rank")
ggsave("figure3_outputs/S19c.OSgene_HPA_expression_rank.pdf",g_osgene_hpa_ranks,height=6,width=12)
knitr::include_graphics("figure3_outputs/S19c.OSgene_HPA_expression_rank.pdf")
```

<embed src="figure3_outputs/S19c.OSgene_HPA_expression_rank.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
dt_hpa_os_zscore <- data.table(t(log2(hpa_tab[unlist(lapply(strsplit(rownames(val_df_res_os_sig[val_df_res_os_sig$log2FoldChange<0,]),"\\."),"[",1)),]+1)))
colnames(dt_hpa_os_zscore) <- val_df_res_os_sig[val_df_res_os_sig$log2FoldChange<0,]$gene
dt_hpa_os_zscore <- data.table(apply(dt_hpa_os_zscore,2,scale))
dt_hpa_os_zscore$cell_type <- unlist(lapply(strsplit(colnames(hpa_tab),"\\.\\.\\."),"[",1))
df_hpa_os_zscore <- data.frame(dt_hpa_os_zscore[ , lapply(.SD, mean) , by=cell_type],stringsAsFactors=F)
rownames(df_hpa_os_zscore) <- gsub("\\."," ",df_hpa_os_zscore$cell_type)
df_hpa_os_zscore$cell_type <- NULL

dev.off()
```

    ## null device 
    ##           1

``` {.r}
pdf("figure3_outputs/S19b.OSgene_HPA_zscore_nobreaks.pdf",height=7,width=8)
pheatmap(df_hpa_os_zscore[,!(colnames(df_hpa_os_zscore) %in% c("C4A","CCL21"))],color=rev(brewer.pal(n=11, name="RdBu")),fontsize_col=8,fontsize_row=8,cellheight=8,cellwidth=8,treeheight_row=15,treeheight_col=15)
dev.off()
```

    ## null device 
    ##           1

``` {.r}
knitr::include_graphics("figure3_outputs/S19b.OSgene_HPA_zscore_nobreaks.pdf")
```

<embed src="figure3_outputs/S19b.OSgene_HPA_zscore_nobreaks.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
temp_colors <- c("#56B4E9","#009E73","#999999","#F0E442","#E69F00")
names(temp_colors) <- c("MITF low","Immune","Intermediate","Keratin high","MITF high")

dt_cpb_os <- data.table(t(fcmeta_tab[rownames(val_df_res_os_sig[val_df_res_os_sig$log2FoldChange>0,]),]))
colnames(dt_cpb_os) <- val_df_res_os_sig[val_df_res_os_sig$log2FoldChange>0,]$gene
dt_cpb_os$plot_cluster <- clin$plot_cluster
tmp_cpb_os_rank <- data.frame(dt_cpb_os[ , lapply(.SD, mean) , by=plot_cluster],stringsAsFactors=F)
tmp_cpb_os_rank <- data.frame(apply(tmp_cpb_os_rank[,2:ncol(tmp_cpb_os_rank)],2,rank),row.names=tmp_cpb_os_rank[,1],stringsAsFactors=F)
m_cpb_os_rank <- reshape2::melt(as.matrix(tmp_cpb_os_rank))
m_cpb_os_rank$Var1 <- factor(as.character(m_cpb_os_rank$Var1),levels=rev(as.character(arrange(data.table(m_cpb_os_rank)[,median(value),by=Var1],V1)$Var1)))
g_osgene_cpb_ranks <- ggplot(m_cpb_os_rank,aes(x=Var1,y=value)) + geom_boxplot(aes(fill=Var1),outlier.shape=NA) + geom_jitter(height=0) + theme_classic(base_size=16) +
theme(legend.position="none",axis.text.x=element_text(angle=90,hjust=1)) + xlab("") + ylab("Melanoma subtype expression rank") +
scale_fill_manual(values=temp_colors[rev(as.character(arrange(data.table(m_cpb_os_rank)[,median(value),by=Var1],V1)$Var1))])
ggsave("figure3_outputs/S19e.OSgene_subtype_expression_rank.pdf",g_osgene_cpb_ranks,height=7,width=12)
knitr::include_graphics("figure3_outputs/S19e.OSgene_subtype_expression_rank.pdf")
```

<embed src="figure3_outputs/S19e.OSgene_subtype_expression_rank.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
dt_cpb_os_zscore <- data.table(t(fcmeta_tab[rownames(val_df_res_os_sig[val_df_res_os_sig$log2FoldChange>0,]),]))
colnames(dt_cpb_os_zscore) <- val_df_res_os_sig[val_df_res_os_sig$log2FoldChange>0,]$gene
dt_cpb_os_zscore <- data.table(apply(dt_cpb_os_zscore,2,scale))
dt_cpb_os_zscore$plot_cluster <- clin$plot_cluster
df_cpb_os_zscore <- data.frame(dt_cpb_os_zscore[ , lapply(.SD, mean) , by=plot_cluster],stringsAsFactors=F)
rownames(df_cpb_os_zscore) <- df_cpb_os_zscore$plot_cluster
df_cpb_os_zscore$plot_cluster <- NULL

pdf("figure3_outputs/S19d.OSgene_subtype_zscore_nobreaks.pdf",height=4,width=6)
pheatmap(df_cpb_os_zscore,color=rev(brewer.pal(n=11, name="RdBu")),fontsize_col=8,fontsize_row=8,cellheight=8,cellwidth=8,treeheight_row=15,treeheight_col=15)
dev.off()
```

    ## null device 
    ##           1

``` {.r}
knitr::include_graphics("figure3_outputs/S19d.OSgene_subtype_zscore_nobreaks.pdf")
```

<embed src="figure3_outputs/S19d.OSgene_subtype_zscore_nobreaks.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
r_nr_genes <- sort(val_df_res_sig$gene)
long_short_os_genes <- sort(val_df_res_os_sig$gene)

print(phyper(length(intersect(r_nr_genes,long_short_os_genes))-1,length(long_short_os_genes),12134,length(r_nr_genes),lower.tail=FALSE,log.p=FALSE))
```

    ## [1] 6.783916e-41

``` {.r}
#6.783916e-41

futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
```

    ## NULL

``` {.r}
pdf("figure3_outputs/S18b.RNR_longshortOS_euler.pdf")
group.venn(list("R and NR genes"=r_nr_genes,"Long and Short OS genes"=long_short_os_genes), label=TRUE,fill=c("red","blue"),lab.cex=0.3,cat.pos = c(0, 0))
dev.off()
```

    ## null device 
    ##           1

``` {.r}
knitr::include_graphics("figure3_outputs/S18b.RNR_longshortOS_euler.pdf")
```

<embed src="figure3_outputs/S18b.RNR_longshortOS_euler.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
val_df_res_sig$CCLE_median_TPM <- as.numeric(apply(pcaexl_ccle_melanoma_tab[rownames(val_df_res_sig),],1,median))
val_df_res_sig$human_protein_atlas_median <- log2(as.numeric(apply(hpa_tab[unlist(lapply(strsplit(rownames(val_df_res_sig),"\\."),"[",1)),],1,median))+1)

g_val_rnr_CCLE_HPA <- ggplot(val_df_res_sig[!is.na(val_df_res_sig$human_protein_atlas_median),],aes(x=human_protein_atlas_median,y=CCLE_median_TPM)) +
geom_point(aes(color=log2FoldChange>0)) + geom_vline(xintercept=1) +  geom_hline(yintercept=1) + theme_classic(base_size=16) + scale_color_manual(values=c("red","blue")) +
theme(legend.position="bottom") + geom_text_repel(aes(label=gene),size=3,max.overlaps=100) + geom_abline(slope=1,intercept=0) + xlim(0,13) + ylim(0,7) +
xlab("HPA Blood median log2(TPM+1)") + ylab("CCLE Melanoma median log2(TPM+1)")
ggsave("figure3_outputs/S18.c.CCLE_HPA_expression_RNRgenes.pdf",g_val_rnr_CCLE_HPA)
knitr::include_graphics("figure3_outputs/S18.c.CCLE_HPA_expression_RNRgenes.pdf")
```

<embed src="figure3_outputs/S18.c.CCLE_HPA_expression_RNRgenes.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
tmp_cor_res <- cor(t(fcmeta_tab[rownames(val_df_res_sig),]))
rownames(tmp_cor_res) <- val_df_res_sig$gene
colnames(tmp_cor_res) <- val_df_res_sig$gene
m_tmp_cor_res <- reshape2::melt(tmp_cor_res)
m_tmp_cor_res$gene_pair <- paste0(m_tmp_cor_res$Var1,"_",m_tmp_cor_res$Var2)

pdf("figure3_outputs/S20a.val_gene_cor.pdf")
corrplot(tmp_cor_res, tl.cex=0.35,type = "upper", order = "hclust",tl.col = "black", tl.srt = 90,col=colorRampPalette(c("blue","white","red"))(200))
dev.off()
```

    ## null device 
    ##           1

``` {.r}
knitr::include_graphics("figure3_outputs/S20a.val_gene_cor.pdf")
```

<embed src="figure3_outputs/S20a.val_gene_cor.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
dt_hpa_rnr <- data.table(t(log2(hpa_tab[unlist(lapply(strsplit(rownames(val_df_res_sig[val_df_res_sig$log2FoldChange>0,]),"\\."),"[",1)),]+1)))
colnames(dt_hpa_rnr) <- val_df_res_sig[val_df_res_sig$log2FoldChange>0,]$gene
dt_hpa_rnr$cell_type <- unlist(lapply(strsplit(colnames(hpa_tab),"\\.\\.\\."),"[",1))
tmp_hpa_rnr_rank <- data.frame(dt_hpa_rnr[ , lapply(.SD, mean) , by=cell_type],stringsAsFactors=F)
tmp_hpa_rnr_rank <- data.frame(apply(tmp_hpa_rnr_rank[,2:ncol(tmp_hpa_rnr_rank)],2,rank),row.names=tmp_hpa_rnr_rank[,1],stringsAsFactors=F)
m_hpa_rnr_rank <- reshape2::melt(as.matrix(tmp_hpa_rnr_rank))
m_hpa_rnr_rank$Var1 <- gsub("\\."," ",m_hpa_rnr_rank$Var1)
m_hpa_rnr_rank$Var1 <- factor(as.character(m_hpa_rnr_rank$Var1),levels=rev(as.character(arrange(data.table(m_hpa_rnr_rank)[,median(value),by=Var1],V1)$Var1)))
g_rnrgene_hpa_ranks <- ggplot(m_hpa_rnr_rank,aes(x=Var1,y=value)) + geom_boxplot(aes(fill=Var1),outlier.shape=NA) + geom_jitter(height=0) + theme_classic(base_size=16) +
theme(legend.position="none",axis.text.x=element_text(angle=90,hjust=1)) + xlab("") + ylab("HPA Blood cell type expression rank")
ggsave("figure3_outputs/S20c.RNRgene_HPA_expression_rank.pdf",g_rnrgene_hpa_ranks,height=6,width=12)
knitr::include_graphics("figure3_outputs/S20c.RNRgene_HPA_expression_rank.pdf")
```

<embed src="figure3_outputs/S20c.RNRgene_HPA_expression_rank.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
dt_hpa_rnr_zscore <- data.table(t(log2(hpa_tab[unlist(lapply(strsplit(rownames(val_df_res_sig[val_df_res_sig$log2FoldChange>0,]),"\\."),"[",1)),]+1)))
colnames(dt_hpa_rnr_zscore) <- val_df_res_sig[val_df_res_sig$log2FoldChange>0,]$gene
dt_hpa_rnr_zscore <- data.table(apply(dt_hpa_rnr_zscore,2,scale))
dt_hpa_rnr_zscore$cell_type <- unlist(lapply(strsplit(colnames(hpa_tab),"\\.\\.\\."),"[",1))
df_hpa_rnr_zscore <- data.frame(dt_hpa_rnr_zscore[ , lapply(.SD, mean) , by=cell_type],stringsAsFactors=F)
rownames(df_hpa_rnr_zscore) <- gsub("\\."," ",df_hpa_rnr_zscore$cell_type)
df_hpa_rnr_zscore$cell_type <- NULL

pdf("figure3_outputs/S20b.RNRgene_HPA_zscore_nobreaks.pdf",height=7,width=9)
pheatmap(df_hpa_rnr_zscore[,!(colnames(df_hpa_rnr_zscore) %in% c("FDCSP","C4A","CCL21"))],color=rev(brewer.pal(n=11, name="RdBu")),fontsize_col=7,fontsize_row=7,cellheight=6,
cellwidth=6,treeheight_row=15,treeheight_col=15)
dev.off()
```

    ## null device 
    ##           1

``` {.r}
knitr::include_graphics("figure3_outputs/S20b.RNRgene_HPA_zscore_nobreaks.pdf")
```

<embed src="figure3_outputs/S20b.RNRgene_HPA_zscore_nobreaks.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
dt_cpb_rnr_zscore <- data.table(t(fcmeta_tab[rownames(val_df_res_sig[val_df_res_sig$log2FoldChange<0,]),]))
colnames(dt_cpb_rnr_zscore) <- val_df_res_sig[val_df_res_sig$log2FoldChange<0,]$gene
dt_cpb_rnr_zscore <- data.table(apply(dt_cpb_rnr_zscore,2,scale))
dt_cpb_rnr_zscore$plot_cluster <- clin$plot_cluster
df_cpb_rnr_zscore <- data.frame(dt_cpb_rnr_zscore[ , lapply(.SD, mean) , by=plot_cluster],stringsAsFactors=F)
rownames(df_cpb_rnr_zscore) <- df_cpb_rnr_zscore$plot_cluster
df_cpb_rnr_zscore$plot_cluster <- NULL

pdf("figure3_outputs/S20d.RNRgene_subtype_zscore_nobreaks.pdf",height=4,width=6)
pheatmap(df_cpb_rnr_zscore,color=rev(brewer.pal(n=11, name="RdBu")),fontsize_col=8,fontsize_row=8,cellheight=8,cellwidth=8,treeheight_row=15,treeheight_col=15)
dev.off()
```

    ## null device 
    ##           1

``` {.r}
knitr::include_graphics("figure3_outputs/S20d.RNRgene_subtype_zscore_nobreaks.pdf")
```

<embed src="figure3_outputs/S20d.RNRgene_subtype_zscore_nobreaks.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
dt_cpb_rnr <- data.table(t(fcmeta_tab[rownames(val_df_res_sig[val_df_res_sig$log2FoldChange<0,]),]))
colnames(dt_cpb_rnr) <- val_df_res_sig[val_df_res_sig$log2FoldChange<0,]$gene
dt_cpb_rnr$plot_cluster <- clin$plot_cluster
tmp_cpb_rnr_rank <- data.frame(dt_cpb_rnr[ , lapply(.SD, mean) , by=plot_cluster],stringsAsFactors=F)
tmp_cpb_rnr_rank <- data.frame(apply(tmp_cpb_rnr_rank[,2:ncol(tmp_cpb_rnr_rank)],2,rank),row.names=tmp_cpb_rnr_rank[,1],stringsAsFactors=F)
m_cpb_rnr_rank <- reshape2::melt(as.matrix(tmp_cpb_rnr_rank))
m_cpb_rnr_rank$Var1 <- factor(as.character(m_cpb_rnr_rank$Var1),levels=rev(as.character(arrange(data.table(m_cpb_rnr_rank)[,median(value),by=Var1],V1)$Var1)))
g_rnrgene_cpb_ranks <- ggplot(m_cpb_rnr_rank,aes(x=Var1,y=value)) + geom_boxplot(aes(fill=Var1),outlier.shape=NA) + geom_jitter(height=0) + theme_classic(base_size=16) +
theme(legend.position="none",axis.text.x=element_text(angle=90,hjust=1)) + xlab("") + ylab("Melanoma subtype expression rank") +
scale_fill_manual(values=temp_colors[rev(as.character(arrange(data.table(m_cpb_rnr_rank)[,median(value),by=Var1],V1)$Var1))])
ggsave("figure3_outputs/S20e.RNRgene_subtype_expression_rank.pdf",g_rnrgene_cpb_ranks,height=7,width=12)
knitr::include_graphics("figure3_outputs/S20e.RNRgene_subtype_expression_rank.pdf")
```

<embed src="figure3_outputs/S20e.RNRgene_subtype_expression_rank.pdf" width="600px" height="600px" type="application/pdf" />

Figure 4

``` {.r}
dir.create("figure4_outputs", showWarnings = FALSE)

mat_m_tmp_cor_res <- as.matrix(tmp_cor_res)
mat_m_tmp_cor_res[upper.tri(mat_m_tmp_cor_res)] <- NA
diag(mat_m_tmp_cor_res) <- NA
fix_tmp_mat_m_tmp_cor_res <- reshape2::melt(mat_m_tmp_cor_res)
fix_tmp_mat_m_tmp_cor_res$gene_pair <- paste0(fix_tmp_mat_m_tmp_cor_res$Var1,"_",fix_tmp_mat_m_tmp_cor_res$Var2)

mat_m_tmp_cor_res_os <- as.matrix(tmp_cor_res_os)
mat_m_tmp_cor_res_os[upper.tri(mat_m_tmp_cor_res_os)] <- NA
diag(mat_m_tmp_cor_res_os) <- NA
fix_tmp_mat_m_tmp_cor_res_os <- reshape2::melt(mat_m_tmp_cor_res_os)
fix_tmp_mat_m_tmp_cor_res_os$gene_pair <- paste0(fix_tmp_mat_m_tmp_cor_res_os$Var1,"_",fix_tmp_mat_m_tmp_cor_res_os$Var2)

half_m_tmp_cor_res <- m_tmp_cor_res[m_tmp_cor_res$gene_pair %in% fix_tmp_mat_m_tmp_cor_res[!is.na(fix_tmp_mat_m_tmp_cor_res$value),]$gene_pair,]
half_m_tmp_cor_res_os <- m_tmp_cor_res_os[m_tmp_cor_res_os$gene_pair %in% fix_tmp_mat_m_tmp_cor_res_os[!is.na(fix_tmp_mat_m_tmp_cor_res_os$value),]$gene_pair,]

tmp_val_df_res_sig <- val_df_res_sig[,c("gene","log2FoldChange")]
rownames(tmp_val_df_res_sig) <- tmp_val_df_res_sig$gene
tmp_val_df_res_os_sig <- val_df_res_os_sig[,c("gene","log2FoldChange")]
rownames(tmp_val_df_res_os_sig) <- tmp_val_df_res_os_sig$gene

half_m_tmp_cor_res$Var1_log2FoldChange <- tmp_val_df_res_sig[half_m_tmp_cor_res$Var1,]$log2FoldChange
half_m_tmp_cor_res$Var2_log2FoldChange <- tmp_val_df_res_sig[half_m_tmp_cor_res$Var2,]$log2FoldChange
half_m_tmp_cor_res_os$Var1_log2FoldChange <- tmp_val_df_res_os_sig[half_m_tmp_cor_res_os$Var1,]$log2FoldChange
half_m_tmp_cor_res_os$Var2_log2FoldChange <- tmp_val_df_res_os_sig[half_m_tmp_cor_res_os$Var2,]$log2FoldChange

half_m_tmp_cor_res$model_type <- ifelse(half_m_tmp_cor_res$Var1_log2FoldChange>0&half_m_tmp_cor_res$Var2_log2FoldChange>0,"R gene/R gene",
ifelse(half_m_tmp_cor_res$Var1_log2FoldChange<0&half_m_tmp_cor_res$Var2_log2FoldChange<0,"NR gene/NR gene","R gene/NR gene"))
half_m_tmp_cor_res$model_type <- factor(half_m_tmp_cor_res$model_type,levels=c("NR gene/NR gene","R gene/R gene","R gene/NR gene"))
half_m_tmp_cor_res_os$model_type <- ifelse(half_m_tmp_cor_res_os$Var1_log2FoldChange<0& half_m_tmp_cor_res_os$Var2_log2FoldChange<0,"High OS gene/High OS gene",
ifelse(half_m_tmp_cor_res_os$Var1_log2FoldChange>0&half_m_tmp_cor_res_os$Var2_log2FoldChange>0,"Low OS gene/Low OS gene","High OS gene/Low OS gene"))
half_m_tmp_cor_res_os$model_type <- factor(half_m_tmp_cor_res_os$model_type,levels=c("Low OS gene/Low OS gene","High OS gene/High OS gene","High OS gene/Low OS gene"))

tmp_half_m_tmp_cor_res <- half_m_tmp_cor_res
tmp_half_m_tmp_cor_res_os <- half_m_tmp_cor_res_os

#add mean or r/nr z and high/low os gene z models to half_m_tmp_cor_res/half_m_tmp_cor_res_os
#mean of r genes z + mean of nr genes z
#mean of highos genes z + mean of lowos genes z
clin[,paste0(val_df_res_sig$gene,"_z")] <- apply(fcmeta_tab[rownames(val_df_res_sig),],1,scale)
clin$mean_r_genes_z <- rowMeans(clin[,paste0(val_df_res_sig[val_df_res_sig$log2FoldChange>0,]$gene,"_z")])
clin$mean_nr_genes_z <- rowMeans(clin[,paste0(val_df_res_sig[val_df_res_sig$log2FoldChange<0,]$gene,"_z")])

clin[,paste0(val_df_res_os_sig$gene,"_z")] <- apply(fcmeta_tab[rownames(val_df_res_os_sig),],1,scale)
clin$mean_lowos_genes_z <- rowMeans(clin[,paste0(val_df_res_os_sig[val_df_res_os_sig$log2FoldChange>0,]$gene,"_z")])
clin$mean_highos_genes_z <- rowMeans(clin[,paste0(val_df_res_os_sig[val_df_res_os_sig$log2FoldChange<0,]$gene,"_z")])

tmpdf <- half_m_tmp_cor_res[half_m_tmp_cor_res$model_type=="R gene/NR gene",][1,]
tmpdf$Var1 <- "mean_r_genes_z"
tmpdf$Var2 <- "mean_nr_genes_z"
tmpdf$value <- cor(clin$mean_r_genes_z,clin$mean_nr_genes_z,method="spearman")
tmpdf$gene_pair <- paste0(tmpdf$Var1,"_",tmpdf$Var2)
tmpdf$Var1_log2FoldChange <- mean(val_df_res_sig[val_df_res_sig$log2FoldChange>0,]$log2FoldChange)
tmpdf$Var2_log2FoldChange <- mean(val_df_res_sig[val_df_res_sig$log2FoldChange<0,]$log2FoldChange)

half_m_tmp_cor_res <- rbind.data.frame(half_m_tmp_cor_res,tmpdf)

tmpdf_os <- half_m_tmp_cor_res_os[half_m_tmp_cor_res_os$model_type=="High OS gene/Low OS gene",][1,]
tmpdf_os$Var1 <- "mean_highos_genes_z"
tmpdf_os$Var2 <- "mean_lowos_genes_z"
tmpdf_os$value <- cor(clin$mean_highos_genes_z,clin$mean_lowos_genes_z,method="spearman")
tmpdf_os$gene_pair <- paste0(tmpdf_os$Var1,"_",tmpdf_os$Var2)
tmpdf_os$Var1_log2FoldChange <- mean(val_df_res_os_sig[val_df_res_os_sig$log2FoldChange<0,]$log2FoldChange)
tmpdf_os$Var2_log2FoldChange <- mean(val_df_res_os_sig[val_df_res_os_sig$log2FoldChange>0,]$log2FoldChange)

half_m_tmp_cor_res_os <- rbind.data.frame(half_m_tmp_cor_res_os,tmpdf_os)

print("Running Cox OS models for R/NR gene pairs and High/Low OS gene pairs")
```

    ## [1] "Running Cox OS models for R/NR gene pairs and High/Low OS gene pairs"

``` {.r}
all_single_strings <- c(paste0("`",unique(c(val_df_res_sig$gene,val_df_res_os_sig$gene)),"`"),"mean_r_genes_z","mean_nr_genes_z","mean_highos_genes_z","mean_lowos_genes_z")
print(paste0("Running ",length(all_single_strings)," single gene Cox OS models"))
```

    ## [1] "Running 159 single gene Cox OS models"

``` {.r}
all_single_formulas <- sapply(all_single_strings,function(x) as.formula(paste('Surv(overall_survival,DOD=="DOD")~',x)))
all_single_models <- lapply(all_single_formulas, function(x){coxph(x,data=clin)})
s_all_single_models <- lapply(all_single_models,summary)

half_m_tmp_cor_res$gene1_continuous_logrank_p <- NA
half_m_tmp_cor_res$gene1_continuous_os_hr <- NA
half_m_tmp_cor_res$gene1_continuous_os_hr_lo <- NA
half_m_tmp_cor_res$gene1_continuous_os_hr_hi <- NA
half_m_tmp_cor_res$gene1_continuous_cindex <- NA
half_m_tmp_cor_res$gene2_continuous_logrank_p <- NA
half_m_tmp_cor_res$gene2_continuous_os_hr <- NA
half_m_tmp_cor_res$gene2_continuous_os_hr_lo <- NA
half_m_tmp_cor_res$gene2_continuous_os_hr_hi <- NA
half_m_tmp_cor_res$gene2_continuous_cindex <- NA

half_m_tmp_cor_res_os$gene1_continuous_logrank_p <- NA
half_m_tmp_cor_res_os$gene1_continuous_os_hr <- NA
half_m_tmp_cor_res_os$gene1_continuous_os_hr_lo <- NA
half_m_tmp_cor_res_os$gene1_continuous_os_hr_hi <- NA
half_m_tmp_cor_res_os$gene1_continuous_cindex <- NA
half_m_tmp_cor_res_os$gene2_continuous_logrank_p <- NA
half_m_tmp_cor_res_os$gene2_continuous_os_hr <- NA
half_m_tmp_cor_res_os$gene2_continuous_os_hr_lo <- NA
half_m_tmp_cor_res_os$gene2_continuous_os_hr_hi <- NA
half_m_tmp_cor_res_os$gene2_continuous_cindex <- NA

single_model_result_list <- lapply(s_all_single_models,function(x){return(cbind(
                                                        gene1_continuous_logrank_p=x$sctest[3],
                                                        gene1_continuous_os_hr=x$conf.int[1,1],
                                                        gene1_continuous_os_hr_lo=x$conf.int[1,3],
                                                        gene1_continuous_os_hr_hi=x$conf.int[1,4],
                                                        gene1_continuous_cindex=x$concordance[1],
                                                        gene2_continuous_logrank_p=x$sctest[3],
                                                        gene2_continuous_os_hr=x$conf.int[1,1],
                                                        gene2_continuous_os_hr_lo=x$conf.int[1,3],
                                                        gene2_continuous_os_hr_hi=x$conf.int[1,4],
                                                        gene2_continuous_cindex=x$concordance[1]
                                                        ))})

single_model_df <- as.data.frame(data.table::transpose(single_model_result_list), col.names = colnames(single_model_result_list[[1]]))
single_model_df$gene <- c(unique(c(val_df_res_sig$gene,val_df_res_os_sig$gene)),"mean_r_genes_z","mean_nr_genes_z","mean_highos_genes_z","mean_lowos_genes_z")
rownames(single_model_df) <- single_model_df$gene

half_m_tmp_cor_res[,c("gene1_continuous_logrank_p","gene1_continuous_os_hr","gene1_continuous_os_hr_lo","gene1_continuous_os_hr_hi","gene1_continuous_cindex")] <- single_model_df[
        as.character(half_m_tmp_cor_res$Var1),c("gene1_continuous_logrank_p","gene1_continuous_os_hr","gene1_continuous_os_hr_lo","gene1_continuous_os_hr_hi","gene1_continuous_cindex")
]
half_m_tmp_cor_res[,c("gene2_continuous_logrank_p","gene2_continuous_os_hr","gene2_continuous_os_hr_lo","gene2_continuous_os_hr_hi","gene2_continuous_cindex")] <- single_model_df[
        as.character(half_m_tmp_cor_res$Var2),c("gene2_continuous_logrank_p","gene2_continuous_os_hr","gene2_continuous_os_hr_lo","gene2_continuous_os_hr_hi","gene2_continuous_cindex")
]

half_m_tmp_cor_res_os[,c("gene1_continuous_logrank_p","gene1_continuous_os_hr","gene1_continuous_os_hr_lo","gene1_continuous_os_hr_hi","gene1_continuous_cindex")] <- single_model_df[
        as.character(half_m_tmp_cor_res_os$Var1),c("gene1_continuous_logrank_p","gene1_continuous_os_hr","gene1_continuous_os_hr_lo","gene1_continuous_os_hr_hi","gene1_continuous_cindex")
]
half_m_tmp_cor_res_os[,c("gene2_continuous_logrank_p","gene2_continuous_os_hr","gene2_continuous_os_hr_lo","gene2_continuous_os_hr_hi","gene2_continuous_cindex")] <- single_model_df[
        as.character(half_m_tmp_cor_res_os$Var2),c("gene2_continuous_logrank_p","gene2_continuous_os_hr","gene2_continuous_os_hr_lo","gene2_continuous_os_hr_hi","gene2_continuous_cindex")
]

all_pairwise_strings <- c(paste0("`",gsub("_","`+`",unique(c(tmp_half_m_tmp_cor_res$gene_pair,tmp_half_m_tmp_cor_res_os$gene_pair))),"`"),
"mean_r_genes_z+mean_nr_genes_z","mean_highos_genes_z+mean_lowos_genes_z")

print(paste0("Running ",length(all_pairwise_strings)," pairwise gene Cox OS models"))
```

    ## [1] "Running 8049 pairwise gene Cox OS models"

``` {.r}
all_pairwise_formulas <- sapply(all_pairwise_strings,function(x)as.formula(paste('Surv(overall_survival,DOD=="DOD")~',x)))
all_pairwise_models <- lapply(all_pairwise_formulas, function(x){coxph(x,data=clin)})
```

    ## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, :
    ## Loglik converged before variable 1 ; beta may be infinite.

``` {.r}
s_all_pairwise_models <- lapply(all_pairwise_models,summary)

half_m_tmp_cor_res$gene1gene2_continuous_logrank_p <- NA
half_m_tmp_cor_res$gene1gene2_continuous_os_hr_1 <- NA
half_m_tmp_cor_res$gene1gene2_continuous_os_hr_1_lo <- NA
half_m_tmp_cor_res$gene1gene2_continuous_os_hr_1_hi <- NA
half_m_tmp_cor_res$gene1gene2_continuous_cox_p_1 <- NA
half_m_tmp_cor_res$gene1gene2_continuous_os_hr_2 <- NA
half_m_tmp_cor_res$gene1gene2_continuous_os_hr_2_lo <- NA
half_m_tmp_cor_res$gene1gene2_continuous_os_hr_2_hi <- NA
half_m_tmp_cor_res$gene1gene2_continuous_cox_p_2 <- NA
half_m_tmp_cor_res$gene1gene2_continuous_cindex <- NA

half_m_tmp_cor_res_os$gene1gene2_continuous_logrank_p <- NA
half_m_tmp_cor_res_os$gene1gene2_continuous_os_hr_1 <- NA
half_m_tmp_cor_res_os$gene1gene2_continuous_os_hr_1_lo <- NA
half_m_tmp_cor_res_os$gene1gene2_continuous_os_hr_1_hi <- NA
half_m_tmp_cor_res_os$gene1gene2_continuous_cox_p_1 <- NA
half_m_tmp_cor_res_os$gene1gene2_continuous_os_hr_2 <- NA
half_m_tmp_cor_res_os$gene1gene2_continuous_os_hr_2_lo <- NA
half_m_tmp_cor_res_os$gene1gene2_continuous_os_hr_2_hi <- NA
half_m_tmp_cor_res_os$gene1gene2_continuous_cox_p_2 <- NA
half_m_tmp_cor_res_os$gene1gene2_continuous_cindex <- NA

pairwise_model_result_list <- lapply(s_all_pairwise_models,function(x){return(cbind(
                                                        gene1gene2_continuous_logrank_p=x$sctest[3],
                                                        gene1gene2_continuous_os_hr_1=x$conf.int[1,1],
                                                        gene1gene2_continuous_os_hr_1_lo=x$conf.int[1,3],
                                                        gene1gene2_continuous_os_hr_1_hi=x$conf.int[1,4],
                                                        gene1gene2_continuous_cox_p_1=x$coef[1,"Pr(>|z|)"],
                                                        gene1gene2_continuous_os_hr_2=x$conf.int[2,1],
                                                        gene1gene2_continuous_os_hr_2_lo=x$conf.int[2,3],
                                                        gene1gene2_continuous_os_hr_2_hi=x$conf.int[2,4],
                                                        gene1gene2_continuous_cox_p_2=x$coef[2,"Pr(>|z|)"],
                                                        gene1gene2_continuous_cindex=x$concordance[1]
                                                        ))})
pairwise_model_df <- as.data.frame(data.table::transpose(pairwise_model_result_list), col.names = colnames(pairwise_model_result_list[[1]]))
pairwise_model_df$gene_pair <- c(unique(c(tmp_half_m_tmp_cor_res$gene_pair,tmp_half_m_tmp_cor_res_os$gene_pair)),"mean_r_genes_z_mean_nr_genes_z","mean_highos_genes_z_mean_lowos_genes_z")
rownames(pairwise_model_df) <- pairwise_model_df$gene_pair

half_m_tmp_cor_res[,c("gene1gene2_continuous_logrank_p","gene1gene2_continuous_os_hr_1","gene1gene2_continuous_os_hr_1_lo","gene1gene2_continuous_os_hr_1_hi",
                                        "gene1gene2_continuous_cox_p_1","gene1gene2_continuous_os_hr_2","gene1gene2_continuous_os_hr_2_lo","gene1gene2_continuous_os_hr_2_hi",
                                        "gene1gene2_continuous_cox_p_2","gene1gene2_continuous_cindex")] <- pairwise_model_df[
                as.character(half_m_tmp_cor_res$gene_pair),
                c("gene1gene2_continuous_logrank_p","gene1gene2_continuous_os_hr_1","gene1gene2_continuous_os_hr_1_lo","gene1gene2_continuous_os_hr_1_hi",
                "gene1gene2_continuous_cox_p_1","gene1gene2_continuous_os_hr_2","gene1gene2_continuous_os_hr_2_lo","gene1gene2_continuous_os_hr_2_hi",
                "gene1gene2_continuous_cox_p_2","gene1gene2_continuous_cindex")
                ]

half_m_tmp_cor_res_os[,c("gene1gene2_continuous_logrank_p","gene1gene2_continuous_os_hr_1","gene1gene2_continuous_os_hr_1_lo","gene1gene2_continuous_os_hr_1_hi",
                                        "gene1gene2_continuous_cox_p_1","gene1gene2_continuous_os_hr_2","gene1gene2_continuous_os_hr_2_lo","gene1gene2_continuous_os_hr_2_hi",
                                        "gene1gene2_continuous_cox_p_2","gene1gene2_continuous_cindex")] <- pairwise_model_df[
                as.character(half_m_tmp_cor_res_os$gene_pair),
                c("gene1gene2_continuous_logrank_p","gene1gene2_continuous_os_hr_1","gene1gene2_continuous_os_hr_1_lo","gene1gene2_continuous_os_hr_1_hi",
                "gene1gene2_continuous_cox_p_1","gene1gene2_continuous_os_hr_2","gene1gene2_continuous_os_hr_2_lo","gene1gene2_continuous_os_hr_2_hi",
                "gene1gene2_continuous_cox_p_2","gene1gene2_continuous_cindex")
                ]

half_m_tmp_cor_res$gene1gene2_continuous_logrank_bonferroni_p <- p.adjust(half_m_tmp_cor_res$gene1gene2_continuous_logrank_p,method="bonferroni")
half_m_tmp_cor_res_os$gene1gene2_continuous_logrank_bonferroni_p <- p.adjust(half_m_tmp_cor_res_os$gene1gene2_continuous_logrank_p,method="bonferroni")

#calculate R/NR AUC values

print("Running R/NR classification models for R/NR gene pairs and High/Low OS gene pairs")
```

    ## [1] "Running R/NR classification models for R/NR gene pairs and High/Low OS gene pairs"

``` {.r}
twosd <<- qnorm(0.975,0,1)

clin$Response_numeric <- as.numeric(as.factor(clin$Response))-1

all_single_formulas_glm <- sapply(all_single_strings,function(x) as.formula(paste('Response_numeric~',x)))
print(paste0("Running ",length(all_single_strings)," single gene R/NR AUC models"))
```

    ## [1] "Running 159 single gene R/NR AUC models"

``` {.r}
all_single_models_glm <- lapply(all_single_formulas_glm, function(x){glm(x,family=binomial(link="logit"),data=clin)})
s_all_single_models_glm <- lapply(all_single_models_glm,summary)
all_single_models_rocobj <- lapply(all_single_models_glm,function(x){ifelse(roc.area(clin$Response_numeric,x$fitted.values)$A<0.5,
                           return(roc.area(clin$Response_numeric,1/x$fitted.values)),
                           return(roc.area(clin$Response_numeric,x$fitted.values)))})
single_model_result_list_glm <- mapply(function(x,y){return(cbind(
                                                        gene1_continuous_rnr_logistic_p=x$coef[2,"Pr(>|z|)"],
                                                        gene1_continuous_rnr_logistic_or=x$coef[2,"Estimate"],
                                                        gene1_continuous_rnr_logistic_or_lo=x$coef[2,"Estimate"]-twosd*x$coef[2,"Std. Error"],
                                                        gene1_continuous_rnr_logistic_or_hi=x$coef[2,"Estimate"]+twosd*x$coef[2,"Std. Error"],
                                                        gene1_continuous_rnr_roc_auc=y$A,
                                                        gene1_continuous_rnr_roc_p=y$p.value,
                                                        gene2_continuous_rnr_logistic_p=x$coef[2,"Pr(>|z|)"],
                                                        gene2_continuous_rnr_logistic_or=x$coef[2,"Estimate"],
                                                        gene2_continuous_rnr_logistic_or_lo=x$coef[2,"Estimate"]-twosd*x$coef[2,"Std. Error"],
                                                        gene2_continuous_rnr_logistic_or_hi=x$coef[2,"Estimate"]+twosd*x$coef[2,"Std. Error"],
                                                        gene2_continuous_rnr_roc_auc=y$A,
                                                        gene2_continuous_rnr_roc_p=y$p.value
                                                        ))},x=s_all_single_models_glm,y=all_single_models_rocobj)

single_model_df_glm <- data.frame(t(single_model_result_list_glm), stringsAsFactors=F)
colnames(single_model_df_glm) <- c("gene1_continuous_rnr_logistic_p","gene1_continuous_rnr_logistic_or","gene1_continuous_rnr_logistic_or_lo","gene1_continuous_rnr_logistic_or_hi",
                         "gene1_continuous_rnr_roc_auc","gene1_continuous_rnr_roc_p","gene2_continuous_rnr_logistic_p","gene2_continuous_rnr_logistic_or","gene2_continuous_rnr_logistic_or_lo",
                         "gene2_continuous_rnr_logistic_or_hi","gene2_continuous_rnr_roc_auc","gene2_continuous_rnr_roc_p")
single_model_df_glm$gene <- c(unique(c(val_df_res_sig$gene,val_df_res_os_sig$gene)),"mean_r_genes_z","mean_nr_genes_z","mean_highos_genes_z","mean_lowos_genes_z")
rownames(single_model_df_glm) <- single_model_df_glm$gene

half_m_tmp_cor_res$gene1_continuous_rnr_logistic_p <- NA
half_m_tmp_cor_res$gene1_continuous_rnr_logistic_or <- NA
half_m_tmp_cor_res$gene1_continuous_rnr_logistic_or_lo <- NA
half_m_tmp_cor_res$gene1_continuous_rnr_logistic_or_hi <- NA
half_m_tmp_cor_res$gene1_continuous_rnr_roc_auc <- NA
half_m_tmp_cor_res$gene1_continuous_rnr_roc_p <- NA
half_m_tmp_cor_res$gene2_continuous_rnr_logistic_p <- NA
half_m_tmp_cor_res$gene2_continuous_rnr_logistic_or <- NA
half_m_tmp_cor_res$gene2_continuous_rnr_logistic_or_lo <- NA
half_m_tmp_cor_res$gene2_continuous_rnr_logistic_or_hi <- NA
half_m_tmp_cor_res$gene2_continuous_rnr_roc_auc <- NA
half_m_tmp_cor_res$gene2_continuous_rnr_roc_p <- NA

half_m_tmp_cor_res_os$gene1_continuous_rnr_logistic_p <- NA
half_m_tmp_cor_res_os$gene1_continuous_rnr_logistic_or <- NA
half_m_tmp_cor_res_os$gene1_continuous_rnr_logistic_or_lo <- NA
half_m_tmp_cor_res_os$gene1_continuous_rnr_logistic_or_hi <- NA
half_m_tmp_cor_res_os$gene1_continuous_rnr_roc_auc <- NA
half_m_tmp_cor_res_os$gene1_continuous_rnr_roc_p <- NA
half_m_tmp_cor_res_os$gene2_continuous_rnr_logistic_p <- NA
half_m_tmp_cor_res_os$gene2_continuous_rnr_logistic_or <- NA
half_m_tmp_cor_res_os$gene2_continuous_rnr_logistic_or_lo <- NA
half_m_tmp_cor_res_os$gene2_continuous_rnr_logistic_or_hi <- NA
half_m_tmp_cor_res_os$gene2_continuous_rnr_roc_auc <- NA
half_m_tmp_cor_res_os$gene2_continuous_rnr_roc_p <- NA

half_m_tmp_cor_res[,c("gene1_continuous_rnr_logistic_p","gene1_continuous_rnr_logistic_or","gene1_continuous_rnr_logistic_or_lo","gene1_continuous_rnr_logistic_or_hi",
                                "gene1_continuous_rnr_roc_auc","gene1_continuous_rnr_roc_p")] <- single_model_df_glm[
                as.character(half_m_tmp_cor_res$Var1),c("gene1_continuous_rnr_logistic_p","gene1_continuous_rnr_logistic_or","gene1_continuous_rnr_logistic_or_lo",
                "gene1_continuous_rnr_logistic_or_hi","gene1_continuous_rnr_roc_auc","gene1_continuous_rnr_roc_p")
                ]

half_m_tmp_cor_res[,c("gene2_continuous_rnr_logistic_p","gene2_continuous_rnr_logistic_or","gene2_continuous_rnr_logistic_or_lo","gene2_continuous_rnr_logistic_or_hi",
                                "gene2_continuous_rnr_roc_auc","gene2_continuous_rnr_roc_p")] <- single_model_df_glm[
                as.character(half_m_tmp_cor_res$Var2),c("gene2_continuous_rnr_logistic_p","gene2_continuous_rnr_logistic_or","gene2_continuous_rnr_logistic_or_lo",
                "gene2_continuous_rnr_logistic_or_hi","gene2_continuous_rnr_roc_auc","gene2_continuous_rnr_roc_p")
                ]

half_m_tmp_cor_res_os[,c("gene1_continuous_rnr_logistic_p","gene1_continuous_rnr_logistic_or","gene1_continuous_rnr_logistic_or_lo","gene1_continuous_rnr_logistic_or_hi",
                                "gene1_continuous_rnr_roc_auc","gene1_continuous_rnr_roc_p")] <- single_model_df_glm[
                as.character(half_m_tmp_cor_res_os$Var1),c("gene1_continuous_rnr_logistic_p","gene1_continuous_rnr_logistic_or","gene1_continuous_rnr_logistic_or_lo",
                "gene1_continuous_rnr_logistic_or_hi","gene1_continuous_rnr_roc_auc","gene1_continuous_rnr_roc_p")
                ]

half_m_tmp_cor_res_os[,c("gene2_continuous_rnr_logistic_p","gene2_continuous_rnr_logistic_or","gene2_continuous_rnr_logistic_or_lo","gene2_continuous_rnr_logistic_or_hi",
                                "gene2_continuous_rnr_roc_auc","gene2_continuous_rnr_roc_p")] <- single_model_df_glm[
                as.character(half_m_tmp_cor_res_os$Var2),c("gene2_continuous_rnr_logistic_p","gene2_continuous_rnr_logistic_or","gene2_continuous_rnr_logistic_or_lo",
                "gene2_continuous_rnr_logistic_or_hi","gene2_continuous_rnr_roc_auc","gene2_continuous_rnr_roc_p")
                ]

all_pairwise_formulas_glm <- sapply(all_pairwise_strings,function(x) as.formula(paste('Response_numeric~',x)))
print(paste0("Running ",length(all_pairwise_strings)," pairwise gene R/NR classification models"))
```

    ## [1] "Running 8049 pairwise gene R/NR classification models"

``` {.r}
all_pairwise_models_glm <- lapply(all_pairwise_formulas_glm, function(x){glm(x,family=binomial(link="logit"),data=clin)})
s_all_pairwise_models_glm <- lapply(all_pairwise_models_glm,summary)
all_pairwise_models_rocobj <- lapply(all_pairwise_models_glm,function(x){return(roc.area(clin$Response_numeric,x$fitted.values))})
pairwise_model_result_list_glm <- mapply(function(x,y){return(cbind(
                                                        gene1gene2_continuous_rnr_logistic_or_1=x$coef[2,"Estimate"],
                                                        gene1gene2_continuous_rnr_logistic_or_1_lo=x$coef[2,"Estimate"]-twosd*x$coef[2,"Std. Error"],
                                                        gene1gene2_continuous_rnr_logistic_or_1_hi=x$coef[2,"Estimate"]+twosd*x$coef[2,"Std. Error"],
                                                        gene1gene2_continuous_rnr_logistic_p_1=x$coef[2,"Pr(>|z|)"],
                                                        gene1gene2_continuous_rnr_logistic_or_2=x$coef[3,"Estimate"],
                                                        gene1gene2_continuous_rnr_logistic_or_2_lo=x$coef[3,"Estimate"]-twosd*x$coef[3,"Std. Error"],
                                                        gene1gene2_continuous_rnr_logistic_or_2_hi=x$coef[3,"Estimate"]+twosd*x$coef[3,"Std. Error"],
                                                        gene1gene2_continuous_rnr_logistic_p_2=x$coef[3,"Pr(>|z|)"],
                                                        gene1gene2_continuous_rnr_roc_auc=y$A,
                                                        gene1gene2_continuous_rnr_roc_p=y$p.value
                                                        ))},x=s_all_pairwise_models_glm,y=all_pairwise_models_rocobj)
pairwise_model_df_glm <- data.frame(t(pairwise_model_result_list_glm), stringsAsFactors=F)
colnames(pairwise_model_df_glm) <- c("gene1gene2_continuous_rnr_logistic_or_1","gene1gene2_continuous_rnr_logistic_or_1_lo","gene1gene2_continuous_rnr_logistic_or_1_hi",
                                   "gene1gene2_continuous_rnr_logistic_p_1","gene1gene2_continuous_rnr_logistic_or_2","gene1gene2_continuous_rnr_logistic_or_2_lo",
                                   "gene1gene2_continuous_rnr_logistic_or_2_hi","gene1gene2_continuous_rnr_logistic_p_2","gene1gene2_continuous_rnr_roc_auc","gene1gene2_continuous_rnr_roc_p")
pairwise_model_df_glm$gene_pair <- c(unique(c(tmp_half_m_tmp_cor_res$gene_pair,tmp_half_m_tmp_cor_res_os$gene_pair)),
                                "mean_r_genes_z_mean_nr_genes_z","mean_highos_genes_z_mean_lowos_genes_z")
rownames(pairwise_model_df_glm) <- pairwise_model_df_glm$gene_pair

half_m_tmp_cor_res$gene1gene2_continuous_rnr_logistic_or_1 <- NA
half_m_tmp_cor_res$gene1gene2_continuous_rnr_logistic_or_1_lo <- NA
half_m_tmp_cor_res$gene1gene2_continuous_rnr_logistic_or_1_hi <- NA
half_m_tmp_cor_res$gene1gene2_continuous_rnr_logistic_p_1 <- NA
half_m_tmp_cor_res$gene1gene2_continuous_rnr_logistic_or_2 <- NA
half_m_tmp_cor_res$gene1gene2_continuous_rnr_logistic_or_2_lo <- NA
half_m_tmp_cor_res$gene1gene2_continuous_rnr_logistic_or_2_hi <- NA
half_m_tmp_cor_res$gene1gene2_continuous_rnr_logistic_p_2 <- NA
half_m_tmp_cor_res$gene1gene2_continuous_rnr_roc_auc <- NA
half_m_tmp_cor_res$gene1gene2_continuous_rnr_roc_p <- NA

half_m_tmp_cor_res_os$gene1gene2_continuous_rnr_logistic_or_1 <- NA
half_m_tmp_cor_res_os$gene1gene2_continuous_rnr_logistic_or_1_lo <- NA
half_m_tmp_cor_res_os$gene1gene2_continuous_rnr_logistic_or_1_hi <- NA
half_m_tmp_cor_res_os$gene1gene2_continuous_rnr_logistic_p_1 <- NA
half_m_tmp_cor_res_os$gene1gene2_continuous_rnr_logistic_or_2 <- NA
half_m_tmp_cor_res_os$gene1gene2_continuous_rnr_logistic_or_2_lo <- NA
half_m_tmp_cor_res_os$gene1gene2_continuous_rnr_logistic_or_2_hi <- NA
half_m_tmp_cor_res_os$gene1gene2_continuous_rnr_logistic_p_2 <- NA
half_m_tmp_cor_res_os$gene1gene2_continuous_rnr_roc_auc <- NA
half_m_tmp_cor_res_os$gene1gene2_continuous_rnr_roc_p <- NA

half_m_tmp_cor_res[,c("gene1gene2_continuous_rnr_logistic_or_1","gene1gene2_continuous_rnr_logistic_or_1_lo","gene1gene2_continuous_rnr_logistic_or_1_hi",
                                        "gene1gene2_continuous_rnr_logistic_p_1","gene1gene2_continuous_rnr_logistic_or_2","gene1gene2_continuous_rnr_logistic_or_2_lo",
                                        "gene1gene2_continuous_rnr_logistic_or_2_hi","gene1gene2_continuous_rnr_logistic_p_2","gene1gene2_continuous_rnr_roc_auc",
                                        "gene1gene2_continuous_rnr_roc_p")] <- pairwise_model_df_glm[as.character(half_m_tmp_cor_res$gene_pair),
                c("gene1gene2_continuous_rnr_logistic_or_1","gene1gene2_continuous_rnr_logistic_or_1_lo","gene1gene2_continuous_rnr_logistic_or_1_hi",
                "gene1gene2_continuous_rnr_logistic_p_1","gene1gene2_continuous_rnr_logistic_or_2","gene1gene2_continuous_rnr_logistic_or_2_lo",
                "gene1gene2_continuous_rnr_logistic_or_2_hi","gene1gene2_continuous_rnr_logistic_p_2","gene1gene2_continuous_rnr_roc_auc","gene1gene2_continuous_rnr_roc_p")
                ]

half_m_tmp_cor_res_os[,c("gene1gene2_continuous_rnr_logistic_or_1","gene1gene2_continuous_rnr_logistic_or_1_lo","gene1gene2_continuous_rnr_logistic_or_1_hi",
                                        "gene1gene2_continuous_rnr_logistic_p_1","gene1gene2_continuous_rnr_logistic_or_2","gene1gene2_continuous_rnr_logistic_or_2_lo",
                                        "gene1gene2_continuous_rnr_logistic_or_2_hi","gene1gene2_continuous_rnr_logistic_p_2","gene1gene2_continuous_rnr_roc_auc",
                                        "gene1gene2_continuous_rnr_roc_p")] <- pairwise_model_df_glm[as.character(half_m_tmp_cor_res_os$gene_pair),
                c("gene1gene2_continuous_rnr_logistic_or_1","gene1gene2_continuous_rnr_logistic_or_1_lo","gene1gene2_continuous_rnr_logistic_or_1_hi",
                "gene1gene2_continuous_rnr_logistic_p_1","gene1gene2_continuous_rnr_logistic_or_2","gene1gene2_continuous_rnr_logistic_or_2_lo",
                "gene1gene2_continuous_rnr_logistic_or_2_hi","gene1gene2_continuous_rnr_logistic_p_2","gene1gene2_continuous_rnr_roc_auc","gene1gene2_continuous_rnr_roc_p")
                ]

half_m_tmp_cor_res$gene1gene2_continuous_rnr_roc_bonferroni_p <- p.adjust(half_m_tmp_cor_res$gene1gene2_continuous_rnr_roc_p,method="bonferroni")
half_m_tmp_cor_res_os$gene1gene2_continuous_rnr_roc_bonferroni_p <- p.adjust(half_m_tmp_cor_res_os$gene1gene2_continuous_rnr_roc_p,method="bonferroni")

tmp_cph <- coxph(surv ~ mean_highos_genes_z + mean_lowos_genes_z,data=clin[!is.na(clin$surv),])
tmp_cpb_ggforest_mean_highos_z_mean_lowos_z <- ggforest(tmp_cph,data=clin[!is.na(clin$surv),],fontsize=0.9)
ggsave("figure4_outputs/S21a.CPB_mean_highos_z_mean_lowos_z_forestplot.pdf",tmp_cpb_ggforest_mean_highos_z_mean_lowos_z,height=3,width=9)
knitr::include_graphics("figure4_outputs/S21a.CPB_mean_highos_z_mean_lowos_z_forestplot.pdf")
```

<embed src="figure4_outputs/S21a.CPB_mean_highos_z_mean_lowos_z_forestplot.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
mod_meanrz_meannrz <- glm(as.numeric(as.factor(clin$Response))-1 ~ clin$mean_r_genes_z + clin$mean_nr_genes_z,family=binomial(link="logit"))
mod_meanhighosz_meanlowosz <- glm(as.numeric(as.factor(clin$Response))-1 ~ clin$mean_highos_genes_z + clin$mean_lowos_genes_z,family=binomial(link="logit"))

rocobj_meanr_z_meannr_z <- roc.area(as.numeric(as.factor(clin$Response))-1,mod_meanrz_meannrz$fitted.values)
rocobj_meanhighos_z_meanlowos_z <- roc.area(as.numeric(as.factor(clin$Response))-1,mod_meanhighosz_meanlowosz$fitted.values)

rocobj_cpb_meanr_z_meannr_z_string <- paste0("AUC: ",formatC(rocobj_meanr_z_meannr_z$A,digits=3,format="f"),", p-value: ",formatC(rocobj_meanr_z_meannr_z$p.value,digits=2,format="e"))
rocobj_cpb_meanhighos_z_meanlowos_z_string <- paste0("AUC: ",formatC(rocobj_meanhighos_z_meanlowos_z$A,digits=3,format="f"),", p-value: ",formatC(rocobj_meanhighos_z_meanlowos_z$p.value,digits=2,
format="e"))

pdf("figure4_outputs/S21b.CPB_aggregatemodels_zscores_ROC_withpvalues.pdf")
pROC::plot.roc(clin$Response,mod_meanrz_meannrz$fitted.values,col="blue",title="CPB cohort response classification")
pROC::lines.roc(clin$Response,mod_meanhighosz_meanlowosz$fitted.values,col="black")
legend("bottomright", legend=c(paste0("mean(R genes Z-score)+mean(NR genes Z-score): ",rocobj_cpb_meanr_z_meannr_z_string),paste0("mean(high OS genes Z-score)+mean(low OS genes Z-score): ",
rocobj_cpb_meanhighos_z_meanlowos_z_string)), col=c("blue","black"), lwd=2,cex=0.7)
dev.off()
```

    ## pdf 
    ##   2

``` {.r}
knitr::include_graphics("figure4_outputs/S21b.CPB_aggregatemodels_zscores_ROC_withpvalues.pdf")
```

<embed src="figure4_outputs/S21b.CPB_aggregatemodels_zscores_ROC_withpvalues.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
g_cindex_roc_os <- ggplot(half_m_tmp_cor_res_os,aes(x=gene1gene2_continuous_cindex,y=gene1gene2_continuous_rnr_roc_auc)) +
geom_point(aes(color=model_type)) + theme_classic() + ylab("CPB Cohort R/NR ROC AUC") + xlab("CPB Cohort Cox model C-index") + geom_abline(slope=1,intercept=0) +
theme(legend.position="bottom",text=element_text(size=16),legend.text=element_text(size=9),legend.title=element_text(size=9)) +
geom_text_repel(data=half_m_tmp_cor_res_os[half_m_tmp_cor_res_os$gene1gene2_continuous_logrank_bonferroni_p<0.05&half_m_tmp_cor_res_os$gene1gene2_continuous_rnr_roc_bonferroni_p<0.05,],
aes(label=gene_pair),size=4)
ggsave("figure4_outputs/S21c.OSgenes_CPB_ROC_vs_Cindex.pdf",g_cindex_roc_os,height=7,width=7)
knitr::include_graphics("figure4_outputs/S21c.OSgenes_CPB_ROC_vs_Cindex.pdf")
```

<embed src="figure4_outputs/S21c.OSgenes_CPB_ROC_vs_Cindex.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
g_cpb_cindex_by_model_type_os <- ggplot(half_m_tmp_cor_res_os,aes(x=model_type,y= gene1gene2_continuous_cindex)) +
geom_boxplot(aes(fill=model_type),outlier.shape=NA) + geom_jitter(alpha=1/10) +
geom_signif(comparisons=list(c("Low OS gene/Low OS gene","High OS gene/High OS gene"),c("High OS gene/High OS gene","High OS gene/Low OS gene"),
c("Low OS gene/Low OS gene","High OS gene/Low OS gene")), step_increase = .05) +
theme_classic(base_size=16) + theme(legend.position="none",axis.text.x=element_text(size=rel(0.8),colour="black"),axis.text.y=element_text(colour="black")) + xlab("") +
ylab("Primary cohort model C-index")
ggsave("figure4_outputs/S21d.OSgene_CPB_cindex_comparison.pdf",g_cpb_cindex_by_model_type_os,height=7,width=7)
knitr::include_graphics("figure4_outputs/S21d.OSgene_CPB_cindex_comparison.pdf")
```

<embed src="figure4_outputs/S21d.OSgene_CPB_cindex_comparison.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
g_cpb_auc_by_model_type_os <- ggplot(half_m_tmp_cor_res_os,aes(x=model_type,y=gene1gene2_continuous_rnr_roc_auc)) +
geom_boxplot(aes(fill=model_type),outlier.shape=NA) + geom_jitter(alpha=1/20) +
geom_signif(comparisons=list(c("Low OS gene/Low OS gene","High OS gene/High OS gene"),c("High OS gene/High OS gene","High OS gene/Low OS gene"),c("Low OS gene/Low OS gene",
"High OS gene/Low OS gene")), step_increase = .05) +
theme_classic(base_size=16) + theme(legend.position="none",axis.text.x=element_text(colour="black",size=10),axis.text.y=element_text(colour="black")) +
xlab("") + ylab("Primary cohort model ROC AUC")
ggsave("figure4_outputs/S21e.OSgene_CPB_roc_auc_comparison.pdf",g_cpb_auc_by_model_type_os,height=7,width=7)
knitr::include_graphics("figure4_outputs/S21e.OSgene_CPB_roc_auc_comparison.pdf")
```

<embed src="figure4_outputs/S21e.OSgene_CPB_roc_auc_comparison.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
g_best_models_comparison_rnr <- ggplot(half_m_tmp_cor_res,aes(x=-1*log10(gene1gene2_continuous_logrank_p),y=-1*log10(gene1gene2_continuous_rnr_roc_p))) + geom_point(aes(color=model_type)) +
#geom_text_repel(aes(label=model_name,color=model_type),size=2) + scale_color_manual(values=c("black","red"))
theme_classic() + xlab("CPB Cohort -log10(Cox log-rank p)") + ylab("CPB Cohort -log10(Response AUC p)") + theme(legend.position="bottom",text=element_text(size=16)) +
geom_vline(xintercept=-1*log10(0.05/nrow(half_m_tmp_cor_res)),color="red") + geom_hline(yintercept=-1*log10(0.05/nrow(half_m_tmp_cor_res)),color="red")
ggsave("figure4_outputs/S21f.RNRgenes_CPB_pairwise_models.pdf",g_best_models_comparison_rnr,height=7,width=7)
knitr::include_graphics("figure4_outputs/S21f.RNRgenes_CPB_pairwise_models.pdf")
```

<embed src="figure4_outputs/S21f.RNRgenes_CPB_pairwise_models.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
g_cindex_roc_rnr <- ggplot(half_m_tmp_cor_res,aes(x=gene1gene2_continuous_cindex,y=gene1gene2_continuous_rnr_roc_auc)) +
geom_point(aes(color=model_type)) + theme_classic() + ylab("CPB Cohort R/NR ROC AUC") + xlab("CPB Cohort Cox model C-index") + geom_abline(slope=1,intercept=0) +
theme(legend.position="bottom",text=element_text(size=16),legend.text=element_text(size=9),legend.title=element_text(size=9))
ggsave("figure4_outputs/S21g.RNRgenes_CPB_ROC_vs_Cindex.pdf",g_cindex_roc_rnr,height=7,width=7)
knitr::include_graphics("figure4_outputs/S21g.RNRgenes_CPB_ROC_vs_Cindex.pdf")
```

<embed src="figure4_outputs/S21g.RNRgenes_CPB_ROC_vs_Cindex.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
g_cpb_cindex_by_model_type_rnr <- ggplot(half_m_tmp_cor_res,aes(x=model_type,y=gene1gene2_continuous_cindex)) +
geom_boxplot(aes(fill=model_type),outlier.shape=NA) + geom_jitter(alpha=1/20) +
geom_signif(comparisons = list(c("NR gene/NR gene","R gene/R gene"),c("R gene/R gene","R gene/NR gene"),c("NR gene/NR gene","R gene/NR gene")), step_increase = .05) +
theme_classic(base_size=16) + theme(legend.position="none",axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black")) +
xlab("") + ylab("Primary cohort model C-index")
ggsave("figure4_outputs/S21h.RNRgene_CPB_cindex_comparison.pdf",g_cpb_cindex_by_model_type_rnr,height=7,width=7)
knitr::include_graphics("figure4_outputs/S21h.RNRgene_CPB_cindex_comparison.pdf")
```

<embed src="figure4_outputs/S21h.RNRgene_CPB_cindex_comparison.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
g_cpb_auc_by_model_type_rnr <- ggplot(half_m_tmp_cor_res,aes(x=model_type,y=gene1gene2_continuous_rnr_roc_auc)) +
geom_boxplot(aes(fill=model_type),outlier.shape=NA) + geom_jitter(alpha=1/20) +
geom_signif(comparisons = list(c("NR gene/NR gene","R gene/R gene"),c("R gene/R gene","R gene/NR gene"),c("NR gene/NR gene","R gene/NR gene")), step_increase = .05) +
theme_classic(base_size=16) + theme(legend.position="none",axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black")) +
xlab("") + ylab("Primary cohort model ROC AUC")
ggsave("figure4_outputs/S21i.RNRgene_CPB_roc_auc_comparison.pdf",g_cpb_auc_by_model_type_rnr,height=7,width=7)
knitr::include_graphics("figure4_outputs/S21i.RNRgene_CPB_roc_auc_comparison.pdf")
```

<embed src="figure4_outputs/S21i.RNRgene_CPB_roc_auc_comparison.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
col_rbp <- c("#D55E00","#0072B2", "#CC79A7")

half_m_tmp_cor_res_os$gene_pair_name <- paste0(half_m_tmp_cor_res_os$Var1,"+",half_m_tmp_cor_res_os$Var2)
half_m_tmp_cor_res_os$gene_pair_name[nrow(half_m_tmp_cor_res_os)] <- "mean(High OS genes Z)+mean(Low OS gene Z)"

g_best_models_comparison_os_withsiglabels_big <- ggplot(half_m_tmp_cor_res_os,aes(x=-1*log10(gene1gene2_continuous_logrank_p),y=-1*log10(gene1gene2_continuous_rnr_roc_p))) +
geom_point(aes(color=model_type)) + theme_classic() + xlab("Primary Cohort -log10(Cox log-rank p)") + ylab("Primary Cohort -log10(Response AUC p)") +
theme(legend.position="bottom",text=element_text(size=16),legend.text=element_text(size=9),legend.title=element_text(size=9)) +
geom_vline(xintercept=-1*log10(0.05/nrow(half_m_tmp_cor_res_os)),color="black",linetype="dashed") +
geom_hline(yintercept=-1*log10(0.05/nrow(half_m_tmp_cor_res_os)),color="black",linetype="dashed") +
geom_text_repel(data=half_m_tmp_cor_res_os[half_m_tmp_cor_res_os$gene1gene2_continuous_logrank_bonferroni_p<0.05&half_m_tmp_cor_res_os$gene1gene2_continuous_rnr_roc_bonferroni_p<0.05,],
aes(label=gene_pair_name),size=5,min.segment.length=unit(0, "lines"),force=2) + scale_color_manual(values=col_rbp,name="Gene pair")
ggsave("figure4_outputs/4a.OSgenes_CPB_pairwise_models_withsiglabels_big.pdf",g_best_models_comparison_os_withsiglabels_big,height=7,width=7)
knitr::include_graphics("figure4_outputs/4a.OSgenes_CPB_pairwise_models_withsiglabels_big.pdf")
```

<embed src="figure4_outputs/4a.OSgenes_CPB_pairwise_models_withsiglabels_big.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
clin$log10_tmb <- log10(clin$nonsilent_snpdnpindel)

tmp_cph <- coxph(surv ~ MAP4K1 + AGER + log10_tmb,data=clin[!is.na(clin$surv)&!is.na(clin$log10_tmb),])
tmp_ggforest_MAP4K1_AGER_log10tmb <- ggforest(tmp_cph,data=clin[!is.na(clin$surv)&!is.na(clin$log10_tmb),],fontsize=0.9)
ggsave("figure4_outputs/S22a.i.CPB_MAP4K1_AGER_log10tmb_forestplot.pdf", tmp_ggforest_MAP4K1_AGER_log10tmb,height=3,width=7)
knitr::include_graphics("figure4_outputs/S22a.i.CPB_MAP4K1_AGER_log10tmb_forestplot.pdf")
```

<embed src="figure4_outputs/S22a.i.CPB_MAP4K1_AGER_log10tmb_forestplot.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
tmp_cph2 <- coxph(surv ~ MAP4K1 + AGER,data=clin[!is.na(clin$surv)&!is.na(clin$log10_tmb),])
print(anova(tmp_cph,tmp_cph2))
```

    ## Analysis of Deviance Table
    ##  Cox model: response is  surv
    ##  Model 1: ~ MAP4K1 + AGER + log10_tmb
    ##  Model 2: ~ MAP4K1 + AGER
    ##    loglik  Chisq Df P(>|Chi|)  
    ## 1 -118.02                      
    ## 2 -120.36 4.6855  1   0.03042 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` {.r}
os_map4k1_ager_plustmb <- anova(tmp_cph,tmp_cph2)[2,4]

rev_rocobj_cpb_glm_MAP4K1AGER <- glm(as.numeric(as.factor(clin[!is.na(clin$log10_tmb),]$Response))-1 ~ clin[!is.na(clin$log10_tmb),]$MAP4K1 + clin[!is.na(clin$log10_tmb),]$AGER,
family=binomial(link="logit"))
rev_rocobj_cpb_MAP4K1AGER <- roc.area(as.numeric(as.factor(clin[!is.na(clin$log10_tmb),]$Response))-1,rev_rocobj_cpb_glm_MAP4K1AGER$fitted.values)
rev_rocobj_cpb_MAP4K1AGER_p <- formatC(rev_rocobj_cpb_MAP4K1AGER$p.value,digits=2,format="e")
rev_rocobj_cpb_MAP4K1AGER_auc <- formatC(rev_rocobj_cpb_MAP4K1AGER$A,digits=3,format="f")
rev_rocobj_cpb_MAP4K1AGER_string <- paste0("AUC: ",formatC(rev_rocobj_cpb_MAP4K1AGER$A,digits=3,format="f"),", p-value: ",formatC(rev_rocobj_cpb_MAP4K1AGER$p.value,digits=2,format="e"))

rev_rocobj_cpb_glm_MAP4K1AGERlog10tmb <- glm(as.numeric(as.factor(clin[!is.na(clin$log10_tmb),]$Response))-1 ~ clin[!is.na(clin$log10_tmb),]$MAP4K1 + clin[!is.na(clin$log10_tmb),]$AGER +
clin[!is.na(clin$log10_tmb),]$log10_tmb,family=binomial(link="logit"))
rev_rocobj_cpb_MAP4K1AGERlog10tmb <- roc.area(as.numeric(as.factor(clin[!is.na(clin$log10_tmb),]$Response))-1,rev_rocobj_cpb_glm_MAP4K1AGERlog10tmb$fitted.values)
rev_rocobj_cpb_MAP4K1AGERlog10tmb_p <- formatC(rev_rocobj_cpb_MAP4K1AGERlog10tmb$p.value,digits=2,format="e")
rev_rocobj_cpb_MAP4K1AGERlog10tmb_auc <- formatC(rev_rocobj_cpb_MAP4K1AGERlog10tmb$A,digits=3,format="f")
rev_rocobj_cpb_MAP4K1AGERlog10tmb_string <- paste0("AUC: ",formatC(rev_rocobj_cpb_MAP4K1AGERlog10tmb$A,digits=3,format="f"),", p-value: ",formatC(rev_rocobj_cpb_MAP4K1AGERlog10tmb$p.value,digits=2,
format="e"))

print(roc.test(roc(as.numeric(as.factor(clin[!is.na(clin$log10_tmb),]$Response))-1,rev_rocobj_cpb_glm_MAP4K1AGER$fitted.values),
roc(as.numeric(as.factor(clin[!is.na(clin$log10_tmb),]$Response))-1,rev_rocobj_cpb_glm_MAP4K1AGERlog10tmb$fitted.values)))
```

    ## 
    ##  DeLong's test for two correlated ROC curves
    ## 
    ## data:  roc(as.numeric(as.factor(clin[!is.na(clin$log10_tmb), ]$Response)) -  and roc(as.numeric(as.factor(clin[!is.na(clin$log10_tmb), ]$Response)) -     1, rev_rocobj_cpb_glm_MAP4K1AGER$fitted.values) and     1, rev_rocobj_cpb_glm_MAP4K1AGERlog10tmb$fitted.values)
    ## Z = -1.5331, p-value = 0.1252
    ## alternative hypothesis: true difference in AUC is not equal to 0
    ## sample estimates:
    ## AUC of roc1 AUC of roc2 
    ##   0.6177945   0.7080201

``` {.r}
rnr_map4k1_ager_plustmb <- roc.test(roc(as.numeric(as.factor(clin[!is.na(clin$log10_tmb),]$Response))-1,rev_rocobj_cpb_glm_MAP4K1AGER$fitted.values),
roc(as.numeric(as.factor(clin[!is.na(clin$log10_tmb),]$Response))-1,rev_rocobj_cpb_glm_MAP4K1AGERlog10tmb$fitted.values))$p.value


pdf("figure4_outputs/S22a.ii.MAP4K1_AGER_log10tmb_ROC.pdf")
pROC::plot.roc(clin[!is.na(clin$log10_tmb),]$Response,rev_rocobj_cpb_glm_MAP4K1AGER$fitted.values,col="blue")
pROC::lines.roc(clin[!is.na(clin$log10_tmb),]$Response,rev_rocobj_cpb_glm_MAP4K1AGERlog10tmb$fitted.values,col="red")
legend("bottomright", legend=c(paste0("MAP4K1+AGER ",rev_rocobj_cpb_MAP4K1AGER_string),paste0("MAP4K1+AGER+log10(TMB) ",rev_rocobj_cpb_MAP4K1AGERlog10tmb_string)),col=c("blue","red"),lwd=2,cex=0.71)
dev.off()
```

    ## pdf 
    ##   2

``` {.r}
knitr::include_graphics("figure4_outputs/S22a.ii.MAP4K1_AGER_log10tmb_ROC.pdf")
```

<embed src="figure4_outputs/S22a.ii.MAP4K1_AGER_log10tmb_ROC.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
tmp_cph <- coxph(surv ~ MAP4K1 + TBX3 + log10_tmb,data=clin[!is.na(clin$surv)&!is.na(clin$log10_tmb),])
tmp_ggforest_MAP4K1_TBX3_log10tmb <- ggforest(tmp_cph,data=clin[!is.na(clin$surv)&!is.na(clin$log10_tmb),],fontsize=0.9)
ggsave("figure4_outputs/S22b.i.MAP4K1_TBX3_log10tmb_forestplot.pdf", tmp_ggforest_MAP4K1_TBX3_log10tmb,height=3,width=7)
knitr::include_graphics("figure4_outputs/S22b.i.MAP4K1_TBX3_log10tmb_forestplot.pdf")
```

<embed src="figure4_outputs/S22b.i.MAP4K1_TBX3_log10tmb_forestplot.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
tmp_cph2 <- coxph(surv ~ MAP4K1 + TBX3,data=clin[!is.na(clin$surv)&!is.na(clin$log10_tmb),])
print(anova(tmp_cph,tmp_cph2))
```

    ## Analysis of Deviance Table
    ##  Cox model: response is  surv
    ##  Model 1: ~ MAP4K1 + TBX3 + log10_tmb
    ##  Model 2: ~ MAP4K1 + TBX3
    ##    loglik  Chisq Df P(>|Chi|)  
    ## 1 -115.81                      
    ## 2 -118.54 5.4607  1   0.01945 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` {.r}
os_map4k1_tbx3_plustmb <- anova(tmp_cph,tmp_cph2)[2,4]


rev_rocobj_cpb_glm_MAP4K1TBX3 <- glm(as.numeric(as.factor(clin[!is.na(clin$log10_tmb),]$Response))-1 ~ clin[!is.na(clin$log10_tmb),]$MAP4K1 + clin[!is.na(clin$log10_tmb),]$TBX3,
family=binomial(link="logit"))
rev_rocobj_cpb_MAP4K1TBX3 <- roc.area(as.numeric(as.factor(clin[!is.na(clin$log10_tmb),]$Response))-1,rev_rocobj_cpb_glm_MAP4K1TBX3$fitted.values)
rev_rocobj_cpb_MAP4K1TBX3_p <- formatC(rev_rocobj_cpb_MAP4K1TBX3$p.value,digits=2,format="e")
rev_rocobj_cpb_MAP4K1TBX3_auc <- formatC(rev_rocobj_cpb_MAP4K1TBX3$A,digits=3,format="f")
rev_rocobj_cpb_MAP4K1TBX3_string <- paste0("AUC: ",formatC(rev_rocobj_cpb_MAP4K1TBX3$A,digits=3,format="f"),", p-value: ",formatC(rev_rocobj_cpb_MAP4K1TBX3$p.value,digits=2,format="e"))

rev_rocobj_cpb_glm_MAP4K1TBX3log10tmb <- glm(as.numeric(as.factor(clin[!is.na(clin$log10_tmb),]$Response))-1 ~ clin[!is.na(clin$log10_tmb),]$MAP4K1 + clin[!is.na(clin$log10_tmb),]$TBX3 +
clin[!is.na(clin$log10_tmb),]$log10_tmb,family=binomial(link="logit"))
rev_rocobj_cpb_MAP4K1TBX3log10tmb <- roc.area(as.numeric(as.factor(clin[!is.na(clin$log10_tmb),]$Response))-1,rev_rocobj_cpb_glm_MAP4K1TBX3log10tmb$fitted.values)
rev_rocobj_cpb_MAP4K1TBX3log10tmb_p <- formatC(rev_rocobj_cpb_MAP4K1TBX3log10tmb$p.value,digits=2,format="e")
rev_rocobj_cpb_MAP4K1TBX3log10tmb_auc <- formatC(rev_rocobj_cpb_MAP4K1TBX3log10tmb$A,digits=3,format="f")
rev_rocobj_cpb_MAP4K1TBX3log10tmb_string <- paste0("AUC: ",formatC(rev_rocobj_cpb_MAP4K1TBX3log10tmb$A,digits=3,format="f"),", p-value: ",formatC(rev_rocobj_cpb_MAP4K1TBX3log10tmb$p.value,digits=2,
format="e"))

print(roc.test(roc(as.numeric(as.factor(clin[!is.na(clin$log10_tmb),]$Response))-1,rev_rocobj_cpb_glm_MAP4K1TBX3$fitted.values),
roc(as.numeric(as.factor(clin[!is.na(clin$log10_tmb),]$Response))-1,rev_rocobj_cpb_glm_MAP4K1TBX3log10tmb$fitted.values)))
```

    ## 
    ##  DeLong's test for two correlated ROC curves
    ## 
    ## data:  roc(as.numeric(as.factor(clin[!is.na(clin$log10_tmb), ]$Response)) -  and roc(as.numeric(as.factor(clin[!is.na(clin$log10_tmb), ]$Response)) -     1, rev_rocobj_cpb_glm_MAP4K1TBX3$fitted.values) and     1, rev_rocobj_cpb_glm_MAP4K1TBX3log10tmb$fitted.values)
    ## Z = -1.2187, p-value = 0.223
    ## alternative hypothesis: true difference in AUC is not equal to 0
    ## sample estimates:
    ## AUC of roc1 AUC of roc2 
    ##   0.7005013   0.7606516

``` {.r}
rnr_map4k1_tbx3_plustmb <- roc.test(roc(as.numeric(as.factor(clin[!is.na(clin$log10_tmb),]$Response))-1,rev_rocobj_cpb_glm_MAP4K1TBX3$fitted.values),
roc(as.numeric(as.factor(clin[!is.na(clin$log10_tmb),]$Response))-1,rev_rocobj_cpb_glm_MAP4K1TBX3log10tmb$fitted.values))$p.value


pdf("figure4_outputs/S22b.ii.MAP4K1_TBX3_log10tmb_ROC.pdf")
pROC::plot.roc(clin[!is.na(clin$log10_tmb),]$Response,rev_rocobj_cpb_glm_MAP4K1TBX3$fitted.values,col="blue")
pROC::lines.roc(clin[!is.na(clin$log10_tmb),]$Response,rev_rocobj_cpb_glm_MAP4K1TBX3log10tmb$fitted.values,col="red")
legend("bottomright", legend=c(paste0("MAP4K1+TBX3 ",rev_rocobj_cpb_MAP4K1TBX3_string),paste0("MAP4K1+TBX3+log10(TMB) ",rev_rocobj_cpb_MAP4K1TBX3log10tmb_string)),col=c("blue","red"),lwd=2,cex=0.71)
dev.off()
```

    ## pdf 
    ##   2

``` {.r}
knitr::include_graphics("figure4_outputs/S22b.ii.MAP4K1_TBX3_log10tmb_ROC.pdf")
```

<embed src="figure4_outputs/S22b.ii.MAP4K1_TBX3_log10tmb_ROC.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
tmp_cph <- coxph(surv ~ mean_highos_genes_z + mean_lowos_genes_z + log10_tmb,data=clin[!is.na(clin$surv)&!is.na(clin$log10_tmb),])
tmp_ggforest_highos_lowos_log10tmb <- ggforest(tmp_cph,data=clin[!is.na(clin$surv)&!is.na(clin$log10_tmb),],fontsize=0.9)
ggsave("figure4_outputs/S22c.i.highos_lowos_log10tmb_forestplot.pdf", tmp_ggforest_highos_lowos_log10tmb,height=3,width=7)
knitr::include_graphics("figure4_outputs/S22c.i.highos_lowos_log10tmb_forestplot.pdf")
```

<embed src="figure4_outputs/S22c.i.highos_lowos_log10tmb_forestplot.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
tmp_cph2 <- coxph(surv ~ mean_highos_genes_z + mean_lowos_genes_z,data=clin[!is.na(clin$surv)&!is.na(clin$log10_tmb),])
print(anova(tmp_cph,tmp_cph2))
```

    ## Analysis of Deviance Table
    ##  Cox model: response is  surv
    ##  Model 1: ~ mean_highos_genes_z + mean_lowos_genes_z + log10_tmb
    ##  Model 2: ~ mean_highos_genes_z + mean_lowos_genes_z
    ##    loglik  Chisq Df P(>|Chi|)
    ## 1 -112.13                    
    ## 2 -113.26 2.2561  1    0.1331

``` {.r}
os_meanhighos_meanlowos_plustmb <- anova(tmp_cph,tmp_cph2)[2,4]

rev_rocobj_cpb_glm_highoslowos <- glm(as.numeric(as.factor(clin[!is.na(clin$log10_tmb),]$Response))-1 ~ clin[!is.na(clin$log10_tmb),]$mean_highos_genes_z + clin[!is.na(clin$log10_tmb),]$mean_lowos_genes_z,
family=binomial(link="logit"))
rev_rocobj_cpb_highoslowos <- roc.area(as.numeric(as.factor(clin[!is.na(clin$log10_tmb),]$Response))-1,rev_rocobj_cpb_glm_highoslowos$fitted.values)
rev_rocobj_cpb_highoslowos_p <- formatC(rev_rocobj_cpb_highoslowos$p.value,digits=2,format="e")
rev_rocobj_cpb_highoslowos_auc <- formatC(rev_rocobj_cpb_highoslowos$A,digits=3,format="f")
rev_rocobj_cpb_highoslowos_string <- paste0("AUC: ",formatC(rev_rocobj_cpb_highoslowos$A,digits=3,format="f"),", p-value: ",formatC(rev_rocobj_cpb_highoslowos$p.value,digits=2,format="e"))

rev_rocobj_cpb_glm_highoslowoslog10tmb <- glm(as.numeric(as.factor(clin[!is.na(clin$log10_tmb),]$Response))-1 ~ clin[!is.na(clin$log10_tmb),]$mean_highos_genes_z +
clin[!is.na(clin$log10_tmb),]$mean_lowos_genes_z + clin[!is.na(clin$log10_tmb),]$log10_tmb,family=binomial(link="logit"))
rev_rocobj_cpb_highoslowoslog10tmb <- roc.area(as.numeric(as.factor(clin[!is.na(clin$log10_tmb),]$Response))-1,rev_rocobj_cpb_glm_highoslowoslog10tmb$fitted.values)
rev_rocobj_cpb_highoslowoslog10tmb_p <- formatC(rev_rocobj_cpb_highoslowoslog10tmb$p.value,digits=2,format="e")
rev_rocobj_cpb_highoslowoslog10tmb_auc <- formatC(rev_rocobj_cpb_highoslowoslog10tmb$A,digits=3,format="f")
rev_rocobj_cpb_highoslowoslog10tmb_string <- paste0("AUC: ",formatC(rev_rocobj_cpb_highoslowoslog10tmb$A,digits=3,format="f"),", p-value: ",formatC(rev_rocobj_cpb_highoslowoslog10tmb$p.value,digits=2,
format="e"))

print(roc.test(roc(as.numeric(as.factor(clin[!is.na(clin$log10_tmb),]$Response))-1,rev_rocobj_cpb_glm_highoslowos$fitted.values),
roc(as.numeric(as.factor(clin[!is.na(clin$log10_tmb),]$Response))-1,rev_rocobj_cpb_glm_highoslowoslog10tmb$fitted.values)))
```

    ## 
    ##  DeLong's test for two correlated ROC curves
    ## 
    ## data:  roc(as.numeric(as.factor(clin[!is.na(clin$log10_tmb), ]$Response)) -  and roc(as.numeric(as.factor(clin[!is.na(clin$log10_tmb), ]$Response)) -     1, rev_rocobj_cpb_glm_highoslowos$fitted.values) and     1, rev_rocobj_cpb_glm_highoslowoslog10tmb$fitted.values)
    ## Z = -0.15234, p-value = 0.8789
    ## alternative hypothesis: true difference in AUC is not equal to 0
    ## sample estimates:
    ## AUC of roc1 AUC of roc2 
    ##   0.7869674   0.7919799

``` {.r}
rnr_meanhighos_meanlowos_plustmb <- roc.test(roc(as.numeric(as.factor(clin[!is.na(clin$log10_tmb),]$Response))-1,rev_rocobj_cpb_glm_highoslowos$fitted.values),
roc(as.numeric(as.factor(clin[!is.na(clin$log10_tmb),]$Response))-1,rev_rocobj_cpb_glm_highoslowoslog10tmb$fitted.values))$p.value

pdf("figure4_outputs/S22c.ii.highos_lowos_log10tmb_ROC.pdf")
pROC::plot.roc(clin[!is.na(clin$log10_tmb),]$Response,rev_rocobj_cpb_glm_highoslowos$fitted.values,col="blue")
pROC::lines.roc(clin[!is.na(clin$log10_tmb),]$Response,rev_rocobj_cpb_glm_highoslowoslog10tmb$fitted.values,col="red")
legend("bottomright", legend=c(paste0("Long OS metagene + Short OS metagene + ",rev_rocobj_cpb_highoslowos_string),paste0("MAP4K1+AGER+log10(TMB) ",rev_rocobj_cpb_highoslowoslog10tmb_string)),
col=c("blue","red"),lwd=2,cex=0.71)
dev.off()
```

    ## pdf 
    ##   2

``` {.r}
knitr::include_graphics("figure4_outputs/S22c.ii.highos_lowos_log10tmb_ROC.pdf")
```

<embed src="figure4_outputs/S22c.ii.highos_lowos_log10tmb_ROC.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
df_tmb_pvals <- data.frame(model=c("MAP4K1&AGER","MAP4K1&TBX3","Long OS metagene & Short OS metagene"),os_pval=c(os_map4k1_ager_plustmb,os_map4k1_tbx3_plustmb,os_meanhighos_meanlowos_plustmb),
rnr_pval=c(rnr_map4k1_ager_plustmb,rnr_map4k1_tbx3_plustmb,rnr_meanhighos_meanlowos_plustmb),stringsAsFactors=F)
df_tmb_pvals$os_pval <- formatC(df_tmb_pvals$os_pval,width=2,digits=2,format="e")
df_tmb_pvals$rnr_pval <- formatC(df_tmb_pvals$rnr_pval,width=2,digits=2,format="e")
colnames(df_tmb_pvals) <- c("model","Likelihood Ratio Test p value\nfor OS model with\n vs. without log10(TMB)","DeLong's Test p value\nfor Response model with\nvs. without log10(TMB)")
pdf("figure4_outputs/S22d.models_with_vs_without_TMB.pdf",height=7,width=8)
grid.table(df_tmb_pvals, rows=NULL)
dev.off()
```

    ## pdf 
    ##   2

``` {.r}
knitr::include_graphics("figure4_outputs/S22d.models_with_vs_without_TMB.pdf")
```

<embed src="figure4_outputs/S22d.models_with_vs_without_TMB.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
m_hpa_map4k1 <- reshape2::melt(hpa_tab[strsplit(rownames(val_df_res_os_sig[val_df_res_os_sig$gene=="MAP4K1",]),"\\.")[[1]][1],])
m_hpa_map4k1$cell_type <- unlist(lapply(strsplit(as.character(m_hpa_map4k1$variable),"\\.\\.\\."),"[",1))
dt_m_hpa_map4k1 <- data.table(m_hpa_map4k1)
m_hpa_map4k1$cell_type <- factor(as.character(m_hpa_map4k1$cell_type),levels=rev(arrange(dt_m_hpa_map4k1[,median(value),by=cell_type],V1)$cell_type))
g_hpa_map4k1 <- ggplot(m_hpa_map4k1,aes(x=cell_type,y=log2(value+1))) + geom_boxplot(aes(fill=cell_type),outlier.shape=NA) + geom_jitter() + theme_classic(base_size=16) +
theme(legend.position="none",axis.text.x=element_text(angle=90,hjust=0.95)) + ylab("log2(MAP4K1 TPM + 1)")
ggsave("figure4_outputs/S23a.HPA_MAP4K1.pdf",g_hpa_map4k1,height=7,width=12)
knitr::include_graphics("figure4_outputs/S23a.HPA_MAP4K1.pdf")
```

<embed src="figure4_outputs/S23a.HPA_MAP4K1.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
cbpalette_reordered <- c("#E69F00","#999999","#56B4E9","#009E73","#F0E442")

kt <- kruskal.test(clin$TBX3,clin$plot_cluster)
kt_p <- formatC(kt$p.value,digits=3,format="e")
g_TBX3_nmf <- ggplot(clin,aes(x=plot_cluster,y=TBX3)) + geom_boxplot(aes(fill=plot_cluster),outlier.shape = NA) + geom_jitter() + theme_classic() + ylab("TBX3") +
scale_fill_manual(values=cbpalette_reordered) + ggtitle(paste0("Kruskal-Wallis p =",kt_p)) + xlab("") + theme(legend.position="none")
ggsave("figure4_outputs/S23b.cpb_nmf_cluster_TBX3.pdf",g_TBX3_nmf)
knitr::include_graphics("figure4_outputs/S23b.cpb_nmf_cluster_TBX3.pdf")
```

<embed src="figure4_outputs/S23b.cpb_nmf_cluster_TBX3.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
kt <- kruskal.test(clin$AGER,clin$plot_cluster)
kt_p <- formatC(kt$p.value,digits=3,format="e")
g_AGER_nmf <- ggplot(clin,aes(x=plot_cluster,y=AGER)) + geom_boxplot(aes(fill=plot_cluster),outlier.shape = NA) + geom_jitter() + theme_classic() + ylab("AGER") +
scale_fill_manual(values=cbpalette_reordered) + ggtitle(paste0("Kruskal-Wallis p =",kt_p)) + xlab("") + theme(legend.position="none")
ggsave("figure4_outputs/S23c.cpb_nmf_cluster_AGER.pdf",g_AGER_nmf)
knitr::include_graphics("figure4_outputs/S23c.cpb_nmf_cluster_AGER.pdf")
```

<embed src="figure4_outputs/S23c.cpb_nmf_cluster_AGER.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
kt <- kruskal.test(clin$mean_lowos_genes_z,clin$plot_cluster)
kt_p <- formatC(kt$p.value,digits=3,format="e")
g_lowosgenes_z_nmf <- ggplot(clin,aes(x=plot_cluster,y=mean_lowos_genes_z)) + geom_boxplot(aes(fill=plot_cluster),outlier.shape = NA) + geom_jitter() + theme_classic() +
ylab("Mean Low OS genes Z-score") + scale_fill_manual(values=cbpalette_reordered) + ggtitle(paste0("Kruskal-Wallis p =",kt_p)) + xlab("") + theme(legend.position="none")
ggsave("figure4_outputs/S23d.cpb_nmf_cluster_lowosmetagene.pdf",g_lowosgenes_z_nmf)
knitr::include_graphics("figure4_outputs/S23d.cpb_nmf_cluster_lowosmetagene.pdf")
```

<embed src="figure4_outputs/S23d.cpb_nmf_cluster_lowosmetagene.pdf" width="600px" height="600px" type="application/pdf" />


``` {.r}
clin$log10_rna_tcb <- log10(clin$rna_tcb)
clin$scale_log10_rna_tcb <- scale(clin$log10_rna_tcb)
clin$log10_rna_bcb <- log10(clin$rna_bcb)
clin$scale_log10_rna_bcb <- scale(clin$log10_rna_bcb)
clin$scale_TIDE_Score <- scale(clin$TIDE_Score)
clin$MHC_II_zscore <- scale(clin$MHC_II)
clin$scale_CYT <- scale(clin$CYT)
clin$scale_GEP <- scale(clin$GEP)
clin$scale_CD274 <- scale(clin$CD274)
clin$scale_IMPRES <- scale(clin$IMPRES)

clin$OS <- ifelse(clin$overall_survival>=365,"OS >= 1yr","OS < 1yr")

rwb <- colorRampPalette(colors = c("red", "white", "blue"))
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}
tmp_rawdata_breaks <- c(quantile_breaks(unlist(clin[,c("scale_log10_rna_tcb","scale_log10_rna_bcb","scale_TIDE_Score","MHC_II_zscore","scale_CYT","scale_GEP","scale_CD274","scale_IMPRES",
"MAP4K1_z","TBX3_z","AGER_z","mean_lowos_genes_z","mean_highos_genes_z")]),10))

raw_heatmap_annotation_colors <- list(
                              OS=c("OS >= 1yr"="blue","OS < 1yr"="red"),
                              Response=c("nonresponder"="red","responder"="blue"),
                              plot_cluster=c("MITF low"="#56B4E9","Immune"="#009E73","Intermediate"="#999999","Keratin high"="#F0E442","MITF high"="#E69F00"),
                              cohort=c("Hugo"="orange","MGH"="Blue","Riaz"="green","VanAllen"="purple")
)

main_heatmap <- pheatmap(t(clin[order(clin$overall_survival),
c("scale_log10_rna_tcb","scale_log10_rna_bcb","scale_TIDE_Score","MHC_II_zscore","scale_CYT","scale_GEP","scale_CD274","scale_IMPRES","MAP4K1_z","TBX3_z","AGER_z","mean_lowos_genes_z",
"mean_highos_genes_z")]),
annotation_col=clin[,c("cohort","Response","OS","plot_cluster")],
cluster_rows=TRUE,cluster_cols=TRUE,clustering_distance_rows="correlation",fontsize_row=5, fontsize_col=3,cellwidth=3,cellheight=5,fontsize=5,na_col="grey",
breaks=tmp_rawdata_breaks,color=rev(brewer.pal(n = 9, name = "RdBu")),treeheight_row=10,treeheight_col=10,annotation_colors=raw_heatmap_annotation_colors,
legend_breaks=c(-4,-3,-2,-1,0,1,2,3,4))

save_pheatmap_pdf <- function(x, filename, width=10, height=7) {
   pdf(filename, width=width, height=height)
   grid::grid.newpage()
   grid::grid.draw(x$gtable)
   dev.off()
}
save_pheatmap_pdf(main_heatmap, "figure4_outputs/4b.OSgenes_rawdata_heatmap_clustered_9col.pdf")
```

    ## pdf 
    ##   2

``` {.r}
knitr::include_graphics("figure4_outputs/4b.OSgenes_rawdata_heatmap_clustered_9col.pdf")
```

<embed src="figure4_outputs/4b.OSgenes_rawdata_heatmap_clustered_9col.pdf" width="700px" height="700px" type="application/pdf" />

Functions

``` {.r}
make_os_p_cindex_plot <- function(clin,pdf_string) {
cph_cohort_MAP4K1TBX3_p <- summary(coxph(surv ~ MAP4K1 + TBX3,data=clin))$sctest[3]
cph_cohort_MAP4K1TBX3_cindex <- summary(coxph(surv ~ MAP4K1 + TBX3,data=clin))$concordance[1]

cph_cohort_MAP4K1AGER_p <- summary(coxph(surv ~ MAP4K1 + AGER,data=clin))$sctest[3]
cph_cohort_MAP4K1AGER_cindex <- summary(coxph(surv ~ MAP4K1 + AGER,data=clin))$concordance[1]

cph_cohort_highosgeneszlowosgenesz_p <- summary(coxph(surv ~ mean_highos_genes_z + mean_lowos_genes_z,data=clin))$sctest[3]
cph_cohort_highosgeneszlowosgenesz_cindex <- summary(coxph(surv ~ mean_highos_genes_z + mean_lowos_genes_z,data=clin))$concordance[1]

cph_cohort_CD274_p <- summary(coxph(surv ~ CD274,data=clin))$sctest[3]
cph_cohort_CD274_cindex <- summary(coxph(surv ~ CD274,data=clin))$concordance[1]
cph_cohort_CYT_p <- summary(coxph(surv ~ CYT,data=clin))$sctest[3]
cph_cohort_CYT_cindex <- summary(coxph(surv ~ CYT,data=clin))$concordance[1]
cph_cohort_GEP_p <- summary(coxph(surv ~ GEP,data=clin))$sctest[3]
cph_cohort_GEP_cindex <- summary(coxph(surv ~ GEP,data=clin))$concordance[1]
cph_cohort_IMPRES_p <- summary(coxph(surv ~ IMPRES,data=clin))$sctest[3]
cph_cohort_IMPRES_cindex <- summary(coxph(surv ~ IMPRES,data=clin))$concordance[1]
cph_cohort_TIDE_p <- summary(coxph(surv ~ TIDE_Score,data=clin))$sctest[3]
cph_cohort_TIDE_cindex <- summary(coxph(surv ~ TIDE_Score,data=clin))$concordance[1]
cph_cohort_MHCII_p <- summary(coxph(surv ~ MHC_II,data=clin))$sctest[3]
cph_cohort_MHCII_cindex <- summary(coxph(surv ~ MHC_II,data=clin))$concordance[1]
cph_cohort_TCB_p <- summary(coxph(surv ~ log10(rna_tcb),data=clin))$sctest[3]
cph_cohort_TCB_cindex <- summary(coxph(surv ~ log10(rna_tcb),data=clin))$concordance[1]
cph_cohort_BCB_p <- summary(coxph(surv ~ log10(rna_bcb),data=clin))$sctest[3]
cph_cohort_BCB_cindex <- summary(coxph(surv ~ log10(rna_bcb),data=clin))$concordance[1]

df_model_comparison_cohort <- data.frame(X1=c("CD274","GEP","IMPRES","TIDE","CYT","TCB","BCB","MHC II","MAP4K1+TBX3","MAP4K1+AGER","High OS genes+Low OS genes"),
X2=c(cph_cohort_CD274_p,cph_cohort_GEP_p,cph_cohort_IMPRES_p,cph_cohort_TIDE_p,cph_cohort_CYT_p,cph_cohort_TCB_p,cph_cohort_BCB_p,cph_cohort_MHCII_p,cph_cohort_MAP4K1TBX3_p,cph_cohort_MAP4K1AGER_p,cph_cohort_highosgeneszlowosgenesz_p),
X3=c(cph_cohort_CD274_cindex,cph_cohort_GEP_cindex,cph_cohort_IMPRES_cindex,cph_cohort_TIDE_cindex,cph_cohort_CYT_cindex,cph_cohort_TCB_cindex,cph_cohort_BCB_cindex,cph_cohort_MHCII_cindex,cph_cohort_MAP4K1TBX3_cindex,cph_cohort_MAP4K1AGER_cindex,
cph_cohort_highosgeneszlowosgenesz_cindex),stringsAsFactors=F)
df_model_comparison_cohort$X1 <- factor(df_model_comparison_cohort$X1,levels=c("CD274","GEP","IMPRES","TIDE","CYT","TCB","BCB","MHC II","MAP4K1+TBX3","MAP4K1+AGER","High OS genes+Low OS genes"))
df_model_comparison_cohort$X2 <- as.numeric(df_model_comparison_cohort$X2)
df_model_comparison_cohort$X3 <- as.numeric(df_model_comparison_cohort$X3)
g_survmodel_compare_cohort_withBF3models <- ggplot(df_model_comparison_cohort,aes(x=as.numeric(X3),y=-1*log10(as.numeric(X2)))) + geom_point(aes(color=grepl("\\+",as.character(X1))),size=2) + geom_text_repel(aes(label=X1),size=6) + theme_classic(base_size=16) + geom_hline(yintercept = -1*log10(0.05),linetype = "dashed") + xlab("Cox model C-index") + ylab("-log10(Cox model Log-Rank p-value)") + scale_color_manual(values=c("black","red")) + theme(legend.position = "none")

ggsave(pdf_string, g_survmodel_compare_cohort_withBF3models,height=7,width=7)
return(df_model_comparison_cohort)
}

make_pfs_p_cindex_plot <- function(clin,pdf_string) {
p_cph_cohort_MAP4K1TBX3_p <- summary(coxph(pfs_surv ~ MAP4K1 + TBX3,data=clin))$sctest[3]
p_cph_cohort_MAP4K1TBX3_cindex <- summary(coxph(pfs_surv ~ MAP4K1 + TBX3,data=clin))$concordance[1]

p_cph_cohort_MAP4K1AGER_p <- summary(coxph(pfs_surv ~ MAP4K1 + AGER,data=clin))$sctest[3]
p_cph_cohort_MAP4K1AGER_cindex <- summary(coxph(pfs_surv ~ MAP4K1 + AGER,data=clin))$concordance[1]

p_cph_cohort_highosgeneszlowosgenesz_p <- summary(coxph(pfs_surv ~ mean_highos_genes_z + mean_lowos_genes_z,data=clin))$sctest[3]
p_cph_cohort_highosgeneszlowosgenesz_cindex <- summary(coxph(pfs_surv ~ mean_highos_genes_z + mean_lowos_genes_z,data=clin))$concordance[1]

p_cph_cohort_CD274_p <- summary(coxph(pfs_surv ~ CD274,data=clin))$sctest[3]
p_cph_cohort_CD274_cindex <- summary(coxph(pfs_surv ~ CD274,data=clin))$concordance[1]
p_cph_cohort_CYT_p <- summary(coxph(pfs_surv ~ CYT,data=clin))$sctest[3]
p_cph_cohort_CYT_cindex <- summary(coxph(pfs_surv ~ CYT,data=clin))$concordance[1]
p_cph_cohort_GEP_p <- summary(coxph(pfs_surv ~ GEP,data=clin))$sctest[3]
p_cph_cohort_GEP_cindex <- summary(coxph(pfs_surv ~ GEP,data=clin))$concordance[1]
p_cph_cohort_IMPRES_p <- summary(coxph(pfs_surv ~ IMPRES,data=clin))$sctest[3]
p_cph_cohort_IMPRES_cindex <- summary(coxph(pfs_surv ~ IMPRES,data=clin))$concordance[1]
p_cph_cohort_TIDE_p <- summary(coxph(pfs_surv ~ TIDE_Score,data=clin))$sctest[3]
p_cph_cohort_TIDE_cindex <- summary(coxph(pfs_surv ~ TIDE_Score,data=clin))$concordance[1]
p_cph_cohort_MHCII_p <- summary(coxph(pfs_surv ~ MHC_II,data=clin))$sctest[3]
p_cph_cohort_MHCII_cindex <- summary(coxph(pfs_surv ~ MHC_II,data=clin))$concordance[1]
p_cph_cohort_TCB_p <- summary(coxph(pfs_surv ~ log10(rna_tcb),data=clin))$sctest[3]
p_cph_cohort_TCB_cindex <- summary(coxph(pfs_surv ~ log10(rna_tcb),data=clin))$concordance[1]
p_cph_cohort_BCB_p <- summary(coxph(pfs_surv ~ log10(rna_bcb),data=clin))$sctest[3]
p_cph_cohort_BCB_cindex <- summary(coxph(pfs_surv ~ log10(rna_bcb),data=clin))$concordance[1]

df_model_comparison_cohort_pfs <- data.frame(X1=c("CD274","GEP","IMPRES","TIDE","CYT","TCB","BCB","MHC II","MAP4K1+TBX3","MAP4K1+AGER","High OS genes+Low OS genes"),X2=c(p_cph_cohort_CD274_p,p_cph_cohort_GEP_p,p_cph_cohort_IMPRES_p,p_cph_cohort_TIDE_p,p_cph_cohort_CYT_p,p_cph_cohort_TCB_p,p_cph_cohort_BCB_p,p_cph_cohort_MHCII_p,p_cph_cohort_MAP4K1TBX3_p,p_cph_cohort_MAP4K1AGER_p,p_cph_cohort_highosgeneszlowosgenesz_p),X3=c(p_cph_cohort_CD274_cindex,p_cph_cohort_GEP_cindex,p_cph_cohort_IMPRES_cindex,p_cph_cohort_TIDE_cindex,p_cph_cohort_CYT_cindex,p_cph_cohort_TCB_cindex,p_cph_cohort_BCB_cindex,p_cph_cohort_MHCII_cindex,p_cph_cohort_MAP4K1TBX3_cindex,p_cph_cohort_MAP4K1AGER_cindex,p_cph_cohort_highosgeneszlowosgenesz_cindex),stringsAsFactors=F)
print(df_model_comparison_cohort_pfs)
df_model_comparison_cohort_pfs$X1 <- factor(df_model_comparison_cohort_pfs$X1,levels=c("CD274","GEP","IMPRES","TIDE","CYT","TCB","BCB","MHC II","MAP4K1+TBX3","MAP4K1+AGER","High OS genes+Low OS genes"))
df_model_comparison_cohort_pfs$X2 <- as.numeric(df_model_comparison_cohort_pfs$X2)
df_model_comparison_cohort_pfs$X3 <- as.numeric(df_model_comparison_cohort_pfs$X3)
g_survmodel_compare_cohort_pfs_withBF3models <- ggplot(df_model_comparison_cohort_pfs,aes(x=as.numeric(X3),y=-1*log10(as.numeric(X2)))) + geom_point(aes(color=grepl("\\+",as.character(X1))),size=2) + geom_text_repel(aes(label=X1),size=6) + theme_classic(base_size=16) + geom_hline(yintercept = -1*log10(0.05),linetype = "dashed") + xlab("Cox model C-index") + ylab("-log10(Cox model Log-Rank p-value)") + scale_color_manual(values=c("black","red")) + theme(legend.position = "none")
ggsave(pdf_string, g_survmodel_compare_cohort_pfs_withBF3models,height=7,width=7)

return(df_model_comparison_cohort_pfs)
}

make_rnr_p_auc_plot <- function(clin,pdf_string_1,pdf_string_2) {
rocobj_cohort_CD274 <- roc.area(as.numeric(as.factor(clin$Response))-1,clin$CD274)
rocobj_cohort_CD274_p <- formatC(rocobj_cohort_CD274$p.value,digits=3,format="f")
rocobj_cohort_CD274_auc <- formatC(rocobj_cohort_CD274$A,digits=3,format="f")
rocobj_cohort_CD274_string <- paste0("AUC: ",formatC(rocobj_cohort_CD274$A,digits=3,format="f"),", p-value: ",formatC(rocobj_cohort_CD274$p.value,digits=3,format="f"))
rocobj_cohort_CYT <- roc.area(as.numeric(as.factor(clin$Response))-1,clin$CYT)
rocobj_cohort_CYT_p <- formatC(rocobj_cohort_CYT$p.value,digits=3,format="f")
rocobj_cohort_CYT_auc <- formatC(rocobj_cohort_CYT$A,digits=3,format="f")
rocobj_cohort_CYT_string <- paste0("AUC: ",formatC(rocobj_cohort_CYT$A,digits=3,format="f"),", p-value: ",formatC(rocobj_cohort_CYT$p.value,digits=3,format="f"))
rocobj_cohort_GEP <- roc.area(as.numeric(as.factor(clin$Response))-1,clin$GEP)
rocobj_cohort_GEP_p <- formatC(rocobj_cohort_GEP$p.value,digits=3,format="f")
rocobj_cohort_GEP_auc <- formatC(rocobj_cohort_GEP$A,digits=3,format="f")
rocobj_cohort_GEP_string <- paste0("AUC: ",formatC(rocobj_cohort_GEP$A,digits=3,format="f"),", p-value: ",formatC(rocobj_cohort_GEP$p.value,digits=3,format="f"))
rocobj_cohort_IMPRES <- roc.area(as.numeric(as.factor(clin$Response))-1,clin$IMPRES)
rocobj_cohort_IMPRES_p <- formatC(rocobj_cohort_IMPRES$p.value,digits=3,format="f")
rocobj_cohort_IMPRES_auc <- formatC(rocobj_cohort_IMPRES$A,digits=3,format="f")
rocobj_cohort_IMPRES_string <- paste0("AUC: ",formatC(rocobj_cohort_IMPRES$A,digits=3,format="f"),", p-value: ",formatC(rocobj_cohort_IMPRES$p.value,digits=3,format="f"))
rocobj_cohort_TIDE <- roc.area(as.numeric(as.factor(clin$Response))-1,-1*clin$TIDE_Score)
rocobj_cohort_TIDE_p <- formatC(rocobj_cohort_TIDE$p.value,digits=3,format="e")
rocobj_cohort_TIDE_auc <- formatC(rocobj_cohort_TIDE$A,digits=3,format="f")
rocobj_cohort_TIDE_string <- paste0("AUC: ",formatC(rocobj_cohort_TIDE$A,digits=3,format="f"),", p-value: ",formatC(rocobj_cohort_TIDE$p.value,digits=2,format="e"))
rocobj_cohort_MHCII <- roc.area(as.numeric(as.factor(clin$Response))-1,clin$MHC_II)
rocobj_cohort_MHCII_p <- formatC(rocobj_cohort_MHCII$p.value,digits=3,format="f")
rocobj_cohort_MHCII_auc <- formatC(rocobj_cohort_MHCII$A,digits=3,format="f")
rocobj_cohort_MHCII_string <- paste0("AUC: ",formatC(rocobj_cohort_MHCII$A,digits=3,format="f"),", p-value: ",formatC(rocobj_cohort_MHCII$p.value,digits=3,format="f"))
rocobj_cohort_TCB <- roc.area(as.numeric(as.factor(clin$Response))-1,log10(clin$rna_tcb))
rocobj_cohort_TCB_p <- formatC(rocobj_cohort_TCB$p.value,digits=3,format="e")
rocobj_cohort_TCB_auc <- formatC(rocobj_cohort_TCB$A,digits=3,format="f")
rocobj_cohort_TCB_string <- paste0("AUC: ",formatC(rocobj_cohort_TCB$A,digits=3,format="f"),", p-value: ",formatC(rocobj_cohort_TCB$p.value,digits=2,format="e"))
rocobj_cohort_BCB <- roc.area(as.numeric(as.factor(clin$Response))-1,log10(clin$rna_bcb))
rocobj_cohort_BCB_p <- formatC(rocobj_cohort_BCB$p.value,digits=3,format="e")
rocobj_cohort_BCB_auc <- formatC(rocobj_cohort_BCB$A,digits=3,format="f")
rocobj_cohort_BCB_string <- paste0("AUC: ",formatC(rocobj_cohort_BCB$A,digits=3,format="f"),", p-value: ",formatC(rocobj_cohort_BCB$p.value,digits=2,format="e"))

rocobj_cohort_glm_MAP4K1TBX3 <- glm(as.numeric(as.factor(clin$Response))-1 ~ clin$MAP4K1 + clin$TBX3,family=binomial(link="logit"))
rocobj_cohort_MAP4K1TBX3 <- roc.area(as.numeric(as.factor(clin$Response))-1,rocobj_cohort_glm_MAP4K1TBX3$fitted.values)
rocobj_cohort_MAP4K1TBX3_p <- formatC(rocobj_cohort_MAP4K1TBX3$p.value,digits=2,format="e")
rocobj_cohort_MAP4K1TBX3_auc <- formatC(rocobj_cohort_MAP4K1TBX3$A,digits=3,format="f")
rocobj_cohort_MAP4K1TBX3_string <- paste0("AUC: ",formatC(rocobj_cohort_MAP4K1TBX3$A,digits=3,format="f"),", p-value: ",formatC(rocobj_cohort_MAP4K1TBX3$p.value,digits=2,format="e"))
rocobj_cohort_glm_MAP4K1AGER <- glm(as.numeric(as.factor(clin$Response))-1 ~ clin$MAP4K1 + clin$AGER,family=binomial(link="logit"))
rocobj_cohort_MAP4K1AGER <- roc.area(as.numeric(as.factor(clin$Response))-1,rocobj_cohort_glm_MAP4K1AGER$fitted.values)
rocobj_cohort_MAP4K1AGER_p <- formatC(rocobj_cohort_MAP4K1AGER$p.value,digits=2,format="e")
rocobj_cohort_MAP4K1AGER_auc <- formatC(rocobj_cohort_MAP4K1AGER$A,digits=3,format="f")
rocobj_cohort_MAP4K1AGER_string <- paste0("AUC: ",formatC(rocobj_cohort_MAP4K1AGER$A,digits=3,format="f"),", p-value: ",formatC(rocobj_cohort_MAP4K1AGER$p.value,digits=2,format="e"))

rocobj_cohort_glm_highosgeneszlowosgenesz <- glm(as.numeric(as.factor(clin$Response))-1 ~ clin$mean_highos_genes_z + clin$mean_lowos_genes_z,family=binomial(link="logit"))
rocobj_cohort_highosgeneszlowosgenesz <- roc.area(as.numeric(as.factor(clin$Response))-1,rocobj_cohort_glm_highosgeneszlowosgenesz$fitted.values)
rocobj_cohort_highosgeneszlowosgenesz_p <- formatC(rocobj_cohort_highosgeneszlowosgenesz$p.value,digits=2,format="e")
rocobj_cohort_highosgeneszlowosgenesz_auc <- formatC(rocobj_cohort_highosgeneszlowosgenesz$A,digits=3,format="f")
rocobj_cohort_highosgeneszlowosgenesz_string <- paste0("AUC: ",formatC(rocobj_cohort_highosgeneszlowosgenesz$A,digits=3,format="f"),", p-value: ",formatC(rocobj_cohort_highosgeneszlowosgenesz$p.value,digits=2,format="e"))

#col_set3_10 <- brewer.pal(n = 10, name = "Set3")
col_set3_11 <- brewer.pal(n = 11, name = "Set3")

pdf(pdf_string_1)
pROC::plot.roc(clin$Response,clin$CD274,col=col_set3_11[1],title="cohort response classification")
pROC::lines.roc(clin$Response,clin$CYT,col=col_set3_11[2])
pROC::lines.roc(clin$Response,clin$GEP,col=col_set3_11[3])
pROC::lines.roc(clin$Response,clin$IMPRES,col=col_set3_11[4])
pROC::lines.roc(clin$Response,clin$TIDE_Score,col=col_set3_11[5])
pROC::lines.roc(clin$Response,clin$MHC_II,col=col_set3_11[6])
pROC::lines.roc(clin$Response,log10(clin$rna_tcb),col=col_set3_11[7])
pROC::lines.roc(clin$Response,log10(clin$rna_bcb),col=col_set3_11[8])
pROC::lines.roc(clin$Response,rocobj_cohort_glm_MAP4K1TBX3$fitted.values,col=col_set3_11[9])
pROC::lines.roc(clin$Response,rocobj_cohort_glm_MAP4K1AGER$fitted.values,col=col_set3_11[10])
pROC::lines.roc(clin$Response,rocobj_cohort_glm_highosgeneszlowosgenesz$fitted.values,col=col_set3_11[11])
legend("bottomright", legend=c(paste0("CD274 ",rocobj_cohort_CD274_string), paste0("CYT ",rocobj_cohort_CYT_string), paste0("GEP ",rocobj_cohort_GEP_string), paste0("IMPRES ",rocobj_cohort_IMPRES_string),paste0("TIDE ",rocobj_cohort_TIDE_string),paste0("MHC II ",rocobj_cohort_MHCII_string),paste0("TCB ",rocobj_cohort_TCB_string),paste0("BCB ",rocobj_cohort_BCB_string),paste0("MAP4K1+TBX3 ",rocobj_cohort_MAP4K1TBX3_string),paste0("MAP4K1+AGER ",rocobj_cohort_MAP4K1AGER_string),paste0("High OS genes+Low OS genes ",rocobj_cohort_highosgeneszlowosgenesz_string)), col=c(col_set3_11), lwd=2,cex=0.71)
dev.off()

df_model_comparison_cohort_rnr <- data.frame(X1=c("CD274","GEP","IMPRES","TIDE","CYT","TCB","BCB","MHC II","MAP4K1+TBX3","MAP4K1+AGER","High OS genes+Low OS genes"),X2=c(rocobj_cohort_CD274_p,rocobj_cohort_GEP_p,rocobj_cohort_IMPRES_p,rocobj_cohort_TIDE_p,rocobj_cohort_CYT_p,rocobj_cohort_TCB_p,rocobj_cohort_BCB_p,rocobj_cohort_MHCII_p,rocobj_cohort_MAP4K1TBX3_p,rocobj_cohort_MAP4K1AGER_p,rocobj_cohort_highosgeneszlowosgenesz_p),X3=c(rocobj_cohort_CD274_auc,rocobj_cohort_GEP_auc,rocobj_cohort_IMPRES_auc,rocobj_cohort_TIDE_auc,rocobj_cohort_CYT_auc,rocobj_cohort_TCB_auc,rocobj_cohort_BCB_auc,rocobj_cohort_MHCII_auc,rocobj_cohort_MAP4K1TBX3_auc,rocobj_cohort_MAP4K1AGER_auc,rocobj_cohort_highosgeneszlowosgenesz_auc),stringsAsFactors=F)
df_model_comparison_cohort_rnr$X1 <- factor(df_model_comparison_cohort_rnr$X1,levels=c("CD274","GEP","IMPRES","TIDE","CYT","TCB","BCB","MHC II","MAP4K1+TBX3","MAP4K1+AGER","High OS genes+Low OS genes"))
df_model_comparison_cohort_rnr$X2 <- as.numeric(df_model_comparison_cohort_rnr$X2)
df_model_comparison_cohort_rnr$X3 <- as.numeric(df_model_comparison_cohort_rnr$X3)

g_rnr_compare_cohort_withBF3models <- ggplot(df_model_comparison_cohort_rnr,aes(x=as.numeric(X3),y=-1*log10(as.numeric(X2)))) + geom_point(aes(color=grepl("\\+",as.character(X1))),size=2) + geom_text_repel(aes(label=X1),size=6) + theme_classic(base_size=16) + geom_hline(yintercept = -1*log10(0.05),linetype = "dashed") + xlab("R/NR classification AUC") + ylab("-log10(R/NR classification AUC p-value)") + scale_color_manual(values=c("black","red")) + theme(legend.position = "none")
ggsave(pdf_string_2, g_rnr_compare_cohort_withBF3models,height=7,width=7)
return(df_model_comparison_cohort_rnr)
}
```

``` {.r}
skcm_tab <- readRDS("figure4_inputs/small_skcm_tab.rds")

skcm_clin <- readRDS("figure4_inputs/secondary_clin.rds")
skcm_clin$surv <- Surv(skcm_clin$overall_survival,skcm_clin$DOD=="DOD")

skcm_clin[,val_df_res_sig$gene] <- t(skcm_tab[rownames(val_df_res_sig),])
skcm_clin[,val_df_res_os_sig$gene] <- t(skcm_tab[rownames(val_df_res_os_sig),])
skcm_clin[,paste0(val_df_res_os_sig$gene,"_z")] <- apply(skcm_tab[rownames(val_df_res_os_sig),],1,scale)
skcm_clin$mean_lowos_genes_z <- rowMeans(skcm_clin[,paste0(val_df_res_os_sig[val_df_res_os_sig$log2FoldChange>0,]$gene,"_z")])
skcm_clin$mean_highos_genes_z <- rowMeans(skcm_clin[,paste0(val_df_res_os_sig[val_df_res_os_sig$log2FoldChange<0,]$gene,"_z")])

clin$pfs_surv <- Surv(clin$PFS,event=clin$progression_status=="Yes")
skcm_clin$pfs_surv <- Surv(skcm_clin$PFS,event=skcm_clin$progression_status=="Yes")

df_model_comparison_cpb <- make_os_p_cindex_plot(clin,"figure4_outputs/S23e.cpb_rna_survival_coxcontinuous_pvaluecindex_withBF3models.pdf")
knitr::include_graphics("figure4_outputs/S23e.cpb_rna_survival_coxcontinuous_pvaluecindex_withBF3models.pdf")
```

<embed src="figure4_outputs/S23e.cpb_rna_survival_coxcontinuous_pvaluecindex_withBF3models.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
df_model_comparison_cpb_pfs <- make_pfs_p_cindex_plot(clin,"figure4_outputs/cpb_rna_survival_coxcontinuous_pvaluecindex_withBF3models_PFS.pdf")
```

    ##                            X1           X2        X3
    ## 1                       CD274 0.7159117012 0.5085514
    ## 2                         GEP 0.0403970291 0.5638767
    ## 3                      IMPRES 0.0351883609 0.5760560
    ## 4                        TIDE 0.0146353875 0.6000259
    ## 5                         CYT 0.0672878646 0.5584348
    ## 6                         TCB 0.0068051414 0.5835709
    ## 7                         BCB 0.0023859674 0.5688002
    ## 8                      MHC II 0.0007576778 0.5506608
    ## 9                 MAP4K1+TBX3 0.0011923014 0.6281420
    ## 10                MAP4K1+AGER 0.0019710587 0.6159627
    ## 11 High OS genes+Low OS genes 0.0006413309 0.6302151

``` {.r}
knitr::include_graphics("figure4_outputs/cpb_rna_survival_coxcontinuous_pvaluecindex_withBF3models_PFS.pdf")
```

<embed src="figure4_outputs/cpb_rna_survival_coxcontinuous_pvaluecindex_withBF3models_PFS.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
df_model_comparison_cpb_rnr <- make_rnr_p_auc_plot(clin,"figure4_outputs/S23g.cpb_rna_ROC_withpvalues_withBF3models.pdf","figure4_outputs/S23f.cpb_rna_RNR_pvalueAUC_withBF3models.pdf")
knitr::include_graphics("figure4_outputs/S23f.cpb_rna_RNR_pvalueAUC_withBF3models.pdf")
```

<embed src="figure4_outputs/S23f.cpb_rna_RNR_pvalueAUC_withBF3models.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
knitr::include_graphics("figure4_outputs/S23g.cpb_rna_ROC_withpvalues_withBF3models.pdf")
```

<embed src="figure4_outputs/S23g.cpb_rna_ROC_withpvalues_withBF3models.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
print(df_model_comparison_cpb_pfs)
```

    ##                            X1           X2        X3
    ## 1                       CD274 0.7159117012 0.5085514
    ## 2                         GEP 0.0403970291 0.5638767
    ## 3                      IMPRES 0.0351883609 0.5760560
    ## 4                        TIDE 0.0146353875 0.6000259
    ## 5                         CYT 0.0672878646 0.5584348
    ## 6                         TCB 0.0068051414 0.5835709
    ## 7                         BCB 0.0023859674 0.5688002
    ## 8                      MHC II 0.0007576778 0.5506608
    ## 9                 MAP4K1+TBX3 0.0011923014 0.6281420
    ## 10                MAP4K1+AGER 0.0019710587 0.6159627
    ## 11 High OS genes+Low OS genes 0.0006413309 0.6302151

``` {.r}
df_model_comparison_cpb_survmodelrnrmodel <- data.frame(X1=df_model_comparison_cpb$X1,X2=df_model_comparison_cpb$X2,X3=df_model_comparison_cpb$X3,
                                                        X4=df_model_comparison_cpb_rnr$X2,X5=df_model_comparison_cpb_rnr$X3,
                                                        stringsAsFactors=F)
df_model_comparison_cpb_survmodelrnrmodel$X1 <- factor(as.character(df_model_comparison_cpb_survmodelrnrmodel$X1),
                                             levels=c("CD274","GEP","IMPRES","TIDE","CYT","TCB","BCB","MHC II","MAP4K1+TBX3","MAP4K1+AGER","High OS genes+Low OS genes"))
df_model_comparison_cpb_survmodelrnrmodel$X2 <- as.numeric(df_model_comparison_cpb_survmodelrnrmodel$X2)
df_model_comparison_cpb_survmodelrnrmodel$X3 <- as.numeric(df_model_comparison_cpb_survmodelrnrmodel$X3)
df_model_comparison_cpb_survmodelrnrmodel$X4 <- as.numeric(df_model_comparison_cpb_survmodelrnrmodel$X4)
df_model_comparison_cpb_survmodelrnrmodel$X5 <- as.numeric(df_model_comparison_cpb_survmodelrnrmodel$X5)

g_survmodelrnrmodel_compare_cpb_withBF3models_pvals <- ggplot(df_model_comparison_cpb_survmodelrnrmodel,aes(x=-1*log10(as.numeric(X2)),y=-1*log10(as.numeric(X4)))) + 
geom_point(aes(color=grepl("\\+",as.character(X1))),size=2) + geom_text_repel(aes(label=X1),size=6) + theme_classic(base_size=16) + 
geom_hline(yintercept = -1*log10(0.05),linetype = "dashed") + geom_vline(xintercept = -1*log10(0.05),linetype = "dashed") +
xlab("-log10(OS Cox model Log-Rank p-value)") + ylab("-log10(R/NR classification AUC p-value)") + 
scale_color_manual(values=c("black","red")) + theme(legend.position = "none")
ggsave("figure4_outputs/4c.cpb_rna_survival_coxcontinuous_pvaluecindex_withBF3models_pvals.pdf", g_survmodelrnrmodel_compare_cpb_withBF3models_pvals,height=7,width=7)
knitr::include_graphics("figure4_outputs/4c.cpb_rna_survival_coxcontinuous_pvaluecindex_withBF3models_pvals.pdf")
```

<embed src="figure4_outputs/4c.cpb_rna_survival_coxcontinuous_pvaluecindex_withBF3models_pvals.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
g_survmodelrnrmodel_compare_cpb_withBF3models_effectsize <- ggplot(df_model_comparison_cpb_survmodelrnrmodel,aes(x=as.numeric(X3),y=as.numeric(X5))) + 
geom_point(aes(color=grepl("\\+",as.character(X1))),size=2) + geom_text_repel(aes(label=X1),size=6) + theme_classic(base_size=16) + xlab("OS Cox model C-index") + ylab("R/NR classification AUC") + 
scale_color_manual(values=c("black","red")) + theme(legend.position = "none")
ggsave("figure4_outputs/4d.cpb_rna_survival_coxcontinuous_pvaluecindex_withBF3models_effectsize.pdf", g_survmodelrnrmodel_compare_cpb_withBF3models_effectsize,height=7,width=7)
knitr::include_graphics("figure4_outputs/4d.cpb_rna_survival_coxcontinuous_pvaluecindex_withBF3models_effectsize.pdf")
```

<embed src="figure4_outputs/4d.cpb_rna_survival_coxcontinuous_pvaluecindex_withBF3models_effectsize.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
df_model_comparison_val <- make_os_p_cindex_plot(skcm_clin,"figure4_outputs/S26b.val_rna_survival_coxcontinuous_pvaluecindex_withBF3models.pdf")
knitr::include_graphics("figure4_outputs/S26b.val_rna_survival_coxcontinuous_pvaluecindex_withBF3models.pdf")
```

<embed src="figure4_outputs/S26b.val_rna_survival_coxcontinuous_pvaluecindex_withBF3models.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
df_model_comparison_val_pfs <- make_pfs_p_cindex_plot(skcm_clin,"figure4_outputs/val_rna_survival_coxcontinuous_pvaluecindex_PFS_withBF3models.pdf")
```

    ##                            X1           X2        X3
    ## 1                       CD274 0.0020885203 0.5956763
    ## 2                         GEP 0.0004041490 0.6043742
    ## 3                      IMPRES 0.0035949049 0.5883860
    ## 4                        TIDE 0.2213872829 0.5410351
    ## 5                         CYT 0.0001241811 0.6069727
    ## 6                         TCB 0.0010298063 0.5937635
    ## 7                         BCB 0.0015503101 0.5931861
    ## 8                      MHC II 0.0005577493 0.6058900
    ## 9                 MAP4K1+TBX3 0.0051765440 0.5880612
    ## 10                MAP4K1+AGER 0.0073283189 0.5975170
    ## 11 High OS genes+Low OS genes 0.0014684244 0.6048795

``` {.r}
knitr::include_graphics("figure4_outputs/val_rna_survival_coxcontinuous_pvaluecindex_PFS_withBF3models.pdf")
```

<embed src="figure4_outputs/val_rna_survival_coxcontinuous_pvaluecindex_PFS_withBF3models.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
df_model_comparison_val_rnr <- make_rnr_p_auc_plot(skcm_clin,"figure4_outputs/S26d.val_rna_ROC_withpvalues_withBF3models.pdf","figure4_outputs/S26c.val_rna_RNR_pvalueAUC_withBF3models.pdf")
knitr::include_graphics("figure4_outputs/S26c.val_rna_RNR_pvalueAUC_withBF3models.pdf")
```

<embed src="figure4_outputs/S26c.val_rna_RNR_pvalueAUC_withBF3models.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
knitr::include_graphics("figure4_outputs/S26d.val_rna_ROC_withpvalues_withBF3models.pdf")
```

<embed src="figure4_outputs/S26d.val_rna_ROC_withpvalues_withBF3models.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
print(df_model_comparison_val_pfs)
```

    ##                            X1           X2        X3
    ## 1                       CD274 0.0020885203 0.5956763
    ## 2                         GEP 0.0004041490 0.6043742
    ## 3                      IMPRES 0.0035949049 0.5883860
    ## 4                        TIDE 0.2213872829 0.5410351
    ## 5                         CYT 0.0001241811 0.6069727
    ## 6                         TCB 0.0010298063 0.5937635
    ## 7                         BCB 0.0015503101 0.5931861
    ## 8                      MHC II 0.0005577493 0.6058900
    ## 9                 MAP4K1+TBX3 0.0051765440 0.5880612
    ## 10                MAP4K1+AGER 0.0073283189 0.5975170
    ## 11 High OS genes+Low OS genes 0.0014684244 0.6048795

``` {.r}
df_model_comparison_val_survmodelrnrmodel <- data.frame(X1=df_model_comparison_val$X1,X2=df_model_comparison_val$X2,X3=df_model_comparison_val$X3,
                                                        X4=df_model_comparison_val_rnr$X2,X5=df_model_comparison_val_rnr$X3,
                                                        stringsAsFactors=F)
df_model_comparison_val_survmodelrnrmodel$X1 <- factor(as.character(df_model_comparison_val_survmodelrnrmodel$X1),
                                             levels=c("CD274","GEP","IMPRES","TIDE","CYT","TCB","BCB","MHC II","MAP4K1+TBX3","MAP4K1+AGER","High OS genes+Low OS genes"))
df_model_comparison_val_survmodelrnrmodel$X2 <- as.numeric(df_model_comparison_val_survmodelrnrmodel$X2)
df_model_comparison_val_survmodelrnrmodel$X3 <- as.numeric(df_model_comparison_val_survmodelrnrmodel$X3)
df_model_comparison_val_survmodelrnrmodel$X4 <- as.numeric(df_model_comparison_val_survmodelrnrmodel$X4)
df_model_comparison_val_survmodelrnrmodel$X5 <- as.numeric(df_model_comparison_val_survmodelrnrmodel$X5)

g_survmodelrnrmodel_compare_val_withBF3models_pvals <- ggplot(df_model_comparison_val_survmodelrnrmodel,aes(x=-1*log10(as.numeric(X2)),y=-1*log10(as.numeric(X4)))) + 
geom_point(aes(color=grepl("\\+",as.character(X1))),size=2) + geom_text_repel(aes(label=X1),size=6) + theme_classic(base_size=16) + 
geom_hline(yintercept = -1*log10(0.05),linetype = "dashed") + geom_vline(xintercept = -1*log10(0.05),linetype = "dashed") +
xlab("-log10(OS Cox model Log-Rank p-value)") + ylab("-log10(R/NR classification AUC p-value)") + 
scale_color_manual(values=c("black","red")) + theme(legend.position = "none")
ggsave("figure4_outputs/4g.val_rna_survival_coxcontinuous_pvaluecindex_withBF3models_pvals.pdf", g_survmodelrnrmodel_compare_val_withBF3models_pvals,height=7,width=7)
knitr::include_graphics("figure4_outputs/4g.val_rna_survival_coxcontinuous_pvaluecindex_withBF3models_pvals.pdf")
```

<embed src="figure4_outputs/4g.val_rna_survival_coxcontinuous_pvaluecindex_withBF3models_pvals.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
g_survmodelrnrmodel_compare_val_withBF3models_effectsize <- ggplot(df_model_comparison_val_survmodelrnrmodel,aes(x=as.numeric(X3),y=as.numeric(X5))) + 
geom_point(aes(color=grepl("\\+",as.character(X1))),size=2) + geom_text_repel(aes(label=X1),size=6) + theme_classic(base_size=16) + xlab("OS Cox model C-index") + ylab("R/NR classification AUC") + 
scale_color_manual(values=c("black","red")) + theme(legend.position = "none")
ggsave("figure4_outputs/4h.val_rna_survival_coxcontinuous_pvaluecindex_withBF3models_effectsize.pdf", g_survmodelrnrmodel_compare_val_withBF3models_effectsize,height=7,width=7)
knitr::include_graphics("figure4_outputs/4h.val_rna_survival_coxcontinuous_pvaluecindex_withBF3models_effectsize.pdf")
```

<embed src="figure4_outputs/4h.val_rna_survival_coxcontinuous_pvaluecindex_withBF3models_effectsize.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
print(df_model_comparison_cpb_pfs)
```

    ##                            X1           X2        X3
    ## 1                       CD274 0.7159117012 0.5085514
    ## 2                         GEP 0.0403970291 0.5638767
    ## 3                      IMPRES 0.0351883609 0.5760560
    ## 4                        TIDE 0.0146353875 0.6000259
    ## 5                         CYT 0.0672878646 0.5584348
    ## 6                         TCB 0.0068051414 0.5835709
    ## 7                         BCB 0.0023859674 0.5688002
    ## 8                      MHC II 0.0007576778 0.5506608
    ## 9                 MAP4K1+TBX3 0.0011923014 0.6281420
    ## 10                MAP4K1+AGER 0.0019710587 0.6159627
    ## 11 High OS genes+Low OS genes 0.0006413309 0.6302151

``` {.r}
out_cpb_model_comparison <- data.frame(model=df_model_comparison_cpb$X1,stringsAsFactors=F)
out_cpb_model_comparison$os_cox_logrank_p <- df_model_comparison_cpb$X2
out_cpb_model_comparison$os_cindex <- df_model_comparison_cpb$X3
out_cpb_model_comparison$pfs_cox_logrank_p <- df_model_comparison_cpb_pfs$X2
out_cpb_model_comparison$pfs_cindex <- df_model_comparison_cpb_pfs$X3
out_cpb_model_comparison$rnr_roc_p <- df_model_comparison_cpb_rnr$X2
out_cpb_model_comparison$rnr_roc <- df_model_comparison_cpb_rnr$X3
write.table(out_cpb_model_comparison,"figure4_outputs/CPB_model_comparison.txt",quote=F,sep="\t")

out_val_model_comparison <- data.frame(model=df_model_comparison_val$X1,stringsAsFactors=F)
out_val_model_comparison$os_cox_logrank_p <- df_model_comparison_val$X2
out_val_model_comparison$os_cindex <- df_model_comparison_val$X3
out_val_model_comparison$pfs_cox_logrank_p <- df_model_comparison_val_pfs$X2
out_val_model_comparison$pfs_cindex <- df_model_comparison_val_pfs$X3
out_val_model_comparison$rnr_roc_p <- df_model_comparison_val_rnr$X2
out_val_model_comparison$rnr_roc <- df_model_comparison_val_rnr$X3
write.table(out_val_model_comparison,"figure4_outputs/Val_model_comparison.txt",quote=F,sep="\t")
```

Figure S24

``` {.r}
load("figure4_inputs/crossvalidation.RData")

all_os_results_meanz <- all_os_results[all_os_results$gene_pair=="mean.highos.genes.z_mean.lowos.genes.z",]
all_os_results_meanz_top25 <- all_os_results[all_os_results$gene_pair=="mean.highos.genes.z.top25_mean.lowos.genes.z.top25",]

all_rnr_results_meanz <- all_rnr_results[all_rnr_results$gene_pair=="mean.r.genes.z_mean.nr.genes.z",]
all_rnr_results_meanz_top25 <- all_rnr_results[all_rnr_results$gene_pair=="mean.r.genes.z.top25_mean.nr.genes.z.top25",]

all_os_bothsig_results_meanz <- all_os_bothsig_results[all_os_bothsig_results$gene_pair=="mean.highos.genes.z_mean.lowos.genes.z",]
all_os_bothsig_results_meanz_top25 <- all_os_bothsig_results[all_os_bothsig_results$gene_pair=="mean.highos.genes.z.top25_mean.lowos.genes.z.top25",]

tmp_all_split_index <- expand.grid(c("80_20","75_25","70_30","66_33","50_50"),c(1:250))
all_split_index <- paste0(tmp_all_split_index$Var1,"_",tmp_all_split_index$Var2)
all_os_results$f_split_index <- factor(all_os_results$split_index,levels=all_split_index)
tmp_table_all_os_results <- table(all_os_results$f_split_index)
osgenes_ossig <- data.frame(cbind(split_index=names(tmp_table_all_os_results),count=tmp_table_all_os_results),stringsAsFactors=F)
osgenes_ossig$count <- as.numeric(osgenes_ossig$count)
osgenes_ossig$split <- paste0(unlist(lapply(strsplit(osgenes_ossig$split_index,"_"),"[",1)),"_",unlist(lapply(strsplit(osgenes_ossig$split_index,"_"),"[",2)))

g_osgenes_ossig <- ggplot(osgenes_ossig,aes(x=count+1,color=split)) + stat_ecdf() + theme_classic() + ylab("Count of High/Low OS Bonferroni p<0.05 genes") + xlab("count + 1") + scale_x_log10()
ggsave("figure4_outputs/S24a.OSgenes_OSsig_count.pdf",g_osgenes_ossig,height=7,width=7)
knitr::include_graphics("figure4_outputs/S24a.OSgenes_OSsig_count.pdf")
```

<embed src="figure4_outputs/S24a.OSgenes_OSsig_count.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
all_rnr_results$f_split_index <- factor(all_rnr_results$split_index,levels=all_split_index)
tmp_table_all_rnr_results <- table(all_rnr_results$f_split_index)
rnrgenes_rnrsig <- data.frame(cbind(split_index=names(tmp_table_all_rnr_results),count=tmp_table_all_rnr_results),stringsAsFactors=F)
rnrgenes_rnrsig$count <- as.numeric(rnrgenes_rnrsig$count)
rnrgenes_rnrsig$split <- paste0(unlist(lapply(strsplit(rnrgenes_rnrsig$split_index,"_"),"[",1)),"_",unlist(lapply(strsplit(rnrgenes_rnrsig$split_index,"_"),"[",2)))

all_os_bothsig_results$f_split_index <- factor(all_os_bothsig_results$split_index,levels=all_split_index)
tmp_table_all_os_bothsig_results <- table(all_os_bothsig_results$f_split_index)
osgenes_bothsig <- data.frame(cbind(split_index=names(tmp_table_all_os_bothsig_results),count=tmp_table_all_os_bothsig_results),stringsAsFactors=F)
osgenes_bothsig$count <- as.numeric(osgenes_bothsig$count)
osgenes_bothsig$split <- paste0(unlist(lapply(strsplit(osgenes_bothsig$split_index,"_"),"[",1)),"_",unlist(lapply(strsplit(osgenes_bothsig$split_index,"_"),"[",2)))

g_osgenes_bothsig <- ggplot(osgenes_bothsig,aes(x=count+1,color=split)) + stat_ecdf() + theme_classic() + ylab("Count of High/Low OS Bonferroni p<0.05 genes for OS and R/NR") +
xlab("count + 1") + scale_x_log10()
ggsave("figure4_outputs/S24b.OSgenes_Bothsig_count.pdf",g_osgenes_bothsig,height=7,width=7)
knitr::include_graphics("figure4_outputs/S24b.OSgenes_Bothsig_count.pdf")
```

<embed src="figure4_outputs/S24b.OSgenes_Bothsig_count.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
g_mean_z_hr_1_bonferronipval <- ggplot(all_os_results_meanz,aes(x=gene1gene2_continuous_os_hr_1,y=test_gene1gene2_continuous_os_hr_1)) +
geom_point(aes(color=test_gene1gene2_continuous_logrank_bonferroni_p<0.05)) +
theme_classic() + xlab("training set mean highOS genes z score HR") + ylab("validation set mean highOS genes z score HR") + scale_x_log10() + scale_y_log10() + facet_grid(~ split) +
geom_smooth(method="lm",se=FALSE,color="red") + annotation_logticks() + coord_fixed() +
scale_color_manual(values=c("black","red"),name="validation set bonferroni p < 0.05") + theme(legend.position="bottom")
ggsave("figure4_outputs/S24c.meanz_train_test_highOS_HR_bonferronipval.pdf",g_mean_z_hr_1_bonferronipval,height=7,width=14)
knitr::include_graphics("figure4_outputs/S24c.meanz_train_test_highOS_HR_bonferronipval.pdf")
```

<embed src="figure4_outputs/S24c.meanz_train_test_highOS_HR_bonferronipval.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
g_mean_z_hr_2_bonferronipval <- ggplot(all_os_results_meanz,aes(x=gene1gene2_continuous_os_hr_2,y=test_gene1gene2_continuous_os_hr_2)) +
geom_point(aes(color=test_gene1gene2_continuous_logrank_bonferroni_p<0.05)) +
theme_classic() + xlab("training set mean lowOS genes z score HR") + ylab("validation set mean lowOS genes z score HR") + scale_x_log10() + scale_y_log10() + facet_grid(~ split) +
geom_smooth(method="lm",se=FALSE,color="red") +  annotation_logticks() + coord_fixed() +
scale_color_manual(values=c("black","red"),name="validation set bonferroni p < 0.05") + theme(legend.position="bottom")
ggsave("figure4_outputs/S24d.meanz_train_test_lowOS_HR_bonferronipval.pdf",g_mean_z_hr_2_bonferronipval,height=7,width=14)
knitr::include_graphics("figure4_outputs/S24d.meanz_train_test_lowOS_HR_bonferronipval.pdf")
```

<embed src="figure4_outputs/S24d.meanz_train_test_lowOS_HR_bonferronipval.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
g_mean_z_hr_1_top25_bonferronipval <- ggplot(all_os_results_meanz_top25,aes(x=gene1gene2_continuous_os_hr_1,y=test_gene1gene2_continuous_os_hr_1)) +
geom_point(aes(color=test_gene1gene2_continuous_logrank_bonferroni_p<0.05)) +
theme_classic() + xlab("training set mean highOS top 25 genes z score HR") + ylab("validation set mean highOS top 25 genes z score HR") + scale_x_log10() + scale_y_log10() + facet_grid(~ split) +
geom_smooth(method="lm",se=FALSE,color="red") + annotation_logticks() + coord_fixed() +
scale_color_manual(values=c("black","red"),name="validation set bonferroni p < 0.05") + theme(legend.position="bottom")
ggsave("figure4_outputs/S24e.meanz_train_test_highOS_HR_top25_bonferronipval.pdf",g_mean_z_hr_1_top25_bonferronipval,height=7,width=14)
knitr::include_graphics("figure4_outputs/S24e.meanz_train_test_highOS_HR_top25_bonferronipval.pdf")
```

<embed src="figure4_outputs/S24e.meanz_train_test_highOS_HR_top25_bonferronipval.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
g_mean_z_hr_2_top25_bonferronipval <- ggplot(all_os_results_meanz_top25,aes(x=gene1gene2_continuous_os_hr_2,y=test_gene1gene2_continuous_os_hr_2)) +
geom_point(aes(color=test_gene1gene2_continuous_logrank_bonferroni_p<0.05)) +
theme_classic() + xlab("training set mean lowOS top 25 genes z score HR") + ylab("validation set mean lowOS top 25 genes z score HR") + scale_x_log10() + scale_y_log10() + facet_grid(~ split) +
geom_smooth(method="lm",se=FALSE,color="red") +  annotation_logticks() + coord_fixed() +
scale_color_manual(values=c("black","red"),name="validation set bonferroni p < 0.05") + theme(legend.position="bottom")
ggsave("figure4_outputs/S24f.meanz_train_test_lowOS_HR_top25_bonferronipval.pdf",g_mean_z_hr_2_top25_bonferronipval,height=7,width=14)
knitr::include_graphics("figure4_outputs/S24f.meanz_train_test_lowOS_HR_top25_bonferronipval.pdf")
```

<embed src="figure4_outputs/S24f.meanz_train_test_lowOS_HR_top25_bonferronipval.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
print(data.frame(cbind(names(sort(table(all_os_results[all_os_results$split=="80_20",]$gene_pair),decreasing=T)[1:10]),
as.numeric(sort(table(all_os_results[all_os_results$split=="80_20",]$gene_pair),decreasing=T)[1:10]),
as.numeric(sort(table(all_os_results[all_os_results$split=="80_20",]$gene_pair),decreasing=T)[1:10]/250),"80_20"),stringsAsFactors=F),row.names=F,justify="left")
```

    ##                                                  X1  X2    X3    X4
    ##              mean.highos.genes.z_mean.lowos.genes.z 249 0.996 80_20
    ##  mean.highos.genes.z.top25_mean.lowos.genes.z.top25 244 0.976 80_20
    ##                                     RAB11FIP4_USP49  81 0.324 80_20
    ##                                        MAP4K1_USP49  78 0.312 80_20
    ##                                    RAB11FIP4_RAD54B  68 0.272 80_20
    ##                                      RAB11FIP4_TBX3  54 0.216 80_20
    ##                                           TBX3_IDO1  51 0.204 80_20
    ##                                        RAD54B_ABHD6  49 0.196 80_20
    ##                                         IKZF3_USP49  46 0.184 80_20
    ##                                         MAP4K1_AGER  43 0.172 80_20

``` {.r}
print(data.frame(cbind(names(sort(table(all_os_results[all_os_results$split=="75_25",]$gene_pair),decreasing=T)[1:10]),
as.numeric(sort(table(all_os_results[all_os_results$split=="75_25",]$gene_pair),decreasing=T)[1:10]),
as.numeric(sort(table(all_os_results[all_os_results$split=="75_25",]$gene_pair),decreasing=T)[1:10]/250),"75_25"),stringsAsFactors=F),row.names=F,justify="left")
```

    ##                                                  X1  X2    X3    X4
    ##  mean.highos.genes.z.top25_mean.lowos.genes.z.top25 243 0.972 75_25
    ##              mean.highos.genes.z_mean.lowos.genes.z 241 0.964 75_25
    ##                                        MAP4K1_USP49  64 0.256 75_25
    ##                                     RAB11FIP4_USP49  63 0.252 75_25
    ##                                    RAB11FIP4_RAD54B  36 0.144 75_25
    ##                                         IKZF3_USP49  34 0.136 75_25
    ##                                         MAP4K1_AGER  34 0.136 75_25
    ##                                      RAB11FIP4_TBX3  34 0.136 75_25
    ##                                         USP49_ABHD6  34 0.136 75_25
    ##                                        RAD54B_ABHD6  33 0.132 75_25

``` {.r}
print(data.frame(cbind(names(sort(table(all_os_results[all_os_results$split=="70_30",]$gene_pair),decreasing=T)[1:10]),
as.numeric(sort(table(all_os_results[all_os_results$split=="70_30",]$gene_pair),decreasing=T)[1:10]),
as.numeric(sort(table(all_os_results[all_os_results$split=="70_30",]$gene_pair),decreasing=T)[1:10]/250),"70_30"),stringsAsFactors=F),row.names=F,justify="left")
```

    ##                                                  X1  X2    X3    X4
    ##  mean.highos.genes.z.top25_mean.lowos.genes.z.top25 242 0.968 70_30
    ##              mean.highos.genes.z_mean.lowos.genes.z 237 0.948 70_30
    ##                                     RAB11FIP4_USP49  42 0.168 70_30
    ##                                        MAP4K1_USP49  40  0.16 70_30
    ##                                    RAB11FIP4_RAD54B  27 0.108 70_30
    ##                                         USP49_ABHD6  26 0.104 70_30
    ##                                         USP49_FAIM3  26 0.104 70_30
    ##                                        RAD54B_ABHD6  24 0.096 70_30
    ##                                           TBX3_IDO1  24 0.096 70_30
    ##                                          SLAIN1_PLN  23 0.092 70_30

``` {.r}
print(data.frame(cbind(names(sort(table(all_os_results[all_os_results$split=="66_33",]$gene_pair),decreasing=T)[1:10]),
as.numeric(sort(table(all_os_results[all_os_results$split=="66_33",]$gene_pair),decreasing=T)[1:10]),
as.numeric(sort(table(all_os_results[all_os_results$split=="66_33",]$gene_pair),decreasing=T)[1:10]/250),"66_33"),stringsAsFactors=F),row.names=F,justify="left")
```

    ##                                                  X1  X2    X3    X4
    ##  mean.highos.genes.z.top25_mean.lowos.genes.z.top25 241 0.964 66_33
    ##              mean.highos.genes.z_mean.lowos.genes.z 228 0.912 66_33
    ##                                     RAB11FIP4_USP49  32 0.128 66_33
    ##                                        MAP4K1_USP49  25   0.1 66_33
    ##                                         USP49_FAIM3  23 0.092 66_33
    ##                                         USP49_FYCO1  21 0.084 66_33
    ##                                    RAB11FIP4_RAD54B  20  0.08 66_33
    ##                                          SLAIN1_PLN  20  0.08 66_33
    ##                                        RAD54B_ABHD6  19 0.076 66_33
    ##                                          IDO1_USP49  18 0.072 66_33

``` {.r}
print(data.frame(cbind(names(sort(table(all_os_results[all_os_results$split=="50_50",]$gene_pair),decreasing=T)[1:10]),
as.numeric(sort(table(all_os_results[all_os_results$split=="50_50",]$gene_pair),decreasing=T)[1:10]),
as.numeric(sort(table(all_os_results[all_os_results$split=="50_50",]$gene_pair),decreasing=T)[1:10]/250),"50_50"),stringsAsFactors=F),row.names=F,justify="left")
```

    ##                                                  X1  X2    X3    X4
    ##  mean.highos.genes.z.top25_mean.lowos.genes.z.top25 243 0.972 50_50
    ##              mean.highos.genes.z_mean.lowos.genes.z 173 0.692 50_50
    ##                                     RAB11FIP4_USP49  12 0.048 50_50
    ##                                               GEP__  11 0.044 50_50
    ##                                               CYT__   9 0.036 50_50
    ##                                        USP49_CXCL13   8 0.032 50_50
    ##                                         USP49_FAIM3   8 0.032 50_50
    ##                                       CCNDBP1_USP49   7 0.028 50_50
    ##                                           PLN_FAIM3   7 0.028 50_50
    ##                                         USP49_FYCO1   7 0.028 50_50

``` {.r}
print("-------------------")
```

    ## [1] "-------------------"

``` {.r}
print(data.frame(cbind(names(sort(table(all_os_results[all_os_results$test_gene1gene2_continuous_logrank_bonferroni_p<0.05&all_os_results$split=="80_20",]$gene_pair),decreasing=T)[1:10]),
as.numeric(sort(table(all_os_results[all_os_results$test_gene1gene2_continuous_logrank_bonferroni_p<0.05&all_os_results$split=="80_20",]$gene_pair),decreasing=T)[1:10]),
as.numeric(sort(table(all_os_results[all_os_results$test_gene1gene2_continuous_logrank_bonferroni_p<0.05&all_os_results$split=="80_20",]$gene_pair),decreasing=T)[1:10]/250),"80_20"),
stringsAsFactors=F),row.names=F,justify="left")
```

    ##                                                  X1   X2    X3    X4
    ##  mean.highos.genes.z.top25_mean.lowos.genes.z.top25   14 0.056 80_20
    ##              mean.highos.genes.z_mean.lowos.genes.z    9 0.036 80_20
    ##                                         MAP4K1_AGER    2 0.008 80_20
    ##                                        MAP4K1_USP49    2 0.008 80_20
    ##                                          AGER_FAIM3    1 0.004 80_20
    ##                                      CCNDBP1_RAD54B    1 0.004 80_20
    ##                                                <NA> <NA>  <NA> 80_20
    ##                                                <NA> <NA>  <NA> 80_20
    ##                                                <NA> <NA>  <NA> 80_20
    ##                                                <NA> <NA>  <NA> 80_20

``` {.r}
print(data.frame(cbind(names(sort(table(all_os_results[all_os_results$test_gene1gene2_continuous_logrank_bonferroni_p<0.05&all_os_results$split=="75_25",]$gene_pair),decreasing=T)[1:10]),
as.numeric(sort(table(all_os_results[all_os_results$test_gene1gene2_continuous_logrank_bonferroni_p<0.05&all_os_results$split=="75_25",]$gene_pair),decreasing=T)[1:10]),
as.numeric(sort(table(all_os_results[all_os_results$test_gene1gene2_continuous_logrank_bonferroni_p<0.05&all_os_results$split=="75_25",]$gene_pair),decreasing=T)[1:10]/250),"75_25"),
stringsAsFactors=F),row.names=F,justify="left")
```

    ##                                                  X1 X2    X3    X4
    ##  mean.highos.genes.z.top25_mean.lowos.genes.z.top25 12 0.048 75_25
    ##              mean.highos.genes.z_mean.lowos.genes.z  5  0.02 75_25
    ##                                     RAB11FIP4_USP49  4 0.016 75_25
    ##                                         MAP4K1_AGER  2 0.008 75_25
    ##                                        MAP4K1_USP49  2 0.008 75_25
    ##                                          ABHD6_RCC2  1 0.004 75_25
    ##                                    RAB11FIP4_ADAM12  1 0.004 75_25
    ##                                        SH3GL1_ABHD6  1 0.004 75_25
    ##                                         SH3GL1_IL16  1 0.004 75_25
    ##                                       SH3GL1_OXNAD1  1 0.004 75_25

``` {.r}
print(data.frame(cbind(names(sort(table(all_os_results[all_os_results$test_gene1gene2_continuous_logrank_bonferroni_p<0.05&all_os_results$split=="70_30",]$gene_pair),decreasing=T)[1:10]),
as.numeric(sort(table(all_os_results[all_os_results$test_gene1gene2_continuous_logrank_bonferroni_p<0.05&all_os_results$split=="70_30",]$gene_pair),decreasing=T)[1:10]),
as.numeric(sort(table(all_os_results[all_os_results$test_gene1gene2_continuous_logrank_bonferroni_p<0.05&all_os_results$split=="70_30",]$gene_pair),decreasing=T)[1:10]/250),"70_30"),
stringsAsFactors=F),row.names=F,justify="left")
```

    ##                                                  X1 X2    X3    X4
    ##  mean.highos.genes.z.top25_mean.lowos.genes.z.top25 15  0.06 70_30
    ##              mean.highos.genes.z_mean.lowos.genes.z  6 0.024 70_30
    ##                                         MAP4K1_AGER  2 0.008 70_30
    ##                                        ABHD6_ATAD3A  1 0.004 70_30
    ##                                        ATP2A3_USP49  1 0.004 70_30
    ##                                     CD79A_RAB11FIP4  1 0.004 70_30
    ##                                         EFNB2_ABHD6  1 0.004 70_30
    ##                                        MAP4K1_BCAS3  1 0.004 70_30
    ##                                        MAP4K1_USP49  1 0.004 70_30
    ##                                          NREP_ABHD6  1 0.004 70_30

``` {.r}
print(data.frame(cbind(names(sort(table(all_os_results[all_os_results$test_gene1gene2_continuous_logrank_bonferroni_p<0.05&all_os_results$split=="66_33",]$gene_pair),decreasing=T)[1:10]),
as.numeric(sort(table(all_os_results[all_os_results$test_gene1gene2_continuous_logrank_bonferroni_p<0.05&all_os_results$split=="66_33",]$gene_pair),decreasing=T)[1:10]),
as.numeric(sort(table(all_os_results[all_os_results$test_gene1gene2_continuous_logrank_bonferroni_p<0.05&all_os_results$split=="66_33",]$gene_pair),decreasing=T)[1:10]/250),"66_33"),
stringsAsFactors=F),row.names=F,justify="left")
```

    ##                                                  X1 X2    X3    X4
    ##  mean.highos.genes.z.top25_mean.lowos.genes.z.top25 20  0.08 66_33
    ##              mean.highos.genes.z_mean.lowos.genes.z  7 0.028 66_33
    ##                                         AGER_CXCL13  1 0.004 66_33
    ##                                         ATP2A3_AGER  1 0.004 66_33
    ##                                     CD79A_RAB11FIP4  1 0.004 66_33
    ##                                         EFNB2_ABHD6  1 0.004 66_33
    ##                                        FAIM3_ATAD3A  1 0.004 66_33
    ##                                            ICK_CD52  1 0.004 66_33
    ##                                          NREP_ABHD6  1 0.004 66_33
    ##                                           TBX3_IDO1  1 0.004 66_33

``` {.r}
print(data.frame(cbind(names(sort(table(all_os_results[all_os_results$test_gene1gene2_continuous_logrank_bonferroni_p<0.05&all_os_results$split=="50_50",]$gene_pair),decreasing=T)[1:10]),
as.numeric(sort(table(all_os_results[all_os_results$test_gene1gene2_continuous_logrank_bonferroni_p<0.05&all_os_results$split=="50_50",]$gene_pair),decreasing=T)[1:10]),
as.numeric(sort(table(all_os_results[all_os_results$test_gene1gene2_continuous_logrank_bonferroni_p<0.05&all_os_results$split=="50_50",]$gene_pair),decreasing=T)[1:10]/250),"50_50"),
stringsAsFactors=F),row.names=F,justify="left")
```

    ##                                                  X1 X2    X3    X4
    ##  mean.highos.genes.z.top25_mean.lowos.genes.z.top25 19 0.076 50_50
    ##              mean.highos.genes.z_mean.lowos.genes.z  6 0.024 50_50
    ##                                     DYNC2H1_PTPRCAP  1 0.004 50_50
    ##                                    FAM118A_CRISPLD1  1 0.004 50_50
    ##                                       FAM118A_FAIM3  1 0.004 50_50
    ##                                      FAM118A_FBXO41  1 0.004 50_50
    ##                                       FAM118A_GFOD1  1 0.004 50_50
    ##                                         FAM118A_KHK  1 0.004 50_50
    ##                                        FAM118A_UBR5  1 0.004 50_50
    ##                                           PLN_FAIM3  1 0.004 50_50

``` {.r}
print("-------------------")
```

    ## [1] "-------------------"

``` {.r}
print(data.frame(cbind(names(sort(table(all_os_bothsig_results[all_os_bothsig_results$split=="80_20",]$gene_pair),decreasing=T)[1:10]),
as.numeric(sort(table(all_os_bothsig_results[all_os_bothsig_results$split=="80_20",]$gene_pair),decreasing=T)[1:10]),
as.numeric(sort(table(all_os_bothsig_results[all_os_bothsig_results$split=="80_20",]$gene_pair),decreasing=T)[1:10]/250),"80_20"),stringsAsFactors=F),row.names=F,justify="left")
```

    ##                                                  X1  X2    X3    X4
    ##              mean.highos.genes.z_mean.lowos.genes.z 119 0.476 80_20
    ##  mean.highos.genes.z.top25_mean.lowos.genes.z.top25  85  0.34 80_20
    ##                                         MAP4K1_AGER  14 0.056 80_20
    ##                                    RAB11FIP4_RAD54B  13 0.052 80_20
    ##                                         TBX3_CXCL13  12 0.048 80_20
    ##                                         MAP4K1_TBX3  11 0.044 80_20
    ##                                         SLAIN1_TBX3  10  0.04 80_20
    ##                                        MAP4K1_USP49   9 0.036 80_20
    ##                                           TBX3_IDO1   9 0.036 80_20
    ##                                        BAHD1_LEPRE1   7 0.028 80_20

``` {.r}
print(data.frame(cbind(names(sort(table(all_os_bothsig_results[all_os_bothsig_results$split=="75_25",]$gene_pair),decreasing=T)[1:10]),
as.numeric(sort(table(all_os_bothsig_results[all_os_bothsig_results$split=="75_25",]$gene_pair),decreasing=T)[1:10]),
as.numeric(sort(table(all_os_bothsig_results[all_os_bothsig_results$split=="75_25",]$gene_pair),decreasing=T)[1:10]/250),"75_25"),stringsAsFactors=F),row.names=F,justify="left")
```

    ##                                                  X1  X2    X3    X4
    ##              mean.highos.genes.z_mean.lowos.genes.z 109 0.436 75_25
    ##  mean.highos.genes.z.top25_mean.lowos.genes.z.top25  75   0.3 75_25
    ##                                         MAP4K1_AGER  18 0.072 75_25
    ##                                         MAP4K1_TBX3  12 0.048 75_25
    ##                                         TBX3_CXCL13  12 0.048 75_25
    ##                                       MAP4K1_DCBLD1   7 0.028 75_25
    ##                                        MAP4K1_USP49   7 0.028 75_25
    ##                                         SLAIN1_TBX3   7 0.028 75_25
    ##                                           TBX3_IDO1   7 0.028 75_25
    ##                                    RAB11FIP4_RAD54B   6 0.024 75_25

``` {.r}
print(data.frame(cbind(names(sort(table(all_os_bothsig_results[all_os_bothsig_results$split=="70_30",]$gene_pair),decreasing=T)[1:10]),
as.numeric(sort(table(all_os_bothsig_results[all_os_bothsig_results$split=="70_30",]$gene_pair),decreasing=T)[1:10]),
as.numeric(sort(table(all_os_bothsig_results[all_os_bothsig_results$split=="70_30",]$gene_pair),decreasing=T)[1:10]/250),"70_30"),stringsAsFactors=F),row.names=F,justify="left")
```

    ##                                                  X1  X2    X3    X4
    ##              mean.highos.genes.z_mean.lowos.genes.z 109 0.436 70_30
    ##  mean.highos.genes.z.top25_mean.lowos.genes.z.top25  95  0.38 70_30
    ##                                         MAP4K1_TBX3   8 0.032 70_30
    ##                                         SLAIN1_TBX3   7 0.028 70_30
    ##                                           TBX3_IDO1   7 0.028 70_30
    ##                                         TBX3_CXCL13   6 0.024 70_30
    ##                                          AGER_FAIM3   5  0.02 70_30
    ##                                          CD79A_AGER   5  0.02 70_30
    ##                                           IDO1_NREP   5  0.02 70_30
    ##                                        MAP4K1_USP49   5  0.02 70_30

``` {.r}
print(data.frame(cbind(names(sort(table(all_os_bothsig_results[all_os_bothsig_results$split=="66_33",]$gene_pair),decreasing=T)[1:10]),
as.numeric(sort(table(all_os_bothsig_results[all_os_bothsig_results$split=="66_33",]$gene_pair),decreasing=T)[1:10]),
as.numeric(sort(table(all_os_bothsig_results[all_os_bothsig_results$split=="66_33",]$gene_pair),decreasing=T)[1:10]/250),"66_33"),stringsAsFactors=F),row.names=F,justify="left")
```

    ##                                                  X1 X2    X3    X4
    ##              mean.highos.genes.z_mean.lowos.genes.z 95  0.38 66_33
    ##  mean.highos.genes.z.top25_mean.lowos.genes.z.top25 88 0.352 66_33
    ##                                         TBX3_CXCL13  7 0.028 66_33
    ##                                         SLAIN1_TBX3  6 0.024 66_33
    ##                                         MAP4K1_TBX3  5  0.02 66_33
    ##                                           TBX3_IDO1  5  0.02 66_33
    ##                                          AGER_FAIM3  4 0.016 66_33
    ##                                          IDO1_USP49  4 0.016 66_33
    ##                                          TBX3_PDZK1  4 0.016 66_33
    ##                                         MAP4K1_AGER  3 0.012 66_33

``` {.r}
print(data.frame(cbind(names(sort(table(all_os_bothsig_results[all_os_bothsig_results$split=="50_50",]$gene_pair),decreasing=T)[1:10]),
as.numeric(sort(table(all_os_bothsig_results[all_os_bothsig_results$split=="50_50",]$gene_pair),decreasing=T)[1:10]),
as.numeric(sort(table(all_os_bothsig_results[all_os_bothsig_results$split=="50_50",]$gene_pair),decreasing=T)[1:10]/250),"50_50"),stringsAsFactors=F),row.names=F,justify="left")
```

    ##                                                  X1  X2    X3    X4
    ##  mean.highos.genes.z.top25_mean.lowos.genes.z.top25 112 0.448 50_50
    ##              mean.highos.genes.z_mean.lowos.genes.z  56 0.224 50_50
    ##                                          AGER_FAIM3   2 0.008 50_50
    ##                                          AGER_RUNX3   2 0.008 50_50
    ##                                          GFOD1_NREP   2 0.008 50_50
    ##                                           GJB2_NREP   2 0.008 50_50
    ##                                          IDO1_USP49   2 0.008 50_50
    ##                                           NREP_CD52   2 0.008 50_50
    ##                                          NREP_FAIM3   2 0.008 50_50
    ##                                        TSPAN9_RUNX3   2 0.008 50_50

``` {.r}
print("-------------------")
```

    ## [1] "-------------------"

``` {.r}
print(data.frame(cbind(names(sort(table(all_os_bothsig_results[all_os_bothsig_results$test_gene1gene2_continuous_logrank_bonferroni_p<0.05&
all_os_bothsig_results$test_gene1gene2_continuous_rnr_roc_bonferroni_p<0.05&all_os_bothsig_results$split=="80_20",]$gene_pair),decreasing=T)[1:10]),
as.numeric(sort(table(all_os_bothsig_results[all_os_bothsig_results$test_gene1gene2_continuous_logrank_bonferroni_p<0.05
&all_os_bothsig_results$test_gene1gene2_continuous_rnr_roc_bonferroni_p<0.05&all_os_bothsig_results$split=="80_20",]$gene_pair),decreasing=T)[1:10]),
as.numeric(sort(table(all_os_bothsig_results[all_os_bothsig_results$test_gene1gene2_continuous_logrank_bonferroni_p<0.05&
all_os_bothsig_results$test_gene1gene2_continuous_rnr_roc_bonferroni_p<0.05&all_os_bothsig_results$split=="80_20",]$gene_pair),decreasing=T)[1:10]/250),"80_20"),
stringsAsFactors=F),row.names=F,justify="left")
```

    ##                                                  X1   X2    X3    X4
    ##  mean.highos.genes.z.top25_mean.lowos.genes.z.top25    7 0.028 80_20
    ##              mean.highos.genes.z_mean.lowos.genes.z    4 0.016 80_20
    ##                                         MAP4K1_AGER    3 0.012 80_20
    ##                                         MAP4K1_TBX3    1 0.004 80_20
    ##                                                <NA> <NA>  <NA> 80_20
    ##                                                <NA> <NA>  <NA> 80_20
    ##                                                <NA> <NA>  <NA> 80_20
    ##                                                <NA> <NA>  <NA> 80_20
    ##                                                <NA> <NA>  <NA> 80_20
    ##                                                <NA> <NA>  <NA> 80_20

``` {.r}
print(data.frame(cbind(names(sort(table(all_os_bothsig_results[all_os_bothsig_results$test_gene1gene2_continuous_logrank_bonferroni_p<0.05&
all_os_bothsig_results$test_gene1gene2_continuous_rnr_roc_bonferroni_p<0.05&all_os_bothsig_results$split=="75_25",]$gene_pair),decreasing=T)[1:10]),
as.numeric(sort(table(all_os_bothsig_results[all_os_bothsig_results$test_gene1gene2_continuous_logrank_bonferroni_p<0.05
&all_os_bothsig_results$test_gene1gene2_continuous_rnr_roc_bonferroni_p<0.05&all_os_bothsig_results$split=="75_25",]$gene_pair),decreasing=T)[1:10]),
as.numeric(sort(table(all_os_bothsig_results[all_os_bothsig_results$test_gene1gene2_continuous_logrank_bonferroni_p<0.05&
all_os_bothsig_results$test_gene1gene2_continuous_rnr_roc_bonferroni_p<0.05&all_os_bothsig_results$split=="75_25",]$gene_pair),decreasing=T)[1:10]/250),"75_25"),
stringsAsFactors=F),row.names=F,justify="left")
```

    ##                                                  X1   X2    X3    X4
    ##              mean.highos.genes.z_mean.lowos.genes.z    9 0.036 75_25
    ##  mean.highos.genes.z.top25_mean.lowos.genes.z.top25    7 0.028 75_25
    ##                                         MAP4K1_AGER    1 0.004 75_25
    ##                                     RAB11FIP4_USP49    1 0.004 75_25
    ##                                                <NA> <NA>  <NA> 75_25
    ##                                                <NA> <NA>  <NA> 75_25
    ##                                                <NA> <NA>  <NA> 75_25
    ##                                                <NA> <NA>  <NA> 75_25
    ##                                                <NA> <NA>  <NA> 75_25
    ##                                                <NA> <NA>  <NA> 75_25

``` {.r}
print(data.frame(cbind(names(sort(table(all_os_bothsig_results[all_os_bothsig_results$test_gene1gene2_continuous_logrank_bonferroni_p<0.05&
all_os_bothsig_results$test_gene1gene2_continuous_rnr_roc_bonferroni_p<0.05&all_os_bothsig_results$split=="70_30",]$gene_pair),decreasing=T)[1:10]),
as.numeric(sort(table(all_os_bothsig_results[all_os_bothsig_results$test_gene1gene2_continuous_logrank_bonferroni_p<0.05
&all_os_bothsig_results$test_gene1gene2_continuous_rnr_roc_bonferroni_p<0.05&all_os_bothsig_results$split=="70_30",]$gene_pair),decreasing=T)[1:10]),
as.numeric(sort(table(all_os_bothsig_results[all_os_bothsig_results$test_gene1gene2_continuous_logrank_bonferroni_p<0.05&
all_os_bothsig_results$test_gene1gene2_continuous_rnr_roc_bonferroni_p<0.05&all_os_bothsig_results$split=="70_30",]$gene_pair),decreasing=T)[1:10]/250),"70_30"),
stringsAsFactors=F),row.names=F,justify="left")
```

    ##                                                  X1   X2    X3    X4
    ##  mean.highos.genes.z.top25_mean.lowos.genes.z.top25    8 0.032 70_30
    ##              mean.highos.genes.z_mean.lowos.genes.z    4 0.016 70_30
    ##                                         MAP4K1_AGER    1 0.004 70_30
    ##                                                <NA> <NA>  <NA> 70_30
    ##                                                <NA> <NA>  <NA> 70_30
    ##                                                <NA> <NA>  <NA> 70_30
    ##                                                <NA> <NA>  <NA> 70_30
    ##                                                <NA> <NA>  <NA> 70_30
    ##                                                <NA> <NA>  <NA> 70_30
    ##                                                <NA> <NA>  <NA> 70_30

``` {.r}
print(data.frame(cbind(names(sort(table(all_os_bothsig_results[all_os_bothsig_results$test_gene1gene2_continuous_logrank_bonferroni_p<0.05&
all_os_bothsig_results$test_gene1gene2_continuous_rnr_roc_bonferroni_p<0.05&all_os_bothsig_results$split=="66_33",]$gene_pair),decreasing=T)[1:10]),
as.numeric(sort(table(all_os_bothsig_results[all_os_bothsig_results$test_gene1gene2_continuous_logrank_bonferroni_p<0.05
&all_os_bothsig_results$test_gene1gene2_continuous_rnr_roc_bonferroni_p<0.05&all_os_bothsig_results$split=="66_33",]$gene_pair),decreasing=T)[1:10]),
as.numeric(sort(table(all_os_bothsig_results[all_os_bothsig_results$test_gene1gene2_continuous_logrank_bonferroni_p<0.05&
all_os_bothsig_results$test_gene1gene2_continuous_rnr_roc_bonferroni_p<0.05&all_os_bothsig_results$split=="66_33",]$gene_pair),decreasing=T)[1:10]/250),"66_33"),
stringsAsFactors=F),row.names=F,justify="left")
```

    ##                                                  X1   X2    X3    X4
    ##  mean.highos.genes.z.top25_mean.lowos.genes.z.top25    9 0.036 66_33
    ##              mean.highos.genes.z_mean.lowos.genes.z    3 0.012 66_33
    ##                                          CIDEB_TBX3    1 0.004 66_33
    ##                                                <NA> <NA>  <NA> 66_33
    ##                                                <NA> <NA>  <NA> 66_33
    ##                                                <NA> <NA>  <NA> 66_33
    ##                                                <NA> <NA>  <NA> 66_33
    ##                                                <NA> <NA>  <NA> 66_33
    ##                                                <NA> <NA>  <NA> 66_33
    ##                                                <NA> <NA>  <NA> 66_33

``` {.r}
print(data.frame(cbind(names(sort(table(all_os_bothsig_results[all_os_bothsig_results$test_gene1gene2_continuous_logrank_bonferroni_p<0.05&
all_os_bothsig_results$test_gene1gene2_continuous_rnr_roc_bonferroni_p<0.05&all_os_bothsig_results$split=="50_50",]$gene_pair),decreasing=T)[1:10]),
as.numeric(sort(table(all_os_bothsig_results[all_os_bothsig_results$test_gene1gene2_continuous_logrank_bonferroni_p<0.05
&all_os_bothsig_results$test_gene1gene2_continuous_rnr_roc_bonferroni_p<0.05&all_os_bothsig_results$split=="50_50",]$gene_pair),decreasing=T)[1:10]),
as.numeric(sort(table(all_os_bothsig_results[all_os_bothsig_results$test_gene1gene2_continuous_logrank_bonferroni_p<0.05&
all_os_bothsig_results$test_gene1gene2_continuous_rnr_roc_bonferroni_p<0.05&all_os_bothsig_results$split=="50_50",]$gene_pair),decreasing=T)[1:10]/250),"50_50"),
stringsAsFactors=F),row.names=F,justify="left")
```

    ##                                                  X1   X2    X3    X4
    ##  mean.highos.genes.z.top25_mean.lowos.genes.z.top25   12 0.048 50_50
    ##              mean.highos.genes.z_mean.lowos.genes.z    1 0.004 50_50
    ##                                                <NA> <NA>  <NA> 50_50
    ##                                                <NA> <NA>  <NA> 50_50
    ##                                                <NA> <NA>  <NA> 50_50
    ##                                                <NA> <NA>  <NA> 50_50
    ##                                                <NA> <NA>  <NA> 50_50
    ##                                                <NA> <NA>  <NA> 50_50
    ##                                                <NA> <NA>  <NA> 50_50
    ##                                                <NA> <NA>  <NA> 50_50

``` {.r}
all_os_results$training_set_size <- unlist(lapply(strsplit(all_os_results$split,"_"),"[",1))
all_os_results$validation_set_size <- unlist(factor(lapply(strsplit(all_os_results$split,"_"),"[",2),levels=c("50","33","30","25","20")))

all_os_bothsig_results$training_set_size <- unlist(lapply(strsplit(all_os_bothsig_results$split,"_"),"[",1))
all_os_bothsig_results$validation_set_size <- unlist(factor(lapply(strsplit(all_os_bothsig_results$split,"_"),"[",2),levels=c("50","33","30","25","20")))

g_bestmod_bar <- ggplot(all_os_results[all_os_results$gene_pair %in% c("MAP4K1_TBX3","MAP4K1_AGER","mean.highos.genes.z.top25_mean.lowos.genes.z.top25","mean.highos.genes.z_mean.lowos.genes.z"),],
aes(x=gene_pair)) + geom_bar(stat="count",fill="lightblue") + facet_grid(~ training_set_size) + theme_classic() +
theme(axis.text.x = element_text(angle = 90)) + stat_count(geom = "text", colour = "black", size = 3.5,aes(label = ..count..),position=position_stack(vjust=0.5)) + ylim(0,250)
ggsave("figure4_outputs/S24g.training_models_barchart.pdf",g_bestmod_bar)
knitr::include_graphics("figure4_outputs/S24g.training_models_barchart.pdf")
```

<embed src="figure4_outputs/S24g.training_models_barchart.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
g_bestmod_bar_bothsig <- ggplot(all_os_bothsig_results[all_os_bothsig_results$gene_pair %in% c("MAP4K1_TBX3","MAP4K1_AGER","mean.highos.genes.z.top25_mean.lowos.genes.z.top25",
"mean.highos.genes.z_mean.lowos.genes.z"),],aes(x=gene_pair)) + geom_bar(stat="count",fill="lightblue") + facet_grid(~ training_set_size) + theme_classic() +
theme(axis.text.x = element_text(angle = 90)) + stat_count(geom = "text", colour = "black", size = 3.5,aes(label = ..count..),position=position_stack(vjust=0.5)) + ylim(0,250)
ggsave("figure4_outputs/S24h.training_models_barchart_bothsig.pdf",g_bestmod_bar_bothsig)
knitr::include_graphics("figure4_outputs/S24h.training_models_barchart_bothsig.pdf")
```

<embed src="figure4_outputs/S24h.training_models_barchart_bothsig.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
g_bestmod_bar_validated <- ggplot(all_os_results[all_os_results$test_gene1gene2_continuous_logrank_bonferroni_p<0.05&all_os_results$gene_pair %in%  c("MAP4K1_TBX3","MAP4K1_AGER",
"mean.highos.genes.z.top25_mean.lowos.genes.z.top25","mean.highos.genes.z_mean.lowos.genes.z"),],aes(x=gene_pair)) + geom_bar(stat="count",fill="lightblue") +
facet_grid(~ validation_set_size) + theme_classic() + theme(axis.text.x = element_text(angle = 90)) +
stat_count(geom = "text", colour = "black", size = 3.5,aes(label = ..count..),position=position_stack(vjust=0.5)) + ylim(0,250) +
scale_x_discrete(limits=c("MAP4K1_AGER","MAP4K1_TBX3","mean.highos.genes.z.top25_mean.lowos.genes.z.top25","mean.highos.genes.z_mean.lowos.genes.z"))
ggsave("figure4_outputs/S24i.validated_models_barchart.pdf",g_bestmod_bar_validated)
knitr::include_graphics("figure4_outputs/S24i.validated_models_barchart.pdf")
```

<embed src="figure4_outputs/S24i.validated_models_barchart.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
g_bestmod_bar_bothsig_validated <- ggplot(all_os_bothsig_results[all_os_bothsig_results$test_gene1gene2_continuous_logrank_bonferroni_p<0.05&
all_os_bothsig_results$test_gene1gene2_continuous_rnr_roc_bonferroni_p<0.05&all_os_bothsig_results$gene_pair %in% c("MAP4K1_TBX3","MAP4K1_AGER","mean.highos.genes.z.top25_mean.lowos.genes.z.top25",
"mean.highos.genes.z_mean.lowos.genes.z"),],aes(x=gene_pair)) + geom_bar(stat="count",fill="lightblue") + facet_grid(~ validation_set_size) + theme_classic() +
theme(axis.text.x = element_text(angle = 90)) + stat_count(geom = "text", colour = "black", size = 3.5,aes(label = ..count..),position=position_stack(vjust=0.5)) + ylim(0,250)
ggsave("figure4_outputs/S24j.validated_models_barchart_bothsig.pdf",g_bestmod_bar_bothsig_validated)
knitr::include_graphics("figure4_outputs/S24j.validated_models_barchart_bothsig.pdf")
```

<embed src="figure4_outputs/S24j.validated_models_barchart_bothsig.pdf" width="600px" height="600px" type="application/pdf" />

Functions

``` {.r}
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

Figure 4 continued

``` {.r}
load("figure4_inputs/val.H.Bayes.RData")
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

gene1.marker <- read.table("figure3_inputs/SKCM.1.markers.txt",as.is=T)$x
gene2.marker <- read.table("figure3_inputs/SKCM.2.markers.txt",as.is=T)$x
gene3.marker <- read.table("figure3_inputs/SKCM.3.markers.txt",as.is=T)$x
gene4.marker <- read.table("figure3_inputs/SKCM.4.markers.txt",as.is=T)$x
gene5.marker <- read.table("figure3_inputs/SKCM.5.markers.txt",as.is=T)$x

cut.NA <- 0.1
OUTPUT <- "figure4_outputs/"
cohort <- "Validation"

master <- skcm_clin
rownames(master) <- master$rna_id
expr <- readRDS("figure4_inputs/val_secondary_180_combat_corrected_log2TPM.rds")
expr[expr==0] <- NA

index.NA <- rowSums(is.na(expr))<= round(cut.NA*ncol(expr))
expr <- expr[index.NA,]
expr.fold <- t(apply(expr,1,function(x) x-median(x,na.rm=T)))

gencode_coding <- read.delim("figure4_inputs/gencode_coding_withentrez.txt",header=T,as.is=T)

x <- c(gene2.marker,gene1.marker,gene3.marker,gene4.marker,gene5.marker) ## marker genes
x <- x[x%in%rownames(expr.fold)]  ### expr.fold=fold-change based expression matrix
x0.new <- expr.fold[match(x,rownames(expr.fold),nomatch=0),ordering.new]
rownames(x0.new) <- gencode_coding$Gene[match(rownames(x0.new),gencode_coding$Gene.id,nomatch=0)]

pdf(file=paste(OUTPUT,paste("S25c",cohort,"marker0.expr.fold.Bayes",0.5,"ordered.full.pdf",sep="."),sep=""),width=12,height=12)
#        p <- plot.expr.fold.heatmap0(x0,1.0,3.0,g.Bayes[ordering])
        p <- plot.expr.fold.heatmap0(x0.new,0.5,3,g.Bayes.new[ordering.new])
        plot(p)
dev.off()
```

    ## pdf 
    ##   2

``` {.r}
knitr::include_graphics("figure4_outputs/S25c.Validation.marker0.expr.fold.Bayes.0.5.ordered.full.pdf")
```

<embed src="figure4_outputs/S25c.Validation.marker0.expr.fold.Bayes.0.5.ordered.full.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
plist_val_proj <- get.sample.association.heatmap(H.Bayes.new,g.Bayes.new,scale=1)
ggsave("figure4_outputs/S25b.Validation_H_norm_newcolors.pdf",plist_val_proj[[2]],height=7,width=7)
knitr::include_graphics("figure4_outputs/S25b.Validation_H_norm_newcolors.pdf")
```

<embed src="figure4_outputs/S25b.Validation_H_norm_newcolors.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
clin$study <- "CPB"

tdf <- reshape2::melt(table(c(as.character(clin$plot_cluster),as.character(skcm_clin$plot_cluster)),c(clin$study,skcm_clin$cohort)))
tmp_ft <- chisq.test(table(c(clin$study,skcm_clin$cohort),c(as.character(clin$plot_cluster),as.character(skcm_clin$plot_cluster))))
tmp_ft_pstring <- paste0("Chi-Sqaured test p=",formatC(tmp_ft$p.value,digits=2))
g_subtype_cohorts <- ggplot(tdf) + geom_mosaic(aes(weight=value,x=product(Var2),fill=Var1),offset=0.01) + scale_fill_manual(values=c(cbpalette_reordered[4],cbpalette_reordered[2],
cbpalette_reordered[5],cbpalette_reordered[1],cbpalette_reordered[3])) + theme_bw() +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.title=element_blank(),axis.text.x=element_text(size=12),axis.text.y=element_text(size=12)) +
xlab("") + ylab("") + annotate(geom="text",x=0.5,y=1.02,label=tmp_ft_pstring) + scale_y_continuous(breaks=seq(0,1,by=0.25))
ggsave("figure4_outputs/S25d.Validation_CPB_subtype_redblue_GideLiu.pdf",g_subtype_cohorts + geom_text(data = ggplot_build(g_subtype_cohorts)$data[[1]], aes(x = (xmin+xmax)/2, y = (ymin+ymax)/2, label=.wt)))
knitr::include_graphics("figure4_outputs/S25d.Validation_CPB_subtype_redblue_GideLiu.pdf")
```

<embed src="figure4_outputs/S25d.Validation_CPB_subtype_redblue_GideLiu.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
tmpkt <- kruskal.test(skcm_clin$rna_tcb,skcm_clin$plot_cluster)
g_tcb_cluster <- ggplot(skcm_clin,aes(x=plot_cluster,y=rna_tcb)) + geom_boxplot(aes(fill= plot_cluster),outlier.shape=NA) + geom_jitter() + theme_classic() + ylab(expression(TCB[RNA])) + xlab("") +
scale_fill_manual(values=cbpalette_reordered) + annotate(geom="text",x=2.5,y=10,label=paste0("Kruskal-Wallis p = ",formatC(tmpkt$p.value,digits=2,format="e"))) + scale_y_log10()
ggsave("figure4_outputs/S25e.Validation_rna_tcb_cluster.pdf", g_tcb_cluster)
knitr::include_graphics("figure4_outputs/S25e.Validation_rna_tcb_cluster.pdf")
```

<embed src="figure4_outputs/S25e.Validation_rna_tcb_cluster.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
tmpkt <- kruskal.test(skcm_clin$rna_bcb,skcm_clin$plot_cluster)
g_bcb_cluster <- ggplot(skcm_clin,aes(x=plot_cluster,y=rna_bcb)) + geom_boxplot(aes(fill= plot_cluster),outlier.shape=NA) + geom_jitter() + theme_classic() + ylab(expression(BCB[RNA])) + xlab("") +
scale_fill_manual(values=cbpalette_reordered) + annotate(geom="text",x=2.5,y=1000,label=paste0("Kruskal-Wallis p = ",formatC(tmpkt$p.value,digits=2,format="e"))) + scale_y_log10()
ggsave("figure4_outputs/S25f.Validation_rna_bcb_cluster.pdf", g_bcb_cluster)
knitr::include_graphics("figure4_outputs/S25f.Validation_rna_bcb_cluster.pdf")
```

<embed src="figure4_outputs/S25f.Validation_rna_bcb_cluster.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
skcm_clin$class <- skcm_clin$plot_cluster
tmp_cph <- summary(coxph(surv ~ class, data = skcm_clin))
skcm_nmf_cluster_p_string <- paste0(" Log-rank p=",as.character(formatC(tmp_cph$sctest["pvalue"],digits=2,format="e")))
pdf("figure4_outputs/S25g.Validation_nmf_cluster_survival.pdf",onefile=F)
fit<- survfit(surv ~ class, data = skcm_clin)
ggsurvplot(fit, data = skcm_clin,pval=TRUE,pval.method=TRUE,legend=c(0.75,0.9),censor=TRUE,palette=cbpalette_reordered[1:5],title=skcm_nmf_cluster_p_string)
dev.off()
```

    ## pdf 
    ##   2

``` {.r}
knitr::include_graphics("figure4_outputs/S25g.Validation_nmf_cluster_survival.pdf")
```

<embed src="figure4_outputs/S25g.Validation_nmf_cluster_survival.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
tdf <- reshape2::melt(table(skcm_clin$Response,skcm_clin$plot_cluster))
tmp_ft <- fisher.test(table(skcm_clin$Response,skcm_clin$plot_cluster))
tmp_ft_pstring <- paste0("Fisher p=",formatC(tmp_ft$p.value,digits=2))
tmp_ft_orstring <- paste0("OR: ",formatC(tmp_ft$estimate,digits=2),", 95% CI ",formatC(tmp_ft$conf.int[1],digits=2),"-",formatC(tmp_ft$conf.int[2],digits=2))
g_rnr_subtype <- ggplot(tdf) + geom_mosaic(aes(weight=value,x=product(Var2),fill=Var1),offset=0.01) + scale_fill_manual(values=c("red","blue")) + theme_bw() +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.title=element_blank(),axis.text.x=element_text(size=12),axis.text.y=element_text(size=12)) +
xlab("") + ylab("") + annotate(geom="text",x=0.5,y=1.02,label=tmp_ft_pstring) + scale_y_continuous(breaks=seq(0,1,by=0.25))
ggsave("figure4_outputs/S25h.subtype_Validation_response_redblue.pdf",g_rnr_subtype + geom_text(data = ggplot_build(g_rnr_subtype)$data[[1]], aes(x = (xmin+xmax)/2, y = (ymin+ymax)/2, label=.wt)))
knitr::include_graphics("figure4_outputs/S25h.subtype_Validation_response_redblue.pdf")
```

<embed src="figure4_outputs/S25h.subtype_Validation_response_redblue.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
skcm_clin$plot_cluster_immuneonly <- ifelse(skcm_clin$plot_cluster=="Immune","Immune","Others")
skcm_immune_hr <- summary(coxph(surv ~ plot_cluster_immuneonly,data=skcm_clin))$coef[1,2]
skcm_immune_hr_lo <- summary(coxph(surv ~ plot_cluster_immuneonly,data=skcm_clin))$conf.int[1,3]
skcm_immune_hr_hi <- summary(coxph(surv ~ plot_cluster_immuneonly,data=skcm_clin))$conf.int[1,4]
skcm_clin$class <- skcm_clin$plot_cluster_immuneonly
class_table <- table(skcm_clin$class)
skcm_immune_hr_string <- paste0("HR=",formatC(skcm_immune_hr,digits=3)," (95% CI, ",formatC(skcm_immune_hr_lo,digits=3),"-",formatC(skcm_immune_hr_hi,digits=3),")\n",
names(class_table)[1],"=",class_table[[1]]," ",names(class_table)[2],"=",class_table[[2]])
pdf("figure4_outputs/4f.Validation_Immuneonly_survival_bestrest.pdf",onefile=F)
fit<- survfit(surv ~ class, data = skcm_clin)
ggsurvplot(fit, data = skcm_clin,pval=TRUE,pval.method=TRUE,legend=c(0.75,0.90),censor=TRUE,palette=c("blue","red"),title=skcm_immune_hr_string)
dev.off()
```

    ## pdf 
    ##   2

``` {.r}
knitr::include_graphics("figure4_outputs/4f.Validation_Immuneonly_survival_bestrest.pdf")
```

<embed src="figure4_outputs/4f.Validation_Immuneonly_survival_bestrest.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
tmp_cph <- coxph(surv ~ MAP4K1 + TBX3,data=clin[!is.na(clin$surv),])
tmp_cpb_ggforest_MAP4K1_TBX3 <- ggforest(tmp_cph,data=clin[!is.na(clin$surv),],fontsize=0.9)
ggsave("figure4_outputs/4i.CPB_MAP4K1TBX3_forestplot.pdf",tmp_cpb_ggforest_MAP4K1_TBX3,height=3,width=6)
knitr::include_graphics("figure4_outputs/4i.CPB_MAP4K1TBX3_forestplot.pdf")
```

<embed src="figure4_outputs/4i.CPB_MAP4K1TBX3_forestplot.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
tmp_cph <- coxph(surv ~ MAP4K1 + TBX3,data=skcm_clin[!is.na(skcm_clin$surv),])
tmp_val_ggforest_MAP4K1_TBX3 <- ggforest(tmp_cph,data=skcm_clin[!is.na(skcm_clin$surv),],fontsize=0.9)
ggsave("figure4_outputs/4j.Val_MAP4K1TBX3_forestplot.pdf",tmp_val_ggforest_MAP4K1_TBX3,height=3,width=6)
knitr::include_graphics("figure4_outputs/4j.Val_MAP4K1TBX3_forestplot.pdf")
```

<embed src="figure4_outputs/4j.Val_MAP4K1TBX3_forestplot.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
clin$tmp_pred_MAP4K1_TBX3 <- NA
clin[!is.na(clin$surv),]$tmp_pred_MAP4K1_TBX3 <- predict(coxph(surv ~ MAP4K1 + TBX3,data=clin[!is.na(clin$surv),]))
clin$tmp_pred_MAP4K1_TBX3_cut2 <- cut(clin$tmp_pred_MAP4K1_TBX3, breaks=c(quantile(clin$tmp_pred_MAP4K1_TBX3, probs = c(0,1/2,1),na.rm=T)),labels=c("low risk","high risk"),na.rm=T,include.lowest=T)

clin$class <- clin$tmp_pred_MAP4K1_TBX3_cut2
tmp_cph <- summary(coxph(surv ~ tmp_pred_MAP4K1_TBX3_cut2,data=clin[!is.na(clin$surv),]))
cpb_pred_hr <- tmp_cph$coef[1,2]
cpb_pred_hr_lo <- tmp_cph$conf.int[1,3]
cpb_pred_hr_hi <- tmp_cph$conf.int[1,4]
class_table <- table(clin[!is.na(clin$surv),]$class)
cpb_map4k1tbx3_hr_string <- paste0("Primary: ","HR=",formatC(cpb_pred_hr,digits=3)," (95% CI, ",formatC(cpb_pred_hr_lo,digits=3),"-",formatC(cpb_pred_hr_hi,digits=3),"\n",
" Log-rank p=",as.character(formatC(tmp_cph$sctest["pvalue"],digits=2,format="e")),names(class_table)[1],"=",class_table[[1]]," ",names(class_table)[2],"=",class_table[[2]])
pdf("figure4_outputs/4k.CPB_MAP4K1_TBX3_cut2.pdf",onefile=F)
fit<- survfit(surv ~ tmp_pred_MAP4K1_TBX3_cut2, data = clin[!is.na(clin$surv),])
ggsurvplot(fit, data = clin[!is.na(clin$surv),],pval=TRUE,pval.method=TRUE,legend=c(0.75,0.90),censor=TRUE,palette=c("blue","red"),legend.labs=c("Low risk", "High risk"),
legend.title="MAP4K1+TBX3 model risk",title=cpb_map4k1tbx3_hr_string,font.legend=12)
dev.off()
```

    ## pdf 
    ##   2

``` {.r}
knitr::include_graphics("figure4_outputs/4k.CPB_MAP4K1_TBX3_cut2.pdf")
```

<embed src="figure4_outputs/4k.CPB_MAP4K1_TBX3_cut2.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
skcm_clin$tmp_pred_MAP4K1_TBX3 <- predict(coxph(surv ~ MAP4K1 + TBX3,data=skcm_clin))
skcm_clin$tmp_pred_MAP4K1_TBX3_cut2 <- cut(skcm_clin$tmp_pred_MAP4K1_TBX3, breaks=c(quantile(skcm_clin$tmp_pred_MAP4K1_TBX3, probs = c(0,1/2,1),na.rm=T)),labels=c("low risk","high risk"),na.rm=T,
include.lowest=T)

skcm_clin$class <- skcm_clin$tmp_pred_MAP4K1_TBX3_cut2
tmp_cph <- summary(coxph(surv ~ tmp_pred_MAP4K1_TBX3_cut2,data=skcm_clin))
val_pred_hr <- tmp_cph$coef[1,2]
val_pred_hr_lo <- tmp_cph$conf.int[1,3]
val_pred_hr_hi <- tmp_cph$conf.int[1,4]
class_table <- table(skcm_clin$class)
val_map4k1tbx3_hr_string <- paste0("Secondary: ","HR=",formatC(val_pred_hr,digits=3)," (95% CI, ",formatC(val_pred_hr_lo,digits=3),"-",formatC(val_pred_hr_hi,digits=3),"\n",
" Log-rank p=",as.character(formatC(tmp_cph$sctest["pvalue"],digits=2,format="e")),names(class_table)[1],"=",class_table[[1]]," ",names(class_table)[2],"=",class_table[[2]])
pdf("figure4_outputs/4l.Val_MAP4K1_TBX3_cut2.pdf",onefile=F)
fit<- survfit(surv ~ tmp_pred_MAP4K1_TBX3_cut2, data = skcm_clin)
ggsurvplot(fit, data = skcm_clin,pval=TRUE,pval.method=TRUE,legend=c(0.75,0.90),censor=TRUE,palette=c("blue","red"),legend.labs=c("Low risk", "High risk"),
legend.title="MAP4K1+TBX3 model risk",title=val_map4k1tbx3_hr_string,font.legend=12)
dev.off()
```

    ## pdf 
    ##   2

``` {.r}
knitr::include_graphics("figure4_outputs/4l.Val_MAP4K1_TBX3_cut2.pdf")
```

<embed src="figure4_outputs/4l.Val_MAP4K1_TBX3_cut2.pdf" width="600px" height="600px" type="application/pdf" />


``` {.r}
val_raw_heatmap_annotation_colors <- list(
                              OS=c("OS >= 1yr"="blue","OS < 1yr"="red"),
                              Response=c("NR"="red","R"="blue"),
                              plot_cluster=c("MITF low"="#56B4E9","Immune"="#009E73","Intermediate"="#999999","Keratin high"="#F0E442","MITF high"="#E69F00"),
                              cohort=c("Gide"="gray","Liu"="black")
)

skcm_clin$OS <- ifelse(skcm_clin$overall_survival>=365,"OS >= 1yr","OS < 1yr")

skcm_clin$log10_rna_tcb <- log10(skcm_clin$rna_tcb)
skcm_clin$scale_log10_rna_tcb <- scale(skcm_clin$log10_rna_tcb)
skcm_clin$log10_rna_bcb <- log10(skcm_clin$rna_bcb)
skcm_clin$scale_log10_rna_bcb <- scale(skcm_clin$log10_rna_bcb)
skcm_clin$scale_TIDE_Score <- scale(skcm_clin$TIDE_Score)
skcm_clin$MHC_II_zscore <- scale(skcm_clin$MHC_II)
skcm_clin$scale_CYT <- scale(skcm_clin$CYT)
skcm_clin$scale_GEP <- scale(skcm_clin$GEP)
skcm_clin$scale_CD274 <- scale(skcm_clin$CD274)
skcm_clin$scale_IMPRES <- scale(skcm_clin$IMPRES)

tmp_val_rawdata_breaks_9 <- c(quantile_breaks(unlist(skcm_clin[,c("scale_log10_rna_tcb","scale_log10_rna_bcb","scale_TIDE_Score","MHC_II_zscore","scale_CYT","scale_GEP","scale_CD274","scale_IMPRES",
"MAP4K1_z","TBX3_z","AGER_z","mean_lowos_genes_z","mean_highos_genes_z")]),10))

secondary_heatmap <- pheatmap(t(skcm_clin[order(skcm_clin$overall_survival),
c("scale_log10_rna_tcb","scale_log10_rna_bcb","scale_TIDE_Score","MHC_II_zscore","scale_CYT","scale_GEP","scale_CD274","scale_IMPRES","MAP4K1_z","TBX3_z","AGER_z",
"mean_lowos_genes_z","mean_highos_genes_z")]),
annotation_col=skcm_clin[,c("cohort","Response","OS","plot_cluster")],
cluster_rows=TRUE,cluster_cols=TRUE,clustering_distance_rows="correlation",fontsize_row=5, fontsize_col=3,cellwidth=3,cellheight=5,fontsize=5,na_col="grey",
breaks=tmp_val_rawdata_breaks_9,color=rev(brewer.pal(n = 9, name = "RdBu")),treeheight_row=10,treeheight_col=10,annotation_colors=val_raw_heatmap_annotation_colors,
legend_breaks=c(-4,-3,-2,-1,0,1,2,3,4))
save_pheatmap_pdf(secondary_heatmap, "figure4_outputs/S26a.Secondary_OSgenes_rawdata_heatmap_clustered_9col.pdf")
```

    ## pdf 
    ##   2

``` {.r}
knitr::include_graphics("figure4_outputs/S26a.Secondary_OSgenes_rawdata_heatmap_clustered_9col.pdf")
```

<embed src="figure4_outputs/S26a.Secondary_OSgenes_rawdata_heatmap_clustered_9col.pdf" width="700px" height="700px" type="application/pdf" />


``` {.r}
tmp_cph <- coxph(surv ~ mean_highos_genes_z + mean_lowos_genes_z,data=skcm_clin[!is.na(skcm_clin$surv),])
tmp_val_ggforest_mean_highos_z_mean_lowos_z <- ggforest(tmp_cph,data=skcm_clin[!is.na(skcm_clin$surv),],fontsize=0.9)
ggsave("figure4_outputs/S26e.Val_mean_highos_z_mean_lowos_z_forestplot.pdf",tmp_val_ggforest_mean_highos_z_mean_lowos_z,height=3,width=9)
knitr::include_graphics("figure4_outputs/S26e.Val_mean_highos_z_mean_lowos_z_forestplot.pdf")
```

<embed src="figure4_outputs/S26e.Val_mean_highos_z_mean_lowos_z_forestplot.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
tmp_cph <- coxph(surv ~ MAP4K1 + AGER,data=skcm_clin[!is.na(skcm_clin$surv),])
tmp_val_ggforest_MAP4K1_AGER <- ggforest(tmp_cph,data=skcm_clin[!is.na(skcm_clin$surv),],fontsize=0.9)
ggsave("figure4_outputs/S26f.Val_MAP4K1AGER_forestplot.pdf",tmp_val_ggforest_MAP4K1_AGER,height=3,width=9)
knitr::include_graphics("figure4_outputs/S26f.Val_MAP4K1AGER_forestplot.pdf")
```

<embed src="figure4_outputs/S26f.Val_MAP4K1AGER_forestplot.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
cpb_MAP4K1_vs_MAP4K1TBX3 <- formatC(anova(coxph(surv ~ MAP4K1,data=clin),coxph(surv ~ MAP4K1 + TBX3,data=clin))$"P(>|Chi|)"[2],format="e",digits=2)
val_MAP4K1_vs_MAP4K1TBX3 <- formatC(anova(coxph(surv ~ MAP4K1 + treatment_type,data=skcm_clin),coxph(surv ~ MAP4K1 + TBX3 + treatment_type,data=skcm_clin))$"P(>|Chi|)"[2],format="e",digits=2)
cpb_TBX3_vs_MAP4K1TBX3 <- formatC(anova(coxph(surv ~ TBX3,data=clin),coxph(surv ~ MAP4K1 + TBX3,data=clin))$"P(>|Chi|)"[2],format="e",digits=2)
val_TBX3_vs_MAP4K1TBX3 <- formatC(anova(coxph(surv ~ TBX3 + treatment_type,data=skcm_clin),coxph(surv ~ MAP4K1 + TBX3 + treatment_type,data=skcm_clin))$"P(>|Chi|)"[2],format="e",digits=2)

df_MAP4K1TBX3_LRT_p <- data.frame(cohort=c("Primary","Secondary"),model_plusTBX3=c(cpb_MAP4K1_vs_MAP4K1TBX3,val_MAP4K1_vs_MAP4K1TBX3),TBX3_plus_model=c(cpb_TBX3_vs_MAP4K1TBX3,val_TBX3_vs_MAP4K1TBX3),
stringsAsFactors=F)
colnames(df_MAP4K1TBX3_LRT_p) <- c("cohort","Likelihood Ratio Test p\nfor MAP4K1 + TBX3 model\nvs. MAP4K1 model","Likelihood Ratio Test p\nfor MAP4K1 + TBX3 model\nvs. TBX3 model")
pdf("figure4_outputs/S26g.Model_comparison_MAP4K1_TBX3_likelihoodratiotests_bothorders.pdf")
grid.table(df_MAP4K1TBX3_LRT_p, rows=NULL)
dev.off()
```

    ## pdf 
    ##   2

``` {.r}
knitr::include_graphics("figure4_outputs/S26g.Model_comparison_MAP4K1_TBX3_likelihoodratiotests_bothorders.pdf")
```

<embed src="figure4_outputs/S26g.Model_comparison_MAP4K1_TBX3_likelihoodratiotests_bothorders.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
clin$MAP4K1_TBX3_categ <- paste0(ifelse(clin$TBX3>median(clin$TBX3),"High TBX3","Low TBX3"),", ",ifelse(clin$MAP4K1>median(clin$MAP4K1),"High MAP4K1","Low MAP4K1"))
skcm_clin$MAP4K1_TBX3_categ <- paste0(ifelse(skcm_clin$TBX3>median(skcm_clin$TBX3),"High TBX3","Low TBX3"),", ",ifelse(skcm_clin$MAP4K1>median(skcm_clin$MAP4K1),"High MAP4K1","Low MAP4K1"))

clin$class <- clin$MAP4K1_TBX3_categ
class_table <- table(clin[!is.na(clin$surv),]$class)
tmp_cph <- summary(coxph(surv ~ class,data=clin[!is.na(clin$surv),]))
cpb_map4k1tbx3_binary_string <- paste0("Log-rank p=",as.character(formatC(tmp_cph$sctest["pvalue"],digits=2,format="e")),"\n",
names(class_table)[1],"=",class_table[[1]]," ",names(class_table)[2],"=",class_table[[2]],"\n",
names(class_table)[3],"=",class_table[[3]]," ",names(class_table)[4],"=",class_table[[4]])
pdf("figure4_outputs/S26h.CPB_MAP4K1_TBX3_highlow.pdf",onefile=F)
fit<- survfit(surv ~ MAP4K1_TBX3_categ, data = clin)
ggsurvplot(fit, data = clin,pval=TRUE,pval.method=TRUE,legend=c(0.75,0.90),censor=TRUE,palette=c("blue","red","green","purple"),legend.title="MAP4K1+TBX3 binary models",
legend.labs=c("High TBX3, High MAP4K1", "High TBX3, Low MAP4K1","Low TBX3, High MAP4K1","Low TBX3, Low MAP4K1"),
title=cpb_map4k1tbx3_binary_string,font.legend=12)
dev.off()
```

    ## pdf 
    ##   2

``` {.r}
knitr::include_graphics("figure4_outputs/S26h.CPB_MAP4K1_TBX3_highlow.pdf")
```

<embed src="figure4_outputs/S26h.CPB_MAP4K1_TBX3_highlow.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
skcm_clin$class <- skcm_clin$MAP4K1_TBX3_categ
class_table <- table(skcm_clin$class)
tmp_cph <- summary(coxph(surv ~ class,data=skcm_clin))
val_map4k1tbx3_binary_string <- paste0("Log-rank p=",as.character(formatC(tmp_cph$sctest["pvalue"],digits=2,format="e")),"\n",
names(class_table)[1],"=",class_table[[1]]," ",names(class_table)[2],"=",class_table[[2]],"\n",
names(class_table)[3],"=",class_table[[3]]," ",names(class_table)[4],"=",class_table[[4]])
pdf("figure4_outputs/S26i.Val_MAP4K1_TBX3_highlow.pdf",onefile=F)
fit<- survfit(surv ~ MAP4K1_TBX3_categ, data = skcm_clin)
ggsurvplot(fit, data = skcm_clin,pval=TRUE,pval.method=TRUE,legend=c(0.75,0.90),censor=TRUE,palette=c("blue","red","green","purple"),legend.title="MAP4K1+TBX3 binary models",
legend.labs=c("High TBX3, High MAP4K1", "High TBX3, Low MAP4K1","Low TBX3, High MAP4K1","Low TBX3, Low MAP4K1"),
title=val_map4k1tbx3_binary_string,font.legend=12)
dev.off()
```

    ## pdf 
    ##   2

``` {.r}
knitr::include_graphics("figure4_outputs/S26i.Val_MAP4K1_TBX3_highlow.pdf")
```

<embed src="figure4_outputs/S26i.Val_MAP4K1_TBX3_highlow.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
skcm_therapy_hr <- summary(coxph(surv ~ treatment_type=="PD1",data=skcm_clin))$coef[1,2] #1.62
skcm_therapy_hr_lo <- summary(coxph(surv ~ treatment_type=="PD1",data=skcm_clin))$conf.int[1,3] #1.12
skcm_therapy_hr_hi <- summary(coxph(surv ~ treatment_type=="PD1",data=skcm_clin))$conf.int[1,4] #2.33
skcm_clin$class <- skcm_clin$treatment_type
class_table <- table(skcm_clin$class)
skcm_therapy_hr_string <- paste0("HR=",formatC(skcm_therapy_hr,digits=3)," (95% CI, ",formatC(skcm_therapy_hr_lo,digits=3),"-",formatC(skcm_therapy_hr_hi,digits=3),")","\n",
names(class_table)[1],"=",class_table[[1]]," ",names(class_table)[2],"=",class_table[[2]])
pdf("figure4_outputs/S26jValidation_therapy_survival.pdf",onefile=F)
fit<- survfit(surv ~ class, data = skcm_clin)
ggsurvplot(fit, data = skcm_clin,pval=TRUE,pval.method=TRUE,legend=c(0.75,0.90),censor=TRUE,palette=c("blue","red"),title=skcm_therapy_hr_string)
dev.off()
```

    ## pdf 
    ##   2

``` {.r}
knitr::include_graphics("figure4_outputs/S26jValidation_therapy_survival.pdf")
```

<embed src="figure4_outputs/S26jValidation_therapy_survival.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
tmp_cph <- coxph(surv ~ MAP4K1 + TBX3 + treatment_type,data=skcm_clin)
tmp_ggforest_MAP4K1_TBX3_val_withtreatment <- ggforest(tmp_cph,data=skcm_clin,fontsize=0.9)
ggsave("figure4_outputs/S26k.Validation_MAP4K1_TBX3_forestplot_withtreatment.pdf",tmp_ggforest_MAP4K1_TBX3_val_withtreatment,height=4,width=7)
knitr::include_graphics("figure4_outputs/S26k.Validation_MAP4K1_TBX3_forestplot_withtreatment.pdf")
```

<embed src="figure4_outputs/S26k.Validation_MAP4K1_TBX3_forestplot_withtreatment.pdf" width="600px" height="600px" type="application/pdf" />

Functions

``` {.r}
plot_multicohort_survival <- function(df1,df2,var1,var1_string,pdf_string,ymax_label_shift,var1_y_breaks) {
#plot_multicohort_survival(clin,skcm_clin,"CD274","CD274","S27a.multicohort_CD274_HR_p.pdf",0.5,c(0.5,0.75,1,1.5,2,3,4,5))

#print(pdf_string)

df1$var1 <- df1[,var1]
df2$var1 <- df2[,var1]

cph_var1_mgh <- summary(coxph(surv ~ var1,data=df1[df1$cohort=="MGH",]))
cph_var1_hug <- summary(coxph(surv ~ var1,data=df1[df1$cohort=="Hugo"&!is.na(df1$surv),]))
cph_var1_ria <- summary(coxph(surv ~ var1,data=df1[df1$cohort=="Riaz",]))
cph_var1_van <- summary(coxph(surv ~ var1,data=df1[df1$cohort=="VanAllen",]))
cph_var1_aus <- summary(coxph(surv ~ var1 + treatment_type,data=df2[df2$cohort=="Gide",]))
cph_var1_liu <- summary(coxph(surv ~ var1,data=df2[df2$cohort=="Liu",]))

cph_var1_mgh_hr <-    cph_var1_mgh$conf.int["var1",1]
cph_var1_mgh_hr_lo <- cph_var1_mgh$conf.int["var1",3]
cph_var1_mgh_hr_hi <- cph_var1_mgh$conf.int["var1",4]
cph_var1_mgh_hr_p <-  cph_var1_mgh$coef["var1",5]
cph_var1_mgh_lr_p <-  cph_var1_mgh$sctest[3]
cph_var1_mgh_hr_se <- cph_var1_mgh$coef["var1",3]

cph_var1_hug_hr <-    cph_var1_hug$conf.int["var1",1]
cph_var1_hug_hr_lo <- cph_var1_hug$conf.int["var1",3]
cph_var1_hug_hr_hi <- cph_var1_hug$conf.int["var1",4]
cph_var1_hug_hr_p <-  cph_var1_hug$coef["var1",5]
cph_var1_hug_lr_p <-  cph_var1_hug$sctest[3]
cph_var1_hug_hr_se <- cph_var1_hug$coef["var1",3]

cph_var1_ria_hr <-    cph_var1_ria$conf.int["var1",1]
cph_var1_ria_hr_lo <- cph_var1_ria$conf.int["var1",3]
cph_var1_ria_hr_hi <- cph_var1_ria$conf.int["var1",4]
cph_var1_ria_hr_p <-  cph_var1_ria$coef["var1",5]
cph_var1_ria_lr_p <-  cph_var1_ria$sctest[3]
cph_var1_ria_hr_se <- cph_var1_ria$coef["var1",3]

cph_var1_van_hr <-    cph_var1_van$conf.int["var1",1]
cph_var1_van_hr_lo <- cph_var1_van$conf.int["var1",3]
cph_var1_van_hr_hi <- cph_var1_van$conf.int["var1",4]
cph_var1_van_hr_p <-  cph_var1_van$coef["var1",5]
cph_var1_van_lr_p <-  cph_var1_van$sctest[3]
cph_var1_van_hr_se <- cph_var1_van$coef["var1",3]

cph_var1_aus_hr <-    cph_var1_aus$conf.int["var1",1]
cph_var1_aus_hr_lo <- cph_var1_aus$conf.int["var1",3]
cph_var1_aus_hr_hi <- cph_var1_aus$conf.int["var1",4]
cph_var1_aus_hr_p <-  cph_var1_aus$coef["var1",5]
cph_var1_aus_lr_p <-  cph_var1_aus$sctest[3]
cph_var1_aus_hr_se <- cph_var1_aus$coef["var1",3]

cph_var1_liu_hr <-    cph_var1_liu$conf.int["var1",1]
cph_var1_liu_hr_lo <- cph_var1_liu$conf.int["var1",3]
cph_var1_liu_hr_hi <- cph_var1_liu$conf.int["var1",4]
cph_var1_liu_hr_p <-  cph_var1_liu$coef["var1",5]
cph_var1_liu_lr_p <-  cph_var1_liu$sctest[3]
cph_var1_liu_hr_se <- cph_var1_liu$coef["var1",3]

n_mgh <- nrow(df1[df1$cohort=="MGH",])
n_hug <- nrow(df1[df1$cohort=="Hugo"&!is.na(df1$surv),])
n_ria <- nrow(df1[df1$cohort=="Riaz",])
n_van <- nrow(df1[df1$cohort=="VanAllen",])
n_aus <- nrow(df2[df2$cohort=="Gide",])
n_liu <- nrow(df2[df2$cohort=="Liu",])

df_cph_var1 <- data.frame(rbind(c("MGH",paste0("MGH ",var1_string),cph_var1_mgh_hr,cph_var1_mgh_hr_lo,cph_var1_mgh_hr_hi,cph_var1_mgh_hr_p,cph_var1_mgh_lr_p,n_mgh),
                                      c("Hugo",paste0("Hugo ",var1_string),cph_var1_hug_hr,cph_var1_hug_hr_lo,cph_var1_hug_hr_hi,cph_var1_hug_hr_p,cph_var1_hug_lr_p,n_hug),
                                      c("Riaz",paste0("Riaz ",var1_string),cph_var1_ria_hr,cph_var1_ria_hr_lo,cph_var1_ria_hr_hi,cph_var1_ria_hr_p,cph_var1_ria_lr_p,n_ria),
                                      c("Van Allen",paste0("Van Allen ",var1_string),cph_var1_van_hr,cph_var1_van_hr_lo,cph_var1_van_hr_hi,cph_var1_van_hr_p,cph_var1_van_lr_p,n_van),
                                      c("Gide",paste0("Gide ",var1_string),cph_var1_aus_hr,cph_var1_aus_hr_lo,cph_var1_aus_hr_hi,cph_var1_aus_hr_p,cph_var1_aus_lr_p,n_aus),
                                      c("Liu",paste0("Liu ",var1_string),cph_var1_liu_hr,cph_var1_liu_hr_lo,cph_var1_liu_hr_hi,cph_var1_liu_hr_p,cph_var1_liu_lr_p,n_liu)),
                                stringsAsFactors=F)


colnames(df_cph_var1) <- c("cohort","name","HR","HR_lo","HR_hi","HR_p","LR_p","cohort_n")

df_cph_var1$print_lrt_p <- paste0("p = ",ifelse(df_cph_var1$LR_p<0.01,formatC(as.numeric(df_cph_var1$LR_p),width=3,format="e",digits=2),formatC(as.numeric(df_cph_var1$LR_p),width=3,digits=2)))
df_cph_var1$print_hr_p <- paste0("p = ",ifelse(df_cph_var1$HR_p<0.01,formatC(as.numeric(df_cph_var1$HR_p),width=3,format="e",digits=2),formatC(as.numeric(df_cph_var1$HR_p),width=3,digits=2)))

df_cph_var1$cohort_with_n <- paste0(df_cph_var1$cohort,"\n","(n=",df_cph_var1$cohort_n,")")
df_cph_var1$cohort_with_n <- factor(df_cph_var1$cohort_with_n,levels=c("Gide\n(n=73)","Hugo\n(n=25)","Liu\n(n=107)","MGH\n(n=38)","Riaz\n(n=51)","Van Allen\n(n=39)"))

g_var1_allmeta_hr_p <- ggplot(df_cph_var1,aes(x=cohort_with_n,y=as.numeric(HR))) +
  geom_point(aes(size=as.numeric(cohort_n))) + geom_pointrange(aes(ymin=as.numeric(HR_lo),ymax=as.numeric(HR_hi))) + theme_classic() +
  xlab("") + ylab("Hazard Ratio") + geom_text(aes(x=cohort_with_n,y=max(as.numeric(df_cph_var1$HR_hi))+ymax_label_shift,label=print_hr_p),size=3.5) +
  geom_hline(yintercept=1,linetype="dashed") + ggtitle(var1_string) +
  theme(legend.position="bottom",plot.title=element_text(hjust=0.5,size=16),axis.text.x=element_text(size=16,color="black"),axis.text.y=element_text(size=14,color="black"),
  axis.title.y=element_text(size=16)) + coord_trans(y="log2") + scale_color_manual(values=c("red","blue")) + scale_y_continuous(breaks=var1_y_breaks)
ggsave(pdf_string,g_var1_allmeta_hr_p,height=7,width=8)

return(df_cph_var1)
}

plot_multicohort_survival_double <- function(df1,df2,var1,var1_string,var2,var2_string,var1var2_string,pdf_string,ymin_label_shift,ymax_label_shift,var1_y_breaks) {
#print(pdf_string)

df1$var1 <- df1[,var1]
df2$var1 <- df2[,var1]
df1$var2 <- df1[,var2]
df2$var2 <- df2[,var2]

cph_var1var2_mgh <- summary(coxph(surv ~ var1 + var2,data=df1[df1$cohort=="MGH",]))
cph_var1var2_hug <- summary(coxph(surv ~ var1 + var2,data=df1[df1$cohort=="Hugo"&!is.na(df1$surv),]))
cph_var1var2_ria <- summary(coxph(surv ~ var1 + var2,data=df1[df1$cohort=="Riaz",]))
cph_var1var2_van <- summary(coxph(surv ~ var1 + var2,data=df1[df1$cohort=="VanAllen",]))
cph_var1var2_aus <- summary(coxph(surv ~ var1 + var2 + treatment_type,data=df2[df2$cohort=="Gide",]))
cph_var1var2_liu <- summary(coxph(surv ~ var1 + var2,data=df2[df2$cohort=="Liu",]))

cph_var1var2_var1_mgh_hr <-    cph_var1var2_mgh$conf.int["var1",1]
cph_var1var2_var1_mgh_hr_lo <- cph_var1var2_mgh$conf.int["var1",3]
cph_var1var2_var1_mgh_hr_hi <- cph_var1var2_mgh$conf.int["var1",4]
cph_var1var2_var1_mgh_hr_p <-  cph_var1var2_mgh$coef["var1",5]
cph_var1var2_var1_mgh_lr_p <-  cph_var1var2_mgh$sctest[3]
cph_var1var2_var1_mgh_hr_se <- cph_var1var2_mgh$coef["var1",3]
cph_var1var2_var2_mgh_hr <-    cph_var1var2_mgh$conf.int["var2",1]
cph_var1var2_var2_mgh_hr_lo <- cph_var1var2_mgh$conf.int["var2",3]
cph_var1var2_var2_mgh_hr_hi <- cph_var1var2_mgh$conf.int["var2",4]
cph_var1var2_var2_mgh_hr_p <-  cph_var1var2_mgh$coef["var2",5]
cph_var1var2_var2_mgh_lr_p <-  cph_var1var2_mgh$sctest[3]
cph_var1var2_var2_mgh_hr_se <- cph_var1var2_mgh$coef["var2",3]

cph_var1var2_var1_hug_hr <-    cph_var1var2_hug$conf.int["var1",1]
cph_var1var2_var1_hug_hr_lo <- cph_var1var2_hug$conf.int["var1",3]
cph_var1var2_var1_hug_hr_hi <- cph_var1var2_hug$conf.int["var1",4]
cph_var1var2_var1_hug_hr_p <-  cph_var1var2_hug$coef["var1",5]
cph_var1var2_var1_hug_lr_p <-  cph_var1var2_hug$sctest[3]
cph_var1var2_var1_hug_hr_se <- cph_var1var2_hug$coef["var1",3]
cph_var1var2_var2_hug_hr <-    cph_var1var2_hug$conf.int["var2",1]
cph_var1var2_var2_hug_hr_lo <- cph_var1var2_hug$conf.int["var2",3]
cph_var1var2_var2_hug_hr_hi <- cph_var1var2_hug$conf.int["var2",4]
cph_var1var2_var2_hug_hr_p <-  cph_var1var2_hug$coef["var2",5]
cph_var1var2_var2_hug_lr_p <-  cph_var1var2_hug$sctest[3]
cph_var1var2_var2_hug_hr_se <- cph_var1var2_hug$coef["var2",3]

cph_var1var2_var1_ria_hr <-    cph_var1var2_ria$conf.int["var1",1]
cph_var1var2_var1_ria_hr_lo <- cph_var1var2_ria$conf.int["var1",3]
cph_var1var2_var1_ria_hr_hi <- cph_var1var2_ria$conf.int["var1",4]
cph_var1var2_var1_ria_hr_p <-  cph_var1var2_ria$coef["var1",5]
cph_var1var2_var1_ria_lr_p <-  cph_var1var2_ria$sctest[3]
cph_var1var2_var1_ria_hr_se <- cph_var1var2_ria$coef["var1",3]
cph_var1var2_var2_ria_hr <-    cph_var1var2_ria$conf.int["var2",1]
cph_var1var2_var2_ria_hr_lo <- cph_var1var2_ria$conf.int["var2",3]
cph_var1var2_var2_ria_hr_hi <- cph_var1var2_ria$conf.int["var2",4]
cph_var1var2_var2_ria_hr_p <-  cph_var1var2_ria$coef["var2",5]
cph_var1var2_var2_ria_lr_p <-  cph_var1var2_ria$sctest[3]
cph_var1var2_var2_ria_hr_se <- cph_var1var2_ria$coef["var2",3]

cph_var1var2_var1_van_hr <-    cph_var1var2_van$conf.int["var1",1]
cph_var1var2_var1_van_hr_lo <- cph_var1var2_van$conf.int["var1",3]
cph_var1var2_var1_van_hr_hi <- cph_var1var2_van$conf.int["var1",4]
cph_var1var2_var1_van_hr_p <-  cph_var1var2_van$coef["var1",5]
cph_var1var2_var1_van_lr_p <-  cph_var1var2_van$sctest[3]
cph_var1var2_var1_van_hr_se <- cph_var1var2_van$coef["var1",3]
cph_var1var2_var2_van_hr <-    cph_var1var2_van$conf.int["var2",1]
cph_var1var2_var2_van_hr_lo <- cph_var1var2_van$conf.int["var2",3]
cph_var1var2_var2_van_hr_hi <- cph_var1var2_van$conf.int["var2",4]
cph_var1var2_var2_van_hr_p <-  cph_var1var2_van$coef["var2",5]
cph_var1var2_var2_van_lr_p <-  cph_var1var2_van$sctest[3]
cph_var1var2_var2_van_hr_se <- cph_var1var2_van$coef["var2",3]

cph_var1var2_var1_aus_hr <-    cph_var1var2_aus$conf.int["var1",1]
cph_var1var2_var1_aus_hr_lo <- cph_var1var2_aus$conf.int["var1",3]
cph_var1var2_var1_aus_hr_hi <- cph_var1var2_aus$conf.int["var1",4]
cph_var1var2_var1_aus_hr_p <-  cph_var1var2_aus$coef["var1",5]
cph_var1var2_var1_aus_lr_p <-  cph_var1var2_aus$sctest[3]
cph_var1var2_var1_aus_hr_se <- cph_var1var2_aus$coef["var1",3]
cph_var1var2_var2_aus_hr <-    cph_var1var2_aus$conf.int["var2",1]
cph_var1var2_var2_aus_hr_lo <- cph_var1var2_aus$conf.int["var2",3]
cph_var1var2_var2_aus_hr_hi <- cph_var1var2_aus$conf.int["var2",4]
cph_var1var2_var2_aus_hr_p <-  cph_var1var2_aus$coef["var2",5]
cph_var1var2_var2_aus_lr_p <-  cph_var1var2_aus$sctest[3]
cph_var1var2_var2_aus_hr_se <- cph_var1var2_aus$coef["var2",3]

cph_var1var2_var1_liu_hr <-    cph_var1var2_liu$conf.int["var1",1]
cph_var1var2_var1_liu_hr_lo <- cph_var1var2_liu$conf.int["var1",3]
cph_var1var2_var1_liu_hr_hi <- cph_var1var2_liu$conf.int["var1",4]
cph_var1var2_var1_liu_hr_p <-  cph_var1var2_liu$coef["var1",5]
cph_var1var2_var1_liu_lr_p <-  cph_var1var2_liu$sctest[3]
cph_var1var2_var1_liu_hr_se <- cph_var1var2_liu$coef["var1",3]
cph_var1var2_var2_liu_hr <-    cph_var1var2_liu$conf.int["var2",1]
cph_var1var2_var2_liu_hr_lo <- cph_var1var2_liu$conf.int["var2",3]
cph_var1var2_var2_liu_hr_hi <- cph_var1var2_liu$conf.int["var2",4]
cph_var1var2_var2_liu_hr_p <-  cph_var1var2_liu$coef["var2",5]
cph_var1var2_var2_liu_lr_p <-  cph_var1var2_liu$sctest[3]
cph_var1var2_var2_liu_hr_se <- cph_var1var2_liu$coef["var2",3]

df_cph_var1var2 <- data.frame(rbind(c("MGH",
                     paste0("MGH ",var1_string),cph_var1var2_var1_mgh_hr,cph_var1var2_var1_mgh_hr_lo,cph_var1var2_var1_mgh_hr_hi,cph_var1var2_var1_mgh_hr_p,cph_var1var2_var1_mgh_lr_p,n_mgh),
             c("MGH",paste0("MGH ",var2_string),cph_var1var2_var2_mgh_hr,cph_var1var2_var2_mgh_hr_lo,cph_var1var2_var2_mgh_hr_hi,cph_var1var2_var2_mgh_hr_p,cph_var1var2_var2_mgh_lr_p,n_mgh),
           c("Hugo",paste0("Hugo ",var1_string),cph_var1var2_var1_hug_hr,cph_var1var2_var1_hug_hr_lo,cph_var1var2_var1_hug_hr_hi,cph_var1var2_var1_hug_hr_p,cph_var1var2_var1_hug_lr_p,n_hug),
           c("Hugo",paste0("Hugo ",var2_string),cph_var1var2_var2_hug_hr,cph_var1var2_var2_hug_hr_lo,cph_var1var2_var2_hug_hr_hi,cph_var1var2_var2_hug_hr_p,cph_var1var2_var2_hug_lr_p,n_hug),
           c("Riaz",paste0("Riaz ",var1_string),cph_var1var2_var1_ria_hr,cph_var1var2_var1_ria_hr_lo,cph_var1var2_var1_ria_hr_hi,cph_var1var2_var1_ria_hr_p,cph_var1var2_var1_ria_lr_p,n_ria),
           c("Riaz",paste0("Riaz ",var2_string),cph_var1var2_var2_ria_hr,cph_var1var2_var2_ria_hr_lo,cph_var1var2_var2_ria_hr_hi,cph_var1var2_var2_ria_hr_p,cph_var1var2_var2_ria_lr_p,n_ria),
 c("Van Allen",paste0("Van Allen ",var1_string),cph_var1var2_var1_van_hr,cph_var1var2_var1_van_hr_lo,cph_var1var2_var1_van_hr_hi,cph_var1var2_var1_van_hr_p,cph_var1var2_var1_van_lr_p,n_van),
 c("Van Allen",paste0("Van Allen ",var2_string),cph_var1var2_var2_van_hr,cph_var1var2_var2_van_hr_lo,cph_var1var2_var2_van_hr_hi,cph_var1var2_var2_van_hr_p,cph_var1var2_var2_van_lr_p,n_van),
           c("Gide",paste0("Gide ",var1_string),cph_var1var2_var1_aus_hr,cph_var1var2_var1_aus_hr_lo,cph_var1var2_var1_aus_hr_hi,cph_var1var2_var1_aus_hr_p,cph_var1var2_var1_aus_lr_p,n_aus),
           c("Gide",paste0("Gide ",var2_string),cph_var1var2_var2_aus_hr,cph_var1var2_var2_aus_hr_lo,cph_var1var2_var2_aus_hr_hi,cph_var1var2_var2_aus_hr_p,cph_var1var2_var2_aus_lr_p,n_aus),
             c("Liu",paste0("Liu ",var1_string),cph_var1var2_var1_liu_hr,cph_var1var2_var1_liu_hr_lo,cph_var1var2_var1_liu_hr_hi,cph_var1var2_var1_liu_hr_p,cph_var1var2_var1_liu_lr_p,n_liu),
             c("Liu",paste0("Liu ",var2_string),cph_var1var2_var2_liu_hr,cph_var1var2_var2_liu_hr_lo,cph_var1var2_var2_liu_hr_hi,cph_var1var2_var2_liu_hr_p,cph_var1var2_var2_liu_lr_p,n_liu)),
                    stringsAsFactors=F)


colnames(df_cph_var1var2) <- c("cohort","name","HR","HR_lo","HR_hi","HR_p","LR_p","cohort_n")
df_cph_var1var2$gene <- ifelse(grepl(var1_string,df_cph_var1var2$name),var1_string,var2_string)

if (var1var2_string == "MAP4K1+AGER") {
df_cph_var1var2$gene <- factor(df_cph_var1var2$gene,levels=c("MAP4K1","AGER"))
}

df_cph_var1var2$print_lrt_p <- paste0("p = ",ifelse(df_cph_var1var2$LR_p<0.01,formatC(as.numeric(df_cph_var1var2$LR_p),width=3,format="e",digits=2),formatC(as.numeric(df_cph_var1var2$LR_p),
width=3,digits=2)))
df_cph_var1var2$print_hr_p <- paste0("p = ",ifelse(df_cph_var1var2$HR_p<0.01,formatC(as.numeric(df_cph_var1var2$HR_p),width=3,format="e",digits=2),formatC(as.numeric(df_cph_var1var2$HR_p),
width=3,digits=2)))

df_cph_var1var2$cohort_with_n <- paste0(df_cph_var1var2$cohort,"\n","(n=",df_cph_var1var2$cohort_n,")")
df_cph_var1var2$cohort_with_n <- factor(df_cph_var1var2$cohort_with_n,levels=c("Gide\n(n=73)","Hugo\n(n=25)","Liu\n(n=107)","MGH\n(n=38)","Riaz\n(n=51)","Van Allen\n(n=39)"))


if ("High OS genes" %in% df_cph_var1var2$gene) {
g_var1var2_allmeta <- ggplot(df_cph_var1var2,aes(x=cohort_with_n,y=as.numeric(HR),color=gene)) +
  geom_point(aes(size=as.numeric(cohort_n)),position=position_dodge(width = 0.25)) + geom_pointrange(aes(ymin=as.numeric(HR_lo),ymax=as.numeric(HR_hi)),position=position_dodge(width = 0.25)) +
  theme_classic() + xlab("") + ylab("Hazard Ratio") +
  geom_text(aes(x=cohort_with_n,y=ifelse(gene==var1_string,min(as.numeric(df_cph_var1var2$HR_lo))-ymin_label_shift,max(as.numeric(df_cph_var1var2$HR_hi))+1),label=print_hr_p),
  position=position_dodge(width = 0.25),size=3.5) + geom_hline(yintercept=1,linetype="dashed") + ggtitle(var1var2_string) +
  theme(legend.position="bottom",plot.title=element_text(hjust=0.5,size=16),axis.text.x=element_text(size=16,color="black"),axis.text.y=element_text(size=14,color="black"),
  axis.title.y=element_text(size=16)) + coord_trans(y="log2") + scale_color_manual(values=c("blue","red")) + scale_y_continuous(breaks=var1_y_breaks)
ggsave(pdf_string,g_var1var2_allmeta,height=7,width=8)
} else {
g_var1var2_allmeta <- ggplot(df_cph_var1var2,aes(x=cohort_with_n,y=as.numeric(HR),color=gene)) +
  geom_point(aes(size=as.numeric(cohort_n)),position=position_dodge(width = 0.25)) + geom_pointrange(aes(ymin=as.numeric(HR_lo),ymax=as.numeric(HR_hi)),position=position_dodge(width = 0.25)) +
  theme_classic() + xlab("") + ylab("Hazard Ratio") +
  geom_text(aes(x=cohort_with_n,y=ifelse(gene==var1_string,min(as.numeric(df_cph_var1var2$HR_hi))-ymin_label_shift,
  max(as.numeric(df_cph_var1var2$HR_hi))+ymax_label_shift),label=print_hr_p),
  position=position_dodge(width = 0.25),size=3.5) + geom_hline(yintercept=1,linetype="dashed") + ggtitle(var1var2_string) +
  theme(legend.position="bottom",plot.title=element_text(hjust=0.5,size=16),axis.text.x=element_text(size=16,color="black"),axis.text.y=element_text(size=14,color="black"),
  axis.title.y=element_text(size=16)) + coord_trans(y="log2") + scale_color_manual(values=c("blue","red")) + scale_y_continuous(breaks=var1_y_breaks)
ggsave(pdf_string,g_var1var2_allmeta,height=7,width=8)
}

return(df_cph_var1var2)
}
```

Figure S27

``` {.r}
clin$log10_rna_bcb <- log10(clin$rna_bcb)
skcm_clin$log10_rna_bcb <- log10(skcm_clin$rna_bcb)

n_mgh <- nrow(clin[clin$cohort=="MGH",])
n_hug <- nrow(clin[clin$cohort=="Hugo"&!is.na(clin$surv),])
n_ria <- nrow(clin[clin$cohort=="Riaz",])
n_van <- nrow(clin[clin$cohort=="VanAllen",])
n_aus <- nrow(skcm_clin[skcm_clin$cohort=="Gide",])
n_liu <- nrow(skcm_clin[skcm_clin$cohort=="Liu",])

df_cph_CD274 <- plot_multicohort_survival(clin,skcm_clin,"CD274","CD274","figure4_outputs/S27a.multicohort_CD274_HR_p.pdf",0.5,c(0.5,0.75,1,1.5,2,3,4,5))
knitr::include_graphics("figure4_outputs/S27a.multicohort_CD274_HR_p.pdf")
```

<embed src="figure4_outputs/S27a.multicohort_CD274_HR_p.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
df_cph_CYT <- plot_multicohort_survival(clin,skcm_clin,"CYT","CYT","figure4_outputs/S27b.multicohort_CYT_HR_p.pdf",0.5,c(0.5,0.75,1,1.5,2,3,4,5))
knitr::include_graphics("figure4_outputs/S27b.multicohort_CYT_HR_p.pdf")
```

<embed src="figure4_outputs/S27b.multicohort_CYT_HR_p.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
df_cph_GEP <- plot_multicohort_survival(clin,skcm_clin,"GEP","GEP","figure4_outputs/S27c.multicohort_GEP_HR_p.pdf",0.5,c(0.5,0.75,1,1.5,2,3,4,5))
knitr::include_graphics("figure4_outputs/S27c.multicohort_GEP_HR_p.pdf")
```

<embed src="figure4_outputs/S27c.multicohort_GEP_HR_p.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
df_cph_IMPRES <- plot_multicohort_survival(clin,skcm_clin,"IMPRES","IMPRES","figure4_outputs/S27d.multicohort_IMPRES_HR_p.pdf",0.5,c(0.5,0.75,1,1.5,2,3,4,5))
knitr::include_graphics("figure4_outputs/S27d.multicohort_IMPRES_HR_p.pdf")
```

<embed src="figure4_outputs/S27d.multicohort_IMPRES_HR_p.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
df_cph_MHC_II_zscore <- plot_multicohort_survival(clin,skcm_clin,"MHC_II_zscore","MHC II z-score","figure4_outputs/S27e.multicohort_MHC_II_zscore_HR_p.pdf",0.5,c(0.5,0.75,1,1.5,2,3,4,5))
knitr::include_graphics("figure4_outputs/S27e.multicohort_MHC_II_zscore_HR_p.pdf")
```

<embed src="figure4_outputs/S27e.multicohort_MHC_II_zscore_HR_p.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
clin$TIDE <- clin$TIDE_Score
skcm_clin$TIDE <- skcm_clin$TIDE_Score
df_cph_TIDE <- plot_multicohort_survival(clin,skcm_clin,"TIDE","TIDE","figure4_outputs/S27f.multicohort_TIDE_HR_p.pdf",0.5,c(0.5,0.75,1,1.5,2,3,4,5))
knitr::include_graphics("figure4_outputs/S27f.multicohort_TIDE_HR_p.pdf")
```

<embed src="figure4_outputs/S27f.multicohort_TIDE_HR_p.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
df_cph_TCB <- plot_multicohort_survival(clin,skcm_clin,"log10_rna_tcb","TCB","figure4_outputs/S27g.multicohort_TCB_HR_p.pdf",0.5,c(0.25,0.5,0.75,1,1.5,2,3,4,5))
knitr::include_graphics("figure4_outputs/S27g.multicohort_TCB_HR_p.pdf")
```

<embed src="figure4_outputs/S27g.multicohort_TCB_HR_p.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
df_cph_BCB <- plot_multicohort_survival(clin,skcm_clin,"log10_rna_bcb","BCB","figure4_outputs/S27h.multicohort_BCB_HR_p.pdf",0.5,c(0.25,0.5,0.75,1,1.5,2,3,4,5))
knitr::include_graphics("figure4_outputs/S27h.multicohort_BCB_HR_p.pdf")
```

<embed src="figure4_outputs/S27h.multicohort_BCB_HR_p.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
df_cph_MAP4K1TBX3 <- plot_multicohort_survival_double(clin,skcm_clin,"MAP4K1","MAP4K1","TBX3","TBX3","MAP4K1+TBX3","figure4_outputs/S27i.multicohort_MAP4K1TBX3.pdf",0.55,0.5,c(0.5,0.75,1,1.5,2,3,4,5))
knitr::include_graphics("figure4_outputs/S27i.multicohort_MAP4K1TBX3.pdf")
```

<embed src="figure4_outputs/S27i.multicohort_MAP4K1TBX3.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
df_cph_MAP4K1AGER <- plot_multicohort_survival_double(clin,skcm_clin,"MAP4K1","MAP4K1","AGER","AGER","MAP4K1+AGER","figure4_outputs/S27j.multicohort_MAP4K1AGER.pdf",0.55,0.5,c(0.5,0.75,1,1.5,2,3,4,5))
knitr::include_graphics("figure4_outputs/S27j.multicohort_MAP4K1AGER.pdf")
```

<embed src="figure4_outputs/S27j.multicohort_MAP4K1AGER.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
df_cph_meanhighOSmeanlowOS <- plot_multicohort_survival_double(clin,skcm_clin,"mean_highos_genes_z","High OS genes",
                                   "mean_lowos_genes_z","Low OS genes","Mean of high OS and low OS genes","figure4_outputs/S27k.multicohort_meanhighOSmeanlowOS.pdf",0.01,1,c(0.25,0.5,0.75,1,1.5,2,3,4,5))
knitr::include_graphics("figure4_outputs/S27k.multicohort_meanhighOSmeanlowOS.pdf")
```

<embed src="figure4_outputs/S27k.multicohort_meanhighOSmeanlowOS.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
multicohort_table <- data.frame(df_cph_CD274$cohort,stringsAsFactors=F)
multicohort_table[,c("CD274_HR","CD274_HR_lo","CD274_HR_hi","CD274_HR_p","CD274_LRT_p")] <- df_cph_CD274[,c("HR","HR_lo","HR_hi","HR_p","LR_p")]
multicohort_table[,c("GEP_HR","GEP_HR_lo","GEP_HR_hi","GEP_HR_p","GEP_LRT_p")] <- df_cph_GEP[,c("HR","HR_lo","HR_hi","HR_p","LR_p")]
multicohort_table[,c("CYT_HR","CYT_HR_lo","CYT_HR_hi","CYT_HR_p","CYT_LRT_p")] <- df_cph_CYT[,c("HR","HR_lo","HR_hi","HR_p","LR_p")]
multicohort_table[,c("IMPRES_HR","IMPRES_HR_lo","IMPRES_HR_hi","IMPRES_HR_p","IMPRES_LRT_p")] <- df_cph_IMPRES[,c("HR","HR_lo","HR_hi","HR_p","LR_p")]
multicohort_table[,c("TIDE_HR","TIDE_HR_lo","TIDE_HR_hi","TIDE_HR_p","TIDE_LRT_p")] <- df_cph_TIDE[,c("HR","HR_lo","HR_hi","HR_p","LR_p")]
multicohort_table[,c("MHC_II_HR","MHC_II_HR_lo","MHC_II_HR_hi","MHC_II_HR_p","MHC_II_LRT_p")] <- df_cph_MHC_II_zscore[,c("HR","HR_lo","HR_hi","HR_p","LR_p")]
multicohort_table[,c("TCB_HR","TCB_HR_lo","TCB_HR_hi","TCB_HR_p","TCB_LRT_p")] <- df_cph_TCB[,c("HR","HR_lo","HR_hi","HR_p","LR_p")]
multicohort_table[,c("BCB_HR","BCB_HR_lo","BCB_HR_hi","BCB_HR_p","BCB_LRT_p")] <- df_cph_BCB[,c("HR","HR_lo","HR_hi","HR_p","LR_p")]
multicohort_table[,c("MAP4K1TBX3_MAP4K1_HR","MAP4K1TBX3_MAP4K1_HR_lo","MAP4K1TBX3_MAP4K1_HR_hi","MAP4K1TBX3_MAP4K1_HR_p","MAP4K1TBX3_MAP4K1_LRT_p")] <-
                                                                                                        df_cph_MAP4K1TBX3[df_cph_MAP4K1TBX3$gene=="MAP4K1",c("HR","HR_lo","HR_hi","HR_p","LR_p")]
multicohort_table[,c("MAP4K1TBX3_TBX3_HR","MAP4K1TBX3_TBX3_HR_lo","MAP4K1TBX3_TBX3_HR_hi","MAP4K1TBX3_TBX3_HR_p","MAP4K1TBX3_TBX3_LRT_p")] <-
                                                                                                        df_cph_MAP4K1TBX3[df_cph_MAP4K1TBX3$gene=="TBX3",c("HR","HR_lo","HR_hi","HR_p","LR_p")]

multicohort_table[,c("MAP4K1AGER_MAP4K1_HR","MAP4K1AGER_MAP4K1_HR_lo","MAP4K1AGER_MAP4K1_HR_hi","MAP4K1AGER_MAP4K1_HR_p","MAP4K1AGER_MAP4K1_LRT_p")] <-
                                                                                                        df_cph_MAP4K1AGER[df_cph_MAP4K1AGER$gene=="MAP4K1",c("HR","HR_lo","HR_hi","HR_p","LR_p")]
multicohort_table[,c("MAP4K1AGER_AGER_HR","MAP4K1AGER_AGER_HR_lo","MAP4K1AGER_AGER_HR_hi","MAP4K1AGER_AGER_HR_p","MAP4K1AGER_AGER_LRT_p")] <-
                                                                                                        df_cph_MAP4K1AGER[df_cph_MAP4K1AGER$gene=="AGER",c("HR","HR_lo","HR_hi","HR_p","LR_p")]

multicohort_table[,c("OSmetagene_LongOS_HR","OSmetagene_LongOS_HR_lo","OSmetagene_LongOS_HR_hi","OSmetagene_LongOS_HR_p","OSmetagene_LongOS_LRT_p")] <-
                                                                        df_cph_meanhighOSmeanlowOS[df_cph_meanhighOSmeanlowOS$gene=="High OS genes",c("HR","HR_lo","HR_hi","HR_p","LR_p")]
multicohort_table[,c("OSmetagene_ShortOS_HR","OSmetagene_ShortOS_HR_lo","OSmetagene_ShortOS_HR_hi","OSmetagene_ShortOS_HR_p","OSmetagene_ShortOS_LRT_p")] <-
                                                                        df_cph_meanhighOSmeanlowOS[df_cph_meanhighOSmeanlowOS$gene=="Low OS genes",c("HR","HR_lo","HR_hi","HR_p","LR_p")]

write.table(multicohort_table,"figure4_outputs/CPBVal_withincohort_performance.txt",quote=F,sep="\t")
```

Functions

``` {.r}
plot_bytreatment_survival <- function(df1,var1,var1_string,pdf_string,ymax_label_shift,var1_y_breaks) {
#plot_bytreatment_survival(comb_clin,"CD274","CD274","S28a.bytreatment_CD274_HR_p.pdf",0.5,c(0.5,0.75,1,1.5,2,3,4,5))

#print(pdf_string)

df1$var1 <- df1[,var1]

cph_var1_ctla4 <- summary(coxph(surv ~ var1,data=df1[df1$therapy_simplified=="CTLA4",]))
cph_var1_pd1 <- summary(coxph(surv ~ var1,data=df1[df1$therapy_simplified=="PD1",]))
cph_var1_pd1ctla4 <- summary(coxph(surv ~ var1,data=df1[df1$therapy_simplified=="CTLA4+PD1",]))
cph_var1_pd1nopriorctla4 <- summary(coxph(surv ~ var1,data=df1[!is.na(df1$prior_CTLA4)&df1$therapy_simplified=="PD1"&df1$prior_CTLA4==FALSE,]))
cph_var1_pd1priorctla4 <- summary(coxph(surv ~ var1,data=df1[!is.na(df1$prior_CTLA4)&df1$therapy_simplified=="PD1"&df1$prior_CTLA4==TRUE,]))

cph_var1_ctla4_hr <-    cph_var1_ctla4$conf.int["var1",1]
cph_var1_ctla4_hr_lo <- cph_var1_ctla4$conf.int["var1",3]
cph_var1_ctla4_hr_hi <- cph_var1_ctla4$conf.int["var1",4]
cph_var1_ctla4_hr_p <-  cph_var1_ctla4$coef["var1",5]
cph_var1_ctla4_lr_p <-  cph_var1_ctla4$sctest[3]
cph_var1_ctla4_hr_se <- cph_var1_ctla4$coef["var1",3]

cph_var1_pd1_hr <-    cph_var1_pd1$conf.int["var1",1]
cph_var1_pd1_hr_lo <- cph_var1_pd1$conf.int["var1",3]
cph_var1_pd1_hr_hi <- cph_var1_pd1$conf.int["var1",4]
cph_var1_pd1_hr_p <-  cph_var1_pd1$coef["var1",5]
cph_var1_pd1_lr_p <-  cph_var1_pd1$sctest[3]
cph_var1_pd1_hr_se <- cph_var1_pd1$coef["var1",3]

cph_var1_pd1ctla4_hr <-    cph_var1_pd1ctla4$conf.int["var1",1]
cph_var1_pd1ctla4_hr_lo <- cph_var1_pd1ctla4$conf.int["var1",3]
cph_var1_pd1ctla4_hr_hi <- cph_var1_pd1ctla4$conf.int["var1",4]
cph_var1_pd1ctla4_hr_p <-  cph_var1_pd1ctla4$coef["var1",5]
cph_var1_pd1ctla4_lr_p <-  cph_var1_pd1ctla4$sctest[3]
cph_var1_pd1ctla4_hr_se <- cph_var1_pd1ctla4$coef["var1",3]

cph_var1_pd1nopriorctla4_hr <-    cph_var1_pd1nopriorctla4$conf.int["var1",1]
cph_var1_pd1nopriorctla4_hr_lo <- cph_var1_pd1nopriorctla4$conf.int["var1",3]
cph_var1_pd1nopriorctla4_hr_hi <- cph_var1_pd1nopriorctla4$conf.int["var1",4]
cph_var1_pd1nopriorctla4_hr_p <-  cph_var1_pd1nopriorctla4$coef["var1",5]
cph_var1_pd1nopriorctla4_lr_p <-  cph_var1_pd1nopriorctla4$sctest[3]
cph_var1_pd1nopriorctla4_hr_se <- cph_var1_pd1nopriorctla4$coef["var1",3]

cph_var1_pd1priorctla4_hr <-    cph_var1_pd1priorctla4$conf.int["var1",1]
cph_var1_pd1priorctla4_hr_lo <- cph_var1_pd1priorctla4$conf.int["var1",3]
cph_var1_pd1priorctla4_hr_hi <- cph_var1_pd1priorctla4$conf.int["var1",4]
cph_var1_pd1priorctla4_hr_p <-  cph_var1_pd1priorctla4$coef["var1",5]
cph_var1_pd1priorctla4_lr_p <-  cph_var1_pd1priorctla4$sctest[3]
cph_var1_pd1priorctla4_hr_se <- cph_var1_pd1priorctla4$coef["var1",3]

n_ctla4 <- nrow(df1[df1$therapy_simplified=="CTLA4",])
n_pd1 <- nrow(df1[df1$therapy_simplified=="PD1",])
n_pd1ctla4 <- nrow(df1[df1$therapy_simplified=="CTLA4+PD1",])
n_pd1nopriorctla4 <- nrow(df1[!is.na(df1$prior_CTLA4)&df1$therapy_simplified=="PD1"&df1$prior_CTLA4==FALSE,])
n_pd1priorctla4 <- nrow(df1[!is.na(df1$prior_CTLA4)&df1$therapy_simplified=="PD1"&df1$prior_CTLA4==TRUE,])

df_cph_var1 <- data.frame(rbind(c("CTLA4",paste0("CTLA4 ",var1_string),cph_var1_ctla4_hr,cph_var1_ctla4_hr_lo,cph_var1_ctla4_hr_hi,cph_var1_ctla4_hr_p,cph_var1_ctla4_lr_p,n_ctla4),
               c("PD1",paste0("PD1 ",var1_string),cph_var1_pd1_hr,cph_var1_pd1_hr_lo,cph_var1_pd1_hr_hi,cph_var1_pd1_hr_p,cph_var1_pd1_lr_p,n_pd1),
               c("PD1+CTLA4",paste0("PD1+CTLA4 ",var1_string),cph_var1_pd1ctla4_hr,cph_var1_pd1ctla4_hr_lo,cph_var1_pd1ctla4_hr_hi,cph_var1_pd1ctla4_hr_p,cph_var1_pd1ctla4_lr_p,n_pd1ctla4),
               c("PD1\nno prior CTLA4",paste0("PD1, no prior CTLA4 ",var1_string),cph_var1_pd1nopriorctla4_hr,cph_var1_pd1nopriorctla4_hr_lo,cph_var1_pd1nopriorctla4_hr_hi,
                       cph_var1_pd1nopriorctla4_hr_p,cph_var1_pd1nopriorctla4_lr_p,n_pd1nopriorctla4),
               c("PD1\nprior CTLA4",paste0("PD1, priorCTLA4 ",var1_string),cph_var1_pd1priorctla4_hr,cph_var1_pd1priorctla4_hr_lo,cph_var1_pd1priorctla4_hr_hi,
                       cph_var1_pd1priorctla4_hr_p,cph_var1_pd1priorctla4_lr_p,n_pd1priorctla4)),
               stringsAsFactors=F)

colnames(df_cph_var1) <- c("cohort","name","HR","HR_lo","HR_hi","HR_p","LR_p","cohort_n")

df_cph_var1$print_lrt_p <- paste0("p = ",ifelse(df_cph_var1$LR_p<0.01,formatC(as.numeric(df_cph_var1$LR_p),width=3,format="e",digits=2),formatC(as.numeric(df_cph_var1$LR_p),width=3,digits=2)))
df_cph_var1$print_hr_p <- paste0("p = ",ifelse(df_cph_var1$HR_p<0.01,formatC(as.numeric(df_cph_var1$HR_p),width=3,format="e",digits=2),formatC(as.numeric(df_cph_var1$HR_p),width=3,digits=2)))

df_cph_var1$cohort_with_n <- paste0(df_cph_var1$cohort,"\n","(n=",df_cph_var1$cohort_n,")")
df_cph_var1$cohort_with_n <- factor(df_cph_var1$cohort_with_n,levels=c("CTLA4\n(n=57)","PD1\n(n=242)","PD1+CTLA4\n(n=34)","PD1\nno prior CTLA4\n(n=129)",
                          "PD1\nprior CTLA4\n(n=72)"))

g_var1_allmeta_hr_p <- ggplot(df_cph_var1,aes(x=cohort_with_n,y=as.numeric(HR))) +
  geom_point(aes(size=as.numeric(cohort_n))) + geom_pointrange(aes(ymin=as.numeric(HR_lo),ymax=as.numeric(HR_hi))) + theme_classic() +
  xlab("") + ylab("Hazard Ratio") + geom_text(aes(x=cohort_with_n,y=max(as.numeric(df_cph_var1$HR_hi))+ymax_label_shift,label=print_hr_p),size=3.5) +
  geom_hline(yintercept=1,linetype="dashed") + ggtitle(var1_string) +
  theme(legend.position="bottom",plot.title=element_text(hjust=0.5,size=16),axis.text.x=element_text(size=16,color="black"),axis.text.y=element_text(size=14,color="black"),
  axis.title.y=element_text(size=16)) + coord_trans(y="log2") + scale_color_manual(values=c("red","blue")) + scale_y_continuous(breaks=var1_y_breaks) +
  geom_vline(xintercept=3.5,linetype="dashed")
ggsave(pdf_string,g_var1_allmeta_hr_p,height=7,width=9)

return(df_cph_var1)
}

plot_bytreatment_survival_double <- function(df1,var1,var1_string,var2,var2_string,var1var2_string,pdf_string,ymin_label_shift,ymax_label_shift,var1_y_breaks) {
#print(pdf_string)

df1$var1 <- df1[,var1]
df1$var2 <- df1[,var2]

cph_var1var2_ctla4 <- summary(coxph(surv ~ var1 + var2,data=df1[df1$therapy_simplified=="CTLA4",]))
cph_var1var2_pd1 <- summary(coxph(surv ~ var1 + var2,data=df1[df1$therapy_simplified=="PD1",]))
cph_var1var2_pd1ctla4 <- summary(coxph(surv ~ var1 + var2,data=df1[df1$therapy_simplified=="CTLA4+PD1",]))
cph_var1var2_pd1nopriorctla4 <- summary(coxph(surv ~ var1 + var2,data=df1[!is.na(df1$prior_CTLA4)&df1$therapy_simplified=="PD1"&df1$prior_CTLA4==FALSE,]))
cph_var1var2_pd1priorctla4 <- summary(coxph(surv ~ var1 + var2,data=df1[!is.na(df1$prior_CTLA4)&df1$therapy_simplified=="PD1"&df1$prior_CTLA4==TRUE,]))

cph_var1var2_var1_ctla4_hr <-    cph_var1var2_ctla4$conf.int["var1",1]
cph_var1var2_var1_ctla4_hr_lo <- cph_var1var2_ctla4$conf.int["var1",3]
cph_var1var2_var1_ctla4_hr_hi <- cph_var1var2_ctla4$conf.int["var1",4]
cph_var1var2_var1_ctla4_hr_p <-  cph_var1var2_ctla4$coef["var1",5]
cph_var1var2_var1_ctla4_lr_p <-  cph_var1var2_ctla4$sctest[3]
cph_var1var2_var1_ctla4_hr_se <- cph_var1var2_ctla4$coef["var1",3]
cph_var1var2_var2_ctla4_hr <-    cph_var1var2_ctla4$conf.int["var2",1]
cph_var1var2_var2_ctla4_hr_lo <- cph_var1var2_ctla4$conf.int["var2",3]
cph_var1var2_var2_ctla4_hr_hi <- cph_var1var2_ctla4$conf.int["var2",4]
cph_var1var2_var2_ctla4_hr_p <-  cph_var1var2_ctla4$coef["var2",5]
cph_var1var2_var2_ctla4_lr_p <-  cph_var1var2_ctla4$sctest[3]
cph_var1var2_var2_ctla4_hr_se <- cph_var1var2_ctla4$coef["var2",3]

cph_var1var2_var1_pd1_hr <-    cph_var1var2_pd1$conf.int["var1",1]
cph_var1var2_var1_pd1_hr_lo <- cph_var1var2_pd1$conf.int["var1",3]
cph_var1var2_var1_pd1_hr_hi <- cph_var1var2_pd1$conf.int["var1",4]
cph_var1var2_var1_pd1_hr_p <-  cph_var1var2_pd1$coef["var1",5]
cph_var1var2_var1_pd1_lr_p <-  cph_var1var2_pd1$sctest[3]
cph_var1var2_var1_pd1_hr_se <- cph_var1var2_pd1$coef["var1",3]
cph_var1var2_var2_pd1_hr <-    cph_var1var2_pd1$conf.int["var2",1]
cph_var1var2_var2_pd1_hr_lo <- cph_var1var2_pd1$conf.int["var2",3]
cph_var1var2_var2_pd1_hr_hi <- cph_var1var2_pd1$conf.int["var2",4]
cph_var1var2_var2_pd1_hr_p <-  cph_var1var2_pd1$coef["var2",5]
cph_var1var2_var2_pd1_lr_p <-  cph_var1var2_pd1$sctest[3]
cph_var1var2_var2_pd1_hr_se <- cph_var1var2_pd1$coef["var2",3]

cph_var1var2_var1_pd1ctla4_hr <-    cph_var1var2_pd1ctla4$conf.int["var1",1]
cph_var1var2_var1_pd1ctla4_hr_lo <- cph_var1var2_pd1ctla4$conf.int["var1",3]
cph_var1var2_var1_pd1ctla4_hr_hi <- cph_var1var2_pd1ctla4$conf.int["var1",4]
cph_var1var2_var1_pd1ctla4_hr_p <-  cph_var1var2_pd1ctla4$coef["var1",5]
cph_var1var2_var1_pd1ctla4_lr_p <-  cph_var1var2_pd1ctla4$sctest[3]
cph_var1var2_var1_pd1ctla4_hr_se <- cph_var1var2_pd1ctla4$coef["var1",3]
cph_var1var2_var2_pd1ctla4_hr <-    cph_var1var2_pd1ctla4$conf.int["var2",1]
cph_var1var2_var2_pd1ctla4_hr_lo <- cph_var1var2_pd1ctla4$conf.int["var2",3]
cph_var1var2_var2_pd1ctla4_hr_hi <- cph_var1var2_pd1ctla4$conf.int["var2",4]
cph_var1var2_var2_pd1ctla4_hr_p <-  cph_var1var2_pd1ctla4$coef["var2",5]
cph_var1var2_var2_pd1ctla4_lr_p <-  cph_var1var2_pd1ctla4$sctest[3]
cph_var1var2_var2_pd1ctla4_hr_se <- cph_var1var2_pd1ctla4$coef["var2",3]
cph_var1var2_var1_pd1nopriorctla4_hr <-    cph_var1var2_pd1nopriorctla4$conf.int["var1",1]
cph_var1var2_var1_pd1nopriorctla4_hr_lo <- cph_var1var2_pd1nopriorctla4$conf.int["var1",3]
cph_var1var2_var1_pd1nopriorctla4_hr_hi <- cph_var1var2_pd1nopriorctla4$conf.int["var1",4]
cph_var1var2_var1_pd1nopriorctla4_hr_p <-  cph_var1var2_pd1nopriorctla4$coef["var1",5]
cph_var1var2_var1_pd1nopriorctla4_lr_p <-  cph_var1var2_pd1nopriorctla4$sctest[3]
cph_var1var2_var1_pd1nopriorctla4_hr_se <- cph_var1var2_pd1nopriorctla4$coef["var1",3]
cph_var1var2_var2_pd1nopriorctla4_hr <-    cph_var1var2_pd1nopriorctla4$conf.int["var2",1]
cph_var1var2_var2_pd1nopriorctla4_hr_lo <- cph_var1var2_pd1nopriorctla4$conf.int["var2",3]
cph_var1var2_var2_pd1nopriorctla4_hr_hi <- cph_var1var2_pd1nopriorctla4$conf.int["var2",4]
cph_var1var2_var2_pd1nopriorctla4_hr_p <-  cph_var1var2_pd1nopriorctla4$coef["var2",5]
cph_var1var2_var2_pd1nopriorctla4_lr_p <-  cph_var1var2_pd1nopriorctla4$sctest[3]
cph_var1var2_var2_pd1nopriorctla4_hr_se <- cph_var1var2_pd1nopriorctla4$coef["var2",3]

cph_var1var2_var1_pd1priorctla4_hr <-    cph_var1var2_pd1priorctla4$conf.int["var1",1]
cph_var1var2_var1_pd1priorctla4_hr_lo <- cph_var1var2_pd1priorctla4$conf.int["var1",3]
cph_var1var2_var1_pd1priorctla4_hr_hi <- cph_var1var2_pd1priorctla4$conf.int["var1",4]
cph_var1var2_var1_pd1priorctla4_hr_p <-  cph_var1var2_pd1priorctla4$coef["var1",5]
cph_var1var2_var1_pd1priorctla4_lr_p <-  cph_var1var2_pd1priorctla4$sctest[3]
cph_var1var2_var1_pd1priorctla4_hr_se <- cph_var1var2_pd1priorctla4$coef["var1",3]
cph_var1var2_var2_pd1priorctla4_hr <-    cph_var1var2_pd1priorctla4$conf.int["var2",1]
cph_var1var2_var2_pd1priorctla4_hr_lo <- cph_var1var2_pd1priorctla4$conf.int["var2",3]
cph_var1var2_var2_pd1priorctla4_hr_hi <- cph_var1var2_pd1priorctla4$conf.int["var2",4]
cph_var1var2_var2_pd1priorctla4_hr_p <-  cph_var1var2_pd1priorctla4$coef["var2",5]
cph_var1var2_var2_pd1priorctla4_lr_p <-  cph_var1var2_pd1priorctla4$sctest[3]
cph_var1var2_var2_pd1priorctla4_hr_se <- cph_var1var2_pd1priorctla4$coef["var2",3]

df_cph_var1var2 <- data.frame(rbind(c("CTLA4",
              paste0("CTLA4 ",var1_string),cph_var1var2_var1_ctla4_hr,cph_var1var2_var1_ctla4_hr_lo,cph_var1var2_var1_ctla4_hr_hi,cph_var1var2_var1_ctla4_hr_p,cph_var1var2_var1_ctla4_lr_p,n_ctla4),
       c("CTLA4",paste0("CTLA4 ",var2_string),cph_var1var2_var2_ctla4_hr,cph_var1var2_var2_ctla4_hr_lo,cph_var1var2_var2_ctla4_hr_hi,cph_var1var2_var2_ctla4_hr_p,cph_var1var2_var2_ctla4_lr_p,n_ctla4),
           c("PD1",paste0("PD1 ",var1_string),cph_var1var2_var1_pd1_hr,cph_var1var2_var1_pd1_hr_lo,cph_var1var2_var1_pd1_hr_hi,cph_var1var2_var1_pd1_hr_p,cph_var1var2_var1_pd1_lr_p,n_pd1),
           c("PD1",paste0("PD1 ",var2_string),cph_var1var2_var2_pd1_hr,cph_var1var2_var2_pd1_hr_lo,cph_var1var2_var2_pd1_hr_hi,cph_var1var2_var2_pd1_hr_p,cph_var1var2_var2_pd1_lr_p,n_pd1),
           c("PD1+CTLA4",paste0("PD1+CTLA4 ",var1_string),cph_var1var2_var1_pd1ctla4_hr,cph_var1var2_var1_pd1ctla4_hr_lo,cph_var1var2_var1_pd1ctla4_hr_hi,cph_var1var2_var1_pd1ctla4_hr_p,
           cph_var1var2_var1_pd1ctla4_lr_p,n_pd1ctla4),
           c("PD1+CTLA4",paste0("PD1+CTLA4 ",var2_string),cph_var1var2_var2_pd1ctla4_hr,cph_var1var2_var2_pd1ctla4_hr_lo,cph_var1var2_var2_pd1ctla4_hr_hi,cph_var1var2_var2_pd1ctla4_hr_p,
           cph_var1var2_var2_pd1ctla4_lr_p,n_pd1ctla4),
           c("PD1\nno prior CTLA4",paste0("PD1\nno prior CTLA4 ",var1_string),cph_var1var2_var1_pd1nopriorctla4_hr,cph_var1var2_var1_pd1nopriorctla4_hr_lo,cph_var1var2_var1_pd1nopriorctla4_hr_hi,
           cph_var1var2_var1_pd1nopriorctla4_hr_p,cph_var1var2_var1_pd1nopriorctla4_lr_p,n_pd1nopriorctla4),
           c("PD1\nno prior CTLA4",paste0("PD1\nno prior CTLA4 ",var2_string),cph_var1var2_var2_pd1nopriorctla4_hr,cph_var1var2_var2_pd1nopriorctla4_hr_lo,cph_var1var2_var2_pd1nopriorctla4_hr_hi,
           cph_var1var2_var2_pd1nopriorctla4_hr_p,cph_var1var2_var2_pd1nopriorctla4_lr_p,n_pd1nopriorctla4),
           c("PD1\nprior CTLA4",paste0("PD1\nprior CTLA4 ",var1_string),cph_var1var2_var1_pd1priorctla4_hr,cph_var1var2_var1_pd1priorctla4_hr_lo,cph_var1var2_var1_pd1priorctla4_hr_hi,
           cph_var1var2_var1_pd1priorctla4_hr_p,cph_var1var2_var1_pd1priorctla4_lr_p,n_pd1priorctla4),
           c("PD1\nprior CTLA4",paste0("PD1\nprior CTLA4 ",var2_string),cph_var1var2_var2_pd1priorctla4_hr,cph_var1var2_var2_pd1priorctla4_hr_lo,cph_var1var2_var2_pd1priorctla4_hr_hi,
           cph_var1var2_var2_pd1priorctla4_hr_p,cph_var1var2_var2_pd1priorctla4_lr_p,n_pd1priorctla4)),
           stringsAsFactors=F)

colnames(df_cph_var1var2) <- c("cohort","name","HR","HR_lo","HR_hi","HR_p","LR_p","cohort_n")
df_cph_var1var2$gene <- ifelse(grepl(var1_string,df_cph_var1var2$name),var1_string,var2_string)

if (var1var2_string == "MAP4K1+AGER") {
df_cph_var1var2$gene <- factor(df_cph_var1var2$gene,levels=c("MAP4K1","AGER"))
}

df_cph_var1var2$print_lrt_p <- paste0("p = ",ifelse(df_cph_var1var2$LR_p<0.01,formatC(as.numeric(df_cph_var1var2$LR_p),width=3,format="e",digits=2),formatC(as.numeric(df_cph_var1var2$LR_p),
width=3,digits=2)))
df_cph_var1var2$print_hr_p <- paste0("p = ",ifelse(df_cph_var1var2$HR_p<0.01,formatC(as.numeric(df_cph_var1var2$HR_p),width=3,format="e",digits=2),formatC(as.numeric(df_cph_var1var2$HR_p),
width=3,digits=2)))

df_cph_var1var2$cohort_with_n <- paste0(df_cph_var1var2$cohort,"\n","(n=",df_cph_var1var2$cohort_n,")")
df_cph_var1var2$cohort_with_n <- factor(df_cph_var1var2$cohort_with_n,levels=c("CTLA4\n(n=57)","PD1\n(n=242)","PD1+CTLA4\n(n=34)","PD1\nno prior CTLA4\n(n=129)",
                          "PD1\nprior CTLA4\n(n=72)"))

if ("High OS genes" %in% df_cph_var1var2$gene) {
g_var1var2_allmeta <- ggplot(df_cph_var1var2,aes(x=cohort_with_n,y=as.numeric(HR),color=gene)) +
  geom_point(aes(size=as.numeric(cohort_n)),position=position_dodge(width = 0.25)) + geom_pointrange(aes(ymin=as.numeric(HR_lo),ymax=as.numeric(HR_hi)),position=position_dodge(width = 0.25)) +
  theme_classic() + xlab("") + ylab("Hazard Ratio") +
  geom_text(aes(x=cohort_with_n,y=ifelse(gene==var1_string,min(as.numeric(df_cph_var1var2$HR_lo))-ymin_label_shift,max(as.numeric(df_cph_var1var2$HR_hi))+1),label=print_hr_p),
  position=position_dodge(width = 0.25),size=3.5) + geom_hline(yintercept=1,linetype="dashed") + ggtitle(var1var2_string) +
  theme(legend.position="bottom",plot.title=element_text(hjust=0.5,size=16),axis.text.x=element_text(size=16,color="black"),axis.text.y=element_text(size=14,color="black"),
  axis.title.y=element_text(size=16)) + coord_trans(y="log2") + scale_color_manual(values=c("blue","red")) + scale_y_continuous(breaks=var1_y_breaks) +
  geom_vline(xintercept=3.5,linetype="dashed")
ggsave(pdf_string,g_var1var2_allmeta,height=7,width=9)
} else {
g_var1var2_allmeta <- ggplot(df_cph_var1var2,aes(x=cohort_with_n,y=as.numeric(HR),color=gene)) +
  geom_point(aes(size=as.numeric(cohort_n)),position=position_dodge(width = 0.25)) + geom_pointrange(aes(ymin=as.numeric(HR_lo),ymax=as.numeric(HR_hi)),position=position_dodge(width = 0.25)) +
  theme_classic() + xlab("") + ylab("Hazard Ratio") +
  geom_text(aes(x=cohort_with_n,y=ifelse(gene==var1_string,min(as.numeric(df_cph_var1var2$HR_hi))-ymin_label_shift,
  max(as.numeric(df_cph_var1var2$HR_hi))+ymax_label_shift),label=print_hr_p),
  position=position_dodge(width = 0.25),size=3.5) + geom_hline(yintercept=1,linetype="dashed") + ggtitle(var1var2_string) +
  theme(legend.position="bottom",plot.title=element_text(hjust=0.5,size=16),axis.text.x=element_text(size=16,color="black"),axis.text.y=element_text(size=14,color="black"),
  axis.title.y=element_text(size=16)) + coord_trans(y="log2") + scale_color_manual(values=c("blue","red")) + scale_y_continuous(breaks=var1_y_breaks) +
  geom_vline(xintercept=3.5,linetype="dashed")
ggsave(pdf_string,g_var1var2_allmeta,height=7,width=9)
}


return(df_cph_var1var2)
}
```

Figure S28

``` {.r}
comb_clin <- readRDS("figure4_inputs/comb_clin.rds")
comb_clin$surv <- Surv(comb_clin$overall_survival,comb_clin$DOD=="DOD")
comb_clin <- comb_clin[!is.na(comb_clin$surv),]

n_ctla4 <- nrow(comb_clin[comb_clin$therapy_simplified=="CTLA4",])
n_pd1 <- nrow(comb_clin[comb_clin$therapy_simplified=="PD1",])
n_pd1ctla4 <- nrow(comb_clin[comb_clin$therapy_simplified=="CTLA4+PD1",])
n_pd1nopriorctla4 <- nrow(comb_clin[!is.na(comb_clin$prior_CTLA4)&comb_clin$therapy_simplified=="PD1"&comb_clin$prior_CTLA4==FALSE,])
n_pd1priorctla4 <- nrow(comb_clin[!is.na(comb_clin$prior_CTLA4)&comb_clin$therapy_simplified=="PD1"&comb_clin$prior_CTLA4==TRUE,])

df_comb_CD274 <- plot_bytreatment_survival(comb_clin,"CD274","CD274","figure4_outputs/S28a.bytreatment_CD274_HR_p.pdf",0.5,c(0.5,0.75,1,1.5,2,3,4,5))
knitr::include_graphics("figure4_outputs/S28a.bytreatment_CD274_HR_p.pdf")
```

<embed src="figure4_outputs/S28a.bytreatment_CD274_HR_p.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
df_comb_CYT <- plot_bytreatment_survival(comb_clin,"CYT","CYT","figure4_outputs/S28b.bytreatment_CYT_HR_p.pdf",0.5,c(0.5,0.75,1,1.5,2,3,4,5))
knitr::include_graphics("figure4_outputs/S28b.bytreatment_CYT_HR_p.pdf")
```

<embed src="figure4_outputs/S28b.bytreatment_CYT_HR_p.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
df_comb_GEP <- plot_bytreatment_survival(comb_clin,"GEP","GEP","figure4_outputs/S28c.bytreatment_GEP_HR_p.pdf",0.5,c(0.5,0.75,1,1.5,2,3,4,5))
knitr::include_graphics("figure4_outputs/S28c.bytreatment_GEP_HR_p.pdf")
```

<embed src="figure4_outputs/S28c.bytreatment_GEP_HR_p.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
df_comb_IMPRES <- plot_bytreatment_survival(comb_clin,"IMPRES","IMPRES","figure4_outputs/S28d.bytreatment_IMPRES_HR_p.pdf",0.5,c(0.5,0.75,1,1.5,2,3,4,5))
knitr::include_graphics("figure4_outputs/S28d.bytreatment_IMPRES_HR_p.pdf")
```

<embed src="figure4_outputs/S28d.bytreatment_IMPRES_HR_p.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
df_comb_MHC_II_zscore <- plot_bytreatment_survival(comb_clin,"MHC_II_zscore","MHC II z-score","figure4_outputs/S28e.bytreatment_MHC_II_zscore_HR_p.pdf",0.5,c(0.5,0.75,1,1.5,2,3,4,5))
knitr::include_graphics("figure4_outputs/S28e.bytreatment_MHC_II_zscore_HR_p.pdf")
```

<embed src="figure4_outputs/S28e.bytreatment_MHC_II_zscore_HR_p.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
df_comb_TIDE <- plot_bytreatment_survival(comb_clin,"TIDE","TIDE","figure4_outputs/S28f.bytreatment_TIDE_HR_p.pdf",0.5,c(0.5,0.75,1,1.5,2,3,4,5))
knitr::include_graphics("figure4_outputs/S28f.bytreatment_TIDE_HR_p.pdf")
```

<embed src="figure4_outputs/S28f.bytreatment_TIDE_HR_p.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
df_comb_TCB <- plot_bytreatment_survival(comb_clin,"log10_rna_tcb","TCB","figure4_outputs/S28g.bytreatment_TCB_HR_p.pdf",0.5,c(0.25,0.5,0.75,1,1.5,2,3,4,5))
knitr::include_graphics("figure4_outputs/S28g.bytreatment_TCB_HR_p.pdf")
```

<embed src="figure4_outputs/S28g.bytreatment_TCB_HR_p.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
df_comb_BCB <- plot_bytreatment_survival(comb_clin,"log10_rna_bcb","BCB","figure4_outputs/S28h.bytreatment_BCB_HR_p.pdf",0.5,c(0.25,0.5,0.75,1,1.5,2,3,4,5))
knitr::include_graphics("figure4_outputs/S28h.bytreatment_BCB_HR_p.pdf")
```

<embed src="figure4_outputs/S28h.bytreatment_BCB_HR_p.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
df_comb_MAP4K1TBX3 <- plot_bytreatment_survival_double(comb_clin,"MAP4K1","MAP4K1","TBX3","TBX3","MAP4K1+TBX3","figure4_outputs/S28i.bytreatment_MAP4K1TBX3.pdf",0.55,0.5,c(0.5,0.75,1,1.5,2,3,4,5))
knitr::include_graphics("figure4_outputs/S28i.bytreatment_MAP4K1TBX3.pdf")
```

<embed src="figure4_outputs/S28i.bytreatment_MAP4K1TBX3.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
df_comb_MAP4K1AGER <- plot_bytreatment_survival_double(comb_clin,"MAP4K1","MAP4K1","AGER","AGER","MAP4K1+AGER","figure4_outputs/S28j.bytreatment_MAP4K1AGER.pdf",0.55,0.5,c(0.5,0.75,1,1.5,2,3,4,5))
knitr::include_graphics("figure4_outputs/S28j.bytreatment_MAP4K1AGER.pdf")
```

<embed src="figure4_outputs/S28j.bytreatment_MAP4K1AGER.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
df_comb_meanhighOSmeanlowOS <- plot_bytreatment_survival_double(comb_clin,"mean_highos_genes_z","High OS genes",
                                   "mean_lowos_genes_z","Low OS genes","Mean of high OS and low OS genes","figure4_outputs/S28k.bytreatment_meanhighOSmeanlowOS.pdf",0.01,1,c(0.25,0.5,0.75,1,1.5,2,3,4,5))
knitr::include_graphics("figure4_outputs/S28k.bytreatment_meanhighOSmeanlowOS.pdf")
```

<embed src="figure4_outputs/S28k.bytreatment_meanhighOSmeanlowOS.pdf" width="600px" height="600px" type="application/pdf" />

``` {.r}
bytreatment_table <- data.frame(df_comb_CD274$cohort,stringsAsFactors=F)
bytreatment_table[,c("CD274_HR","CD274_HR_lo","CD274_HR_hi","CD274_HR_p","CD274_LRT_p")] <- df_comb_CD274[,c("HR","HR_lo","HR_hi","HR_p","LR_p")]
bytreatment_table[,c("GEP_HR","GEP_HR_lo","GEP_HR_hi","GEP_HR_p","GEP_LRT_p")] <- df_comb_GEP[,c("HR","HR_lo","HR_hi","HR_p","LR_p")]
bytreatment_table[,c("CYT_HR","CYT_HR_lo","CYT_HR_hi","CYT_HR_p","CYT_LRT_p")] <- df_comb_CYT[,c("HR","HR_lo","HR_hi","HR_p","LR_p")]
bytreatment_table[,c("IMPRES_HR","IMPRES_HR_lo","IMPRES_HR_hi","IMPRES_HR_p","IMPRES_LRT_p")] <- df_comb_IMPRES[,c("HR","HR_lo","HR_hi","HR_p","LR_p")]
bytreatment_table[,c("TIDE_HR","TIDE_HR_lo","TIDE_HR_hi","TIDE_HR_p","TIDE_LRT_p")] <- df_comb_TIDE[,c("HR","HR_lo","HR_hi","HR_p","LR_p")]
bytreatment_table[,c("MHC_II_HR","MHC_II_HR_lo","MHC_II_HR_hi","MHC_II_HR_p","MHC_II_LRT_p")] <- df_comb_MHC_II_zscore[,c("HR","HR_lo","HR_hi","HR_p","LR_p")]
bytreatment_table[,c("TCB_HR","TCB_HR_lo","TCB_HR_hi","TCB_HR_p","TCB_LRT_p")] <- df_comb_TCB[,c("HR","HR_lo","HR_hi","HR_p","LR_p")]
bytreatment_table[,c("BCB_HR","BCB_HR_lo","BCB_HR_hi","BCB_HR_p","BCB_LRT_p")] <- df_comb_BCB[,c("HR","HR_lo","HR_hi","HR_p","LR_p")]
bytreatment_table[,c("MAP4K1TBX3_MAP4K1_HR","MAP4K1TBX3_MAP4K1_HR_lo","MAP4K1TBX3_MAP4K1_HR_hi","MAP4K1TBX3_MAP4K1_HR_p","MAP4K1TBX3_MAP4K1_LRT_p")] <-
                                                                                                        df_comb_MAP4K1TBX3[df_comb_MAP4K1TBX3$gene=="MAP4K1",c("HR","HR_lo","HR_hi","HR_p","LR_p")]
bytreatment_table[,c("MAP4K1TBX3_TBX3_HR","MAP4K1TBX3_TBX3_HR_lo","MAP4K1TBX3_TBX3_HR_hi","MAP4K1TBX3_TBX3_HR_p","MAP4K1TBX3_TBX3_LRT_p")] <-
                                                                                                        df_comb_MAP4K1TBX3[df_comb_MAP4K1TBX3$gene=="TBX3",c("HR","HR_lo","HR_hi","HR_p","LR_p")]

bytreatment_table[,c("MAP4K1AGER_MAP4K1_HR","MAP4K1AGER_MAP4K1_HR_lo","MAP4K1AGER_MAP4K1_HR_hi","MAP4K1AGER_MAP4K1_HR_p","MAP4K1AGER_MAP4K1_LRT_p")] <-
                                                                                                        df_comb_MAP4K1AGER[df_comb_MAP4K1AGER$gene=="MAP4K1",c("HR","HR_lo","HR_hi","HR_p","LR_p")]
bytreatment_table[,c("MAP4K1AGER_AGER_HR","MAP4K1AGER_AGER_HR_lo","MAP4K1AGER_AGER_HR_hi","MAP4K1AGER_AGER_HR_p","MAP4K1AGER_AGER_LRT_p")] <-
                                                                                                        df_comb_MAP4K1AGER[df_comb_MAP4K1AGER$gene=="AGER",c("HR","HR_lo","HR_hi","HR_p","LR_p")]

bytreatment_table[,c("OSmetagene_LongOS_HR","OSmetagene_LongOS_HR_lo","OSmetagene_LongOS_HR_hi","OSmetagene_LongOS_HR_p","OSmetagene_LongOS_LRT_p")] <-
                                                                        df_comb_meanhighOSmeanlowOS[df_comb_meanhighOSmeanlowOS$gene=="High OS genes",c("HR","HR_lo","HR_hi","HR_p","LR_p")]
bytreatment_table[,c("OSmetagene_ShortOS_HR","OSmetagene_ShortOS_HR_lo","OSmetagene_ShortOS_HR_hi","OSmetagene_ShortOS_HR_p","OSmetagene_ShortOS_LRT_p")] <-
                                                                        df_comb_meanhighOSmeanlowOS[df_comb_meanhighOSmeanlowOS$gene=="Low OS genes",c("HR","HR_lo","HR_hi","HR_p","LR_p")]

write.table(bytreatment_table,"figure4_outputs/CPBVal_bytreatment_performance.txt",quote=F,sep="\t")
```

Figure S29

``` {.r}
library(fgsea)

pcaexl_ccle_melanoma_tab <- readRDS("figure4_inputs/ccle_melanoma_tab.rds")

tbx3_corr_res <- data.frame(gene=gencode_coding$Gene,ensg=gencode_coding$Gene.id,chr=gencode_coding$Chromosome,
                                start=gencode_coding$Start,end=gencode_coding$End,pcaexcl_ccle_tbx3_cor=NA)
tbx3_corr_res <- tbx3_corr_res[tbx3_corr_res$gene!="TBX3",]

tbx3_ensg <- as.character(gencode_coding[gencode_coding$Gene=="TBX3",]$Gene.id)
tbx3_ccle_i <- which(rownames(pcaexl_ccle_melanoma_tab)==tbx3_ensg)

tbx3_corr_res$pcaexcl_ccle_tbx3_cor <- cor(t(as.matrix(pcaexl_ccle_melanoma_tab[-tbx3_ccle_i,])),as.numeric(pcaexl_ccle_melanoma_tab[tbx3_ensg,]))[,1]
```

    ## Warning in cor(t(as.matrix(pcaexl_ccle_melanoma_tab[-tbx3_ccle_i, ])),
    ## as.numeric(pcaexl_ccle_melanoma_tab[tbx3_ensg, : the standard deviation is
    ## zero

``` {.r}
tbx3_corr_res$pcaexcl_ccle_tbx3_cor_z <- scale(tbx3_corr_res$pcaexcl_ccle_tbx3_cor)
tmp_tbx3_corr_res <- tbx3_corr_res[!is.na(tbx3_corr_res$pcaexcl_ccle_tbx3_cor_z)&tbx3_corr_res$gene!="TBX3",]

s_de_res <- tmp_tbx3_corr_res[order(tmp_tbx3_corr_res$pcaexcl_ccle_tbx3_cor_z),]
s_de_res <- s_de_res[!is.na(s_de_res$pcaexcl_ccle_tbx3_cor_z),]
rownames(s_de_res) <- s_de_res$ensg
s_de_res_tmpRanks <- s_de_res[!is.na(gencode_coding[rownames(s_de_res),]$entrez_id),]$pcaexcl_ccle_tbx3_cor_z
names(s_de_res_tmpRanks) <- gencode_coding[rownames(s_de_res[!is.na(gencode_coding[rownames(s_de_res),]$entrez_id),]),]$entrez_id

pathways.go <- gmtPathways("figure4_inputs/c5.all.v6.2.entrez.gmt")
s_de_res_tmpRanks_go <- fgsea(pathways = pathways.go, stats = s_de_res_tmpRanks,minSize=10, maxSize=500,  nperm=100000)
#write.table(data.frame(arrange(s_de_res_tmpRanks_go[s_de_res_tmpRanks_go$padj<0.01,],NES))[,c(1:7)],"figure4_outputs/fcmeta_pcaexcl_ccle.TBX3corr_GOterms.txt",quote=F,sep="\t")
topPathwaysUp <- arrange(s_de_res_tmpRanks_go[padj<0.01&NES<0,],NES)[1:20]$pathway
topPathwaysDown <- arrange(s_de_res_tmpRanks_go[padj<0.01&NES>0,],-NES)[1:20]$pathway
pdf("figure4_outputs/S29a.i.ccle.gsea_GOterms.highTBX3_enriched_pathways.pdf",height=7,width=16)
plotGseaTable(pathways.go[topPathwaysDown[!is.na(topPathwaysDown)]], s_de_res_tmpRanks, s_de_res_tmpRanks_go, gseaParam = 0.5)
dev.off()
```

    ## pdf 
    ##   2

``` {.r}
knitr::include_graphics("figure4_outputs/S29a.i.ccle.gsea_GOterms.highTBX3_enriched_pathways.pdf")
```

<embed src="figure4_outputs/S29a.i.ccle.gsea_GOterms.highTBX3_enriched_pathways.pdf" width="700px" height="700px" type="application/pdf" />

``` {.r}
pdf("figure4_outputs/S29a.ii.ccle.gsea_GOterms.lowTBX3_enriched_pathways.pdf",height=7,width=16)
plotGseaTable(pathways.go[topPathwaysUp[!is.na(topPathwaysUp)]], s_de_res_tmpRanks, s_de_res_tmpRanks_go, gseaParam = 0.5)
dev.off()
```

    ## pdf 
    ##   2

``` {.r}
knitr::include_graphics("figure4_outputs/S29a.ii.ccle.gsea_GOterms.lowTBX3_enriched_pathways.pdf")
```

<embed src="figure4_outputs/S29a.ii.ccle.gsea_GOterms.lowTBX3_enriched_pathways.pdf" width="700px" height="700px" type="application/pdf" />


``` {.r}
color20 <- c("#e6194b","#3cb44b","#ffe119","#4363d8","#f58231","#911eb4","#46f0f0","#f032e6","#bcf60c","#fabebe","#008080",
             "#e6beff","#9a6324","#fffac8","#800000","#aaffc3","#808000","#ffd8b1","#000075","#808080","#ffffff","#000000")

cbbPalette <- c("#000000","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7")

jerby_tumor_tsne <- readRDS("figure4_inputs/jerby_tumor_tsne.rds")
jerby_normal_tsne <- readRDS("figure4_inputs/jerby_normal_tsne.rds")
jerby_tpm <- readRDS("figure4_inputs/jerby_tpm_mini.rds")

jerby_tumor_tsne$TBX3 <- as.numeric(jerby_tpm[which(jerby_tpm$GENE=="TBX3"),jerby_tumor_tsne$NAME])
jerby_tumor_tsne$NGFR <- as.numeric(jerby_tpm[which(jerby_tpm$GENE=="NGFR"),jerby_tumor_tsne$NAME])
jerby_normal_tsne$TBX3 <- as.numeric(jerby_tpm[which(jerby_tpm$GENE=="TBX3"),jerby_normal_tsne$NAME])
jerby_normal_tsne$NGFR <- as.numeric(jerby_tpm[which(jerby_tpm$GENE=="NGFR"),jerby_normal_tsne$NAME])

g_jerby_normal_celltype <- ggplot(jerby_normal_tsne,aes(x=X,y=Y,color=cell.type)) + geom_point(size=1) + xlab("tSNE1") + ylab("tSNE2") + theme_classic() + scale_color_manual(values=cbbPalette) + theme(legend.position="bottom")
ggsave("figure4_outputs/S29b.jerby_normal_celltype.pdf",g_jerby_normal_celltype,height=7,width=7)
knitr::include_graphics("figure4_outputs/S29b.jerby_normal_celltype.pdf")
```

<embed src="figure4_outputs/S29b.jerby_normal_celltype.pdf" width="500px" height="500px" type="application/pdf" />

``` {.r}
g_jerby_normal_TBX3_zero <- ggplot(jerby_normal_tsne,aes(x=X,y=Y,color=TBX3>0)) + geom_point(size=1) + xlab("tSNE1") + ylab("tSNE2") + theme_classic() + scale_color_manual(values=c("black","red")) + theme(legend.position="bottom")
ggsave("figure4_outputs/S29c.jerby_normal_TBX3_zero_redblack.pdf",g_jerby_normal_TBX3_zero,height=7,width=7)
knitr::include_graphics("figure4_outputs/S29c.jerby_normal_TBX3_zero_redblack.pdf")
```

<embed src="figure4_outputs/S29c.jerby_normal_TBX3_zero_redblack.pdf" width="500px" height="500px" type="application/pdf" />

``` {.r}
g_jerby_tumor_bypatient <- ggplot(jerby_tumor_tsne,aes(x=X,y=Y,color=LABELS)) + geom_point(size=1) + xlab("tSNE1") + ylab("tSNE2") + theme_classic() + scale_color_manual(values=color20) + theme(legend.position="bottom")
ggsave("figure4_outputs/S29d.jerby_tumor_bypatient.pdf",g_jerby_tumor_bypatient,height=7,width=7)
knitr::include_graphics("figure4_outputs/S29d.jerby_tumor_bypatient.pdf")
```

<embed src="figure4_outputs/S29d.jerby_tumor_bypatient.pdf" width="500px" height="500px" type="application/pdf" />

``` {.r}
g_jerby_tumor_TBX3_zero <- ggplot(jerby_tumor_tsne,aes(x=X,y=Y,color=TBX3>0)) + geom_point(size=1) + xlab("tSNE1") + ylab("tSNE2") + theme_classic() + scale_color_manual(values=c("black","red")) + theme(legend.position="bottom")
ggsave("figure4_outputs/S29e.jerby_tumor_TBX3_zero_redblack.pdf",g_jerby_tumor_TBX3_zero,height=7,width=7)
knitr::include_graphics("figure4_outputs/S29e.jerby_tumor_TBX3_zero_redblack.pdf")
```

<embed src="figure4_outputs/S29e.jerby_tumor_TBX3_zero_redblack.pdf" width="500px" height="500px" type="application/pdf" />

``` {.r}
jerby_tumor_tsne$cell.type <- "tumor"
jerby_tn_tbx3 <- rbindlist(list(jerby_normal_tsne[,c("NAME","tumor","cell.type","TBX3","NGFR")],jerby_tumor_tsne[,c("NAME","tumor","cell.type","TBX3","NGFR")]))
jerby_tn_tbx3$cell.type_TBX3 <- factor(jerby_tn_tbx3$cell.type,levels=arrange(jerby_tn_tbx3[,mean(TBX3),by=cell.type],V1)$cell.type)
jerby_tn_tbx3$cell_type_TBX3NGFR <- ifelse(jerby_tn_tbx3$cell.type=="tumor",ifelse(jerby_tn_tbx3$NGFR==0,"tumor_NGFRzero","tumor_NGFRexpressed"),jerby_tn_tbx3$cell.type)
jerby_tn_tbx3$cell_type_TBX3NGFR <- factor(jerby_tn_tbx3$cell_type_TBX3NGFR,levels=arrange(jerby_tn_tbx3[,mean(TBX3),by=cell_type_TBX3NGFR],V1)$cell_type_TBX3NGFR)


tmp_wt <- wilcox.test(jerby_tn_tbx3[jerby_tn_tbx3$cell.type=="tumor",]$TBX3,jerby_tn_tbx3[jerby_tn_tbx3$cell.type!="tumor",]$TBX3)
tmp_wt_p <- formatC(tmp_wt$p.value,format="e",digits=2)
tmp_wt_2 <- wilcox.test(jerby_tn_tbx3[jerby_tn_tbx3$cell.type=="tumor"&jerby_tn_tbx3$NGFR>0,]$TBX3,jerby_tn_tbx3[jerby_tn_tbx3$cell.type=="tumor"&jerby_tn_tbx3$NGFR==0,]$TBX3)
tmp_wt_2_p <- formatC(tmp_wt_2$p.value,format="e",digits=2)

jerby_tumornormal_TBX3_NGFR_celltype_boxplot <- ggplot(jerby_tn_tbx3,aes(x=cell_type_TBX3NGFR,y=TBX3)) + geom_boxplot(outlier.shape = NA) + geom_jitter() + theme_classic() + xlab("") + theme(legend.position="none") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + stat_summary(fun.y=mean, geom="point", shape=18,size=3, color="red") + ggtitle(paste0("Tumor vs. Normal Wilcox p =",tmp_wt_p,"\n","Tumor NGFR positive vs. Tumor NGFR negative Wilcox p =",tmp_wt_2_p))
ggsave("figure4_outputs/S29f.jerby_tumornormal_TBX3_withNGFR_boxplot_withmeans.pdf",jerby_tumornormal_TBX3_NGFR_celltype_boxplot,height=7,width=7)
knitr::include_graphics("figure4_outputs/S29f.jerby_tumornormal_TBX3_withNGFR_boxplot_withmeans.pdf")
```

<embed src="figure4_outputs/S29f.jerby_tumornormal_TBX3_withNGFR_boxplot_withmeans.pdf" width="500px" height="500px" type="application/pdf" />


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
    ##  [1] fgsea_1.2.1        Rcpp_1.0.7         RAM_1.2.1.3       
    ##  [4] vegan_2.4-3        lattice_0.20-35    permute_0.9-4     
    ##  [7] gridExtra_2.3      data.table_1.14.0  psych_1.7.5       
    ## [10] RColorBrewer_1.1-2 pROC_1.9.1         verification_1.42 
    ## [13] dtw_1.18-1         proxy_0.4-22       CircStats_0.2-6   
    ## [16] MASS_7.3-50        boot_1.3-20        fields_8.10       
    ## [19] maps_3.3.0         spam_2.2-0         dotCall64_0.9-5.2 
    ## [22] pheatmap_1.0.10    ggrepel_0.9.1      plyr_1.8.6        
    ## [25] meta_4.9-2         ggsignif_0.2.0     corrplot_0.84     
    ## [28] reshape2_1.4.3     ggmosaic_0.2.0     survminer_0.4.3   
    ## [31] ggpubr_0.2         magrittr_2.0.1     survival_2.42-6   
    ## [34] ggplot2_3.2.1     
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] backports_1.1.2         fastmatch_1.1-0        
    ##   [3] igraph_1.2.1            lazyeval_0.2.2         
    ##   [5] sp_1.3-1                splines_3.4.0          
    ##   [7] BiocParallel_1.10.1     digest_0.6.27          
    ##   [9] htmltools_0.5.2         gdata_2.18.0           
    ##  [11] productplots_0.1.1      fansi_0.5.0            
    ##  [13] phytools_0.6-00         cluster_2.0.7-1        
    ##  [15] jpeg_0.1-8              colorspace_1.4-1       
    ##  [17] dplyr_0.8.4             crayon_1.4.1           
    ##  [19] jsonlite_1.7.2          zoo_1.8-3              
    ##  [21] phangorn_2.2.0          ape_5.1                
    ##  [23] glue_1.4.2              gtable_0.3.0           
    ##  [25] abind_1.4-5             scales_0.5.0           
    ##  [27] futile.options_1.0.1    mvtnorm_1.0-8          
    ##  [29] plotrix_3.7-2           viridisLite_0.3.0      
    ##  [31] xtable_1.8-2            cmprsk_2.2-7           
    ##  [33] magic_1.5-8             foreign_0.8-70         
    ##  [35] mapproj_1.2.6           km.ci_0.5-2            
    ##  [37] animation_2.5           htmlwidgets_1.5.4      
    ##  [39] httr_1.4.2              gplots_3.0.1           
    ##  [41] geosphere_1.5-5         ellipsis_0.3.2         
    ##  [43] pkgconfig_2.0.3         reshape_0.8.8          
    ##  [45] utf8_1.2.2              tidyselect_1.0.0       
    ##  [47] labeling_0.3            rlang_0.4.11           
    ##  [49] munsell_0.5.0           tools_3.4.0            
    ##  [51] ggmap_2.6.1             ade4_1.7-11            
    ##  [53] broom_0.5.0             evaluate_0.11          
    ##  [55] geometry_0.3-6          stringr_1.3.1          
    ##  [57] fastmap_1.1.0           FD_1.0-12              
    ##  [59] yaml_2.2.1              knitr_1.16             
    ##  [61] survMisc_0.5.4          caTools_1.17.1         
    ##  [63] purrr_0.3.3             RgoogleMaps_1.4.2      
    ##  [65] nlme_3.1-137            formatR_1.5            
    ##  [67] compiler_3.4.0          plotly_4.9.4.1         
    ##  [69] png_0.1-7               clusterGeneration_1.3.4
    ##  [71] tibble_2.1.3            stringi_1.2.3          
    ##  [73] futile.logger_1.4.3     Matrix_1.2-14          
    ##  [75] KMsurv_0.1-5            vctrs_0.3.8            
    ##  [77] labdsv_1.8-0            msm_1.6.4              
    ##  [79] pillar_1.6.2            lifecycle_1.0.0        
    ##  [81] combinat_0.0-8          cowplot_0.9.3          
    ##  [83] bitops_1.0-6            R6_2.5.1               
    ##  [85] KernSmooth_2.23-15      lambda.r_1.2.3         
    ##  [87] gtools_3.8.1            assertthat_0.2.0       
    ##  [89] proto_1.0.0             rprojroot_1.3-2        
    ##  [91] rjson_0.2.20            withr_2.4.2            
    ##  [93] mnormt_1.5-5            mgcv_1.8-17            
    ##  [95] expm_0.999-2            parallel_3.4.0         
    ##  [97] quadprog_1.5-5          VennDiagram_1.6.20     
    ##  [99] tidyr_1.0.2             coda_0.19-1            
    ## [101] rmarkdown_1.5           numDeriv_2016.8-1      
    ## [103] scatterplot3d_0.3-41
