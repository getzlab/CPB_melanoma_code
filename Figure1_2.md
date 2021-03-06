Figure1\_2
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

makeKMplot_2group <- function(df,var1,png_name) {
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
png(png_name,height=700,width=700,res=100)
#fit<- do.call(survfit,list(formula = surv ~ class, data = df))
fit<- surv_fit(surv ~ class, data = df)
print(ggsurvplot(fit, data = df,pval=TRUE,pval.method=TRUE,legend=c(0.75,0.75),censor=TRUE,palette=c("blue","red"),title=hr_string))
dev.off()
}

makeKMplot_2group_withp <- function(df,var1,png_name) {
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
hr_string <- paste0("HR=",formatSF(hr,digits=3)," (95% CI, ",formatSF(hr_lo,digits=3),"-",formatSF(hr_hi,digits=3),")",
" Log-rank p=",as.character(formatSF(s_coxph$sctest["pvalue"],digits=3,format="e")),"\n",
names(class_table)[1],"=",class_table[[1]]," ",names(class_table)[2],"=",class_table[[2]])
png(png_name,height=700,width=700,res=100)
#fit<- do.call(survfit,list(formula = surv ~ class, data = df))
fit<- surv_fit(surv ~ class, data = df)
print(ggsurvplot(fit, data = df,pval=TRUE,pval.method=TRUE,legend=c(0.75,0.75),censor=TRUE,palette=c("blue","red"),title=hr_string))
dev.off()
}

makeKMplot_4group <- function(df,var1,png_name) {
#note: need to pass in clin[!is.na(clin$surv),] or class counts will be wrong
df$class <- df[,var1]
class_table <- table(df$class)
count_string <- paste0(names(class_table)[1],"=",class_table[[1]]," ",names(class_table)[2],"=",class_table[[2]],"\n",
names(class_table)[3],"=",class_table[[3]]," ",names(class_table)[4],"=",class_table[[4]])
png(png_name,height=700,width=700,res=100)
#fit<- do.call(survfit,list(formula = surv ~ class, data = df))
fit<- surv_fit(surv ~ class, data = df)
print(ggsurvplot(fit, data = df,pval=TRUE,pval.method=TRUE,legend=c(0.75,0.75),censor=TRUE,palette=cbpalette[1:4],title=count_string))
dev.off()
}

makeKMplot_4group_withsmalltitle <- function(df,var1,png_name) {
#note: need to pass in clin[!is.na(clin$surv),] or class counts will be wrong
df$class <- df[,var1]
class_table <- table(df$class)
count_string <- paste0(names(class_table)[1],"=",class_table[[1]]," ",names(class_table)[2],"=",class_table[[2]],"\n",
names(class_table)[3],"=",class_table[[3]]," ",names(class_table)[4],"=",class_table[[4]])
png(png_name,height=700,width=700,res=100)
#fit<- do.call(survfit,list(formula = surv ~ class, data = df))
fit<- surv_fit(surv ~ class, data = df)
print(ggsurvplot(fit, data = df,pval=TRUE,pval.method=TRUE,legend=c(0.75,0.75),censor=TRUE,palette=cbpalette[1:4],title=count_string,font.title=c(10,"plain","black")))
dev.off()
}

makeResponseMosaic <- function(df,var1,png_name) {
f_ft <- fisher.test(table(df$Response,df[,var1]))
f_ft_pstring <- paste0("Fisher p=",formatSF(f_ft$p.value,digits=3))
tdf <- melt(table(df$Response,df[,var1]))
g_rnr_temp <- ggplot(tdf) + geom_mosaic(aes(weight=value,x=product(Var2),fill=Var1),offset=0.01) + 
scale_fill_manual(values=c("red","blue")) + theme_bw() +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.title=element_blank(),axis.text.x=element_text(size=12),axis.text.y=element_text(size=12)) +
xlab("") + ylab("") + annotate(geom="text",x=0.5,y=1.02,label=f_ft_pstring) + scale_y_continuous(breaks=seq(0,1,by=0.25))
ggsave(png_name,g_rnr_temp + geom_text(data = ggplot_build(g_rnr_temp)$data[[1]], aes(x = (xmin+xmax)/2, y = (ymin+ymax)/2, label=.wt)),height=7,width=7)
}

makeResponseMosaic_rotatelabels <- function(df,var1,png_name) {
f_ft <- fisher.test(table(df$Response,df[,var1]))
f_ft_pstring <- paste0("Fisher p=",formatSF(f_ft$p.value,digits=3))
tdf <- melt(table(df$Response,df[,var1]))
g_rnr_temp <- ggplot(tdf) + geom_mosaic(aes(weight=value,x=product(Var2),fill=Var1),offset=0.01) + scale_fill_manual(values=c("red","blue")) + theme_bw() +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.title=element_blank(),axis.text.x=element_text(size=12,angle=90),axis.text.y=element_text(size=12)) +
xlab("") + ylab("") + annotate(geom="text",x=0.5,y=1.02,label=f_ft_pstring) + scale_y_continuous(breaks=seq(0,1,by=0.25))
ggsave(png_name,g_rnr_temp + geom_text(data = ggplot_build(g_rnr_temp)$data[[1]], aes(x = (xmin+xmax)/2, y = (ymin+ymax)/2, label=.wt)),height=7,width=7)
}

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

Figure 1

``` {.r}
dna_clin$tmb_binary <- ifelse(dna_clin$nonsilent_snpdnpindel>median(dna_clin$nonsilent_snpdnpindel),"high_tmb","low_tmb")
dna_clin$plot_tmb_binary <- ifelse(dna_clin$nonsilent_snpdnpindel>median(dna_clin$nonsilent_snpdnpindel),"High TMB","Low TMB")

dir.create("figure12_outputs", showWarnings = FALSE)
makeKMplot_2group(dna_clin,"plot_tmb_binary","figure12_outputs/1b.1.dna_tmb_survival.png")
```

    ## png 
    ##   2

``` {.r}
knitr::include_graphics("figure12_outputs/1b.1.dna_tmb_survival.png")
```

<img src="figure12_outputs/1b.1.dna_tmb_survival.png" width="75%" height="75%" />

``` {.r}
tmpwt <- wilcox.test(dna_clin[dna_clin$Response=="responder",]$nonsilent_snpdnpindel,dna_clin[dna_clin$Response=="nonresponder",]$nonsilent_snpdnpindel)
g_tmb_resp <- ggplot(dna_clin,aes(x=Response,y=nonsilent_snpdnpindel)) + geom_boxplot(aes(fill=Response),outlier.shape=NA) + geom_jitter() + theme_classic() + ylab("log10(TMB)") +
xlab("") + scale_fill_manual(values=c("red","blue")) + annotate("text",x=1.5,y=20,label=paste0("Wilcox p = ",formatSF(tmpwt$p.value,digits=3,format="e"))) + scale_y_log10() +
ggtitle(paste0("NR=",nrow(dna_clin[dna_clin$Response=="nonresponder",])," R=",nrow(dna_clin[dna_clin$Response=="responder",])))
ggsave("figure12_outputs/1b.2.dna_tmb_response.png",g_tmb_resp,height=7,width=7)
knitr::include_graphics("figure12_outputs/1b.2.dna_tmb_response.png")
```

<img src="figure12_outputs/1b.2.dna_tmb_response.png" width="75%" height="75%" />

``` {.r}
dna_clin$plot_purity_binary <- NA
dna_clin[dna_clin$not_cell_line,]$plot_purity_binary <- ifelse(dna_clin[dna_clin$not_cell_line,]$purity>median(dna_clin[dna_clin$not_cell_line,]$purity),"High tumor purity","Low tumor purity")
makeKMplot_2group(dna_clin[dna_clin$not_cell_line,],"plot_purity_binary","figure12_outputs/1c.1.dna_purity_survival.png")
```

    ## png 
    ##   2

``` {.r}
knitr::include_graphics("figure12_outputs/1c.1.dna_purity_survival.png")
```

<img src="figure12_outputs/1c.1.dna_purity_survival.png" width="75%" height="75%" />

``` {.r}
tmpwt <- wilcox.test(dna_clin[dna_clin$not_cell_line&dna_clin$Response=="responder",]$purity,dna_clin[dna_clin$not_cell_line&dna_clin$Response=="nonresponder",]$purity)
g_purity_resp <- ggplot(dna_clin[dna_clin$not_cell_line,],aes(x=Response,y=purity)) + geom_boxplot(aes(fill=Response),outlier.shape=NA) + geom_jitter() + theme_classic() +
ylab("Tumor Purity") + xlab("") + scale_fill_manual(values=c("red","blue")) + annotate("text",x=1.5,y=0.9,label=paste0("Wilcox p = ",formatSF(tmpwt$p.value,digits=3,format="e"))) +
ggtitle(paste0("NR=",nrow(dna_clin[dna_clin$not_cell_line&dna_clin$Response=="nonresponder",])," R=",nrow(dna_clin[dna_clin$not_cell_line&dna_clin$Response=="responder",])))
ggsave("figure12_outputs/1c.2.dna_purity_response.png",g_purity_resp,height=7,width=7)
knitr::include_graphics("figure12_outputs/1c.2.dna_purity_response.png")
```

<img src="figure12_outputs/1c.2.dna_purity_response.png" width="75%" height="75%" />

``` {.r}
ct_r <- formatSF(cor.test(dna_clin[dna_clin$not_cell_line,]$purity,log10(dna_clin[dna_clin$not_cell_line,]$nonsilent_snpdnpindel),method="spearman")$estimate,format="f",digits=3)
```

    ## Warning in cor.test.default(dna_clin[dna_clin$not_cell_line, ]$purity,
    ## log10(dna_clin[dna_clin$not_cell_line, : Cannot compute exact p-value with
    ## ties

``` {.r}
ct_p <- formatSF(cor.test(dna_clin[dna_clin$not_cell_line,]$purity,log10(dna_clin[dna_clin$not_cell_line,]$nonsilent_snpdnpindel),method="spearman")$p.value,format="e",digits=3)
```

    ## Warning in cor.test.default(dna_clin[dna_clin$not_cell_line, ]$purity,
    ## log10(dna_clin[dna_clin$not_cell_line, : Cannot compute exact p-value with
    ## ties

``` {.r}
g_dna_purity_tmb <- ggplot(dna_clin[dna_clin$not_cell_line,],aes(x=purity,y=nonsilent_snpdnpindel)) + geom_point(size=3) + geom_smooth(method="lm",se=FALSE) + theme_classic() + scale_y_log10() +
xlab("Tumor Purity") + ylab("TMB (log10)") + annotate(geom="text",x=0.15,y=3000,label=paste0("rho = ",ct_r,"\np = ",ct_p,"\nn = ",nrow(dna_clin[dna_clin$not_cell_line,])))
ggsave("figure12_outputs/1d.dna_purity_tmb.png",g_dna_purity_tmb,height=7,width=7)
knitr::include_graphics("figure12_outputs/1d.dna_purity_tmb.png")
```

<img src="figure12_outputs/1d.dna_purity_tmb.png" width="75%" height="75%" />

``` {.r}
dna_clin$plot_combined_tmb_purity <- NA
dna_clin[dna_clin$not_cell_line,]$plot_combined_tmb_purity <- paste0(dna_clin[dna_clin$not_cell_line,]$plot_tmb_binary,", ",dna_clin[dna_clin$not_cell_line,]$plot_purity_binary)
dna_clin$plot_combined_tmb_purity_bestrest <- NA
dna_clin[dna_clin$not_cell_line,]$plot_combined_tmb_purity_bestrest <- ifelse(dna_clin[dna_clin$not_cell_line,]$plot_combined_tmb_purity=="High TMB, Low tumor purity",
"High TMB, Low tumor purity","Others")

makeKMplot_2group(dna_clin[dna_clin$not_cell_line,],"plot_combined_tmb_purity_bestrest","figure12_outputs/1e.dna_tmb_purity_survival.png")
```

    ## png 
    ##   2

``` {.r}
knitr::include_graphics("figure12_outputs/1e.dna_tmb_purity_survival.png")
```

<img src="figure12_outputs/1e.dna_tmb_purity_survival.png" width="75%" height="75%" />

``` {.r}
makeResponseMosaic(dna_clin[dna_clin$not_cell_line,],"plot_combined_tmb_purity_bestrest","figure12_outputs/1f.dna_tmbpurity_bestrest_response.png")
knitr::include_graphics("figure12_outputs/1f.dna_tmbpurity_bestrest_response.png")
```

<img src="figure12_outputs/1f.dna_tmbpurity_bestrest_response.png" width="75%" height="75%" />

``` {.r}
makeKMplot_2group_withp(dna_clin,"Response","figure12_outputs/S1b.dna_response_survival.png")
```

    ## png 
    ##   2

``` {.r}
knitr::include_graphics("figure12_outputs/S1b.dna_response_survival.png")
```

<img src="figure12_outputs/S1b.dna_response_survival.png" width="75%" height="75%" />

``` {.r}
makeKMplot_2group_withp(clin[!is.na(clin$surv),],"Response","figure12_outputs/S1c.rna_response_survival.png")
```

    ## png 
    ##   2

``` {.r}
knitr::include_graphics("figure12_outputs/S1c.rna_response_survival.png")
```

<img src="figure12_outputs/S1c.rna_response_survival.png" width="75%" height="75%" />

``` {.r}
cbpalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442")
cbpalette_reordered <- c("#E69F00","#999999","#56B4E9","#009E73","#F0E442")

dna_clin$plot_tmb_over10 <- ifelse(10^6*dna_clin$nonsilent_snpdnpindel/28665628>10,"TMB > 10 mutations/Mb","TMB <= 10 mutations/Mb")

makeKMplot_2group(dna_clin,"plot_tmb_over10","figure12_outputs/S5a.dna_tmbover10_survival.png")
```

    ## png 
    ##   2

``` {.r}
knitr::include_graphics("figure12_outputs/S5a.dna_tmbover10_survival.png")
```

<img src="figure12_outputs/S5a.dna_tmbover10_survival.png" width="75%" height="75%" />

``` {.r}
makeResponseMosaic(dna_clin,"plot_tmb_over10","figure12_outputs/S5b.dna_tmbover10_response.png")
knitr::include_graphics("figure12_outputs/S5b.dna_tmbover10_response.png")
```

<img src="figure12_outputs/S5b.dna_tmbover10_response.png" width="75%" height="75%" />

``` {.r}
makeKMplot_4group(dna_clin[dna_clin$not_cell_line,],"plot_combined_tmb_purity","figure12_outputs/S5c.dna_purity_tmb_survival.png")
```

    ## png 
    ##   2

``` {.r}
knitr::include_graphics("figure12_outputs/S5c.dna_purity_tmb_survival.png")
```

<img src="figure12_outputs/S5c.dna_purity_tmb_survival.png" width="75%" height="75%" />

``` {.r}
ct_r <- formatSF(cor.test(log10(dna_clin$n_neoag_underrank2_snp_9mer10mer),log10(dna_clin$nonsilent_snpdnpindel),method="spearman")$estimate,format="f",digits=3)
```

    ## Warning in
    ## cor.test.default(log10(dna_clin$n_neoag_underrank2_snp_9mer10mer), : Cannot
    ## compute exact p-value with ties

``` {.r}
ct_p <- formatSF(cor.test(log10(dna_clin$n_neoag_underrank2_snp_9mer10mer),log10(dna_clin$nonsilent_snpdnpindel),method="spearman")$p.value,format="e",digits=3)
```

    ## Warning in
    ## cor.test.default(log10(dna_clin$n_neoag_underrank2_snp_9mer10mer), : Cannot
    ## compute exact p-value with ties

``` {.r}
g_dna_tmb_neoantigens <- ggplot(dna_clin,aes(x=nonsilent_snpdnpindel,y=n_neoag_underrank2_snp_9mer10mer)) + geom_point(size=3) + geom_smooth(method="lm",se=FALSE) + theme_classic() +
scale_x_log10() + scale_y_log10() +  xlab("TMB") + ylab("neoantigens") +
annotate(geom="text",x=100,y=3000,label=paste0("rho = ",ct_r,"\np = ",ct_p,"\nn = ",nrow(dna_clin)))
ggsave("figure12_outputs/S5d.dna_tmb_neoantigens.png",g_dna_tmb_neoantigens,height=7,width=7)
knitr::include_graphics("figure12_outputs/S5d.dna_tmb_neoantigens.png")
```

<img src="figure12_outputs/S5d.dna_tmb_neoantigens.png" width="75%" height="75%" />

``` {.r}
dna_clin$neoag_binary <- ifelse(dna_clin$n_neoag_underrank2_snp_9mer10mer>median(dna_clin$n_neoag_underrank2_snp_9mer10mer),"High neoantigen burden","Low neoantigen burden")
makeKMplot_2group(dna_clin,"neoag_binary","figure12_outputs/S5e.dna_neoantigen_survival.png")
```

    ## png 
    ##   2

``` {.r}
knitr::include_graphics("figure12_outputs/S5e.dna_neoantigen_survival.png")
```

<img src="figure12_outputs/S5e.dna_neoantigen_survival.png" width="75%" height="75%" />

``` {.r}
ct_r <- formatSF(cor.test(log10(dna_clin$nonsilent_clonal_snpdnpindel),log10(dna_clin$nonsilent_snpdnpindel),method="spearman")$estimate,format="f",digits=3)
```

    ## Warning in cor.test.default(log10(dna_clin$nonsilent_clonal_snpdnpindel), :
    ## Cannot compute exact p-value with ties

``` {.r}
ct_p <- formatSF(cor.test(log10(dna_clin$nonsilent_clonal_snpdnpindel),log10(dna_clin$nonsilent_snpdnpindel),method="spearman")$p.value,format="e",digits=3)
```

    ## Warning in cor.test.default(log10(dna_clin$nonsilent_clonal_snpdnpindel), :
    ## Cannot compute exact p-value with ties

``` {.r}
g_dna_tmb_clonaltmb <- ggplot(dna_clin,aes(x=nonsilent_snpdnpindel,y=nonsilent_clonal_snpdnpindel)) + geom_point(size=3) + geom_smooth(method="lm",se=FALSE) + theme_classic() +
scale_x_log10() + scale_y_log10() +  xlab("TMB") + ylab("Clonal TMB") +
annotate(geom="text",x=100,y=3000,label=paste0("rho = ",ct_r,"\np = ",ct_p,"\nn = ",nrow(dna_clin)))
ggsave("figure12_outputs/S5f.dna_tmb_clonaltmb.png",g_dna_tmb_clonaltmb,height=7,width=7)
knitr::include_graphics("figure12_outputs/S5f.dna_tmb_clonaltmb.png")
```

<img src="figure12_outputs/S5f.dna_tmb_clonaltmb.png" width="75%" height="75%" />

``` {.r}
dna_clin$clonal_tmb_binary <- ifelse(dna_clin$nonsilent_clonal_snpdnpindel>median(dna_clin$nonsilent_clonal_snpdnpindel),"High Clonal TMB","Low Clonal TMB")
makeKMplot_2group(dna_clin,"clonal_tmb_binary","figure12_outputs/S5g.dna_clonal_tmb_survival.png")
```

    ## png 
    ##   2

``` {.r}
knitr::include_graphics("figure12_outputs/S5g.dna_clonal_tmb_survival.png")
```

<img src="figure12_outputs/S5g.dna_clonal_tmb_survival.png" width="75%" height="75%" />

Figure 2

``` {.r}
clin$plot_tcb_binary <- ifelse(clin$rna_tcb>median(clin$rna_tcb),"High TCB","Low TCB")
makeKMplot_2group(clin[!is.na(clin$surv),],"plot_tcb_binary","figure12_outputs/2a.1.rna_tcb_survival.png")
```

    ## png 
    ##   2

``` {.r}
knitr::include_graphics("figure12_outputs/2a.1.rna_tcb_survival.png")
```

<img src="figure12_outputs/2a.1.rna_tcb_survival.png" width="75%" height="75%" />

``` {.r}
tmpwt <- wilcox.test(log10(clin[clin$Response=="responder",]$rna_tcb),log10(clin[clin$Response=="nonresponder",]$rna_tcb))
g_rna_tcb_resp <- ggplot(clin,aes(x=Response,y=rna_tcb)) + geom_boxplot(aes(fill=Response),outlier.shape=NA) + geom_jitter() + theme_classic() + ylab(expression(TCB[RNA])) +
xlab("") + scale_fill_manual(values=c("red","blue")) + annotate("text",x=1.5,y=20,label=paste0("Wilcox p = ",formatSF(tmpwt$p.value,digits=3,format="e"))) + scale_y_log10() +
ggtitle(paste0("NR=",nrow(clin[clin$Response=="nonresponder",])," R=",nrow(clin[clin$Response=="responder",])))
ggsave("figure12_outputs/2a.2.rna_tcb_response.png",g_rna_tcb_resp,height=7,width=7)
knitr::include_graphics("figure12_outputs/2a.2.rna_tcb_response.png")
```

<img src="figure12_outputs/2a.2.rna_tcb_response.png" width="75%" height="75%" />

``` {.r}
clin$plot_bcb_binary <- ifelse(clin$rna_bcb>median(clin$rna_bcb),"High BCB","Low BCB")
makeKMplot_2group(clin[!is.na(clin$surv),],"plot_bcb_binary","figure12_outputs/2b.1.rna_bcb_survival.png")
```

    ## png 
    ##   2

``` {.r}
knitr::include_graphics("figure12_outputs/2b.1.rna_bcb_survival.png")
```

<img src="figure12_outputs/2b.1.rna_bcb_survival.png" width="75%" height="75%" />

``` {.r}
tmpwt <- wilcox.test(log10(clin[clin$Response=="responder",]$rna_bcb),log10(clin[clin$Response=="nonresponder",]$rna_bcb))
g_rna_bcb_resp <- ggplot(clin,aes(x=Response,y=rna_bcb)) + geom_boxplot(aes(fill=Response),outlier.shape=NA) + geom_jitter() + theme_classic() + ylab(expression(BCB[RNA])) +
xlab("") + scale_fill_manual(values=c("red","blue")) + annotate("text",x=1.5,y=2000,label=paste0("Wilcox p = ",formatSF(tmpwt$p.value,digits=3,format="e"))) + scale_y_log10() +
ggtitle(paste0("NR=",nrow(clin[clin$Response=="nonresponder",])," R=",nrow(clin[clin$Response=="responder",])))
ggsave("figure12_outputs/2b.2.rna_bcb_response.png",g_rna_bcb_resp,height=7,width=7)
knitr::include_graphics("figure12_outputs/2b.2.rna_bcb_response.png")
```

<img src="figure12_outputs/2b.2.rna_bcb_response.png" width="75%" height="75%" />

``` {.r}
ct_r <- formatSF(cor.test(log10(clin$rna_tcb),clin$mean_CD2_CD3D_CD3E_CD3G,method="spearman")$estimate,format="f",digits=3)
ct_p <- formatSF(cor.test(log10(clin$rna_tcb),clin$mean_CD2_CD3D_CD3E_CD3G,method="spearman")$p.value,format="e",digits=3)
g_rna_tcb_CD2CD3 <- ggplot(clin,aes(x=mean_CD2_CD3D_CD3E_CD3G,y=rna_tcb)) + geom_point(size=3) + geom_smooth(method="lm",se=FALSE) + theme_classic() + scale_y_log10() +
xlab("mean CD2 CD3D CD3E CD3G") + ylab(expression(TCB[RNA])) + annotate(geom="text",x=0.5,y=10,label=paste0("rho = ",ct_r,"\np = ",ct_p,"\nn = ",nrow(clin)))
ggsave("figure12_outputs/2c.rna_CD2CD3_rna_tcb.png",g_rna_tcb_CD2CD3,height=7,width=7)
knitr::include_graphics("figure12_outputs/2c.rna_CD2CD3_rna_tcb.png")
```

<img src="figure12_outputs/2c.rna_CD2CD3_rna_tcb.png" width="75%" height="75%" />

``` {.r}
ct_r <- formatSF(cor.test(log10(clin$rna_bcb),clin$mean_IGJ_MZB1,method="spearman")$estimate,format="f",digits=3)
```

    ## Warning in cor.test.default(log10(clin$rna_bcb), clin$mean_IGJ_MZB1, method
    ## = "spearman"): Cannot compute exact p-value with ties

``` {.r}
ct_p <- formatSF(cor.test(log10(clin$rna_bcb),clin$mean_IGJ_MZB1,method="spearman")$p.value,format="e",digits=3)
```

    ## Warning in cor.test.default(log10(clin$rna_bcb), clin$mean_IGJ_MZB1, method
    ## = "spearman"): Cannot compute exact p-value with ties

``` {.r}
g_rna_bcb_IGJMZB1 <- ggplot(clin,aes(x=mean_IGJ_MZB1,y=rna_bcb)) + geom_point(size=3) + geom_smooth(method="lm",se=FALSE) + theme_classic() + scale_y_log10() +
xlab("mean IGJ MZB1") + ylab(expression(BCB[RNA])) + annotate(geom="text",x=0.5,y=1000,label=paste0("rho = ",ct_r,"\np = ",ct_p,"\nn = ",nrow(clin)))
ggsave("figure12_outputs/2d.rna_IGJMZB1_rna_tcb.png",g_rna_bcb_IGJMZB1,height=7,width=7)
knitr::include_graphics("figure12_outputs/2d.rna_IGJMZB1_rna_tcb.png")
```

<img src="figure12_outputs/2d.rna_IGJMZB1_rna_tcb.png" width="75%" height="75%" />

``` {.r}
clin$plot_tcb_bcb_binary <- paste0(clin$plot_tcb_binary,", ",clin$plot_bcb_binary)
clin$plot_tcb_bcb_binary_bestrest <- ifelse(clin$plot_tcb_bcb_binary=="High TCB, High BCB",clin$plot_tcb_bcb_binary,"Others")
makeKMplot_2group(clin[!is.na(clin$surv),],"plot_tcb_bcb_binary_bestrest","figure12_outputs/2e.rna_tcbbcb_survival.png")
```

    ## png 
    ##   2

``` {.r}
knitr::include_graphics("figure12_outputs/2e.rna_tcbbcb_survival.png")
```

<img src="figure12_outputs/2e.rna_tcbbcb_survival.png" width="75%" height="75%" />

``` {.r}
ct_r <- formatSF(cor.test(log10(clin[!is.na(clin$same_location)&clin$same_location==TRUE,]$rna_tcb),log10(clin[!is.na(clin$same_location)&clin$same_location==TRUE,]$dna_tcb))$estimate,
format="f",digits=3)
ct_p <- formatSF(cor.test(log10(clin[!is.na(clin$same_location)&clin$same_location==TRUE,]$rna_tcb),log10(clin[!is.na(clin$same_location)&clin$same_location==TRUE,]$dna_tcb))$p.value,
format="e",digits=3)
g_rna_tcb_dna_tcb <- ggplot(clin[!is.na(clin$same_location)&clin$same_location==TRUE,],aes(x=rna_tcb,y=dna_tcb)) + geom_point(size=3) + geom_smooth(method="lm",se=FALSE) +
theme_classic() + annotate(geom="text",x=0.02,y=0.2,label=paste0("r = ",ct_r,"\np = ",ct_p,"\nn = ",nrow(clin[!is.na(clin$same_location)&clin$same_location==TRUE,]))) +
xlab(expression(TCB[RNA])) + ylab(expression(TCB[DNA])) + scale_x_log10() + scale_y_log10()
ggsave("figure12_outputs/2f.rna_tcb_dna_tcb.png",g_rna_tcb_dna_tcb,height=7,width=7)
knitr::include_graphics("figure12_outputs/2f.rna_tcb_dna_tcb.png")
```

<img src="figure12_outputs/2f.rna_tcb_dna_tcb.png" width="75%" height="75%" />

``` {.r}
ct_r <- formatSF(cor.test(log10(clin[!is.na(clin$same_location)&clin$same_location==TRUE,]$rna_bcb),log10(clin[!is.na(clin$same_location)&clin$same_location==TRUE,]$dna_bcb))$estimate,
format="f",digits=3)
ct_p <- formatSF(cor.test(log10(clin[!is.na(clin$same_location)&clin$same_location==TRUE,]$rna_bcb),log10(clin[!is.na(clin$same_location)&clin$same_location==TRUE,]$dna_bcb))$p.value,
format="e",digits=3)
g_rna_bcb_dna_bcb <- ggplot(clin[!is.na(clin$same_location)&clin$same_location==TRUE,],aes(x=rna_bcb,y=dna_bcb)) + geom_point(size=3) + geom_smooth(method="lm",se=FALSE) +
theme_classic() + annotate(geom="text",x=0.5,y=0.2,label=paste0("r = ",ct_r,"\np = ",ct_p,"\nn = ",nrow(clin[!is.na(clin$same_location)&clin$same_location==TRUE,]))) +
xlab(expression(BCB[RNA])) + ylab(expression(BCB[DNA])) + scale_x_log10() + scale_y_log10()
ggsave("figure12_outputs/2g.rna_bcb_dna_bcb.png",g_rna_bcb_dna_bcb,height=7,width=7)
knitr::include_graphics("figure12_outputs/2g.rna_bcb_dna_bcb.png")
```

<img src="figure12_outputs/2g.rna_bcb_dna_bcb.png" width="75%" height="75%" />

``` {.r}
clin$tcr_overlap_binary <- NA
clin[!is.na(clin$same_location)&clin$same_location==TRUE,]$tcr_overlap_binary <- ifelse(clin[!is.na(clin$same_location)&clin$same_location==TRUE,]$tcr_dna_overlap_reads==0,
"No overlap","TCR overlap")
tdf <- melt(table(clin[!is.na(clin$same_location)&clin$same_location==TRUE,]$plot_tcb_binary,clin[!is.na(clin$same_location)&clin$same_location==TRUE,]$tcr_overlap_binary))
g_rnr_temp <- ggplot(tdf) + geom_mosaic(aes(weight=value,x=product(Var2),fill=Var1),offset=0.01) + scale_fill_manual(values=c("red","blue")) + theme_bw() +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.title=element_blank(),axis.text.x=element_text(size=12),axis.text.y=element_text(size=12)) +
xlab("") + ylab("") + scale_y_continuous(breaks=seq(0,1,by=0.25))
ggsave("figure12_outputs/2h.1.dna_rna_tcr_overlap.png",g_rnr_temp + geom_text(data = ggplot_build(g_rnr_temp)$data[[1]], aes(x = (xmin+xmax)/2, y = (ymin+ymax)/2, label=.wt)))
knitr::include_graphics("figure12_outputs/2h.1.dna_rna_tcr_overlap.png")
```

<img src="figure12_outputs/2h.1.dna_rna_tcr_overlap.png" width="75%" height="75%" />

``` {.r}
clin$ig_overlap_binary <- NA
clin[!is.na(clin$same_location)&clin$same_location==TRUE,]$ig_overlap_binary <- ifelse(clin[!is.na(clin$same_location)&clin$same_location==TRUE,]$ig_dna_overlap_reads==0,
"No overlap","Ig overlap")
tdf <- melt(table(clin[!is.na(clin$same_location)&clin$same_location==TRUE,]$plot_bcb_binary,clin[!is.na(clin$same_location)&clin$same_location==TRUE,]$ig_overlap_binary))
g_rnr_temp <- ggplot(tdf) + geom_mosaic(aes(weight=value,x=product(Var2),fill=Var1),offset=0.01) + scale_fill_manual(values=c("red","blue")) + theme_bw() +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.title=element_blank(),axis.text.x=element_text(size=12),axis.text.y=element_text(size=12)) +
xlab("") + ylab("") + scale_y_continuous(breaks=seq(0,1,by=0.25))
ggsave("figure12_outputs/2h.2.dna_rna_ig_overlap.png",g_rnr_temp + geom_text(data = ggplot_build(g_rnr_temp)$data[[1]], aes(x = (xmin+xmax)/2, y = (ymin+ymax)/2, label=.wt)))
knitr::include_graphics("figure12_outputs/2h.2.dna_rna_ig_overlap.png")
```

<img src="figure12_outputs/2h.2.dna_rna_ig_overlap.png" width="75%" height="75%" />

``` {.r}
dna_clin$plot_cohort_dna_tcb_binary <- NA
dna_clin[dna_clin$cohort=="MGH",]$plot_cohort_dna_tcb_binary <- ifelse(dna_clin[dna_clin$cohort=="MGH",]$dna_tcb>median(dna_clin[dna_clin$cohort=="MGH",]$dna_tcb),"High TCB","Low TCB")
dna_clin[dna_clin$cohort=="VanAllen",]$plot_cohort_dna_tcb_binary <- ifelse(dna_clin[dna_clin$cohort=="VanAllen",]$dna_tcb>median(dna_clin[dna_clin$cohort=="VanAllen",]$dna_tcb),
"High TCB","Low TCB")
dna_clin[dna_clin$cohort=="Wargo",]$plot_cohort_dna_tcb_binary <- ifelse(dna_clin[dna_clin$cohort=="Wargo",]$dna_tcb>median(dna_clin[dna_clin$cohort=="Wargo",]$dna_tcb),
"High TCB","Low TCB")
dna_clin[dna_clin$not_cell_line&dna_clin$cohort=="Ribas",]$plot_cohort_dna_tcb_binary <-
ifelse(dna_clin[dna_clin$not_cell_line&dna_clin$cohort=="Ribas",]$dna_tcb>median(dna_clin[dna_clin$not_cell_line&dna_clin$cohort=="Ribas",]$dna_tcb),"High TCB","Low TCB")

dna_clin$plot_cohort_dna_bcb_binary <- NA
dna_clin[dna_clin$cohort=="MGH",]$plot_cohort_dna_bcb_binary <- ifelse(dna_clin[dna_clin$cohort=="MGH",]$dna_bcb>median(dna_clin[dna_clin$cohort=="MGH",]$dna_bcb),"High BCB","Low BCB")
dna_clin[dna_clin$cohort=="VanAllen",]$plot_cohort_dna_bcb_binary <- ifelse(dna_clin[dna_clin$cohort=="VanAllen",]$dna_bcb>median(dna_clin[dna_clin$cohort=="VanAllen",]$dna_bcb),
"High BCB","Low BCB")
dna_clin[dna_clin$cohort=="Wargo",]$plot_cohort_dna_bcb_binary <- ifelse(dna_clin[dna_clin$cohort=="Wargo",]$dna_bcb>median(dna_clin[dna_clin$cohort=="Wargo",]$dna_bcb),
"High BCB","Low BCB")
dna_clin[dna_clin$not_cell_line&dna_clin$cohort=="Ribas",]$plot_cohort_dna_bcb_binary <-
ifelse(dna_clin[dna_clin$not_cell_line&dna_clin$cohort=="Ribas",]$dna_bcb>median(dna_clin[dna_clin$not_cell_line&dna_clin$cohort=="Ribas",]$dna_bcb),"High BCB","Low BCB")

makeKMplot_2group(dna_clin[dna_clin$not_cell_line,],"plot_cohort_dna_tcb_binary","figure12_outputs/2i.dna_tcb_survival.png")
```

    ## png 
    ##   2

``` {.r}
knitr::include_graphics("figure12_outputs/2i.dna_tcb_survival.png")
```

<img src="figure12_outputs/2i.dna_tcb_survival.png" width="75%" height="75%" />

``` {.r}
makeKMplot_2group(dna_clin[dna_clin$not_cell_line,],"plot_cohort_dna_bcb_binary","figure12_outputs/2j.dna_bcb_survival.png")
```

    ## png 
    ##   2

``` {.r}
knitr::include_graphics("figure12_outputs/2j.dna_bcb_survival.png")
```

<img src="figure12_outputs/2j.dna_bcb_survival.png" width="75%" height="75%" />

``` {.r}
dna_clin$plot_combined_cohort_dna_tcb_cohort_dna_bcb <- NA
dna_clin[dna_clin$not_cell_line,]$plot_combined_cohort_dna_tcb_cohort_dna_bcb <- paste0(dna_clin[dna_clin$not_cell_line,]$plot_cohort_dna_tcb_binary,", ",
dna_clin[dna_clin$not_cell_line,]$plot_cohort_dna_bcb_binary)
dna_clin$plot_combined_cohort_dna_tcb_cohort_dna_bcb_bestrest <- NA
dna_clin[dna_clin$not_cell_line,]$plot_combined_cohort_dna_tcb_cohort_dna_bcb_bestrest <-
ifelse(dna_clin[dna_clin$not_cell_line,]$plot_combined_cohort_dna_tcb_cohort_dna_bcb=="High TCB, High BCB","High TCB, High BCB","Others")

makeKMplot_2group(dna_clin[dna_clin$not_cell_line,],"plot_combined_cohort_dna_tcb_cohort_dna_bcb_bestrest","figure12_outputs/2k.dna_tcbbcb_survival.png")
```

    ## png 
    ##   2

``` {.r}
knitr::include_graphics("figure12_outputs/2k.dna_tcbbcb_survival.png")
```

<img src="figure12_outputs/2k.dna_tcbbcb_survival.png" width="75%" height="75%" />

``` {.r}
dna_clin$plot_combined_tmb_cohort_dna_tcb <- NA
dna_clin[dna_clin$not_cell_line,]$plot_combined_tmb_cohort_dna_tcb <- paste0(dna_clin[dna_clin$not_cell_line,]$plot_tmb_binary,", ",
dna_clin[dna_clin$not_cell_line,]$plot_cohort_dna_tcb_binary)
dna_clin$plot_combined_tmb_cohort_dna_tcb_bestrest <- NA
dna_clin[dna_clin$not_cell_line,]$plot_combined_tmb_cohort_dna_tcb_bestrest <- ifelse(dna_clin[dna_clin$not_cell_line,]$plot_combined_tmb_cohort_dna_tcb=="High TMB, High TCB",
"High TMB, High TCB","Others")

makeKMplot_2group(dna_clin[dna_clin$not_cell_line,],"plot_combined_tmb_cohort_dna_tcb_bestrest","figure12_outputs/2l.dna_tmbtcb_survival_bestrest.png")
```

    ## png 
    ##   2

``` {.r}
knitr::include_graphics("figure12_outputs/2l.dna_tmbtcb_survival_bestrest.png")
```

<img src="figure12_outputs/2l.dna_tmbtcb_survival_bestrest.png" width="75%" height="75%" />

``` {.r}
makeKMplot_4group(dna_clin[dna_clin$not_cell_line,],"plot_combined_tmb_cohort_dna_tcb","figure12_outputs/2m.dna_tmbtcb_survival.png")
```

    ## png 
    ##   2

``` {.r}
knitr::include_graphics("figure12_outputs/2m.dna_tmbtcb_survival.png")
```

<img src="figure12_outputs/2m.dna_tmbtcb_survival.png" width="75%" height="75%" />

``` {.r}
tmp_ft <- fisher.test(table(dna_clin[dna_clin$not_cell_line,]$Response,factor(dna_clin[dna_clin$not_cell_line,]$plot_combined_tmb_cohort_dna_tcb_bestrest,
levels=c("Others","High TMB, High TCB"))))
tmp_ft_pstring <- paste0("Fisher p=",formatSF(tmp_ft$p.value,digits=3)," OR=",formatSF(tmp_ft$estimate,digits=3))
tdf <- melt(table(dna_clin[dna_clin$not_cell_line,]$Response,factor(dna_clin[dna_clin$not_cell_line,]$plot_combined_tmb_cohort_dna_tcb_bestrest,levels=c("High TMB, High TCB","Others"))))
g_rnr_temp <- ggplot(tdf) + geom_mosaic(aes(weight=value,x=product(Var2),fill=Var1),offset=0.01) + scale_fill_manual(values=c("red","blue")) + theme_bw() +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.title=element_blank(),axis.text.x=element_text(size=12),axis.text.y=element_text(size=12)) +
xlab("") + ylab("") + annotate(geom="text",x=0.5,y=1.02,label=tmp_ft_pstring) + scale_y_continuous(breaks=seq(0,1,by=0.25))
ggsave("figure12_outputs/2n.dna_tmbtcb_bestrest_response.png",g_rnr_temp + geom_text(data = ggplot_build(g_rnr_temp)$data[[1]], aes(x = (xmin+xmax)/2, y = (ymin+ymax)/2, label=.wt)),height=7,width=7)
knitr::include_graphics("figure12_outputs/2n.dna_tmbtcb_bestrest_response.png")
```

<img src="figure12_outputs/2n.dna_tmbtcb_bestrest_response.png" width="75%" height="75%" />

Figure S8

``` {.r}
dna_clin$plot_combined_tmb_cohort_dna_bcb <- NA
dna_clin[dna_clin$not_cell_line,]$plot_combined_tmb_cohort_dna_bcb <- paste0(dna_clin[dna_clin$not_cell_line,]$plot_tmb_binary,", ",
dna_clin[dna_clin$not_cell_line,]$plot_cohort_dna_bcb_binary)
dna_clin$plot_combined_tmb_cohort_dna_bcb_bestrest <- NA
dna_clin[dna_clin$not_cell_line,]$plot_combined_tmb_cohort_dna_bcb_bestrest <- ifelse(dna_clin[dna_clin$not_cell_line,]$plot_combined_tmb_cohort_dna_bcb=="High TMB, High BCB",
"High TMB, High BCB","Others")

tmp_kw <- kruskal.test(log10(clin$rna_tcb),factor(clin$cohort))
tmp_kw_pstring <- paste0("Kruskal-Wallis p=",formatSF(tmp_kw$p.value,digits=3))
g_rna_tcb_cohort <- ggplot(clin,aes(x=cohort,y=rna_tcb)) + geom_boxplot(aes(fill=cohort)) + geom_jitter() + theme_classic() + ylab(expression(TCB[RNA])) + xlab ("") +
scale_y_log10() + annotate(geom="text",x=2.5,y=10,label=tmp_kw_pstring) + theme(legend.position="none")
ggsave("figure12_outputs/S8a.1.rna_tcb_cohort.png",g_rna_tcb_cohort)
knitr::include_graphics("figure12_outputs/S8a.1.rna_tcb_cohort.png")
```

<img src="figure12_outputs/S8a.1.rna_tcb_cohort.png" width="75%" height="75%" />

``` {.r}
tmp_kw <- kruskal.test(log10(clin$rna_bcb),factor(clin$cohort))
tmp_kw_pstring <- paste0("Kruskal-Wallis p=",formatSF(tmp_kw$p.value,digits=3))
g_rna_bcb_cohort <- ggplot(clin,aes(x=cohort,y=rna_bcb)) + geom_boxplot(aes(fill=cohort)) + geom_jitter() + theme_classic() + ylab(expression(BCB[RNA])) + xlab ("") + scale_y_log10() +
annotate(geom="text",x=2.5,y=5000,label=tmp_kw_pstring) + theme(legend.position="none")
ggsave("figure12_outputs/S8a.2.rna_bcb_cohort.png",g_rna_bcb_cohort)
knitr::include_graphics("figure12_outputs/S8a.2.rna_bcb_cohort.png")
```

<img src="figure12_outputs/S8a.2.rna_bcb_cohort.png" width="75%" height="75%" />

``` {.r}
ct_r <- formatSF(cor.test(log10(clin$rna_tcb),log10(clin$rna_bcb),method="spearman")$estimate,format="f",digits=3)
ct_p <- formatSF(cor.test(log10(clin$rna_tcb),log10(clin$rna_bcb),method="spearman")$p.value,format="e",digits=3)
g_rna_tcb_rna_bcb <- ggplot(clin,aes(x=rna_tcb,y=rna_bcb)) + geom_point(size=3) + geom_smooth(method="lm",se=FALSE) + theme_classic() + scale_x_log10() + xlab(expression(TCB[RNA])) +
scale_y_log10() + ylab(expression(BCB[RNA])) + ggtitle(paste0("rho = ",ct_r," p = ",ct_p," n=",nrow(clin)))
ggsave("figure12_outputs/S8b.rna_tcb_rna_bcb.png",g_rna_tcb_rna_bcb,height=7,width=7)
knitr::include_graphics("figure12_outputs/S8b.rna_tcb_rna_bcb.png")
```

<img src="figure12_outputs/S8b.rna_tcb_rna_bcb.png" width="75%" height="75%" />

``` {.r}
makeKMplot_4group(clin[!is.na(clin$surv),],"plot_tcb_bcb_binary","figure12_outputs/S8c.rna_tcbbcb_survival.png")
```

    ## png 
    ##   2

``` {.r}
knitr::include_graphics("figure12_outputs/S8c.rna_tcbbcb_survival.png")
```

<img src="figure12_outputs/S8c.rna_tcbbcb_survival.png" width="75%" height="75%" />

``` {.r}
makeResponseMosaic(clin,"plot_tcb_bcb_binary_bestrest","figure12_outputs/S8d.dna_tcbbcb_bestrest_response.png")
knitr::include_graphics("figure12_outputs/S8d.dna_tcbbcb_bestrest_response.png")
```

<img src="figure12_outputs/S8d.dna_tcbbcb_bestrest_response.png" width="75%" height="75%" />

``` {.r}
makeResponseMosaic_rotatelabels(clin,"plot_tcb_bcb_binary","figure12_outputs/S8e.dna_tcbbcb_response.png")
knitr::include_graphics("figure12_outputs/S8e.dna_tcbbcb_response.png")
```

<img src="figure12_outputs/S8e.dna_tcbbcb_response.png" width="75%" height="75%" />

``` {.r}
ct_r <- formatSF(cor.test(skcm_clin$GEP,log10(skcm_clin$rna_tcb),method="spearman")$estimate,format="f",digits=3)
ct_p <- formatSF(cor.test(skcm_clin$GEP,log10(skcm_clin$rna_tcb),method="spearman")$p.value,format="e",digits=3)
g_tcga_gep_tcb <- ggplot(skcm_clin,aes(x=GEP,y=rna_tcb)) + geom_point(size=3) + geom_smooth(method="lm",se=FALSE) + theme_classic() + scale_y_log10() +
xlab("GEP") + ylab(expression(TCB[RNA])) + annotate(geom="text",x=2,y=10,label=paste0("rho = ",ct_r,"\np = ",ct_p,"\nn=",nrow(skcm_clin)))
ggsave("figure12_outputs/S8f.tcga_gep_tcb.png",g_tcga_gep_tcb,height=7,width=7)
knitr::include_graphics("figure12_outputs/S8f.tcga_gep_tcb.png")
```

<img src="figure12_outputs/S8f.tcga_gep_tcb.png" width="75%" height="75%" />

``` {.r}
ct_r <- formatSF(cor.test(clin$GEP,log10(clin$rna_tcb),method="spearman")$estimate,format="f",digits=3)
ct_p <- formatSF(cor.test(clin$GEP,log10(clin$rna_tcb),method="spearman")$p.value,format="e",digits=3)
g_rna_gep_tcb <- ggplot(clin,aes(x=GEP,y=rna_tcb)) + geom_point(size=3) + geom_smooth(method="lm",se=FALSE) + theme_classic() + scale_y_log10() +
xlab("GEP") + ylab(expression(TCB[RNA])) + annotate(geom="text",x=2,y=10,label=paste0("rho = ",ct_r,"\np = ",ct_p,"\nn=",nrow(clin)))
ggsave("figure12_outputs/S8g.rna_gep_tcb.png",g_rna_gep_tcb,height=7,width=7)
knitr::include_graphics("figure12_outputs/S8g.rna_gep_tcb.png")
```

<img src="figure12_outputs/S8g.rna_gep_tcb.png" width="75%" height="75%" />

Figure S9

``` {.r}
tmp_kw <- kruskal.test(log10(dna_clin[dna_clin$not_cell_line,]$dna_tcb),factor(dna_clin[dna_clin$not_cell_line,]$cohort))
tmp_kw_pstring <- paste0("Kruskal-Wallis p=",formatSF(tmp_kw$p.value,format="e",digits=3))
g_dna_tcb_cohort <- ggplot(dna_clin[dna_clin$not_cell_line,],aes(x=cohort,y=dna_tcb)) + geom_boxplot(aes(fill=cohort)) + geom_jitter() + theme_classic() + ylab(expression(TCB[DNA])) +
xlab ("") + scale_y_log10() + annotate(geom="text",x=2.5,y=1,label=tmp_kw_pstring) + theme(legend.position="none")
ggsave("figure12_outputs/S9a.1.dna_tcb_cohort.png",g_dna_tcb_cohort)
knitr::include_graphics("figure12_outputs/S9a.1.dna_tcb_cohort.png")
```

<img src="figure12_outputs/S9a.1.dna_tcb_cohort.png" width="75%" height="75%" />

``` {.r}
tmp_kw <- kruskal.test(log10(dna_clin[dna_clin$not_cell_line,]$dna_bcb),factor(dna_clin[dna_clin$not_cell_line,]$cohort))
tmp_kw_pstring <- paste0("Kruskal-Wallis p=",formatSF(tmp_kw$p.value,digits=3))
g_dna_bcb_cohort <- ggplot(dna_clin[dna_clin$not_cell_line,],aes(x=cohort,y=dna_bcb)) + geom_boxplot(aes(fill=cohort)) + geom_jitter() + theme_classic() + ylab(expression(BCB[DNA])) +
xlab ("") + scale_y_log10() + annotate(geom="text",x=2.5,y=0.25,label=tmp_kw_pstring) + theme(legend.position="none")
ggsave("figure12_outputs/S9a.2.dna_bcb_cohort.png",g_dna_bcb_cohort)
knitr::include_graphics("figure12_outputs/S9a.2.dna_bcb_cohort.png")
```

<img src="figure12_outputs/S9a.2.dna_bcb_cohort.png" width="75%" height="75%" />

``` {.r}
ct_r <- formatSF(cor.test(log10(dna_clin[dna_clin$not_cell_line,]$dna_tcb),log10(dna_clin[dna_clin$not_cell_line,]$dna_bcb),method="spearman")$estimate,format="f",digits=3)
ct_p <- formatSF(cor.test(log10(dna_clin[dna_clin$not_cell_line,]$dna_tcb),log10(dna_clin[dna_clin$not_cell_line,]$dna_bcb),method="spearman")$p.value,format="e",digits=3)
g_dna_tcb_dna_bcb <- ggplot(dna_clin[dna_clin$not_cell_line,],aes(x=dna_tcb,y=dna_bcb)) + geom_point(size=3) + geom_smooth(method="lm",se=FALSE) + theme_classic() + scale_y_log10() +
scale_x_log10() + xlab(expression(TCB[DNA])) + ylab(expression(BCB[DNA])) +
annotate(geom="text",x=0.01,y=0.75,label=paste0("rho = ",ct_r,"\np = ",ct_p,"\nn=",nrow(dna_clin[dna_clin$not_cell_line,])))
ggsave("figure12_outputs/S9b.dna_tcb_dna_bcb.png",g_dna_tcb_dna_bcb,height=7,width=7)
knitr::include_graphics("figure12_outputs/S9b.dna_tcb_dna_bcb.png")
```

<img src="figure12_outputs/S9b.dna_tcb_dna_bcb.png" width="75%" height="75%" />

``` {.r}
makeResponseMosaic(dna_clin[dna_clin$not_cell_line,],"plot_cohort_dna_tcb_binary","figure12_outputs/S9c.dna_tcb_response.png")
knitr::include_graphics("figure12_outputs/S9c.dna_tcb_response.png")
```

<img src="figure12_outputs/S9c.dna_tcb_response.png" width="75%" height="75%" />

``` {.r}
makeResponseMosaic(dna_clin[dna_clin$not_cell_line,],"plot_cohort_dna_bcb_binary","figure12_outputs/S9d.dna_bcb_response.png")
knitr::include_graphics("figure12_outputs/S9d.dna_bcb_response.png")
```

<img src="figure12_outputs/S9d.dna_bcb_response.png" width="75%" height="75%" />

``` {.r}
makeResponseMosaic(dna_clin[dna_clin$not_cell_line,],"plot_combined_cohort_dna_tcb_cohort_dna_bcb_bestrest","figure12_outputs/S9e.dna_tcbbcb_bestrest_response.png")
knitr::include_graphics("figure12_outputs/S9e.dna_tcbbcb_bestrest_response.png")
```

<img src="figure12_outputs/S9e.dna_tcbbcb_bestrest_response.png" width="75%" height="75%" />

``` {.r}
makeKMplot_4group(dna_clin[dna_clin$not_cell_line,],"plot_combined_cohort_dna_tcb_cohort_dna_bcb","figure12_outputs/S9f.dna_tcbbcb_survival.png")
```

    ## png 
    ##   2

``` {.r}
knitr::include_graphics("figure12_outputs/S9f.dna_tcbbcb_survival.png")
```

<img src="figure12_outputs/S9f.dna_tcbbcb_survival.png" width="75%" height="75%" />

``` {.r}
makeResponseMosaic_rotatelabels(dna_clin[dna_clin$not_cell_line,],"plot_combined_cohort_dna_tcb_cohort_dna_bcb","figure12_outputs/S9g.dna_tcbbcb_response.png")
knitr::include_graphics("figure12_outputs/S9g.dna_tcbbcb_response.png")
```

<img src="figure12_outputs/S9g.dna_tcbbcb_response.png" width="75%" height="75%" />

``` {.r}
ct_r <- formatSF(cor.test(skcm_clin$GEP,log10(skcm_clin$dna_tcb),method="spearman")$estimate,format="f",digits=3)
ct_p <- formatSF(cor.test(skcm_clin$GEP,log10(skcm_clin$dna_tcb),method="spearman")$p.value,format="e",digits=3)
g_tcga_gep_dna_tcb <- ggplot(skcm_clin,aes(x=GEP,y=dna_tcb)) + geom_point(size=3) + geom_smooth(method="lm",se=FALSE) + theme_classic() + scale_y_log10() +
xlab("GEP") + ylab(expression(TCB[DNA])) + annotate(geom="text",x=2,y=0.2,label=paste0("rho = ",ct_r,"\np = ",ct_p,"\nn=",nrow(skcm_clin)))
ggsave("figure12_outputs/S9h.tcga_gep_dna_tcb.png",g_tcga_gep_dna_tcb,height=7,width=7)
knitr::include_graphics("figure12_outputs/S9h.tcga_gep_dna_tcb.png")
```

<img src="figure12_outputs/S9h.tcga_gep_dna_tcb.png" width="75%" height="75%" />

``` {.r}
ct_r <- formatSF(cor.test(clin[!is.na(clin$same_location)&clin$same_location==TRUE,]$GEP,log10(clin[!is.na(clin$same_location)&clin$same_location==TRUE,]$dna_tcb),method="spearman")$estimate,
format="f",digits=3)
ct_p <- formatSF(cor.test(clin[!is.na(clin$same_location)&clin$same_location==TRUE,]$GEP,log10(clin[!is.na(clin$same_location)&clin$same_location==TRUE,]$dna_tcb),method="spearman")$p.value,
format="e",digits=3)
g_gep_dna_tcb <- ggplot(clin[!is.na(clin$same_location)&clin$same_location==TRUE,],aes(x=GEP,y=dna_tcb)) + geom_point(size=3) + geom_smooth(method="lm",se=FALSE) + theme_classic() +
scale_y_log10() + xlab("GEP") + ylab(expression(TCB[RNA])) +
annotate(geom="text",x=2,y=0.2,label=paste0("rho = ",ct_r,"\np = ",ct_p,"\nn=",nrow(clin[!is.na(clin$same_location)&clin$same_location==TRUE,])))
ggsave("figure12_outputs/S9i.rna_gep_dna_tcb.png",g_gep_dna_tcb,height=7,width=7)
knitr::include_graphics("figure12_outputs/S9i.rna_gep_dna_tcb.png")
```

<img src="figure12_outputs/S9i.rna_gep_dna_tcb.png" width="75%" height="75%" />

``` {.r}
tmpwt <- wilcox.test(dna_clin[dna_clin$not_cell_line,]$dna_tcb,dna_clin[dna_clin$not_cell_line,]$dna_bcb)
m_dna_tcbbcb <- melt(dna_clin[dna_clin$not_cell_line,c("pair_id","dna_tcb","dna_bcb")])
g_dnatcb_vs_dnabcb <- ggplot(m_dna_tcbbcb,aes(x=variable,y=value)) + geom_boxplot(aes(fill=variable),outlier.shape=NA) + geom_jitter() + theme_classic() + ylab("Normalized DNA Read Count") +
xlab("") + scale_fill_manual(values=c("#0072B2","#D55E00")) + annotate("text",x=1.5,y=1,label=paste0("Wilcox p = ",formatSF(tmpwt$p.value,digits=3,format="e"))) + scale_y_log10() +
theme(legend.position="none")
ggsave("figure12_outputs/S9j.dna_tcb_vs_dna_bcb.png",g_dnatcb_vs_dnabcb,width=4,height=7)
knitr::include_graphics("figure12_outputs/S9j.dna_tcb_vs_dna_bcb.png")
```

<img src="figure12_outputs/S9j.dna_tcb_vs_dna_bcb.png" width="75%" height="75%" />

``` {.r}
tmpwt <- wilcox.test(clin$rna_tcb,clin$rna_bcb)
m_rna_tcbbcb <- melt(clin[,c("rna_id","rna_tcb","rna_bcb")])
g_rnatcb_vs_rnabcb <- ggplot(m_rna_tcbbcb,aes(x=variable,y=value)) + geom_boxplot(aes(fill=variable),outlier.shape=NA) + geom_jitter() + theme_classic() + ylab("Normalized RNA Read Count") +
xlab("") + scale_fill_manual(values=c("#0072B2","#D55E00")) + annotate("text",x=1.5,y=5000,label=paste0("Wilcox p = ",formatSF(tmpwt$p.value,digits=3,format="e"))) + scale_y_log10() +
theme(legend.position="none")
ggsave("figure12_outputs/S9k.rna_tcb_vs_rna_bcb.png",g_rnatcb_vs_rnabcb,width=4,height=7)
knitr::include_graphics("figure12_outputs/S9k.rna_tcb_vs_rna_bcb.png")
```

<img src="figure12_outputs/S9k.rna_tcb_vs_rna_bcb.png" width="75%" height="75%" />

``` {.r}
clin$T.cells.bmode <- clin$T.cells.CD8.bmode+clin$T.cells.CD4.bmode
clin$T.cells.bmode.norm <- clin$T.cells.bmode/clin$Absolute.score..sig.score..bmode
clin$B.cells.bmode.norm <- clin$B.cells.bmode/clin$Absolute.score..sig.score..bmode

m_cibersort_all_fractions <- melt(clin[,c("rna_id","T.cells.bmode.norm","B.cells.bmode.norm")])

tmpwt <- wilcox.test(clin$T.cells.bmode.norm,clin$B.cells.bmode.norm)
g_cibersort_tcell_bcell <- ggplot(m_cibersort_all_fractions[m_cibersort_all_fractions$variable %in% c("T.cells.bmode.norm","B.cells.bmode.norm"),],aes(x=variable,y=value)) +
geom_boxplot(outlier.shape=NA,aes(fill=variable)) + geom_jitter() + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="none") +
scale_fill_manual(values=c("#0072B2","#D55E00")) + ylab("Cell type fraction") + xlab("") + annotate("text",x=1.5,y=0.25,label=paste0("Wilcox p = ",formatC(tmpwt$p.value,digits=3,format="e")))
ggsave("figure12_outputs/S9l.cibersort_tcell_bcell_boxplot.png",g_cibersort_tcell_bcell,height=7,width=4)
knitr::include_graphics("figure12_outputs/S9l.cibersort_tcell_bcell_boxplot.png")
```

<img src="figure12_outputs/S9l.cibersort_tcell_bcell_boxplot.png" width="75%" height="75%" />

``` {.r}
tmp_ct <- cor.test(clin$sc_naiveb,clin$rna_bcb,method="spearman")
tmp_ct_rho <- formatC(as.numeric(tmp_ct$estimate),digits=2)
tmp_ct_p <- formatC(as.numeric(tmp_ct$p.value),format="e",digits=2)
tmp_ct_string <- paste0("rho = ",tmp_ct_rho,"\n"," p = ",tmp_ct_p,"\n n = ",nrow(clin))
g_bcb_rna_vs_scnaiveb <- ggplot(clin,aes(x=sc_naiveb,y=rna_bcb)) + geom_point() + geom_smooth(method="lm",se=FALSE) +
theme_classic(base_size=16) + theme(axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black")) +
xlab("naive B cell scRNA signature") + ylab("BCB_RNA")+ scale_y_log10() + annotate(geom="text",x=0.01,y=0.1,label=tmp_ct_string)
ggsave("figure12_outputs/S9m.naiveb_vs_rna_bcb.png",g_bcb_rna_vs_scnaiveb,height=7,width=7)
knitr::include_graphics("figure12_outputs/S9m.naiveb_vs_rna_bcb.png")
```

<img src="figure12_outputs/S9m.naiveb_vs_rna_bcb.png" width="75%" height="75%" />

``` {.r}
tmp_ct <- cor.test(clin$sc_plasmab,clin$rna_bcb,method="spearman")
tmp_ct_rho <- formatC(as.numeric(tmp_ct$estimate),digits=2)
tmp_ct_p <- formatC(as.numeric(tmp_ct$p.value),format="e",digits=2)
tmp_ct_string <- paste0("rho = ",tmp_ct_rho,"\n"," p = ",tmp_ct_p,"\n n = ",nrow(clin))
g_bcb_rna_vs_scplasmab <- ggplot(clin,aes(x=sc_plasmab,y=rna_bcb)) + geom_point() + geom_smooth(method="lm",se=FALSE) +
theme_classic(base_size=16) + theme(axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black")) +
xlab("plasma cell scRNA signature") + ylab("BCB_RNA")+ scale_y_log10() + annotate(geom="text",x=0.01,y=0.1,label=tmp_ct_string)
ggsave("figure12_outputs/S9n.scplasmapb_vs_rna_bcb.png",g_bcb_rna_vs_scplasmab,height=7,width=7)
knitr::include_graphics("figure12_outputs/S9n.scplasmapb_vs_rna_bcb.png")
```

<img src="figure12_outputs/S9n.scplasmapb_vs_rna_bcb.png" width="75%" height="75%" />

Figure S10

``` {.r}
ct_r <- formatSF(cor.test(log10(dna_clin[dna_clin$not_cell_line,]$nonsilent_snpdnpindel),log10(dna_clin[dna_clin$not_cell_line,]$dna_tcb),method="spearman")$estimate,format="f",digits=3)
```

    ## Warning in cor.test.default(log10(dna_clin[dna_clin$not_cell_line, ]
    ## $nonsilent_snpdnpindel), : Cannot compute exact p-value with ties

``` {.r}
ct_p <- formatSF(cor.test(log10(dna_clin[dna_clin$not_cell_line,]$nonsilent_snpdnpindel),log10(dna_clin[dna_clin$not_cell_line,]$dna_tcb),method="spearman")$p.value,format="e",digits=3)
```

    ## Warning in cor.test.default(log10(dna_clin[dna_clin$not_cell_line, ]
    ## $nonsilent_snpdnpindel), : Cannot compute exact p-value with ties

``` {.r}
g_dna_tmb_dna_tcb <- ggplot(dna_clin[dna_clin$not_cell_line,],aes(x=nonsilent_snpdnpindel,y=dna_tcb)) + geom_point(size=3) + geom_smooth(method="lm",se=FALSE) + theme_classic() +
scale_y_log10() + scale_x_log10() + xlab("TMB") + ylab(expression(TCB[DNA])) +
annotate(geom="text",x=30,y=1,label=paste0("rho = ",ct_r,"\np = ",ct_p,"\nn=",nrow(dna_clin[dna_clin$not_cell_line,])))
ggsave("figure12_outputs/S10a.dna_tmb_dna_tcb.png",g_dna_tmb_dna_tcb,width=7,height=7)
knitr::include_graphics("figure12_outputs/S10a.dna_tmb_dna_tcb.png")
```

<img src="figure12_outputs/S10a.dna_tmb_dna_tcb.png" width="75%" height="75%" />

``` {.r}
ct_r <- formatSF(cor.test(log10(dna_clin[dna_clin$not_cell_line,]$nonsilent_snpdnpindel),log10(dna_clin[dna_clin$not_cell_line,]$dna_bcb),method="spearman")$estimate,format="f",digits=3)
```

    ## Warning in cor.test.default(log10(dna_clin[dna_clin$not_cell_line, ]
    ## $nonsilent_snpdnpindel), : Cannot compute exact p-value with ties

``` {.r}
ct_p <- formatSF(cor.test(log10(dna_clin[dna_clin$not_cell_line,]$nonsilent_snpdnpindel),log10(dna_clin[dna_clin$not_cell_line,]$dna_bcb),method="spearman")$p.value,format="e",digits=3)
```

    ## Warning in cor.test.default(log10(dna_clin[dna_clin$not_cell_line, ]
    ## $nonsilent_snpdnpindel), : Cannot compute exact p-value with ties

``` {.r}
g_dna_tmb_dna_bcb <- ggplot(dna_clin[dna_clin$not_cell_line,],aes(x=nonsilent_snpdnpindel,y=dna_bcb)) + geom_point(size=3) + geom_smooth(method="lm",se=FALSE) + theme_classic() +
scale_y_log10() + scale_x_log10() + xlab("TMB") + ylab(expression(BCB[DNA])) +
annotate(geom="text",x=25,y=0.3,label=paste0("rho = ",ct_r,"\np = ",ct_p,"\nn=",nrow(dna_clin[dna_clin$not_cell_line,])))
ggsave("figure12_outputs/S10b.dna_tmb_dna_bcb.png",g_dna_tmb_dna_bcb,height=7,width=7)
knitr::include_graphics("figure12_outputs/S10b.dna_tmb_dna_bcb.png")
```

<img src="figure12_outputs/S10b.dna_tmb_dna_bcb.png" width="75%" height="75%" />

``` {.r}
makeKMplot_2group(dna_clin[dna_clin$not_cell_line,],"plot_combined_tmb_cohort_dna_bcb_bestrest","figure12_outputs/S10c.dna_tmbbcb_survival_bestrest.png")
```

    ## png 
    ##   2

``` {.r}
knitr::include_graphics("figure12_outputs/S10c.dna_tmbbcb_survival_bestrest.png")
```

<img src="figure12_outputs/S10c.dna_tmbbcb_survival_bestrest.png" width="75%" height="75%" />

``` {.r}
tmp_ft <- fisher.test(table(dna_clin[dna_clin$not_cell_line,]$Response,factor(dna_clin[dna_clin$not_cell_line,]$plot_combined_tmb_cohort_dna_bcb_bestrest,
levels=c("Others","High TMB, High BCB"))))
tmp_ft_pstring <- paste0("Fisher p=",formatSF(tmp_ft$p.value,digits=3)," OR=",formatSF(tmp_ft$estimate,digits=3))
tdf <- melt(table(dna_clin[dna_clin$not_cell_line,]$Response,factor(dna_clin[dna_clin$not_cell_line,]$plot_combined_tmb_cohort_dna_bcb_bestrest,levels=c("High TMB, High BCB","Others"))))
g_rnr_temp <- ggplot(tdf) + geom_mosaic(aes(weight=value,x=product(Var2),fill=Var1),offset=0.01) + scale_fill_manual(values=c("red","blue")) + theme_bw() +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.title=element_blank(),axis.text.x=element_text(size=12),axis.text.y=element_text(size=12)) +
xlab("") + ylab("") + annotate(geom="text",x=0.5,y=1.02,label=tmp_ft_pstring) + scale_y_continuous(breaks=seq(0,1,by=0.25))
ggsave("figure12_outputs/S10d.dna_tmbbcb_bestrest_response.png",g_rnr_temp + geom_text(data = ggplot_build(g_rnr_temp)$data[[1]], aes(x = (xmin+xmax)/2, y = (ymin+ymax)/2, label=.wt)),height=7,width=7)
knitr::include_graphics("figure12_outputs/S10d.dna_tmbbcb_bestrest_response.png")
```

<img src="figure12_outputs/S10d.dna_tmbbcb_bestrest_response.png" width="75%" height="75%" />

``` {.r}
makeKMplot_4group(dna_clin[dna_clin$not_cell_line,],"plot_combined_tmb_cohort_dna_bcb","figure12_outputs/S10e.dna_tmbbcb_survival.png")
```

    ## png 
    ##   2

``` {.r}
knitr::include_graphics("figure12_outputs/S10e.dna_tmbbcb_survival.png")
```

<img src="figure12_outputs/S10e.dna_tmbbcb_survival.png" width="75%" height="75%" />

``` {.r}
dna_clin$plot_combined_tmb_cohort_dna_tcb_cohort_dna_bcb <- NA
dna_clin[dna_clin$not_cell_line,]$plot_combined_tmb_cohort_dna_tcb_cohort_dna_bcb <- paste0(dna_clin[dna_clin$not_cell_line,]$plot_tmb_binary,", ",
dna_clin[dna_clin$not_cell_line,]$plot_cohort_dna_tcb_binary,", ",dna_clin[dna_clin$not_cell_line,]$plot_cohort_dna_bcb_binary)
dna_clin$plot_combined_tmb_cohort_dna_tcb_cohort_dna_bcb_bestrest <- NA
dna_clin[dna_clin$not_cell_line,]$plot_combined_tmb_cohort_dna_tcb_cohort_dna_bcb_bestrest <-
ifelse(dna_clin[dna_clin$not_cell_line,]$plot_combined_tmb_cohort_dna_tcb_cohort_dna_bcb=="High TMB, High TCB, High BCB",
"High TMB, High TCB, High BCB","Others")

makeKMplot_2group(dna_clin[dna_clin$not_cell_line,],"plot_combined_tmb_cohort_dna_tcb_cohort_dna_bcb_bestrest","figure12_outputs/S10f.dna_tmbtcbbcb_survival_bestrest.png")
```

    ## png 
    ##   2

``` {.r}
knitr::include_graphics("figure12_outputs/S10f.dna_tmbtcbbcb_survival_bestrest.png")
```

<img src="figure12_outputs/S10f.dna_tmbtcbbcb_survival_bestrest.png" width="75%" height="75%" />

``` {.r}
tmp_ft <- fisher.test(table(dna_clin[dna_clin$not_cell_line,]$Response,factor(dna_clin[dna_clin$not_cell_line,]$plot_combined_tmb_cohort_dna_tcb_cohort_dna_bcb_bestrest,
levels=c("Others","High TMB, High TCB, High BCB"))))
tmp_ft_pstring <- paste0("Fisher p=",formatSF(tmp_ft$p.value,digits=3)," OR=",formatSF(tmp_ft$estimate,digits=3))
tdf <- melt(table(dna_clin[dna_clin$not_cell_line,]$Response,
factor(dna_clin[dna_clin$not_cell_line,]$plot_combined_tmb_cohort_dna_tcb_cohort_dna_bcb_bestrest,levels=c("High TMB, High TCB, High BCB","Others"))))
g_rnr_temp <- ggplot(tdf) + geom_mosaic(aes(weight=value,x=product(Var2),fill=Var1),offset=0.01) + scale_fill_manual(values=c("red","blue")) + theme_bw() +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.title=element_blank(),axis.text.x=element_text(size=12),axis.text.y=element_text(size=12)) +
xlab("") + ylab("") + annotate(geom="text",x=0.5,y=1.02,label=tmp_ft_pstring) + scale_y_continuous(breaks=seq(0,1,by=0.25))
ggsave("figure12_outputs/S10g.dna_tmbtcbbcb_bestrest_response.png",g_rnr_temp + geom_text(data = ggplot_build(g_rnr_temp)$data[[1]], aes(x = (xmin+xmax)/2, y = (ymin+ymax)/2, label=.wt)),height=7,width=7)
knitr::include_graphics("figure12_outputs/S10g.dna_tmbtcbbcb_bestrest_response.png")
```

<img src="figure12_outputs/S10g.dna_tmbtcbbcb_bestrest_response.png" width="75%" height="75%" />

``` {.r}
long_cbbpalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
dna_clin$class <- dna_clin$plot_combined_tmb_cohort_dna_tcb_cohort_dna_bcb
class_table <- table(dna_clin[dna_clin$not_cell_line,]$class)
dna_tmbtcbbcb_string <- paste0(names(class_table)[1],"=",class_table[[1]]," ",names(class_table)[2],"=",class_table[[2]],"\n",
names(class_table)[3],"=",class_table[[3]]," ",names(class_table)[4],"=",class_table[[4]],"\n",
names(class_table)[5],"=",class_table[[5]]," ",names(class_table)[6],"=",class_table[[6]],"\n",
names(class_table)[7],"=",class_table[[7]]," ",names(class_table)[8],"=",class_table[[8]])
png("figure12_outputs/S10h.dna_tmbtcbbcb_survival.png",height=700,width=700,res=100)
fit<- survfit(surv ~ class, data = dna_clin[dna_clin$not_cell_line,])
ggsurvplot(fit, data = dna_clin[dna_clin$not_cell_line,],pval=TRUE,pval.method=TRUE,legend=c(0.75,0.8),censor=TRUE,palette=long_cbbpalette,title=dna_tmbtcbbcb_string,
font.title=c(10,"plain","black"),font.legend=c(8,"plain","black"))
dev.off()
```

    ## png 
    ##   2

``` {.r}
knitr::include_graphics("figure12_outputs/S10h.dna_tmbtcbbcb_survival.png")
```

<img src="figure12_outputs/S10h.dna_tmbtcbbcb_survival.png" width="75%" height="75%" />

``` {.r}
tmp_chisq <- chisq.test(table(dna_clin[dna_clin$not_cell_line,]$Response,dna_clin[dna_clin$not_cell_line,]$plot_combined_tmb_cohort_dna_tcb_cohort_dna_bcb),simulate.p.value=TRUE)
tmp_ft_pstring <- paste0("Chi squared p=",formatSF(tmp_chisq$p.value,digits=3))
tdf <- melt(table(dna_clin[dna_clin$not_cell_line,]$Response,dna_clin[dna_clin$not_cell_line,]$plot_combined_tmb_cohort_dna_tcb_cohort_dna_bcb))
g_rnr_temp <- ggplot(tdf) + geom_mosaic(aes(weight=value,x=product(Var2),fill=Var1),offset=0.01) + scale_fill_manual(values=c("red","blue")) + theme_bw() +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.title=element_blank(),axis.text.x=element_text(size=12,angle=90),axis.text.y=element_text(size=12)) +
xlab("") + ylab("") + annotate(geom="text",x=0.5,y=1.02,label=tmp_ft_pstring) + scale_y_continuous(breaks=seq(0,1,by=0.25))
ggsave("figure12_outputs/S10i.dna_tmbtcbbcb_response.png",g_rnr_temp + geom_text(data = ggplot_build(g_rnr_temp)$data[[1]], aes(x = (xmin+xmax)/2, y = (ymin+ymax)/2, label=.wt)),height=7,width=7)
knitr::include_graphics("figure12_outputs/S10i.dna_tmbtcbbcb_response.png")
```

<img src="figure12_outputs/S10i.dna_tmbtcbbcb_response.png" width="75%" height="75%" />

``` {.r}
fit1_tmb <- coxph(surv ~ plot_tmb_binary,data=dna_clin[dna_clin$not_cell_line,])
fit1_tcb <- coxph(surv ~ plot_cohort_dna_tcb_binary,data=dna_clin[dna_clin$not_cell_line,])
fit1_bcb <- coxph(surv ~ plot_cohort_dna_bcb_binary,data=dna_clin[dna_clin$not_cell_line,])
fit2_tmbtcb <- coxph(surv ~ plot_tmb_binary + plot_cohort_dna_tcb_binary,data=dna_clin[dna_clin$not_cell_line,])
fit2_tmbbcb <- coxph(surv ~ plot_tmb_binary + plot_cohort_dna_bcb_binary,data=dna_clin[dna_clin$not_cell_line,])
print("TMB vs. TMB+TCB")
```

    ## [1] "TMB vs. TMB+TCB"

``` {.r}
print(anova(fit1_tmb,fit2_tmbtcb))
```

    ## Analysis of Deviance Table
    ##  Cox model: response is  surv
    ##  Model 1: ~ plot_tmb_binary
    ##  Model 2: ~ plot_tmb_binary + plot_cohort_dna_tcb_binary
    ##    loglik  Chisq Df P(>|Chi|)   
    ## 1 -547.15                       
    ## 2 -542.90 8.4893  1  0.003572 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` {.r}
print("TCB vs. TMB+TCB")
```

    ## [1] "TCB vs. TMB+TCB"

``` {.r}
anova(fit1_tcb,fit2_tmbtcb)
```

    ## Analysis of Deviance Table
    ##  Cox model: response is  surv
    ##  Model 1: ~ plot_cohort_dna_tcb_binary
    ##  Model 2: ~ plot_tmb_binary + plot_cohort_dna_tcb_binary
    ##    loglik  Chisq Df P(>|Chi|)  
    ## 1 -545.06                      
    ## 2 -542.90 4.3123  1   0.03784 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` {.r}
print("TMB vs. TMB+BCB")
```

    ## [1] "TMB vs. TMB+BCB"

``` {.r}
anova(fit1_tmb,fit2_tmbbcb)
```

    ## Analysis of Deviance Table
    ##  Cox model: response is  surv
    ##  Model 1: ~ plot_tmb_binary
    ##  Model 2: ~ plot_tmb_binary + plot_cohort_dna_bcb_binary
    ##    loglik  Chisq Df P(>|Chi|)   
    ## 1 -547.15                       
    ## 2 -543.47 7.3498  1  0.006707 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` {.r}
print("BCB vs. TMB+BCB")
```

    ## [1] "BCB vs. TMB+BCB"

``` {.r}
anova(fit1_bcb,fit2_tmbbcb)
```

    ## Analysis of Deviance Table
    ##  Cox model: response is  surv
    ##  Model 1: ~ plot_cohort_dna_bcb_binary
    ##  Model 2: ~ plot_tmb_binary + plot_cohort_dna_bcb_binary
    ##    loglik Chisq Df P(>|Chi|)  
    ## 1 -546.20                     
    ## 2 -543.47 5.459  1   0.01947 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Figure S11

``` {.r}
clin$Malignant.bmode.norm <- clin$Malignant.bmode/clin$Absolute.score..sig.score..bmode
clin$Total.Immune.cells.bmode <- clin$T.cells.CD8.bmode + clin$NK.cells.bmode + clin$Macrophages.bmode + clin$T.cells.CD4.bmode + clin$B.cells.bmode
clin$Total.Immune.cells.bmode.norm <- clin$Total.Immune.cells.bmode/clin$Absolute.score..sig.score..bmode
clin$Total.Stromal.cells.bmode <- clin$Endothelial.cells.bmode + clin$CAF.bmode
clin$Total.Stromal.cells.bmode.norm <- clin$Total.Stromal.cells.bmode/clin$Absolute.score..sig.score..bmode

clin$T.cells.CD8.bmode.norm <- clin$T.cells.CD8.bmode/clin$Absolute.score..sig.score..bmode
clin$T.cells.CD4.bmode.norm <- clin$T.cells.CD4.bmode/clin$Absolute.score..sig.score..bmode
clin$NK.cells.bmode.norm <- clin$NK.cells.bmode/clin$Absolute.score..sig.score..bmode
clin$Macrophages.bmode.norm <- clin$Macrophages.bmode/clin$Absolute.score..sig.score..bmode
clin$Endothelial.cells.bmode.norm <- clin$Endothelial.cells.bmode/clin$Absolute.score..sig.score..bmode
clin$CAF.bmode.norm <- clin$CAF.bmode/clin$Absolute.score..sig.score..bmode

m_cibersort_fractions <- melt(clin[,c("rna_id","Total.Immune.cells.bmode.norm","Total.Stromal.cells.bmode.norm","Malignant.bmode.norm")])
m_cibersort_fractions$rna_id <- factor(m_cibersort_fractions$rna_id,levels=arrange(clin[,c("rna_id","Malignant.bmode.norm")],-Malignant.bmode.norm)$rna_id)

cibersort_pal <- c("#009E73", "#F0E442", "#0072B2")
g_ciberort_bar <- ggplot(m_cibersort_fractions,aes(x=rna_id,y=value)) + geom_bar(stat="identity",aes(fill=variable)) + theme_classic() +
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=4),legend.position="bottom") + scale_fill_manual(values= cibersort_pal)
ggsave("figure12_outputs/S11a.CPB_cibersort_bar.png",g_ciberort_bar,height=7,width=7)
knitr::include_graphics("figure12_outputs/S11a.CPB_cibersort_bar.png")
```

<img src="figure12_outputs/S11a.CPB_cibersort_bar.png" width="75%" height="75%" />

``` {.r}
g_cibersort_cell_fraction <- ggplot(m_cibersort_fractions,aes(x=variable,y=value)) + geom_boxplot(outlier.shape=NA,aes(fill=variable)) + geom_jitter() + theme_classic() +
theme(legend.position="bottom",axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black")) + scale_fill_manual(values=cibersort_pal) +
geom_signif(comparisons=list(c("Malignant.bmode.norm","Total.Immune.cells.bmode.norm"),c("Malignant.bmode.norm","Total.Stromal.cells.bmode.norm"),c("Total.Immune.cells.bmode.norm","Total.Stromal.cells.bmode.norm")), step_increase= .1) + scale_y_continuous(limits=c(0,1.25),breaks=c(0,0.25,0.5,0.75,1))
ggsave("figure12_outputs/S11b.cibersort_fraction_boxplot.png",g_cibersort_cell_fraction,height=7,width=7)
knitr::include_graphics("figure12_outputs/S11b.cibersort_fraction_boxplot.png")
```

<img src="figure12_outputs/S11b.cibersort_fraction_boxplot.png" width="75%" height="75%" />

``` {.r}
tmp_ct <- cor.test(clin[clin$has_WES==TRUE&clin$same_location==TRUE,]$purity,
clin[clin$has_WES==TRUE&clin$same_location==TRUE,]$Malignant.bmode.norm,method="spearman")
```

    ## Warning in cor.test.default(clin[clin$has_WES == TRUE & clin$same_location
    ## == : Cannot compute exact p-value with ties

``` {.r}
tmp_ct_rho <- formatC(as.numeric(tmp_ct$estimate),digits=2)
tmp_ct_p <- formatC(as.numeric(tmp_ct$p.value),format="e",digits=2)
tmp_ct_string <- paste0("rho = ",tmp_ct_rho,"\n"," p = ",tmp_ct_p,"\n n = ",nrow(clin[clin$has_WES==TRUE&clin$same_location==TRUE,]))
g_cibersort_malignant_DNApurity <- ggplot(clin[clin$has_WES==TRUE&clin$same_location==TRUE,],aes(x=purity,y=Malignant.bmode.norm)) + geom_point() +
geom_smooth(method="lm",se=FALSE) + theme_classic(base_size=16) + theme(axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black")) +
xlab("DNA WES tumor purity") + ylab("CIBERSORTx malignant cell fraction") + annotate(geom="text",x=0.25,y=0.8,label=tmp_ct_string)
ggsave("figure12_outputs/S11c.cibersort_malignant_vs_DNApurity.png",g_cibersort_malignant_DNApurity,height=7,width=7)
knitr::include_graphics("figure12_outputs/S11c.cibersort_malignant_vs_DNApurity.png")
```

<img src="figure12_outputs/S11c.cibersort_malignant_vs_DNApurity.png" width="75%" height="75%" />

``` {.r}
tmp_ct <- cor.test(clin$Total.Immune.cells.bmode.norm,clin$Malignant.bmode.norm,method="spearman")
tmp_ct_rho <- formatC(as.numeric(tmp_ct$estimate),digits=2)
tmp_ct_p <- formatC(as.numeric(tmp_ct$p.value),format="e",digits=2)
tmp_ct_string <- paste0("rho = ",tmp_ct_rho,"\n"," p = ",tmp_ct_p,"\n n = ",nrow(clin))
g_cibersort_malignant_immune <- ggplot(clin,aes(x=Total.Immune.cells.bmode.norm,y=Malignant.bmode.norm))+ geom_point() + geom_smooth(method="lm",se=FALSE) +
theme_classic(base_size=16) + theme(axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black")) +
xlab("CIBERSORTx immune cell fraction") + ylab("CIBERSORTx malignant cell fraction")  + annotate(geom="text",x=0.6,y=0.8,label=tmp_ct_string)
ggsave("figure12_outputs/S11d.cibersort_malignant_vs_immune.png",g_cibersort_malignant_immune,height=7,width=7)
knitr::include_graphics("figure12_outputs/S11d.cibersort_malignant_vs_immune.png")
```

<img src="figure12_outputs/S11d.cibersort_malignant_vs_immune.png" width="75%" height="75%" />

``` {.r}
tmp_ct <- cor.test(clin$Total.Stromal.cells.bmode.norm,clin$Malignant.bmode.norm,method="spearman")
tmp_ct_rho <- formatC(as.numeric(tmp_ct$estimate),digits=2)
tmp_ct_p <- formatC(as.numeric(tmp_ct$p.value),format="e",digits=2)
tmp_ct_string <- paste0("rho = ",tmp_ct_rho,"\n"," p = ",tmp_ct_p,"\n n = ",nrow(clin))
g_cibersort_malignant_stromal <- ggplot(clin,aes(x=Total.Stromal.cells.bmode.norm,y=Malignant.bmode.norm)) + geom_point() + geom_smooth(method="lm",se=FALSE) +
theme_classic(base_size=16) + theme(axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black")) +
xlab("CIBERSORTx stromal cell fraction") + ylab("CIBERSORTx malignant cell fraction") + annotate(geom="text",x=0.6,y=0.8,label=tmp_ct_string)
ggsave("figure12_outputs/S11e.cibersort_malignant_vs_stromal.png",g_cibersort_malignant_stromal,height=7,width=7)
knitr::include_graphics("figure12_outputs/S11e.cibersort_malignant_vs_stromal.png")
```

<img src="figure12_outputs/S11e.cibersort_malignant_vs_stromal.png" width="75%" height="75%" />

``` {.r}
tmp_ct <- cor.test(clin$Total.Stromal.cells.bmode.norm,clin$Total.Immune.cells.bmode.norm,method="spearman")
tmp_ct_rho <- formatC(as.numeric(tmp_ct$estimate),digits=2)
tmp_ct_p <- formatC(as.numeric(tmp_ct$p.value),format="e",digits=2)
tmp_ct_string <- paste0("rho = ",tmp_ct_rho,"\n"," p = ",tmp_ct_p,"\n n = ",nrow(clin))
g_cibersort_immune_stromal <- ggplot(clin,aes(x=Total.Immune.cells.bmode.norm,y=Total.Stromal.cells.bmode.norm)) + geom_point() + geom_smooth(method="lm",se=FALSE) +
theme_classic(base_size=16) + theme(axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black")) +
xlab("CIBERSORTx immune cell fraction") + ylab("CIBERSORTx stromal cell fraction") + annotate(geom="text",x=0.6,y=0.6,label=tmp_ct_string)
ggsave("figure12_outputs/S11f.cibersort_immune_vs_stromal.png",g_cibersort_immune_stromal,height=7,width=7)
knitr::include_graphics("figure12_outputs/S11f.cibersort_immune_vs_stromal.png")
```

<img src="figure12_outputs/S11f.cibersort_immune_vs_stromal.png" width="75%" height="75%" />

Figure S12

``` {.r}
skcm_clin_34$plot_tmb_binary <- ifelse(skcm_clin_34$nonsilent_snpdnpindel>median(skcm_clin_34$nonsilent_snpdnpindel),"High TMB","Low TMB")
skcm_clin_34$plot_dna_tcb_binary <- ifelse(skcm_clin_34$dna_tcb>median(skcm_clin_34$dna_tcb),"High TCB DNA","Low TCB DNA")
skcm_clin_34$plot_dna_bcb_binary <- ifelse(skcm_clin_34$dna_bcb>median(skcm_clin_34$dna_bcb),"High BCB DNA","Low BCB DNA")
skcm_clin_34$plot_rna_tcb_binary <- ifelse(skcm_clin_34$rna_tcb>median(skcm_clin_34$rna_tcb),"High TCB RNA","Low TCB RNA")
skcm_clin_34$plot_rna_bcb_binary <- ifelse(skcm_clin_34$rna_bcb>median(skcm_clin_34$rna_bcb),"High BCB RNA","Low BCB RNA")
skcm_clin_34$plot_dna_tcb_dna_bcb_binary <- paste0(skcm_clin_34$plot_dna_tcb_binary,", ",skcm_clin_34$plot_dna_bcb_binary)
skcm_clin_34$plot_rna_tcb_rna_bcb_binary <- paste0(skcm_clin_34$plot_rna_tcb_binary,", ",skcm_clin_34$plot_rna_bcb_binary)


makeKMplot_2group(skcm_clin_34[!is.na(skcm_clin_34$surv),],"plot_tmb_binary","figure12_outputs/S12a.skcm_dna_tmb_survival.png")
```

    ## png 
    ##   2

``` {.r}
knitr::include_graphics("figure12_outputs/S12a.skcm_dna_tmb_survival.png")
```

<img src="figure12_outputs/S12a.skcm_dna_tmb_survival.png" width="75%" height="75%" />

``` {.r}
makeKMplot_2group(skcm_clin_34[!is.na(skcm_clin_34$surv),],"plot_dna_tcb_binary","figure12_outputs/S12b.skcm_dna_tcb_survival.png")
```

    ## png 
    ##   2

``` {.r}
knitr::include_graphics("figure12_outputs/S12b.skcm_dna_tcb_survival.png")
```

<img src="figure12_outputs/S12b.skcm_dna_tcb_survival.png" width="75%" height="75%" />

``` {.r}
makeKMplot_2group(skcm_clin_34[!is.na(skcm_clin_34$surv),],"plot_dna_bcb_binary","figure12_outputs/S12c.skcm_dna_bcb_survival.png")
```

    ## png 
    ##   2

``` {.r}
knitr::include_graphics("figure12_outputs/S12c.skcm_dna_bcb_survival.png")
```

<img src="figure12_outputs/S12c.skcm_dna_bcb_survival.png" width="75%" height="75%" />

``` {.r}
skcm_clin_34$plot_combined_tmb_dna_tcb <- NA
skcm_clin_34$plot_combined_tmb_dna_tcb <- paste0(skcm_clin_34$plot_tmb_binary,", ",skcm_clin_34$plot_dna_tcb_binary)
skcm_clin_34$plot_combined_tmb_dna_tcb_bestrest <- NA
skcm_clin_34$plot_combined_tmb_dna_tcb_bestrest <- ifelse(skcm_clin_34$plot_combined_tmb_dna_tcb=="High TMB, High TCB DNA","High TMB, High TCB DNA","Others")

makeKMplot_2group(skcm_clin_34[!is.na(skcm_clin_34$surv),],"plot_combined_tmb_dna_tcb_bestrest","figure12_outputs/S12d.skcm_dna_tmbtcb_survival_bestrest.png")
```

    ## png 
    ##   2

``` {.r}
knitr::include_graphics("figure12_outputs/S12d.skcm_dna_tmbtcb_survival_bestrest.png")
```

<img src="figure12_outputs/S12d.skcm_dna_tmbtcb_survival_bestrest.png" width="75%" height="75%" />

``` {.r}
makeKMplot_4group(skcm_clin_34[!is.na(skcm_clin_34$surv),],"plot_combined_tmb_dna_tcb","figure12_outputs/S12e.skcm_dna_tmbtcb_survival.png")
```

    ## png 
    ##   2

``` {.r}
knitr::include_graphics("figure12_outputs/S12e.skcm_dna_tmbtcb_survival.png")
```

<img src="figure12_outputs/S12e.skcm_dna_tmbtcb_survival.png" width="75%" height="75%" />

``` {.r}
skcm_clin_34$plot_combined_tmb_dna_bcb <- NA
skcm_clin_34$plot_combined_tmb_dna_bcb <- paste0(skcm_clin_34$plot_tmb_binary,", ",skcm_clin_34$plot_dna_bcb_binary)
skcm_clin_34$plot_combined_tmb_dna_bcb_bestrest <- NA
skcm_clin_34$plot_combined_tmb_dna_bcb_bestrest <- ifelse(skcm_clin_34$plot_combined_tmb_dna_bcb=="High TMB, High BCB DNA","High TMB, High BCB DNA","Others")

makeKMplot_2group(skcm_clin_34[!is.na(skcm_clin_34$surv),],"plot_combined_tmb_dna_bcb_bestrest","figure12_outputs/S12f.skcm_dna_tmbbcb_survival_bestrest.png")
```

    ## png 
    ##   2

``` {.r}
knitr::include_graphics("figure12_outputs/S12f.skcm_dna_tmbbcb_survival_bestrest.png")
```

<img src="figure12_outputs/S12f.skcm_dna_tmbbcb_survival_bestrest.png" width="75%" height="75%" />

``` {.r}
makeKMplot_4group(skcm_clin_34[!is.na(skcm_clin_34$surv),],"plot_combined_tmb_dna_bcb","figure12_outputs/S12g.skcm_dna_tmbbcb_survival.png")
```

    ## png 
    ##   2

``` {.r}
knitr::include_graphics("figure12_outputs/S12g.skcm_dna_tmbbcb_survival.png")
```

<img src="figure12_outputs/S12g.skcm_dna_tmbbcb_survival.png" width="75%" height="75%" />

``` {.r}
skcm_clin_34$plot_combined_dna_tcb_dna_bcb <- NA
skcm_clin_34$plot_combined_dna_tcb_dna_bcb <- paste0(skcm_clin_34$plot_dna_tcb_binary,", ",skcm_clin_34$plot_dna_bcb_binary)
skcm_clin_34$plot_combined_dna_tcb_dna_bcb_bestrest <- NA
skcm_clin_34$plot_combined_dna_tcb_dna_bcb_bestrest <-ifelse(skcm_clin_34$plot_combined_dna_tcb_dna_bcb=="High TCB DNA, High BCB DNA","High TCB DNA, High BCB DNA","Others")

makeKMplot_2group(skcm_clin_34[!is.na(skcm_clin_34$surv),],"plot_combined_dna_tcb_dna_bcb_bestrest","figure12_outputs/S12h.skcm_dna_tcbbcb_survival_bestrest.png")
```

    ## png 
    ##   2

``` {.r}
knitr::include_graphics("figure12_outputs/S12h.skcm_dna_tcbbcb_survival_bestrest.png")
```

<img src="figure12_outputs/S12h.skcm_dna_tcbbcb_survival_bestrest.png" width="75%" height="75%" />

``` {.r}
makeKMplot_4group_withsmalltitle(skcm_clin_34[!is.na(skcm_clin_34$surv),],"plot_combined_dna_tcb_dna_bcb","figure12_outputs/S12i.skcm_dna_tcbbcb_survival.png")
```

    ## png 
    ##   2

``` {.r}
knitr::include_graphics("figure12_outputs/S12i.skcm_dna_tcbbcb_survival.png")
```

<img src="figure12_outputs/S12i.skcm_dna_tcbbcb_survival.png" width="75%" height="75%" />

``` {.r}
skcm_clin_34$plot_combined_tmb_dna_tcb_dna_bcb <- NA
skcm_clin_34$plot_combined_tmb_dna_tcb_dna_bcb <- paste0(skcm_clin_34$plot_tmb_binary,", ",skcm_clin_34$plot_dna_tcb_binary,", ",skcm_clin_34$plot_dna_bcb_binary)
skcm_clin_34$plot_combined_tmb_dna_tcb_dna_bcb_bestrest <- NA
skcm_clin_34$plot_combined_tmb_dna_tcb_dna_bcb_bestrest <-ifelse(skcm_clin_34$plot_combined_tmb_dna_tcb_dna_bcb=="High TMB, High TCB DNA, High BCB DNA",
"High TMB, High TCB DNA, High BCB DNA","Others")

makeKMplot_2group(skcm_clin_34[!is.na(skcm_clin_34$surv),],"plot_combined_tmb_dna_tcb_dna_bcb_bestrest","figure12_outputs/S12j.skcm_dna_tmbtcbbcb_survival_bestrest.png")
```

    ## png 
    ##   2

``` {.r}
knitr::include_graphics("figure12_outputs/S12j.skcm_dna_tmbtcbbcb_survival_bestrest.png")
```

<img src="figure12_outputs/S12j.skcm_dna_tmbtcbbcb_survival_bestrest.png" width="75%" height="75%" />

``` {.r}
skcm_clin_34$class <- skcm_clin_34$plot_combined_tmb_dna_tcb_dna_bcb
class_table <- table(skcm_clin_34[!is.na(skcm_clin_34$surv),]$class)
skcm_dna_tmbtcbbcb_string <- paste0(names(class_table)[1],"=",class_table[[1]]," ",names(class_table)[2],"=",class_table[[2]],"\n",
names(class_table)[3],"=",class_table[[3]]," ",names(class_table)[4],"=",class_table[[4]],"\n",
names(class_table)[5],"=",class_table[[5]]," ",names(class_table)[6],"=",class_table[[6]],"\n",
names(class_table)[7],"=",class_table[[7]]," ",names(class_table)[8],"=",class_table[[8]])
png("figure12_outputs/S12k.skcm_dna_tmbtcbbcb_survival.png",height=700,width=700,res=100)
fit<- survfit(surv ~ class, data = skcm_clin_34[!is.na(skcm_clin_34$surv),])
ggsurvplot(fit, data = skcm_clin_34[!is.na(skcm_clin_34$surv),],pval=TRUE,pval.method=TRUE,legend=c(0.75,0.8),censor=TRUE,palette=long_cbbpalette,title=skcm_dna_tmbtcbbcb_string,
font.title=c(10,"plain","black"),font.legend=c(8,"plain","black"))
dev.off()
```

    ## png 
    ##   2

``` {.r}
knitr::include_graphics("figure12_outputs/S12k.skcm_dna_tmbtcbbcb_survival.png")
```

<img src="figure12_outputs/S12k.skcm_dna_tmbtcbbcb_survival.png" width="75%" height="75%" />

``` {.r}
write.table(skcm_clin_34,"figure12_outputs/skcm_clin_34_table.txt",quote=F,sep="\t")
```

Figure S13

``` {.r}
comb_tcr_tab <- readRDS("figure12_inputs/comb_tcr_tab.rds")
tab_prepost_AB <- readRDS("figure12_inputs/tab_prepost_AB.rds")

tab_prepost_AB$tcr_norm_A <- (comb_tcr_tab[tab_prepost_AB$sample_id_A,]$total_tcr+1)/comb_tcr_tab[tab_prepost_AB$sample_id_A,]$Mapped.Reads*10^6
tab_prepost_AB$tcr_norm_B <- (comb_tcr_tab[tab_prepost_AB$sample_id_B,]$total_tcr+1)/comb_tcr_tab[tab_prepost_AB$sample_id_B,]$Mapped.Reads*10^6
tab_prepost_AB$ig_norm_A <- (comb_tcr_tab[tab_prepost_AB$sample_id_A,]$total_ig+1)/comb_tcr_tab[tab_prepost_AB$sample_id_A,]$Mapped.Reads*10^6
tab_prepost_AB$ig_norm_B <- (comb_tcr_tab[tab_prepost_AB$sample_id_B,]$total_ig+1)/comb_tcr_tab[tab_prepost_AB$sample_id_B,]$Mapped.Reads*10^6

tab_prepost_AB$tcr_change_AB <- ifelse(tab_prepost_AB$tcr_norm_A<tab_prepost_AB$tcr_norm_B,"TCB increase","TCB decrease")
tab_prepost_AB$ig_change_AB <- ifelse(tab_prepost_AB$ig_norm_A<tab_prepost_AB$ig_norm_B,"BCB increase","BCB decrease")

tmp_tab_1 <- tab_prepost_AB[,c("patient","sample_id_A","sample_order_A","cohort","tcr_norm_A","ig_norm_A","priorCTLA4","tcr_change_AB","ig_change_AB","warn")]
tmp_tab_2 <- tab_prepost_AB[,c("patient","sample_id_B","sample_order_B","cohort","tcr_norm_B","ig_norm_B","priorCTLA4","tcr_change_AB","ig_change_AB","warn")]
colnames(tmp_tab_1) <- c("patient","sample_id","sample_order","cohort","tcr_norm","ig_norm","priorCTLA4","tcr_change_AB","ig_change_AB","warn")
colnames(tmp_tab_2) <- c("patient","sample_id","sample_order","cohort","tcr_norm","ig_norm","priorCTLA4","tcr_change_AB","ig_change_AB","warn")
comb_tab_prepost_AB <- rbind(tmp_tab_1,tmp_tab_2)
comb_tab_prepost_AB$sample_order <- ifelse(comb_tab_prepost_AB$sample_order=="A","pre","post")
comb_tab_prepost_AB$sample_order <- factor(comb_tab_prepost_AB$sample_order,levels=c("pre","post"))

dna_tcb <- readRDS("figure12_inputs/dna_tcb.rds")

dna_tab_prepost_good <- readRDS("figure12_inputs/dna_tcb_prepost.rds")

dna_tab_prepost_good$tcr_norm_A <- (dna_tcb[dna_tab_prepost_good$pair_id_pre,]$dna_tcb)
dna_tab_prepost_good$tcr_norm_B <- (dna_tcb[dna_tab_prepost_good$pair_id_post,]$dna_tcb)
dna_tab_prepost_good$ig_norm_A <- (dna_tcb[dna_tab_prepost_good$pair_id_pre,]$dna_bcb)
dna_tab_prepost_good$ig_norm_B <- (dna_tcb[dna_tab_prepost_good$pair_id_post,]$dna_bcb)

dna_tab_prepost_good$tcr_change_AB <- ifelse(dna_tab_prepost_good$tcr_norm_A<dna_tab_prepost_good$tcr_norm_B,"TCB increase","TCB decrease")
dna_tab_prepost_good$ig_change_AB <- ifelse(dna_tab_prepost_good$ig_norm_A<dna_tab_prepost_good$ig_norm_B,"BCB increase","BCB decrease")

dna_tmp_tab_1 <- dna_tab_prepost_good[,c("patient","sample_id_pre","sample_order_pre","cohort","tcr_norm_A","ig_norm_A","prior_CTLA4_pre","tcr_change_AB","ig_change_AB")]
dna_tmp_tab_1$sample_order_pre <- "pre"
dna_tmp_tab_2 <- dna_tab_prepost_good[,c("patient","sample_id_post","sample_order_post","cohort","tcr_norm_B","ig_norm_B","prior_CTLA4_pre","tcr_change_AB","ig_change_AB")]
dna_tmp_tab_2$sample_order_post <- "post"
colnames(dna_tmp_tab_1) <- c("patient","sample_id","sample_order","cohort","tcr_norm","ig_norm","prior_CTLA4_pre","tcr_change_AB","ig_change_AB")
colnames(dna_tmp_tab_2) <- c("patient","sample_id","sample_order","cohort","tcr_norm","ig_norm","prior_CTLA4_pre","tcr_change_AB","ig_change_AB")

dna_comb_tab_prepost_AB <- rbind(dna_tmp_tab_1,dna_tmp_tab_2)
dna_comb_tab_prepost_AB$sample_order <- factor(dna_comb_tab_prepost_AB$sample_order,levels=c("pre","post"))


tmp_wt <- wilcox.test(tab_prepost_AB$tcr_norm_A,tab_prepost_AB$tcr_norm_B,paired=T)
tmp_wt_string <- paste0("Paired Wilcox test\np=",signif(tmp_wt$p.value,3))
g_paired_boxplot_tcr <- ggplot(comb_tab_prepost_AB,aes(x=sample_order,y=tcr_norm)) + geom_boxplot(outlier.shape=NA) + geom_point(aes(color=cohort),alpha=0.3) +
geom_line(aes(group=patient,color=cohort),alpha=0.3) + theme_classic(base_size=16) + scale_color_manual(values=c("red","blue","black")) + scale_y_log10() + ylab("TCB") +
annotate(x=1.5,y=20,geom="text",label=tmp_wt_string,size=5) + theme(legend.position="bottom")
ggsave("figure12_outputs/S13a.i.tcr_prepost_paired_boxplot.png",g_paired_boxplot_tcr,height=7,width=4)
knitr::include_graphics("figure12_outputs/S13a.i.tcr_prepost_paired_boxplot.png")
```

<img src="figure12_outputs/S13a.i.tcr_prepost_paired_boxplot.png" width="75%" height="75%" />

``` {.r}
tmp_wt <- wilcox.test(tab_prepost_AB$tcr_norm_A,tab_prepost_AB$tcr_norm_B,paired=T)
tmp_wt_string <- paste0("Paired Wilcox test\np=",signif(tmp_wt$p.value,3))
g_paired_boxplot_tcr_increasedecrease <- ggplot(comb_tab_prepost_AB,aes(x=sample_order,y=tcr_norm)) + geom_boxplot(outlier.shape=NA) + geom_point(aes(color=tcr_change_AB),alpha=0.3) +
geom_line(aes(group=patient,color=tcr_change_AB),alpha=0.3) + theme_classic(base_size=16) + scale_color_manual(values=c("red","blue")) + scale_y_log10() + ylab("TCB") +
annotate(x=1.5,y=20,geom="text",label=tmp_wt_string,size=5) + theme(legend.position="bottom")
ggsave("figure12_outputs/S13a.ii.tcr_prepost_paired_boxplot_increasedecrease.png",g_paired_boxplot_tcr_increasedecrease,height=7,width=4)
knitr::include_graphics("figure12_outputs/S13a.ii.tcr_prepost_paired_boxplot_increasedecrease.png")
```

<img src="figure12_outputs/S13a.ii.tcr_prepost_paired_boxplot_increasedecrease.png" width="75%" height="75%" />

``` {.r}
tmp_wt <- wilcox.test(tab_prepost_AB$ig_norm_A,tab_prepost_AB$ig_norm_B,paired=T)
tmp_wt_string <- paste0("Paired Wilcox test\np=",signif(tmp_wt$p.value,3))
g_paired_boxplot_ig <- ggplot(comb_tab_prepost_AB,aes(x=sample_order,y=ig_norm)) + geom_boxplot(outlier.shape=NA) + geom_point(aes(color=cohort),alpha=0.3) +
geom_line(aes(group=patient,color=cohort),alpha=0.3) + theme_classic(base_size=16) + scale_color_manual(values=c("red","blue","black")) + scale_y_log10() + ylab("BCB") +
annotate(x=1.5,y=10000,geom="text",label=tmp_wt_string,size=5) + theme(legend.position="bottom")
ggsave("figure12_outputs/S13b.i.ig_prepost_paired_boxplot.png",g_paired_boxplot_ig,height=7,width=4)
knitr::include_graphics("figure12_outputs/S13b.i.ig_prepost_paired_boxplot.png")
```

<img src="figure12_outputs/S13b.i.ig_prepost_paired_boxplot.png" width="75%" height="75%" />

``` {.r}
tmp_wt <- wilcox.test(tab_prepost_AB$ig_norm_A,tab_prepost_AB$ig_norm_B,paired=T)
tmp_wt_string <- paste0("Paired Wilcox test\np=",signif(tmp_wt$p.value,3))
g_paired_boxplot_ig_increasedecrease <- ggplot(comb_tab_prepost_AB,aes(x=sample_order,y=ig_norm)) + geom_boxplot(outlier.shape=NA) + geom_point(aes(color=ig_change_AB),alpha=0.3) +
geom_line(aes(group=patient,color=ig_change_AB),alpha=0.3) + theme_classic(base_size=16) + scale_color_manual(values=c("red","blue")) + scale_y_log10() + ylab("BCB") +
annotate(x=1.5,y=10000,geom="text",label=tmp_wt_string,size=5) + theme(legend.position="bottom")
ggsave("figure12_outputs/S13b.ii.ig_prepost_paired_boxplot_increasedecrease.png",g_paired_boxplot_ig_increasedecrease,height=7,width=4)
knitr::include_graphics("figure12_outputs/S13b.ii.ig_prepost_paired_boxplot_increasedecrease.png")
```

<img src="figure12_outputs/S13b.ii.ig_prepost_paired_boxplot_increasedecrease.png" width="75%" height="75%" />

``` {.r}
tmp_wt <- wilcox.test(dna_tab_prepost_good$tcr_norm_A,dna_tab_prepost_good$tcr_norm_B,paired=T)
tmp_wt_string <- paste0("Paired Wilcox test\np=",signif(tmp_wt$p.value,3))
g_dna_paired_boxplot_tcr <- ggplot(dna_comb_tab_prepost_AB,aes(x=sample_order,y=tcr_norm)) + geom_boxplot(outlier.shape=NA) + geom_point(aes(color=cohort),alpha=0.3) +
geom_line(aes(group=patient,color=cohort),alpha=0.3) + theme_classic(base_size=16) + scale_color_manual(values=c("red","blue","black")) + scale_y_log10() + ylab("TCB") +
annotate(x=1.5,y=0.5,geom="text",label=tmp_wt_string,size=5) + theme(legend.position="bottom")
ggsave("figure12_outputs/S13c.i.dna_tcr_prepost_paired_boxplot.png",g_dna_paired_boxplot_tcr,height=7,width=4.5)
knitr::include_graphics("figure12_outputs/S13c.i.dna_tcr_prepost_paired_boxplot.png")
```

<img src="figure12_outputs/S13c.i.dna_tcr_prepost_paired_boxplot.png" width="75%" height="75%" />

``` {.r}
tmp_wt <- wilcox.test(dna_tab_prepost_good$tcr_norm_A,dna_tab_prepost_good$tcr_norm_B,paired=T)
tmp_wt_string <- paste0("Paired Wilcox test\np=",signif(tmp_wt$p.value,3))
g_dna_paired_boxplot_tcr_increasedecrease <- ggplot(dna_comb_tab_prepost_AB,aes(x=sample_order,y=tcr_norm)) + geom_boxplot(outlier.shape=NA) + geom_point(aes(color=tcr_change_AB),alpha=0.3) +
geom_line(aes(group=patient,color=tcr_change_AB),alpha=0.3) + theme_classic(base_size=16) + scale_color_manual(values=c("red","blue")) + scale_y_log10() + ylab("TCB") +
annotate(x=1.5,y=0.5,geom="text",label=tmp_wt_string,size=5) + theme(legend.position="bottom")
ggsave("figure12_outputs/S13c.ii.dna_tcr_prepost_paired_boxplot_increasedecrease.png",g_dna_paired_boxplot_tcr_increasedecrease,height=7,width=4.5)
knitr::include_graphics("figure12_outputs/S13c.ii.dna_tcr_prepost_paired_boxplot_increasedecrease.png")
```

<img src="figure12_outputs/S13c.ii.dna_tcr_prepost_paired_boxplot_increasedecrease.png" width="75%" height="75%" />

``` {.r}
tmp_wt <- wilcox.test(dna_tab_prepost_good$ig_norm_A,dna_tab_prepost_good$ig_norm_B,paired=T)
tmp_wt_string <- paste0("Paired Wilcox test\np=",signif(tmp_wt$p.value,3))
g_dna_paired_boxplot_ig <- ggplot(dna_comb_tab_prepost_AB,aes(x=sample_order,y=ig_norm)) + geom_boxplot(outlier.shape=NA) + geom_point(aes(color=cohort),alpha=0.3) +
geom_line(aes(group=patient,color=cohort),alpha=0.3) + theme_classic(base_size=16) + scale_color_manual(values=c("red","blue","black")) + scale_y_log10() + ylab("BCB") +
annotate(x=1.5,y=0.5,geom="text",label=tmp_wt_string,size=5) + theme(legend.position="bottom")
ggsave("figure12_outputs/S13d.i.dna_ig_prepost_paired_boxplot.png",g_dna_paired_boxplot_ig,height=7,width=4.5)
knitr::include_graphics("figure12_outputs/S13d.i.dna_ig_prepost_paired_boxplot.png")
```

<img src="figure12_outputs/S13d.i.dna_ig_prepost_paired_boxplot.png" width="75%" height="75%" />

``` {.r}
tmp_wt <- wilcox.test(dna_tab_prepost_good$ig_norm_A,dna_tab_prepost_good$ig_norm_B,paired=T)
tmp_wt_string <- paste0("Paired Wilcox test\np=",signif(tmp_wt$p.value,3))
g_dna_paired_boxplot_ig_increasedecrease <- ggplot(dna_comb_tab_prepost_AB,aes(x=sample_order,y=ig_norm)) + geom_boxplot(outlier.shape=NA) + geom_point(aes(color=ig_change_AB),alpha=0.3) +
geom_line(aes(group=patient,color=ig_change_AB),alpha=0.3) + theme_classic(base_size=16) + scale_color_manual(values=c("red","blue")) + scale_y_log10() + ylab("BCB") +
annotate(x=1.5,y=0.5,geom="text",label=tmp_wt_string,size=5) + theme(legend.position="bottom")
ggsave("figure12_outputs/S13d.ii.dna_ig_prepost_paired_boxplot_increasedecrease.png",g_dna_paired_boxplot_ig_increasedecrease,height=7,width=4.5)
knitr::include_graphics("figure12_outputs/S13d.ii.dna_ig_prepost_paired_boxplot_increasedecrease.png")
```

<img src="figure12_outputs/S13d.ii.dna_ig_prepost_paired_boxplot_increasedecrease.png" width="75%" height="75%" />

``` {.r}
rna_tab_outcomes <- readRDS("figure12_inputs/rna_tab_outcomes.rds")

comb_tab_prepost_AB[comb_tab_prepost_AB$cohort=="MGH",]$priorCTLA4 <- ifelse(rna_tab_outcomes[comb_tab_prepost_AB[comb_tab_prepost_AB$cohort=="MGH",]$patient,]$prior_CTLA4_A=="TRUE","prior_CTLA4","no_CTLA4")
tab_prepost_AB[tab_prepost_AB$cohort=="MGH",]$priorCTLA4 <- ifelse(rna_tab_outcomes[tab_prepost_AB[tab_prepost_AB$cohort=="MGH",]$patient,]$prior_CTLA4_A=="TRUE","prior_CTLA4","no_CTLA4")

comb_tab_prepost_AB$Response_A <- rna_tab_outcomes[comb_tab_prepost_AB$patient,]$Response_A
tab_prepost_AB$Response_A <- rna_tab_outcomes[tab_prepost_AB$patient,]$Response_A
tab_prepost_AB$overall_survival <- rna_tab_outcomes[tab_prepost_AB$patient,]$overall_survival
tab_prepost_AB$DOD <- rna_tab_outcomes[tab_prepost_AB$patient,]$DOD
tab_prepost_AB$surv <- Surv(tab_prepost_AB$overall_survival,tab_prepost_AB$DOD=="DOD")

tmp_wt_1 <- wilcox.test(tab_prepost_AB[tab_prepost_AB$cohort %in% c("MGH","Riaz")&tab_prepost_AB$priorCTLA4=="prior_CTLA4",]$tcr_norm_A,tab_prepost_AB[tab_prepost_AB$cohort %in% c("MGH","Riaz")&tab_prepost_AB$priorCTLA4=="prior_CTLA4",]$tcr_norm_B,paired=T)
tmp_wt_1_string <- paste0("Paired Wilcox test\np=",signif(tmp_wt_1$p.value,3))
tmp_wt_2 <- wilcox.test(tab_prepost_AB[tab_prepost_AB$cohort %in% c("MGH","Riaz")&tab_prepost_AB$priorCTLA4=="no_CTLA4",]$tcr_norm_A,tab_prepost_AB[tab_prepost_AB$cohort %in% c("MGH","Riaz")&tab_prepost_AB$priorCTLA4=="no_CTLA4",]$tcr_norm_B,paired=T)
tmp_wt_2_string <- paste0("Paired Wilcox test\np=",signif(tmp_wt_2$p.value,3))
tmp_label_df <- data.frame(priorCTLA4=c("prior_CTLA4","no_CTLA4"),temp_label=c(tmp_wt_1_string,tmp_wt_2_string),stringsAsFactors=F)

g_paired_boxplot_tcr_mghriaz <- ggplot(comb_tab_prepost_AB[comb_tab_prepost_AB$cohort %in% c("MGH","Riaz"),],aes(x=sample_order,y=tcr_norm)) + geom_boxplot(outlier.shape=NA) +
geom_point(aes(color=cohort),alpha=0.3) + geom_line(aes(group=patient,color=cohort),alpha=0.3) + theme_classic(base_size=16) + scale_color_manual(values=c("blue","black")) +
scale_y_log10() + ylab("TCB") + geom_text(data=tmp_label_df,aes(x=1.5,y=20,label=temp_label),size=4) + facet_wrap(~priorCTLA4) + theme(legend.position="bottom")
ggsave("figure12_outputs/S13e.tcr_prepost_paired_boxplot_MGHRiaz_priorCTLA4.png",g_paired_boxplot_tcr_mghriaz,height=7,width=5)
knitr::include_graphics("figure12_outputs/S13e.tcr_prepost_paired_boxplot_MGHRiaz_priorCTLA4.png")
```

<img src="figure12_outputs/S13e.tcr_prepost_paired_boxplot_MGHRiaz_priorCTLA4.png" width="75%" height="75%" />

``` {.r}
tmp_wt_1 <- wilcox.test(tab_prepost_AB[tab_prepost_AB$cohort %in% c("MGH","Riaz")&tab_prepost_AB$priorCTLA4=="prior_CTLA4",]$ig_norm_A,tab_prepost_AB[tab_prepost_AB$cohort %in% c("MGH","Riaz")&tab_prepost_AB$priorCTLA4=="prior_CTLA4",]$ig_norm_B,paired=T)
tmp_wt_1_string <- paste0("Paired Wilcox test\np=",signif(tmp_wt_1$p.value,3))
tmp_wt_2 <- wilcox.test(tab_prepost_AB[tab_prepost_AB$cohort %in% c("MGH","Riaz")&tab_prepost_AB$priorCTLA4=="no_CTLA4",]$ig_norm_A,tab_prepost_AB[tab_prepost_AB$cohort %in% c("MGH","Riaz")&tab_prepost_AB$priorCTLA4=="no_CTLA4",]$ig_norm_B,paired=T)
tmp_wt_2_string <- paste0("Paired Wilcox test\np=",signif(tmp_wt_2$p.value,3))
tmp_label_df <- data.frame(priorCTLA4=c("prior_CTLA4","no_CTLA4"),temp_label=c(tmp_wt_1_string,tmp_wt_2_string),stringsAsFactors=F)

g_paired_boxplot_ig_mghriaz <- ggplot(comb_tab_prepost_AB[comb_tab_prepost_AB$cohort %in% c("MGH","Riaz"),],aes(x=sample_order,y=ig_norm)) + geom_boxplot(outlier.shape=NA) +
geom_point(aes(color=cohort),alpha=0.3) + geom_line(aes(group=patient,color=cohort),alpha=0.3) + theme_classic(base_size=16) + scale_color_manual(values=c("blue","black")) +
scale_y_log10() + ylab("BCB") + geom_text(data=tmp_label_df,aes(x=1.5,y=10000,label=temp_label),size=4) + facet_wrap(~priorCTLA4) + theme(legend.position="bottom")
ggsave("figure12_outputs/S13f.ig_prepost_paired_boxplot_MGHRiaz_priorCTLA4.png",g_paired_boxplot_ig_mghriaz,height=7,width=5)
knitr::include_graphics("figure12_outputs/S13f.ig_prepost_paired_boxplot_MGHRiaz_priorCTLA4.png")
```

<img src="figure12_outputs/S13f.ig_prepost_paired_boxplot_MGHRiaz_priorCTLA4.png" width="75%" height="75%" />

``` {.r}
tmp_wt_1 <- wilcox.test(dna_tab_prepost_good[dna_tab_prepost_good$prior_CTLA4_pre==TRUE,]$tcr_norm_A,dna_tab_prepost_good[dna_tab_prepost_good$prior_CTLA4_pre==TRUE,]$tcr_norm_B,paired=T)
tmp_wt_1_string <- paste0("Paired Wilcox test\np=",signif(tmp_wt_1$p.value,3))
tmp_wt_2 <- wilcox.test(dna_tab_prepost_good[dna_tab_prepost_good$prior_CTLA4_pre==FALSE,]$tcr_norm_A,dna_tab_prepost_good[dna_tab_prepost_good$prior_CTLA4_pre==FALSE,]$tcr_norm_B,paired=T)
tmp_wt_2_string <- paste0("Paired Wilcox test\np=",signif(tmp_wt_2$p.value,3))
tmp_label_df <- data.frame(prior_CTLA4_pre=c(TRUE,FALSE),temp_label=c(tmp_wt_1_string,tmp_wt_2_string),stringsAsFactors=F)

g_paired_boxplot_tcr_dna <- ggplot(dna_comb_tab_prepost_AB,aes(x=sample_order,y=tcr_norm)) + geom_boxplot(outlier.shape=NA) +
geom_point(aes(color=cohort),alpha=0.3) + geom_line(aes(group=patient,color=cohort),alpha=0.3) + theme_classic(base_size=16) + scale_color_manual(values=c("red","blue","black")) +
scale_y_log10() + ylab("TCB") + geom_text(data=tmp_label_df,aes(x=1.5,y=0.5,label=temp_label),size=4) + facet_wrap(~prior_CTLA4_pre) + theme(legend.position="bottom")
ggsave("figure12_outputs/S13g.dna_tcr_prepost_paired_boxplot_priorCTLA4.png",g_paired_boxplot_tcr_dna,height=7,width=5)
knitr::include_graphics("figure12_outputs/S13g.dna_tcr_prepost_paired_boxplot_priorCTLA4.png")
```

<img src="figure12_outputs/S13g.dna_tcr_prepost_paired_boxplot_priorCTLA4.png" width="75%" height="75%" />

``` {.r}
tmp_wt_1 <- wilcox.test(dna_tab_prepost_good[dna_tab_prepost_good$prior_CTLA4_pre==TRUE,]$ig_norm_A,dna_tab_prepost_good[dna_tab_prepost_good$prior_CTLA4_pre==TRUE,]$ig_norm_B,paired=T)
tmp_wt_1_string <- paste0("Paired Wilcox test\np=",signif(tmp_wt_1$p.value,3))
tmp_wt_2 <- wilcox.test(dna_tab_prepost_good[dna_tab_prepost_good$prior_CTLA4_pre==FALSE,]$ig_norm_A,dna_tab_prepost_good[dna_tab_prepost_good$prior_CTLA4_pre==FALSE,]$ig_norm_B,paired=T)
tmp_wt_2_string <- paste0("Paired Wilcox test\np=",signif(tmp_wt_2$p.value,3))
tmp_label_df <- data.frame(prior_CTLA4_pre=c(TRUE,FALSE),temp_label=c(tmp_wt_1_string,tmp_wt_2_string),stringsAsFactors=F)

g_paired_boxplot_ig_dna <- ggplot(dna_comb_tab_prepost_AB,aes(x=sample_order,y=ig_norm)) + geom_boxplot(outlier.shape=NA) +
geom_point(aes(color=cohort),alpha=0.3) + geom_line(aes(group=patient,color=cohort),alpha=0.3) + theme_classic(base_size=16) + scale_color_manual(values=c("red","blue","black")) +
scale_y_log10() + ylab("BCB") + geom_text(data=tmp_label_df,aes(x=1.5,y=0.5,label=temp_label),size=4) + facet_wrap(~prior_CTLA4_pre) + theme(legend.position="bottom")
ggsave("figure12_outputs/S13h.dna_ig_prepost_paired_boxplot_priorCTLA4.png",g_paired_boxplot_ig_dna,height=7,width=5)
knitr::include_graphics("figure12_outputs/S13h.dna_ig_prepost_paired_boxplot_priorCTLA4.png")
```

<img src="figure12_outputs/S13h.dna_ig_prepost_paired_boxplot_priorCTLA4.png" width="75%" height="75%" />

Figure S14

``` {.r}
makeResponseMosaic_responsevar <- function(df,var1,response_var,png_name) {
f_ft <- fisher.test(table(df[,response_var],df[,var1]))
f_ft_pstring <- paste0("Fisher p=",formatSF(f_ft$p.value,digits=3))
tdf <- melt(table(df[,response_var],df[,var1]))
g_rnr_temp <- ggplot(tdf) + geom_mosaic(aes(weight=value,x=product(Var2),fill=Var1),offset=0.01) + 
scale_fill_manual(values=c("red","blue")) + theme_bw() +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.title=element_blank(),axis.text.x=element_text(size=12),axis.text.y=element_text(size=12)) +
xlab("") + ylab("") + annotate(geom="text",x=0.5,y=1.02,label=f_ft_pstring) + scale_y_continuous(breaks=seq(0,1,by=0.25))
ggsave(png_name,g_rnr_temp + geom_text(data = ggplot_build(g_rnr_temp)$data[[1]], aes(x = (xmin+xmax)/2, y = (ymin+ymax)/2, label=.wt)),height=7,width=7)
}

makeResponseMosaic_responsevar(tab_prepost_AB,"tcr_change_AB","Response_A","figure12_outputs/S14a.tcr_prepost_response_redblue.png")
knitr::include_graphics("figure12_outputs/S14a.tcr_prepost_response_redblue.png")
```

<img src="figure12_outputs/S14a.tcr_prepost_response_redblue.png" width="75%" height="75%" />

``` {.r}
makeResponseMosaic_responsevar(tab_prepost_AB,"ig_change_AB","Response_A","figure12_outputs/S14b.ig_prepost_response_redblue.png")
knitr::include_graphics("figure12_outputs/S14b.ig_prepost_response_redblue.png")
```

<img src="figure12_outputs/S14b.ig_prepost_response_redblue.png" width="75%" height="75%" />

``` {.r}
makeResponseMosaic_responsevar(dna_tab_prepost_good,"tcr_change_AB","Response_pre","figure12_outputs/S14c.dna_tcr_prepost_response_redblue.png")
knitr::include_graphics("figure12_outputs/S14c.dna_tcr_prepost_response_redblue.png")
```

<img src="figure12_outputs/S14c.dna_tcr_prepost_response_redblue.png" width="75%" height="75%" />

``` {.r}
makeResponseMosaic_responsevar(dna_tab_prepost_good,"ig_change_AB","Response_pre","figure12_outputs/S14d.dna_ig_prepost_response_redblue.png")
knitr::include_graphics("figure12_outputs/S14d.dna_ig_prepost_response_redblue.png")
```

<img src="figure12_outputs/S14d.dna_ig_prepost_response_redblue.png" width="75%" height="75%" />

``` {.r}
makeKMplot_2group(tab_prepost_AB,"tcr_change_AB","figure12_outputs/S14e.rna_tcr_prepost_survival.png")
```

    ## png 
    ##   2

``` {.r}
knitr::include_graphics("figure12_outputs/S14e.rna_tcr_prepost_survival.png")
```

<img src="figure12_outputs/S14e.rna_tcr_prepost_survival.png" width="75%" height="75%" />

``` {.r}
makeKMplot_2group(tab_prepost_AB,"ig_change_AB","figure12_outputs/S14f.rna_ig_prepost_survival.png")
```

    ## png 
    ##   2

``` {.r}
knitr::include_graphics("figure12_outputs/S14f.rna_ig_prepost_survival.png")
```

<img src="figure12_outputs/S14f.rna_ig_prepost_survival.png" width="75%" height="75%" />

``` {.r}
dna_tab_prepost_good$surv <- Surv(dna_tab_prepost_good$overall_survival,dna_tab_prepost_good$DOD=="DOD")

makeKMplot_2group(dna_tab_prepost_good,"tcr_change_AB","figure12_outputs/S14g.dna_tcr_prepost_survival.png")
```

    ## png 
    ##   2

``` {.r}
knitr::include_graphics("figure12_outputs/S14g.dna_tcr_prepost_survival.png")
```

<img src="figure12_outputs/S14g.dna_tcr_prepost_survival.png" width="75%" height="75%" />

``` {.r}
makeKMplot_2group(dna_tab_prepost_good,"ig_change_AB","figure12_outputs/S14h.dna_ig_prepost_survival.png")
```

    ## png 
    ##   2

``` {.r}
knitr::include_graphics("figure12_outputs/S14h.dna_ig_prepost_survival.png")
```

<img src="figure12_outputs/S14h.dna_ig_prepost_survival.png" width="75%" height="75%" />

``` {.r}
write.table(tab_prepost_AB[,c("patient","sample_id_A","sample_id_B","cohort","warn","priorCTLA4","tcr_norm_A","tcr_norm_B","ig_norm_A","ig_norm_B","tcr_change_AB","ig_change_AB")],"figure12_outputs/RNA_prepost_tcr_table.txt",quote=F,sep="\t")

write.table(dna_tab_prepost_good,"figure12_outputs/DNA_prepost_tcr_table.txt",quote=F,sep="\t")
```

``` {.r}
library(xtable)
library(data.table)

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
rownames(dna_clin) <- dna_clin$tumor_id
dna_clin$mutrate <- log10(dna_clin$nonsilent_snpdnpindel)

dna_clin$Response_binary <- as.numeric(factor(dna_clin$Response,levels=c("responder","nonresponder")))*1-1

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
df_u_lof$RNR <- dna_clin[df_u_lof$tumor_id,]$Response
df_lof_res <- data.frame(cbind(gene=over3_lof_genes,n_NR=-1,n_R=-1,fisher_p=-1,fisher_or=-1),stringsAsFactors=F)
tot_n_NR <- as.numeric(table(dna_clin$Response)["nonresponder"])
tot_n_R <- as.numeric(table(dna_clin$Response)["responder"])

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
df_u_ns$RNR <- dna_clin[df_u_ns$tumor_id,]$Response
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

png("figure12_outputs/S6b.RNR_lof_qq.png",height=700,width=700,res=100)
qq_pval(as.numeric(df_lof_res$fisher_p))
dev.off()
```

    ## png 
    ##   2

``` {.r}
knitr::include_graphics("figure12_outputs/S6b.RNR_lof_qq.png")
```

<img src="figure12_outputs/S6b.RNR_lof_qq.png" width="75%" height="75%" />

``` {.r}
png("figure12_outputs/S6d.RNR_ns_qq.png",height=700,width=700,res=100)
qq_pval(as.numeric(df_ns_res$fisher_p))
dev.off()
```

    ## png 
    ##   2

``` {.r}
knitr::include_graphics("figure12_outputs/S6d.RNR_ns_qq.png")
```

<img src="figure12_outputs/S6d.RNR_ns_qq.png" width="75%" height="75%" />

``` {.r}
df_lof_res_cox <- data.frame(cbind(gene=over3_lof_genes,n_NR=-1,n_R=-1,cox_p=-1,cox_hr=-1,notmb_cox_p=-1,notmb_cox_hr=-1,logistic_p=-1,logistic_or=-1,notmb_logistic_p=-1,notmb_logistic_or=-1),stringsAsFactors=F)

print("lof cox")
```

    ## [1] "lof cox"

``` {.r}
dna_clin$mutstatus <- 0
for (i in 1:nrow(df_lof_res_cox)) {
#    if (i %% 100 == 0) {print(i)}
    dna_clin$mutstatus <- 0
    tmp_g <- df_lof_res_cox$gene[i]
    n_R_mut <- nrow(df_u_lof[df_u_lof$gene==tmp_g & df_u_lof$RNR=="responder",])
    n_NR_mut <- nrow(df_u_lof[df_u_lof$gene==tmp_g & df_u_lof$RNR=="nonresponder",])
    df_lof_res_cox[i,]$n_NR <- n_NR_mut
    df_lof_res_cox[i,]$n_R <- n_R_mut
    dna_clin[df_u_lof[df_u_lof$gene==tmp_g,]$tumor_id,]$mutstatus <- 1
    tmp_cox <- summary(coxph(dna_clin$surv ~ dna_clin$mutrate + dna_clin$mutstatus))
    df_lof_res_cox[i,]$cox_p <- tmp_cox$coef[2,5]
    df_lof_res_cox[i,]$cox_hr <- tmp_cox$coef[2,2]
    tmp_logistic <- summary(glm(dna_clin$Response_binary ~ dna_clin$mutrate + dna_clin$mutstatus,family=binomial(link="logit")))
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
df_lof_res_cox$cox_p <- as.numeric(df_lof_res_cox$cox_p)
df_lof_res_cox$logistic_p <- as.numeric(df_lof_res_cox$logistic_p)

df_lof_res_cox$cox_q <- p.adjust(df_lof_res_cox$cox_p,method="BH")
df_lof_res_cox$logistic_q <- p.adjust(df_lof_res_cox$logistic_p,method="BH")

png("figure12_outputs/S6a.OS_lof_qq.png",height=700,width=700,res=100)
qq_pval(as.numeric(df_lof_res_cox$cox_p))
dev.off()
```

    ## png 
    ##   2

``` {.r}
knitr::include_graphics("figure12_outputs/S6a.OS_lof_qq.png")
```

<img src="figure12_outputs/S6a.OS_lof_qq.png" width="75%" height="75%" />

``` {.r}
write.table(df_lof_res_cox[,c("gene","cox_hr","cox_p","cox_q","logistic_or","logistic_p","logistic_q")],gzfile("figure12_outputs/lof_table.txt.gz"),quote=F,sep="\t")

print("nsSNV cox")
```

    ## [1] "nsSNV cox"

``` {.r}
df_ns_res_cox <- data.frame(cbind(gene=over3_ns_genes,n_NR=-1,n_R=-1,cox_p=-1,cox_hr=-1,notmb_cox_p=-1,notmb_cox_hr=-1,logistic_p=-1,logistic_or=-1,notmb_logistic_p=-1,notmb_logistic_or=-1),stringsAsFactors=F)

for (i in 1:nrow(df_ns_res)) {
#    if (i %% 1000 == 0) {print(i)}
    dna_clin$mutstatus <- 0
    tmp_g <- df_ns_res_cox$gene[i]
    n_R_mut <- nrow(df_u_ns[df_u_ns$gene==tmp_g & df_u_ns$RNR=="responder",])
    n_NR_mut <- nrow(df_u_ns[df_u_ns$gene==tmp_g & df_u_ns$RNR=="nonresponder",])
    df_ns_res_cox[i,]$n_NR <- n_NR_mut
    df_ns_res_cox[i,]$n_R <- n_R_mut
    dna_clin[df_u_ns[df_u_ns$gene==tmp_g,]$tumor_id,]$mutstatus <- 1
    tmp_cox <- summary(coxph(dna_clin$surv ~ dna_clin$mutrate + dna_clin$mutstatus))
    df_ns_res_cox[i,]$cox_p <- tmp_cox$coef[2,5]
    df_ns_res_cox[i,]$cox_hr <- tmp_cox$coef[2,2]
    tmp_logistic <- summary(glm(dna_clin$Response_binary ~ dna_clin$mutrate + dna_clin$mutstatus,family=binomial(link="logit")))
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


png("figure12_outputs/S6c.OS_ns_qq.png",height=700,width=700,res=100)
qq_pval(as.numeric(df_ns_res_cox$cox_p))
dev.off()
```

    ## png 
    ##   2

``` {.r}
knitr::include_graphics("figure12_outputs/S6c.OS_ns_qq.png")
```

<img src="figure12_outputs/S6c.OS_ns_qq.png" width="75%" height="75%" />

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
    ## [1] grid      methods   stats     graphics  grDevices utils     datasets 
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] data.table_1.14.0 xtable_1.8-2      ggsignif_0.2.0   
    ##  [4] knitr_1.16        reshape2_1.4.3    ggmosaic_0.2.0   
    ##  [7] gridExtra_2.3     psych_1.7.5       rms_5.1-1        
    ## [10] SparseM_1.77      Hmisc_4.0-3       Formula_1.2-3    
    ## [13] lattice_0.20-35   survminer_0.4.3   ggpubr_0.2       
    ## [16] magrittr_2.0.1    survival_2.42-6   ggplot2_3.2.1    
    ## [19] plyr_1.8.6       
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
    ## [25] mvtnorm_1.0-8       scales_0.5.0        km.ci_0.5-2        
    ## [28] KMsurv_0.1-5        MatrixModels_0.4-1  htmlTable_1.9      
    ## [31] tibble_2.1.3        ellipsis_0.3.2      TH.data_1.0-9      
    ## [34] withr_2.4.2         nnet_7.3-12         lazyeval_0.2.2     
    ## [37] mnormt_1.5-5        crayon_1.4.1        polspline_1.1.13   
    ## [40] evaluate_0.11       fansi_0.5.0         nlme_3.1-137       
    ## [43] MASS_7.3-50         foreign_0.8-70      tools_3.4.0        
    ## [46] lifecycle_1.0.0     multcomp_1.4-6      stringr_1.3.1      
    ## [49] plotly_4.9.4.1      munsell_0.5.0       cluster_2.0.7-1    
    ## [52] compiler_3.4.0      rlang_0.4.11        htmlwidgets_1.5.4  
    ## [55] labeling_0.3        base64enc_0.1-3     rmarkdown_1.5      
    ## [58] gtable_0.3.0        codetools_0.2-15    R6_2.5.1           
    ## [61] zoo_1.8-3           dplyr_0.8.4         fastmap_1.1.0      
    ## [64] survMisc_0.5.4      utf8_1.2.2          rprojroot_1.3-2    
    ## [67] stringi_1.2.3       parallel_3.4.0      Rcpp_1.0.7         
    ## [70] vctrs_0.3.8         rpart_4.1-13        acepack_1.4.1      
    ## [73] tidyselect_1.0.0
