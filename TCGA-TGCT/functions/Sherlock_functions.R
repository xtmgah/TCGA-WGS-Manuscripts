

set_wd <- function() {
  library(rstudioapi) # make sure you have it installed
  current_path <- getSourceEditorContext()$path
  setwd(dirname(current_path ))
  print( getwd() )
  print("ok")
}


# Load library ------------------------------------------------------------
libztw <- function(...){
  library(tidyverse)
  #library(tidylog)
  library(ggsci)
  library(ggpubr)
  library(hrbrthemes)
  library(scales)
  library(cowplot)
  library(inspectdf)
  library(broom)
  library(janitor)
  library(clipr)
  library(gtsummary)
  library(hablar)
  
  library(conflicted)
  conflict_prefer("fct",winner = 'hablar')
  conflict_prefer("filter",winner = 'dplyr')
  conflict_prefer("select",winner = 'dplyr')
  conflict_prefer("fisher.test",winner = 'janitor')
  
  #fct <- hablar::fct
  
}

myggstyle <- function(myfont='Roboto Condensed',...){
  library(dunnr)
  theme_set(theme_td(base_family = myfont))
  set_geom_fonts(family = myfont,ggrepel = TRUE )
  
}



# Color sets --------------------------------------------------------------
ncicolpal <- c('#BB0E3D','#ff7f00','#4daf4a','#007BBD','#984ea3','#467E84','#14315C','#3D4551','#947100','#4BBFC6','#FACE00','#a65628','#f781bf','#71767A','#b2df8a','#cab2d6','#ffff33','#542788','#bababa','#01665e')




# PDF plotting function ------------------------------------------------------------

pdfhr <- function(...){
  hrbrthemes::import_roboto_condensed()
  d <- read.csv(extrafont:::fonttable_file(), stringsAsFactors = FALSE)
  d[grepl("Light", d$FontName),]$FamilyName <- font_rc_light
  write.csv(d,extrafont:::fonttable_file(), row.names = FALSE)
  extrafont::loadfonts()
}

pdfhr2 <- function(...){
  require(showtext)
  #font_add_google(name = "Roboto Condensed", family = "Roboto Condensed",regular.wt = 400, bold.wt = 700)
  font_add(family = "Roboto Condensed",regular = '/Users/zhangt8/Library/Fonts/RobotoCondensed-Regular.ttf',bold = '/Users/zhangt8/Library/Fonts/RobotoCondensed-Bold.ttf',italic = '/Users/zhangt8/Library/Fonts/RobotoCondensed-Italic.ttf',bolditalic = "/Users/zhangt8/Library/Fonts/RobotoCondensed-BlackItalic.ttf")
  showtext_auto()
  showtext_opts(dpi = 300)
}



# uniform funciton for data visualization


first_cap <- function(x) {
  s <- tolower(x)
  paste0(toupper(substring(s, 1, 1)), substring(s, 2))
} 




barplot_fisher <- function(mdata0,var1name,var2name,samplelist=NULL,filename=NULL, var1lab=NULL, var2lab=NULL,width = 5,height = 8,fisher=TRUE,textcol='white',Pcol='black',var1_num=FALSE){
  #pdfhr2()
  
  if(is.null(var1lab)){ var1lab = var1name}
  if(is.null(var2lab)){ var1lab = var2name}
  print(var1lab)
  
  mdata <- mdata0
  
  if(!is.null(samplelist)){
    mdata <- mdata0 %>% filter(Tumor_Barcode %in% samplelist)
  }
  
  mdata <- mdata %>% select(Tumor_Barcode,one_of(c(var1name,var2name))) %>% drop_na() 
  
  if(fisher){
    tresult <- mdata %>% select(-Tumor_Barcode) %>% table() %>% fisher.test() %>% tidy()
    print(tresult)
    titletext <- paste0("OR = ",round(tresult$estimate,2),"; P = ",scientific(tresult$p.value,3))
  }else{
    titletext <- ""
  }
  
  
  colnames(mdata)[2] <- 'Var1'
  colnames(mdata)[3] <- 'Var2'
  
  if(var1_num){
    vartmp1 <- mdata %>% group_by(Var1) %>% tally() %>% dplyr::mutate(percent=percent_format(accuracy = 0.1)(n/sum(n))) %>% mutate(Var1_Lab=paste0(Var1,' (N=',n,')')) %>% select(Var1,Var1_Lab)
  }else{
    vartmp1 <- mdata %>% group_by(Var1) %>% tally() %>% dplyr::mutate(percent=percent_format(accuracy = 0.1)(n/sum(n))) %>% mutate(Var1_Lab=paste0(Var1,' (',percent,')')) %>% select(Var1,Var1_Lab)
  }
  
  vartmp2 <- mdata %>% group_by(Var2) %>% tally() %>% dplyr::mutate(percent=percent_format(accuracy = 0.1)(n/sum(n))) %>% mutate(Var2_Lab=paste0(Var2,' (',percent,')')) %>% select(Var2,Var2_Lab)
  
  p <-  mdata %>% 
    group_by(Var1,Var2) %>% 
    dplyr::tally()%>%
    dplyr::mutate(percent=n/sum(n)) %>% 
    ungroup() %>% 
    arrange(Var1,Var2) %>% 
    left_join(vartmp1) %>% 
    left_join(vartmp2) %>%
    mutate(Var1_Lab = fct_inorder(Var1_Lab)) %>% 
    mutate(Var2_Lab = fct_inorder(Var2_Lab)) %>% 
    ggplot(aes(x=Var1_Lab, y=n, fill=Var2_Lab))+
    geom_bar(stat="identity", position ="fill",width = 0.8)+
    geom_text(aes(label=paste0(n,'\n',sprintf("%1.1f", percent*100),"%")), position=position_fill(vjust=0.5), colour=textcol,family='Roboto Condensed')+
    scale_y_continuous(breaks = pretty_breaks(),labels = percent_format(),expand = c(0,0))+
    scale_fill_jama()+
    theme_ipsum_rc(base_size = 12,axis_title_just = 'm',axis_title_size = 14,axis_text_size = 12,grid = 'Yy',ticks = FALSE)+
    labs(title = titletext,fill=var2lab,x=var1lab, y = 'Percentage')+
    theme(text = element_text(family = 'Roboto Condensed'),
          plot.title = element_text(size = 15,hjust = 0.5,face = 'plain',family = 'Roboto Condensed',colour =Pcol),
          axis.text.x = element_text(angle = 30,hjust = 1,vjust = 1))+
    coord_cartesian(clip = 'off')
  
  if(!is.null(filename)){
    ggsave(filename = filename,width = width,height = height,device = cairo_pdf )
  }else{
    return(p)
  }
  
}




# plot_group_compare ------------------------------------------------------
# tdata
# A data frame containing at least:
#   
#   one numeric column with the measurement of interest (e.g. age_at_diagnosis, MEDIAN_COVERAGE, BB_Purity),
# 
# one categorical grouping column (default: "Study", but you can change with group=),
# 
# optionally one facet column (e.g. "Type") if you want to facet with facet=.
# value:  An unquoted column name inside tdata that holds the numeric values to compare.

plot_group_compare <- function(
    tdata,
    value,                    # unquoted column name for y
    group      = "Study",     # string column name for group on x
    facet      = NULL,        # optional string col name for facet (e.g. "Type")
    FDR_facet = FALSE,
    ylab       = NULL,        # y-axis label
    palette    = NULL,        # named vector for fill colors; defaults to Set2 if NULL
    pcol_sig   = NULL,        # color for significant brackets; defaults to ncicolpal[1] if available else "black"
    p_sig_cut  = 0.05,        # significance threshold for coloring brackets
    hide_ns    = TRUE,        # hide non-significant brackets
    output_file = NULL,       # path to save (extension decides device); if NULL just returns plot
    width      = NULL,        # plot width (inches); if NULL auto by #groups
    height     = NULL,        # plot height (inches); if NULL defaults to 5
    facet_scales = "fixed",   # scales for facet_wrap
    facet_ncol   = NULL,      # ncol for facet_wrap (optional)
    y_pos_adjust = NULL,      #
    seed       = 1,            # for quasirandom jitter reproducibility
    hline_y        = NULL,
    hline_linetype = "dashed",
    hline_size     = 0.5,
    hline_color    = "gray40",
    point_alpha = 0.6,
    point_size = 2,
    point_stroke = 0.2,
    point_color = 'white',
    point_width = 0.4,
    p_override     = NULL,
    legend_position = NULL,
    p_off = FALSE
) {
  # --- deps ---
  req_pkgs <- c("dplyr","ggplot2","ggbeeswarm","rstatix","ggpubr","scales","rlang","hrbrthemes")
  missing  <- req_pkgs[!vapply(req_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Please install packages: ", paste(missing, collapse=", "))
  
  # cowplot::panel_border is used if available
  has_cow  <- requireNamespace("cowplot", quietly = TRUE)
  panel_border_layer <- if (has_cow) cowplot::panel_border(color = "black",size = 0.7) else ggplot2::theme()
  
  # --- tidy-eval setup ---
  value_sym <- rlang::ensym(value)
  group_chr <- group
  facet_chr <- facet
  
  df <- tdata %>% dplyr::filter(!is.na(!!value_sym), !is.na(.data[[group_chr]]))
  
  # number of groups (per facet if present, but we use total unique for sizing)
  n_groups <- df %>% dplyr::pull(.data[[group_chr]]) %>% unique() %>% length()
  
  # label chooses "P" for <=2 groups, "FDR" for >2 groups
  label_prefix <- if (n_groups > 2 | (!is.null(facet_chr) & FDR_facet & n_groups == 2) ) "FDR" else "P"
  
  # stats: pairwise Wilcoxon (BH)
  if (!is.null(facet_chr)) {
    
    if(FDR_facet & n_groups == 2 ){
      stat.test <- df %>%
        dplyr::group_by(.data[[facet_chr]]) %>%
        rstatix::pairwise_wilcox_test(stats::as.formula(paste0(rlang::as_name(value_sym), " ~ ", group_chr)),
                                      p.adjust.method = "BH") %>%
        ungroup() %>% 
        rstatix::add_xy_position(x = group_chr) %>% 
        mutate(p.adj = p.adjust(p,method = 'BH'))
      
    }else{
      
      stat.test <- df %>%
        dplyr::group_by(.data[[facet_chr]]) %>%
        rstatix::pairwise_wilcox_test(stats::as.formula(paste0(rlang::as_name(value_sym), " ~ ", group_chr)),
                                      p.adjust.method = "BH") %>%
        ungroup() %>% 
        rstatix::add_xy_position(x = group_chr)
    }
  } else {
    stat.test <- df %>%
      rstatix::pairwise_wilcox_test(stats::as.formula(paste0(rlang::as_name(value_sym), " ~ ", group_chr)),
                                    p.adjust.method = "BH") %>%
      rstatix::add_xy_position(x = group_chr)
  }
  
  
  
  # ---- attach medians (and n) per group (and facet if present)
  med_tbl <- df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c(facet_chr, group_chr)))) %>%
    dplyr::summarise(
      median = stats::median(!!value_sym, na.rm = TRUE),
      n = dplyr::n(),
      .groups = "drop"
    )
  
  if (!is.null(facet_chr)) {
    stat.test <- stat.test %>%
      dplyr::left_join(
        med_tbl %>% dplyr::rename(median1 = median, n1 = n, group1 = !!sym(group_chr)),
        by = c("group1", facet_chr)
      ) %>%
      dplyr::left_join(
        med_tbl %>% dplyr::rename(median2 = median, n2 = n, group2 = !!sym(group_chr)),
        by = c("group2", facet_chr)
      )
  } else {
    stat.test <- stat.test %>%
      dplyr::left_join(
        med_tbl %>% dplyr::rename(median1 = median, n1 = n, group1 = !!sym(group_chr)),
        by = "group1"
      ) %>%
      dplyr::left_join(
        med_tbl %>% dplyr::rename(median2 = median, n2 = n, group2 = !!sym(group_chr)),
        by = "group2"
      )
  }
  
  stat.test <- stat.test %>%
    dplyr::mutate(
      median_diff = median1 - median2,
      abs_diff    = abs(median_diff)
    )
  
  if(!p_off){
    
    # ---- OPTIONAL: manual override of adjusted p-values
    if (!is.null(p_override) && nrow(stat.test) > 0) {
      if (is.numeric(p_override) && length(p_override) == 1) {
        stat.test$p.adj <- as.numeric(p_override)
      } else if (is.numeric(p_override) && length(p_override) == nrow(stat.test)) {
        stat.test$p.adj <- as.numeric(p_override)
      } else if (is.function(p_override)) {
        vec <- p_override(stat.test)
        stopifnot(is.numeric(vec), length(vec) == nrow(stat.test))
        stat.test$p.adj <- as.numeric(vec)
      } else if (is.data.frame(p_override)) {
        # try to find a value column
        val_col <- intersect(c("p_override","p.adj","p"), names(p_overrideerride))
        if (length(val_col) == 0) stop("p_override data.frame must contain one of: p_override, p.adj, p")
        val_col <- val_col[1]
        
        join_keys <- c("group1","group2")
        if (!is.null(facet_chr) && facet_chr %in% names(stat.test) && facet_chr %in% names(p_override)) {
          join_keys <- c(join_keys, facet_chr)
        }
        stat.test <- stat.test %>%
          dplyr::left_join(p_override %>%
                             dplyr::select(dplyr::all_of(c(join_keys, val_col))) %>%
                             dplyr::rename(p_adj_override = !!rlang::sym(val_col)),
                           by = join_keys) %>%
          dplyr::mutate(p.adj = dplyr::coalesce(p_adj_override, p.adj)) %>%
          dplyr::select(-dplyr::any_of("p_adj_override"))
      } else {
        stop("Unsupported p_override type. Use numeric (length 1 or nrow), function, or data.frame.")
      }
    }
    
    # p/FDR label formatting and color
    fmt_p <- function(p) ifelse(p < 0.001, sprintf("%s = %.2e", label_prefix, p), sprintf("%s = %.3f", label_prefix, p))
    stat.test <- stat.test %>%
      dplyr::mutate(
        p.label = fmt_p(p.adj),
        col     = ifelse(p.adj < p_sig_cut,
                         if (!is.null(pcol_sig)) pcol_sig else if (exists("ncicolpal", inherits = TRUE)) get("ncicolpal", inherits = TRUE)[1] else "black",
                         "black")
      )
    
    
    # --- y.position adjustment ---
    if (!is.null(y_pos_adjust)) {
      if (is.numeric(y_pos_adjust) && length(y_pos_adjust) != nrow(stat.test)) {
        stat.test$y.position <- stat.test$y.position + y_pos_adjust
      } else if (is.numeric(y_pos_adjust) && length(y_pos_adjust) == nrow(stat.test)) {
        stat.test$y.position <- y_pos_adjust
      } else if (is.function(y_pos_adjust)) {
        stat.test$y.position <- y_pos_adjust(stat.test)
      } else {
        warning("y_pos_adjust must be NULL, a single number, a numeric vector of length nrow(stat.test), or a function.")
      }
    }
    
  }else{
    stat.test <-  data.frame()
  }
  print(stat.test$y.position)
  
  # base plot
  set.seed(seed)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[group_chr]], y = !!value_sym, fill = .data[[group_chr]])) +
    ggbeeswarm::geom_quasirandom(pch = 21, size = point_size, width = point_width, color = point_color, stroke = point_stroke, alpha = point_alpha) +
    ggplot2::geom_boxplot(width = 0.5, fill = NA, color = "gray20", outlier.shape = NA, size = 0.6) 
  
  # NEW: optional horizontal reference line
  if (!is.null(hline_y)) {
    p <- p + ggplot2::geom_hline(yintercept = hline_y,
                                 linetype = hline_linetype,
                                 linewidth = hline_size,
                                 color = hline_color)
  }
  
  # Add brackets if tests exist
  if (nrow(stat.test) > 0) {
    p <- p + ggpubr::stat_pvalue_manual(
      stat.test,
      label        = "p.label",
      y.position   = "y.position",
      xmin         = "xmin",
      xmax         = "xmax",
      color        = "col",
      tip.length   = 0.01,
      step.increase= 0.05,
      hide.ns      = hide_ns
    )
  }
  
  p <- p +
    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 7)) +
    scale_color_identity() +     # <--- forces ggplot to use your provided hex/named colors
    hrbrthemes::theme_ipsum_rc(axis_text_size = 14, axis_title_just = "m",
                               axis_title_size = 16, grid = "XY", ticks = TRUE) +
    panel_border_layer +
    ggplot2::theme(
      axis.text.x      = ggplot2::element_blank(),
      panel.spacing.x  = grid::unit(0.2, "cm"),
      axis.ticks.x     = ggplot2::element_blank(),
      plot.margin      = ggplot2::margin(4, 4, 4, 4),
      legend.position  = "right",
      legend.text      = ggplot2::element_text(size = 14),
      legend.title     = ggplot2::element_text(size = 14),
      strip.text.x     = ggplot2::element_text(hjust = 0.5, face = "plain", size = 14)
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size = 4))) +
    ggplot2::labs(x = NULL, y = if (is.null(ylab)) rlang::as_name(value_sym) else ylab)
  
  
  
  
  # fill scale
  if (!is.null(palette)) {
    p <- p + ggplot2::scale_fill_manual(values = palette)
  } else if (exists("study_color", inherits = TRUE) && is.null(facet_chr) && group_chr == "Study") {
    p <- p + ggplot2::scale_fill_manual(values = get("study_color", inherits = TRUE))
  } else {
    p <- p + ggplot2::scale_fill_brewer(palette = "Set2")
  }
  
  # facet (optional)
  if (!is.null(facet_chr)) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste0("~", facet_chr)), scales = facet_scales, ncol = facet_ncol)
  }
  
  # legend position
  if(!is.null(legend_position)){
    p  <- p + theme(legend.position = legend_position)
  }
  
  # --- saving logic ---
  if (!is.null(output_file)) {
    # Auto height/width if not provided
    # Example rule: width ~ 2.0 * #groups, minimum 4"; height default 5"
    auto_w <- max(4, 2.0 * n_groups)
    w <- if (is.null(width)) auto_w else width
    h <- if (is.null(height)) 5 else height
    
    # choose device by extension
    ext <- tools::file_ext(output_file)
    dev <- tolower(ext)
    # Prefer cairo for PDF if available
    if (dev == "pdf") {
      if (capabilities("cairo")) {
        ggplot2::ggsave(output_file, plot = p, width = w, height = h, device = cairo_pdf)
      } else {
        ggplot2::ggsave(output_file, plot = p, width = w, height = h)
      }
    } else {
      ggplot2::ggsave(output_file, plot = p, width = w, height = h)
    }
  }
  
  return(invisible(list(plot = p, stats = stat.test)))
}




# Define a function for hierarchical clustering
# Example usage
# Assuming `sbs288_activity` is your data frame and "Tumor_Barcode" is the ID column
#ordered_barcodes <- get_ordered_barcodes(sbs288_activity, "Tumor_Barcode")

# Print the result
#print(ordered_barcodes)
get_ordered_barcodes <- function(data, id_column,dmethod='euclidean',hmethod='complete') {
  # Ensure the input is a data frame or tibble
  if (!is.data.frame(data)) {
    stop("Input data must be a data frame or tibble.")
  }
  
  # Check if the specified ID column exists
  if (!id_column %in% colnames(data)) {
    stop("Specified ID column does not exist in the data.")
  }
  
  # Separate the ID column and numeric data
  id_values <- data[[id_column]]  # Extract the ID column
  numeric_data <- data %>%
    select(-all_of(id_column)) %>%  # Remove the ID column
    as.matrix()                    # Convert to a matrix for clustering
  
  # Ensure numeric data is valid
  if (!is.numeric(numeric_data)) {
    stop("All non-ID columns must be numeric for clustering.")
  }
  
  # Perform hierarchical clustering
  dist_matrix <- dist(numeric_data,method = dmethod)   # Compute the distance matrix
  hclust_result <- hclust(dist_matrix,method = hmethod )  # Hierarchical clustering
  
  # Get the ordered IDs based on the clustering
  ordered_ids <- id_values[hclust_result$order]
  
  return(ordered_ids)
}




# PieDonut_ztw ------------------------------------------------------------
require(webr)
require(ggrepel)
PieDonut_ztw <- function (data, mapping, start = getOption("PieDonut.start", 0), addPieLabel = TRUE, addDonutLabel = TRUE, showRatioDonut = TRUE, showRatioPie = TRUE, ratioByGroup = TRUE, showRatioThreshold = getOption("PieDonut.showRatioThreshold", 0.02), labelposition = getOption("PieDonut.labelposition", 2), labelpositionThreshold = 0.1, r0 = getOption("PieDonut.r0", 0.3), r1 = getOption("PieDonut.r1", 1), r2 = getOption("PieDonut.r2", 1.2), explode = NULL, selected = NULL, explodePos = 0.1, ThresholdToLable=FALSE,
                          color = "white", pieAlpha = 0.8, donutAlpha = 1, maxx = NULL, 
                          showPieName = TRUE, showDonutName = FALSE, title = NULL, 
                          pieLabelSize = 4, donutLabelSize = 3, titlesize = 5, explodePie = TRUE, 
                          explodeDonut = FALSE, use.label = TRUE, use.labels = TRUE, mainCol = NULL,showNum = FALSE,returnp=FALSE,hjustvalue=NULL, vjustvalue=NULL,
                          family = getOption("PieDonut.family", "")) 
{
  (cols = colnames(data))
  if (use.labels) 
    data = addLabelDf(data, mapping)
  count <- NULL
  if ("count" %in% names(mapping)) 
    count <- getMapping(mapping, "count")
  count
  pies <- donuts <- NULL
  (pies = getMapping(mapping, "pies"))
  if (is.null(pies)) 
    (pies = getMapping(mapping, "pie"))
  if (is.null(pies)) 
    (pies = getMapping(mapping, "x"))
  (donuts = getMapping(mapping, "donuts"))
  if (is.null(donuts)) 
    (donuts = getMapping(mapping, "donut"))
  if (is.null(donuts)) 
    (donuts = getMapping(mapping, "y"))
  if (!is.null(count)) {
    df <- data %>% group_by(.data[[pies]]) %>% dplyr::summarize(Freq = sum(.data[[count]]))
    df
  }
  else {
    df = data.frame(table(data[[pies]]))
  }
  colnames(df)[1] = pies
  df$end = cumsum(df$Freq)
  df$start = dplyr::lag(df$end)
  df$start[1] = 0
  total = sum(df$Freq)
  df$start1 = df$start * 2 * pi/total
  df$end1 = df$end * 2 * pi/total
  df$start1 = df$start1 + start
  df$end1 = df$end1 + start
  df$focus = 0
  if (explodePie) 
    df$focus[explode] = explodePos
  df$mid = (df$start1 + df$end1)/2
  df$x = ifelse(df$focus == 0, 0, df$focus * sin(df$mid))
  df$y = ifelse(df$focus == 0, 0, df$focus * cos(df$mid))
  df$label = df[[pies]]
  df$ratio = df$Freq/sum(df$Freq)
  
  if(showRatioPie & showNum){
    df$label = ifelse(df$ratio >= showRatioThreshold, paste0(df$label, "\n(",df$Freq, ", ",scales::percent(accuracy = 0.1,df$ratio), ")"), as.character(df$label))
  }else {
    if (showRatioPie) {
      df$label = ifelse(df$ratio >= showRatioThreshold, paste0(df$label, "\n(", scales::percent(accuracy = 0.1,df$ratio), ")"), as.character(df$label))
    }
    
    if (showNum) {
      df$label = ifelse(df$ratio >= showRatioThreshold, paste0(df$label, "\n(", df$Freq, ")"), as.character(df$label))
    }
    
  }
  
  
  df$labelx = (r0 + r1)/2 * sin(df$mid) + df$x
  df$labely = (r0 + r1)/2 * cos(df$mid) + df$y
  if (!is.factor(df[[pies]])) 
    df[[pies]] <- factor(df[[pies]])
  df
  if(is.null(mainCol)){
    mainCol = gg_color_hue(nrow(df)) 
  }
  df$radius = r1
  df$radius[df$focus != 0] = df$radius[df$focus != 0] + df$focus[df$focus != 
                                                                   0]
  if(is.null(hjustvalue)){
    df$hjust = ifelse((df$mid%%(2 * pi)) > pi, 1, 0)
  }else{
    df$hjust = hjustvalue
  }
  
  if(is.null(vjustvalue)){
    df$vjust = ifelse(((df$mid%%(2 * pi)) < (pi/2)) | (df$mid%%(2 * pi) > (pi * 3/2)), 0, 1)
  }else{
    df$vjust = vjustvalue
  }
  
  df$segx = df$radius * sin(df$mid)
  df$segy = df$radius * cos(df$mid)
  df$segxend = (df$radius + 0.05) * sin(df$mid)
  df$segyend = (df$radius + 0.05) * cos(df$mid)
  df
  if (!is.null(donuts)) {
    subColor = makeSubColor(mainCol, no = length(unique(data[[donuts]])))
    subColor
    data
    if (!is.null(count)) {
      df3 <- as.data.frame(data[c(donuts, pies, count)])
      colnames(df3) = c("donut", "pie", "Freq")
      df3
      df3 <- eval(parse(text = "complete(df3,donut,pie)"))
      df3$Freq[is.na(df3$Freq)] = 0
      if (!is.factor(df3[[1]])) 
        df3[[1]] = factor(df3[[1]])
      if (!is.factor(df3[[2]])) 
        df3[[2]] = factor(df3[[2]])
      df3 <- df3 %>% arrange(.data$pie, .data$donut)
      a <- df3 %>% spread(.data$pie, value = .data$Freq)
      a = as.data.frame(a)
      a
      rownames(a) = a[[1]]
      a = a[-1]
      a
      colnames(df3)[1:2] = c(donuts, pies)
    }
    else {
      df3 = data.frame(table(data[[donuts]], data[[pies]]), 
                       stringsAsFactors = FALSE)
      colnames(df3)[1:2] = c(donuts, pies)
      a = table(data[[donuts]], data[[pies]])
      a
    }
    a
    df3
    df3$group = rep(colSums(a), each = nrow(a))
    df3$pie = rep(1:ncol(a), each = nrow(a))
    total = sum(df3$Freq)
    total
    df3$ratio1 = df3$Freq/total
    df3
    if (ratioByGroup) {
      df3$ratio = scales::percent(accuracy = 0.1,df3$Freq/df3$group)
    }
    else {
      df3$ratio <- scales::percent(accuracy = 0.1,df3$ratio1)
    }
    df3$end = cumsum(df3$Freq)
    df3
    df3$start = dplyr::lag(df3$end)
    df3$start[1] = 0
    df3$start1 = df3$start * 2 * pi/total
    df3$end1 = df3$end * 2 * pi/total
    df3$start1 = df3$start1 + start
    df3$end1 = df3$end1 + start
    df3$mid = (df3$start1 + df3$end1)/2
    df3$focus = 0
    if (!is.null(selected)) {
      df3$focus[selected] = explodePos
    }
    else if (!is.null(explode)) {
      selected = c()
      for (i in 1:length(explode)) {
        start = 1 + nrow(a) * (explode[i] - 1)
        selected = c(selected, start:(start + nrow(a) - 
                                        1))
      }
      selected
      df3$focus[selected] = explodePos
    }
    df3
    df3$x = 0
    df3$y = 0
    df
    if (!is.null(explode)) {
      explode
      for (i in 1:length(explode)) {
        xpos = df$focus[explode[i]] * sin(df$mid[explode[i]])
        ypos = df$focus[explode[i]] * cos(df$mid[explode[i]])
        df3$x[df3$pie == explode[i]] = xpos
        df3$y[df3$pie == explode[i]] = ypos
      }
    }
    df3$no = 1:nrow(df3)
    df3$label = df3[[donuts]]
    df3$label0= df3$label
    
    if(showRatioDonut & showNum){
      if (max(nchar(levels(df3$label))) <= 2) 
        df3$label = paste0(df3$label, "(", df3$Freq, ", ",df3$ratio, ")")
      else df3$label = paste0(df3$label, "\n(",df3$Freq,", ", df3$ratio,  ")")
    }else{
      
      if (showRatioDonut) {
        if (max(nchar(levels(df3$label))) <= 2) 
          df3$label = paste0(df3$label, "(", df3$ratio, ")")
        else df3$label = paste0(df3$label, "\n(", df3$ratio,  ")")
      }
      
      if (showNum) {
        df3$label <- df3$label0
        if (max(nchar(levels(df3$label))) <= 2) 
          df3$label = paste0(df3$label, "(", df3$Freq, ")")
        else df3$label = paste0(df3$label, "\n(", df3$Freq,  ")")
      }
      
    }
    
    
    df3$label[df3$ratio1 == 0] = ""
    df3$label[df3$ratio1 < showRatioThreshold] = ""
    
    if(is.null(hjustvalue)){
      df3$hjust = ifelse((df3$mid%%(2 * pi)) > pi, 1, 0)
    }else{
      df3$hjust = hjustvalue
    }
    
    if(is.null(vjustvalue)){
      df3$vjust = ifelse(((df3$mid%%(2 * pi)) < (pi/2)) | (df3$mid%%(2 * pi) > (pi * 3/2)), 0, 1)
    }else{
      df3$vjust=vjustvalue
    }
    df3$no = factor(df3$no)
    df3
    labelposition
    if (labelposition > 0) {
      df3$radius = r2
      if (explodeDonut) 
        df3$radius[df3$focus != 0] = df3$radius[df3$focus != 
                                                  0] + df3$focus[df3$focus != 0]
      df3$segx = df3$radius * sin(df3$mid) + df3$x
      df3$segy = df3$radius * cos(df3$mid) + df3$y
      df3$segxend = (df3$radius + 0.05) * sin(df3$mid) + 
        df3$x
      df3$segyend = (df3$radius + 0.05) * cos(df3$mid) + 
        df3$y
      if (labelposition == 2) 
        df3$radius = (r1 + r2)/2
      df3$labelx = (df3$radius) * sin(df3$mid) + df3$x
      df3$labely = (df3$radius) * cos(df3$mid) + df3$y
    }
    else {
      df3$radius = (r1 + r2)/2
      if (explodeDonut) 
        df3$radius[df3$focus != 0] = df3$radius[df3$focus != 
                                                  0] + df3$focus[df3$focus != 0]
      df3$labelx = df3$radius * sin(df3$mid) + df3$x
      df3$labely = df3$radius * cos(df3$mid) + df3$y
    }
    df3$segx[df3$ratio1 == 0] = 0
    df3$segxend[df3$ratio1 == 0] = 0
    df3$segy[df3$ratio1 == 0] = 0
    df3$segyend[df3$ratio1 == 0] = 0
    if (labelposition == 0) {
      df3$segx[df3$ratio1 < showRatioThreshold] = 0
      df3$segxend[df3$ratio1 < showRatioThreshold] = 0
      df3$segy[df3$ratio1 < showRatioThreshold] = 0
      df3$segyend[df3$ratio1 < showRatioThreshold] = 0
    }
    df3
    del = which(df3$Freq == 0)
    del
    if (length(del) > 0) 
      subColor <- subColor[-del]
    subColor
  }
  p <- ggplot() + theme_no_axes() + coord_fixed()+theme(plot.background = element_blank(),panel.border = element_blank(),plot.margin = margin(-30,-30,-30,-30))
  if (is.null(maxx)) {
    r3 = r2 + 0.3
  }
  else {
    r3 = maxx
  }
  ## change df
  
  if(ThresholdToLable){
    df$label[df$ratio < showRatioThreshold] <- ''
  }
  
  
  p1 <- p + geom_arc_bar(aes_string(x0 = "x", y0 = "y", r0 = as.character(r0), 
                                    r = as.character(r1), start = "start1", end = "end1", 
                                    fill = pies), alpha = pieAlpha, color = color,size=0.1,data = df) + 
    transparent() + scale_fill_manual(values = mainCol) + 
    xlim(r3 * c(-1, 1)) + ylim(r3 * c(-1, 1)) + guides(fill = FALSE)
  if ((labelposition == 1) & (is.null(donuts))) {
    p1 <- p1 + geom_segment(aes_string(x = "segx", y = "segy", 
                                       xend = "segxend", yend = "segyend"), data = df) + 
      geom_text(aes_string(x = "segxend", y = "segyend", 
                           label = "label", hjust = "hjust", vjust = "vjust"), 
                size = pieLabelSize, data = df, family = family)
  }
  else if ((labelposition == 2) & (is.null(donuts))) {
    p1 <- p1 + geom_segment(aes_string(x = "segx", y = "segy", 
                                       xend = "segxend", yend = "segyend"), data = df[df$ratio <  labelpositionThreshold, ]) + geom_text(aes_string(x = "segxend",  y = "segyend", label = "label", hjust = "hjust",  vjust = "vjust"), size = pieLabelSize, data = df[df$ratio < 
                                                                                                                                                                                                                                                                         labelpositionThreshold, ], family = family) + geom_text(aes_string(x = "labelx", y = "labely", label = "label"), size = pieLabelSize,  data = df[df$ratio >= labelpositionThreshold, ],  family = family)
  }
  else {
    p1 <- p1 + geom_text(aes_string(x = "labelx", y = "labely", 
                                    label = "label"), size = pieLabelSize, data = df, 
                         family = family)
  }
  if (showPieName) 
    p1 <- p1 + annotate("text", x = 0, y = 0, label = pies, 
                        size = titlesize, family = family)
  p1 <- p1 + theme(text = element_text(family = family))
  if (!is.null(donuts)) {
    if (explodeDonut) {
      p3 <- p + geom_arc_bar(aes_string(x0 = "x", y0 = "y", 
                                        r0 = as.character(r1), r = as.character(r2), 
                                        start = "start1", end = "end1", fill = "no", 
                                        explode = "focus"), alpha = donutAlpha, color = color, 
                             data = df3)
    }
    else {
      p3 <- p + geom_arc_bar(aes_string(x0 = "x", y0 = "y", 
                                        r0 = as.character(r1), r = as.character(r2), 
                                        start = "start1", end = "end1", fill = "no"), 
                             alpha = donutAlpha, color = color, data = df3)
    }
    p3 <- p3 + transparent() + scale_fill_manual(values = subColor) + 
      xlim(r3 * c(-1, 1)) + ylim(r3 * c(-1, 1)) + guides(fill = FALSE)
    #p3
    if (labelposition == 1) {
      p3 <- p3 + geom_segment(aes_string(x = "segx", y = "segy", 
                                         xend = "segxend", yend = "segyend"), data = df3) + 
        geom_text(aes_string(x = "segxend", y = "segyend", 
                             label = "label", hjust = "hjust", vjust = "vjust"), 
                  size = donutLabelSize, data = df3, family = family)
    }
    else if (labelposition == 0) {
      p3 <- p3 + geom_text(aes_string(x = "labelx", y = "labely", 
                                      label = "label"), size = donutLabelSize, data = df3, 
                           family = family)
    }
    else {
      p3 <- p3 + geom_segment(aes_string(x = "segx", y = "segy", 
                                         xend = "segxend", yend = "segyend"), data = df3[df3$ratio1 < 
                                                                                           labelpositionThreshold, ]) + geom_text(aes_string(x = "segxend", 
                                                                                                                                             y = "segyend", label = "label", hjust = "hjust", 
                                                                                                                                             vjust = "vjust"), size = donutLabelSize, data = df3[df3$ratio1 < 
                                                                                                                                                                                                   labelpositionThreshold, ], family = family) + 
        geom_text(aes_string(x = "labelx", y = "labely", 
                             label = "label"), size = donutLabelSize, data = df3[df3$ratio1 >= 
                                                                                   labelpositionThreshold, ], family = family)
    }
    if (!is.null(title)) 
      p3 <- p3 + annotate("text", x = 0, y = r3, label = title, 
                          size = titlesize, family = family)
    else if (showDonutName) 
      p3 <- p3 + annotate("text", x = (-1) * r3, y = r3, 
                          label = donuts, hjust = 0, size = titlesize, 
                          family = family)
    p3 <- p3 + theme(text = element_text(family = family))
    
    if(returnp){
      return(list(p1,p3))
    }else{
      grid.newpage()
      print(p1, vp = viewport(height = 1, width = 1))
      print(p3, vp = viewport(height = 1, width = 1))
    }
    
    
    
  } else {
    if(returnp){
      return(p1)
    }else{
      p1
      #df
    }
  }
}

environment(PieDonut_ztw) <- asNamespace("webr")



## PieDonut_ztw2



# PieDonut_ztw ------------------------------------------------------------
PieDonut_ztw2<- function (data, mapping, start = getOption("PieDonut.start", 0), addPieLabel = TRUE, addDonutLabel = TRUE, showRatioDonut = TRUE, showRatioPie = TRUE, ratioByGroup = TRUE, showRatioThreshold = getOption("PieDonut.showRatioThreshold", 0.02), labelposition = getOption("PieDonut.labelposition", 2), labelpositionThreshold = 0.1, r0 = getOption("PieDonut.r0", 0.3), r1 = getOption("PieDonut.r1", 1), r2 = getOption("PieDonut.r2", 1.2), explode = NULL, selected = NULL, explodePos = 0.1, 
                          color = "white", pieAlpha = 0.8, donutAlpha = 1, maxx = NULL, 
                          showPieName = TRUE, showDonutName = FALSE, title = NULL, 
                          pieLabelSize = 4, donutLabelSize = 3, titlesize = 5, explodePie = TRUE, 
                          explodeDonut = FALSE, use.label = TRUE, use.labels = TRUE, mainCol = NULL,showNum = FALSE,returnp=FALSE,
                          family = getOption("PieDonut.family", "")) 
{
  (cols = colnames(data))
  if (use.labels) 
    data = addLabelDf(data, mapping)
  count <- NULL
  if ("count" %in% names(mapping)) 
    count <- getMapping(mapping, "count")
  count
  pies <- donuts <- NULL
  (pies = getMapping(mapping, "pies"))
  if (is.null(pies)) 
    (pies = getMapping(mapping, "pie"))
  if (is.null(pies)) 
    (pies = getMapping(mapping, "x"))
  (donuts = getMapping(mapping, "donuts"))
  if (is.null(donuts)) 
    (donuts = getMapping(mapping, "donut"))
  if (is.null(donuts)) 
    (donuts = getMapping(mapping, "y"))
  if (!is.null(count)) {
    df <- data %>% group_by(.data[[pies]]) %>% dplyr::summarize(Freq = sum(.data[[count]]))
    df
  }
  else {
    df = data.frame(table(data[[pies]]))
  }
  colnames(df)[1] = pies
  df$end = cumsum(df$Freq)
  df$start = dplyr::lag(df$end)
  df$start[1] = 0
  total = sum(df$Freq)
  df$start1 = df$start * 2 * pi/total
  df$end1 = df$end * 2 * pi/total
  df$start1 = df$start1 + start
  df$end1 = df$end1 + start
  df$focus = 0
  if (explodePie) 
    df$focus[explode] = explodePos
  df$mid = (df$start1 + df$end1)/2
  df$x = ifelse(df$focus == 0, 0, df$focus * sin(df$mid))
  df$y = ifelse(df$focus == 0, 0, df$focus * cos(df$mid))
  df$label = df[[pies]]
  df$ratio = df$Freq/sum(df$Freq)
  
  if(showRatioPie & showNum){
    df$label = ifelse(df$ratio >= showRatioThreshold, paste0(df$label, "\n(",df$Freq, ", ",scales::percent(accuracy = 0.1,df$ratio), ")"), as.character(df$label))
  }else {
    if (showRatioPie) {
      df$label = ifelse(df$ratio >= showRatioThreshold, paste0(df$label, "\n(", scales::percent(accuracy = 0.1,df$ratio), ")"), as.character(df$label))
    }
    
    if (showNum) {
      df$label = ifelse(df$ratio >= showRatioThreshold, paste0(df$label, "\n(", df$Freq, ")"), as.character(df$label))
    }
    
  }
  
  
  df$labelx = (r0 + r1)/2 * sin(df$mid) + df$x
  df$labely = (r0 + r1)/2 * cos(df$mid) + df$y
  if (!is.factor(df[[pies]])) 
    df[[pies]] <- factor(df[[pies]])
  df
  if(is.null(mainCol)){
    mainCol = gg_color_hue(nrow(df)) 
  }
  df$radius = r1
  df$radius[df$focus != 0] = df$radius[df$focus != 0] + df$focus[df$focus != 
                                                                   0]
  df$hjust = ifelse((df$mid%%(2 * pi)) > pi, 1, 0)
  df$vjust = ifelse(((df$mid%%(2 * pi)) < (pi/2)) | (df$mid%%(2 * 
                                                                pi) > (pi * 3/2)), 0, 1)
  df$segx = df$radius * sin(df$mid)
  df$segy = df$radius * cos(df$mid)
  df$segxend = (df$radius + 0.05) * sin(df$mid)
  df$segyend = (df$radius + 0.05) * cos(df$mid)
  df
  if (!is.null(donuts)) {
    subColor = makeSubColor(mainCol, no = length(unique(data[[donuts]])))
    subColor
    data
    if (!is.null(count)) {
      df3 <- as.data.frame(data[c(donuts, pies, count)])
      colnames(df3) = c("donut", "pie", "Freq")
      df3
      df3 <- eval(parse(text = "complete(df3,donut,pie)"))
      df3$Freq[is.na(df3$Freq)] = 0
      if (!is.factor(df3[[1]])) 
        df3[[1]] = factor(df3[[1]])
      if (!is.factor(df3[[2]])) 
        df3[[2]] = factor(df3[[2]])
      df3 <- df3 %>% arrange(.data$pie, .data$donut)
      a <- df3 %>% spread(.data$pie, value = .data$Freq)
      a = as.data.frame(a)
      a
      rownames(a) = a[[1]]
      a = a[-1]
      a
      colnames(df3)[1:2] = c(donuts, pies)
    }
    else {
      df3 = data.frame(table(data[[donuts]], data[[pies]]), 
                       stringsAsFactors = FALSE)
      colnames(df3)[1:2] = c(donuts, pies)
      a = table(data[[donuts]], data[[pies]])
      a
    }
    a
    df3
    df3$group = rep(colSums(a), each = nrow(a))
    df3$pie = rep(1:ncol(a), each = nrow(a))
    total = sum(df3$Freq)
    total
    df3$ratio1 = df3$Freq/total
    df3
    if (ratioByGroup) {
      df3$ratio = scales::percent(accuracy = 0.1,df3$Freq/df3$group)
    }
    else {
      df3$ratio <- scales::percent(accuracy = 0.1,df3$ratio1)
    }
    df3$end = cumsum(df3$Freq)
    df3
    df3$start = dplyr::lag(df3$end)
    df3$start[1] = 0
    df3$start1 = df3$start * 2 * pi/total
    df3$end1 = df3$end * 2 * pi/total
    df3$start1 = df3$start1 + start
    df3$end1 = df3$end1 + start
    df3$mid = (df3$start1 + df3$end1)/2
    df3$focus = 0
    if (!is.null(selected)) {
      df3$focus[selected] = explodePos
    }
    else if (!is.null(explode)) {
      selected = c()
      for (i in 1:length(explode)) {
        start = 1 + nrow(a) * (explode[i] - 1)
        selected = c(selected, start:(start + nrow(a) - 
                                        1))
      }
      selected
      df3$focus[selected] = explodePos
    }
    df3
    df3$x = 0
    df3$y = 0
    df
    if (!is.null(explode)) {
      explode
      for (i in 1:length(explode)) {
        xpos = df$focus[explode[i]] * sin(df$mid[explode[i]])
        ypos = df$focus[explode[i]] * cos(df$mid[explode[i]])
        df3$x[df3$pie == explode[i]] = xpos
        df3$y[df3$pie == explode[i]] = ypos
      }
    }
    df3$no = 1:nrow(df3)
    df3$label = df3[[donuts]]
    df3$label0= df3$label
    
    if(showRatioDonut & showNum){
      if (max(nchar(levels(df3$label))) <= 2) 
        df3$label = paste0(df3$label, "(", df3$Freq, ", ",df3$ratio, ")")
      else df3$label = paste0(df3$label, "\n(",df3$Freq,", ", df3$ratio,  ")")
    }else{
      
      if (showRatioDonut) {
        if (max(nchar(levels(df3$label))) <= 2) 
          df3$label = paste0(df3$label, "(", df3$ratio, ")")
        else df3$label = paste0(df3$label, "\n(", df3$ratio,  ")")
      }
      
      if (showNum) {
        df3$label <- df3$label0
        if (max(nchar(levels(df3$label))) <= 2) 
          df3$label = paste0(df3$label, "(", df3$Freq, ")")
        else df3$label = paste0(df3$label, "\n(", df3$Freq,  ")")
      }
      
    }
    
    
    df3$label[df3$ratio1 == 0] = ""
    df3$label[df3$ratio1 < showRatioThreshold] = ""
    df3$hjust = ifelse((df3$mid%%(2 * pi)) > pi, 1, 0)
    df3$vjust = ifelse(((df3$mid%%(2 * pi)) < (pi/2)) | 
                         (df3$mid%%(2 * pi) > (pi * 3/2)), 0, 1)
    df3$no = factor(df3$no)
    df3
    labelposition
    if (labelposition > 0) {
      df3$radius = r2
      if (explodeDonut) 
        df3$radius[df3$focus != 0] = df3$radius[df3$focus != 
                                                  0] + df3$focus[df3$focus != 0]
      df3$segx = df3$radius * sin(df3$mid) + df3$x
      df3$segy = df3$radius * cos(df3$mid) + df3$y
      df3$segxend = (df3$radius + 0.05) * sin(df3$mid) + 
        df3$x
      df3$segyend = (df3$radius + 0.05) * cos(df3$mid) + 
        df3$y
      if (labelposition == 2) 
        df3$radius = (r1 + r2)/2
      df3$labelx = (df3$radius) * sin(df3$mid) + df3$x
      df3$labely = (df3$radius) * cos(df3$mid) + df3$y
    }
    else {
      df3$radius = (r1 + r2)/2
      if (explodeDonut) 
        df3$radius[df3$focus != 0] = df3$radius[df3$focus != 
                                                  0] + df3$focus[df3$focus != 0]
      df3$labelx = df3$radius * sin(df3$mid) + df3$x
      df3$labely = df3$radius * cos(df3$mid) + df3$y
    }
    df3$segx[df3$ratio1 == 0] = 0
    df3$segxend[df3$ratio1 == 0] = 0
    df3$segy[df3$ratio1 == 0] = 0
    df3$segyend[df3$ratio1 == 0] = 0
    if (labelposition == 0) {
      df3$segx[df3$ratio1 < showRatioThreshold] = 0
      df3$segxend[df3$ratio1 < showRatioThreshold] = 0
      df3$segy[df3$ratio1 < showRatioThreshold] = 0
      df3$segyend[df3$ratio1 < showRatioThreshold] = 0
    }
    df3
    del = which(df3$Freq == 0)
    del
    if (length(del) > 0) 
      subColor <- subColor[-del]
    subColor
  }
  p <- ggplot() + theme_no_axes() + coord_fixed()
  if (is.null(maxx)) {
    r3 = r2 + 0.3
  }
  else {
    r3 = maxx
  }
  p1 <- p + geom_arc_bar(aes_string(x0 = "x", y0 = "y", r0 = as.character(r0), 
                                    r = as.character(r1), start = "start1", end = "end1", 
                                    fill = pies), alpha = pieAlpha, color = color, data = df) + 
    transparent() + scale_fill_manual(values = mainCol) + 
    xlim(r3 * c(-1, 1)) + ylim(r3 * c(-1, 1)) + guides(fill = FALSE)
  if ((labelposition == 1) & (is.null(donuts))) {
    p1 <- p1 + geom_segment(aes_string(x = "segx", y = "segy", 
                                       xend = "segxend", yend = "segyend"), data = df) + 
      geom_text_repel(aes_string(x = "segxend", y = "segyend", 
                                 label = "label", hjust = "hjust", vjust = "vjust"), 
                      size = pieLabelSize, data = df, family = family)
  }
  else if ((labelposition == 2) & (is.null(donuts))) {
    p1 <- p1 + geom_segment(aes_string(x = "segx", y = "segy", 
                                       xend = "segxend", yend = "segyend"), data = df[df$ratio <  labelpositionThreshold, ]) + geom_text_repel(aes_string(x = "segxend",  y = "segyend", label = "label", hjust = "hjust",  vjust = "vjust"), size = pieLabelSize, data = df[df$ratio < 
                                                                                                                                                                                                                                                                               labelpositionThreshold, ], family = family) + geom_text_repel(aes_string(x = "labelx", y = "labely", label = "label"), size = pieLabelSize,  data = df[df$ratio >= labelpositionThreshold, ],  family = family)
  }
  else {
    p1 <- p1 + geom_text_repel(aes_string(x = "labelx", y = "labely", 
                                          label = "label"), size = pieLabelSize, data = df, 
                               family = family)
  }
  if (showPieName) 
    p1 <- p1 + annotate("text", x = 0, y = 0, label = pies, 
                        size = titlesize, family = family)
  p1 <- p1 + theme(text = element_text(family = family))
  if (!is.null(donuts)) {
    if (explodeDonut) {
      p3 <- p + geom_arc_bar(aes_string(x0 = "x", y0 = "y", 
                                        r0 = as.character(r1), r = as.character(r2), 
                                        start = "start1", end = "end1", fill = "no", 
                                        explode = "focus"), alpha = donutAlpha, color = color, 
                             data = df3)
    }
    else {
      p3 <- p + geom_arc_bar(aes_string(x0 = "x", y0 = "y", 
                                        r0 = as.character(r1), r = as.character(r2), 
                                        start = "start1", end = "end1", fill = "no"), 
                             alpha = donutAlpha, color = color, data = df3)
    }
    p3 <- p3 + transparent() + scale_fill_manual(values = subColor) + 
      xlim(r3 * c(-1, 1)) + ylim(r3 * c(-1, 1)) + guides(fill = FALSE)
    #p3
    if (labelposition == 1) {
      p3 <- p3 + geom_segment(aes_string(x = "segx", y = "segy", 
                                         xend = "segxend", yend = "segyend"), data = df3) + 
        geom_text_repel(aes_string(x = "segxend", y = "segyend", 
                                   label = "label", hjust = "hjust", vjust = "vjust"), 
                        size = donutLabelSize, data = df3, family = family)
    }
    else if (labelposition == 0) {
      p3 <- p3 + geom_text_repel(aes_string(x = "labelx", y = "labely", 
                                            label = "label"), size = donutLabelSize, data = df3, 
                                 family = family)
    }
    else {
      p3 <- p3 + geom_segment(aes_string(x = "segx", y = "segy", 
                                         xend = "segxend", yend = "segyend"), data = df3[df3$ratio1 < 
                                                                                           labelpositionThreshold, ]) + geom_text_repel(aes_string(x = "segxend", 
                                                                                                                                                   y = "segyend", label = "label", hjust = "hjust", 
                                                                                                                                                   vjust = "vjust"), size = donutLabelSize, data = df3[df3$ratio1 < 
                                                                                                                                                                                                         labelpositionThreshold, ], family = family) + 
        geom_text_repel(aes_string(x = "labelx", y = "labely", 
                                   label = "label"), size = donutLabelSize, data = df3[df3$ratio1 >= 
                                                                                         labelpositionThreshold, ], family = family)
    }
    if (!is.null(title)) 
      p3 <- p3 + annotate("text", x = 0, y = r3, label = title, 
                          size = titlesize, family = family)
    else if (showDonutName) 
      p3 <- p3 + annotate("text", x = (-1) * r3, y = r3, 
                          label = donuts, hjust = 0, size = titlesize, 
                          family = family)
    p3 <- p3 + theme(text = element_text(family = family))
    
    if(returnp){
      return(list(p1,p3))
    }else{
      grid.newpage()
      print(p1, vp = viewport(height = 1, width = 1))
      print(p3, vp = viewport(height = 1, width = 1))
    }
    
    
    
  } else {
    if(returnp){
      return(p1)
    }else{
      p1
    }
  }
}

environment(PieDonut_ztw2) <- asNamespace("webr")








plot_mutational_signature_landscape <- function(
    tumor_activity_all,
    tumor_activity_all_ratio,
    tumor_decompsite,
    sigcol,
    group_data = NULL,
    samlevs = NULL,
    sortby_signatures = NULL,
    sortby_cutoff = 0.1,
    output_file = NULL,
    width = 18,
    height = 6,
    legend_sig_nrow =1
) {
  # req pkgs
  require(dplyr)
  require(tidyr)
  require(ggplot2)
  require(scales)
  require(ggnewscale)
  require(gtools)
  require(hrbrthemes)   # for theme_ipsum
  require(rlang)
  require(tibble)
  require(viridisLite)  # for viridis scale (via viridis_c)
  
  # ---- helpers ----
  round_up_multiple <- function(x, bases = c(2,3,4,5)) {
    if (is.infinite(x) || is.na(x)) return(1)
    candidates <- sapply(bases, function(b) ceiling(x / b) * b)
    min(candidates)
  }
  
  # Ensure required columns exist
  stopifnot("Tumor_Barcode" %in% colnames(tumor_activity_all))
  stopifnot("Tumor_Barcode" %in% colnames(tumor_activity_all_ratio))
  stopifnot("Tumor_Barcode" %in% colnames(tumor_decompsite))
  
  # Clean up group_data if provided
  if (!is.null(group_data)) {
    # If first two columns are barcode+group (as in your snippet), coerce their names
    if (!all(c("Tumor_Barcode","Group") %in% colnames(group_data))) {
      if (ncol(group_data) >= 2) {
        colnames(group_data)[1:2] <- c("Tumor_Barcode","Group")
      } else {
        stop("group_data must have at least two columns: Tumor_Barcode and Group")
      }
    }
    # Keep only needed cols, drop duplicates
    group_data <- group_data %>%
      select(Tumor_Barcode, Group) %>%
      distinct()
  }
  
  # ---- long-format combined data ----
  t_ratio <- tumor_activity_all_ratio %>%
    tidyr::pivot_longer(-Tumor_Barcode) %>%
    mutate(type = "Ratio")
  
  t_count <- tumor_activity_all %>%
    tidyr::pivot_longer(-Tumor_Barcode) %>%
    group_by(Tumor_Barcode) %>%
    summarise(value = sum(value), .groups = "drop") %>%
    mutate(value = log2(value),
           type = "Count")
  
  tdata <- bind_rows(t_ratio, t_count)
  
  # Join group info only if provided
  if (!is.null(group_data)) {
    tdata <- tdata %>% left_join(group_data, by = "Tumor_Barcode")
  }
  
  # Signature ordering for facet bars
  tdata <- tdata %>%
    mutate(
      name = factor(
        name,
        levels = gtools::mixedsort(names(tumor_activity_all)[names(tumor_activity_all) != "Tumor_Barcode"])
      )
    )
  
  
  
  if(!is.null(samlevs)){
    if(!("Group" %in% colnames(samlevs))){
      samlevs <-  samlevs %>% left_join(group_data)
    }
  }else{
    
    # ---- sample ordering (clustering-based) ----
    if (is.null(group_data)) {
      ordered_barcodes <- get_ordered_barcodes(tumor_activity_all_ratio, "Tumor_Barcode")
      samlev_data <- tibble(Tumor_Barcode = ordered_barcodes) %>%
        mutate(seq_old = seq_along(Tumor_Barcode),
               Group  = 1)
    } else {
      group_values <- sort(unique(group_data$Group))
      ordered_barcodes <- NULL
      for (gvalue in group_values) {
        group_sample <- group_data %>%
          filter(Group == gvalue) %>%
          pull(Tumor_Barcode)
        ordered_barcodes_tmp <- get_ordered_barcodes(
          tumor_activity_all_ratio %>% filter(Tumor_Barcode %in% group_sample),
          "Tumor_Barcode"
        )
        ordered_barcodes <- c(ordered_barcodes, ordered_barcodes_tmp)
      }
      samlev_data <- tibble(Tumor_Barcode = ordered_barcodes) %>%
        mutate(seq_old = seq_along(Tumor_Barcode)) %>%
        left_join(group_data, by = "Tumor_Barcode")
    }
    
    # Optionally re-sort within groups by select signatures
    if (is.null(sortby_signatures)) {
      samlevs <- samlev_data %>%
        group_by(Group) %>%
        mutate(seq = seq_along(Tumor_Barcode)) %>%
        ungroup()
    } else {
      tmpsig <- tumor_activity_all_ratio %>%
        select(Tumor_Barcode, tidyselect::any_of(sortby_signatures)) %>%
        rowwise() %>%
        mutate(Summary = sum(c_across(-Tumor_Barcode), na.rm = TRUE)) %>%
        ungroup() %>%
        mutate(Signature = Summary > sortby_cutoff) %>%
        select(Tumor_Barcode, Signature)
      
      samlevs <- samlev_data %>%
        left_join(tmpsig, by = "Tumor_Barcode") %>%
        arrange(Signature, seq_old) %>%
        group_by(Group) %>%
        mutate(seq = seq_along(Tumor_Barcode)) %>%
        ungroup() %>%
        select(-Signature, -seq_old)
    }
  }
  
  # Merge for plotting
  if(!is.null(group_data)){
    pdata <- tdata %>% left_join(samlevs, by = c("Tumor_Barcode","Group"))
  }else{
    pdata <- tdata %>% left_join(samlevs, by = c("Tumor_Barcode"))
  }
  
  
  # scale for dual-axis alignment
  max_val <- suppressWarnings(max(pdata$value, na.rm = TRUE))
  scale_y <- round_up_multiple(max_val)
  if (!is.finite(scale_y) || is.na(scale_y)) scale_y <- 1
  
  pdata_ratio <- pdata %>%
    filter(type == "Ratio") %>%
    mutate(value = value * scale_y)
  
  pdata_count <- pdata %>% filter(type == "Count")
  
  # cosine similarity strip
  pdata_cos <- tumor_decompsite %>% 
    dplyr::select(Tumor_Barcode, Cosine_Similarity) %>% 
    dplyr::left_join(samlevs, by = "Tumor_Barcode") %>%
    { 
      # ensure a Group column exists (for the no-group case)
      if (!"Group" %in% names(.)) dplyr::mutate(., Group = 1) else .
    }
  
  
  # ---- plot ----
  p <- ggplot() +
    geom_bar(
      data = pdata_ratio,
      aes(x = seq, y = value, fill = name),
      stat = "identity", width = 1
    ) +
    scale_fill_manual(values = sigcol, breaks = names(sigcol)) +
    geom_point(
      data = pdata_count,
      aes(x = seq, y = value),
      inherit.aes = FALSE, col = "#6baed6", size = 2
    ) +
    geom_point(
      data = pdata_count,
      aes(x = seq, y = value),
      inherit.aes = FALSE, col = "black", size = 0.1
    ) +
    geom_line(
      data = pdata_count,
      aes(x = seq, y = value),
      inherit.aes = FALSE, col = "white", linewidth = 0.15
    ) +
    scale_y_continuous(
      name = "Mutational signature proportion",
      expand = c(0, 0),
      breaks = seq(0, scale_y, length.out = 6),
      labels = seq(0, 1, length.out = 6),
      sec.axis = sec_axis(~ ., breaks = seq(0, scale_y, length.out = 6),
                          name = "Number of mutations (log2)")
    ) +
    scale_x_continuous(breaks = scales::pretty_breaks(), expand = c(0, 0)) +
    theme_ipsum(base_size = 12, grid = FALSE, ticks = TRUE,
                axis_title_just = "m", axis_title_size = 14) +
    theme(
      axis.title.y.right = element_text(color = "#08519c"),
      axis.text.y.right  = element_text(color = "#08519c"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      strip.text.x = element_text(hjust = 0.5, face = 2),
      panel.spacing = unit(0.1, "cm"),
      legend.position = "bottom",
      text = element_text(family = "Roboto Condensed"),
      plot.margin = ggplot2::margin(4, 4, 4, 4),
      legend.margin = margin(t = -4)
    ) +
    labs(x = NULL, fill = "Signature") +
    guides(fill = guide_legend(nrow = legend_sig_nrow, direction = "horizontal")) +
    ggnewscale::new_scale_fill() +
    geom_tile(
      data = pdata_cos,
      aes(x = seq, y = -0.4, fill = Cosine_Similarity),
      height = 0.5
    ) +
    scale_fill_viridis_c(limits = c(0, 1)) +
    guides(fill = guide_colorbar(
      title = "Cosine similarity\n",
      keywidth = unit(8, "cm"),
      keyheight = unit(0.4, "cm")
    ))
  
  if (!is.null(group_data)) {
    p <- p + facet_grid(~ Group, scales = "free_x",space = 'free_x')
  }
  
  # ---- optional PDF output ----
  if (!is.null(output_file)) {
    ggsave(filename = output_file, plot = p, width = width, height = height, units = "in")
  }
  
  # Return both objects in either case
  return(list(p = p, samlevs = samlevs))
}



# plot_scatter ------------------------------------------------------------
plot_scatter <- function(
    data,
    x, y,
    group = NULL,
    facet = FALSE,
    facet_scales = "free",
    facet_nrow = NULL,
    point_color = "black",
    point_fill  = "#007BBD",
    point_size  = 3,
    point_stroke = 0.4,
    x_lab = NULL, y_lab = NULL,
    # --- NEW: relative label positions in [0,1] ---
    label_x_rel = 0.03,              # 3% from left by default
    label_y_rel = 0.95,              # 5% down from top by default
    alpha_sig = 0.05,
    label_size = 5,
    p_override = NULL,
    label_override = NULL,
    label_color_override = NULL,
    add_reg_eq = FALSE,
    reg_eq_label_y_rel = NULL,       # optional: equation label Y (relative)
    palette = NULL
) {
  require(rlang)
  x <- enquo(x); y <- enquo(y); grp_quosure <- enquo(group)
  has_group <- !quo_is_null(grp_quosure)
  
  df <- data %>% mutate(.row_id = row_number())
  
  aes_base <- aes(x = !!x, y = !!y)
  aes_points <- if (has_group & !is.null(palette) & length(palette)>1) aes(fill = !!grp_quosure) else NULL
  
  if(has_group & !is.null(palette) & length(palette)>1){
    p <- ggplot(df, aes_base) +
      geom_point(
        mapping = aes_points,
        pch = 21, size = point_size, stroke = point_stroke,
        color  = point_color
      )
  }else{
    p <- ggplot(df, aes_base) +
      geom_point(
        pch = 21, size = point_size, stroke = point_stroke,
        color  = point_color,
        fill  = point_fill
      )
  }
  
  p <- p + geom_smooth(method = "lm", se = TRUE) +
    scale_x_continuous(breaks = pretty_breaks(n = 7)) +
    scale_y_continuous(breaks = pretty_breaks(n = 7)) +
    labs(x = x_lab, y = y_lab) +
    theme_ipsum_rc(axis_title_just = "m", axis_title_size = 14, ticks = TRUE, grid = "Y") +
    theme(
      panel.spacing = unit(0.4, "lines"),
      strip.background = element_blank(),
      strip.text = element_text(hjust = 0.5, face = "bold"),
      plot.margin = margin(4, 4, 0, 4),
      legend.position = "none"
    ) +
    ggpubr::border(color = "black", size = 0.35) +
    coord_cartesian(clip = "off")
  
  if (facet && has_group) p <- p + facet_wrap(vars(!!grp_quosure), scales = facet_scales,nrow = facet_nrow)
  if (has_group && !is.null(palette)) p <- p + scale_fill_manual(values = palette)
  
  # ---- per-group stats & ranges ----
  if (has_group) {
    stats_df <- df %>%
      group_by(!!grp_quosure) %>%
      summarise(
        r = cor(!!x, !!y, use = "complete.obs", method = "pearson"),
        pval = tryCatch(cor.test(!!x, !!y, method = "pearson")$p.value, error = function(e) NA_real_),
        ymin = min(!!y, na.rm = TRUE), ymax = max(!!y, na.rm = TRUE),
        xmin = min(!!x, na.rm = TRUE), xmax = max(!!x, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(group_label = as.character(!!grp_quosure))
  } else {
    stats_df <- tibble(
      r    = cor(df %>% pull(!!x), df %>% pull(!!y), use = "complete.obs", method = "pearson"),
      pval = tryCatch(cor.test(df %>% pull(!!x), df %>% pull(!!y), method = "pearson")$p.value, error = function(e) NA_real_),
      ymin = min(df %>% pull(!!y), na.rm = TRUE), ymax = max(df %>% pull(!!y), na.rm = TRUE),
      xmin = min(df %>% pull(!!x), na.rm = TRUE), xmax = max(df %>% pull(!!x), na.rm = TRUE),
      group_label = "All"
    )
  }
  
  # optional p-value override
  if (!is.null(p_override)) {
    if (length(p_override) == 1L) {
      stats_df$pval <- as.numeric(p_override)
    } else if (!is.null(names(p_override))) {
      stats_df$pval <- as.numeric(p_override[stats_df$group_label])
    }
  }
  
  stats_df <- stats_df %>%
    mutate(
      sig_color = ifelse(!is.na(pval) & pval < alpha_sig, "#BB0E3D", "black"),
      r_txt = sprintf("R = %.2f", r),
      p_txt = paste0("P = ", ifelse(pval < .001, scientific_format(digits = 3)(pval), formatC(pval, format = "f", digits = 3))),
      label = paste(r_txt, p_txt, sep = ", ")
    )
  
  if (!is.null(label_override)) {
    if (length(label_override) == 1L) {
      stats_df$label <- as.character(label_override)
    } else if (!is.null(names(label_override))) {
      stats_df$label <- as.character(label_override[stats_df$group_label])
    }
  }
  
  if (!is.null(label_color_override)) {
    if (length(label_color_override) == 1L) {
      stats_df$sig_color <- as.character(label_color_override)
    } else if (!is.null(names(label_color_override))) {
      stats_df$sig_color <- as.character(label_color_override[stats_df$group_label])
    }
  }
  
  
  # --- helpers to resolve relative values to each facet ---
  resolve_rel <- function(rel, key, n) {
    # rel: NULL / scalar / vector / named vector keyed by facet label
    if (is.null(rel)) return(rep(NA_real_, n))
    if (!is.null(names(rel))) {
      out <- as.numeric(rel[key])
    } else if (length(rel) == 1L) {
      out <- rep(rel, n)
    } else {
      out <- rep(rel, length.out = n)
    }
    # clamp to [0,1] just in case
    pmax(0, pmin(1, out))
  }
  
  # positions in RELATIVE [0,1]
  rel_x <- resolve_rel(label_x_rel, stats_df$group_label, nrow(stats_df))
  rel_y <- resolve_rel(label_y_rel, stats_df$group_label, nrow(stats_df))
  
  # convert to DATA coordinates per panel
  stats_df$label_x <- stats_df$xmin + rel_x * (stats_df$xmax - stats_df$xmin)
  stats_df$label_y <- stats_df$ymin + rel_y * (stats_df$ymax - stats_df$ymin)
  
  # ---- add correlation text ----
  p <- p + geom_text(
    data = stats_df,
    aes(x = label_x, y = label_y, label = label),
    color = stats_df$sig_color,
    hjust = 0, size = label_size,
    inherit.aes = FALSE
  )
  
  print(stats_df)
  
  # ---- optional regression equation ----
  if (add_reg_eq) {
    eq_rel_y <- if (is.null(reg_eq_label_y_rel)) {
      # 10% below the main label
      pmax(0, rel_y - 0.10)
    } else if (length(reg_eq_label_y_rel) == 1L) {
      rep(pmax(0, pmin(1, reg_eq_label_y_rel)), nrow(stats_df))
    } else if (!is.null(names(reg_eq_label_y_rel))) {
      as.numeric(pmax(0, pmin(1, reg_eq_label_y_rel[stats_df$group_label])))
    } else {
      pmax(0, pmin(1, rep(reg_eq_label_y_rel, length.out = nrow(stats_df))))
    }
    
    eq_y <- stats_df$ymin + eq_rel_y * (stats_df$ymax - stats_df$ymin)
    
    if (has_group) {
      for (i in seq_len(nrow(stats_df))) {
        p <- p + ggpubr::stat_regline_equation(
          data = df %>% filter(!!grp_quosure == stats_df$group_label[i]),
          mapping  = aes(x = !!x, y = !!y),     # <-- supply x & y
          label.y = eq_y[i],
          label.x = stats_df$label_x[i],
          inherit.aes = FALSE,
          size=5
        )
      }
    } else {
      p <- p + ggpubr::stat_regline_equation(
        mapping  = aes(x = !!x, y = !!y),       # <-- supply x & y
        label.y = eq_y[1],
        label.x = stats_df$label_x[1],
        inherit.aes = FALSE,
        size=5
      )
    }
  }
  
  return(p)
}


# plot_scatter(
#   data = tdata %>% filter(!is.na(age_at_exposure), name =='TL_Ratio'),
#   x = age_at_exposure, y = value,
#   #group = name,
#   #facet = TRUE, facet_scales = "free", palette = ncicolpal[1],
#   #point_color = "black",
#   #point_fill = "#007BBD",
#   x_lab = "Age at exposure", y_lab = "Telomere length",
#   label_x_rel = 0.4,
#   label_y_rel=0,               # one per facet (recycled if lengths differ)
#   alpha_sig = 0.05,
#   label_override = ("β = 0.01, P = 6.70e-04"),
#   label_color_override = ncicolpal[1]
#   # p_override = c(GroupA = 0.012, GroupB = 0.18),  # optional override by group
#   #add_reg_eq = FALSE                      # set TRUE to show the equation
# )
