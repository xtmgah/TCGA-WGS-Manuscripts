# Pairwise Wilcoxon comparisons with BH adjustment and optional facet adjustment.
# Quasirandom point placement uses seed 1 by default.
# Required libraries: dplyr, ggplot2 and rlang; other dependencies are namespaced.
plot_group_compare <- function(
    tdata,
    value,                    # unquoted column name for y
    group      = "Study",     # string column name for group on x
    group_name  = NULL,
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
    p_off = FALSE,
    point_off = FALSE
) {
  req_pkgs <- c("dplyr","ggplot2","ggbeeswarm","rstatix","ggpubr","scales","rlang","hrbrthemes")
  missing  <- req_pkgs[!vapply(req_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Please install packages: ", paste(missing, collapse=", "))

  if(is.null(group_name)){ group_name = group}

  has_cow  <- requireNamespace("cowplot", quietly = TRUE)
  panel_border_layer <- if (has_cow) cowplot::panel_border(color = "black",size = 0.3) else ggplot2::theme()

  value_sym <- rlang::ensym(value)
  group_chr <- group
  facet_chr <- facet

  df <- tdata %>% dplyr::filter(!is.na(!!value_sym), !is.na(.data[[group_chr]]))

  n_groups <- df %>% dplyr::pull(.data[[group_chr]]) %>% unique() %>% length()

  label_prefix <- if (n_groups > 2 | (!is.null(facet_chr) & FDR_facet & n_groups == 2) ) "FDR" else "P"

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
        val_col <- intersect(c("p_override","p.adj","p"), names(p_override))
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

    fmt_p <- function(p) ifelse(p < 0.001, sprintf("%s = %.2e", label_prefix, p), sprintf("%s = %.3f", label_prefix, p))
    stat.test <- stat.test %>%
      dplyr::mutate(
        p.label = fmt_p(p.adj),
        col     = ifelse(p.adj < p_sig_cut,
                         if (!is.null(pcol_sig)) pcol_sig else if (exists("ncicolpal", inherits = TRUE)) get("ncicolpal", inherits = TRUE)[1] else "black",
                         "black")
      )


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

  set.seed(seed)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[group_chr]], y = !!value_sym, fill = .data[[group_chr]]))

  if(!point_off){
    p <-  p + ggbeeswarm::geom_quasirandom(pch = 21, size = point_size, width = point_width, color = point_color, stroke = point_stroke, alpha = point_alpha)
  }

  p <- p + ggplot2::geom_boxplot(width = 0.5, fill = NA, color = "gray20", outlier.shape = NA, size = 0.6)

  if (!is.null(hline_y)) {
    p <- p + ggplot2::geom_hline(yintercept = hline_y,
                                 linetype = hline_linetype,
                                 linewidth = hline_size,
                                 color = hline_color)
  }


  if (nrow(stat.test) > 0) {
    if(nrow(stat.test)==1){
      p <- p + ggpubr::stat_pvalue_manual(
        stat.test,
        label        = "p.label",
        y.position   = "y.position",
        xmin         = "xmin",
        xmax         = "xmax",
        color        = "col",
        tip.length   = 0.01,
        label.size   = 4.5, family = "Roboto Condensed",
        hide.ns      = hide_ns
      )
    }else{
      p <- p + ggpubr::stat_pvalue_manual(
        stat.test,
        label        = "p.label",
        y.position   = "y.position",
        xmin         = "xmin",
        xmax         = "xmax",
        color        = "col",
        tip.length   = 0.01,
        label.size   = 4.5, family = "Roboto Condensed",
        step.increase= 0.05,
        hide.ns      = hide_ns
      )


    }
  }

  p <- p +
    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 7)) +
    scale_color_identity() +     # <--- forces ggplot to use your provided hex/named colors
    ptc_theme() +
    panel_border_layer +
    ggplot2::theme(
      axis.text.x      = ggplot2::element_blank(),
      panel.spacing.x  = grid::unit(0.2, "cm"),
      axis.ticks.x     = ggplot2::element_blank(),
      plot.margin      = ggplot2::margin(4, 4, 4, 4),
      legend.position  = "right",
      legend.text      = ggplot2::element_text(size = 14),
      legend.title     = ggplot2::element_text(size = 14),
      legend.box.spacing = unit(x = 0.05,"cm"),
      strip.text.x     = ggplot2::element_text(hjust = 0.5, face = "plain", size = 14)
    ) +
    coord_cartesian(clip = "off")+
    ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size = 4))) +
    ggplot2::labs(x = NULL, y = if (is.null(ylab)) rlang::as_name(value_sym) else ylab,fill=group_name)




  if (!is.null(palette)) {
    p <- p + ggplot2::scale_fill_manual(values = palette)
  } else if (exists("study_color", inherits = TRUE) && is.null(facet_chr) && group_chr == "Study") {
    p <- p + ggplot2::scale_fill_manual(values = get("study_color", inherits = TRUE))
  } else {
    p <- p + ggplot2::scale_fill_brewer(palette = "Set2")
  }

  if (!is.null(facet_chr)) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste0("~", facet_chr)), scales = facet_scales, ncol = facet_ncol)
  }

  if(!is.null(legend_position)){
    p  <- p + theme(legend.position = legend_position)
  }

  if (!is.null(output_file)) {
    output_file <- file.path(ptc_external_directory(dirname(output_file),'Output directory',must_exist=FALSE,create=TRUE),basename(output_file))
    auto_w <- max(4, 2.0 * n_groups)
    w <- if (is.null(width)) auto_w else width
    h <- if (is.null(height)) 5 else height

    ext <- tools::file_ext(output_file)
    dev <- tolower(ext)
    if (dev == "pdf") {
      if (capabilities("cairo")) {
        ggplot2::ggsave(output_file, plot = p, width = w, height = h, device = grDevices::cairo_pdf)
      } else {
        ggplot2::ggsave(output_file, plot = p, width = w, height = h)
      }
    } else {
      ggplot2::ggsave(output_file, plot = p, width = w, height = h)
    }
  }

  return(invisible(list(plot = p, stats = stat.test)))
}
