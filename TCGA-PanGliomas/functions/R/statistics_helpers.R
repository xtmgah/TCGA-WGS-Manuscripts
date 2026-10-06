# Small, deterministic estimators and strict regression comparison helpers.
# P/q values have relative tolerance only: zero must never stand in for a small P.
numeric_agreement <- function(actual, expected, probability = FALSE) {
  if (length(actual) != length(expected) || !identical(is.na(actual), is.na(expected))) return(FALSE)
  keep <- !is.na(expected); a <- actual[keep]; e <- expected[keep]
  if (any(!is.finite(a)) || any(!is.finite(e))) return(FALSE)
  if (probability && (any(a < 0 | a > 1) || any(e < 0 | e > 1))) return(FALSE)
  all(abs(a-e) <= if (probability) 1e-8 * abs(e) else 1e-10 + 1e-8 * abs(e))
}
probability_column <- function(x) grepl('(^[Pp]$|^p[._]|^pvalue$|^padj|^q_|^fdr$|^fdr_|_pvalue$|_p_value$|adjusted_pvalue$)', x)

pairwise_metrics <- function(d, metric_col, group_col, pairs, metrics) {
  do.call(rbind, lapply(metrics, function(m) {
    rows <- lapply(seq_along(pairs), function(i) {
      g <- pairs[[i]]; x <- d$value[d[[metric_col]] == m & d[[group_col]] == g[1] & is.finite(d$value)]
      y <- d$value[d[[metric_col]] == m & d[[group_col]] == g[2] & is.finite(d$value)]
      stopifnot(length(x) > 0, length(y) > 0)
      data.frame(metric=m, group_1=g[1], group_2=g[2], comparison_order=i,
                 n_group_1=length(x), n_group_2=length(y), pvalue=wilcox.test(x,y,exact=FALSE)$p.value)
    })
    z <- do.call(rbind, rows); z$padj_bh <- p.adjust(z$pvalue, 'BH'); z
  }))
}

# Specimen-level sufficient statistics retain segment-length weighting while
# avoiding the expensive timing model and repeated segment-level joins.
weighted_summary <- function(d, group_col) {
  groups <- unique(as.character(d[[group_col]]))
  do.call(rbind, lapply(groups, function(g) {
    x <- d[d[[group_col]] == g, ]; stopifnot(!anyDuplicated(x$Tumor_Barcode), all(x$segment_weight_sum > 0))
    data.frame(group=g, weighted_mean_time=sum(x$weighted_time_sum)/sum(x$segment_weight_sum),
               n_samples=nrow(x), n_segments=sum(x$n_segments), total_gain_mb=sum(x$segment_weight_sum)/1e6)
  }))
}
weighted_permutation <- function(d, group_col, seed, n_perm=5000L, plus_one=FALSE, two_group=FALSE) {
  # Pin both sampling algorithm and RNG, independent of callers' global options.
  RNGkind('Mersenne-Twister', 'Inversion', 'Rejection'); set.seed(seed)
  stopifnot(!anyDuplicated(d$Tumor_Barcode), all(is.finite(d$weighted_time_sum)), all(d$segment_weight_sum > 0))
  labels <- as.character(d[[group_col]])
  score <- function(g) {
    means <- tapply(d$weighted_time_sum, g, sum)/tapply(d$segment_weight_sum, g, sum)
    if (two_group) unname(means['ASTRO_Group1']-means['ASTRO_Group2']) else diff(range(means))
  }
  observed <- score(labels)
  null <- replicate(n_perm, {
    if (two_group) {
      chosen <- sample(d$Tumor_Barcode, sum(labels == 'ASTRO_Group1'))
      score(ifelse(d$Tumor_Barcode %in% chosen, 'ASTRO_Group1', 'ASTRO_Group2'))
    } else score(sample(labels, length(labels), replace=FALSE))
  })
  exceed <- sum(if (two_group) abs(null) >= abs(observed) else null >= observed)
  list(observed=observed, pvalue=(exceed + as.integer(plus_one))/(n_perm + as.integer(plus_one)),
       n_permutations=n_perm, exceedances=exceed, seed=seed, plus_one=plus_one)
}
