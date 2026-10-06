# Stepwise Kaplan–Meier confidence-band helper extracted from the original project.
step_ribbon_df <- function(df) {
    df <- arrange(df, .data$time)
    out <- vector("list", nrow(df))
    prev_lower <- df$lower[1]
    prev_upper <- df$upper[1]
    prev_surv <- df$surv[1]
    out[[1]] <- df[1, ]
    for (i in seq(2, nrow(df))) {
        before <- df[i, ]
        before$lower <- prev_lower
        before$upper <- prev_upper
        before$surv <- prev_surv
        after <- df[i, ]
        out[[i]] <- bind_rows(before, after)
        prev_lower <- df$lower[i]
        prev_upper <- df$upper[i]
        prev_surv <- df$surv[i]
    }
    bind_rows(out)
}
