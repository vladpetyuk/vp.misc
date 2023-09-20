#' Plot AUC
#'
#' Plot AUC after LOOCV Model Evaluation
#'
#' @param modelingResult output from either
#'          \code{lr_modeling} or \code{rf_modeling}
#' @param CI (logical) Plot confidence interval or not. Default is False.
#' @param ... further arguments passes to the \code{plot} method of the
#'             \code{\link[pROC]{ci.se}} object of the \code{pROC} package.
#'
#'
#' @importFrom ROCR performance
#' @importFrom graphics abline
#' @importFrom pROC roc ci.se
#'
#' @export plotAUC
#'
plotAUC <- function(modelingResult, CI = FALSE, ...) {
   if (!CI) {
      perf <- performance(modelingResult$pred, "tpr", "fpr")
      # x=1-spec, y=sens
      plot(perf,
         main = sprintf("AUC: %s", round(modelingResult$auc, 2)),
         col = 2, lwd = 2
      )
      abline(a = 0, b = 1, lwd = 2, lty = 2, col = "gray")
   } else {
      old_par <- par()
      par(pty = "s")
      pROC_obj <- roc(modelingResult$pred@labels[[1]],
         modelingResult$pred@predictions[[1]],
         direction = "<",
         smoothed = T,
         # arguments for ci
         ci = TRUE,
         ci.alpha = 0.95,
         stratified = FALSE,
         # arguments for plot
         plot = TRUE,
         auc.polygon = F,
         max.auc.polygon = F,
         grid = F,
         print.auc = TRUE,
         print.auc.y = 0.1,
         print.auc.x = 0.8,
         show.thres = TRUE
      )
      sens.ci <- ci.se(pROC_obj)
      plot(sens.ci, type = "shape", col = "#FF8888", conf = 95, ...)
      abline(a = 1, b = -1, lty = 2, lwd = 2, col = "grey50")
      par(old_par)
   }
}

#' Plot AUC (with ggplot2)
#'
#' Plot AUC after LOOCV Model Evaluation (with ggplot2)
#'
#' @param modelingResult output from either
#'          \code{lr_modeling} or \code{rf_modeling}
#' @param CI (logical) Plot confidence interval or not. Default is False.
#' @param rectilinear (logical) Prevents diagonal lines being formed from jumps in TPR from CI boundaries. (Technical -- transforms segment y = mx + b between fpr1 and fpr2 to the line x = avg(fpr1, fpr2) between y = m*fpr1 + b and y = m*fpr2 + b. Surrounding horizontal segments are extended to this new vertical segment.)
#' @param ... further arguments passes to the \code{plot} method of the
#'             \code{\link[pROC]{ci.se}} object of the \code{pROC} package.
#' @param no_numeric_policy. One of "plot_blank" or "warning" or "error". Defaults to "warning." What to do if \code{modelingResult} does not contain any numeric values, e.g., if \code{compute_rf} is set to \code{FALSE}. "plot_blank" will plot a blank ROC curve. "warning" will emit a warning but then will perform "plot_blank" behavior. "error" throw an error.
#'
#'
#' @importFrom ROCR performance
#' @importFrom graphics abline
#' @importFrom pROC roc ci.se
#'
#' @export plotAUC_gg
#'
plotAUC_gg <- function(
    modelingResult,
    CI = TRUE,
    rectilinear = FALSE,
    no_numeric_policy = c("warning", "plot_blank", "error"),
    seed = 0) {
   if (identical(no_numeric_policy, c("warning", "plot_blank", "error"))) {
      no_numeric_policy <- "warning"
   } else {
      if (!no_numeric_policy %in% c("warning", "plot_blank", "error")) {
         stop("no_numeric_policy must be one of warning, plot_blank, or error.")
      }
   }
   fpr <- NULL
   lo <- NULL
   hi <- NULL

   random_line_col <- "#888888"
   ci <- CI

   plot_blank <- function() {
      p <- ggplot2::ggplot() +
         ggplot2::geom_abline(
            mapping = aes(color = "", slope = 1, intercept = 0),
            linetype = "dashed"
         ) +
         ggplot2::scale_color_manual(values = c(random_line_col), name = "Random Model") +
         ggplot2::xlab("False Positive Rate") +
         ggplot2::ylab("True Positive Rate") +
         ggplot2::labs(
            title = "ROC Curve",
            subtitle = "AUC: NA"
         ) +
         ggplot2::theme(
            aspect.ratio = 1,
            plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
            plot.subtitle = element_text(hjust = 0.5, family = "mono")
         ) +
         ggplot2::scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
         ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2))
      return(p)
   }

   make_rectilinear_roc_ci <- function(ci_df) {
      intermediate_points_lo <- list()
      intermediate_points_hi <- list()
      prev_y_lo <- 0
      prev_y_hi <- 0
      for (i in 2:nrow(ci_df)) {
         mid_val <- (ci_df$fpr[i - 1] + ci_df$fpr[i]) / 2
         if (ci_df$lo[i] != prev_y_lo) {
            intermediate_points_lo <- append(
               intermediate_points_lo,
               list(
                  c(mid_val, ci_df$lo[i - 1]),
                  c(mid_val, ci_df$lo[i])
               )
            )
         }

         if (ci_df$hi[i] != prev_y_hi) {
            intermediate_points_hi <- append(
               intermediate_points_hi,
               list(
                  c(mid_val, ci_df$hi[i - 1]),
                  c(mid_val, ci_df$hi[i])
               )
            )
         }

         prev_y_lo <- ci_df$lo[i]
         prev_y_hi <- ci_df$hi[i]
      }

      lo_df <- as.data.frame(t(do.call(cbind, intermediate_points_lo))) %>%
         mutate(hi = NA) %>%
         `colnames<-`(c("fpr", "lo", "hi"))
      hi_df <- as.data.frame(t(do.call(cbind, intermediate_points_hi))) %>%
         mutate(lo = NA) %>%
         `colnames<-`(c("fpr", "hi", "lo"))
      hi_df <- hi_df[, c("fpr", "lo", "hi")]
      ci_df <- rbind(ci_df, lo_df, hi_df) %>% arrange(fpr)

      for (i in seq_len(nrow(ci_df))) {
         if (is.na(ci_df[i, "lo"])) {
            ci_df[i, "lo"] <- ci_df[i - 1, "lo"]
         }
         if (is.na(ci_df[i, "hi"])) {
            ci_df[i, "hi"] <- ci_df[i - 1, "hi"]
         }
      }

      return(ci_df)
   }

   plot_main <- function(seed = seed) {
      perf <- performance(modelingResult$pred, "tpr", "fpr")
      set.seed(seed)
      pROC_obj <- suppressMessages(pROC::roc(modelingResult$pred@labels[[1]],
         modelingResult$pred@predictions[[1]],
         direction = "<",
         smoothed = TRUE,
         # arguments for ci
         ci = TRUE,
         ci.alpha = 0.95,
         stratified = FALSE,
         # arguments for plot
         plot = FALSE,
         auc.polygon = FALSE,
         max.auc.polygon = FALSE,
         grid = FALSE
      ))
      df <- data.frame(x = perf@x.values[[1]], y = perf@y.values[[1]])
      the_ci <- pROC::ci.se(
         pROC_obj,
         specificities = seq(0, 1, 0.025),
         progress = "none"
      )
      the_ci <- as.data.frame(the_ci) %>%
         `colnames<-`(c("lo", "mid", "hi")) %>%
         tibble::rownames_to_column("fpr") %>%
         mutate(fpr = as.double(sub("X", "", fpr))) %>%
         select(-c("mid"))
      the_ci$fpr <- 1 - the_ci$fpr
      the_ci <- the_ci %>% arrange(fpr)
      if (rectilinear) {
         the_ci <- make_rectilinear_roc_ci(the_ci)
      }

      auc.ci.lo <- pROC_obj$ci[1]
      auc.ci.hi <- pROC_obj$ci[3]
      ribbon_col <- "#c96f6f"
      ci_col <- "red"
      p <- ggplot2::ggplot(data = df) +
         ggplot2::geom_ribbon(
            data = the_ci,
            mapping = aes(x = fpr, ymin = lo, ymax = hi, fill = ""),
            color = ci_col, alpha = 0.5
         ) +
         ggplot2::scale_fill_manual(values = c(ribbon_col), name = "95% CI") +
         ggplot2::geom_line(mapping = aes(x = x, y = y), lwd = 1) +
         ggplot2::geom_abline(
            mapping = aes(color = "", slope = 1, intercept = 0),
            linetype = "dashed"
         ) +
         ggplot2::scale_color_manual(values = c(random_line_col), name = "Random Model") +
         ggplot2::xlab("False Positive Rate") +
         ggplot2::ylab("True Positive Rate") +
         ggplot2::labs(
            title = "ROC Curve",
            subtitle = sprintf(
               "AUC: %.3f\nAUC CI: [%.3f, %.3f]",
               modelingResult$auc,
               auc.ci.lo,
               auc.ci.hi
            )
         ) +
         ggplot2::theme(
            aspect.ratio = 1,
            plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
            plot.subtitle = element_text(hjust = 0.5, family = "mono")
         ) +
         ggplot2::scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
         ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2))
      return(p)
   }


   if (is.na(modelingResult$auc)) {
      if (no_numeric_policy == "warning") {
         warning("modelingResult$auc is NA. Plotting blank ROC Curve.")
         return(plot_blank())
      } else {
         if (no_numeric_policy == "plot_blank") {
            return(plot_blank())
         } else {
            stop("modelingResult$auc is NA and `modelingResult` == \"error\".")
         }
      }
   }

   return(plot_main())
}
