#' Plot QC Filters for Gene Expression Data
#'
#' This function creates a quality control (QC) plot for a `SummarizedExperiment` object, highlighting samples that pass or fail specified QC thresholds.
#'
#' @param gexp A `SummarizedExperiment` object containing gene expression data. The `colData` must include `nfeats` (number of detected features) and `ncounts` (total counts per sample).
#' @param fname A character string specifying the file path to save the plot as a PDF. Default is `"results/qc/plot_QC_filters.pdf"`. Set to `NULL` to skip saving.
#' @param min_nfeats An integer specifying the minimum number of features for a sample to pass QC. Default is `5000`.
#' @param min_ncounts A numeric value specifying the minimum total counts for a sample to pass QC. Default is `5e6`.
#' @param out A character string indicating the output type: `"plotly"` (interactive Plotly plot) or `"ggplot"` (static ggplot). Default is `"plotly"`.
#'
#' @return A QC plot, either as a `plotly` interactive object or a `ggplot` static object, depending on the `out` parameter.
#'
#' @details
#' The function identifies samples that fail QC based on thresholds for the number of detected features (`min_nfeats`) and total counts (`min_ncounts`). It creates a scatter plot with:
#' - The x-axis representing the number of features (`nfeats`).
#' - The y-axis representing the total counts (`ncounts`).
#' - Points color-coded based on QC status (`pass` or `fail`).
#' - Dashed lines indicating the QC thresholds.
#'
#' If `out = "plotly"`, the plot is interactive, allowing for sample tooltips. If a file name (`fname`) is provided, the static ggplot is saved as a PDF.
#'
#' @export
#'
#' @importFrom ggplot2 ggplot geom_vline geom_hline geom_point ggsave scale_x_continuous scale_y_continuous scale_color_manual labs expansion
#' @importFrom plotly ggplotly
#' @importFrom dplyr mutate case_when if_else
#' @importFrom ggrepel geom_text_repel
#' @importFrom SummarizedExperiment colData

plot_qc_filters <- function(gexp,
                            fname = "results/qc/plot_QC_filters.pdf",
                            min_nfeats = 5000, min_ncounts = 5e6,
                            out = c("plotly", "ggplot")[1]) {

  # Ensure gexp is a SummarizedExperiment object
  if (!"SummarizedExperiment" %in% class(gexp)) {
    stop("gexp must be a SummarizedExperiment object")
  }

  # Extract the QC data from colData
  qc_table <- SummarizedExperiment::colData(gexp) |>
    as.data.frame() |>
    dplyr::mutate(qc_status = dplyr::case_when(
      nfeats < min_nfeats ~ "QC fail",
      ncounts < min_ncounts ~ "QC fail",
      TRUE ~ "QC pass"
    ))

  # Create the plot
  plt_qc1 <- ggplot2::ggplot(qc_table, aes(nfeats, ncounts, label = sample_id)) +
    ggplot2::geom_vline(xintercept = min_nfeats, linetype = "dashed") +
    ggplot2::geom_hline(yintercept = min_ncounts, linetype = "dashed") +
    ggplot2::geom_point(aes(color = qc_status)) +
    ggrepel::geom_text_repel(aes(label = dplyr::if_else(qc_status == "QC fail", sample_id, NA)),
                    min.segment.length = 0.1, size = 3) +
    ggplot2::scale_x_continuous(limits = c(0, max(colData(gexp)$nfeats, 6000)),
                       expand = ggplot2::expansion(c(0,0.1))) +
    ggplot2::scale_y_continuous(limits = c(0, max(max(colData(gexp)$ncounts), 5.1e6)),
                       expand = ggplot2::expansion(c(0,0.1))) +
    ggplot2::scale_color_manual(values = c("QC pass" = "black", "QC fail" = "red")) +
    ggplot2::labs(x = "Number of genes with non-zero counts", y = "Total number of counts",
         color = "QC status")

  # Save the plot if fname is not null
  if (!is.null(fname)) {
    dir.create(dirname(fname), showWarnings = FALSE, recursive = TRUE)
    ggplot2::ggsave(fname, plt_qc1, create.dir = TRUE)
  }

  # Convert to plotly if desired
  plt_qc_int <- plotly::ggplotly(plt_qc1, tooltip = c("sample_id", "qc_status"))

  # Return the desired plot type
  if(out == "plotly") {
    return(plt_qc_int)
  } else {
    return(plt_qc1)
  }
}
