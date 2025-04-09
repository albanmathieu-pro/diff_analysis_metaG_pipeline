#' Produce a PCA plot from count_matrix
#'
#' The PCA is produced using the count_matrix by default (can also use ruvg, combat, or clr normalization;
#' see the \code{use_normalisation} param).
#'
#' @param count_matrix A matrix of counts where rows are features (e.g., genes) and columns are samples.
#' @param use_normalisation What kind of normalisation should be used? Can be either "none", "ruvg", "combat", or "clr". Default: "none".
#' @param min_counts The minimal mean number of counts required to keep the feature for the PCA.
#' @param metadata The metadata to add to the coordinate table to add color or shape with the results and the \code{plot_pca} function.
#' If there is no metadata available, keep the default value of \code{NULL}. Value can either be a csv file or a \code{data.frame}. Default: \code{NULL}.
#' @param id_metadata The colname to use to join the metadata with the coordinate table. Must contain all the sample names found in the count_matrix.
#' If \code{NULL} and a metadata table is provided, the first column will be used to join the tables. Default: \code{NULL}.
#' @param ncp Number of components to include in the graph. Default: 2.
#'
#' @return Produce the PCA and silently returns a \code{list} with the coordinates and the x and y labels.
#'
#' @examples
#' count_matrix <- matrix(runif(100), nrow = 10, ncol = 10)
#' colnames(count_matrix) <- paste0("Sample", 1:10)
#' rownames(count_matrix) <- paste0("Gene", 1:10)
#' res_pca <- produce_pca_df(count_matrix)
#'
#' @export
produce_pca_df <- function(count_matrix, use_normalisation = "none", min_counts = 5,
                           metadata = NULL, id_metadata = NULL, ncp = 2) {

    stopifnot(is(count_matrix, "matrix"))
    stopifnot(use_normalisation %in% c("none", "ruvg", "combat", "clr"))
    stopifnot(is(ncp, "numeric"))
    stopifnot(identical(ncp, round(ncp)))
    stopifnot(ncp > 1)

    if (!is.null(metadata)) {
        id_metadata <- validate_metadata(metadata, id_metadata, count_matrix)
    }

    if (use_normalisation == "none") {
        tpm <- as.data.frame(count_matrix)
    } else if (use_normalisation == "ruvg") {
        stopifnot("ruvg_counts" %in% names(attributes(count_matrix)))
        tpm <- as.data.frame(attr(count_matrix, "ruvg_counts"))
    } else if (use_normalisation == "combat") {
        stopifnot("combat_counts" %in% names(attributes(count_matrix)))
        tpm <- as.data.frame(attr(count_matrix, "combat_counts"))
    } else if (use_normalisation == "clr") {
        stopifnot("clr" %in% names(attributes(count_matrix)))
        tpm <- as.data.frame(attr(count_matrix, "clr"))
    }

    tpm <- tpm %>%
             dplyr::mutate(feature_s = rownames(tpm)) %>%
             tidyr::gather(sample, tpm, -feature_s)

    min_tpm <- dplyr::group_by(tpm, feature_s) %>%
        dplyr::summarize(tpm = sum(tpm)) %>%
        dplyr::filter(tpm >= min_counts) %>%
        dplyr::pull(feature_s)

    tpm_filter <- dplyr::filter(tpm, feature_s %in% min_tpm)

    df <- tidyr::spread(tpm_filter, sample, tpm) %>%
        as.data.frame
    rownames(df) <- df$feature_s
    m <- dplyr::select(df, -feature_s) %>%
        as.matrix %>%
        t

    pca <- FactoMineR::PCA(m, graph = FALSE)
    coord <- pca$ind$coord
    dims <- stringr::str_extract(colnames(coord), "[0-9]*$") %>%
        as.numeric
    stopifnot(ncp <= tail(dims, 1))
    df <- data.frame(Dim1 = coord[,1])
    for (i in seq(2, ncp)) {
        df[[paste0("Dim", i)]] <- coord[,i]
    }
    df <- df %>%
        dplyr::mutate(sample = rownames(df))

    if (!is.null(metadata)) {
        if (is.null(id_metadata)) {
            id_metadata <- colnames(metadata)[1]
        }
        df <- dplyr::left_join(df, metadata, by = c("sample" = id_metadata))
    }

    xlab <- paste0("Dim1 (", pca$eig[1,2] %>% round(2), "%)")
    ylab <- paste0("Dim2 (", pca$eig[2,2] %>% round(2), "%)")

    list(coord = df, xlab = xlab, ylab = ylab)
}

validate_metadata <- function(metadata, id_metadata, count_matrix) {
    if (!is.null(metadata)) {
        stopifnot(is(metadata, "data.frame") | is(metadata, "character"))
        if (is.character(metadata)) {
            stopifnot(file.exists(metadata))
            metadata <- readr::read_csv(metadata, show_col_types = FALSE)
        }
        if (is.null(id_metadata)) {
            id_metadata <- colnames(metadata)[1]
        }
        stopifnot(is(id_metadata, "character"))
        stopifnot(id_metadata %in% colnames(metadata))
        stopifnot(all(colnames(count_matrix) %in% metadata[[id_metadata]]))
    }
    id_metadata
}

#' Produce a PCA plot from produce_pca_df results
#'
#' @param res_pca The \code{list} returned from the \code{produce_pca} function
#' @param size The size of the points. Default: 3.
#' @param color The name of the column in \code{res_pca$coord} to use to color the
#' points in the PCA. If \code{NULL}, won't add color. Default: \code{NULL}.
#' @param shape The name of the column in \code{res_pca$coord} to use to define
#' the shape of the points in the PCA. If \code{NULL}, won't add change the shape.
#' Default: \code{NULL}.
#' @param show_names Try to add the labels of each point? Will only add the
#' labels if the points are not too closely clustered together in the PCA.
#' Default: \code{TRUE}.
#' @param show_ellipses Add ellipses around groups? Default: \code{FALSE}.
#' @param title Add a title to the graph? If \code{NULL}, no title is added to
#' the graph. Otherwise this param must be a \code{character}. Default:
#' \code{NULL}.
#' @param graph Produce the graph. \code{TRUE} or \code{FALSE}. Default:
#' \code{TRUE}.
#' @param legend.position Value for the \code{legend.position} param from
#' ggplot2. Default: "right".
#' @param legend.box Value for the \code{legend.box} param from ggplot2.
#' Default: "vertical".
#' @param manual_colors A named vector of colors to manually set group colors. Default: \code{NULL}.
#'
#' @return Returns the \code{ggplot} object
#'
#' @examples
#' txi <- get_demo_txi()
#' res_pca <- produce_pca_df(txi)
#' p <- plot_pca(res_pca, graph = FALSE)
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 scale_shape_manual
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 ggtitle
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggnewscale new_scale_colour
#' @importFrom ggplot2 stat_ellipse
#'
#' @export
plot_pca <- function(res_pca, size = 3, color = NULL, shape = NULL,
                     show_names = TRUE, show_ellipses = FALSE, title = NULL, graph = TRUE,
                     legend.position = "right",
                     legend.box = "vertical",
                     manual_colors = NULL) {

    # Validate the params
    stopifnot(is(res_pca, "list"))
    stopifnot(all(c("coord", "xlab", "ylab") %in% names(res_pca)))
    stopifnot(is(res_pca$coord, "data.frame"))
    stopifnot(nrow(res_pca$coord) > 0)
    stopifnot(all(c("Dim1", "Dim2", "sample") %in% colnames(res_pca$coord)))

    stopifnot(is(res_pca$xlab, "character"))
    stopifnot(is(res_pca$ylab, "character"))

    stopifnot(is(size, "numeric"))
    stopifnot(identical(size, round(size)))
    stopifnot(size > 0)

    if (!is.null(color)) {
        stopifnot(is(color, "character"))
        stopifnot(color %in% colnames(res_pca$coord))
    }
    if (!is.null(shape)) {
        stopifnot(is(shape, "character"))
        stopifnot(shape %in% colnames(res_pca$coord))
    }
    stopifnot(is(show_names, "logical"))
    stopifnot(is(show_ellipses, "logical"))
    if (!is.null(title)) {
        stopifnot(is(title, "character"))
    }
    stopifnot(is(graph, "logical"))

    # To avoid 6 shapes limitation
    if (!is.null(shape)) {
        shape_vec <- c(15,16,17,18,19,20,1,2,3,4,5,6,7,8,9,10,11,12,13,14)
        shape_vec <- shape_vec[1:length(unique(res_pca$coord[[shape]]))]
    }

    # Produce the graph
    gg <- ggplot2::ggplot(res_pca$coord, ggplot2::aes(x = Dim1, y = Dim2)) +
        ggplot2::theme_bw() +
        ggplot2::labs(x = res_pca$xlab, y = res_pca$ylab)

    ## Add color and shape
    if (!is.null(color) & !is.null(shape)) {
        gg <- gg +
            ggplot2::geom_point(ggplot2::aes_string(col = color, shape = shape),
                                size = size) +
            ggplot2::scale_shape_manual(values = shape_vec)
    } else if (!is.null(color) & is.null(shape)) {
        gg <- gg +
            ggplot2::geom_point(ggplot2::aes_string(col = color), size = size)
    } else if (is.null(color) & !is.null(shape)) {
        gg <- gg +
            ggplot2::geom_point(ggplot2::aes_string(shape = shape), size = size) +
            ggplot2::scale_shape_manual(values = shape_vec)
    } else {
        gg <- gg + ggplot2::geom_point(size = size)
    }

    # Apply manual colors if provided
    if (!is.null(manual_colors) & !is.null(color)) {
        gg <- gg + ggplot2::scale_color_manual(values = manual_colors)
    }

    # Add ellipses if show_ellipses is TRUE
    if (show_ellipses) {
        gg <- gg +
            ggplot2::stat_ellipse(ggplot2::aes(color = group_val,fill = group_val), 
                                  geom = "polygon", alpha = 0.1 ) +
            ggplot2::scale_color_manual(values = manual_colors)+
            ggplot2::scale_fill_manual(values = manual_colors)
    }

    # legend
    if (!is.null(color) | !is.null(shape)) {
        gg <- gg + ggplot2::theme(legend.position = legend.position,
                                  legend.box = legend.box)
    }

    # show_names
    if (show_names) {
        gg <- gg + ggrepel::geom_text_repel(ggplot2::aes(label = sample),
                                            color = "black", force = 10)
    }

    # title
    if (!is.null(title)) {
        gg <- gg + ggplot2::ggtitle(title)
    }

    # graph
    if (graph) {
        print(gg)
    }

    invisible(gg)
}

#' Deprecated
#'
#' @param txi The \code{txi} object returned by the \code{import_kallisto}
#' function.
#' @param graph Produce the graph. \code{TRUE} or \code{FALSE}. Default:
#' \code{TRUE}.
#' @param use_ruvg Use RUVg normalization? Needs to be pre-computed using the
#' \code{ruvg_normalization} function. Default: \code{FALSE}.
#'
#' @export
produce_pca <- function(txi, graph = TRUE, use_ruvg = FALSE) {
    msg <- paste0("This function is now deprecated, please use the ",
                  "`produce_pca_df` and the `plot_pca` functions to produce ",
                  "your PCA.")
    message(msg)
}
