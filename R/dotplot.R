#' Draw sequence dotplot
#'
#' Draws a sequence dotplot with forward and reverse complement matching. The
#' longer sequence will be represented on the X-axis.
#'
#' @param s1 the first sequence.
#' @param s2 the second sequence.
#' @param width the width of the window for comparing matches.
#' @param step the step size to move the window forward.
#' @param n_mismatches the allowed number of mismatches within window.
#' @param geom the geom used to draw the dots.
#' @param longer_as_x if TRUE, the longer sequence will be represented on the
#'  x-axis. If FALSE, the first sequence will be represented on the x-axis.
#' @param threads the number of threads to use.
#'
#' @return ggplot2::ggplot object containing the sequence dotplot
#' @export
#'
#' @examples
seqdotplot <- function(
    s1, s2,
    width = 10, step = 1, n_mismatches = 0,
    geom = c("auto", "point", "tile"),
    longer_as_x = FALSE,
    match_reverse = TRUE,
    threads = getOption("mc.cores", 1L)
) {
    stopifnot(n_mismatches < width)
    geom <- match.arg(geom)

    xlab <- deparse(substitute(s1))
    ylab <- deparse(substitute(s2))

    s1 <- as.character(s1)
    s2 <- as.character(s2)

    s1_longer <- nchar(s1) > nchar(s2)
    if (longer_as_x && !s1_longer) {
        temp <- xlab
        xlab <- ylab
        ylab <- temp
    }

    chunk_width <- 1000
    str_break <- function(str) {
        width <- chunk_width
        substring(
            str,
            seq(1, nchar(str), width),
            seq(width, nchar(str) + width - 1, width)
        )
    }

    if (threads == 1 || (nchar(s1) <= chunk_width && nchar(s2) <= chunk_width)) {
        # if single threaded or both strings are short
        dotplot_data_fwd <- compute_dotplot_data(s1, s2, width, step, n_mismatches)
        if (match_reverse) {
            s2 <- s2 |>
                Biostrings::DNAString() |>
                Biostrings::reverseComplement() |>
                as.character()
            dotplot_data_rev <- compute_dotplot_data(s1, s2, width, step, n_mismatches)
            dotplot_data <- dplyr::bind_rows(dotplot_data_fwd, dotplot_data_rev) |>
                dplyr::arrange(.data$x)
        } else {
            dotplot_data <- dotplot_data_fwd
        }
        
    } else {
        # otherwise, chunk the longer string and compute dotplot data in parallel
        compute_chunk <- function(s1, s2, width, step, n_mismatches) {
            compute_dotplot_data(s1, s2, width, step, n_mismatches)
        }

        # chunk along the longer string
        chunked_str <- str_break(if (s1_longer) s1 else s2)

        # compute dotplot data for forward matches
        dotplot_data_fwd <- parallel::mclapply(
            chunked_str,
            compute_chunk,
            # args
            s2 = if (s1_longer) s2 else s1,
            width = width,
            step = step,
            n_mismatches = n_mismatches,
            mc.cores = threads
        )

        dotplot_data_fwd <- purrr::map2(
            dotplot_data_fwd,
            seq_along(dotplot_data_fwd),
            function(df, i) {
                df$x <- df$x + (i-1) * chunk_width
                df
            }
        )

        dotplot_data_fwd <- dplyr::bind_rows(dotplot_data_fwd) |>
            dplyr::mutate(orientation = "forward")

        if (match_reverse) {
            # compute dotplot data for reverse complement matches
            s2 <- s2 |>
                Biostrings::DNAString() |>
                Biostrings::reverseComplement() |>
                as.character()

            # chunk along the longer string
            chunked_str <- str_break(if (s1_longer) s1 else s2)
            
            dotplot_data_rev <- parallel::mclapply(
                chunked_str,
                compute_chunk,
                # args
                s2 = if (s1_longer) s2 else s1,
                width = width,
                step = step,
                n_mismatches = n_mismatches,
                mc.cores = threads
            )

            dotplot_data_rev <- purrr::map2(
                dotplot_data_rev,
                seq_along(dotplot_data_rev),
                function(df, i) {
                    df$x <- df$x + (i-1) * chunk_width
                    df
                }
            )

            offset <- if (s1_longer) nchar(s1) else nchar(s2)
            dotplot_data_rev <- dplyr::bind_rows(dotplot_data_rev) |>
                dplyr::mutate(x = offset - x) |>
                dplyr::mutate(orientation = "reverse")

            dotplot_data <- dplyr::bind_rows(dotplot_data_fwd, dotplot_data_rev) |>
            dplyr::arrange(.data$x)
        } else {
            dotplot_data <- dotplot_data_fwd
        }
    }

    if (geom == "auto") {
        if (nchar(s1) > 2000 || nchar(s2) > 2000) {
            geom <- "point"
        } else {
            geom <- "tile"
        }
    }

    p <- ggplot2::ggplot(dotplot_data, ggplot2::aes(x = .data$x, y = .data$y, col = orientation))

    if (geom == "point") {
        p <- p + ggplot2::geom_point(shape = 15, size = 0.1)
    } else if (geom == "tile") {
        p <- p + ggplot2::geom_tile()
    }

    if (match_reverse) {
        p + ggplot2::scale_x_continuous(expand = c(0, 0), position = "top") +
            ggplot2::scale_y_reverse(expand = c(0, 0)) +
            ggplot2::theme_classic() +
            ggplot2::scale_colour_manual(values = c("forward" = "black", "reverse" = "red")) +
            ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=3))) +
            ggplot2::theme(panel.border = ggplot2::element_rect(colour = "black", fill=NA)) +
            xlab(xlab) +
            ylab(ylab)
    } else {
        p + ggplot2::scale_x_continuous(expand = c(0, 0), position = "top") +
            ggplot2::scale_y_reverse(expand = c(0, 0)) +
            ggplot2::theme_classic() +
            ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=3))) +
            ggplot2::theme(panel.border = ggplot2::element_rect(colour = "black", fill=NA)) +
            xlab(xlab) +
            ylab(ylab)
    }
}
