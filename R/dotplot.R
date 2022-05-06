#' Title
#'
#' @param s1
#' @param s2
#' @param width
#' @param step
#'
#' @return
#' @export
#'
#' @examples
seqdotplot <- function(s1, s2, width = 10, step = 1, n_mismatches = 0, geom = "auto", threads = 1) {
    chunk_width <- 500
    str_break <- function(str) {
        width <- chunk_width
        substring(
            str,
            seq(1, nchar(str), width),
            seq(width, nchar(str) + width - 1, width)
        )
    }

    if (nchar(s1) <= chunk_width || threads == 1) {
        dotplot_data <- compute_dotplot_data(s1, s2, width, step, n_mismatches)
    } else {
        s1_chunked <- str_break(s1)
        compute_chunk <- function(s1, s2, width, step, n_mismatches) {
            compute_dotplot_data(s1, s2, width, step, n_mismatches)
        }
        dotplot_data <- parallel::mclapply(
            s1_chunked,
            compute_chunk,
            # args
            s2 = s2,
            width = width,
            step = step,
            n_mismatches = n_mismatches,
            mc.cores = threads
        )
        dotplot_data <- map2(
            dotplot_data,
            seq_along(dotplot_data),
            function(df, i) {
                df$x <- df$x + (i-1) * chunk_width
                df
            }
        )

        dotplot_data <- bind_rows(dotplot_data)
    }

    if (geom == "auto") {
        if (nchar(s1) > 2000 || nchar(s2) > 2000) {
            geom <- "point"
        } else {
            geom <- "tile"
        }
    }

    p <- ggplot(dotplot_data, aes(x = x, y = y))

    if (geom == "point") {
        p <- p + geom_point(shape = 15, size = 0.1)
    } else if (geom == "tile") {
        p <- p + geom_tile()
    }

    p + scale_x_continuous(expand = c(0, 0), position = "top") +
        scale_y_reverse(expand = c(0, 0)) +
        theme_classic() +
        theme(panel.border = element_rect(colour = "black", fill=NA)) +
        xlab("Sequence 1") +
        ylab("Sequence 2")
}
