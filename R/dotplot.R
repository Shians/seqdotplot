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
#' @param threads the number of threads to use.
#'
#' @return ggplot object containing the sequence dotplot
#' @export
#'
#' @examples
seqdotplot <- function(
    s1, s2,
    width = 10, step = 1, n_mismatches = 0,
    geom = c("auto", "point", "tile"),
    threads = getOption("mc.cores", 1L)
) {
    stopifnot(n_mismatches < width)
    geom <- match.arg(geom)

    xlab <- deparse(substitute(s1))
    ylab <- deparse(substitute(s2))

    s1 <- as.character(s1)
    s2 <- as.character(s2)

    s1_longer <- nchar(s1) > nchar(s2)
    if (!s1_longer) {
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

    if (nchar(s1) <= chunk_width || threads == 1) {
        dotplot_data <- compute_dotplot_data(s1, s2, width, step, n_mismatches)
    } else {
        compute_chunk <- function(s1, s2, width, step, n_mismatches) {
            compute_dotplot_data(s1, s2, width, step, n_mismatches)
        }

        chunked_str <- str_break(if (s1_longer) s1 else s2)

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

        dotplot_data_fwd <- map2(
            dotplot_data_fwd,
            seq_along(dotplot_data_fwd),
            function(df, i) {
                df$x <- df$x + (i-1) * chunk_width
                df
            }
        )

        dotplot_data_fwd <- bind_rows(dotplot_data_fwd) %>%
            mutate(orientation = "forward")

        s2 <- s2 %>%
            Biostrings::DNAString() %>%
            Biostrings::reverseComplement() %>%
            as.character

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

        dotplot_data_rev <- map2(
            dotplot_data_rev,
            seq_along(dotplot_data_rev),
            function(df, i) {
                df$x <- df$x + (i-1) * chunk_width
                df
            }
        )

        offset <- if (s1_longer) nchar(s1) else nchar(s2)
        dotplot_data_rev <- bind_rows(dotplot_data_rev) %>%
            mutate(x = offset - x) %>%
            mutate(orientation = "reverse")

        dotplot_data <- bind_rows(dotplot_data_fwd, dotplot_data_rev) %>%
            arrange(x)
    }

    if (geom == "auto") {
        if (nchar(s1) > 2000 || nchar(s2) > 2000) {
            geom <- "point"
        } else {
            geom <- "tile"
        }
    }

    p <- ggplot(dotplot_data, aes(x = x, y = y, col = orientation))

    if (geom == "point") {
        p <- p + geom_point(shape = 15, size = 0.1)
    } else if (geom == "tile") {
        p <- p + geom_tile()
    }

    p + scale_x_continuous(expand = c(0, 0), position = "top") +
        scale_y_reverse(expand = c(0, 0)) +
        theme_classic() +
        scale_colour_manual(values = c("forward" = "black", "reverse" = "red")) +
        guides(colour = guide_legend(override.aes = list(size=3))) +
        theme(panel.border = element_rect(colour = "black", fill=NA)) +
        xlab(xlab) +
        ylab(ylab)
}
