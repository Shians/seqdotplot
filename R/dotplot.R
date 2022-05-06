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
seqdotplot <- function(s1, s2, width = 10, step = 1, n_mismatches = 0, geom = "auto") {
    dotplot_data <- compute_dotplot_data(s1, s2, width, step, n_mismatches)

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
