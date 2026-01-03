#' Create a ggplot2 scale bar grob
#'
#' Constructs a white-background, black-rectangle scale bar with a text label
#' (in µm or mm) that can be added to a `ggplot2` object using
#' `annotation_custom()`. The bar length can be given explicitly in microns
#' or automatically chosen based on the x-range of the data.
#'
#' @param x_vals Numeric vector of x coordinates (in pixels) used to determine
#'   the horizontal extent and cropping of the scale bar.
#' @param y_vals Numeric vector of y coordinates (in pixels) used to determine
#'   the vertical extent and cropping of the scale bar.
#' @param microns_per_pixel Numeric scalar giving the conversion factor from
#'   pixels to microns. Defaults to 0.12028.
#' @param scale_length_um Optional numeric scalar: desired scale bar length in
#'   microns. If `NULL`, a "nice" value is chosen automatically based on the
#'   x-range (approximately one quarter of the field of view).
#' @param x_low_cutoff,x_high_cutoff Numeric in \[0, 1\]; quantiles used to
#'   clamp extreme values of `x_vals` before computing the range.
#' @param y_low_cutoff,y_high_cutoff Numeric in \[0, 1\]; quantiles used to
#'   clamp extreme values of `y_vals` before computing the range.
#' @param position Character string giving the position of the scale bar
#'   within the plot area. One of `"bottom-right"`, `"bottom-left"`,
#'   `"top-right"`, or `"top-left"`.
#' @param label_size Numeric scalar giving the font size (points) of the scale
#'   bar label text.
#'
#' @importFrom ggplot2 annotation_custom
#' @importFrom grid rectGrob textGrob gpar
#' @export
#' @return A named list of three `annotation_custom()` layers:
#'   - `bg`: white background rectangle behind the bar and label
#'   - `rect`: black rectangle representing the scale bar
#'   - `label`: text label showing the bar length (e.g. `"50 µm"` or `"1 mm"`)
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#'
#' scale_bar <- make_scale_bar(
#'     x_vals = cell_meta$CenterX_global_px,
#'     y_vals = cell_meta$CenterY_global_px
#' )
#'
#' ggplot() +
#'     scale_bar$bg +
#'     scale_bar$rect +
#'     scale_bar$label
#' }
#'
make_scale_bar <- function(
    x_vals,
    y_vals,
    microns_per_pixel = 0.12028,
    scale_length_um = NULL,
    x_low_cutoff = 0,
    x_high_cutoff = 1,
    y_low_cutoff = 0,
    y_high_cutoff = 1,
    position = c("bottom-right", "bottom-left", "top-right", "top-left"),
    label_size = 15) {
    # Resolve position
    position <- match.arg(position)

    # Quantile cutoff
    x_low_cutoff <- stats::quantile(x_vals, x_low_cutoff, na.rm = TRUE)
    x_high_cutoff <- stats::quantile(x_vals, x_high_cutoff, na.rm = TRUE)
    y_low_cutoff <- stats::quantile(y_vals, y_low_cutoff, na.rm = TRUE)
    y_high_cutoff <- stats::quantile(y_vals, y_high_cutoff, na.rm = TRUE)
    x_vals[x_vals < x_low_cutoff] <- x_low_cutoff
    x_vals[x_vals > x_high_cutoff] <- x_high_cutoff
    y_vals[y_vals < y_low_cutoff] <- y_low_cutoff
    y_vals[y_vals > y_high_cutoff] <- y_high_cutoff

    # Ranges
    x_range <- range(x_vals, na.rm = TRUE)
    y_range <- range(y_vals, na.rm = TRUE)
    x_length <- diff(x_range)
    x_length_um <- x_length * microns_per_pixel

    # Auto-determine scale length if not provided
    if (is.null(scale_length_um)) {
        target <- x_length_um / 4
        order <- 10^floor(log10(target))
        mantissa <- target / order
        nice_mantissa <- if (mantissa < 1.5) {
            1
        } else if (mantissa < 3.5) {
            2
        } else if (mantissa < 7.5) {
            5
        } else {
            10
        }
        scale_length_um <- nice_mantissa * order
    }

    # Convert to pixels
    scale_length_px <- scale_length_um / microns_per_pixel

    # Format label
    scale_label <- if (scale_length_um >= 1000) {
        paste0(scale_length_um / 1000, " mm")
    } else {
        paste0(scale_length_um, " µm")
    }

    # Margin proportional to bar height/width for clean padding
    pad_x <- scale_length_px * 0.1
    pad_y <- scale_length_px * 0.1
    bar_height_px <- scale_length_px * 0.05
    label_offset_y <- bar_height_px * 1.5

    # Determine positions
    if (position %in% c("bottom-right", "top-right")) {
        x_end <- x_range[2] - pad_x
        x_start <- x_end - scale_length_px
    } else {
        x_start <- x_range[1] + pad_x
        x_end <- x_start + scale_length_px
    }

    if (position %in% c("bottom-right", "bottom-left")) {
        y_pos <- y_range[1] + pad_y
    } else {
        y_pos <- y_range[2] - pad_y - bar_height_px
    }

    # Background rectangle
    bg_rect <- ggplot2::annotation_custom(
        grob = grid::rectGrob(gp = grid::gpar(fill = "white", alpha = 0.8, col = NA)),
        xmin = x_start - pad_x * 0.5,
        xmax = x_end + pad_x * 0.5,
        ymin = y_pos - pad_y * 0.5,
        ymax = y_pos + bar_height_px + label_offset_y * 5 + pad_y * 0.5
    )

    # Black scale bar
    bar_rect <- ggplot2::annotation_custom(
        grob = grid::rectGrob(gp = grid::gpar(fill = "black")),
        xmin = x_start,
        xmax = x_end,
        ymin = y_pos,
        ymax = y_pos + bar_height_px
    )

    # Label — uses label_size
    label_text <- ggplot2::annotation_custom(
        grob = grid::textGrob(
            scale_label,
            gp = grid::gpar(col = "black", fontsize = label_size),
            just = "center",
            vjust = 0
        ),
        xmin = (x_start + x_end) / 2,
        xmax = (x_start + x_end) / 2,
        ymin = y_pos + bar_height_px + label_offset_y,
        ymax = y_pos + bar_height_px + label_offset_y
    )

    list(bg = bg_rect, rect = bar_rect, label = label_text)
}
