#' Plot CosMx cell polygons and gene transcripts with full + zoom panels
#'
#' @param polygons A data.frame containing polygon vertices. Must include columns
#'   \code{fov}, \code{x_global_px}, \code{y_global_px}, and \code{cellID} (or \code{cell_ID}).
#' @param tx A data.frame of transcript coordinates. Must include \code{target}, \code{x_global_px}, \code{y_global_px}.
#' @param genes Character vector of genes to plot.
#' @param colors Colors for each gene in \code{genes}.
#' @param celltype_df Optional. data.frame mapping cells to types. Must include \code{fov}, \code{cellID} (or \code{cell_ID}) and \code{celltype_col}.
#' @param celltype_col Column name in \code{celltype_df} storing cell types.
#' @param celltype_palette Named vector mapping cell types to colors.
#' @param ... Other parameters controlling zoom, scales, and point sizes.
#'
#' @return A list with two ggplot objects: \code{$full} and \code{$zoom}.
#' @export
plotCosMxGenes <- function(
    polygons,
    tx,
    genes,
    colors = NULL,
    
    # ---- cell type overlay ----
    celltype_df = NULL,                 # data.frame with fov + cellID (or cell_ID) + cell type column
    celltype_col = "cell_type",          # column name in celltype_df storing cell type
    celltype_palette = NULL,             # named vector: names=cell types, values=colors (RECOMMENDED)
    celltype_alpha_full = 0.35,          # FULL: fill alpha
    celltype_alpha_zoom = 1,             # ZOOM: (outline mode, not used; keep for compatibility)
    show_celltype_legend_full = TRUE,    # show legend on FULL
    show_celltype_legend_zoom = FALSE,   # show legend on ZOOM (usually FALSE to avoid duplication)
    # ---------------------------
    
    # ---- which layers to show ----
    show_genes_full = FALSE,             # FULL: show gene points? (you want FALSE)
    show_genes_zoom = TRUE,              # ZOOM: show gene points? (you want TRUE)
    full_celltype_mode = c("fill", "outline", "both"), # you want "fill"
    zoom_celltype_mode = c("outline", "fill", "both"), # you want "outline"
    # -----------------------------
    
    zoom.x = NULL,
    zoom.y = NULL,
    zoom.width.x = 1000,
    zoom.width.y = 1000,
    microns_per_pixel = 0.12028,
    scale_length_um_full = 50,
    scale_length_um_zoom = 25,
    label_size_full = 22,
    label_size_zoom = 25,
    point_size_full = 0.3,
    point_size_zoom = 2,
    point_alpha_full = 1,
    point_alpha_zoom = 1,
    polygon_linewidth_full = 0.3,
    polygon_linewidth_zoom = 1,
    show_zoom_rect = TRUE,
    zoom.rect.color = "white",
    zoom.rect.alpha = 0,
    zoom.rect.linewidth = 1,
    zoom.rect.linetype = "dashed"
) {
  # Dependencies
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Please install dplyr")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install ggplot2")
  if (!requireNamespace("rlang", quietly = TRUE)) stop("Please install rlang")
  
  full_celltype_mode <- match.arg(full_celltype_mode)
  zoom_celltype_mode <- match.arg(zoom_celltype_mode)
  
  # Helper: detect which column name to use for cell id
  detect_cellid_col <- function(df) {
    if ("cellID" %in% colnames(df)) return("cellID")
    if ("cell_ID" %in% colnames(df)) return("cell_ID")
    stop("Cannot find cell id column. Expect 'cellID' or 'cell_ID'.")
  }
  
  # Basic checks
  if (is.null(colors)) colors <- grDevices::rainbow(length(genes))
  if (length(colors) != length(genes)) stop("Length of colors must match length of genes")
  
  # Ensure required columns exist
  stopifnot(all(c("fov", "x_global_px", "y_global_px") %in% colnames(polygons)))
  stopifnot(all(c("x_global_px", "y_global_px") %in% colnames(tx)))
  
  poly_cellid_col <- detect_cellid_col(polygons)
  
  # --- Join cell type info into polygons (optional) ---
  celltype_palette_resolved <- NULL
  if (!is.null(celltype_df)) {
    if (!("fov" %in% colnames(celltype_df))) stop("celltype_df must contain column: fov")
    ct_cellid_col <- detect_cellid_col(celltype_df)
    if (!(celltype_col %in% colnames(celltype_df))) {
      stop(paste0("celltype_df must contain column: ", celltype_col))
    }
    
    # Standardize join keys to: fov + cellID
    polygons2 <- polygons
    if (poly_cellid_col != "cellID") {
      polygons2 <- dplyr::rename(polygons2, cellID = !!rlang::sym(poly_cellid_col))
    }
    
    celltype_df2 <- celltype_df %>%
      dplyr::select(
        fov,
        cellID = !!rlang::sym(ct_cellid_col),
        cell_type_joined = !!rlang::sym(celltype_col)
      ) %>%
      dplyr::distinct()
    
    polygons <- polygons2 %>%
      dplyr::left_join(celltype_df2, by = c("fov", "cellID"))
    
    # Resolve palette (named vector)
    if (is.null(celltype_palette)) {
      levs <- sort(unique(polygons$cell_type_joined))
      levs <- levs[!is.na(levs)]
      if (length(levs) > 0) {
        celltype_palette_resolved <- stats::setNames(grDevices::hcl.colors(length(levs), "Dark 3"), levs)
      } else {
        celltype_palette_resolved <- NULL
      }
    } else {
      # IMPORTANT: user-provided palette should be a NAMED vector
      if (is.null(names(celltype_palette)) || any(names(celltype_palette) == "")) {
        stop("celltype_palette must be a NAMED vector: c('TypeA'='#RRGGBB', 'TypeB'='#RRGGBB', ...)")
      }
      celltype_palette_resolved <- celltype_palette
    }
    
    # Warn if palette missing some levels
    types_in_plot <- unique(polygons$cell_type_joined)
    types_in_plot <- types_in_plot[!is.na(types_in_plot)]
    missing_types <- setdiff(types_in_plot, names(celltype_palette_resolved))
    if (length(missing_types) > 0) {
      warning(
        "celltype_palette is missing these levels: ",
        paste(missing_types, collapse = ", "),
        ". They will be shown as grey30."
      )
    }
  } else {
    # If no celltype_df, standardize cell id for grouping if needed
    if (poly_cellid_col != "cellID") {
      polygons <- dplyr::rename(polygons, cellID = !!rlang::sym(poly_cellid_col))
    }
  }
  
  # Filter transcripts for genes (CosMx tx file commonly uses 'target' for gene)
  if (!("target" %in% colnames(tx))) stop("tx must contain column: target (gene name)")
  tx.data_list <- lapply(genes, function(g) tx[tx$target == g, , drop = FALSE])
  names(tx.data_list) <- genes
  
  # Full-range bounds
  x_range <- range(polygons$x_global_px, na.rm = TRUE)
  y_range <- range(polygons$y_global_px, na.rm = TRUE)
  x_min <- x_range[1]; x_max <- x_range[2]
  y_min <- y_range[1]; y_max <- y_range[2]
  x_span <- x_max - x_min
  y_span <- y_max - y_min
  
  if (is.null(zoom.x)) zoom.x <- 0.5
  if (is.null(zoom.y)) zoom.y <- 0.5
  if (zoom.x < 0 || zoom.x > 1 || zoom.y < 0 || zoom.y > 1) stop("zoom.x and zoom.y must be in [0,1]")
  
  zoom.x.start <- x_min + zoom.x * x_span
  zoom.y.start <- y_min + zoom.y * y_span
  if (zoom.x.start + zoom.width.x > x_max) zoom.x.start <- x_max - zoom.width.x
  if (zoom.y.start + zoom.width.y > y_max) zoom.y.start <- y_max - zoom.width.y
  zoom.window <- c(zoom.x.start, zoom.x.start + zoom.width.x, zoom.y.start, zoom.y.start + zoom.width.y)
  
  # Scale bars (requires sputilsR::make_scale_bar)
  if (!exists("make_scale_bar", where = asNamespace("sputilsR"), inherits = FALSE) &&
      !("make_scale_bar" %in% getNamespaceExports("sputilsR"))) {
    stop("Cannot find sputilsR::make_scale_bar. Please ensure sputilsR is loaded and updated.")
  }
  
  scale.bar.full <- sputilsR::make_scale_bar(
    x_vals = polygons$x_global_px,
    y_vals = polygons$y_global_px,
    microns_per_pixel = microns_per_pixel,
    scale_length_um = scale_length_um_full,
    label_size = label_size_full
  )
  
  # ------------------------
  # FULL plot: celltype FILLED, no genes (per your request)
  # ------------------------
  p_full <- ggplot2::ggplot()
  
  if (!is.null(celltype_df)) {
    if (full_celltype_mode == "fill") {
      p_full <- p_full +
        ggplot2::geom_polygon(
          data = polygons,
          mapping = ggplot2::aes(
            x = x_global_px, y = y_global_px,
            group = interaction(fov, cellID),
            fill = cell_type_joined
          ),
          color = NA,
          linewidth = 0,
          alpha = 0.8
        ) +
        ggplot2::scale_fill_manual(values = celltype_palette_resolved, na.value = "grey30") +
        ggplot2::guides(
          fill = if (show_celltype_legend_full) ggplot2::guide_legend(override.aes = list(alpha = 1)) else "none"
        )
    } else if (full_celltype_mode == "outline") {
      p_full <- p_full +
        ggplot2::geom_polygon(
          data = polygons,
          mapping = ggplot2::aes(
            x = x_global_px, y = y_global_px,
            group = interaction(fov, cellID),
            color = cell_type_joined
          ),
          fill = NA,
          linewidth = polygon_linewidth_full,
          alpha = 1
        ) +
        ggplot2::scale_color_manual(values = celltype_palette_resolved, na.value = "grey30") +
        ggplot2::guides(color = if (show_celltype_legend_full) ggplot2::guide_legend(override.aes = list(linewidth = 3)) else "none")
    } else { # both
      p_full <- p_full +
        ggplot2::geom_polygon(
          data = polygons,
          mapping = ggplot2::aes(
            x = x_global_px, y = y_global_px,
            group = interaction(fov, cellID),
            fill = cell_type_joined,
            color = cell_type_joined
          ),
          linewidth = polygon_linewidth_full,
          alpha = celltype_alpha_full
        ) +
        ggplot2::scale_fill_manual(values = celltype_palette_resolved, na.value = "grey30") +
        ggplot2::scale_color_manual(values = celltype_palette_resolved, na.value = "grey30") +
        ggplot2::guides(
          color = "none",
          fill = if (show_celltype_legend_full) ggplot2::guide_legend(override.aes = list(alpha = 1)) else "none"
        )
    }
  } else {
    # fallback: draw polygons in grey
    p_full <- p_full +
      ggplot2::geom_polygon(
        data = polygons,
        mapping = ggplot2::aes(x = x_global_px, y = y_global_px, group = interaction(fov, cellID)),
        fill = NA,
        color = "darkgrey",
        linewidth = polygon_linewidth_full,
        alpha = 1
      )
  }
  
  # IMPORTANT: genes in FULL only if requested
  if (isTRUE(show_genes_full)) {
    for (i in seq_along(genes)) {
      p_full <- p_full +
        ggplot2::geom_point(
          data = tx.data_list[[i]],
          mapping = ggplot2::aes(x = x_global_px, y = y_global_px),
          color = colors[i],
          size = point_size_full,
          alpha = point_alpha_full
        )
    }
  }
  
  # Zoom rectangle on full plot
  if (isTRUE(show_zoom_rect)) {
    p_full <- p_full +
      ggplot2::annotate(
        "rect",
        xmin = zoom.window[1], xmax = zoom.window[2],
        ymin = zoom.window[3], ymax = zoom.window[4],
        color = zoom.rect.color,
        alpha = zoom.rect.alpha,
        linewidth = zoom.rect.linewidth,
        linetype = zoom.rect.linetype
      )
  }
  
  p_full <- p_full +
    scale.bar.full$bg + scale.bar.full$rect + scale.bar.full$label +
    ggplot2::coord_equal(expand = FALSE) +
    ggplot2::theme_void() +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "black")) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0))
  
  # ------------------------
  # ZOOM plot: celltype OUTLINE + genes (per your request)
  # ------------------------
  polygons_zoom <- polygons %>%
    dplyr::filter(
      x_global_px >= zoom.window[1], x_global_px <= zoom.window[2],
      y_global_px >= zoom.window[3], y_global_px <= zoom.window[4]
    )
  
  tx_zoom_list <- lapply(tx.data_list, function(df) {
    df %>%
      dplyr::filter(
        x_global_px >= zoom.window[1], x_global_px <= zoom.window[2],
        y_global_px >= zoom.window[3], y_global_px <= zoom.window[4]
      )
  })
  
  scale.bar.zoom <- sputilsR::make_scale_bar(
    x_vals = zoom.window[1:2],
    y_vals = zoom.window[3:4],
    microns_per_pixel = microns_per_pixel,
    scale_length_um = scale_length_um_zoom,
    label_size = label_size_zoom
  )
  
  p_zoom <- ggplot2::ggplot()
  
  if (!is.null(celltype_df)) {
    if (zoom_celltype_mode == "outline") {
      p_zoom <- p_zoom +
        ggplot2::geom_polygon(
          data = polygons_zoom,
          mapping = ggplot2::aes(
            x = x_global_px, y = y_global_px,
            group = interaction(fov, cellID),
            color = cell_type_joined
          ),
          fill = NA,
          linewidth = polygon_linewidth_zoom,
          alpha = 1
        ) +
        ggplot2::scale_color_manual(values = celltype_palette_resolved, na.value = "grey30") +
        ggplot2::guides(
          color = if (show_celltype_legend_zoom) ggplot2::guide_legend(override.aes = list(linewidth = 3)) else "none"
        )
    } else if (zoom_celltype_mode == "fill") {
      p_zoom <- p_zoom +
        ggplot2::geom_polygon(
          data = polygons_zoom,
          mapping = ggplot2::aes(
            x = x_global_px, y = y_global_px,
            group = interaction(fov, cellID),
            fill = cell_type_joined
          ),
          color = NA,
          linewidth = polygon_linewidth_zoom,
          alpha = celltype_alpha_full
        ) +
        ggplot2::scale_fill_manual(values = celltype_palette_resolved, na.value = "grey30") +
        ggplot2::guides(
          fill = if (show_celltype_legend_zoom) ggplot2::guide_legend(override.aes = list(alpha = 1)) else "none"
        )
    } else { # both
      p_zoom <- p_zoom +
        ggplot2::geom_polygon(
          data = polygons_zoom,
          mapping = ggplot2::aes(
            x = x_global_px, y = y_global_px,
            group = interaction(fov, cellID),
            fill = cell_type_joined,
            color = cell_type_joined
          ),
          linewidth = polygon_linewidth_zoom,
          alpha = celltype_alpha_full
        ) +
        ggplot2::scale_fill_manual(values = celltype_palette_resolved, na.value = "grey30") +
        ggplot2::scale_color_manual(values = celltype_palette_resolved, na.value = "grey30") +
        ggplot2::guides(
          color = "none",
          fill = if (show_celltype_legend_zoom) ggplot2::guide_legend(override.aes = list(alpha = 1)) else "none"
        )
    }
  } else {
    p_zoom <- p_zoom +
      ggplot2::geom_polygon(
        data = polygons_zoom,
        mapping = ggplot2::aes(x = x_global_px, y = y_global_px, group = interaction(fov, cellID)),
        fill = NA,
        color = "darkgrey",
        linewidth = polygon_linewidth_zoom,
        alpha = 1
      )
  }
  
  # Genes in ZOOM only if requested
  if (isTRUE(show_genes_zoom)) {
    for (i in seq_along(genes)) {
      p_zoom <- p_zoom +
        ggplot2::geom_point(
          data = tx_zoom_list[[i]],
          mapping = ggplot2::aes(x = x_global_px, y = y_global_px),
          color = colors[i],
          size = point_size_zoom,
          alpha = point_alpha_zoom
        )
    }
  }
  
  p_zoom <- p_zoom +
    scale.bar.zoom$bg + scale.bar.zoom$rect + scale.bar.zoom$label +
    ggplot2::coord_equal(
      xlim = zoom.window[1:2],
      ylim = zoom.window[3:4],
      expand = FALSE
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "black"))
  
  return(list(full = p_full, zoom = p_zoom))
}
