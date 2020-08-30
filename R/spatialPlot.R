#' Plot cluster labels in spatial context
#'
#' After running \code{spatialCluster()} or \code{spatialEnhance()}, we can
#' visualize the cluster assignments of each spot using \code{clusterPlot()} or
#' \code{enhancePlot()}, respectively.
#' 
#' @param sce SingleCellExperiment containing cluster assignments in
#'   \code{sce$spatial.cluster}
#' @param platform Spatial sequencing platform. If "Visium", the hex spot layout
#'   will be used, otherwise square spots will be plotted.
#' 
#' @return Both functions return a \code{ggplot} object.
#' 
#' @examples
#' set.seed(149)
#' sce <- exampleSCE()
#' sce$spatial.cluster <- floor(runif(ncol(sce), 1, 4))
#' clusterPlot(sce)
#'
#' @name spatialPlot
NULL

palette <- c("#0173b2", "#de8f05", "#029e73", "#d55e00", "#cc78bc",
             "#ca9161", "#fbafe4", "#949494", "#ece133", "#56b4e9")

#' Plot spatial cluster assignments.
#'
#' @importFrom ggplot2 ggplot aes_ geom_polygon scale_fill_manual coord_equal labs theme_void
#'
#' @export
#' @rdname spatialPlot
clusterPlot <- function(sce, platform=c("Visium", "ST")) {
    # TODO: add user-specified palette
    # TODO: add platform/lattice to sce metadata instead of passing
    platform <- match.arg(platform)

    cdata <- data.frame(colData(sce))
    if (platform == "Visium") {
        vertices <- .make_hex_spots(cdata)
    } else {
        vertices <- .make_square_spots(cdata)
    }

    splot <- ggplot(data=vertices, 
                    aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~factor(fill))) + 
        geom_polygon() +
        # scale_fill_manual() +
        labs(fill="Cluster") +
        coord_equal() +
        theme_void()

    splot
}

#' Plot enhanced spatial cluster assignments.
#'
#' @importFrom ggplot2 ggplot aes_ geom_polygon scale_fill_manual coord_equal labs theme_void
#'
#' @export
#' @rdname spatialPlot
enhancePlot <- function(sce, platform=c("Visium", "ST")) {
    # TODO: add user-specified palette
    # TODO: add platform/lattice to sce metadata instead of passing
    platform <- match.arg(platform)

    cdata <- data.frame(colData(sce))
    ## TODO: add "enhanced" to metadata and put this logic in clusterPlot()
    if (platform == "Visium") {
        vertices <- .make_triangle_subspots(cdata)
    } else {
        vertices <- .make_square_spots(cdata, scale.factor=(1/3))
    }

    ## TODO: add color (edge color) parameter
    splot <- ggplot(data=vertices,
                    aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~factor(fill))) +
        geom_polygon() +
        # scale_fill_manual() +
        labs(fill="Cluster") +
        coord_equal() +
        theme_void()

    splot
}

#' Helper to extract x, y, fill ID
#' @keywords internal
.select_spot_positions <- function(cdata, x="col", y="row", fill.col="spatial.cluster") {
    spot_positions <- cdata[, c(x, y, fill.col)]
    colnames(spot_positions) <- c("x.pos", "y.pos", "fill")
    spot_positions$spot <- rownames(spot_positions)

    spot_positions
}

#' Compute vertex coordinates for each spot in frame of plot
#'
#' @param spot_positions Center for hex, top left for square
#' @param vertex_offsets Data frame of (x, y) offsets wrt spot position for each
#'   vertex of spot
#' 
#' @return Cartesian product of positions and offsets, with coordinates
#'   computed as (pos + offset)
#'
#' @keywords internal
.make_spot_vertices <- function(spot_positions, vertex_offsets) {
    spot_vertices <- merge(spot_positions, vertex_offsets)
    spot_vertices$x.vertex <- spot_vertices$x.pos + spot_vertices$x.offset
    spot_vertices$y.vertex <- spot_vertices$y.pos + spot_vertices$y.offset

    as.data.frame(spot_vertices)
}

#' Make vertices for each hex spot
#' 
#' @keywords internal
.make_hex_spots <- function(cdata) {
    ## R = circumradius, distance from center to vertex
    ## r = inradius, distance from center to edge midpoint
    r <- 1/2
    R <- (2 / sqrt(3)) * r

    spot_positions <- .select_spot_positions(cdata)
    spot_positions <- .adjust_hex_centers(spot_positions)

    ## vertices of each hex (with respect to center coordinates)
    ## start at top center, loop clockwise
    vertex_offsets <- data.frame(x.offset=c(0, r, r, 0, -r, -r),
                                 y.offset=c(-R, -R/2, R/2, R, R/2, -R/2))

    spot_vertices <- .make_spot_vertices(spot_positions, vertex_offsets)

    ## Flip to match image orientation
    spot_vertices$y.vertex <- -spot_vertices$y.vertex

    spot_vertices
}

#' Adjust hex spot positions so hexagons are adjacent to each other in plot
#'
#' Spots are regular hexagons with one unit of horizontal distance
#' between centers
#' 
#' @keywords internal
.adjust_hex_centers <- function(spot_positions) {
    ## R = circumradius, distance from center to vertex
    ## r = inradius, distance from center to edge midpoint
    r <- 1/2
    R <- (2 / sqrt(3)) * r
    
    ## Start at (1-indexed origin)
    spot_positions$x.pos <- spot_positions$x.pos - min(spot_positions$x.pos) + 1
    spot_positions$y.pos <- spot_positions$y.pos - min(spot_positions$y.pos) + 1
    
    ## Shift centers up so rows are adjacent
    spot_positions$y.pos <- spot_positions$y.pos * R * (3/2)
    
    ## Spot columns are offset by row
    ## (i.e. odd rows have odd numbered columns, even rows have even)
    ## Shift centers to the left so columns are adjacent (but hexes stay offset)
    spot_positions$x.pos <- (spot_positions$x.pos + 1) / 2
    
    spot_positions
}

#' Make vertices for each square spot
#'
#' Squares are simple, just mae a unit square at each array coordinate
#' 
#' @keywords internal
.make_square_spots <- function(cdata, fill.col="spatial.cluster", scale.factor=1) {
    spot_positions <- .select_spot_positions(cdata)
    
    vertex_offsets <- data.frame(x.offset=c(0, 1, 1, 0),
                                  y.offset=c(0, 0, 1, 1))
    vertex_offsets <- vertex_offsets * scale.factor

    .make_spot_vertices(spot_positions, vertex_offsets)
}

#' Helper to pull out subspot position columns
#' Probably redundant with select_spot_positions above, but we need subspot.idx
#' 
#' @keywords internal
.select_subspot_positions <- function(cdata, x="spot.col", y="spot.row", fill.col="spatial.cluster") {
    spot_positions <- cdata[, c(x, y, fill.col, "subspot.idx")]
    colnames(spot_positions) <- c("x.pos", "y.pos", "fill", "subspot.idx")
    spot_positions$spot <- rownames(spot_positions)
    
    spot_positions
}

#' Make vertices for each triangle subspot of a hex
#' 
#' @keywords internal
.make_triangle_subspots <- function(cdata, fill.col="spatial.cluster") {
    spot_positions <- .select_subspot_positions(cdata, x="spot.col", y="spot.row")
    spot_positions <- .adjust_hex_centers(spot_positions)
    
    ## R = circumradius, distance from center to vertex
    ## r = inradius, distance from center to edge midpoint
    r <- 1/2
    R <- (2 / sqrt(3)) * r
    
    ## Make lists of triangle vertices (with respect to hex center)
    ## subspot.idx is same ordering as `shift` in spatialEnhance
    ## that is, beginning in top right and proceeding clockwise, (1, 5, 3, 4, 6, 2)
    vertex_offsets <- do.call(rbind, list(
        data.frame(x.offset=c(0, 0, r), y.offset=c(0, -R, -R/2), subspot.idx=1),
        data.frame(x.offset=c(0, r, r), y.offset=c(0, -R/2, R/2), subspot.idx=5),
        data.frame(x.offset=c(0, r, 0), y.offset=c(0, R/2, R), subspot.idx=3),
        data.frame(x.offset=c(0, 0, -r), y.offset=c(0, R, R/2), subspot.idx=4),
        data.frame(x.offset=c(0, -r, -r), y.offset=c(0, R/2, -R/2), subspot.idx=6),
        data.frame(x.offset=c(0, -r, 0), y.offset=c(0, -R/2, -R), subspot.idx=2)
    ))
    
    ## note that instead of cartesian product, `merge()` does an outer join
    ## on subspot.idx here
    spot_vertices <- .make_spot_vertices(spot_positions, vertex_offsets)
    spot_vertices$y.vertex <- -spot_vertices$y.vertex
    
    spot_vertices
}