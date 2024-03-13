#' Spatial plotting functions
#'
#' @param color Optional hex code to set color of borders around spots. Set to
#'   \code{NA} to remove borders.
#' @param ... Additional arguments for \code{geom_polygon()}. \code{size}, to
#'   specify the linewidth of these borders, is likely the most useful.
#' @param platform Spatial sequencing platform. If "Visium", the hex spot layout
#'   will be used, otherwise square spots will be plotted.\cr
#'   NOTE: specifying this argument is only necessary if \code{sce} was not
#'   created by \code{spatialCluster()} or \code{spatialEnhance()}.
#' @param is.enhanced True if \code{sce} contains subspot-level data instead of
#'   spots. Spatial sequencing platform. If true, the respective subspot lattice
#'   for each platform will be plotted.\cr
#'   NOTE: specifying this argument is only necessary if \code{sce} was not
#'   created by \code{spatialCluster()} or \code{spatialEnhance()}.
#'
#' @keywords internal
#' @name spatialPlot
NULL

## Use Seaborn colorblind palette as default
palette <- c("#0173b2", "#de8f05", "#029e73", "#d55e00", "#cc78bc",
             "#ca9161", "#fbafe4", "#949494", "#ece133", "#56b4e9")

#' Plot spatial cluster assignments. 
#' 
#' @param sce SingleCellExperiment. If \code{fill} is specified and is a string,
#'   it must exist as a column in \code{colData(sce)}.
#' @param label Labels used to color each spot. May be the name of a column in
#'   \code{colData(sce)}, or a vector of discrete values.
#' @param palette Optional vector of hex codes to use for discrete spot values.
#' @inheritParams spatialPlot
#' 
#' @return Returns a ggplot object.
#' 
#' @examples
#' sce <- exampleSCE()
#' clusterPlot(sce)
#'
#' @family spatial plotting functions
#'
#' @importFrom ggplot2 ggplot aes_ geom_polygon scale_fill_manual coord_equal labs theme_void
#' @export
clusterPlot <- function(sce, label="spatial.cluster",
                        palette=NULL, color=NULL,
                        platform=NULL, is.enhanced=NULL,
                        ...) {
    
    if (is.null(platform))
        platform <- .bsData(sce, "platform", "Visium")
    if (is.null(is.enhanced))
        is.enhanced <- .bsData(sce, "is.enhanced", FALSE)
    
    vertices <- .make_vertices(sce, label, platform, is.enhanced)
    
    ## No borders around subspots by default
    if (is.null(color)) {
        color <- ifelse(is.enhanced, NA, "#d8dcd6")
    }

    splot <- ggplot(data=vertices, 
                    aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~factor(fill))) +
        geom_polygon(color=color, ...) +
        labs(fill="Cluster") +
        coord_equal() +
        theme_void()

    if (!is.null(palette))
        splot <- splot + scale_fill_manual(values=palette)

    splot
}

#' Plot spatial gene expression.
#' 
#' @param sce SingleCellExperiment. If \code{feature} is specified and is a 
#'   string, it must exist as a row in the specified assay of \code{sce}.
#' @param feature Feature vector used to color each spot. May be the name of a
#'   gene/row in an assay of \code{sce}, or a vector of continuous values.
#' @param assay.type String indicating which assay in \code{sce} the expression
#'   vector should be taken from.
#' @param low,mid,high Optional hex codes for low, mid, and high values of the
#'   color gradient used for continuous spot values.
#' @param diverging If true, use a diverging color gradient in
#'   \code{featurePlot()} (e.g. when plotting a fold change) instead of a
#'   sequential gradient (e.g. when plotting expression).
#' @inheritParams spatialPlot
#' 
#' @return Returns a ggplot object.
#' 
#' @examples
#' sce <- exampleSCE()
#' featurePlot(sce, "gene_2")
#' 
#' @family spatial plotting functions
#'
#' @importFrom ggplot2 ggplot aes_ geom_polygon scale_fill_gradient scale_fill_gradient2 coord_equal labs theme_void
#' @importFrom scales muted
#' @importFrom assertthat assert_that
#' @export
featurePlot <- function(sce, feature,
                        assay.type="logcounts", 
                        diverging=FALSE,
                        low=NULL, high=NULL, mid=NULL,
                        color=NULL,
                        platform=NULL, is.enhanced=NULL,
                        ...) {
    
    if (is.null(platform))
        platform <- .bsData(sce, "platform", "Visium")
    if (is.null(is.enhanced))
        is.enhanced <- .bsData(sce, "is.enhanced", FALSE)
    
    ## extract expression from logcounts if a gene name is passed.
    ## otherwise, assume a vector of counts was passed and let
    ## .make_vertices helpers check validity
    if (is.character(feature)) {
        assert_that(feature %in% rownames(sce),
                    msg=sprintf("Feature %s not in SCE.", feature))
        fill <- assay(sce, assay.type)[feature, ]
        fill.name <- feature
    } else {
        fill <- feature
        
        ## this could be an argument, but it's easily overwritten with labs()
        ## and we should encourage composing ggplot functions instead
        fill.name <- "Expression"
    }
    
    vertices <- .make_vertices(sce, fill, platform, is.enhanced)
    
    ## No borders around subspots by default
    if (is.null(color)) {
        color <- ifelse(is.enhanced, NA, "#d8dcd6")
    }
    
    splot <- ggplot(data=vertices, 
                    aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~fill)) +
        geom_polygon(color=color, ...) +
        labs(fill=fill.name) +
        coord_equal() +
        theme_void()
    
    if (diverging) {
        low <- ifelse(is.null(low), "#F0F0F0", low)
        high <- ifelse(is.null(high), muted("red"), high)
        splot <- splot + scale_fill_gradient(low=low, high=high)
    } else {
        low <- ifelse(is.null(low), muted("blue"), low)
        mid <- ifelse(is.null(mid), "#F0F0F0", mid)
        high <- ifelse(is.null(high), muted("red"), high)
        splot <- splot + scale_fill_gradient2(low=low, mid=mid, high=high)
    }
    
    splot
}



#' Make vertices outlining spots/subspots for geom_polygon()
#' 
#' @param sce SingleCellExperiment with row/col in colData
#' @param fill Name of a column in \code{colData(sce)} or a vector of values to
#'   use as fill for each spot
#' @param platform "Visium", "VisiumHD" or"ST", used to determine spot layout
#' @param is.enhanced If true, \code{sce} contains enhanced subspot data instead
#'   of spot-level expression. Used to determine spot layout.
#'   
#' @return Table of (x.pos, y.pos, spot, fill); where \code{spot} groups the
#'   vertices outlining the spot's border
#' 
#' @keywords internal
.make_vertices <- function(sce, fill, platform, is.enhanced) {
    cdata <- data.frame(colData(sce))
    
    if (platform == "Visium") {
        if (is.enhanced) {
            vertices <- .make_triangle_subspots(cdata, fill)
        } else {
            vertices <- .make_hex_spots(cdata, fill)
        }
    } else if (platform %in% c("VisiumHD", "ST")) {
        if (is.enhanced) {
            vertices <- .make_square_spots(cdata, fill, scale.factor = 1/3, offset = -1/6)
        } else {
            vertices <- .make_square_spots(cdata, fill)
        }
    } else {
        stop("Unsupported platform: \"", platform, "\". Cannot create spot layout.")
    }
    
    vertices
}

#' Helper to extract x, y, fill ID from colData
#' 
#' @return Dataframe of (x.pos, y.pos, fill) for each spot
#' 
#' @keywords internal
#' @importFrom assertthat assert_that
.select_spot_positions <- function(cdata, x="array_col", y="array_row", fill="spatial.cluster") {
    ## Provide either a column name or vector of labels/values
    assert_that(is.vector(fill) || is.character(fill) || is.factor(fill))
    
    ## I think this is the best way to check if something is a string
    if (is.character(fill) && length(fill) == 1) {
        spot_positions <- cdata[, c(x, y, fill)]
        colnames(spot_positions) <- c("x.pos", "y.pos", "fill")    
    } else if (is.vector(fill) || is.factor(fill)) {
        assert_that(nrow(cdata) == length(fill))
        spot_positions <- cdata[, c(x, y)]
        colnames(spot_positions) <- c("x.pos", "y.pos")    
        spot_positions$fill <- fill
    }
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
#' @return Table of (x.pos, y.pos, spot, fill); where \code{spot} groups the
#'   vertices outlining the spot's border
#' 
#' @keywords internal
.make_hex_spots <- function(cdata, fill) {
    ## R = circumradius, distance from center to vertex
    ## r = inradius, distance from center to edge midpoint
    r <- 1/2
    R <- (2 / sqrt(3)) * r

    spot_positions <- .select_spot_positions(cdata, fill=fill)
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
#' @return Shifted spot centers
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
#' Squares are simple, just make a unit square at each array coordinate
#' 
#' @return Table of (x.pos, y.pos, spot, fill); where \code{spot} groups the
#'   vertices outlining the spot's border
#' 
#' @keywords internal
.make_square_spots <- function(cdata, fill="spatial.cluster", scale.factor = 1, offset = 0) {
    spot_positions <- .select_spot_positions(cdata, fill=fill)
    
    vertex_offsets <- data.frame(x.offset=c(0, 1, 1, 0),
                                  y.offset=c(0, 0, 1, 1))
    vertex_offsets <- vertex_offsets * scale.factor + offset

    .make_spot_vertices(spot_positions, vertex_offsets)
}

#' Helper to pull out subspot position columns
#' Probably redundant with select_spot_positions above, but we need subspot.idx
#' 
#' @return Dataframe of (x.pos, y.pos, fill) for each spot
#' 
#' @keywords internal
.select_subspot_positions <- function(cdata, x="spot.col", y="spot.row", fill="spatial.cluster") {
    ## Provide either a column name or vector of labels/values
    assert_that(is.vector(fill) || is.character(fill) || is.factor(fill))
    
    if (is.character(fill) && length(fill) == 1) {
        spot_positions <- cdata[, c(x, y, "subspot.idx", fill)]
        colnames(spot_positions) <- c("x.pos", "y.pos", "subspot.idx", "fill")
    } else if (is.vector(fill) || is.factor(fill)) {
        assert_that(nrow(cdata) == length(fill))
        spot_positions <- cdata[, c(x, y, "subspot.idx")]
        colnames(spot_positions) <- c("x.pos", "y.pos", "subspot.idx")    
        spot_positions$fill <- fill
    }

    spot_positions$spot <- rownames(spot_positions)
    
    spot_positions
}

#' Make vertices for each triangle subspot of a hex
#' 
#' @return Table of (x.pos, y.pos, spot, fill); where \code{spot} groups the
#'   vertices outlining the spot's border
#'
#' @keywords internal
.make_triangle_subspots <- function(cdata, fill="spatial.cluster") {
    spot_positions <- .select_subspot_positions(cdata, x="spot.col", y="spot.row", fill=fill)
    spot_positions <- .adjust_hex_centers(spot_positions)
    
    ## R = circumradius, distance from center to vertex
    ## r = inradius, distance from center to edge midpoint
    r <- 1/2
    R <- (2 / sqrt(3)) * r
    
    ## Make lists of triangle vertices (with respect to hex center)
    ## subspot.idx is same ordering as `shift` in spatialEnhance
    ## that is, beginning in top right and proceeding clockwise, (1, 5, 3, 4, 6, 2)
    ## NOTE: however, we need to reflect over x-axis to match global negation of y-coordinate
    vertex_offsets <- do.call(rbind, list(
        data.frame(x.offset=c(0, 0, r), y.offset=c(0, -R, -R/2), subspot.idx=3),
        data.frame(x.offset=c(0, r, r), y.offset=c(0, -R/2, R/2), subspot.idx=5),
        data.frame(x.offset=c(0, r, 0), y.offset=c(0, R/2, R), subspot.idx=1),
        data.frame(x.offset=c(0, 0, -r), y.offset=c(0, R, R/2), subspot.idx=2),
        data.frame(x.offset=c(0, -r, -r), y.offset=c(0, R/2, -R/2), subspot.idx=6),
        data.frame(x.offset=c(0, -r, 0), y.offset=c(0, -R/2, -R), subspot.idx=4)
    ))
    
    ## note that instead of cartesian product, `merge()` does an outer join
    ## on subspot.idx here
    spot_vertices <- .make_spot_vertices(spot_positions, vertex_offsets)
    spot_vertices$y.vertex <- -spot_vertices$y.vertex
    
    spot_vertices
}
